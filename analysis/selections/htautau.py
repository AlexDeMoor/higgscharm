"""Helpers for the H->tautau + charm selection (mu tau_h, e tau_h, tau_h tau_h).

Everything here returns awkward arrays; the ObjectSelector methods in
object_selections.py are thin wrappers that store the results in `objects`.
"""
import numpy as np
import awkward as ak
from coffea.nanoevents.methods import candidate

# gen-level decay channel codes (see gen_htautau_channel)
GEN_CHANNEL = {"none": 0, "mutau": 1, "etau": 2, "tautau": 3, "leptonic": 4}

# pair-ranking recipes: (key, ascending) from the most to the least significant
# TODO (mu tau_h): RANK_MUTAU, the ranking of the (muon, tau_h) candidates
RANK_ETAU = [
    (lambda p: p.l1.pfRelIso03_all, True),
    (lambda p: p.l1.pt, False),
    (lambda p: p.l2.rawDeepTau2018v2p5VSjet, False),
    (lambda p: p.l2.pt, False),
]
RANK_TAUTAU = [
    (lambda p: p.l1.rawDeepTau2018v2p5VSjet, False),
    (lambda p: p.l1.pt, False),
    (lambda p: p.l2.rawDeepTau2018v2p5VSjet, False),
    (lambda p: p.l2.pt, False),
]


def pt_sorted(collection):
    return collection[ak.argsort(collection.pt, ascending=False, axis=1)]


def to_candidates(collection, extra_fields=()):
    """re-zip a collection as PtEtaPhiMCandidate keeping only the listed extra fields"""
    fields = {
        "pt": collection.pt,
        "eta": collection.eta,
        "phi": collection.phi,
        "mass": collection.mass,
    }
    for f in extra_fields:
        if f in collection.fields:
            fields[f] = collection[f]
    return ak.zip(fields, with_name="PtEtaPhiMCandidate", behavior=candidate.behavior)


# ----------------------------------------------------------------------------------
# generator level
# ----------------------------------------------------------------------------------
def gen_htautau_channel(events):
    """Classify the tau-pair decay from GenPart: count electrons/muons whose DIRECT
    mother is a tau (GenVisTau is NOT used: it has a pT > 10 GeV builder cut).

    Returns (channel, gen_leptons, hard_taus) with channel: 1 mu tau_h, 2 e tau_h,
    3 tau_h tau_h, 4 fully leptonic (ee, emu, mumu), 0 none (fewer than two hard-process taus).
    """
    gp = events.GenPart
    parent_pdg = ak.fill_none(gp.parent.pdgId, 0)
    is_lep = (np.abs(gp.pdgId) == 11) | (np.abs(gp.pdgId) == 13)
    # only the decay products of the HARD-PROCESS taus (H -> tautau, Z -> tautau): taus from
    # B/D-hadron decays in the bH / cH samples would otherwise contaminate the classification
    direct = gp.hasFlags(["isDirectHardProcessTauDecayProduct"])
    leptons_from_tau = gp[is_lep & (np.abs(parent_pdg) == 15) & direct]
    n_mu = ak.sum(np.abs(leptons_from_tau.pdgId) == 13, axis=1)
    n_e = ak.sum(np.abs(leptons_from_tau.pdgId) == 11, axis=1)
    hard_taus = gp[(np.abs(gp.pdgId) == 15) & gp.hasFlags(["isHardProcess"])]
    n_tau = ak.num(hard_taus)
    channel = ak.where(
        n_tau < 2,
        GEN_CHANNEL["none"],
        ak.where(
            (n_mu == 0) & (n_e == 0),
            GEN_CHANNEL["tautau"],
            ak.where(
                (n_mu == 1) & (n_e == 0),
                GEN_CHANNEL["mutau"],
                ak.where((n_e == 1) & (n_mu == 0), GEN_CHANNEL["etau"], GEN_CHANNEL["leptonic"]),
            ),
        ),
    )
    return ak.values_astype(channel, np.int32), leptons_from_tau, hard_taus


# ----------------------------------------------------------------------------------
# di-tau pair
# ----------------------------------------------------------------------------------
def build_htt_pair(l1coll, l2coll, same_flavour, rank, min_dr=0.5):
    """All (l1, l2) candidates with dR(l1, l2) > min_dr, ranked with successive stable
    argsorts (least significant key first => lexicographic order)."""
    if same_flavour:
        pairs = ak.combinations(l2coll, 2, fields=["l1", "l2"])
    else:
        pairs = ak.cartesian({"l1": l1coll, "l2": l2coll}, axis=1)
    pairs = pairs[pairs.l1.delta_r(pairs.l2) > min_dr]
    for key, ascending in reversed(rank):
        order = ak.argsort(key(pairs), axis=1, ascending=ascending, stable=True)
        pairs = pairs[order]
    return pairs


def transverse_mass(a, b):
    return np.sqrt(2.0 * a.pt * b.pt * (1.0 - np.cos(a.delta_phi(b))))


def dzeta(l1, l2, met, alpha=0.85):
    """D_zeta = p_zeta^miss - alpha * p_zeta^vis (CMS H->tautau convention)"""
    zx = np.cos(l1.phi) + np.cos(l2.phi)
    zy = np.sin(l1.phi) + np.sin(l2.phi)
    zn = np.sqrt(zx**2 + zy**2)
    zx, zy = zx / zn, zy / zn
    vis_px = l1.pt * np.cos(l1.phi) + l2.pt * np.cos(l2.phi)
    vis_py = l1.pt * np.sin(l1.phi) + l2.pt * np.sin(l2.phi)
    pz_vis = vis_px * zx + vis_py * zy
    pz_miss = met.pt * np.cos(met.phi) * zx + met.pt * np.sin(met.phi) * zy
    return pz_miss - alpha * pz_vis


def collinear_mass(l1, l2, met):
    """collinear-approximation di-tau mass (NaN when unphysical) and the validity flag"""
    mx = met.pt * np.cos(met.phi)
    my = met.pt * np.sin(met.phi)
    p1x, p1y = l1.pt * np.cos(l1.phi), l1.pt * np.sin(l1.phi)
    p2x, p2y = l2.pt * np.cos(l2.phi), l2.pt * np.sin(l2.phi)
    det = p1x * p2y - p1y * p2x
    det = ak.where(det == 0, np.nan, det)
    a = (mx * p2y - my * p2x) / det
    b = (p1x * my - p1y * mx) / det
    x1, x2 = 1.0 / (1.0 + a), 1.0 / (1.0 + b)
    mvis = (l1 + l2).mass
    physical = ak.fill_none((x1 > 0) & (x1 <= 1) & (x2 > 0) & (x2 <= 1), False)
    arg = ak.where(physical, x1 * x2, 1.0)
    return ak.where(physical, mvis / np.sqrt(arg), np.nan), physical
