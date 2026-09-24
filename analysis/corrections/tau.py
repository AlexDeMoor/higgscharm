"""Tau corrections (TAU POG / CAT correctionlib, DeepTau2018v2p5):

* TauWeights            -- tau ID scale factors (VSjet, VSe, VSmu) and tau-leg trigger
                           scale factors as event weights
* apply_tau_es_corrections -- tau energy scale (object correction), propagated to MET

The tau ID and trigger SFs are evaluated on the tau_h legs of the selected di-tau pair
(`events.selected_htt_taus`), not on every tau in the event.

genmatch convention (Tau_genPartFlav): 1 prompt e, 2 prompt mu, 3 tau->e, 4 tau->mu,
5 genuine tau_h, 0 jet/unmatched.
"""
import copy
import warnings
import numpy as np
import awkward as ak
import correctionlib
from typing import Type
from coffea.analysis_tools import Weights
from analysis.corrections.utils import unflat_sf
from analysis.filesets.utils import get_nano_version
from analysis.corrections.correctionlib_files import correction_files, load_correction_set

ALLOWED_DECAY_MODES = (0, 1, 10, 11)
# lower pT edge of the 2024 `tau_trigger` fits (the same for every trigtype); the correction
# has flow=error, so a leg below it raises instead of being clamped
TAU_TRIGGER_PT_MIN = 24.6
_clamp_warned = set()


def _dm_mask(dm):
    mask = ak.zeros_like(dm, dtype=bool)
    for allowed in ALLOWED_DECAY_MODES:
        mask = mask | (dm == allowed)
    return mask


class TauWeights:
    """
    Tau ID scale factors

    Parameters:
    -----------
        events:
            events collection (must carry `selected_htt_taus`)
        weights:
            Weights container from coffea.analysis_tools
        year:
            dataset year {2024}
        shift:
            if None add 'nominal', 'up' and 'down' variations, else only 'nominal'
    """

    def __init__(
        self,
        events: ak.Array,
        weights: Type[Weights],
        year: str,
        shift: str,
        taus_field: str = "selected_htt_taus",
    ) -> None:
        self.events = events
        self.taus = events[taus_field]
        self.weights = weights
        self.year = year
        self.shift = shift
        self.nano_version = get_nano_version(year)

        self.flat_taus = ak.flatten(self.taus)
        self.taus_counts = ak.num(self.taus)
        self.cset = load_correction_set(correction_files["tau"][year])

    # ---- helpers --------------------------------------------------------------
    def _inputs(self):
        t = self.flat_taus
        pt = ak.to_numpy(ak.fill_none(t.pt, 30.0)).astype(np.float64)
        eta = np.abs(ak.to_numpy(ak.fill_none(t.eta, 0.0)).astype(np.float64))
        dm = ak.to_numpy(ak.fill_none(t.decayMode, 0)).astype(np.int64)
        gm = ak.to_numpy(ak.fill_none(t.genPartFlav, 0)).astype(np.int64)
        return pt, eta, dm, gm

    def _evaluate(self, name, mask, args, syst):
        # the ID corrections clamp inputs outside their domain (tau_trigger does NOT, see
        # add_trigger_weights); objects outside `mask` get SF = 1 through unflat_sf
        sf = self.cset[name].evaluate(*args)
        return unflat_sf(ak.Array(sf), ak.Array(mask), self.taus_counts)

    def get_vsjet_weights(self, wp_vsjet, wp_vse, syst):
        pt, eta, dm, gm = self._inputs()
        mask = (gm == 5) & _dm_mask(dm) & (pt > 20.0)
        dm_safe = np.where(_dm_mask(dm), dm, 0)
        gm_safe = np.where(mask, gm, 5)
        return self._evaluate(
            "DeepTau2018v2p5VSjet",
            mask,
            [pt, dm_safe, gm_safe, wp_vsjet, wp_vse, syst, "dm"],
            syst,
        )

    def get_vse_weights(self, wp_vse, syst):
        pt, eta, dm, gm = self._inputs()
        mask = ((gm == 1) | (gm == 3)) & _dm_mask(dm)
        dm_safe = np.where(_dm_mask(dm), dm, 0)
        gm_safe = np.where(mask, gm, 1)
        return self._evaluate(
            "DeepTau2018v2p5VSe", mask, [eta, dm_safe, gm_safe, wp_vse, syst], syst
        )

    def get_vsmu_weights(self, wp_vsmu, wp_vse, wp_vsjet, syst):
        pt, eta, dm, gm = self._inputs()
        mask = (gm == 2) | (gm == 4)
        gm_safe = np.where(mask, gm, 2)
        return self._evaluate(
            "DeepTau2018v2p5VSmu",
            mask,
            [eta, gm_safe, wp_vsmu, wp_vse, wp_vsjet, syst],
            syst,
        )

    def _add(self, name, func, args):
        nominal = func(*args, syst="nom")
        if self.shift is None:
            up = func(*args, syst="up")
            down = func(*args, syst="down")
            self.weights.add(name=name, weight=nominal, weightUp=up, weightDown=down)
        else:
            self.weights.add(name=name, weight=nominal)

    # ---- public ---------------------------------------------------------------
    def add_id_weights(self, vsjet_wp="Medium", vse_wp="VVLoose", vsmu_wp="VLoose"):
        """add the three DeepTau ID scale factors to the weights container"""
        yr = self.year[:4]
        self._add(f"CMS_eff_t_id_vsjet_{yr}", self.get_vsjet_weights, [vsjet_wp, vse_wp])
        self._add(f"CMS_eff_t_id_vse_{yr}", self.get_vse_weights, [vse_wp])
        self._add(
            f"CMS_eff_t_id_vsmu_{yr}", self.get_vsmu_weights, [vsmu_wp, vse_wp, vsjet_wp]
        )

    def add_trigger_weights(self, region_masks, wp="Medium", dm_binned=True):
        """tau-leg trigger scale factors (`tau_trigger`), one factor per tau_h leg.

        region_masks: {trigtype: per-event bool mask} with trigtype a `tau_trigger` key
        ('ditau', 'ditaujet', 'mutau', 'etau'). The regions must be disjoint; an event in a
        region gets the product of the SFs of all its tau_h legs for that trigtype, every
        other event gets 1. The 2024 SFs describe the PNet HLT paths.

        Nominal: DM-binned SFs (dm_binned=False: inclusive, dm = -1). Uncertainties: one
        nuisance per trigtype and decay mode, CMS_eff_trigger_t_<trigtype>_dm<dm>_<year>,
        moving the legs of that decay mode only (fully correlated between the two legs).
        """
        yr = self.year[:4]
        n_events = len(self.taus_counts)
        counts = ak.to_numpy(self.taus_counts)
        event_index = np.repeat(np.arange(n_events), counts)

        leg_trigtype = np.full(len(event_index), "", dtype=object)
        in_any_region = np.zeros(n_events, dtype=bool)
        for trigtype, mask in region_masks.items():
            mask = np.asarray(mask, dtype=bool)
            if (mask & in_any_region).any():
                raise ValueError(f"tau trigger regions overlap (trigtype {trigtype!r}): the SF product needs disjoint regions")
            in_any_region |= mask
            leg_trigtype[mask[event_index]] = trigtype

        t = self.flat_taus
        pt = ak.to_numpy(t.pt).astype(np.float64)
        dm = ak.to_numpy(t.decayMode).astype(np.int64)
        in_region = leg_trigtype != ""
        if not np.isin(dm[in_region], ALLOWED_DECAY_MODES).all():
            raise ValueError(f"tau trigger SFs exist for decay modes {ALLOWED_DECAY_MODES} only")
        below = in_region & (pt < TAU_TRIGGER_PT_MIN)
        if below.any() and "tau_trigger" not in _clamp_warned:
            _clamp_warned.add("tau_trigger")
            warnings.warn(
                f"{below.sum()} tau_h trigger legs below {TAU_TRIGGER_PT_MIN} GeV (outside the "
                f"tau_trigger fit range): SF evaluated at {TAU_TRIGGER_PT_MIN} GeV",
                RuntimeWarning,
            )
        pt = np.maximum(pt, TAU_TRIGGER_PT_MIN)
        dm_input = dm if dm_binned else np.full(len(dm), -1)

        def leg_sfs(syst):
            sf = np.ones(len(pt))
            for trigtype in region_masks:
                sel = leg_trigtype == trigtype
                if sel.any():
                    sf[sel] = self.cset["tau_trigger"].evaluate(
                        pt[sel], dm_input[sel], trigtype, wp, "sf", syst
                    )
            return sf

        def event_product(sf):
            return ak.to_numpy(ak.prod(ak.unflatten(sf, counts), axis=1))

        nominal_legs = leg_sfs("nom")
        nominal = event_product(nominal_legs)
        self.weights.add(name=f"CMS_eff_trigger_t_{yr}", weight=nominal)
        if self.shift is not None:
            return
        # variations as ratios to the nominal product (weight = 1 for these entries)
        up_legs, down_legs = leg_sfs("up"), leg_sfs("down")
        dm_groups = ALLOWED_DECAY_MODES if dm_binned else ("all",)
        ones = np.ones(n_events)
        for trigtype in region_masks:
            for dm_group in dm_groups:
                sel = leg_trigtype == trigtype
                if dm_group != "all":
                    sel = sel & (dm == dm_group)
                varied = {}
                for direction, legs in (("up", up_legs), ("down", down_legs)):
                    shifted = np.where(sel, legs, nominal_legs)
                    product = event_product(shifted)
                    varied[direction] = np.divide(product, nominal, out=np.ones(n_events), where=nominal != 0)
                self.weights.add(
                    name=f"CMS_eff_trigger_t_{trigtype}_dm{dm_group}_{yr}",
                    weight=ones,
                    weightUp=varied["up"],
                    weightDown=varied["down"],
                )


# -------------------------------------------------------------------------------
# tau energy scale
# -------------------------------------------------------------------------------
def _propagate_to_met(met, pt_before, pt_after, phi):
    """return a copy of `met` with the per-event sum of (pt_after - pt_before) vectors
    subtracted (objects gain pT -> MET loses it)"""
    dpx = ak.sum((pt_after - pt_before) * np.cos(phi), axis=1)
    dpy = ak.sum((pt_after - pt_before) * np.sin(phi), axis=1)
    met_px = met.pt * np.cos(met.phi) - dpx
    met_py = met.pt * np.sin(met.phi) - dpy
    new_met = copy.copy(met)
    new_met["pt"] = ak.values_astype(np.sqrt(met_px**2 + met_py**2), np.float32)
    new_met["phi"] = ak.values_astype(np.arctan2(met_py, met_px), np.float32)
    return new_met


def get_tau_es_factors(events, year, wp_vsjet, wp_vse, syst="nom"):
    """per-tau energy-scale factor (1 for taus outside the correction domain)"""
    cset = load_correction_set(correction_files["tau"][year])
    taus = events.Tau
    flat, counts = ak.flatten(taus), ak.num(taus)
    pt = ak.to_numpy(flat.pt).astype(np.float64)
    eta = np.abs(ak.to_numpy(flat.eta).astype(np.float64))
    dm = ak.to_numpy(flat.decayMode).astype(np.int64)
    gm = ak.to_numpy(flat.genPartFlav).astype(np.int64)
    # genuine tau_h (5), e->tau_h fakes (1, 3) and mu->tau_h fakes (2, 4) are covered;
    # jets (0) are not corrected
    valid = _dm_mask(dm) & (gm >= 1) & (gm <= 5)
    factor = cset["tau_energy_scale"].evaluate(
        pt,
        eta,
        np.where(_dm_mask(dm), dm, 0),
        np.where(valid, gm, 5),
        "DeepTau2018v2p5",
        wp_vsjet,
        wp_vse,
        syst,
    )
    factor = np.where(valid, factor, 1.0)
    return ak.unflatten(ak.Array(factor), counts)


def apply_tau_es_corrections(events, year, shifts, corrections_config):
    """apply the tau energy scale to MC taus (pt and mass), propagate the change to
    the MET collection of every entry of `shifts`, and (optionally) add the
    CMS_scale_t up/down object shifts"""
    met_key = "PuppiMET" if int(get_nano_version(year)) >= 12 else "MET"
    if not hasattr(events, "genWeight"):
        # data: no energy scale, but attach the same bookkeeping fields as for MC so that
        # the histogram expressions (tau_pt_raw, ...) resolve
        taus = copy.copy(events.Tau)
        taus["pt_raw"] = events.Tau.pt
        taus["mass_raw"] = events.Tau.mass
        for collections, _ in shifts:
            collections["Tau"] = taus
        return shifts

    tes_cfg = corrections_config.get("tau_es", {}) or {}
    wp_vsjet = tes_cfg.get("vsjet", "Medium")
    wp_vse = tes_cfg.get("vse", "VVLoose")

    def corrected_taus(factor):
        taus = copy.copy(events.Tau)
        taus["pt_raw"] = events.Tau.pt
        taus["mass_raw"] = events.Tau.mass
        taus["pt"] = ak.values_astype(events.Tau.pt * factor, np.float32)
        taus["mass"] = ak.values_astype(events.Tau.mass * factor, np.float32)
        return taus

    nominal_factor = get_tau_es_factors(events, year, wp_vsjet, wp_vse, "nom")
    nominal_taus = corrected_taus(nominal_factor)
    for collections, _ in shifts:
        collections["Tau"] = nominal_taus
        if met_key in collections:
            collections[met_key] = _propagate_to_met(
                collections[met_key], nominal_taus.pt_raw, nominal_taus.pt, nominal_taus.phi
            )

    if corrections_config["object_shifts"]:
        nominal_collections = shifts[0][0]
        for var in ["up", "down"]:
            taus_var = corrected_taus(get_tau_es_factors(events, year, wp_vsjet, wp_vse, var))
            collections = dict(nominal_collections)
            collections["Tau"] = taus_var
            if met_key in nominal_collections:
                # propagate relative to the nominal-TES MET
                collections[met_key] = _propagate_to_met(
                    nominal_collections[met_key], nominal_taus.pt, taus_var.pt, taus_var.phi
                )
            shifts.append((collections, f"CMS_scale_t_{year[:4]}{var.capitalize()}"))
    return shifts
