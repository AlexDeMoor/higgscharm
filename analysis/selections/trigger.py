import yaml
import warnings
import numpy as np
import awkward as ak
import importlib.resources
from analysis.filesets.utils import get_dataset_name


# flags registered at run time (register_trigger_flags), on top of trigger_flags.yaml:
# {year (int): {flag: [HLT paths]}}
_REGISTERED_FLAGS = {}


def register_trigger_flags(year, flags):
    """Add (or override) trigger flags for `year` without editing trigger_flags.yaml, e.g.
    from a workflow plugin module (`plugins:` in the workflow yaml).
    flags: {flag: [HLT path names without the 'HLT_' prefix]}"""
    _REGISTERED_FLAGS.setdefault(int(str(year)[:4]), {}).update(
        {flag: list(paths) for flag, paths in flags.items()}
    )


def get_hltpaths_from_flag(flag, year):
    registered = _REGISTERED_FLAGS.get(int(year[:4]), {})
    if flag in registered:
        return registered[flag]
    with importlib.resources.open_text(
        f"analysis.selections", f"trigger_flags.yaml"
    ) as file:
        hlt_paths = yaml.safe_load(file)
    trigger_flags = hlt_paths[int(year[:4])]
    if flag in trigger_flags:
        return trigger_flags[flag]
    else:
        print(f"Trigger flag {flag} is not defined for year {year}")


_missing_hlt_warned = set()

# TrigObj filter bits of the 2024 NanoAODv15 files (from the TrigObj_filterBits branch
# title). Values are lists of bits: an object matches a leg if ANY of the bits is set.
# select_trigger_legs (object_selections.py) produces l1_match_<leg> / l2_match_<leg> for
# every leg listed under the kind of the pair leg (muon, electron or tau).
TRIGOBJ_BITS = {
    "muon": {
        # TODO (mu tau_h): the muon legs of the triggers you use (read the bit meanings from the
        # TrigObj_filterBits branch title of a 2024 NanoAODv15 file)
    },
    "electron": {
        # dedicated hltEle30WPTightGsfTrackIsoFilter bit; bit 1 (hltEle*WPTight*TrackIsoFilter*) is a
        # superset that also fires for the Ele24er cross-trigger leg and Ele32, so it is not used
        "singleele": [18],
        "etau": [17],           # 1e-1tau PNet cross-trigger electron leg (bit 7 = the DeepTau cross path, not used)
    },
    "tau": {
        # TODO (mu tau_h): the tau leg of a mu-tau cross trigger, if you use one
        "etau": [12],           # e-tau cross-trigger tau leg
        "ditau": [11],          # di-tau leg
        "ditaujet": [14],       # di-tau + PFJet tau leg
    },
    "jet": {
        "ditaujet": [17],       # DiTau+Jet (Jet) Signal
    },
}

# offline legs matched to TrigObj of ANOTHER type: kind -> {leg: (TrigObj id, bits, dR)}.
# e tau_h events fire the di-tau paths with the electron reconstructed as an HLT tau (no anti-e at HLT)
TRIGOBJ_CROSS_TYPE = {
    "electron": {
        "ditau": (15, [11], 0.5),
    },
}


def match_trigger_objects(objects, trigobjs, obj_id, filter_bits, dr_max):
    """Per-object bool: at least one TrigObj with |id| == obj_id and any of `filter_bits`
    set lies within dr_max of the object. Works on jagged collections (returns jagged bools)."""
    mask = ak.zeros_like(trigobjs.pt, dtype=bool)
    for bit in filter_bits:
        mask = mask | ((trigobjs.filterBits & (1 << bit)) != 0)
    cands = trigobjs[(np.abs(trigobjs.id) == obj_id) & mask]
    dr = objects.metric_table(cands)
    return ak.any(dr < dr_max, axis=-1)


def trigger_flags_or(events, flags, year):
    """plain OR of a list of trigger flags (MC and data alike; no PD exclusivity)"""
    mask = np.zeros(len(events), dtype="bool")
    for flag in flags:
        mask = mask | trigger_from_flag(events, flag, year)
    return mask


def trigger_region_mask(events, hlt_paths, dataset, year, flags):
    """OR of `flags` for a trigger REGION, keeping the primary-dataset ownership and
    exclusivity encoded in `hlt_paths` (PD -> flags, dict order = PD priority) for data:
      * MC: plain OR of the flags;
      * data: 0 unless this dataset's PD is the owner of the flags; events firing any flag
        of a higher-priority PD are removed (so that every event is taken from one PD only)."""
    dataset_name = get_dataset_name(dataset)
    if hasattr(events, "genWeight"):
        return trigger_flags_or(events, flags, year)
    if dataset_name == "mc":
        # data without a recognised primary-dataset prefix: refuse instead of silently
        # returning an empty (or, in trigger_mask, an inclusive) mask
        raise KeyError(f"data dataset {dataset!r} has no recognised primary dataset name (get_dataset_name)")
    owners = [pd for pd, pd_flags in hlt_paths.items() if any(f in pd_flags for f in flags)]
    if not owners:
        raise KeyError(f"trigger region flags {flags} are not listed under any primary dataset in hlt_paths {hlt_paths}")
    owner = owners[0]
    if dataset_name != owner:
        return np.zeros(len(events), dtype="bool")
    mask = trigger_flags_or(events, flags, year)
    for pd, pd_flags in hlt_paths.items():
        if pd == owner:
            break
        mask = mask & ~trigger_flags_or(events, pd_flags, year)
    return mask


def trigger_from_flag(events, flag, year):
    """OR of the HLT paths behind `flag`. Paths absent from the input file (e.g. a
    path added mid-year, or a menu change) are skipped with a one-time warning
    instead of raising a KeyError, so the mask is always a plain bool array."""
    trigger_mask = np.zeros(len(events), dtype="bool")
    hlt_paths = get_hltpaths_from_flag(flag, year)
    if hlt_paths:
        for hlt_path in hlt_paths:
            if hlt_path not in events.HLT.fields:
                if hlt_path not in _missing_hlt_warned:
                    _missing_hlt_warned.add(hlt_path)
                    warnings.warn(
                        f"HLT path {hlt_path} (flag {flag}, {year}) not in the input "
                        f"file; treated as not fired",
                        RuntimeWarning,
                    )
                continue
            trigger_mask = trigger_mask | ak.to_numpy(events.HLT[hlt_path])
    return trigger_mask


def zzto4l_trigger(events, hlt_paths, dataset, year):
    nevents = len(events)
    dataset_name = get_dataset_name(dataset)
    # compute all trigger masks based on the flags in hlt_paths
    trigger_flags = {}
    for flag in hlt_paths:
        trigger_flags[flag] = trigger_from_flag(events, flag, year)

    if hasattr(events, "genWeight"):
        # compute the combined OR of all flags (for MC datasets)
        mc_trigger_mask = np.zeros(nevents, dtype="bool")
        for flag in trigger_flags:
            mc_trigger_mask = mc_trigger_mask | trigger_flags[flag]
        return mc_trigger_mask

    trigger_flags["TriEle"] = np.zeros(nevents, dtype="bool")

    # ensure each event is taken only from a single PD
    pd_trigger_mask = (
        (np.zeros(nevents, dtype="bool"))
        | (
            ((dataset_name == "EGamma") * np.ones(nevents, dtype=bool))
            & (trigger_flags["DiEle"] | trigger_flags["TriEle"])
        )
        | (
            (
                ((dataset_name == "Muon") or (dataset_name == "DoubleMuon"))
                * np.ones(nevents, dtype=bool)
            )
            & (trigger_flags["DiMu"] | trigger_flags["TriMu"])
            & ~trigger_flags["DiEle"]
            & ~trigger_flags["TriEle"]
        )
        | (
            ((dataset_name == "MuonEG") * np.ones(nevents, dtype=bool))
            & trigger_flags["MuEle"]
            & ~trigger_flags["DiMu"]
            & ~trigger_flags["TriMu"]
            & ~trigger_flags["DiEle"]
            & ~trigger_flags["TriEle"]
        )
        | (
            ((dataset_name == "EGamma") * np.ones(nevents, dtype=bool))
            & trigger_flags["SingleEle"]
            & ~trigger_flags["MuEle"]
            & ~trigger_flags["DiMu"]
            & ~trigger_flags["TriMu"]
            & ~trigger_flags["DiEle"]
            & ~trigger_flags["TriEle"]
        )
        | (
            (
                ((dataset_name == "Muon") or (dataset_name == "SingleMuon"))
                * np.ones(nevents, dtype=bool)
            )
            & trigger_flags["SingleMu"]
            & ~trigger_flags["SingleEle"]
            & ~trigger_flags["MuEle"]
            & ~trigger_flags["DiMu"]
            & ~trigger_flags["TriMu"]
            & ~trigger_flags["DiEle"]
            & ~trigger_flags["TriEle"]
        )
    )
    return pd_trigger_mask


def trigger_mask(events, hlt_paths, dataset, year):
    dataset_name = get_dataset_name(dataset)
    # compute all trigger masks based on the flags in hlt_paths
    trigger_flags = {}
    for pd_flags in hlt_paths.values():
        for flag in pd_flags:
            trigger_flags[flag] = trigger_from_flag(events, flag, year)

    # compute the combined OR of all flags (for MC)
    all_combined_mask = np.zeros(len(events), dtype="bool")
    for flag in trigger_flags:
        all_combined_mask = all_combined_mask | trigger_flags[flag]

    dataset_masks = {}
    for pd_name, flags in hlt_paths.items():
        # OR of all flags for this primary dataset
        mask = np.zeros(len(events), dtype="bool")
        for flag in flags:
            mask = mask | trigger_flags[flag]

        if dataset_masks:
            # exclude the flags of the higher-priority primary datasets
            other_masks = np.zeros(len(events), dtype="bool")
            for other_pd in dataset_masks:
                for flag in hlt_paths[other_pd]:
                    other_masks = other_masks | trigger_flags[flag]
            mask = mask & ~other_masks
        dataset_masks[pd_name] = mask

    if dataset_name == "mc":
        if not hasattr(events, "genWeight"):
            raise KeyError(f"data dataset {dataset!r} has no recognised primary dataset name (get_dataset_name)")
        return all_combined_mask
    if dataset_name not in dataset_masks:
        # a data primary dataset that is not listed in hlt_paths must never fall back to the
        # inclusive OR (it would double count against the listed PDs)
        raise KeyError(
            f"primary dataset {dataset_name!r} (from {dataset!r}) is not listed in hlt_paths "
            f"{list(hlt_paths)}: add it to the workflow or do not run this dataset"
        )
    return dataset_masks[dataset_name]


def trigger_match(leptons: ak.Array, trigobjs: ak.Array, hlt_path: str, year: str):
    """
    Returns DeltaR matched trigger objects

    leptons:
        electrons or muons arrays
    trigobjs:
        trigger objects array
    hlt_path:
        trigger to match
    year:
        dataset year {2016preVFP, 2016postVFP, 2017, 2018, 2022preEE, 2022postEE, 2023preBPix, 2023postBPix}

    https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaNanoAOD#Trigger_bits_how_to
    """
    match_configs = {
        "run2": {
            "IsoMu24": {
                "pt": trigobjs.pt > 22,
                "id": abs(trigobjs.id) == 13,
                "filterbit": (trigobjs.filterBits & 8) > 0,
            },
            "IsoMu27": {
                "pt": trigobjs.pt > 25,
                "filterbit": (trigobjs.filterBits & 8) > 0,
                "id": abs(trigobjs.id) == 13,
            },
            "Mu50": {
                "pt": trigobjs.pt > 45,
                "filterbit": (trigobjs.filterBits & 1024) > 0,
                "id": abs(trigobjs.id) == 13,
            },
            "OldMu100": {
                "pt": trigobjs.pt > 95,
                "filterbit": (trigobjs.filterBits & 2048) > 0,
                "id": abs(trigobjs.id) == 13,
            },
            "TkMu100": {
                "pt": trigobjs.pt > 95,
                "filterbit": (trigobjs.filterBits & 2048) > 0,
                "id": abs(trigobjs.id) == 13,
            },
            "Ele35_WPTight_Gsf": {
                "pt": trigobjs.pt > 33,
                "filterbit": (trigobjs.filterBits & 2) > 0,
                "id": abs(trigobjs.id) == 11,
            },
            "Ele32_WPTight_Gsf": {
                "pt": trigobjs.pt > 30,
                "filterbit": (trigobjs.filterBits & 2) > 0,
                "id": abs(trigobjs.id) == 11,
            },
            "Ele27_WPTight_Gsf": {
                "pt": trigobjs.pt > 25,
                "filterbit": (trigobjs.filterBits & 2) > 0,
                "id": abs(trigobjs.id) == 11,
            },
            "Photon175": {
                "pt": trigobjs.pt > 25,
                "filterbit": (trigobjs.filterBits & 8192) > 0,
                "id": abs(trigobjs.id) == 11,
            },
            "Photon200": {
                "pt": trigobjs.pt > 25,
                "filterbit": (trigobjs.filterBits & 8192) > 0,
                "id": abs(trigobjs.id) == 11,
            },
            "IsoTkMu24": {
                "pt": trigobjs.pt > 22,
                "filterbit": (trigobjs.filterBits & 8) > 0,
                "id": abs(trigobjs.id) == 13,
            },
        },
        "run3": {
            # filterbit: 3 => 1mu
            # id: 13 => mu
            "IsoMu24": {
                "pt": trigobjs.pt > 23,
                "id": abs(trigobjs.id) == 13,
                "filterbit": trigobjs.filterBits & (0x1 << 3) > 0,
            },
            # filterbit: 0 => TrkIsoVVL
            # id: 13 => mu
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8": {
                "pt": trigobjs.pt > 7,
                "id": abs(trigobjs.id) == 13,
                "filterbit": trigobjs.filterBits & (0x1 << 0) > 0,
            },
            # filterbit: 1 => 1e (WPTight)
            # id: 11 => ele
            "Ele30_WPTight_Gsf": {
                "pt": trigobjs.pt > 28,
                "id": abs(trigobjs.id) == 11,
                "filterbit": trigobjs.filterBits & (0x1 << 1) > 0,
            },
        },
    }
    run_key = "run2" if year.startswith("201") else "run3"
    pass_pt = match_configs[run_key][hlt_path]["pt"]
    pass_id = match_configs[run_key][hlt_path]["id"]
    pass_filterbit = match_configs[run_key][hlt_path]["filterbit"]
    trigger_cands = trigobjs[pass_pt & pass_id & pass_filterbit]
    delta_r = leptons.metric_table(trigger_cands)
    pass_delta_r = delta_r < 0.1
    n_of_trigger_matches = ak.sum(pass_delta_r, axis=2)
    trig_matched_locs = n_of_trigger_matches >= 1
    return trig_matched_locs


def trigger_match_mask(events, hlt_paths, year, leptons):
    trigger_match_mask = np.zeros(len(events), dtype="bool")
    for dataset_flags in hlt_paths.values():
        for flag in dataset_flags:
            for hlt_path in get_hltpaths_from_flag(flag, year):
                trig_obj_mask = trigger_match(
                    leptons=leptons,
                    trigobjs=events.TrigObj,
                    hlt_path=hlt_path,
                    year=year,
                )
                trigger_match_mask = trigger_match_mask | trig_obj_mask
    return trigger_match_mask
