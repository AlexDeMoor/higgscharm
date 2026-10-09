import yaml
import json
import numpy as np
import awkward as ak
import importlib.resources
from coffea.lumi_tools import LumiMask
from analysis.selections.trigger import (
    trigger_mask, trigger_match_mask, zzto4l_trigger, trigger_from_flag, trigger_flags_or, trigger_region_mask,
)


def get_lumi_mask(events, year):
    year_map = {
        "2016": "analysis/data/lumimask/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
        "2017": "analysis/data/lumimask/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
        "2018": "analysis/data/lumimask/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
        "2022": "analysis/data/lumimask/Cert_Collisions2022_355100_362760_Golden.txt",
        "2023": "analysis/data/lumimask/Cert_Collisions2023_366442_370790_Golden.txt",
        "2024": "analysis/data/lumimask/Cert_Collisions2024_378981_386951_Golden.txt",
    }
    for key in year_map:
        if year.startswith(key):
            goldenjson = year_map[key]
            break
    else:
        raise ValueError(f"Unrecognized year format: '{year}'")

    if hasattr(events, "genWeight"):
        lumi_mask = np.ones(len(events), dtype="bool")
    else:
        lumi_info = LumiMask(goldenjson)
        lumi_mask = lumi_info(events.run, events.luminosityBlock)
    return lumi_mask == 1


def get_zzto4l_trigger_mask(events, hlt_paths, dataset_key, year):
    return zzto4l_trigger(events, hlt_paths, dataset_key, year)


def get_trigger_mask(events, hlt_paths, dataset_key, year):
    return trigger_mask(events, hlt_paths, dataset_key, year)


def get_trigger_flag_mask(events, flag, year):
    """bool mask of a single trigger flag (OR of its HLT paths), for pass/fail axes"""
    return trigger_from_flag(events, flag, year)


def get_trigger_flags_or(events, flags, year):
    """plain OR of several trigger flags (no primary-dataset exclusivity; for study axes)"""
    return trigger_flags_or(events, flags, year)


def get_trigger_region_mask(events, hlt_paths, dataset, year, flags):
    """OR of the flags of one trigger region; PD-exclusive for data (see trigger_region_mask)"""
    return trigger_region_mask(events, hlt_paths, dataset, year, flags)


def get_trigger_match_mask(events, hlt_paths, year, leptons):
    mask = trigger_match_mask(events, hlt_paths, year, leptons)
    return ak.sum(mask, axis=-1) > 0


def get_metfilters_mask(events, year):
    """AND of the MET filters of analysis/data/met/metfilters.json. A filter (or year)
    missing from the json / the input file raises: skipping it would silently loosen the
    selection."""
    with importlib.resources.path("analysis.data.met", "metfilters.json") as path:
        with open(path, "r") as handle:
            metfilters_years = json.load(handle)
    if year not in metfilters_years:
        raise KeyError(f"no MET filters defined for year {year} in metfilters.json")
    metfilters = metfilters_years[year]
    metfilters_mask = np.ones(len(events), dtype="bool")
    metfilterkey = "mc" if hasattr(events, "genWeight") else "data"
    missing = [mf for mf in metfilters[metfilterkey] if mf not in events.Flag.fields]
    if missing:
        raise KeyError(
            f"MET filters {missing} ({year}, {metfilterkey}) are not in the input file "
            f"{events.metadata.get('filename', '')!r}"
        )
    for mf in metfilters[metfilterkey]:
        metfilters_mask = metfilters_mask & events.Flag[mf]
    return metfilters_mask


def get_stitching_mask(events, dataset, dataset_key, ht_value):
    stitching_mask = np.ones(len(events), dtype="bool")
    if dataset.startswith(dataset_key):
        stitching_mask = events.LHE.HT < ht_value
    return stitching_mask
