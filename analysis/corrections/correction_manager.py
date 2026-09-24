import fnmatch
import numpy as np
from coffea.analysis_tools import Weights
from analysis.filesets.utils import get_nano_version
from analysis.corrections.muon import MuonWeights
from analysis.corrections.ctag import CTagCorrector
from analysis.corrections.ctag2d import CTag2DCorrector
from analysis.corrections.jetvetomaps import jet_veto
from analysis.corrections.pileup import add_pileup_weight
from analysis.corrections.nnlops import add_nnlops_weight
from analysis.corrections.lhepdf import add_lhepdf_weight
from analysis.corrections.electron import ElectronWeights
from analysis.corrections.jerc import apply_jerc_corrections
from analysis.corrections.lhescale import add_scalevar_weight
from analysis.corrections.met import apply_met_phi_corrections
from analysis.corrections.muon_ss import apply_muon_ss_corrections
from analysis.corrections.partonshower import add_partonshower_weight
from analysis.corrections.toppt import add_toppt_weight
from analysis.corrections.higgs_hf import add_higgs_hf_weight
from analysis.corrections.electron_ss import apply_electron_ss_corrections
from analysis.corrections.tau import TauWeights, apply_tau_es_corrections
from analysis.utils.helpers import update


def object_corrector_manager(events, year, corrections_config):
    """apply object level corrections"""
    objcorr_config = corrections_config["object"] or []
    met_field_key = "PuppiMET" if int(get_nano_version(year)) >= 12 else "MET"

    vetoed_events = events
    # nominal entry with the uncorrected jets / MET: used only when no jet correction is
    # configured (apply_jerc_corrections inserts its own corrected nominal entry, so it
    # must receive an EMPTY list -- passing a pre-filled one processed every chunk twice)
    shifts = [({"Jet": events.Jet, met_field_key: events[met_field_key]}, None)]
    if objcorr_config:
        # apply jet energy scale and smearing corrections
        if "jet" in objcorr_config:
            shifts = apply_jerc_corrections(
                events=events,
                year=year,
                shifts=[],
                corrections_config=corrections_config,
            )

        if "muon" in objcorr_config:
            # apply muon scale and smearing corrections
            shifts = apply_muon_ss_corrections(
                events=events,
                year=year,
                shifts=shifts,
                corrections_config=corrections_config,
            )
        else:
            for i in range(len(shifts)):
                shifts[i][0]["Muon"] = events.Muon

        if "electron" in objcorr_config:
            # electron scale and smearing corrections are applied only to run3 datasets
            # they are already applied in nanov9 datasets
            if year.startswith("202"):
                shifts = apply_electron_ss_corrections(
                    events=events,
                    year=year,
                    shifts=shifts,
                    corrections_config=corrections_config,
                )
        else:
            for i in range(len(shifts)):
                shifts[i][0]["Electron"] = events.Electron

        if "tau" in objcorr_config:
            # tau energy scale (MC only); propagated to the MET of every shift
            shifts = apply_tau_es_corrections(
                events=events,
                year=year,
                shifts=shifts,
                corrections_config=corrections_config,
            )

        # apply jet veto (on the corrected jets of the nominal shift; the veto is geometric,
        # so it is applied identically to every shift)
        if "jetveto" in objcorr_config:
            nominal_events = update(events, {"Jet": shifts[0][0]["Jet"]})
            event_veto = jet_veto(nominal_events, year)
            vetoed_events = events[~event_veto]
            for collections, _ in shifts:
                for key in collections:
                    collections[key] = collections[key][~event_veto]

    return vetoed_events, shifts


def _trigger_regions(trigger_config, category, selection_masks):
    """region-aware trigger SF config {regions: {label: selection name}, categories: [fnmatch
    patterns]} -> {label: bool mask over the pruned events}, or None when `category` is not
    listed (the SF then does not apply there)"""
    patterns = trigger_config.get("categories", ["*"])
    if not any(fnmatch.fnmatch(str(category), pat) for pat in patterns):
        return None
    regions = trigger_config["regions"]
    missing = [name for name in regions.values() if name not in (selection_masks or {})]
    if missing:
        raise KeyError(f"trigger SF regions {missing} are not defined event selections")
    return {label: selection_masks[name] for label, name in regions.items()}


def weight_manager(pruned_ev, year, dataset, workflow_config, category, shift, selection_masks=None):
    """apply event level corrections (weights)

    selection_masks: {selection name: bool mask over `pruned_ev`}, needed by weights that
    depend on an event-selection region (tau-leg trigger SFs)."""
    # get weights config info
    weights_config = workflow_config.corrections_config["event_weights"]
    # initialize weights container
    weights_container = Weights(len(pruned_ev), storeIndividual=True)
    # get nano version
    nano_version = get_nano_version(year)
    # add weights
    if hasattr(pruned_ev, "genWeight"):
        if weights_config["genWeight"]:
            weights_container.add("genweight", pruned_ev.genWeight)

        if "pileupWeight" in weights_config:
            if weights_config["pileupWeight"]:
                add_pileup_weight(
                    events=pruned_ev,
                    year=year,
                    shift=shift,
                    weights_container=weights_container,
                )
        if "partonshowerWeight" in weights_config:
            if weights_config["partonshowerWeight"]:
                if "PSWeight" in pruned_ev.fields:
                    add_partonshower_weight(
                        events=pruned_ev,
                        weights_container=weights_container,
                        shift=shift,
                    )
        if "lhepdfWeight" in weights_config:
            if weights_config["lhepdfWeight"]:
                if "LHEPdfWeight" in pruned_ev.fields:
                    add_lhepdf_weight(
                        events=pruned_ev,
                        weights_container=weights_container,
                        shift=shift,
                    )
        if "lhescaleWeight" in weights_config:
            if weights_config["lhescaleWeight"]:
                if "LHEScaleWeight" in pruned_ev.fields:
                    add_scalevar_weight(
                        events=pruned_ev,
                        weights_container=weights_container,
                        shift=shift,
                    )
        if "nnlopsWeight" in weights_config:
            if weights_config["nnlopsWeight"]:
                if dataset.startswith("GluGluH"):
                    add_nnlops_weight(
                        events=pruned_ev,
                        weights_container=weights_container,
                    )
        if "toppTWeight" in weights_config:
            if weights_config["toppTWeight"]:
                # tt only (AN-23-102 line 566); the module no-ops on other datasets
                add_toppt_weight(
                    events=pruned_ev,
                    weights_container=weights_container,
                    shift=shift,
                    dataset=dataset,
                )
        if "higgsHFWeight" in weights_config:
            if weights_config["higgsHFWeight"]:
                # ggH/VBF only; per-event, keyed on gen-jet heavy flavour, so it
                # replaces the mis-scoped flat lnN on the merged higgsbkg group
                add_higgs_hf_weight(
                    events=pruned_ev,
                    weights_container=weights_container,
                    shift=shift,
                    dataset=dataset,
                    # charm: AN-23-102 scopes this to the ggH+HF composition and our
                    # signal is H+c, so the c-flavour variant is the relevant one.
                    flav="c",
                )
        # object-level scale factors (lepton / tau ID) only in the categories listed in
        # `event_weights.object_sf_categories` (fnmatch patterns, default: all). Validation
        # regions such as gen-level categories should not be scaled by SFs of whatever
        # reconstructed objects happen to be around.
        sf_patterns = weights_config.get("object_sf_categories") or ["*"]
        if isinstance(sf_patterns, str):
            sf_patterns = [sf_patterns]
        apply_object_sf = any(fnmatch.fnmatch(str(category), pat) for pat in sf_patterns)

        if "muon" in weights_config and apply_object_sf:
            if weights_config["muon"]:
                if "selected_muons" in pruned_ev.fields:
                    muon_weights = MuonWeights(
                        events=pruned_ev,
                        year=year,
                        shift=shift,
                        weights=weights_container,
                    )
                    if "id" in weights_config["muon"]:
                        if weights_config["muon"]["id"]:
                            muon_weights.add_id_weights(
                                id_wp=weights_config["muon"]["id"]
                            )
                    if "iso" in weights_config["muon"]:
                        if weights_config["muon"]["iso"]:
                            muon_weights.add_iso_weights(
                                id_wp=weights_config["muon"]["id"],
                                iso_wp=weights_config["muon"]["iso"],
                            )
                    trigger_config = weights_config["muon"].get("trigger")
                    if isinstance(trigger_config, dict):
                        # region-aware single-muon trigger SF ({regions: {<path>: <selection>}})
                        regions = _trigger_regions(trigger_config, category, selection_masks)
                        for path, mask in (regions or {}).items():
                            muon_weights.add_region_trigger_weights(
                                mask,
                                id_wp=weights_config["muon"]["id"],
                                iso_wp=weights_config["muon"]["iso"],
                                path=path,
                            )
                    elif trigger_config:
                        muon_weights.add_trigger_weights(
                            dataset=dataset,
                            id_wp=weights_config["muon"]["id"],
                            iso_wp=weights_config["muon"]["iso"],
                            hlt_paths=workflow_config.event_selection["hlt_paths"],
                        )
        if "electron" in weights_config and apply_object_sf:
            if weights_config["electron"]:
                if "selected_electrons" in pruned_ev.fields:
                    electron_weights = ElectronWeights(
                        events=pruned_ev,
                        year=year,
                        shift=shift,
                        weights=weights_container,
                    )
                    if "id" in weights_config["electron"]:
                        if weights_config["electron"]["id"]:
                            electron_weights.add_id_weights(
                                id_wp=weights_config["electron"]["id"]
                            )
                    if "reco" in weights_config["electron"]:
                        if weights_config["electron"]["reco"]:
                            if nano_version == "9":
                                electron_weights.add_reco_weights("RecoAbove20")
                                electron_weights.add_reco_weights("RecoBelow20")
                            elif nano_version == "12":
                                electron_weights.add_reco_weights("RecoBelow20")
                                electron_weights.add_reco_weights("Reco20to75")
                                electron_weights.add_reco_weights("RecoAbove75")
                            elif nano_version == "15":
                                electron_weights.add_reco_weights("Reco20to75")
                                electron_weights.add_reco_weights("RecoAbove75")
                    trigger_config = weights_config["electron"].get("trigger")
                    if isinstance(trigger_config, dict):
                        # region-aware single-electron trigger SF ({regions: {Ele30: <selection>}})
                        regions = _trigger_regions(trigger_config, category, selection_masks)
                        for path, mask in (regions or {}).items():
                            electron_weights.add_region_trigger_weights(
                                mask, id_wp=weights_config["electron"]["id"], path=path
                            )
                    elif trigger_config:
                        # legacy path (all selected electrons); 2024 only via the region-aware form
                        if nano_version != "15":
                            electron_weights.add_hlt_weights(
                                id_wp=weights_config["electron"]["id"],
                            )

        if "tau" in weights_config and apply_object_sf:
            if weights_config["tau"]:
                # DeepTau ID SFs on the tau_h legs of the selected di-tau pair
                if "selected_htt_taus" in pruned_ev.fields:
                    tau_weights = TauWeights(
                        events=pruned_ev,
                        weights=weights_container,
                        year=year,
                        shift=shift,
                    )
                    tau_weights.add_id_weights(
                        vsjet_wp=weights_config["tau"].get("vsjet", "Medium"),
                        vse_wp=weights_config["tau"].get("vse", "VVLoose"),
                        vsmu_wp=weights_config["tau"].get("vsmu", "VLoose"),
                    )
                    # tau-leg trigger SFs: {regions: {trigtype: selection name}, categories: [patterns]},
                    # applied only in the categories that require those trigger regions
                    trigger_config = weights_config["tau"].get("trigger")
                    if trigger_config:
                        regions = _trigger_regions(trigger_config, category, selection_masks)
                        if regions is not None:
                            tau_weights.add_trigger_weights(regions)
        if "ctagging" in weights_config:
            if weights_config["ctagging"]:
                working_points = ["loose", "medium", "tight"]
                if weights_config["ctagging"]["wp"] not in working_points:
                    raise ValueError(
                        f"There are no available c-tag SFs for the working point. Please specify {working_points}"
                    )
                ctag_corrector = CTagCorrector(
                    events=pruned_ev,
                    weights=weights_container,
                    worging_point=weights_config["ctagging"]["wp"],
                    category=category,
                    year=year,
                    shift=shift,
                )
                ctag_corrector.add_ctag_weights(flavor="b")
                ctag_corrector.add_ctag_weights(flavor="c")
                ctag_corrector.add_ctag_weights(flavor="light")

        if "ctagging_2d" in weights_config:
            if weights_config["ctagging_2d"]:
                ctag2d_corrector = CTag2DCorrector(
                    events=pruned_ev,
                    weights=weights_container,
                    year=year,
                    shift=shift,
                )
                ctag2d_corrector.add_weights()
    else:
        weights_container.add("weight", np.ones(len(pruned_ev)))
    return weights_container
