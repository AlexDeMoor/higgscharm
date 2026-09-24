import correctionlib
import numpy as np
import awkward as ak
from analysis.working_points import working_points
from analysis.filesets.utils import get_nano_version
from analysis.corrections.correctionlib_files import correction_files, load_correction_set


def jet_veto(events, year: str, mapname: str = "jetvetomap"):
    """
    From: https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/JME_2022_Summer22EE_jetvetomaps.html
        These are the jet veto maps showing regions with an excess of jets (hot zones) and lack of jets
        (cold zones). Using the phi-symmetry of the CMS detector, these areas with detector and or
        calibration issues can be pinpointed. Non-zero value indicates that the region is vetoed

    From https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps:
        The nominal “loose selection” would be:
            Run2:
                - jet pT > 15 GeV
                - tight jet ID
                - PU jet ID for CHS jets with pT < 50 GeV
                - (jet charged EM fraction + jet neutral EM fraction) < 0.9
                - jets that don’t overlap with PF muon (dR < 0.2)
            Run3:
                - jet pT > 15 GeV
                - tightLepVeto jet ID
                - (jet charged EM fraction + jet neutral EM fraction) < 0.9
    """
    vetomap_names = {
        "2016preVFP": "Summer19UL16_V1",
        "2016postVFP": "Summer19UL16_V1",
        "2017": "Summer19UL17_V1",
        "2018": "Summer19UL18_V1",
        "2022preEE": "Summer22_23Sep2023_RunCD_V1",
        "2022postEE": "Summer22EE_23Sep2023_RunEFG_V1",
        "2023preBPix": "Summer23Prompt23_RunC_V1",
        "2023postBPix": "Summer23BPixPrompt23_RunD_V1",
        "2024": "Summer24Prompt24_RunBCDEFGHI_V1",
    }
    cset = load_correction_set(correction_files["jetvetomaps"][year])

    # nominal loose selection
    nano_version = get_nano_version(year)
    jet_pt = events.Jet.pt > 15
    jets_pu_id = working_points.jet_pileup_id(events, year, "tight")
    jets_id = working_points.jet_id(
        events, year, wp="tight" if nano_version == "9" else "tightlepveto"
    )
    jets_em_fraction = (events.Jet.neEmEF + events.Jet.chEmEF) < 0.9
    events["Jet", "loose_mask"] = jet_pt & jets_pu_id & jets_id & jets_em_fraction

    # get jet veto mask
    jets = events.Jet
    j, n = ak.flatten(jets), ak.num(jets)
    jet_loose_mask = j.loose_mask
    jet_eta_mask = np.abs(j.eta) < 5.19
    jet_phi_mask = np.abs(j.phi) < 3.1415

    in_jet_mask = jet_loose_mask & jet_eta_mask & jet_phi_mask
    in_jets = j.mask[in_jet_mask]

    jets_eta = ak.fill_none(in_jets.eta, 0.0)
    jets_phi = ak.fill_none(in_jets.phi, 0.0)

    vetomap = cset[vetomap_names[year]].evaluate(mapname, jets_eta, jets_phi)
    vetomap = ak.where(in_jet_mask, vetomap, 0.0)
    return ak.any(ak.unflatten(vetomap, n) > 0, axis=1)
