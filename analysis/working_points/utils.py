import correctionlib
import awkward as ak
from analysis.filesets.utils import get_nano_version
from analysis.corrections.correctionlib_files import correction_files


def get_ctag_mask(jets: ak.Array, year: str, wp: str):
    """
    Parameters:
    -----------
      jets: Jet collection
      year: {2016preVFP, 2016postVFP, 2017, 2018, 2022preEE, 2022postEE, 2023preBPix, 2023postBPix, 2024}
      wp: {loose, medium, tight}
    """
    nano_version = get_nano_version(year)
    if nano_version == "9":
        tagger = "deepjet"
    elif nano_version == "12":
        tagger = "pnet"
    elif nano_version == "15":
        tagger = "upart"

    wp_map = {"loose": "L", "medium": "M", "tight": "T"}
    tagger_map = {
        "deepjet": "deepJet_wp_values",
        "pnet": "particleNet_wp_values",
        "upart": "UParTAK4_wp_values",
    }

    cset = correctionlib.CorrectionSet.from_file(correction_files["ctagging"][year])
    ctag_wps_evaluator = cset[tagger_map[tagger]]
    cvsb_wp = ctag_wps_evaluator.evaluate(wp_map[wp], "CvB")
    cvsl_wp = ctag_wps_evaluator.evaluate(wp_map[wp], "CvL")

    if tagger == "deepjet":
        pass_ctag_wp = (jets.btagDeepFlavCvB > cvsb_wp) & (
            jets.btagDeepFlavCvL > cvsl_wp
        )
    elif tagger == "pnet":
        pass_ctag_wp = (jets.btagPNetCvB > cvsb_wp) & (jets.btagPNetCvL > cvsl_wp)
    elif tagger == "upart":
        pass_ctag_wp = (jets.btagUParTAK4CvB > cvsb_wp) & (
            jets.btagUParTAK4CvL > cvsl_wp
        )

    return pass_ctag_wp


def get_btag_mask(jets: ak.Array, year: str, wp: str):
    """
    b-tagging working point mask (UParT for nanov15, ParticleNet for nanov12, DeepJet for nanov9)

    Parameters:
    -----------
      jets: Jet collection
      year: {2022preEE, 2022postEE, 2023preBPix, 2023postBPix, 2024}
      wp: {loose, medium, tight, xtight, xxtight}
    """
    nano_version = get_nano_version(year)
    wp_map = {"loose": "L", "medium": "M", "tight": "T", "xtight": "XT", "xxtight": "XXT"}
    if year not in correction_files.get("btagging", {}):
        raise KeyError(f"no b-tagging working points registered for {year}")
    cset = correctionlib.CorrectionSet.from_file(correction_files["btagging"][year])
    if nano_version == "15":
        threshold = cset["UParTAK4_wp_values"].evaluate(wp_map[wp])
        return jets.btagUParTAK4B > threshold
    elif nano_version == "12":
        threshold = cset["particleNet_wp_values"].evaluate(wp_map[wp])
        return jets.btagPNetB > threshold
    threshold = cset["deepJet_wp_values"].evaluate(wp_map[wp])
    return jets.btagDeepFlavB > threshold
