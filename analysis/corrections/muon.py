import json
import correctionlib
import numpy as np
import awkward as ak
from typing import Type
from coffea.analysis_tools import Weights
from analysis.filesets.utils import get_nano_version
from analysis.corrections.utils import unflat_sf
from analysis.corrections.correctionlib_files import correction_files, load_correction_set


class MuonWeights:
    """
    Muon weights class

    Parameters:
    -----------
        events:
            pruned events
        weights:
            Weights container
        year:
            Year of the dataset {2016preVFP, 2016postVFP, 2017, 2018, 2022preEE, 2022postEE, 2023preBPix, 2023postBPix, 2024}
        syst:
            syst syst
        id_wp:
            ID working point {loose, medium, tight}
        iso_wp:
            Iso working point {loose, medium, tight}

    more info at: https://cms-analysis-corrections.docs.cern.ch/corrections/MUO/
    """

    def __init__(
        self,
        events: ak.Array,
        weights: Type[Weights],
        year: str,
        shift: str,
    ) -> None:
        self.events = events
        self.muons = events.selected_muons
        self.weights = weights
        self.year = year
        self.shift = shift
        self.nano_version = get_nano_version(year)

        self.flat_muons = ak.flatten(self.muons)
        self.muons_counts = ak.num(self.muons)

        # get muon id/iso/HLT correction set
        self.cset = load_correction_set(
            correction_files["muon"][year]
        )
        # get muon HLT eff correction set
        self.cset_eff = load_correction_set(
            correction_files["muon_hlt"][year]
        )

    def add_id_weights(self, id_wp):
        """
        add muon ID weights to weights container
        """
        nominal_weights = self.get_id_weights(id_wp, syst="nominal")
        if self.shift is None:
            # get 'up' and 'down' weights
            up_weights = self.get_id_weights(id_wp, syst="systup")
            down_weights = self.get_id_weights(id_wp, syst="systdown")
            # add scale factors to weights container
            self.weights.add(
                name=f"CMS_eff_m_id_{self.year[:4]}",
                weight=nominal_weights,
                weightUp=up_weights,
                weightDown=down_weights,
            )
        else:
            self.weights.add(
                name=f"CMS_eff_m_id_{self.year[:4]}",
                weight=nominal_weights,
            )

    def add_iso_weights(self, id_wp, iso_wp):
        """
        add muon iso weights to weights container
        """
        # get nominal scale factors
        nominal_weights = self.get_iso_weights(id_wp, iso_wp, syst="nominal")
        if self.shift is None:
            # get 'up' and 'down' weights
            up_weights = self.get_iso_weights(id_wp, iso_wp, syst="systup")
            down_weights = self.get_iso_weights(id_wp, iso_wp, syst="systdown")
            # add nominal, up and down weights to weights container
            self.weights.add(
                name=f"CMS_eff_m_iso_{self.year[:4]}",
                weight=nominal_weights,
                weightUp=up_weights,
                weightDown=down_weights,
            )
        else:
            # add nominal weights to weights container
            self.weights.add(
                name=f"CMS_eff_m_iso_{self.year[:4]}",
                weight=nominal_weights,
            )

    def add_trigger_weights(self, id_wp, iso_wp, hlt_paths, dataset):
        """
        add muon iso weights to weights container
        """
        # get nominal scale factors
        nominal_weights = self.get_hlt_weights(
            id_wp=id_wp,
            iso_wp=iso_wp,
            syst="nominal",
        )
        if self.shift is None:
            """
            # get 'up' and 'down' weights
            up_weights = self.get_hlt_weights(
                id_wp=id_wp,
                iso_wp=iso_wp,
                syst="systup",
            )
            down_weights = self.get_hlt_weights(
                id_wp=id_wp,
                iso_wp=iso_wp,
                syst="systdown",
            )
            """
            # add nominal, up and down weights to weights container
            self.weights.add(
                name=f"CMS_eff_m_trigger_{self.year[:4]}",
                weight=nominal_weights,
                # weightUp=up_weights,
                # weightDown=down_weights,
            )
        else:
            # add nominal weights to weights container
            self.weights.add(
                name=f"CMS_eff_m_trigger_{self.year[:4]}",
                weight=nominal_weights,
            )

    def add_region_trigger_weights(self, region_mask, id_wp, iso_wp, path):
        """Single-muon trigger scale factor for the events of one trigger region.

        Called by weight_manager for `event_weights.muon.trigger: {regions: {<path>: <selection>},
        categories: [...]}` (region_mask = that event selection; see analysis/workflows/README.md).

        TODO (mu tau_h): implement it, following ElectronWeights.add_region_trigger_weights
        (analysis/corrections/electron.py): pick the MUO correction for your trigger path and
        (id_wp, iso_wp) denominator from self.cset, evaluate it on the leading selected muon of
        the events in region_mask (SF = 1 elsewhere, mind the pT / eta range of the map) and add
        it to self.weights with its up / down variations."""
        raise NotImplementedError(
            "MuonWeights.add_region_trigger_weights is not implemented yet "
            "(analysis/corrections/muon.py): set `trigger: false` in event_weights.muon "
            "until it is"
        )

    def get_id_weights(self, id_wp, syst):
        """Compute muon ID weights"""
        id_corrections = {
            "loose": "NUM_LooseID_DEN_TrackerMuons",
            "medium": "NUM_MediumID_DEN_TrackerMuons",
            "tight": "NUM_TightID_DEN_TrackerMuons",
        }
        # get muons that pass the id wp, and within SF binning
        pt_lower_limit = 10 if self.nano_version == "15" else 15
        muon_pt_mask = (self.flat_muons.pt > pt_lower_limit) & (
            self.flat_muons.pt < 500.0
        )
        muon_eta_mask = np.abs(self.flat_muons.eta) < 2.399
        in_muon_mask = muon_pt_mask & muon_eta_mask
        in_muons = self.flat_muons.mask[in_muon_mask]

        # get muons pT and eta (replace None values with some 'in-limit' value).
        # Run 2 (nanov9) maps are binned in |eta|; the Run 3 maps (2022-2024) are
        # binned in signed eta [-2.4, 2.4], so only take |eta| for nanov9.
        muon_pt = ak.fill_none(in_muons.pt, 15.0)
        muon_eta = ak.fill_none(in_muons.eta, 0.0)
        if self.nano_version == "9":
            muon_eta = np.abs(muon_eta)

        sfs = self.cset[id_corrections[id_wp]].evaluate(muon_eta, muon_pt, syst)
        weights = unflat_sf(
            sfs,
            in_muon_mask,
            self.muons_counts,
        )
        return weights

    def get_iso_weights(self, id_wp, iso_wp, syst):
        """Compute muon iso weights"""
        iso_corrections = {
            "9": {
                "loose": {
                    "loose": "NUM_LooseRelIso_DEN_LooseID",
                    "medium": None,
                    "tight": None,
                },
                "medium": {
                    "loose": "NUM_LooseRelIso_DEN_MediumID",
                    "medium": None,
                    "tight": "NUM_TightRelIso_DEN_MediumID",
                },
                "tight": {
                    "loose": "NUM_LooseRelIso_DEN_TightIDandIPCut",
                    "medium": None,
                    "tight": "NUM_TightRelIso_DEN_TightIDandIPCut",
                },
            },
            "12": {
                "loose": {
                    "loose": "NUM_LoosePFIso_DEN_LooseID",
                    "medium": "NUM_LoosePFIso_DEN_MediumID",
                    "tight": "NUM_LoosePFIso_DEN_TightID",
                },
                "medium": {
                    "loose": None,
                    "medium": None,
                    "tight": None,
                },
                "tight": {
                    "loose": None,
                    "medium": "NUM_TightPFIso_DEN_MediumID",
                    "tight": "NUM_TightPFIso_DEN_TightID",
                },
            },
            "15": {
                "loose": {
                    "loose": "NUM_LoosePFIso_DEN_LooseID",
                    "medium": "NUM_LoosePFIso_DEN_MediumID",
                    "tight": "NUM_LoosePFIso_DEN_TightID",
                },
                "medium": {
                    "loose": None,
                    "medium": None,
                    "tight": None,
                },
                "tight": {
                    "loose": None,
                    "medium": "NUM_TightPFIso_DEN_MediumID",
                    "tight": "NUM_TightPFIso_DEN_TightID",
                },
            },
        }
        correction_name = iso_corrections[self.nano_version][iso_wp][id_wp]
        if correction_name is None:
            raise ValueError(
                f"There are no muon ISO weights for id wp '{id_wp}' and iso wp '{iso_wp}' combination"
            )
        # get 'in-limits' muons
        pt_lower_limit = 10 if self.nano_version == "15" else 15
        muon_pt_mask = (self.flat_muons.pt > pt_lower_limit) & (
            self.flat_muons.pt < 500.0
        )
        muon_eta_mask = np.abs(self.flat_muons.eta) < 2.399
        in_muon_mask = muon_pt_mask & muon_eta_mask
        in_muons = self.flat_muons.mask[in_muon_mask]

        # get muons pT and abseta (replace None values with some 'in-limit' value)
        muon_pt = ak.fill_none(in_muons.pt, 15)
        muon_eta = ak.fill_none(in_muons.eta, 0.0)
        if self.nano_version == "9":
            muon_eta = np.abs(muon_eta)

        weights = unflat_sf(
            self.cset[correction_name].evaluate(
                muon_eta,
                muon_pt,
                syst,
            ),
            in_muon_mask,
            self.muons_counts,
        )
        return weights

    def get_hlt_weights(self, id_wp, iso_wp, syst):
        """Compute muon HLT weights"""

        hlt_path_id_map = {
            "2016preVFP": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2016postVFP": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2017": "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight",
            "2018": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2022preEE": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2022postEE": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2023preBPix": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2023postBPix": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2024": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
        }

        if not ((id_wp == "tight") & (iso_wp == "tight")):
            raise ValueError(
                f"There are no muon HLT weights for id wp '{self.id_wp}' and iso wp '{self.iso_wp}' combination"
            )

        # get muons within SF binning
        pt_upper_limit = 29.0 if self.year == "2017" else 26.0
        muon_pt_mask = (self.flat_muons.pt > pt_upper_limit) & (
            self.flat_muons.pt < 199.99
        )
        muon_eta_mask = np.abs(self.flat_muons.eta) < 2.399

        in_muons_mask = muon_pt_mask & muon_eta_mask
        in_muons = self.flat_muons.mask[in_muons_mask]

        # get muons pT and abseta (replace None values with some 'in-limit' value)
        muon_pt = ak.fill_none(in_muons.pt, 30.0)
        muon_eta = ak.fill_none(in_muons.eta, 0)
        muon_abseta = ak.fill_none(np.abs(in_muons.eta), 0)

        # for single muon events, compute SF from POG SF
        # for double muon events, compute SF from data/mc muon hlt efficiencies
        kind = "single" if ak.all(ak.num(self.muons) == 1) else "double"

        if kind == "single":
            sf = self.cset[hlt_path_id_map[self.year]].evaluate(
                muon_abseta, muon_pt, syst
            )
            nominal_sf = unflat_sf(sf, in_muons_mask, self.muons_counts)

        elif kind == "double":
            data_eff = self.cset_eff["Muon-HLT-DataEff"].evaluate(
                syst,
                hlt_path_id_map[self.year],
                muon_abseta if self.nano_version == "9" else muon_eta,
                muon_pt,
            )
            data_eff = ak.where(in_muons_mask, data_eff, ak.ones_like(data_eff))
            data_eff = ak.unflatten(data_eff, self.muons_counts)
            data_eff_leading = ak.firsts(data_eff)
            data_eff_subleading = ak.pad_none(data_eff, target=2)[:, 1]
            full_data_eff = (
                data_eff_leading
                + data_eff_subleading
                - data_eff_leading * data_eff_subleading
            )
            full_data_eff = ak.fill_none(full_data_eff, 1)

            mc_eff = self.cset_eff["Muon-HLT-McEff"].evaluate(
                syst,
                hlt_path_id_map[self.year],
                muon_abseta if self.nano_version == "9" else muon_eta,
                muon_pt,
            )
            mc_eff = ak.where(in_muons_mask, mc_eff, ak.ones_like(mc_eff))
            mc_eff = ak.unflatten(mc_eff, self.muons_counts)
            mc_eff_leading = ak.firsts(mc_eff)
            mc_eff_subleading = ak.pad_none(mc_eff, target=2)[:, 1]
            full_mc_eff = (
                mc_eff_leading + mc_eff_subleading - mc_eff_leading * mc_eff_subleading
            )
            full_mc_eff = ak.fill_none(full_mc_eff, 1)

            nominal_sf = full_data_eff / full_mc_eff

        return nominal_sf
