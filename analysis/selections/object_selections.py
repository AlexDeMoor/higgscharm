import warnings
import numpy as np
import awkward as ak
from coffea.nanoevents.methods import candidate
from coffea.nanoevents.methods import vector
from coffea.nanoevents.methods.vector import LorentzVector
from analysis.working_points import working_points
from analysis.filesets.utils import get_nano_version
from analysis.selections.htautau import (
    GEN_CHANNEL,
    RANK_ETAU,
    RANK_TAUTAU,
    pt_sorted,
    to_candidates,
    gen_htautau_channel,
    build_htt_pair,
    dzeta,
    collinear_mass,
)
from analysis.selections.trigger import TRIGOBJ_BITS, TRIGOBJ_CROSS_TYPE, match_trigger_objects
from analysis.selections import (
    delta_r_higher,
    delta_r_lower,
    select_dileptons,
    transverse_mass,
    get_closest_lepton,
    assign_lepton_fsr_idx,
    fourlepcand,
    make_cand,
    select_best_zzcandidate,
    select_candidate_mass,
    select_candidate_cjet_dphi,
)


# select_htt_pair_mutau placeholder: warn once per process
_MUTAU_PLACEHOLDER_WARNED = False


class ObjectSelector:

    def __init__(self, object_selection_config, year):
        self.object_selection_config = object_selection_config
        self.year = year

    def select_objects(self, events):
        self.objects = {}
        self.events = events

        for obj_name, obj_config in self.object_selection_config.items():
            # check if object is defined from events or user defined function
            if "events" in obj_config["field"]:
                self.objects[obj_name] = eval(obj_config["field"])
            else:
                selection_function = getattr(self, obj_config["field"])
                selection_function(obj_name)
            if "add_cut" in obj_config:
                for field_to_add in obj_config["add_cut"]:
                    selection_mask = self.get_selection_mask(
                        events=events,
                        obj_name=obj_name,
                        cuts=obj_config["add_cut"][field_to_add],
                    )
                    self.objects[obj_name][field_to_add] = selection_mask
            if "add_field" in obj_config:
                # bring objects and year to local scope for the eval'd expressions
                objects = self.objects
                year = self.year
                for field_name, field_to_add in obj_config["add_field"].items():
                    self.objects[obj_name][field_name] = eval(field_to_add)
            if "cuts" in obj_config:
                selection_mask = self.get_selection_mask(
                    events=events, obj_name=obj_name, cuts=obj_config["cuts"]
                )
                self.objects[obj_name] = self.objects[obj_name][selection_mask]
        return self.objects

    def get_selection_mask(self, events, obj_name, cuts):
        # bring objects and year to local scope
        objects = self.objects
        year = self.year
        # initialize selection mask
        selection_mask = ak.ones_like(self.objects[obj_name].pt, dtype=bool)
        # iterate over all cuts
        for str_mask in cuts:
            mask = eval(str_mask)
            selection_mask = np.logical_and(selection_mask, mask)
        return selection_mask

    # --------------------------------------------------------------------------------
    # ZToLL
    # --------------------------------------------------------------------------------
    def select_dimuons(self, obj_name):
        if "muons" not in self.objects:
            raise ValueError(f"'muons' object has not been defined!")
        self.objects[obj_name] = select_dileptons(self.objects, "muons")

    def select_dielectrons(self, obj_name):
        if "electrons" not in self.objects:
            raise ValueError(f"'electrons' object has not been defined!")
        self.objects[obj_name] = select_dileptons(self.objects, "electrons")

    # --------------------------------------------------------------------------------
    # ZZTo4L (Z+L, ZLL)
    # --------------------------------------------------------------------------------
    def select_zzto4l_leptons(self, obj_name):
        muons = self.objects["muons"]
        muons["lostHits"] = ak.zeros_like(muons.pt)
        electrons = self.objects["electrons"]
        # leptons before FSR recovery/iso correction
        helper_leptons = ak.concatenate([muons, electrons], axis=1)
        helper_leptons["idx"] = ak.local_index(helper_leptons, axis=1)
        helper_leptons = ak.zip(
            {
                "pt": helper_leptons.pt,
                "eta": helper_leptons.eta,
                "phi": helper_leptons.phi,
                "mass": helper_leptons.mass,
                "charge": helper_leptons.charge,
                "pdgId": helper_leptons.pdgId,
                "idx": helper_leptons.idx,
                "is_loose": helper_leptons.is_loose,
                "is_relaxed": helper_leptons.is_relaxed,
                "is_tight": helper_leptons.is_tight,
                "lostHits": helper_leptons.lostHits,
            },
            with_name="PtEtaPhiMCandidate",
            behavior=candidate.behavior,
        )
        # select FSR photons
        fsr_photons = self.objects["fsr_photons"]
        fsr_photons["mass"] = 0
        fsr_photons["charge"] = 0
        fsr_photons["idx"] = ak.local_index(fsr_photons, axis=1)
        # find closest lepton
        closest_lepton, dr_et2 = get_closest_lepton(fsr_photons, helper_leptons)
        fsr_photons["helper_idx"] = ak.local_index(closest_lepton, axis=1)
        fsr_photons["lepton_idx"] = ak.fill_none(closest_lepton.idx, -1)
        # if the FSR photon is matched to multiple leptons, choose the one with lowest dREt2
        fsr_cartesian = ak.cartesian(
            {"fsr1": fsr_photons, "fsr2": fsr_photons}, nested=True, axis=1
        )
        is_duplicate = (
            ak.sum(
                fsr_cartesian.fsr1.lepton_idx == fsr_cartesian.fsr2.lepton_idx, axis=-1
            )
            > 1
        )
        min_dr_et2 = ak.argmin(dr_et2.mask[is_duplicate], axis=-1, keepdims=True)
        is_min = ak.fill_none(
            ak.flatten(min_dr_et2[:, None] == fsr_photons.helper_idx, axis=-1), False
        )
        new_lepton_idx = ak.where(
            ~is_duplicate | (is_duplicate & is_min), fsr_photons.lepton_idx, -1
        )
        fsr_photons["lepton_idx"] = new_lepton_idx

        # add fsr_idx field to leptons
        helper_leptons = assign_lepton_fsr_idx(fsr_photons, helper_leptons)

        # for each FSR photon that was selected, we exclude that photon from the isolation sum of all the leptons in the event
        # This concerns the photons that are in the isolation cone and outside the isolation veto of said leptons dR < 0.4 AND dR > 0.01
        selected_fsr_photons = fsr_photons[fsr_photons.lepton_idx > -1]
        muons_matchedfsr_cartesian = ak.cartesian(
            {"muon": muons, "fsr": selected_fsr_photons},
            nested=True,
            axis=1,
        )
        muons_matchedfsr_cartesian_dr = muons_matchedfsr_cartesian.muon.delta_r(
            muons_matchedfsr_cartesian.fsr
        )
        muons_matchedfsr_cartesian_dr_valid = (muons_matchedfsr_cartesian_dr > 0.01) & (
            muons_matchedfsr_cartesian_dr < 0.4
        )
        valid_matchedfsr_pt = ak.sum(
            muons_matchedfsr_cartesian.fsr.mask[muons_matchedfsr_cartesian_dr_valid].pt,
            axis=-1,
        )
        muon_corrected_iso = muons.pfRelIso03_all - (valid_matchedfsr_pt / muons.pt)
        # update pfRelIso03_all field with corrected isolation
        muons["pfRelIso03_all"] = ak.where(
            muon_corrected_iso > 0, muon_corrected_iso, 0.0
        )
        # update 'is_tight' selection for muons
        muons["is_tight"] = muons.is_tight & (muons.pfRelIso03_all < 0.35)

        electrons["pfRelIso03_all"] = ak.zeros_like(electrons.pt)
        # concatenate muons and electrons with corrected iso
        leptons = ak.concatenate([muons, electrons], axis=1)
        leptons = leptons[ak.argsort(leptons.pt, axis=1, ascending=False)]
        leptons = ak.zip(
            {
                "pt": leptons.pt,
                "eta": leptons.eta,
                "phi": leptons.phi,
                "mass": leptons.mass,
                "charge": leptons.charge,
                "pdgId": leptons.pdgId,
                "pfRelIso03_all": leptons.pfRelIso03_all,
                "idx": helper_leptons.idx,
                "fsr_idx": helper_leptons.fsr_idx,
                "is_loose": leptons.is_loose,
                "is_relaxed": leptons.is_relaxed,
                "is_tight": leptons.is_tight,
                "lostHits": leptons.lostHits,
            },
            with_name="PtEtaPhiMCandidate",
            behavior=candidate.behavior,
        )
        leptons["idx"] = ak.local_index(leptons, axis=1)
        # assign -1 to FSR lepton_idx associated with the excluded leptons
        index_still_present = ak.any(
            fsr_photons.idx == leptons.fsr_idx[:, None], axis=-1
        )
        new_lepton_idx = ak.where(index_still_present, fsr_photons.lepton_idx, -1)
        fsr_photons["lepton_idx"] = new_lepton_idx

        # select leptons with and without matched FSR photons
        has_matched_fsr_photons = leptons.fsr_idx > -1
        leptons_with_matched_fsrphotons = ak.pad_none(
            leptons[has_matched_fsr_photons], 1
        )
        leptons_without_matched_fsrphotons = leptons[~has_matched_fsr_photons]
        # select FSR photons with matched leptons
        has_matched_leptons = fsr_photons.lepton_idx > -1
        fsr_with_matched_leptons = ak.pad_none(fsr_photons[has_matched_leptons], 1)

        # add p4 field to leptons adding the matched FSR photons
        leptons_with_matched_fsrphotons["p4"] = ak.zip(
            {
                "pt": (leptons_with_matched_fsrphotons + fsr_with_matched_leptons).pt,
                "eta": (leptons_with_matched_fsrphotons + fsr_with_matched_leptons).eta,
                "phi": (leptons_with_matched_fsrphotons + fsr_with_matched_leptons).phi,
                "mass": (
                    leptons_with_matched_fsrphotons + fsr_with_matched_leptons
                ).mass,
                "charge": (
                    leptons_with_matched_fsrphotons + fsr_with_matched_leptons
                ).charge,
                "lostHits": leptons_with_matched_fsrphotons.lostHits,
            },
            with_name="PtEtaPhiMCandidate",
            behavior=candidate.behavior,
        )
        leptons_without_matched_fsrphotons["p4"] = ak.zip(
            {
                "pt": leptons_without_matched_fsrphotons.pt,
                "eta": leptons_without_matched_fsrphotons.eta,
                "phi": leptons_without_matched_fsrphotons.phi,
                "mass": leptons_without_matched_fsrphotons.mass,
                "charge": leptons_without_matched_fsrphotons.charge,
                "lostHits": leptons_without_matched_fsrphotons.lostHits,
            },
            with_name="PtEtaPhiMCandidate",
            behavior=candidate.behavior,
        )
        # concatenate leptons with and without matched FSR photons after FSR recovery
        leptons = ak.where(
            # concatenate only events with matched FSR photons
            # to avoid including None values from 'leptons_with_matched_fsrphotons'
            ak.sum(leptons.fsr_idx > -1, axis=1) == 0,
            leptons_without_matched_fsrphotons,
            ak.concatenate(
                [leptons_with_matched_fsrphotons, leptons_without_matched_fsrphotons],
                axis=1,
            ),
        )
        self.objects[obj_name] = leptons

    def select_zcandidates(self, obj_name):
        """selects Z candidates for SR and all CRS"""
        # get Z candidates
        zcand = ak.combinations(self.objects["leptons"], 2, fields=["l1", "l2"])
        # check that Z candidates leptons pass the loose id
        zcand = zcand[zcand.l1.is_loose & zcand.l2.is_loose]
        # add Z candidate, p4, pT and idx fields
        zcand["p4"] = zcand.l1.p4 + zcand.l2.p4
        zcand["pt"] = zcand.p4.pt
        zcand["idx"] = ak.local_index(zcand, axis=1)
        # add the Z candidates to objects
        self.objects[obj_name] = zcand

    def select_best_zcandidate(self, obj_name):
        """selects best Z candidate as the one closest to the nominal Z mass"""
        zmass = 91.1876
        best_zcand_idx = ak.argmin(
            np.abs(self.objects["zcandidates"].p4.mass - zmass), axis=1
        )
        best_zcand = self.objects["zcandidates"][
            best_zcand_idx == self.objects["zcandidates"].idx
        ]
        self.objects[obj_name] = best_zcand

    def select_other_relaxed_leptons(self, obj_name):
        """
        selects additional relaxed leptons in the Z+L CR. Adds the mask 'pass_selection' to additional relaxed leptons that pass the analysis selection
        """
        # select best Z candidates
        best_zcands = self.objects["best_zcandidates"]
        # select relaxed leptons (whose idx are different to best Z candidate lepton's idx)
        relaxed_leptons = self.objects["leptons"][self.objects["leptons"].is_relaxed]
        relaxed_leptons_bestzl1_idx = relaxed_leptons.idx != best_zcands.l1.idx[:, None]
        relaxed_leptons_bestzl2_idx = relaxed_leptons.idx != best_zcands.l2.idx[:, None]
        relaxed_leptons_bestz_idx_mask = ak.flatten(
            relaxed_leptons_bestzl1_idx & relaxed_leptons_bestzl2_idx, axis=-1
        )
        relaxed_leptons = relaxed_leptons[relaxed_leptons_bestz_idx_mask]

        # add the 'pass_selection' attribute to 'relaxed_leptons'. It flags relaxed leptons passing the analysis selection
        is_tight = relaxed_leptons.is_tight == ak.ones_like(
            best_zcands.l1.idx[:, None], dtype=bool
        )
        is_tight = ak.flatten(is_tight, axis=-1)
        # ghost removal: ∆R(η, φ) > 0.02 between each of the leptons (to protect against split tracks)
        relaxed_leptons_bestzl1_dr = relaxed_leptons.metric_table(best_zcands.l1)
        relaxed_leptons_bestzl2_dr = relaxed_leptons.metric_table(best_zcands.l2)
        relaxed_leptons_bestz_dr_mask = ak.flatten(
            (relaxed_leptons_bestzl1_dr > 0.02) & (relaxed_leptons_bestzl2_dr > 0.02),
            axis=-1,
        )
        # QCD suppression cut: invariant mass of relaxed lepton and the opposite sign tight lepton from the best Z candidate should satisfy m2l > 4 GeV
        relaxed_leptons_bestzl1_opposite_charge = ak.flatten(
            relaxed_leptons.charge != best_zcands.l1.charge[:, None], axis=-1
        )
        relaxed_leptons_bestzl2_opposite_charge = ak.flatten(
            relaxed_leptons.charge != best_zcands.l2.charge[:, None], axis=-1
        )
        relaxed_leptons_bestzl1_cartesian = ak.cartesian(
            {"lepton": relaxed_leptons.p4, "zl1": best_zcands.l1.p4},
            nested=True,
            axis=1,
        )
        relaxed_leptons_bestzl1_mass = ak.flatten(
            (
                relaxed_leptons_bestzl1_cartesian.lepton
                + relaxed_leptons_bestzl1_cartesian.zl1
            ).mass,
            axis=-1,
        )
        relaxed_leptons_bestzl1_mass_mask = relaxed_leptons_bestzl1_mass > 4
        relaxed_leptons_bestzl2_cartesian = ak.cartesian(
            {"lepton": relaxed_leptons.p4, "zl2": best_zcands.l2.p4},
            nested=True,
            axis=1,
        )
        relaxed_leptons_bestzl2_mass = ak.flatten(
            (
                relaxed_leptons_bestzl2_cartesian.lepton
                + relaxed_leptons_bestzl2_cartesian.zl2
            ).mass,
            axis=-1,
        )
        relaxed_leptons_bestzl2_mass_mask = relaxed_leptons_bestzl2_mass > 4
        qcd_suppression_mask = (
            relaxed_leptons_bestzl1_opposite_charge & relaxed_leptons_bestzl1_mass_mask
        ) | (
            relaxed_leptons_bestzl2_opposite_charge & relaxed_leptons_bestzl2_mass_mask
        )

        # get full pass selection mask
        pass_selection = is_tight & relaxed_leptons_bestz_dr_mask & qcd_suppression_mask
        relaxed_leptons["pass_selection"] = pass_selection

        # add relaxed leptons to objects
        self.objects[obj_name] = relaxed_leptons

    def select_trilepton(self, obj_name):
        self.objects[obj_name] = (
            self.objects["best_zcandidates"].l1
            + self.objects["best_zcandidates"].l2
            + ak.firsts(self.objects["other_relaxed_leptons"])
        )

    def select_zzcandidates(self, obj_name):
        """selects ZZ candidates for SR and CRs"""
        self.objects[obj_name] = make_cand(
            self.objects["zcandidates"], kind="zz", sort_by_mass=True
        )

    def select_zllcandidates_os(self, obj_name):
        """selects Zll candidates for CRs"""
        self.objects[obj_name] = make_cand(
            self.objects["zcandidates"], kind="zll", sort_by_mass=False, os_method=True
        )

    def select_zllcandidates_ss(self, obj_name):
        """selects Zll candidates for CRs"""
        self.objects[obj_name] = make_cand(
            self.objects["zcandidates"], kind="zll", sort_by_mass=False, os_method=False
        )

    def select_best_zzcandidate(self, obj_name):
        """selects best ZZ candidates for SR"""
        self.objects[obj_name] = select_best_zzcandidate(self.objects["zzcandidates"])

    def select_best_1fcr_zllcandidate(self, obj_name):
        """selects best Zll candidates for 3P1F CR"""
        self.objects[obj_name] = select_best_zzcandidate(
            self.objects["zllcandidates"], "is_1fcr"
        )

    def select_best_2fcr_zllcandidate(self, obj_name):
        """selects best Zll candidates for 2P2F CR"""
        self.objects[obj_name] = select_best_zzcandidate(
            self.objects["zllcandidates"], "is_2fcr"
        )

    def select_best_sscr_zllcandidate(self, obj_name):
        """selects best Zll candidates for SS CR"""
        self.objects[obj_name] = select_best_zzcandidate(
            self.objects["zllcandidates"], "is_sscr"
        )

    def select_mass_4e_1fcr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zllcandidate_1fcr"], "4e"
        )

    def select_mass_2e2mu_1fcr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zllcandidate_1fcr"], "2e2mu"
        )

    def select_mass_2mu2e_1fcr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zllcandidate_1fcr"], "2mu2e"
        )

    def select_mass_4mu_1fcr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zllcandidate_1fcr"], "4mu"
        )

    def select_mass_4e_2fcr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zllcandidate_2fcr"], "4e"
        )

    def select_mass_2e2mu_2fcr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zllcandidate_2fcr"], "2e2mu"
        )

    def select_mass_2mu2e_2fcr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zllcandidate_2fcr"], "2mu2e"
        )

    def select_mass_4mu_2fcr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zllcandidate_2fcr"], "4mu"
        )

    def select_mass_4e_sscr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zllcandidate_sscr"], "4e"
        )

    def select_mass_2e2mu_sscr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zllcandidate_sscr"], "2e2mu"
        )

    def select_mass_2mu2e_sscr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zllcandidate_sscr"], "2mu2e"
        )

    def select_mass_4mu_sscr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zllcandidate_sscr"], "4mu"
        )

    def select_mass_4mu_sr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zzcandidate"], "4mu"
        )

    def select_mass_4e_sr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zzcandidate"], "4e"
        )

    def select_mass_2e2mu_sr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zzcandidate"], "2e2mu"
        )

    def select_mass_2mu2e_sr(self, obj_name):
        self.objects[obj_name] = select_candidate_mass(
            self.objects["best_zzcandidate"], "2mu2e"
        )

    def select_dphi_4e(self, obj_name):
        self.objects[obj_name] = select_candidate_cjet_dphi(
            self.objects["best_zzcandidate"], "4e", self.objects["cjets"]
        )

    def select_dphi_4mu(self, obj_name):
        self.objects[obj_name] = select_candidate_cjet_dphi(
            self.objects["best_zzcandidate"], "4mu", self.objects["cjets"]
        )

    def select_dphi_2mu2e(self, obj_name):
        self.objects[obj_name] = select_candidate_cjet_dphi(
            self.objects["best_zzcandidate"], "2mu2e", self.objects["cjets"]
        )

    def select_dphi_2e2mu(self, obj_name):
        self.objects[obj_name] = select_candidate_cjet_dphi(
            self.objects["best_zzcandidate"], "2e2mu", self.objects["cjets"]
        )

    def select_dphi_inclusive(self, obj_name):
        self.objects[obj_name] = select_candidate_cjet_dphi(
            self.objects["best_zzcandidate"], "inclusive", self.objects["cjets"]
        )

    # --------------------------------------------------------------------------------
    # Jets (shared by the H+c workflows)
    # --------------------------------------------------------------------------------
    def _cvl_field(self):
        return "btagUParTAK4CvL" if get_nano_version(self.year) == "15" else "btagPNetCvL"

    def select_candidate_cjet(self, obj_name):
        """the c-jet candidate is the c-tagged jet with the highest CvL score
        (UParT for NanoAODv15, ParticleNet otherwise)"""
        cjets = self.objects["cjets"]
        score = cjets[self._cvl_field()]
        self.objects["candidate_cjet"] = cjets[
            ak.argmax(score, axis=1) == ak.local_index(cjets, axis=1)
        ]
        self.objects["candidate_cjet_lorentzvector"] = to_candidates(
            self.objects["candidate_cjet"]
        )

    JET_FIELDS = (
        "hadronFlavour",
        "partonFlavour",
        "btagUParTAK4CvL",
        "btagUParTAK4CvB",
        "btagUParTAK4B",
        "btagPNetCvL",
        "btagPNetCvB",
        "btagPNetB",
        "btagDeepFlavCvL",
        "btagDeepFlavCvB",
        "btagDeepFlavB",
    )

    def select_leading_jets(self, obj_name):
        """jet1/jet2/jet3 = the three highest-pT selected jets (padded with None) as
        per-event candidates carrying the flavour and tagger fields"""
        jets = pt_sorted(self.objects["jets"])
        padded = ak.pad_none(jets, target=3)
        for i in range(3):
            self.objects[f"jet{i + 1}"] = to_candidates(padded[:, i], self.JET_FIELDS)

    # --------------------------------------------------------------------------------
    # H -> tautau: generator level
    # --------------------------------------------------------------------------------
    def select_gen_htautau(self, obj_name):
        """gen-level decay channel + generator objects of the H -> tautau system.
        Data (no GenPart) gets channel 0 and empty collections."""
        ev = self.events
        if "GenPart" not in ev.fields:
            # data: channel 0, empty collections, and per-event placeholders carrying the SAME
            # fields as the MC branch (all None) so that every axis expression still resolves
            nan = ak.values_astype(ak.zeros_like(ev.event), np.float32) * np.nan
            zeros = ak.values_astype(ak.zeros_like(ev.event), np.int32)
            empty = ev.Jet[:, :0]
            self.objects["gen_channel"] = zeros
            for name in ("gen_leptons", "gen_muons", "gen_electrons", "gen_taus",
                         "gen_vistaus", "gen_higgs", "gen_hfquarks", "gen_jets"):
                self.objects[name] = empty

            def placeholder(fields):
                return ak.zip(
                    {f: ak.mask(nan, ak.zeros_like(ev.event, dtype=bool)) for f in ("pt", "eta", "phi", "mass") + tuple(fields)},
                    with_name="PtEtaPhiMCandidate", behavior=candidate.behavior,
                )
            for i in range(3):
                self.objects[f"genjet{i + 1}"] = placeholder(("hadronFlavour", "partonFlavour"))
            self.objects["gen_vistau1"] = placeholder(("status", "charge"))
            self.objects["gen_vistau2"] = placeholder(("status", "charge"))
            self.objects["gen_lepton"] = placeholder(("pdgId",))
            self.objects["gen_met"] = ak.zip({"pt": nan, "phi": nan})
            self.objects["gen_higgs_mass_vis"] = nan
            self.objects["gen_higgs_pt_vis"] = nan
            self.objects["gen_higgs_dr_legs"] = nan
            return

        gp = ev.GenPart
        channel, leptons_from_tau, hard_taus = gen_htautau_channel(ev)
        self.objects["gen_channel"] = channel
        self.objects["gen_leptons"] = pt_sorted(leptons_from_tau)
        self.objects["gen_muons"] = pt_sorted(leptons_from_tau[np.abs(leptons_from_tau.pdgId) == 13])
        self.objects["gen_electrons"] = pt_sorted(leptons_from_tau[np.abs(leptons_from_tau.pdgId) == 11])
        self.objects["gen_taus"] = hard_taus
        # visible hadronic taus of the hard-process taus only (GenVisTau also lists tau_h from
        # B/D decays); the mother of a GenVisTau is the tau
        vistaus = ev.GenVisTau
        vis_from_hard = ak.fill_none(vistaus.parent.hasFlags(["isHardProcess"]), False)
        self.objects["gen_vistaus"] = pt_sorted(vistaus[vis_from_hard])
        self.objects["gen_higgs"] = gp[(gp.pdgId == 25) & gp.hasFlags(["isLastCopy"])]
        # OUTGOING hard-process heavy-flavour quark(s) of the associated production
        # (status 21 = incoming beam partons, which also carry isHardProcess)
        self.objects["gen_hfquarks"] = pt_sorted(
            gp[((np.abs(gp.pdgId) == 4) | (np.abs(gp.pdgId) == 5)) & gp.hasFlags(["isHardProcess"]) & (gp.status != 21)]
        )
        # gen jets cleaned against the visible tau decay products (otherwise the leading gen
        # jets are simply the tau_h / the lepton)
        genjets = ev.GenJet
        genjets = genjets[
            delta_r_higher(genjets, self.objects["gen_vistaus"], 0.4)
            & delta_r_higher(genjets, self.objects["gen_leptons"], 0.4)
        ]
        genjets = pt_sorted(genjets)
        self.objects["gen_jets"] = genjets
        padded = ak.pad_none(genjets, target=3)
        for i in range(3):
            self.objects[f"genjet{i + 1}"] = to_candidates(padded[:, i], ("hadronFlavour", "partonFlavour"))
        vis = ak.pad_none(self.objects["gen_vistaus"], target=2)
        self.objects["gen_vistau1"] = to_candidates(vis[:, 0], ("status", "charge"))
        self.objects["gen_vistau2"] = to_candidates(vis[:, 1], ("status", "charge"))
        self.objects["gen_lepton"] = to_candidates(
            ak.pad_none(self.objects["gen_leptons"], target=1)[:, 0], ("pdgId",)
        )
        self.objects["gen_met"] = ev.GenMET
        # gen-level visible di-tau system: the two visible legs of the hard-process taus
        # (lepton + tau_h, or the two tau_h), in the full phase space
        def to_vectors(coll):
            # plain Lorentz vectors (GenPart carries no charge -> no PtEtaPhiMCandidate)
            return ak.zip(
                {"pt": coll.pt, "eta": coll.eta, "phi": coll.phi, "mass": coll.mass},
                with_name="PtEtaPhiMLorentzVector", behavior=vector.behavior,
            )
        gen_legs = pt_sorted(ak.concatenate(
            [to_vectors(self.objects["gen_leptons"]), to_vectors(self.objects["gen_vistaus"])], axis=1
        ))
        gen_legs = ak.pad_none(gen_legs, target=2)
        gen_vis = gen_legs[:, 0] + gen_legs[:, 1]
        self.objects["gen_higgs_mass_vis"] = gen_vis.mass
        self.objects["gen_higgs_pt_vis"] = gen_vis.pt
        self.objects["gen_higgs_dr_legs"] = gen_legs[:, 0].delta_r(gen_legs[:, 1])

    # --------------------------------------------------------------------------------
    # H -> tautau: reconstructed di-tau pair
    # --------------------------------------------------------------------------------
    def _store_pair(self, pairs):
        # keep the best-ranked candidate only (jagged with 0/1 entries so that the
        # yaml `cuts` (opposite sign, ...) still apply)
        self.objects["htt_pair"] = pairs[:, :1]

    def select_htt_pair_mutau(self, obj_name):
        """TODO (mu tau_h): build the (muon, tau_h) pair candidates and keep the best-ranked one,
        as select_htt_pair_etau does below (build_htt_pair + a ranking recipe in htautau.py).

        Placeholder: stores an EMPTY pair collection (with the right l1 = muon / l2 = tau_h record
        structure), so the workflow runs but every pair-based selection is empty."""
        global _MUTAU_PLACEHOLDER_WARNED
        if not _MUTAU_PLACEHOLDER_WARNED:
            _MUTAU_PLACEHOLDER_WARNED = True
            warnings.warn(
                "select_htt_pair_mutau is a placeholder: no mu tau_h pair is built "
                "(analysis/selections/object_selections.py)",
                RuntimeWarning,
            )
        pairs = ak.cartesian({"l1": self.objects["muons"], "l2": self.objects["taus"]}, axis=1)
        self._store_pair(pairs[:, :0])

    def select_htt_pair_etau(self, obj_name):
        self._store_pair(build_htt_pair(self.objects["electrons"], self.objects["taus"], False, RANK_ETAU))

    def select_htt_pair_tautau(self, obj_name):
        # combinations preserve the input order, so pT-sorting the taus first
        # guarantees l1 = the leading-pT tau_h of every candidate pair
        taus = pt_sorted(self.objects["taus"])
        self._store_pair(build_htt_pair(None, taus, True, RANK_TAUTAU))

    def select_htt_legs(self, obj_name):
        """per-event legs of the selected pair (option records) + the tau_h legs as a
        jagged collection (for the tau ID scale factors) + all pair-forming leptons
        (for the jet cleaning)"""
        pair = self.objects["htt_pair"]
        best = ak.firsts(pair)
        self.objects["l1"] = best.l1
        self.objects["l2"] = best.l2
        is_tau = lambda leg: "decayMode" in leg.fields
        legs = []
        if is_tau(pair.l1):
            legs.append(pair.l1)
        if is_tau(pair.l2):
            legs.append(pair.l2)
        self.objects["htt_taus"] = ak.concatenate(legs, axis=1)
        self.objects["htt_leptons"] = ak.concatenate(
            [to_candidates(pair.l1), to_candidates(pair.l2)], axis=1
        )
        # tighter tau_h working points on the hadronic leg(s) (sensitivity studies)
        for leg_name in ("l1", "l2"):
            leg = self.objects[leg_name]
            is_tau_leg = "decayMode" in leg.fields
            for wp_name, field, idx in (("vse", "idDeepTau2018v2p5VSe", 6), ("vsmu", "idDeepTau2018v2p5VSmu", 4), ("vsjet", "idDeepTau2018v2p5VSjet", 6)):
                self.objects[f"{leg_name}_{wp_name}_tight"] = (
                    ak.fill_none(leg[field] >= idx, False) if is_tau_leg
                    else ak.zeros_like(self.events.event, dtype=bool)
                )
        # pT-sorted leading / sub-leading objects (per-event option records) for the
        # trigger turn-on axes -- the collections are not pT-ordered after the energy
        # corrections, so ak.firsts(objects[...]) would not be the leading object
        for name in ("muons", "electrons", "taus", "muons_study", "electrons_study", "taus_study"):
            if name in self.objects:
                padded = ak.pad_none(pt_sorted(self.objects[name]), target=2)
                self.objects[f"{name}_lead"] = padded[:, 0]
                self.objects[f"{name}_sublead"] = padded[:, 1]

    def select_trigger_legs(self, obj_name):
        """Trigger-object matching of the selected pair legs (2024 NanoAODv15 filter bits, see
        TRIGOBJ_BITS in analysis/selections/trigger.py). Per-event booleans:
          l1_match_<leg>, l2_match_<leg> for leg in the channel-relevant legs (plus the
          TRIGOBJ_CROSS_TYPE legs, e.g. l1_match_ditau = the electron fired an HLT di-tau leg), and
          jet_match_ditaujet = an offline jet (pT > 65, |eta| < 4.7, dR > 0.5 from both
          legs) matched to the di-tau+jet HLT jet leg."""
        best = self.objects["htt_pair"][:, :1]   # best-ranked pair only (jagged, 0 or 1 entry)
        trigobjs = self.events.TrigObj
        for leg_name in ("l1", "l2"):
            leg = best[leg_name]
            if "decayMode" in leg.fields:
                kind, obj_id, dr_max = "tau", 15, 0.5
            elif "deltaEtaSC" in leg.fields:
                kind, obj_id, dr_max = "electron", 11, 0.1
            else:
                kind, obj_id, dr_max = "muon", 13, 0.1
            for spec_name, bits in TRIGOBJ_BITS.get(kind, {}).items():
                matched = match_trigger_objects(leg, trigobjs, obj_id, bits, dr_max)
                self.objects[f"{leg_name}_match_{spec_name}"] = ak.fill_none(ak.firsts(matched), False)
            for spec_name, (cross_id, bits, cross_dr) in TRIGOBJ_CROSS_TYPE.get(kind, {}).items():
                matched = match_trigger_objects(leg, trigobjs, cross_id, bits, cross_dr)
                self.objects[f"{leg_name}_match_{spec_name}"] = ak.fill_none(ak.firsts(matched), False)
        # offline jet matched to the di-tau+jet jet leg (id 1, bit 17), away from both tau legs
        jets = self.events.Jet
        jets = jets[(jets.pt > 65) & (np.abs(jets.eta) < 4.7)]
        jets = jets[delta_r_higher(jets, self.objects["htt_leptons"], 0.5)]
        jet_matched = match_trigger_objects(jets, trigobjs, 1, TRIGOBJ_BITS["jet"]["ditaujet"], 0.4)
        self.objects["jet_match_ditaujet"] = ak.fill_none(ak.any(jet_matched, axis=1), False)

    def select_htt_met_vars(self, obj_name):
        """di-tau + MET observables (per-event option floats)"""
        l1, l2, met = self.objects["l1"], self.objects["l2"], self.objects["met"]
        vis = l1 + l2
        self.objects["higgs_mass_vis"] = vis.mass
        self.objects["higgs_pt"] = vis.pt
        self.objects["higgs_eta"] = vis.eta
        self.objects["higgs_phi"] = vis.phi
        self.objects["higgs_dr_legs"] = l1.delta_r(l2)
        self.objects["higgs_dphi_legs"] = l1.delta_phi(l2)
        self.objects["higgs_deta_legs"] = l1.eta - l2.eta
        self.objects["higgs_charge_product"] = l1.charge * l2.charge
        self.objects["mt_lepton_met"] = transverse_mass(l1, met)
        self.objects["mt_tau_met"] = transverse_mass(l2, met)
        self.objects["higgs_mt_legs"] = np.sqrt(2.0 * l1.pt * l2.pt * (1.0 - np.cos(l1.delta_phi(l2))))
        self.objects["higgs_mt_tot"] = np.sqrt(
            self.objects["higgs_mt_legs"] ** 2 + self.objects["mt_lepton_met"] ** 2 + self.objects["mt_tau_met"] ** 2
        )
        self.objects["dphi_lepton_met"] = l1.delta_phi(met)
        self.objects["dphi_tau_met"] = l2.delta_phi(met)
        self.objects["dphi_higgs_met"] = vis.delta_phi(met)
        self.objects["higgs_dzeta"] = dzeta(l1, l2, met)
        self.objects["higgs_mass_col"], self.objects["higgs_mass_col_valid"] = collinear_mass(l1, l2, met)
        self.objects["higgs_met_pt"] = (vis + met).pt

    def select_acceptance_flags(self, obj_name):
        """gen-level kinematic acceptance of each channel (trigger-study denominators
        evaluated on generator objects); always plain booleans"""
        el = self.objects["gen_electrons"]
        vt = self.objects["gen_vistaus"]
        el_in = ((el.pt > 25) & (np.abs(el.eta) < 2.1)) | ((el.pt > 31) & (np.abs(el.eta) < 2.5))
        tau30 = (vt.pt > 30) & (np.abs(vt.eta) < 2.3)
        tau30_25 = (vt.pt > 30) & (np.abs(vt.eta) < 2.5)
        fill = lambda mask: ak.fill_none(ak.values_astype(mask, bool), False)
        # TODO (mu tau_h): in_gen_acc_mutau, the gen-level image of your mu tau_h object selection
        self.objects["in_gen_acc_etau"] = fill(ak.any(el_in, axis=1) & ak.any(tau30, axis=1))
        self.objects["in_gen_acc_tautau"] = fill(ak.sum(tau30_25, axis=1) >= 2)
