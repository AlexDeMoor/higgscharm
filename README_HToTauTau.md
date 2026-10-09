# cH, H→ττ + charm — analysis guide (2024 / Summer24, NanoAODv15)

This file describes the H→ττ layer built on the framework: samples, objects, the di-tau pair, triggers,
scale factors, categories, histograms and plots of the three `htautau_*` workflows. It ends with
**[the μτh channel: your starting point](#the-μτh-channel-your-starting-point)**. See
[`README.md`](README.md) for setup and how to run things, and
[`analysis/workflows/README.md`](analysis/workflows/README.md) for the YAML reference.

The three workflows share `analysis/workflows/includes/htautau_common.yaml` (`extends:`):

| workflow | channel | pair `(l1, l2)` | status |
|---|---|---|---|
| `htautau_etau` | e τh | (electron, τh) | complete |
| `htautau_tautau` | τh τh | (leading τh, sub-leading τh) | complete |
| `htautau_mutau` | μ τh | (muon, τh) | starting point, to be built |

---

## Samples

| dataset (fileset key) | process | files | events | xsec [pb] | note |
|---|---|---|---|---|---|
| `HplusCharm_HtoTauTau` (`hc_tautau`) | `H+c` | 26 | 10.3 M | 0.006082 | 4FS cH, σ = 97.0 fb × BR(H→ττ) 0.0627 |
| `HplusBottom_HtoTauTau` (`hb_tautau`) | `H+b` | 24 | 10.0 M | 0.045194 | 5FS bH, σ = 720.8 fb × 0.0627 |

Private FullSim productions under `/eos/user/a/ademoor/HcTT/`, H→ττ inclusive (all decay channels, no
generator filter). Both are **signals** (workflow group `datasets.signal`): the analysis measures the cH and
bH cross sections and their ratio, so bH is not part of the background stack or of "Total Background"
in the postprocess cutflows. The generator weights are ±constant (aMC@NLO, 21 % / 30 % negative in cH / bH):
normalise with Σw (the job metadata), never with event counts, and expect small categories to
have non-monotonic cutflows. The pilot file `bH_fullsim_Summer24_v1_merged/tranche0_test/` (5011
events) is kept out of the fileset for smoke tests.

**Backgrounds and data (2026-10-08)**, all `RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2`, in
`analysis/filesets/2024_nanov15.yaml` (the fileset JSON holds 41 entries):

| process | fileset key(s) | samples | xsec [pb] |
|---|---|---|---|
| `DYll` | `dy_ll` | DYto2E / DYto2Mu-2Jets, MLL-50 and MLL-10to50 (amcatnloFXFX) | MLL-50: 2124.08 per flavour (NNLO QCD + NLO EW, FEWZ); MLL-10to50: 7076.7 / 7020.0 (XSDB, NLO) |
| `DYtautau` | `dy_tautau` | DYto2Tau-2Jets, MLL-50 and MLL-10to50 | 2124.08; 7100.0 |
| `Wjets` | `wjets_nlo` | WtoENu / WtoMuNu / WtoTauNu-2Jets (inclusive NLO; no inclusive WtoLNu exists in this campaign) | 21132.0 / 21132.0 / 21449.0 (NNLO QCD × NNLO EW, MATRIX, PDG 2025) |
| `tt`, `Single Top`, `Diboson` | `tt`, `singletop`, `diboson` | as before | as before |
| `ggH_htt`, `qqH_htt`, `VH_htt`, `ttH_nonbb` | `ggh_tautau`, `qqh_tautau`, `vh_tautau`, `tth_nonbb` | GluGluHto2Tau, VBFHto2Tau, W±H / ZH Hto2Tau, TTH HtoNon2B (no ttH→ττ sample) | LHCHWG 13.6 TeV × BR(ττ) 0.0627 (ttH: × (1 − BR(bb)), to confirm) |

Z→ττ (irreducible) and Z→ℓℓ (ℓ→τh fakes) are separate processes; the MLL-50 entries also keep the
key `dy_nlo` for the Z workflows. Data: `Tau` (eras C–I + I v2), `Muon0/1`, `EGamma0/1`. The channel
yamls list their PDs under **`data_blinded:`**, which no runner / postprocess step reads: the data are
defined but **not processed** (rename to `data:` to unblind). `datasets.mc` (common include) makes a
default `runner.py` submission include every background above.

---

## Objects (`includes/htautau_common.yaml`)

| object | selection |
|---|---|
| `electrons` | mvaIso WP90, supercluster gap veto (1.44 < \|η_SC\| < 1.57), \|dxy\| < 0.045, \|dz\| < 0.2, two-leg kinematics (a looser (pT, η) branch for the cross-trigger leg, a harder one for the single-electron leg; eτh overrides the thresholds) |
| `muons` | **defined by the μτh channel**; the common file only holds a loose placeholder (= the veto muon), so the other channels just count them |
| `taus` | DeepTau2018v2p5 **VSjet Medium, VSe VVLoose, VSmu VLoose**; DM ∈ {0, 1, 10, 11}; \|dz\| < 0.2; \|q\| = 1; ΔR > 0.5 from the veto e/μ; generic baseline pT > 20, \|η\| < 2.3 (eτh: pT > 24.6; τhτh: pT > 24.6, \|η\| < 2.1, see Triggers) |
| `veto_muons`, `veto_electrons` | pT > 10; μ: loose ID + PFIso04 < 0.3; e: (mvaNoIso WP90 + PFIso03 < 0.3) **or** mvaIso WP90, so the veto is a superset of the eτh signal electron (since 2026-10-08; with mvaNoIso alone ~4 % of the signal electrons were not veto electrons and `n_veto_e == 1` lost them); **no impact-parameter cuts** (a second lepton must veto even if displaced or non-isolated) |
| `jets`, `cjets`, `bjets` | pT > 20, \|η\| < 2.5, tight-lepton-veto ID (recomputed from the energy fractions: `Jet_jetId` does not exist in v15), ΔR > 0.5 from the pair legs; c / b tagging = UParT medium working points |
| `met` | `PuppiMET`, with the JEC/JER (Type-1 re-propagated relative to the NanoAOD jet pT) and tau-energy-scale corrections propagated (the correction manager writes them back into `PuppiMET`). No recoil or xy correction. NanoAODv15 also has `PuppiMET_cov*` (for SVfit / FastMTT), `PFMET_*` (monitored as `pfmet_*`) and `DeepMET*` |
| `*_study` | validation objects: the object-level denominators of the trigger studies (`muons_study`: channel-defined, placeholder in the common file) |

Working points live in `analysis/working_points/working_points.py` (`muon_id`, `muon_iso`,
`electron_id`, `tau_id_vsjet/vse/vsmu`, `tau_decay_mode`, `jet_id`, `jet_ctagging`, `jet_btagging`).

## The di-tau pair and the derived objects

The pair chain is a sequence of `ObjectSelector` methods (`analysis/selections/object_selections.py`,
helpers in `analysis/selections/htautau.py`), in this order in the common include:

| object (`field:`) | produces (in `objects`) |
|---|---|
| `htt_pair` (`select_htt_pair_<channel>`) | `htt_pair`: all (l1, l2) candidates with ΔR > 0.5, ranked, **best one kept** (jagged, 0 or 1 entry per event; the charge requirement is an event selection, `os_pair`) |
| `legs` (`select_htt_legs`) | `l1`, `l2` (per-event option records, `None` when there is no pair), `htt_taus` (the τh legs, used by the τ scale factors), `htt_leptons` (both legs, used for the jet cleaning), `l{1,2}_vs{e,mu,jet}_tight`, and `<collection>_lead` / `_sublead` for muons, electrons, taus and the `*_study` collections |
| `trigger_legs` (`select_trigger_legs`) | `l{1,2}_match_<leg>` for every trigger leg of `TRIGOBJ_BITS[<kind of the leg>]`, `jet_match_ditaujet` |
| `leading_jets` (`select_leading_jets`) | `jet1`, `jet2`, `jet3` |
| `candidate_cjet` (`select_candidate_cjet`) | the c-jet candidate |
| `htt_met_vars` (`select_htt_met_vars`) | `higgs_mass_vis`, `higgs_pt/eta/phi`, `higgs_dr_legs`, `higgs_dphi_legs`, `higgs_deta_legs`, `higgs_charge_product`, `mt_lepton_met`, `mt_tau_met`, `higgs_mt_legs`, `higgs_mt_tot`, `dphi_*_met`, `higgs_dzeta`, `higgs_mass_col` (+ `higgs_mass_col_valid`), `higgs_met_pt` |
| `gen` (`select_gen_htautau`) | `gen_channel`, `gen_leptons`, `gen_muons`, `gen_electrons`, `gen_vistaus`, `gen_higgs`, `gen_hfquarks`, `gen_jets`, `genjet1..3`, `gen_vistau1/2`, `gen_lepton`, `gen_met`, `gen_higgs_mass_vis`, ... |
| `acc` (`select_acceptance_flags`) | `in_gen_acc_<channel>`: the gen-level image of the **reco kinematic cuts** of the channel (eτh: e (pT > 25, \|η\| < 2.1) or (pT > 32, \|η\| < 2.5), τ_vis pT > 24.6, \|η\| < 2.3; τhτh: two τ_vis with pT > 24.6, \|η\| < 2.1) |
| `best_cvl_jet` (`select_best_cvl_jet`) | the selected jet with the highest CvL (all `jets`, not only the c-tagged ones) |
| `cjet_truth` (`select_cjet_truth`) | ΔR(c-jet candidate, nearest hard-process gen c), the pT and CvL rank of the leading true c-jet (hadronFlavour 4) among the selected jets (MC) |
| `tau_response` (`select_tau_response`) | `l{1,2}_tau_response(_raw)`: reco τh pT (after / before TES) / gen visible pT, genuine τh only (genPartFlav 5; `Tau_genPartIdx` indexes `GenVisTau`, resolved by hand because coffea cross-references it to `GenPart`) |

Pair ranking (`build_htt_pair` + a recipe, applied with successive stable sorts): eτh (`RANK_ETAU`)
takes the most isolated electron, then the highest electron pT, then the highest τh raw VSjet score,
then the highest τh pT. τhτh (`RANK_TAUTAU`) ranks by raw VSjet score and pT of both legs, with
l1 = the leading-pT τh.

**Gen classification** (`gen_htautau_channel`): 1 μτh, 2 eτh, 3 τhτh, 4 fully leptonic, 0 none (and
data). Leptons count only if their direct mother is a τ **and** they carry
`isDirectHardProcessTauDecayProduct`: τ leptons from B/D-hadron decays, frequent in the cH/bH
samples, must not change the classification. `gen_vistaus` are the `GenVisTau` whose mother (the
last-copy τ) carries `fromHardProcess` and `isLastCopy` (pT > 10 builder cut). Until 2026-10-08 the
mother was required to carry `isHardProcess`, which only the first copy has: about 2/3 of the
genuine visible τ were dropped (bH pilot: 2137 → 5960 of 5972 kept; events with two visible τ
628 → 1806). Gen jets are cleaned against the visible τ decay products, and `gen_hfquarks` are the
outgoing hard-process c/b quarks.

**Trigger-study and acceptance numbers produced before 2026-10-08** (`trigger_eff.md`, `in_gen_acc`,
the `trig_gen*` layouts, gen m_vis / gen-jet histograms) used the old `gen_vistaus` and the old
acceptance cuts: they are stale until the samples are rerun.

## Categories

The complete channels follow one naming scheme (one output / plot folder each):

| category | meaning |
|---|---|
| `gen_<ch>` | truth ττ decay = the channel, no reco requirement |
| `den_<ch>_study` | validation: the object-level trigger-study denominator |
| `sel_<ch>` | full offline selection **without** trigger |
| `sel_<ch>_trig` | + `trigger` (the analysis baseline) |
| `sel_<ch>_trig_gen<ch>` | + the truth channel (signal purity) |

Selections defined in the common include: `goodvertex`, `lumimask`, `met_filters`, `one_pair`,
`os_pair`, `ss_pair`, `gen_mutau` / `gen_etau` / `gen_tautau`, plus the optional
`mt_lepton_met_lt70`, `bveto` and `atleast_one_cjet`. Scale factors are applied only in the
categories matching `event_weights.object_sf_categories` (`sel_*`). Cutflows: `cutflow_<cat>.csv`
holds the Σ genWeight (no scale factors) and `cutflow_raw_<cat>.csv` the event counts. Since
2026-10-08 the job metadata cutflow starts with `initial` (before the jet veto, = Σw used for the
normalisation) and `jet_veto` (after it); the per-category CSVs list the category cuts only.

`event_selection.mc_only_categories` lists the categories that are meaningless or double counted on
data (`gen_<ch>`, `den_<ch>_study`, `sel_<ch>` without trigger, `sel_<ch>_trig_gen<ch>`): they are
not filled (no histograms, no cutflow) when running on data.

## Triggers

**Flags.** A flag is a name for one or more HLT paths, defined per year in
`analysis/selections/trigger_flags.yaml` (one flag per path for H→ττ). A plugin can also register
flags at run time with `register_trigger_flags` (see the workflows README).

**Primary datasets.** `event_selection.hlt_paths` maps a primary dataset to its flags. On data an
event is taken from the first PD (in dict order) whose flags fired; on MC it is the OR of the flags.

**Gated trigger regions.** In the complete channels `trigger` is **not** the flat OR of the paths. It
is the OR of disjoint *trigger regions*, and each region requires:
1. the HLT bit(s): `get_trigger_region_mask(events, hlt_paths, dataset, year, [flags])`;
2. offline cuts on the selected pair legs, on the efficiency plateau and inside the pT/η range of the
   trigger scale factor;
3. the pair legs matched to the trigger objects of that path: `objects['l1_match_<leg>']`, with
   ΔR < 0.1 for e/μ and < 0.5 for τh (`match_trigger_objects`). The filter bits are in
   `TRIGOBJ_BITS` in `analysis/selections/trigger.py`, read from the title of the
   `TrigObj_filterBits` branch of a 2024 NanoAODv15 file.

Because the regions are disjoint, a per-region trigger scale factor can be applied (next section).
The flat OR and every single path remain available as study axes (`pass_<flag>`, `pass_union`, via
`get_trigger_flag_mask` / `get_trigger_flags_or`), which is how the trigger studies compare choices.

**eτh** (see the header of `htautau_etau.yaml` for the numbers):

| region | HLT flag | offline gate on the pair | scale factor |
|---|---|---|---|
| `trig_singleele` | `SingleEle` (Ele30 WPTight) | e pT > 32, e matched (bit 18) | EGM `HLT_SF_Ele30_MVAiso90ID` on the electron |
| `trig_crossele` | `EleTau_PNet` | 25 < e pT ≤ 32, \|η_e\| < 2.1, \|η_τ\| < 2.1, both legs matched (e bit 17, τ bit 12) | TAU `tau_trigger` `etau` on the τh; no SF exists for the Ele24 leg (SF = 1, a normalisation uncertainty is needed) |
| `trig_ditau_extra` | di-tau flags, no EGamma flag fired | monitoring only, **not** in `trigger` (an electron firing an HLT τ leg has no SF) | — |

The single-electron gate is 32, not 30: the 30–32 GeV bin is still on the turn-on (SF ≈ 0.70),
so those electrons go to the cross region. Only the PNet cross path is used, because the 2024 τ-leg
SFs are measured for the PNet paths.

**τhτh** (see the header of `htautau_tautau.yaml`):

| region | HLT flag | offline requirement | scale factor |
|---|---|---|---|
| `trig_ditau` | `DiTau_PNet` | both τh (pT > 24.6, \|η\| < 2.1) matched to the di-tau leg (bit 11) | `ditau`(τ1) × `ditau`(τ2) |
| `trig_ditaujet` | `DiTauJet_PNet`, not in `trig_ditau` | both τh matched to the di-tau+jet leg (bit 14) + a jet (pT > 65, \|η\| < 4.7) matched to the jet leg (bit 17) | `ditaujet`(τ1) × `ditaujet`(τ2), no jet-leg SF |

24.6 GeV is the lower edge of the 2024 `tau_trigger` fit range (below it the correction raises), not
a plateau. The DeepTau di-tau paths have no 2024 SF and are study axes only.

`trig_ditaujet` is exclusive w.r.t. the **matched** `trig_ditau` region, not w.r.t. the DiTau_PNet
decision (checked 2026-10-07): the TAU POG per-leg SFs are unconditional efficiencies from μτh
tag-and-probe (the di-tau path is never vetoed), so neither exclusion makes the SF product exact, and
decision-level exclusion moved ~1 % of the region. Public Run-3 PNet frameworks (uhh-cms/hh2bbtautau)
also use matched-trigger indicators. Possible upgrade: their OR formula with the `eff_data` / `eff_mc`
of `tau_trigger`.

**Before running on data**: check the L1 seed prescales. HLT branches were checked (2026-10-08) in
one file per era (C–I + I v2) of the Tau / EGamma0 / Muon0 PDs: every adopted path exists. A flag
whose paths are **all** missing from a file, an unknown flag and a missing MET filter now raise
(a partially missing flag still only warns).

## Scale factors and calibrations (2024)

| correction | where | applied |
|---|---|---|
| JEC (Summer24Prompt24_V5) + JER (JRV2), Type-1 PuppiMET | `analysis/corrections/jerc.py` (re-propagated relative to the NanoAOD jet pT, since the NanoAOD PuppiMET is already Type-1 corrected) | yes |
| jet veto map (Summer24Prompt24_RunBCDEFGHI_V1) | `jetvetomaps.py` (~8 % of the events) | yes |
| muon scale / smearing | `muon_ss.py` (2024 corrections for 26 < pT < 200 GeV; other muons keep their NanoAOD pT; unseeded random numbers) | yes |
| electron scale / smearing | `electron_ss.py` (supercluster η) | yes |
| tau energy scale (propagated to PuppiMET) | `tau.py`; `corrections.tau_es` must match the τh selection WPs | yes (MC) |
| pileup (official LUM 2026-04-15 `puWeights_CDEFGHI.json.gz`, `Collisions24_CDEFGHI_goldenJSON`, pinned; replaces the private `Collisions2024_preliminary` since 2026-10-08), PS / PDF / scale weights | `pileup.py`, `correction_manager.py` | yes (variations filled only with `fill_variations: true`) |
| muon ID / isolation SFs | `MuonWeights` (MUO `muon_Z.json.gz`, signed-η maps) | in the channel that selects muons (`event_weights.muon`) |
| electron ID (wp90iso) / reco SFs | `ElectronWeights` (supercluster η) | eτh |
| τh ID SFs, DeepTau2018v2p5 VSjet / VSe / VSmu | `TauWeights` on the pair legs (`htt_taus`) | all channels |
| single-electron trigger SF (Ele30) | `ElectronWeights.add_region_trigger_weights`, per trigger region | eτh |
| single-muon trigger SF | `MuonWeights.add_region_trigger_weights`: **a stub to implement** (μτh) | — |
| τ-leg trigger SFs (`tau_trigger`: `ditau`, `ditaujet`, `mutau`, `etau`; PNet paths; DM-binned; from 24.6 GeV) | `TauWeights.add_trigger_weights`, per trigger region | eτh, τhτh |
| c-tag / b-tag SFs | no 2024 set for this selection | **no** |
| object shifts (JES/JER, μ/e scale+res, τ ES Up/Down, unclustered MET) | `object_shifts: true` | off |
| MET recoil corrections (DY / W+jets / H MC), MET xy (φ-modulation) correction, μ / e scale-smearing → MET | not implemented (`update_met` and `apply_met_phi_corrections` exist but are not called) | **no** |

The correction files per year are in `analysis/corrections/correctionlib_files.py` (plus
`jec_params_correctionlib.yaml` for the JEC tags). When jobs die with a JEC `KeyError`, compare the
tags with `/cvmfs/.../latest`.

## Histograms

The common include defines ~120 single-axis histograms: τh leg (kinematics, ID scores, DM,
genPartFlav, pT before TES), di-tau system (m_vis, m_col, m_T^tot, D_ζ, ...), MET family, the three
leading jets (flavours, UParT and PNet scores), the c-jet candidate, multiplicities and gen-level
objects. The channel files add their lepton axes and multi-axis layouts:

- `trig_reco*`, `trig_gen*`, `trig_study`: a turn-on axis (a `<collection>_lead` pT / η) × one Boolean
  axis per HLT flag × `pass_union` × `pass_trigger` × `pass_trig_<region>` (× `gen_channel` or the
  acceptance flags). `plot_trigger_study.py` derives efficiencies, exclusive gains and acceptance
  losses from them.
- `*_by_tauwp`: a few key observables × the tighter DeepTau WP flags (`l2_vs*_tight`).
- `higgs_mass_col` is physical (0 < x₁, x₂ ≤ 1) in only ~25–45 % of the events and NaN otherwise:
  use `higgs_mt_tot` as the default mass proxy for now.
- Axes with `mc_only: true` (gen-level quantities, flavours) are filled with NaN on data.

Added 2026-10-08 (plot folders `ctag/` and `detector/`, placed before `jets/` in `plotting.groups`):

| group | axes / layouts |
|---|---|
| `ctag/` | best-CvL jet among **all** selected jets: `bestcvl_pt/eta/hadronflavour`, UParT `bestcvl_cvsl/cvsb/b_upart`; 2D (CvL, CvB) `ctag2d_bestcvl`, `ctag2d_jet1` (20 × 20); `hf_tag_category` (0 = no c/b tag, 1 = ≥ 1 c-tag and no b-tag, 2 = ≥ 1 b-tag; b has priority); c-jet truth (MC): `cjet_partonflavour`, `cjet_dr_genc` (overflow when there is no hard-process c), `truecjet_pt_rank`, `truecjet_cvl_rank` (0 = highest; overflow = no true c-jet) |
| `detector/` | `rho` (`Rho_fixedGridRhoFastjetAll`), `pu_ntrueint` (MC), η–φ maps `etaphi_tau`, `etaphi_jet1` (+ `etaphi_ele` in eτh, `etaphi_tau1` in τhτh) |
| `tau/` | τ response (MC, genuine τh): `tau_response`, `tau_response_raw` and the per-DM layouts `tau_response_by_dm`, `tau_response_raw_by_dm` (τhτh: also `tau1_response*_by_dm`) |

`plot_shapes.py` draws two-axis layouts (other than `trig_*` / `*_by_tauwp`) as one panel per process,
saved under the layout name: a 2D map (`--logz` for a log colour scale), or, when the second axis is
a category (DM), one unit-area curve per category with the weighted median in the legend.

## Output layout

`outputs/<workflow>/2024/` holds the framework files (per-part `.coffea`, merged `<dataset>.coffea`,
`<process>.coffea`, `2024_processed_histograms.coffea`, `<category>/cutflow_*.csv`). Everything meant to
be looked at is under `outputs/plots/<workflow>/2024/`:

```
outputs/plots/htautau_etau/2024/
├── README.md, index.html                 glossary / categories, browsable thumbnails (open from CERNBox)
├── cutflows/<category folder>.csv
├── trigger_study/trigger_eff.{md,json}, trigger_study/<category folder>/turnon_*.png
├── 1_gen_full_region/                    gen_etau
├── 2_trigger_study_denominator/          den_etau_study
├── 3_selection_no_trigger/               sel_etau
├── 4_selection_trigger/                  sel_etau_trig  <- analysis baseline
│   └── higgs_candidate/ electron/ tau/ met/ jets/ multiplicities/ gen/ tau_wp_study/
└── 5_selection_trigger_truthmatched/     sel_etau_trig_genetau
```

Folder names and variable groups come from the `plotting:` block of the workflow (`category_labels`,
`group_labels`, `groups` = fnmatch patterns on the axis names). "lepton" is pair leg 1 (μ, e or
the leading τh) and "tau" is the τh leg.

## First MC campaign (2026-10-08/09): eτh and τhτh

All MC (31 background samples + cH + bH) ran on condor for both workflows (5 files per job); data stay
blinded. Reports with the full numbers: `DeepTau2PNET/eTau_analysis_feedback_2026-10-08.md` and
`DeepTau2PNET/tauTau_analysis_feedback_2026-10-09.md` (§7: the limit method).

| (108.96 fb⁻¹, MC only, no QCD) | eτh | τhτh |
|---|---|---|
| SR yield H+c / H+b | 5.57 / 42.4 | 7.32 / 56.0 |
| total background | 2.14M (DYll e→τh 56 %, W 21 %, Z→ττ 14 %) | 0.40M (Z→ττ 85 %) |
| S/√B H+c / H+b | 0.0038 / 0.029 | 0.0116 / 0.089 |
| stat-only expected limit (× SM, m_vis × 0/c/b) | 439 / 48.4 (best variant 262 / 32.3) | 92 / 10.9 |

- eτh + τhτh: 90 / 10.6 × SM (86 / 10.3 with eτh vsE Tight + vsJet Tight + mT < 70); 120 / 14.2 with a QCD
  placeholder. Statistics only: S/B is 1e-6 – 1e-4, systematics will dominate unless the backgrounds are
  constrained in situ.
- The fit variable is m_vis (10 GeV bins) split into the 0/c/b tag categories by factorisation: there is no
  joint histogram of two physics variables yet, and no m_ττ (SVfit / FastMTT).
- Open items: QCD multijet (needs data sidebands), the τhτh Z→ττ yield (above published CMS numbers), the
  τhτh trigger turn-on (24.6 GeV threshold), vsE Tight in eτh, the cH generator card (~70 % of cH events
  have no hard c quark), the c-tag light mistag rate.

## Known traps

- Σw, not event counts, for normalisation (negative weights). `run_postprocess.py` refuses to
  normalise when a part listed in `local_parts.txt` is missing.
- Boolean selections / axes must never contain `None` (NaN casts to `True`): `ak.fill_none(..., False)`.
- All axes of a multi-axis layout must be per-event scalars (`objects['l1']`, `ak.firsts(...)`).
- Collections are not pT-ordered after the energy corrections: use `<collection>_lead` / `_sublead`.
- Every axis in `axes` is evaluated on every chunk, even without a layout.
- Scale factors apply only in `object_sf_categories` (`sel_*`); a trigger SF block also needs its own
  `categories`. The pT lower edges of the SF maps (Ele30 25 GeV, `tau_trigger` 24.6 GeV) raise or
  are clamped with a warning: the object selection must not admit lower-pT legs into a region.
- `run_postprocess.py --plot` (data/MC stacks) needs an MC process; for the signal-only stage use
  `scripts/htautau/plot_shapes.py`.
- Status of the 2026-10-08 changes (electron veto, pileup, gen visible τ, acceptance, new variables,
  backgrounds, `data_blinded`, `mc_only_categories`, cutflow rows): implemented, validated on the
  pilot / one file per sample, and used by the 2026-10-08 condor campaign (signals included). The
  older outputs were moved to `outputs/<wf>/2024_oldsignal_20260916`. Do not mix per-part outputs from
  before and after (the cutflow `initial` row changed meaning). `plot_shapes.py` on outputs produced
  before 2026-10-08 skips the new layouts with one warning.
- Condor: a job that fails at the same job number in several workflows is a bad file replica, not a
  transient error. Find the file in the `.err`, test-open it through a redirector
  (`root://xrootd-cms.infn.it/` or `root://cms-xrd-global.cern.ch/`), reroute it in the fileset JSON
  **and** in every workflow's `condor/<wf>/2024/<dataset>/partitions.json`, remove the queued copies and
  resubmit only the missing job numbers with the same partitions (never `jobs_status.py` for this: it
  re-fetches and renumbers the parts). A redirector can still send the job back to the bad replica:
  then pin a good site URL.
- The stacked plotter (`run_postprocess.py --plot`) skips the x-range zoom for an all-empty histogram
  (e.g. `gen_lepton_*` in the τhτh gen category) instead of crashing.
- Do not edit processor-side code while `run_signal_local.sh` or condor jobs run: every part/job
  imports the current code.

---

## The μτh channel: your starting point

`analysis/workflows/htautau_mutau.yaml` is a **runnable skeleton**. Everything that is not written
in it comes from the common include. Right now it:

- runs on the pilot file and fills `gen_mutau`, the truth μτh events. Look at them first: gen muon
  and τh spectra, the visible di-tau mass, the heavy-flavour quark;
- builds **no pair**: `select_htt_pair_mutau` is a placeholder that stores an empty pair
  collection and warns once, so `sel_mutau` / `sel_mutau_trig` stay empty;
- uses the single-muon flag as a placeholder `trigger`, with no offline gate and no matching;
- applies no muon scale factors.

`htautau_etau.yaml` is a complete channel built on exactly the same machinery: read it next to this
guide. The skeleton's numbered TODOs correspond to these steps (a sensible order):

| step | what | where |
|---|---|---|
| 1 | **Objects.** The signal muon: kinematics consistent with the triggers you will use, ID, isolation, impact parameter. Check that `taus` suits you. Override `muons` (and later `muons_study` / `taus_study`) in the yaml. | `htautau_mutau.yaml`; `analysis/working_points/working_points.py` |
| 2 | **The pair.** Implement `select_htt_pair_mutau` and a ranking recipe `RANK_MUTAU`. How do you pick one (μ, τh) pair when there are several? | `analysis/selections/object_selections.py`, `analysis/selections/htautau.py` (`build_htt_pair`, `RANK_ETAU`) |
| 3 | **Triggers.** Which single-muon and μ–τh cross triggers exist in the 2024 menu, with which online thresholds? Add their flags. Find the TrigObj filter bits of their muon and τ legs (`TRIGOBJ_BITS`). Define gated trigger regions and make `trigger` their OR. Add `pass_<flag>` / turn-on axes to study the efficiencies. | `analysis/selections/trigger_flags.yaml`, `analysis/selections/trigger.py`, `htautau_mutau.yaml` |
| 4 | **Event selection.** Exactly one signal muon, the extra-lepton vetoes, ≥ 1 τh, the categories of the naming scheme above (incl. a trigger-study denominator). | `htautau_mutau.yaml` |
| 5 | **Scale factors.** Muon ID / isolation SFs matching your working points (`event_weights.muon`); the single-muon trigger SF (implement `MuonWeights.add_region_trigger_weights`, modelled on the electron one); the τ-leg trigger SF of a cross region (`tau.trigger`, trigtype `mutau`). Which (pT, η) ranges do the maps cover? | `htautau_mutau.yaml`, `analysis/corrections/muon.py`, `analysis/corrections/tau.py` |
| 6 | **Histograms and acceptance.** Muon axes (isolation, IDs, genPartFlav, ...), the `trig_*` layouts, `in_gen_acc_mutau` in `select_acceptance_flags`, the `plotting:` labels. | `htautau_mutau.yaml`, `analysis/selections/object_selections.py` |

After every step: config parse → smoke test on the pilot file (check the printed cutflows) → full
local run → `run_postprocess.py` → `plot_shapes.py` / `plot_trigger_study.py`, and compare with the
gen-level category. Keep a copy under another workflow name (e.g. `htautau_mutau_v2.yaml`) when
you want to compare two versions, since a rerun overwrites the outputs of the same workflow name.
Documenting *why* each threshold / path / working point was chosen (turn-on plots, SF ranges,
expected yields) matters as much as the code.
