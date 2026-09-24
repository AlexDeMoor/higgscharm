# higgscharm — cH, H→ττ analysis framework

[![Codestyle](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Columnar CMS analysis package ([coffea](https://coffeateam.github.io/coffea/) 0.7 +
[scikit-hep](https://scikit-hep.org)) that turns NanoAOD into histograms. The analysis on this branch
is **cH, H→ττ + charm** (Higgs boson produced with a charm quark) on 2024 data / Summer24 MC
(NanoAODv15), in three channels: μτh, eτh and τhτh. The signal (cH, H→ττ) and the main
background (bH, H→ττ) are private FullSim productions; the DAS backgrounds and the data come later.

| read | for |
|---|---|
| this file | setup, repository map, the pipeline, how to run everything |
| [`README_HToTauTau.md`](README_HToTauTau.md) | the H→ττ analysis itself: samples, objects, pair building, triggers, scale factors, categories, histograms, plots, **and the μτh starting point** |
| [`analysis/workflows/README.md`](analysis/workflows/README.md) | the workflow YAML reference (every key, the expression scopes, `extends`, `plugins`, weights) |

- [1. Getting started on lxplus](#1-getting-started-on-lxplus)
- [2. Repository map](#2-repository-map)
- [3. The pipeline](#3-the-pipeline)
- [4. Input filesets](#4-input-filesets)
- [5. Workflows](#5-workflows)
- [6. Running](#6-running)
- [7. Postprocessing and plots](#7-postprocessing-and-plots)
- [8. Development loop and traps](#8-development-loop-and-traps)
- [9. Other workflows and legacy tooling](#9-other-workflows-and-legacy-tooling)

---

## 1. Getting started on lxplus

```bash
ssh <user>@lxplus.cern.ch
git clone https://github.com/AlexDeMoor/higgscharm.git && cd higgscharm
git switch Htautau-analysis              # the cH, H->tautau branch
git switch -c <your-branch>              # develop on your own branch, open PRs against Htautau-analysis
```

**Python environment.** The code needs coffea **0.7.x** (0.7.30), not the coffea 2025 of a recent
`python3`. The simplest environment is the image the condor jobs use, through apptainer:

```bash
export IMG=/cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base-almalinux9:0.7.30-py3.10
export PY="apptainer exec -B /afs -B /eos -B /cvmfs -B /run/user/$(id -u) $IMG python3"   # add to your ~/.bashrc
$PY -c "import coffea; print(coffea.__version__)"                    # 0.7.30
```

(`-B /run/user/<uid>` passes your Kerberos ticket into the container; without it EOS paths give
`Permission denied` inside the container.)

A conda environment with coffea 0.7.30, dask, hist, mplhep, correctionlib and XRootD works too (set
`PY` to its `python3`). Every command below is run **from the repository root** (the scripts use
paths relative to the working directory).

**Output area.** All outputs go to your EOS space; `outputs/` in the repo must point there:

```bash
mkdir -p /eos/user/${USER:0:1}/$USER/higgscharm/outputs
ln -sfn /eos/user/${USER:0:1}/$USER/higgscharm/outputs outputs
```

**Grid proxy**: only needed for DAS datasets and condor (`voms-proxy-init --voms cms --valid 192:00`).
The private signal samples are read from EOS directly.

**First run** (a few minutes):

```bash
$PY fetch.py -y 2024 --onlysignal              # fileset JSON of the private signal samples (no proxy needed)
PILOT=/eos/user/a/ademoor/HcTT/bH_fullsim_Summer24_v1_merged/tranche0_test/HbTT_Summer24_NanoAODv15_merged_0000_0001.root
$PY run_local.py -w htautau_etau  -y 2024 -d HplusBottom_HtoTauTau --files $PILOT --iterative
$PY run_local.py -w htautau_mutau -y 2024 -d HplusBottom_HtoTauTau --files $PILOT --iterative
```

The pilot file (5011 bH events, not part of the fileset) is the standard smoke-test input. The
first command runs the complete eτh channel; the second one runs the μτh starting point (it warns
that the μτh pair is not built yet, see [`README_HToTauTau.md`](README_HToTauTau.md#the-μτh-channel-your-starting-point)).

---

## 2. Repository map

| path | content |
|---|---|
| `analysis/workflows/*.yaml` | **the analysis definitions** (one YAML per workflow; `includes/` holds shared bodies) |
| `analysis/processors/base.py` | `BaseProcessor`: the single coffea processor that executes any workflow YAML |
| `analysis/selections/` | `object_selections.py` (`ObjectSelector`: every non-trivial object builder), `htautau.py` (H→ττ helpers: pair building, gen classification, m_col, D_ζ), `trigger.py` + `trigger_flags.yaml` (HLT flags, trigger regions, trigger-object matching), `event_selections.py` (lumi mask, MET filters, ...) |
| `analysis/working_points/` | ID / isolation / tagging working points (`working_points.muon_id(events, 'medium')`, `tau_id_vsjet`, `jet_ctagging`, ...) |
| `analysis/corrections/` | object corrections (JEC/JER, jet veto, lepton scale/smearing, tau energy scale) and event weights (pileup, PS/PDF/scale, lepton and tau ID / trigger scale factors); `correctionlib_files.py` lists the correction files per year |
| `analysis/histograms/` | histogram building / filling from the `histogram_config` block |
| `analysis/filesets/` | dataset definitions (`<year>_<nano>.yaml`) and the generated fileset JSONs |
| `analysis/postprocess/` | merging, normalisation (lumi × xsec / Σw), cutflows, data/MC plots |
| `run_local.py` | run a workflow locally (one dataset or explicit files) |
| `runner.py`, `submit_condor.py`, `submit.py`, `condor/`, `jobs_status.py` | HTCondor submission and monitoring |
| `run_postprocess.py` | merge job outputs, normalise, cutflows |
| `scripts/htautau/` | local production of the signal samples (`run_signal_local.sh`) and the H→ττ plots (trigger study, shapes, index pages) |
| `fetch.py` | builds the fileset JSONs (DAS/Rucio, or the `path:` globs of private samples) |

---

## 3. The pipeline

```
analysis/filesets/2024_nanov15.yaml          dataset -> {era, process, key, xsec, query | path}
   |  fetch.py                                -> analysis/filesets/fileset_2024_nanov15_lxplus.json (generated, not in git)
   v
run_local.py  |  scripts/htautau/run_signal_local.sh  |  runner.py -> condor
   |  BaseProcessor + analysis/workflows/<workflow>.yaml
   v
outputs/<workflow>/2024/<dataset>/<dataset>[_<part>].coffea     histograms (with a `category` axis) + metadata (Σw, cutflows)
   |  run_postprocess.py --postprocess --output_format coffea     (normalised to lumi × xsec / Σw)
   v
outputs/<workflow>/2024/<process>.coffea, 2024_processed_histograms.coffea, <category>/cutflow_*.csv
   |  scripts/htautau/plot_trigger_study.py, plot_shapes.py, make_index.py
   v
outputs/plots/<workflow>/2024/...   PNGs per category / variable group, trigger-study plots, README.md, index.html
```

Inside `BaseProcessor.process`, for every chunk of events:

1. **object corrections** (`corrections.object`): JEC/JER (+ Type-1 PuppiMET), jet veto map, muon /
   electron scale and smearing, tau energy scale;
2. **object selection** (`object_selection`): the objects are built one after the other in YAML order
   into a dictionary `objects` (plain NanoAOD collections with cuts, or `ObjectSelector` methods);
3. **event selection** (`event_selection.selections`): named boolean masks, combined into
   **categories** (analysis regions);
4. **weights** (`corrections.event_weights`): genWeight, pileup, PS/PDF/scale, scale factors;
5. **histograms** (`histogram_config`): every axis expression is evaluated once, then the
   histograms are filled per category.

---

## 4. Input filesets

Each campaign has a dataset file `analysis/filesets/<year>_<nano version>.yaml`
(2024: `2024_nanov15.yaml`):

```yaml
HplusCharm_HtoTauTau:              # dataset name (used with run_local.py -d)
  era: signal                      # data: run letter (C, D, ...); MC: mc; private signal: signal
  query: null                      # DAS query (DAS datasets)
  path: /eos/user/a/ademoor/HcTT/fullsim_Summer24_v1_merged/tranche[12]/*.root   # glob(s) (private samples)
  process: H+c                     # grouping at postprocessing (and plot colours)
  key: hc_tautau                   # referenced by the `datasets:` block of a workflow
  xsec: 0.006082                   # pb (null for data)
```

`fetch.py` turns these into the fileset JSON read by `run_local.py` / condor:

```bash
$PY fetch.py -y 2024 --onlysignal      # only the `path:` (private) samples: no proxy, no container
$PY fetch.py -y 2024                   # DAS datasets (needs a VOMS proxy; runs in the coffea-dask image, --image)
```

Options: `-y {2016preVFP, ..., 2023postBPix, 2024}`, `--samples <names>` (a subset), `-s/--signal`
(also the signal), `--site` (redirector for private samples). The JSON is generated, not versioned.

---

## 5. Workflows

A workflow is one YAML file in `analysis/workflows/` that defines the whole analysis of a region:
`datasets` (fileset keys), `object_selection`, `event_selection` (`hlt_paths`, `selections`,
`categories`), `corrections` (object corrections + event weights), `histogram_config` (axes +
layout), optionally `plotting` (plot folder names), `plugins` (extra Python modules) and
`extends` (inherit a shared body from `includes/`). Cut and axis strings are Python expressions
evaluated by the processor. The full reference is [`analysis/workflows/README.md`](analysis/workflows/README.md).

H→ττ workflows (all `extends: includes/htautau_common`):

| workflow | channel | status |
|---|---|---|
| `htautau_etau` | e τh | complete: objects, gated trigger regions with trigger-object matching, trigger / ID scale factors, trigger-study histograms |
| `htautau_tautau` | τh τh | complete |
| `htautau_mutau` | μ τh | **starting point** (runnable skeleton with TODOs), see [`README_HToTauTau.md`](README_HToTauTau.md#the-μτh-channel-your-starting-point) |

The other workflows (Z→ℓℓ, ZZ→4ℓ, H(→ZZ)+c, c-tag efficiencies) are listed in [section 9](#9-other-workflows-and-legacy-tooling).

---

## 6. Running

### Smoke test (always first)

```bash
$PY run_local.py -w <workflow> -y 2024 -d HplusBottom_HtoTauTau --files $PILOT --max-chunks 1 --iterative
```

`--iterative` runs in one process (full Python tracebacks); `--max-chunks 1` stops after one
chunk of 50k events. The output goes to `outputs/<workflow>/2024/HplusBottom_HtoTauTau/` and the
cutflow of every category is printed at the end.

### Full signal samples, locally

```bash
bash scripts/htautau/run_signal_local.sh "htautau_etau htautau_mutau" "HplusCharm_HtoTauTau HplusBottom_HtoTauTau"
```

The samples (≈20M events) are copied in parts to `/dev/shm` and the workflows run on each part in
parallel (~45 min for three workflows). Reading full samples straight from EOS/xrootd on lxplus is
latency bound (~10 MB/s per node, `Operation expired` timeouts), which is why the script stages
the files. It writes one `<dataset>_<part>.coffea` per part plus `local_parts.txt`; the
`<dataset>_*.coffea` of a previous run of the same workflow are **deleted** first. Environment
knobs: `PY`, `NPART` (files per part), `WORKERS`. Run long jobs inside `tmux`.

`run_local.py` options: `-d <dataset>` (from the fileset JSON), `--files` (explicit files),
`--max-files`, `--max-chunks`, `--workers`, `--chunksize`, `--stage-dir` (copy the inputs locally
first), `--eos` (write to your EOS outputs; same place as the symlink), `--tag` (suffix of the output
file name).

### Condor (DAS backgrounds / data)

```bash
voms-proxy-init --voms cms --valid 192:00
$PY runner.py -w <workflow> -y 2024 --submit --eos --output_format coffea -m 4000   # without --submit: only builds the job files
$PY jobs_status.py -w <workflow> -y 2024 --eos --output_format coffea               # missing / failed jobs, resubmission
```

Jobs run in the `coffea-base-almalinux9:0.7.30-py3.10` image and `cd` into this repository on EOS,
so they run the code as it is **at execution time**: do not edit it while jobs are queued, and
validate one job before a campaign. `submit_condor.py` submits a single dataset.

---

## 7. Postprocessing and plots

```bash
$PY run_postprocess.py -w <workflow> -y 2024 --postprocess --output_format coffea
$PY scripts/htautau/plot_shapes.py        -w <workflow> -y 2024     # normalised shapes, per category / variable group
$PY scripts/htautau/plot_trigger_study.py -w <workflow> -y 2024     # trigger efficiencies / turn-ons from the trig_* histograms
$PY scripts/htautau/make_index.py         -w <workflow> -y 2024     # README.md + index.html of the plot folder
```

- `run_postprocess.py` merges the job outputs per dataset, normalises MC to `lumi × xsec / Σw`
  (Σw = the generator weight sum stored in the job metadata; never event counts, the signal
  genWeights are ±constant with 20–30 % negative) and groups datasets into processes. Outputs:
  `outputs/<wf>/2024/<process>.coffea`, `2024_processed_histograms.coffea`, and per category
  `cutflow_<category>.csv` (weighted) / `cutflow_raw_<category>.csv` (event counts).
- The plots go to `outputs/plots/<workflow>/2024/`, with one numbered folder per category (names from
  the `plotting:` block of the workflow) and one sub-folder per variable group. Open `index.html`
  from CERNBox to browse them.
- `run_postprocess.py --plot` makes data/MC stack plots; it is meant for when backgrounds and data
  are available.

---

## 8. Development loop and traps

There is no test suite. The loop is:

1. edit the workflow YAML (or the Python it calls);
2. check that the config builds:
   `$PY -c "from analysis.workflows.config import WorkflowConfigBuilder as B; B('<workflow>').build_workflow_config()"`;
3. smoke test on the pilot file (section 6) and read the printed cutflows;
4. full local run, postprocess, plots.

Traps worth knowing early (more in [`README_HToTauTau.md`](README_HToTauTau.md#known-traps)):

- Objects are built **in YAML order**; an object can only use objects defined above it. With
  `extends`, an overridden key keeps its place and new keys are appended at the end.
- **Every axis in `histogram_config.axes` is evaluated on every chunk**, whether or not a layout
  uses it: an axis that refers to a missing object or selection breaks the run.
- Boolean selections / axes must never contain `None` (use `ak.fill_none(..., False)`); multi-axis
  histograms need per-event scalars (`objects['l1'].pt`, `ak.firsts(...)`), not jagged arrays.
- Collections are not pT-ordered after the energy corrections: use the `<collection>_lead` /
  `_sublead` objects for "leading object" quantities.
- The muon scale / smearing correction uses unseeded random numbers, so two runs of a μ workflow
  differ slightly in the muon-dependent histograms (the selections and cutflows do not change).
- Rerunning a workflow overwrites its `outputs/<workflow>/2024/<dataset>/` files; use a different
  workflow name (e.g. a copy `htautau_mutau_test.yaml`) to keep an old result.

---

## 9. Other workflows and legacy tooling

These predate the H→ττ analysis and share the same processor (they are useful as examples):

* Z→ℓℓ calibration regions: `ztomumu`, `ztoee` (automated postprocessing: `run_full_postprocess.py`)
* H→ZZ→4ℓ: `zzto4l`; Z+X estimation: `zplusl_os`, `zplusll_os`, `zplusl_ss*`, `zplusll_ss`
* H(→ZZ→4ℓ)+c: `hplusc`; c-tagging efficiency maps: `ctag_eff`

`scripts/mva` (BDT training / inference) and `scripts/combine` (datacards, fits, impacts) are
generic tools inherited from an earlier H→WW analysis; their usage examples still name H→WW
workflows that no longer exist in this branch.
