## Workflow config

A workflow is one yaml file `analysis/workflows/<name>.yaml` (run with `-w <name>`). It defines the
selection, variables, output histograms, triggers, corrections and weights of an analysis region. The
top-level keys are `extends`, `plugins`, `datasets`, `object_selection`, `event_selection`,
`corrections`, `histogram_config` and `plotting` (plus `mva` / `inference` / `combine` for the legacy
MVA and datacard tools). The strings in `cuts`, `selections` and axis `expression`s are **Python
expressions** evaluated by the processor (see [Expression scopes](#expression-scopes)).

### `extends` (shared bodies)

A workflow may start with `extends: includes/<name>`: the yaml `analysis/workflows/includes/<name>.yaml` is loaded first and the workflow is deep-merged on top of it (`WorkflowConfigBuilder`):
- keys with dict values are merged recursively, any other value (lists included) is replaced;
- a `null` value deletes the base key; base key order is preserved and new keys are appended (this matters for `object_selection`, whose objects are built sequentially);
- files under `includes/` are not runnable workflows (they are not listed by `runner.py` / `submit.py`).

The H→ττ workflows (`htautau_mutau`, `htautau_etau`, `htautau_tautau`) share `includes/htautau_common.yaml` this way.

### `plugins` (extra Python modules)

```yaml
plugins:
  - analysis.selections.my_module     # imported by the processor before any event is processed
```

Every listed module is imported (`importlib.import_module`) when `BaseProcessor` is created and again
in each worker, before the first chunk. At import time a plugin can extend the framework without
touching its files, for example:
- add object builders: `ObjectSelector.select_my_object = my_function`, where `my_function(self, obj_name)`
  stores its result in `self.objects[...]`, and then use `field: select_my_object`;
- register trigger flags: `analysis.selections.trigger.register_trigger_flags(2024, {"MyFlag": ["HLT_path_without_prefix"]})`.

An import error stops the run. For code that belongs to the framework, edit the framework modules
directly; plugins are meant for work-in-progress or workflow-specific code.

### `datasets`

Data and MC datasets to process

```yaml
datasets:
  data:               # List of dataset keys corresponding to real data
    - electron
  mc:                 # List of dataset keys corresponding to Monte Carlo samples
    - dy_inclusive
    - singletop
    - tt
    - wjets_ht
    - diboson
  signal:             # List of dataset keys corresponding to signal samples
    - hplusc
    - hplusb
    - hpluslight
```
- Each item refers to the `key` field defined in the fileset YAML
- This section is used by the job submission script `runner.py` to determine which samples are to be run
- `data_blinded` (optional group): data keys that are defined and validated but **never processed**: `runner.py`, `jobs_status.py` and the postprocess only read `data`, `mc` and `signal`, and `run_postprocess.py` forces `--blind` when there is no `data` group. The H→ττ channels list their primary datasets here; rename the group to `data:` to unblind.
- `get_datasets_to_run_over` raises on a key that no fileset yaml of any year defines (typo) and on an unknown group; keys that exist only for other years are skipped with one summary warning.
- With `extends`, the groups are deep-merged: the H→ττ common include holds `signal` and `mc`, each channel adds its `data_blinded`.


* `object_selection`: Contains the information required for object selection:
```yaml
object_selection:
  muons:
    field: events.Muon
    cuts:
      - events.Muon.pt > 10
      - np.abs(events.Muon.eta) < 2.4
      - working_points.muon_iso(events, 'tight')
      - working_points.muon_id(events, 'tight')
  electrons:
    field: events.Electron
    cuts:
      - events.Electron.pt > 10
      - np.abs(events.Electron.eta) < 2.5
      - working_points.electron_id(events, 'wp80iso')
      - delta_r_higher(events.Electron, objects['muons'], 0.4)
  dimuons:
    field: select_dimuons
    cuts:
      - objects['dimuons'].l1.delta_r(objects['dimuons'].l2) > 0.02
      - objects['dimuons'].l1.charge * objects['dimuons'].l2.charge < 0
      - (objects['dimuons'].p4.mass > 60.0) & (objects['dimuons'].p4.mass < 120.0)
```
With `field` you define how to select the object, either through a NanoAOD field (`events.Muon`) or a custom object-selection function (`select_dimuons`) defined as a method of the `ObjectSelector` class located at `analysis/selections/object_selections.py`. A `field` that contains the string `events` is evaluated as an expression; anything else is taken as a method name. Each object is added sequentially (in yaml order) to a dictionary called `objects`, which can later be used to access the objects: an object can only use the objects defined above it. Methods may also write several entries at once. The H→ττ methods (`select_htt_pair_<channel>`, `select_htt_legs`, `select_trigger_legs`, `select_htt_met_vars`, `select_gen_htautau`, `select_acceptance_flags`, ...) and what they produce are listed in [`README_HToTauTau.md`](../../README_HToTauTau.md#the-di-tau-pair-and-the-derived-objects).

`cuts` defines the set of object-level cuts to apply to the object. Similarly, you can use NanoAOD fields (`events.Muon.pt > 24`) to define a cut or any valid expression using the already defined objects (`objects['dimuons'].z.mass < 120.0`). Alternatively, you can also use a working point function (`working_points.muon_iso(events, 'tight')`) defined in `analysis/working_points/working_points.py)`

`add_field` adds computed fields to an object (evaluated before its `cuts`), e.g.
`add_field: {rel_iso: 'events.Muon.pfRelIso04_all'}`; the expressions see `events`, `objects` and `year`.

You can also use `add_cut` to define masks that will be added to the object collection as a new field and can be accessed later in the workflow:

```yaml
muons:
    field: events.Muon
    add_cut:
      is_loose:
        - events.Muon.pt > 5
        - np.abs(events.Muon.eta) < 2.4
        - events.Muon.isGlobal | (events.Muon.isTracker & (events.Muon.nStations > 0))
      is_relaxed:
        - objects['muons'].is_loose
        - np.abs(events.Muon.sip3d) < 4
      is_tight:
        - objects['muons'].is_loose
        - objects['muons'].is_relaxed
        - events.Muon.isPFcand | ((events.Muon.highPtId > 0) & (events.Muon.pt > 200))
zcandidates:
    field: select_zcandidates 
    add_cut:
      is_ossf:
        - objects['zcandidates'].l1.pdgId == -objects['zcandidates'].l2.pdgId
      is_ss:
        - objects['zcandidates'].l1.pdgId == objects['zcandidates'].l2.pdgId
      is_sr:
        - objects['zcandidates'].is_ossf
        - (1*objects['zcandidates'].l1.is_tight + 1*objects['zcandidates'].l2.is_tight) == 2
      is_sscr:
        - objects['zcandidates'].is_ss
        - objects['zcandidates'].l1.is_relaxed
        - objects['zcandidates'].l2.is_relaxed
```

### `event_selection`

Defines event-level selections and analysis regions (categories)
```yaml
event_selection:
  hlt_paths:
    Muon:                        # Primary dataset name (Muon, EGamma, Tau, ...)
      - SingleMu                 # Trigger flag (defined in analysis/selections/trigger_flags.yaml)
  selections:                    # Event-level selections
    trigger: get_trigger_mask(events, hlt_paths, dataset, year)
    trigger_match: get_trigger_match_mask(events, hlt_paths, year, events.Muon)
    lumi: get_lumi_mask(events, year)
    goodvertex: events.PV.npvsGood > 0
    two_muons: ak.num(objects['muons']) == 2
    one_dimuon: ak.num(objects['dimuons']) == 1
    leading_muon_pt: ak.firsts(objects['muons'].pt) > 30
    subleading_muon_pt: ak.pad_none(objects['muons'], target=2)[:, 1].pt > 15
  categories:
    base:                        # Named selection region
      - goodvertex
      - lumi
      - trigger
      - trigger_match
      - two_muons
      - leading_muon_pt
      - subleading_muon_pt
      - one_dimuon
```

- `hlt_paths`: Maps each primary dataset (the prefix of the data dataset names: `Muon`, `EGamma`, `Tau`, ... see `get_dataset_name` in `analysis/filesets/utils.py`) to the list of trigger flags taken from it
    - Trigger flags are defined per year in [`analysis/selections/trigger_flags.yaml`](../selections/trigger_flags.yaml) (a plugin can add more with `register_trigger_flags`)
    - For data: an event is taken from the first primary dataset (in dict order) whose flags fired, so that no event is counted twice; a data primary dataset missing from `hlt_paths` raises.
    - For MC: all flags are combined with a logical OR via `get_trigger_mask()` ([`trigger.py`](../selections/trigger.py)). Alternatively, you can write your own trigger function to be used.
    - The older workflows enforce lepton–trigger-object matching with `get_trigger_match_mask()`; the H→ττ workflows match the pair legs themselves (`select_trigger_legs`, see below).

- `selections`: Defines event-level cuts. Similarly to object selection, you can use any valid expression from a NanoAOD field or a custom event-selection function defined at [`analysis/selections/event_selections.py`](../selections/event_selections.py)

- `categories`: Named groups of selections that define analysis regions. Each category defines a region for histogram filling and postprocessing (a workflow must define at least one). The per-category masks are exposed to the histogram expressions as `category_masks[<category>]` (e.g. a Boolean "passes the reco selection" axis inside a gen-level category).
- Selections are evaluated in order (all of them, whether a category uses them or not) and every selection may refer to the ones defined before it as `selection_masks['<name>']` (e.g. `trigger: selection_masks['trig_singleele'] | selection_masks['trig_crossele']`); the same dictionary is in scope for the histogram-axis expressions. A selection must be a per-event boolean without `None` (`ak.fill_none(..., False)`); at most 32 selections per workflow.
- Trigger helpers in the event-selection scope: `get_trigger_flag_mask(events, '<flag>', year)` (one HLT path), `get_trigger_flags_or(events, [flags], year)` (pure OR, no primary-dataset exclusivity — for study axes) and `get_trigger_region_mask(events, hlt_paths, dataset, year, [flags])` (OR of the flags of one trigger region, PD-exclusive on data following the `hlt_paths` order). `get_trigger_mask` raises for a data primary dataset that is not listed in `hlt_paths`. Trigger-object matching flags come from the `select_trigger_legs` object method (`analysis/selections/trigger.py: match_trigger_objects`, `TRIGOBJ_BITS`).
- Fail-loud checks: an undefined trigger flag (or a year without flags) raises; a flag set to `null` / `[]` for a year (e.g. `TriMu` 2018) is an all-False mask. HLT paths of a flag that are absent from the input file are skipped with a one-time warning as long as at least one path of the flag is present; if all of them are absent the job raises. `get_metfilters_mask` raises when a MET filter of `analysis/data/met/metfilters.json` (or the year) is missing.

- `mc_only_categories` (optional): list of category names or `fnmatch` patterns (`gen_*`) that are not filled for data — no histograms, no parquets and no cutflow in data outputs (their bins stay empty in the data histograms). Use it for the categories that are meaningless on data (gen-level, object-level study denominators without lumimask / MET filters) or double counted on data (a selection without a trigger requirement is filled once per primary dataset). Every entry must match at least one category (checked at config parse time).
```yaml
event_selection:
  mc_only_categories: [gen_etau, den_etau_study, sel_etau, sel_etau_trig_genetau]
```

- Cutflow (`metadata[<category>]["cutflow"]` genWeight-weighted, `["cutflow_raw"]` unweighted; nominal shift only): `initial` = the chunk before the jet veto map (the same events as `metadata["sumw"]`), `jet_veto` = after the jet veto map (equal to `initial` when `jetveto` is not in `corrections.object`), then one row per cut of the category. Data: weight 1.



### `corrections`

Contains the object-level corrections and event-level weights to apply:

```yaml
corrections:
  object:
    - jet                 # Jet scale and resolution (+ Type-1 PuppiMET)
    - jetveto             # Jet veto
    - muon                # Muon scale and resolution  
    - electron            # Electron scale and resolution 
    - tau                 # Tau energy scale (MC; propagated to PuppiMET)
  object_shifts: false    # include object-level shifts 
  tau_es:                 # working points assumed by the tau energy scale (object correction `tau`)
    vsjet: Medium
    vse: VVLoose
  event_weights:
    genWeight: true
    pileupWeight: true
    partonshowerWeight: true
    lhepdfWeight: true
    lhescaleWeight: true
    nnlopsWeight: true
    muon:
      - id: tight
      - iso: tight
      - trigger: true
    electron:
      - id: false
      - reco: false
      - trigger: false
    tau:                  # DeepTau2018v2p5 ID scale factors on the tau_h legs of the selected pair (`selected_htt_taus`)
      - vsjet: Medium
      - vse: VVLoose
      - vsmu: VLoose
      - trigger:          # optional: TAU POG `tau_trigger` SFs on the same legs (see below)
          regions: {ditau: trig_ditau, ditaujet: trig_ditaujet}
          categories: ["sel_tautau_trig*"]
```

The lepton scale factors act on the collections named exactly `muons` / `electrons`, the tau ones on the
τh legs of the selected pair (`htt_taus`).

`muon.trigger` / `electron.trigger` accept either `true` (legacy: the SF on every selected lepton, not for
2024 electrons) or the same region form as `tau.trigger`, e.g.
`trigger: {regions: {Ele30: trig_singleele}, categories: ["sel_etau_trig", "sel_etau_trig_genetau"]}`:
the single-lepton SF of that path (EGM Ele30 with the `id` WP) is applied to the leading selected lepton of
the events in that selection only (`CMS_eff_e_trigger_Ele30_<year>`); pT is clamped to the lower map edge
(25 GeV) with a warning. The electron version is `ElectronWeights.add_region_trigger_weights`
(`analysis/corrections/electron.py`); the muon one, `MuonWeights.add_region_trigger_weights`
(`analysis/corrections/muon.py`), is still a stub to implement.

`tau.trigger` maps a `tau_trigger` trigtype (`ditau`, `ditaujet`, `mutau`, `etau`) to the name of an
event selection (the trigger region). The regions must be disjoint: an event in a region is weighted by
the product of the DM-binned SFs of its τh legs for that trigtype, every other event by 1. The weight is
added only in the categories matching `categories` (fnmatch patterns; they must also match
`object_sf_categories`), as `CMS_eff_trigger_t_<year>` with one nuisance per trigtype and decay mode
(`CMS_eff_trigger_t_<trigtype>_dm<dm>_<year>`). The 2024 correction is binned from 24.6 GeV and raises
below: the τh object selection must not admit lower-pT legs into a region.

- `object_shifts`: If true, object-level systematic shifts (e.g., JES_up, JES_down, JER_up, JER_down, etc) are automatically applied and propagated to downstream event-level calculations.
- Corrections are implemented through a set of utilities managed by two functions:
    - `object_corrector_manager` ([`correction_manager.py`](../corrections/correction_manager.py)): applies object-level corrections (jet energy scale and resolution, jet veto maps, lepton scale and smearing, tau energy scale).
    - `weight_manager` (same file): applies event-level weights (generator weight, pileup, parton shower / PDF / scale, lepton and tau scale factors).

**Note**: Ensure that the working points used for object selection and event-level corrections (ID, isolation, etc.) are consistent!

### `histogram_config`

Use to define processor's output histograms (more info on Hist histograms [here](https://hist.readthedocs.io/en/latest/)). Here you define the histogram axes associated with the variables you want to include in the analysis. 
```yaml
histogram_config:
  add_syst_axis: true
  add_weight: true
  axes:
    muon_pt:
      type: Regular
      bins: 50
      start: 30
      stop: 300
      label: $p_T(\mu)$ [GeV]
      expression: objects['muons'].pt
    muon_eta:
      type: Regular
      bins: 50
      start: -2.5
      stop: 2.5
      label: $\eta(\mu)$
      expression: objects['muons'].eta
    muon_phi:
      type: Regular
      bins: 50
      start: -3.14159
      stop: 3.14159
      label: $\phi(\mu)$
      expression: objects['muons'].phi
    dimuon_mass:
      type: Regular
      bins: 100
      start: 10
      stop: 150
      label: $m(\mu\mu)$ [GeV]
      expression: objects['dimuons'].z.p4.mass
  layout:
    muon:
      - muon_pt
      - muon_eta
      - muon_phi
    zcandidate:
      - dimuon_mass
```
`fill_variations: false` fills only the nominal weight (no Up/Down templates of the weight systematics), which makes shape/efficiency studies much faster; it is forced to `false` when `add_syst_axis` is `false`. An axis may carry `mc_only: true`: on data it is filled with NaN (flow bin) instead of being evaluated (gen-level quantities, flavours). `event_weights.object_sf_categories` (list of fnmatch patterns on the category names, default all) restricts the lepton / tau ID scale factors to those categories.

Note that the variable associated with the axis must be included through the `expression` field using the `objects` dictionary. Output histogram's layout is defined with the `layout` field. In the example above, our output dictionary will contain two histograms labelled `muon` and `zcandidate`, the first with the `muon_pt`, `muon_eta` and `muon_phi` axes, and the second only with the `dimuon_mass` axis (make sure to include axis with the same dimensions within a histogram). If you set `layout: individual` then the output dictionary will contain a histogram for each axis. Note that if you set `add_syst_axis: true`, a StrCategory axis `{"variable_name": {"type": "StrCategory", "categories": [], "growth": True}}` to store systematic variations will be added to each histogram.

**Every axis in `axes` is evaluated on every chunk**, whether or not a `layout` uses it: an axis whose
expression needs a missing object or selection breaks the run (with `extends`, remove an inherited
axis with `<axis>: null`). A Boolean axis must never receive `None`, and all the axes of a multi-axis
histogram must be per-event scalars (jagged arrays are flattened and the per-event weights are
broadcast from the last axis).

### `plotting` (used by `scripts/htautau/*` only)

```yaml
plotting:
  category_labels:                 # <category>: output folder name under outputs/plots/<workflow>/<year>/
    gen_etau: 1_gen_full_region
    sel_etau_trig: 4_selection_trigger
  group_labels:                    # folder name of a variable group
    lepton: electron
  groups:                          # variable groups = fnmatch patterns on the axis names
    lepton: ["ele_*"]
```

### Expression scopes

| where | names available |
|---|---|
| `object_selection` (`cuts`, `add_cut`, `add_field`, a `field` containing `events`) | `events`, `objects`, `year`, `np`, `ak`, `working_points`, the helpers of `analysis/selections/__init__.py` (`delta_r_higher`, `delta_r_lower`, ...), `candidate` |
| `event_selection.selections`, `histogram_config.axes[*].expression` | `events`, `objects`, `year`, `dataset`, `hlt_paths`, `selection_masks` (the selections defined before), `category_masks` (axes only), `np`, `ak`, and the mask helpers `get_trigger_mask`, `get_trigger_flag_mask(events, flag, year)`, `get_trigger_flags_or(events, [flags], year)`, `get_trigger_region_mask(events, hlt_paths, dataset, year, [flags])`, `get_lumi_mask(events, year)`, `get_metfilters_mask(events, year)`, `get_trigger_match_mask`, ... — **not** `working_points` |

