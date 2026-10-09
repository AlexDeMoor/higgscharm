import numpy as np
import awkward as ak
from copy import deepcopy
from coffea import processor
from coffea.nanoevents import NanoAODSchema
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods.vector import LorentzVector
from analysis.utils import dump_lumi, update, add_cutflow, dump_parquet
from analysis.utils.parquet_writer import dump_chunk_sumw
from analysis.workflows.config import WorkflowConfigBuilder
from analysis.histograms import HistBuilder, fill_histograms

from analysis.corrections.correction_manager import (
    object_corrector_manager,
    weight_manager,
)
from analysis.selections import (
    ObjectSelector,
    get_lumi_mask,
    get_trigger_mask,
    get_trigger_flag_mask,
    get_trigger_flags_or,
    get_trigger_region_mask,
    get_zzto4l_trigger_mask,
    get_metfilters_mask,
    get_trigger_match_mask,
    get_stitching_mask,
)

NanoAODSchema.warn_missing_crossrefs = False


class BaseProcessor(processor.ProcessorABC):
    def __init__(
        self,
        workflow: str,
        year: str,
        output_format: str,
        output_location: str,
    ):
        self.year = year
        self.workflow = workflow
        self.output_format = output_format
        self.output_location = output_location

        config_builder = WorkflowConfigBuilder(workflow)
        self.workflow_config = config_builder.build_workflow_config()
        self.workflow_config.load_plugins()
        self.histogram_config = self.workflow_config.histogram_config
        self.histograms = HistBuilder(self.workflow_config).build_histogram()

    def process(self, events):
        # again in the worker process (no-op when already imported, e.g. after a fork)
        self.workflow_config.load_plugins()
        self.is_mc = hasattr(events, "genWeight")
        # Full generator sumw of this chunk BEFORE any selection/veto (the jet veto map
        # removes events before process_shift, so it must be recorded here). Used for
        # the .coffea metadata; the parquet path additionally writes sumw_records.
        self._sumw = float(ak.sum(events.genWeight)) if self.is_mc else len(events)
        # raw event count of the chunk before the jet veto ("initial" cutflow row)
        self._nevents = len(events)
        if self.output_format == "parquet" and self.is_mc:
            dump_chunk_sumw(events, self.workflow, self.year, self.output_location)
        vetoed_events, shifts = object_corrector_manager(
            events=events,
            year=self.year,
            corrections_config=self.workflow_config.corrections_config,
        )
        return processor.accumulate(
            self.process_shift(update(vetoed_events, collections), shift)
            for collections, shift in shifts
        )

    def process_shift(self, events, shift):
        year = self.year
        dataset = events.metadata["dataset"]
        histograms = deepcopy(self.histograms)

        # initialize output dictionary to store histograms/arrays and metadata
        output = {}
        output["metadata"] = {}
        if shift is None:
            # add sum of genweights (before selection and before the jet veto) to metadata
            output["metadata"].update({"sumw": self._sumw})

        if not self.is_mc:
            # back-fill the MC-only object fields used by histogram expressions
            events["Jet", "hadronFlavour"] = ak.zeros_like(events.Jet.pt)
            events["Jet", "partonFlavour"] = ak.zeros_like(events.Jet.pt)
            for collection in ("Muon", "Electron", "Tau"):
                if collection in events.fields and "genPartFlav" not in events[collection].fields:
                    events[collection, "genPartFlav"] = ak.zeros_like(events[collection].pt)

        # --------------------------------------------------------------
        # Object selection
        # --------------------------------------------------------------
        object_selection = self.workflow_config.object_selection
        object_selector = ObjectSelector(object_selection, year)
        objects = object_selector.select_objects(events)

        # --------------------------------------------------------------
        # Event selection
        # --------------------------------------------------------------
        event_selection = self.workflow_config.event_selection
        hlt_paths = event_selection["hlt_paths"]
        if not self.is_mc and shift is None:
            # save (run, luminosityBlock) pairs to metadata (nominal only: the output of
            # every object shift is accumulated)
            lumi_mask = eval(event_selection["selections"]["lumimask"])
            dump_lumi(events[lumi_mask], output)

        selection_manager = PackedSelection()
        # named selection masks, filled incrementally: a selection may refer to the ones
        # defined BEFORE it as selection_masks['<name>'] (e.g. trigger = OR of regions);
        # the dict is also in scope for the histogram-axis expressions
        selection_masks = {}
        for selection, mask in event_selection["selections"].items():
            selection_manager.add(selection, eval(mask))
            selection_masks[selection] = selection_manager.all(selection)

        # categories that are not filled for data (event_selection.mc_only_categories:
        # gen-level, study and no-trigger categories, meaningless or double counted
        # across primary datasets on data): no histograms, parquets or cutflow
        skip_categories = (
            set() if self.is_mc else set(event_selection.get("mc_only_categories", []))
        )
        if shift is None:
            add_cutflow(
                events,
                output,
                selection_manager,
                self.workflow_config,
                initial=(self._sumw, self._nevents),
                skip_categories=skip_categories,
            )

        # per-category masks. Computed up-front so that histogram expressions can refer
        # to them through `category_masks[<category>]` (e.g. an "in acceptance" Boolean
        # axis in a gen-level category = passes the reco-level selection category).
        categories = event_selection["categories"]
        category_masks = {
            category: selection_manager.all(*category_cuts)
            for category, category_cuts in categories.items()
        }
        # evaluate every histogram-axis expression ONCE per chunk (they are the same
        # for all categories); each category then only applies its mask
        full_variables = {}
        nan_column = np.full(len(events), np.nan, dtype=np.float32)
        for variable, axis in self.histogram_config.axes.items():
            # explicit loop: a comprehension would hide the local scope from eval()
            if getattr(axis, "mc_only", False) and not self.is_mc:
                full_variables[variable] = nan_column
                continue
            full_variables[variable] = eval(axis.expression)
        # --------------------------------------------------------------
        # Histogram filling / array dumping
        # --------------------------------------------------------------
        for category, category_mask in category_masks.items():
            if category in skip_categories:
                continue
            nevents_after = ak.sum(category_mask)
            if nevents_after > 0:
                # get pruned events
                pruned_ev = events[category_mask]
                # add each selected object to 'pruned_ev' as a new field
                for obj in objects:
                    pruned_ev[f"selected_{obj}"] = objects[obj][category_mask]
                # get weights container
                weights_container = weight_manager(
                    pruned_ev=pruned_ev,
                    year=year,
                    dataset=dataset,
                    workflow_config=self.workflow_config,
                    category=category,
                    shift=shift,
                    selection_masks={
                        name: mask[category_mask] for name, mask in selection_masks.items()
                    },
                )
                # get analysis variables map
                variables_map = {
                    variable: array[category_mask]
                    for variable, array in full_variables.items()
                }
                # NanoAOD event id — used by scripts/mva/prep_training_inputs.py
                # for the deterministic train/test split
                variables_map["event"] = events.event[category_mask]

                if self.output_format == "coffea":
                    fill_histograms(
                        histogram_config=self.histogram_config,
                        weights_container=weights_container,
                        variables_map=variables_map,
                        histograms=histograms,
                        category=category,
                        is_mc=self.is_mc,
                        shift=shift,
                        flow=True,
                    )
                elif self.output_format == "parquet":
                    dump_parquet(
                        events=events,
                        weights_container=weights_container,
                        variables_map=variables_map,
                        workflow=self.workflow,
                        year=year,
                        category=category,
                        output_location=self.output_location,
                        shift=shift,
                    )

        # add histograms to output dictionary
        if self.output_format == "coffea":
            output["histograms"] = histograms
        return output

    def postprocess(self, accumulator):
        pass
