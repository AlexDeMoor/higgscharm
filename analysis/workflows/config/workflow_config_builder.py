import copy
import fnmatch
import yaml
from pathlib import Path
from analysis.histograms import HistogramConfig
from .workflow_config import WorkflowConfig


WORKFLOWS_DIR = Path(__file__).resolve().parents[1]


def _load_yaml(name: str) -> dict:
    """Load `analysis/workflows/<name>.yaml` (`name` may contain a sub-directory,
    e.g. 'includes/htautau_common')."""
    path = WORKFLOWS_DIR / f"{name}.yaml"
    if not path.exists():
        raise FileNotFoundError(f"workflow config {path} does not exist")
    with open(path, "r") as f:
        return yaml.safe_load(f) or {}


def deep_merge(base: dict, child: dict) -> dict:
    """Recursively merge `child` into a copy of `base`.

    - keys present in both with dict values are merged recursively;
    - any other key in `child` replaces the base value (lists are replaced, not
      concatenated);
    - a `null` value in `child` deletes the key (kept as a None tombstone while an
      `extends` chain is being resolved, stripped by strip_tombstones at the end);
    - base key order is preserved, new keys are appended (this matters for
      `object_selection`, whose objects are built sequentially).
    """
    out = copy.deepcopy(base)
    for key, value in child.items():
        if value is None:
            out[key] = None
        elif isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def strip_tombstones(config):
    """remove the None-valued keys (deletions) recursively"""
    if isinstance(config, dict):
        return {k: strip_tombstones(v) for k, v in config.items() if v is not None}
    return config


class WorkflowConfigBuilder:
    def __init__(self, workflow: str):
        config = _load_yaml(workflow)
        # `extends: <name>` chains: the named yaml is loaded (recursively) and this
        # config is merged on top of it (see deep_merge for the rules)
        seen = [workflow]
        while "extends" in config:
            base_name = config.pop("extends")
            if base_name in seen:
                raise ValueError(f"circular 'extends' chain: {seen + [base_name]}")
            seen.append(base_name)
            base = _load_yaml(base_name)
            base_extends = base.pop("extends", None)
            config = deep_merge(base, config)
            if base_extends is not None:
                config["extends"] = base_extends
        self.config = strip_tombstones(config)
        if not (self.config.get("event_selection") or {}).get("categories"):
            raise ValueError(f"workflow {workflow!r} defines no event_selection.categories")

    def build_workflow_config(self):
        return WorkflowConfig(
            object_selection=self.parse_object_selection(),
            event_selection=self.parse_event_selection(),
            corrections_config=self.parse_corrections_config(),
            histogram_config=self.parse_histogram_config(),
            datasets=self.parse_datasets_config(),
            mva=self.config.get("mva"),
            inference=self.config.get("inference"),
            combine=self.config.get("combine"),
            plotting=self.config.get("plotting"),
            plugins=self.config.get("plugins"),
        )

    def parse_object_selection(self):
        object_selection = {}
        for object_name in self.config["object_selection"]:
            object_selection[object_name] = {
                "field": self.config["object_selection"][object_name]["field"]
            }
            if "cuts" in self.config["object_selection"][object_name]:
                object_selection[object_name]["cuts"] = self.config["object_selection"][
                    object_name
                ]["cuts"]
            if "add_cut" in self.config["object_selection"][object_name]:
                cuts_to_add = self.config["object_selection"][object_name]["add_cut"]
                object_selection[object_name]["add_cut"] = {}
                for cut_name, cuts in cuts_to_add.items():
                    object_selection[object_name]["add_cut"][cut_name] = cuts
            if "add_field" in self.config["object_selection"][object_name]:
                fields_to_add = self.config["object_selection"][object_name][
                    "add_field"
                ]
                object_selection[object_name]["add_field"] = {}
                for field_name, fields in fields_to_add.items():
                    object_selection[object_name]["add_field"][field_name] = fields
        return object_selection

    def parse_event_selection(self):
        event_selection = {}
        for cut_name, cut in self.config["event_selection"].items():
            event_selection[cut_name] = cut
        if "mc_only_categories" in event_selection:
            event_selection["mc_only_categories"] = self.resolve_mc_only_categories(
                event_selection["mc_only_categories"], event_selection["categories"]
            )
        return event_selection

    @staticmethod
    def resolve_mc_only_categories(entries, categories):
        """`event_selection.mc_only_categories`: category names or fnmatch patterns of the
        categories that are not filled for data (BaseProcessor skips them: no histograms,
        parquets or cutflow). Returns the list of matched category names (yaml order); an
        entry that matches no category is an error (typo / renamed category)."""
        if entries is None:
            return []
        if isinstance(entries, str) or not isinstance(entries, list):
            raise ValueError(
                f"event_selection.mc_only_categories must be a list of category names or "
                f"fnmatch patterns, got {entries!r}"
            )
        unmatched = [
            entry
            for entry in entries
            if not any(fnmatch.fnmatchcase(category, str(entry)) for category in categories)
        ]
        if unmatched:
            raise ValueError(
                f"event_selection.mc_only_categories {unmatched} match no category "
                f"(categories: {list(categories)})"
            )
        return [
            category
            for category in categories
            if any(fnmatch.fnmatchcase(category, str(entry)) for entry in entries)
        ]

    def parse_histogram_config(self):
        hist_config = HistogramConfig(**copy.deepcopy(self.config["histogram_config"]))
        hist_config.categories = list(self.parse_event_selection()["categories"].keys())
        return hist_config

    def parse_corrections_config(self):
        corrections = {}
        corrections["object"] = self.config["corrections"]["object"]
        corrections["object_shifts"] = self.config["corrections"]["object_shifts"]
        # working points used by the tau energy scale (must match the tau selection)
        corrections["tau_es"] = self.config["corrections"].get("tau_es", {}) or {}
        corrections["event_weights"] = {}
        for name, vals in self.config["corrections"]["event_weights"].items():
            if isinstance(vals, bool):
                corrections["event_weights"][name] = vals
            elif name == "object_sf_categories":
                corrections["event_weights"][name] = vals
            elif isinstance(vals, list):
                corrections["event_weights"][name] = {}
                for val in vals:
                    for corr, wp in val.items():
                        corrections["event_weights"][name][corr] = wp
        return corrections

    def parse_datasets_config(self):
        return {key: self.config["datasets"][key] for key in self.config["datasets"]}
