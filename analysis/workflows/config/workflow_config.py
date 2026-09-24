import yaml
import importlib


class WorkflowConfig:
    """
    Attributes:
    -----------
        object_selection:
        event_selection:
        corrections_config:
        histogram_config:
    """

    def __init__(
        self,
        object_selection,
        event_selection,
        corrections_config,
        histogram_config,
        datasets,
        mva=None,
        inference=None,
        combine=None,
        plotting=None,
        plugins=None,
    ):
        self.object_selection = object_selection
        self.event_selection = event_selection
        self.corrections_config = corrections_config
        self.histogram_config = histogram_config
        self.datasets = datasets
        self.mva = mva
        self.inference = inference
        self.combine = combine
        # plotting-only settings (output folder names / variable groups), see scripts/htautau
        self.plotting = plotting or {}
        # python modules imported before processing (see load_plugins)
        self.plugins = list(plugins or [])

    def load_plugins(self):
        """Import the modules listed under `plugins:` in the workflow yaml. A plugin module
        extends the framework at import time (e.g. adds ObjectSelector methods or registers
        trigger flags with analysis.selections.trigger.register_trigger_flags). Import errors
        are raised, never swallowed: a workflow must not run without its plugins."""
        for module in self.plugins:
            importlib.import_module(module)

    def to_dict(self):
        """Convert WorkflowConfig to a dictionary."""
        d = {
            "object_selection": self.object_selection,
            "event_selection": self.event_selection,
            "corrections_config": self.corrections_config,
            "histogram_config": self.histogram_config.to_dict(),
            "datasets": self.datasets,
        }
        for k in ("mva", "inference", "combine", "plotting"):
            if getattr(self, k) is not None:
                d[k] = getattr(self, k)
        if self.plugins:
            d["plugins"] = self.plugins
        return d

    def to_yaml(self):
        """Convert WorkflowConfig to a YAML string."""
        return yaml.dump(self.to_dict(), sort_keys=False, default_flow_style=False)
