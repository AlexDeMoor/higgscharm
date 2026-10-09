"""Shared helpers for the H->tautau plotting scripts (trigger study, shapes, index).

Human-facing products live under outputs/plots/<workflow>/<year>/ (the framework's own
outputs/<workflow>/<year>/ tree is left untouched). Folder names and variable groups come
from the workflow's `plotting:` block (category_labels, group_labels, groups).
"""
import fnmatch
import logging
from pathlib import Path

import numpy as np
import yaml
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mplhep as hep  # noqa: E402
from coffea.util import load  # noqa: E402

from analysis.workflows.config import WorkflowConfigBuilder  # noqa: E402

OUTPUT_DIR = Path.cwd() / "outputs"
PLOTS_ROOT = OUTPUT_DIR / "plots"
POSTPROCESS_DIR = Path.cwd() / "analysis" / "postprocess"

PROCESS_LABELS = {"H+c": r"cH, H$\rightarrow\tau\tau$", "H+b": r"bH, H$\rightarrow\tau\tau$"}
CHANNEL_LABELS = {
    "mutau": r"$\mu\tau_h$",
    "etau": r"$e\tau_h$",
    "tautau": r"$\tau_h\tau_h$",
}
CHANNEL_TEXT = {"mutau": "mu tau_h", "etau": "e tau_h", "tautau": "tau_h tau_h"}


def channel_of(workflow: str) -> str:
    """'htautau_mutau', 'htautau_mutau_test', ... -> 'mutau'"""
    name = workflow.replace("htautau_", "")
    for channel in CHANNEL_LABELS:
        if name == channel or name.startswith(channel + "_"):
            return channel
    raise ValueError(f"cannot infer the channel from workflow name {workflow!r}")


def setup_style():
    hep.style.use(hep.style.CMS)
    plt.rcParams.update({"font.size": 14, "legend.fontsize": 11, "axes.labelsize": 16})


def plots_root(workflow, year) -> Path:
    return PLOTS_ROOT / workflow / year


def load_inputs(workflow, year):
    """workflow config + the per-process histograms produced by run_postprocess.py"""
    config = WorkflowConfigBuilder(workflow).build_workflow_config()
    output_dir = OUTPUT_DIR / workflow / year
    hist_file = output_dir / f"{year}_processed_histograms.coffea"
    if not hist_file.exists():
        raise FileNotFoundError(
            f"{hist_file} not found: run 'python3 run_postprocess.py -w {workflow} -y {year} "
            f"--postprocess --output_format coffea' first"
        )
    histograms = load(hist_file)
    with open(POSTPROCESS_DIR / "luminosity.yaml") as f:
        lumi = yaml.safe_load(f)[year]
    with open(POSTPROCESS_DIR / "process_color_map.yaml") as f:
        colors = yaml.safe_load(f)
    return config, histograms, output_dir, lumi, colors


# ---------------------------------------------------------------------------------
# folder / group resolution from the `plotting:` block
# ---------------------------------------------------------------------------------
def category_folder(config, category: str) -> str:
    return (config.plotting.get("category_labels") or {}).get(category, category)


def group_of(config, variable: str) -> str:
    """group (folder) of an axis / layout name: first matching fnmatch pattern; 'other' otherwise"""
    groups = config.plotting.get("groups") or {}
    labels = config.plotting.get("group_labels") or {}
    for group, patterns in groups.items():
        if any(fnmatch.fnmatch(variable, pattern) for pattern in patterns or []):
            return labels.get(group, group)
    return "other"


def is_2d_layout(layout, variables) -> bool:
    """two-axis layout drawn as a 2D map / per-category curves (the trig_* and *_by_tauwp
    multi-axis layouts have their own code paths)"""
    return len(variables) == 2 and not layout.startswith("trig_") and not layout.endswith("_by_tauwp")


def ordered_groups(config):
    groups = config.plotting.get("groups") or {}
    labels = config.plotting.get("group_labels") or {}
    return [labels.get(g, g) for g in groups] + ["other"]


def nominal(hist, category):
    """slice the nominal variation of one category (keeps the remaining axes)"""
    selector = {}
    if "variation" in hist.axes.name:
        selector["variation"] = "nominal"
    if "category" in hist.axes.name:
        selector["category"] = category
    return hist[selector] if selector else hist


def efficiency(num_values, num_variances, den_values, den_variances, min_neff=10.0):
    """weighted efficiency eps = sum(w_pass) / sum(w_all) with the weighted-binomial error

        var(eps) = [ (1 - eps)^2 * sum(w_pass^2) + eps^2 * sum(w_fail^2) ] / sum(w_all)^2

    (num_variances = sum(w_pass^2), den_variances = sum(w_all^2) as stored by hist's Weight
    storage). With negative generator weights eps can leave [0, 1] in poorly populated bins:
    such bins, bins with n_eff = sum(w)^2 / sum(w^2) < min_neff and empty denominators are
    returned as NaN. Returns (eff, err, n_eff), all arrays."""
    num_values = np.asarray(num_values, dtype=float)
    num_variances = np.asarray(num_variances, dtype=float)
    den_values = np.asarray(den_values, dtype=float)
    den_variances = np.asarray(den_variances, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        eff = num_values / den_values
        neff = den_values**2 / den_variances
        fail_var = np.clip(den_variances - num_variances, 0, None)
        err = np.sqrt((1 - eff) ** 2 * num_variances + eff**2 * fail_var) / np.abs(den_values)
    bad = ~np.isfinite(eff) | (neff < min_neff) | (den_values <= 0) | (eff < 0) | (eff > 1)
    eff = np.where(bad, np.nan, eff)
    err = np.where(bad, np.nan, err)
    return eff, err, neff


def category_is_triggered(config, category):
    """True when the category applies a trigger requirement (its pass fractions are then
    P(path | trigger), not efficiencies)"""
    cuts = config.event_selection["categories"].get(category, [])
    return any(c == "trigger" or c.startswith("trig_") for c in cuts)


def cms_label(ax, lumi, extra="Preliminary"):
    # data=False already prepends "Simulation" to the label
    hep.cms.label(extra, data=False, lumi=f"{lumi * 1e-3:.1f}", year="2024", com="13.6", ax=ax, fontsize=13)


def savefig(fig, path, extensions=("png",)):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    for ext in extensions:
        fig.savefig(path.with_suffix(f".{ext}"), bbox_inches="tight", dpi=120)
    plt.close(fig)
    logging.info(f"saved {path.with_suffix('.' + extensions[0])}")
