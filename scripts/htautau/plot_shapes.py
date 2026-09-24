"""cH vs bH shape comparison for one H->tautau channel workflow.

Writes, for every category (folder from `plotting.category_labels`) and every variable
(folder from `plotting.groups`):

    outputs/plots/<workflow>/<year>/<category folder>/<group>/<variable>.png

with two panels: lumi-scaled yields (left) and unit-area shapes (right). The `*_by_tauwp`
layouts (variable split by the tau_h working-point Boolean axes) go to `tau_wp_study/`. The
framework cutflow CSVs are copied to `cutflows/<category folder>.csv`, and the README /
index.html are (re)generated at the end (make_index.py).

    python3 scripts/htautau/plot_shapes.py -w htautau_mutau -y 2024 [--categories sel_mutau ...]
"""
import sys
import shutil
import logging
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repo root
from scripts.htautau.common import (  # noqa: E402
    load_inputs, nominal, cms_label, savefig, setup_style, plots_root, category_folder, group_of,
    PROCESS_LABELS, CHANNEL_LABELS, channel_of,
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-w", "--workflow", required=True, help="htautau_<channel>[_variant]")
    p.add_argument("-y", "--year", default="2024")
    p.add_argument("--categories", nargs="*", default=None)
    p.add_argument("--variables", nargs="*", default=None)
    p.add_argument("--processes", nargs="*", default=None)
    p.add_argument("--extension", nargs="*", default=["png"])
    p.add_argument("--log", action="store_true")
    p.add_argument("--no-index", action="store_true", help="do not regenerate README.md / index.html")
    return p.parse_args()


def draw(ax, h, label, color, density, log):
    values = h.values()
    if density:
        integral = values.sum()
        h = h / integral if integral > 0 else h
    hep.histplot(h, ax=ax, label=label, color=color, histtype="step", linewidth=1.8, flow="none")
    if log:
        ax.set_yscale("log")


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config, histograms, output_dir, lumi, colors = load_inputs(args.workflow, args.year)
    hc = config.histogram_config
    root = plots_root(args.workflow, args.year)
    processes = args.processes or [p for p in histograms if p != "Data"]
    categories = args.categories or hc.categories
    title = CHANNEL_LABELS[channel_of(args.workflow)]

    for category in categories:
        cat_dir = root / category_folder(config, category)
        # cutflow copy (framework CSV -> descriptive folder name)
        for prefix, suffix in (("cutflow", ""), ("cutflow_raw", "_raw")):
            src = output_dir / category / f"{prefix}_{category}.csv"
            if src.exists():
                (root / "cutflows").mkdir(parents=True, exist_ok=True)
                shutil.copy(src, root / "cutflows" / f"{category_folder(config, category)}{suffix}.csv")
        for layout, variables in hc.layout.items():
            if layout.startswith("trig_"):
                continue
            if args.variables and not any(v in args.variables for v in variables):
                continue
            variable = variables[0]
            xlabel = hc.axes[variable].label
            hists = {}
            for process in processes:
                h = nominal(histograms[process][layout], category)
                if h.sum(flow=True).value > 0:
                    hists[process] = h
            if not hists:
                continue
            setup_style()
            if layout.endswith("_by_tauwp"):
                # variable split by each tau_h WP Boolean axis, one panel per process
                bool_axes = [n for n in variables[1:]]
                fig, axes = plt.subplots(1, len(hists), figsize=(8 * len(hists), 6.5), squeeze=False)
                for ax, (process, h) in zip(axes[0], hists.items()):
                    base = h.project(variable)
                    draw(ax, base, f"{PROCESS_LABELS.get(process, process)}: baseline", "k", False, args.log)
                    for i, bax in enumerate(bool_axes):
                        tight = h[{bax: True}].project(variable)
                        draw(ax, tight, f"+ {hc.axes[bax].label}", plt.cm.tab10(i), False, args.log)
                    ax.set_xlabel(xlabel); ax.set_ylabel("Events / bin"); ax.grid(alpha=0.3)
                    ax.legend(fontsize=8, title=f"{title} [{category_folder(config, category)}]", title_fontsize=9)
                    cms_label(ax, lumi)
                savefig(fig, cat_dir / (config.plotting.get("group_labels") or {}).get("tau_wp_study", "tau_wp_study") / layout, args.extension)
                plt.close(fig)
                continue
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6.5))
            try:
                for process, h in hists.items():
                    proj = h.project(variable)
                    label = PROCESS_LABELS.get(process, process)
                    draw(ax1, proj, f"{label} ({proj.sum().value:.2f} ev)", colors.get(process), False, args.log)
                    draw(ax2, proj, label, colors.get(process), True, args.log)
                ax1.set_ylabel(f"Events / bin ({lumi * 1e-3:.1f} fb$^{{-1}}$)"); ax2.set_ylabel("a.u.")
                for ax in (ax1, ax2):
                    ax.set_xlabel(xlabel); ax.grid(alpha=0.3)
                    ax.legend(fontsize=10, title=f"{title} [{category_folder(config, category)}]", title_fontsize=10)
                    cms_label(ax, lumi)
                savefig(fig, cat_dir / group_of(config, variable) / variable, args.extension)
            except Exception as exc:   # one bad axis must not stop the whole run
                logging.warning(f"{category}/{variable}: not plotted ({exc})")
            finally:
                plt.close(fig)

    if not args.no_index:
        from scripts.htautau.make_index import make_index
        make_index(args.workflow, args.year)


if __name__ == "__main__":
    main()
