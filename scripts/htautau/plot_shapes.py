"""cH vs bH shape comparison for one H->tautau channel workflow.

Writes, for every category (folder from `plotting.category_labels`) and every variable
(folder from `plotting.groups`):

    outputs/plots/<workflow>/<year>/<category folder>/<group>/<variable>.png

with two panels: lumi-scaled yields (left) and unit-area shapes (right). The `*_by_tauwp`
layouts (variable split by the tau_h working-point Boolean axes) go to `tau_wp_study/`. Other
two-axis layouts (ctag2d_*, etaphi_*, *_response_by_dm) are drawn one panel per process, under
`<group of the layout name>/<layout>.png`: a 2D map of the yields (--logz: log colour scale), or,
when the second axis is an IntCategory, the unit-area shape per category (+ inclusive). The
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
    PROCESS_LABELS, CHANNEL_LABELS, channel_of, is_2d_layout,
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
    p.add_argument("--logz", action="store_true", help="log colour scale of the 2D maps")
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


# ---------------------------------------------------------------------------------
# two-axis layouts (not trig_* / *_by_tauwp): (x, IntCategory) -> one curve per category,
# (x, y) -> a 2D map; one panel per process
# ---------------------------------------------------------------------------------
def _is_category(axis):
    return type(axis).__name__ in ("IntCategory", "StrCategory")


def _edges_ticks(axis):
    """bin edges of an axis + (tick positions, labels) for categorical axes (None otherwise)"""
    if _is_category(axis):
        n = len(axis)
        return np.arange(n + 1, dtype=float), (np.arange(n) + 0.5, [str(c) for c in axis])
    return np.asarray(axis.edges, dtype=float), None


def weighted_median(edges, values):
    values = np.clip(np.asarray(values, dtype=float), 0, None)   # negative genWeights
    cum = np.cumsum(values)
    if cum.size == 0 or cum[-1] <= 0:
        return np.nan
    i = min(np.searchsorted(cum, cum[-1] / 2), len(values) - 1)
    return 0.5 * (edges[i] + edges[i + 1])


def draw_2d_layout(hists, hc, variables, title, legend_title, lumi, logz=False, log=False):
    """one panel per process. (x, IntCategory y): unit-area x shape of every populated y category
    plus the inclusive one (weighted median in the legend). (x, y): map of the yields, optional
    log colour scale (bins <= 0 left blank; linear when nothing is positive)."""
    from matplotlib.colors import LogNorm
    xname, yname = variables
    fig, axes = plt.subplots(1, len(hists), figsize=(8.5 * len(hists), 6.5), squeeze=False)
    for ax, (process, h) in zip(axes[0], hists.items()):
        h2 = h.project(xname, yname)
        label = PROCESS_LABELS.get(process, process)
        yaxis = h2.axes[yname]
        if _is_category(yaxis):
            xedges = np.asarray(h2.axes[xname].edges, dtype=float)
            values = h2.values()   # (x, y) without flow
            curves = [(f"all {hc.axes[yname].label}", values.sum(axis=1), "k")]
            curves += [(f"{hc.axes[yname].label} = {cat}", values[:, i], plt.cm.tab10(i % 10))
                       for i, cat in enumerate(yaxis)]
            for name, v, color in curves:
                total = v.sum()
                if total <= 0:
                    continue   # empty / net-negative category
                ax.stairs(v / total, xedges, color=color, linewidth=1.8,
                          label=f"{name} (median {weighted_median(xedges, v):.3f})")
            ax.set_xlabel(hc.axes[xname].label); ax.set_ylabel("a.u.")
            if log:
                ax.set_yscale("log")
            ax.grid(alpha=0.3)
            if ax.get_legend_handles_labels()[0]:
                ax.legend(fontsize=8, title=f"{title} {label} {legend_title}", title_fontsize=9)
            else:
                ax.set_title(f"{title} {label} {legend_title}: empty", fontsize=11, pad=28)
        else:
            xedges, xticks = _edges_ticks(h2.axes[xname])
            yedges, yticks = _edges_ticks(yaxis)
            values = h2.values().T   # (y, x) for pcolormesh
            norm, shown = None, values
            positive = values > 0
            if logz and positive.any():
                vmin, vmax = values[positive].min(), values[positive].max()
                norm = LogNorm(vmin=vmin, vmax=max(vmax, vmin * 1.0001))
                shown = np.where(positive, values, np.nan)
            mesh = ax.pcolormesh(xedges, yedges, shown, norm=norm, cmap="viridis", shading="flat")
            fig.colorbar(mesh, ax=ax, label=f"Events / bin ({lumi * 1e-3:.1f} fb$^{{-1}}$)")
            for ticks, setter in ((xticks, (ax.set_xticks, ax.set_xticklabels)), (yticks, (ax.set_yticks, ax.set_yticklabels))):
                if ticks is not None:
                    setter[0](ticks[0]); setter[1](ticks[1])
            ax.set_xlabel(hc.axes[xname].label); ax.set_ylabel(hc.axes[yname].label)
            ax.set_title(f"{title} {label} ({values.sum():.2f} ev) {legend_title}", fontsize=11, pad=28)
        cms_label(ax, lumi)
    return fig


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config, histograms, output_dir, lumi, colors = load_inputs(args.workflow, args.year)
    hc = config.histogram_config
    root = plots_root(args.workflow, args.year)
    processes = args.processes or [p for p in histograms if p != "Data"]
    categories = args.categories or hc.categories
    title = CHANNEL_LABELS[channel_of(args.workflow)]
    # layouts the yaml defines but the processed outputs lack (outputs produced before an axis was
    # added): skipped with one warning instead of a KeyError -- rerun the processing to get them
    missing = sorted({layout for layout in hc.layout for process in processes if layout not in histograms[process]})
    if missing:
        print(f"WARNING: {len(missing)} layouts not in the processed histograms (older outputs, rerun "
              f"the processing to get them), skipped: {', '.join(missing)}")

    for category in categories:
        cat_dir = root / category_folder(config, category)
        # cutflow copy (framework CSV -> descriptive folder name)
        for prefix, suffix in (("cutflow", ""), ("cutflow_raw", "_raw")):
            src = output_dir / category / f"{prefix}_{category}.csv"
            if src.exists():
                (root / "cutflows").mkdir(parents=True, exist_ok=True)
                shutil.copy(src, root / "cutflows" / f"{category_folder(config, category)}{suffix}.csv")
        for layout, variables in hc.layout.items():
            if layout.startswith("trig_") or layout in missing:
                continue
            if args.variables and not any(v in args.variables for v in [*variables, layout]):
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
            if is_2d_layout(layout, variables):
                # saved under the LAYOUT name: the axes may be shared with 1D layouts
                # (etaphi_tau = [tau_eta, tau_phi]) whose plots must not be overwritten
                try:
                    fig = draw_2d_layout(hists, hc, variables, title, f"[{category_folder(config, category)}]",
                                         lumi, logz=args.logz, log=args.log)
                    savefig(fig, cat_dir / group_of(config, layout) / layout, args.extension)
                except Exception as exc:   # one bad layout must not stop the whole run
                    logging.warning(f"{category}/{layout}: not plotted ({exc})")
                finally:
                    plt.close("all")
                continue
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
