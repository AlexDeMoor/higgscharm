"""Trigger study for one H->tautau channel workflow.

Outputs (outputs/plots/<workflow>/<year>/trigger_study/):
    trigger_eff.md, trigger_eff.json                         integrated efficiencies (+ errors), exclusive gains, acceptance
    <category folder>/<turn-on axis>_all_paths.png           every path + union (flat OR) + gated trigger, cH vs bH
    <category folder>/<turn-on axis>_union.png               flat OR, gated trigger and the trigger regions
    <category folder>/<turn-on axis>_acceptance.png          union inside / outside the reco acceptance (gen category only)

Boolean axes of the `trig_*` layouts:
  * `pass_<flag>`      one per HLT path (adopted paths + bookkeeping paths),
  * `pass_union`       pure OR of the adopted paths (what quick_trigg_study measured),
  * `pass_trigger`     the gated `trigger` selection: HLT bit AND offline plateau cuts on the
                       selected pair AND TrigObj matching (PD exclusive on data),
  * `pass_trig_<r>`    the individual gated trigger regions.
Per category the integrated efficiency of each of them, the exclusive gain of each adopted
path (union minus the union without it) and, in the gen-level category, the acceptance
block (efficiency inside / outside the reco acceptance, triggered events lost by the
selection) are tabulated ONCE (they do not depend on the turn-on axis); turn-on curves are
then drawn per layout. Categories that apply a trigger cut are skipped (their pass fractions
are P(path | trigger)). Weighted efficiencies with weighted-binomial errors (negative
generator weights!); bins with n_eff < 10 or eff outside [0, 1] are suppressed.

    python3 scripts/htautau/plot_trigger_study.py -w htautau_mutau -y 2024
"""
import json
import sys
import logging
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repo root
from scripts.htautau.common import (  # noqa: E402
    load_inputs, nominal, efficiency, cms_label, savefig, setup_style, plots_root, category_folder,
    category_is_triggered, PROCESS_LABELS, CHANNEL_LABELS, channel_of,
)

MIN_NEFF = 10.0
PROCESS_STYLE = {"H+c": dict(linestyle="-", marker="o"), "H+b": dict(linestyle="--", marker="s")}
PALETTE = plt.cm.tab10.colors + plt.cm.Dark2.colors


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-w", "--workflow", required=True, help="htautau_<channel>[_variant]")
    p.add_argument("-y", "--year", default="2024")
    p.add_argument("--processes", nargs="*", default=None, help="processes to include (default: all but Data)")
    p.add_argument("--extension", nargs="*", default=["png"])
    return p.parse_args()


def classify_axes(hist):
    """(paths, combined, regions) Boolean axes of a trig_* histogram"""
    names = list(hist.axes.name)
    paths = [n for n in names if n.startswith("pass_") and n not in ("pass_union", "pass_trigger") and not n.startswith("pass_trig_")]
    combined = [n for n in ("pass_union", "pass_trigger") if n in names]
    regions = [n for n in names if n.startswith("pass_trig_")]
    return paths, combined, regions


def bool_efficiency(hist, axis):
    """integrated efficiency (+ error, n_eff) of one Boolean axis"""
    total = hist.sum(flow=True)
    num = hist.project(axis)[True]
    eff, err, neff = efficiency([num.value], [num.variance], [total.value], [total.variance], min_neff=0.0)
    return float(eff[0]), float(err[0]), float(neff[0])


def or_efficiency(h_joint, flags, all_flags):
    """P(OR of `flags`) from the joint projection over `all_flags`: 1 - P(all False)"""
    total = h_joint.sum(flow=True).value
    if total <= 0 or not flags:
        return np.nan if total <= 0 else 0.0
    sel = {f: 0 for f in flags}
    for f in all_flags:
        if f not in flags:
            sel[f] = sum
    none = h_joint[sel]
    none_val = none.sum(flow=True).value if hasattr(none, "sum") else float(none.value)
    return 1.0 - none_val / total


def integrated_table(hist, adopted):
    """per-path / union / trigger / region efficiencies and exclusive gains of one category"""
    paths, combined, regions = classify_axes(hist)
    total = hist.sum(flow=True)
    neff = total.value**2 / total.variance if total.variance > 0 else 0.0
    rows = {}
    for b in paths + combined + regions:
        eff, err, _ = bool_efficiency(hist, b)
        rows[b] = {"eff": eff, "err": err, "kind": "path" if b in paths else "combined" if b in combined else "region",
                   "adopted": b in adopted}
    joint = hist.project(*paths)
    adopted_present = [f for f in adopted if f in paths]
    union_from_flags = or_efficiency(joint, adopted_present, paths)
    for f in paths:
        if f in adopted_present:
            others = [g for g in adopted_present if g != f]
            rows[f]["exclusive_gain"] = float(union_from_flags - or_efficiency(joint, others, paths))
        else:
            rows[f]["exclusive_gain"] = float(or_efficiency(joint, adopted_present + [f], paths) - union_from_flags)
    return {"sumw": float(total.value), "n_eff": float(neff), "rows": rows, "union_from_adopted_flags": float(union_from_flags)}


def acceptance_block(hist):
    """gen category: union efficiency inside / outside the reco acceptance"""
    out = {}
    tot = hist.sum(flow=True).value
    if "in_reco_acc" in hist.axes.name and tot > 0:
        j = hist.project("pass_union", "in_reco_acc")
        n_in, n_out = j[:, True].sum(flow=True).value, j[:, False].sum(flow=True).value
        out["frac_in_reco_acc"] = float(n_in / tot)
        out["eff_union_in_acc"] = float(j[True, True].value / n_in) if n_in else np.nan
        out["eff_union_outside_acc"] = float(j[True, False].value / n_out) if n_out else np.nan
        n_trig = j[True, :].sum(flow=True).value
        out["frac_triggered_lost_by_selection"] = float(j[True, False].value / n_trig) if n_trig else np.nan
        if "pass_trigger" in hist.axes.name:
            jt = hist.project("pass_trigger", "in_reco_acc")
            out["eff_trigger_in_acc"] = float(jt[True, True].value / n_in) if n_in else np.nan
    if "in_gen_acc" in hist.axes.name and tot > 0:
        jg = hist.project("pass_union", "in_gen_acc")
        n_gen = jg[:, True].sum(flow=True).value
        out["frac_in_gen_acc"] = float(n_gen / tot)
        out["eff_union_in_gen_acc"] = float(jg[True, True].value / n_gen) if n_gen else np.nan
    return out


def turnon_curves(hist, xname, flags, extra_sel=None):
    """{flag: (centers, eff, err, edges)} vs the turn-on axis; all-NaN series are dropped"""
    h = hist[extra_sel] if extra_sel else hist
    den = h.project(xname)
    edges = den.axes[0].edges
    centers = 0.5 * (edges[1:] + edges[:-1])
    out = {}
    for f in flags:
        num = h.project(xname, f)[:, True]
        eff, err, _ = efficiency(num.values(), num.variances(), den.values(), den.variances(), MIN_NEFF)
        if np.all(np.isnan(eff)):
            continue
        out[f] = (centers, eff, err, edges)
    return out


def plot_turnons(curves_by_process, xlabel, title, path, lumi, extensions, flags, ylabel="Trigger efficiency"):
    """colour = path / region, marker + linestyle = process"""
    setup_style()
    fig, ax = plt.subplots(figsize=(9, 7))
    colors = {f: PALETTE[i % len(PALETTE)] for i, f in enumerate(flags)}
    edges = None
    for process, curves in curves_by_process.items():
        pstyle = PROCESS_STYLE.get(process, dict(linestyle=":", marker="^"))
        for flag, (centers, eff, err, ed) in curves.items():
            if flag not in flags:
                continue
            edges = ed
            emph = flag in ("pass_union", "pass_trigger")
            ax.errorbar(centers, eff, yerr=err, xerr=(ed[1:] - ed[:-1]) / 2, capsize=2, markersize=5,
                        linewidth=1.8 if emph else 1.0, color=colors[flag],
                        label=f"{PROCESS_LABELS.get(process, process)}: {flag.replace('pass_', '')}", **pstyle)
    if edges is None:
        plt.close(fig)
        return
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_ylim(0, 1.15); ax.set_xlim(edges[0], edges[-1]); ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8, ncol=2 if len(flags) > 3 else 1, title=title, title_fontsize=10)
    cms_label(ax, lumi)
    savefig(fig, path, extensions)
    plt.close(fig)


def fmt(row):
    return f"{row['eff']:.4f} ± {row['err']:.4f}" if np.isfinite(row["eff"]) else "n/a"


def write_tables(summary, config, study_dir, adopted):
    with open(study_dir / "trigger_eff.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    L = [f"# Trigger study: {summary['workflow']} ({summary['year']})", "",
         f"Adopted paths (union = pure OR): {', '.join(a.replace('pass_', '') for a in adopted)}",
         "`trigger` = OR of the gated regions (HLT bit + offline plateau cuts on the selected pair + TrigObj matching). "
         "Efficiencies are weighted (negative generator weights) with weighted-binomial errors; "
         "exclusive gain = union minus the union without the path (for non-adopted paths: the gain it would add).", ""]
    for category, entry in summary["categories"].items():
        per_process = entry["integrated"]
        if not per_process:
            continue
        procs = list(per_process)
        L += [f"## {category_folder(config, category)}  (category `{category}`)", "",
              "| | " + " | ".join(f"eff {p} | excl. gain {p}" for p in procs) + " |", "|---|" + "---|---|" * len(procs)]
        rows = next(iter(per_process.values()))["rows"]
        for kind, header in (("path", "paths"), ("combined", "combined"), ("region", "gated regions")):
            names = [n for n, r in rows.items() if r["kind"] == kind]
            if not names:
                continue
            L.append(f"| *{header}* | " + " | ".join(" | " for _ in procs) + " |")
            for n in names:
                cells = []
                for p in procs:
                    r = per_process[p]["rows"][n]
                    gain = f"{r['exclusive_gain']:+.4f}" if "exclusive_gain" in r else ""
                    tag = "" if r["adopted"] or kind != "path" else " (not adopted)"
                    cells.append(f"{fmt(r)} | {gain}{tag}")
                label = f"**{n.replace('pass_', '')}**" if kind == "combined" else n.replace("pass_", "")
                L.append(f"| {label} | " + " | ".join(cells) + " |")
        L.append("| sumw (n_eff) | " + " | ".join(f"{per_process[p]['sumw']:.1f} ({per_process[p]['n_eff']:.0f}) | " for p in procs) + " |")
        acc = entry.get("acceptance") or {}
        if acc:
            L += ["", "acceptance (gen-level category; union = pure OR, trigger = gated selection):", "",
                  "| quantity | " + " | ".join(acc) + " |", "|---|" + "---|" * len(acc)]
            keys = ["frac_in_gen_acc", "frac_in_reco_acc", "eff_union_in_gen_acc", "eff_union_in_acc", "eff_trigger_in_acc",
                    "eff_union_outside_acc", "frac_triggered_lost_by_selection"]
            for k in keys:
                L.append(f"| {k} | " + " | ".join(f"{acc[p].get(k, float('nan')):.4f}" for p in acc) + " |")
        L.append("")
    (study_dir / "trigger_eff.md").write_text("\n".join(L))
    logging.info(f"tables written to {study_dir}/trigger_eff.{{json,md}}")


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    config, histograms, output_dir, lumi, colors = load_inputs(args.workflow, args.year)
    channel = channel_of(args.workflow)
    study_dir = plots_root(args.workflow, args.year) / "trigger_study"
    study_dir.mkdir(parents=True, exist_ok=True)
    processes = args.processes or [p for p in histograms if p != "Data"]
    hlt_paths = config.event_selection["hlt_paths"]
    adopted = [f"pass_{flag}" for flags in hlt_paths.values() for flag in flags]
    layouts = [name for name in config.histogram_config.layout if name.startswith("trig_")]
    categories = [c for c in config.histogram_config.categories if not category_is_triggered(config, c)]
    skipped = [c for c in config.histogram_config.categories if category_is_triggered(config, c)]
    if skipped:
        logging.info(f"trigger-applied categories skipped (pass fractions are conditional): {skipped}")
    labels = {name: config.histogram_config.axes[name].label for name in config.histogram_config.axes}

    # ---- 1. tables (once per category, from the first layout that has events) ----
    summary = {"workflow": args.workflow, "year": args.year, "adopted_paths": adopted, "categories": {}}
    for category in categories:
        entry = {"integrated": {}, "acceptance": {}}
        for process in processes:
            for layout in layouts:
                if layout not in histograms[process]:
                    continue
                h = nominal(histograms[process][layout], category)
                if h.sum(flow=True).value == 0:
                    continue
                entry["integrated"][process] = integrated_table(h, adopted)
                break
            if category.startswith("gen_"):
                for layout in layouts:
                    if layout not in histograms[process]:
                        continue
                    h = nominal(histograms[process][layout], category)
                    if "in_reco_acc" in h.axes.name and h.sum(flow=True).value > 0:
                        entry["acceptance"][process] = acceptance_block(h)
                        break
        summary["categories"][category] = entry
    write_tables(summary, config, study_dir, adopted)

    # ---- 2. turn-on plots ----
    for category in categories:
        cat_dir = study_dir / category_folder(config, category)
        title = f"{CHANNEL_LABELS[channel]}  [{category_folder(config, category)}]"
        for layout in layouts:
            if layout.startswith("trig_gen") and not category.startswith("gen_"):
                continue    # gen-level turn-on axes / acceptance flags: only meaningful in the gen category
            curves, acc_curves, xname, flag_sets = {}, {}, None, None
            for process in processes:
                if layout not in histograms[process]:
                    continue
                h = nominal(histograms[process][layout], category)
                if h.sum(flow=True).value == 0:
                    continue
                xname = h.axes.name[0]
                paths, combined, regions = classify_axes(h)
                flag_sets = (paths, combined, regions)
                c = turnon_curves(h, xname, paths + combined + regions)
                if c:
                    curves[process] = c
                if "in_reco_acc" in h.axes.name and category.startswith("gen_"):
                    inside = turnon_curves(h, xname, ["pass_union", "pass_trigger"], {"in_reco_acc": True})
                    outside = turnon_curves(h, xname, ["pass_union"], {"in_reco_acc": False})
                    acc = {}
                    if "pass_union" in inside:
                        acc["pass_union_in_acc"] = inside["pass_union"]
                    if "pass_trigger" in inside:
                        acc["pass_trigger_in_acc"] = inside["pass_trigger"]
                    if "pass_union" in outside:
                        acc["pass_union_outside_acc"] = outside["pass_union"]
                    if "pass_union" in curves.get(process, {}):
                        acc["pass_union_all"] = curves[process]["pass_union"]
                    if acc:
                        acc_curves[process] = acc
            if not curves:
                continue
            paths, combined, regions = flag_sets
            xlabel = labels.get(xname, xname)
            plot_turnons(curves, xlabel, title, cat_dir / f"{xname}_all_paths", lumi, args.extension, paths + combined)
            plot_turnons(curves, xlabel, title, cat_dir / f"{xname}_union", lumi, args.extension, combined + regions)
            if acc_curves:
                plot_turnons(acc_curves, xlabel, title, cat_dir / f"{xname}_acceptance", lumi, args.extension,
                             ["pass_union_all", "pass_union_in_acc", "pass_trigger_in_acc", "pass_union_outside_acc"],
                             ylabel="Trigger efficiency (union / gated)")
    logging.info(f"plots written under {study_dir}")


if __name__ == "__main__":
    main()
