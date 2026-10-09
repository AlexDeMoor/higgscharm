"""README.md + index.html for the plots of one H->tautau channel (and the top-level index).

    python3 scripts/htautau/make_index.py -w htautau_mutau -y 2024

Reads the workflow config (categories, selections, axes, `plotting:` block) and the files
present under outputs/plots/<workflow>/<year>/; safe to re-run at any time.
"""
import re
import sys
import csv
import html
import logging
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # repo root
from analysis.workflows.config import WorkflowConfigBuilder  # noqa: E402
from scripts.htautau.common import (  # noqa: E402
    PLOTS_ROOT, plots_root, category_folder, group_of, ordered_groups, channel_of, CHANNEL_TEXT,
    is_2d_layout,
)

DEFINITIONS = {
    "higgs_mass_vis": "visible mass of the di-tau (Higgs candidate) system m(l1 + l2)",
    "higgs_mass_col": "collinear-approximation di-tau mass (NaN/overflow when unphysical)",
    "higgs_mt_tot": "total transverse mass sqrt(mT(l1,l2)^2 + mT(l1,MET)^2 + mT(l2,MET)^2)",
    "higgs_dzeta": "D_zeta = p_zeta^miss - 0.85 p_zeta^vis (zeta = bisector of the two legs)",
    "higgs_met_pt": "pT of the di-tau + MET system",
    "higgs_mt_legs": "mT of the two visible legs",
    "mt_lepton_met": "mT(leg 1, MET); leg 1 = mu / e / leading tau_h",
    "mt_tau_met": "mT(tau_h leg, MET)",
    "tau_idx_vsjet": "DeepTau2018v2p5 VSjet working point passed (cumulative index 1..8)",
    "tau_idx_vse": "DeepTau2018v2p5 VSe working point passed (1..8)",
    "tau_idx_vsmu": "DeepTau2018v2p5 VSmu working point passed (1..4)",
    "gen_channel": "generator tau-tau decay: 1 mu tau_h, 2 e tau_h, 3 tau_h tau_h, 4 leptonic, 0 none",
    "met_pt": "PuppiMET (Type-1 JEC re-propagation w.r.t. the nano jets + tau energy scale propagated)",
    "higgs_mass_col_valid": "collinear-mass solution physical (0 < x1, x2 <= 1); m_col is NaN otherwise",
    "gen_higgs_mass_vis": "gen-level visible di-tau mass (hard-process tau decay products: lepton / visible tau_h)",
    "n_veto_muons": "loose muons (pT>10, looseId, PFIso04<0.3, no IP cut) used for the extra-lepton veto",
    "n_veto_electrons": "loose electrons (pT>10, mvaNoIso WP90, PFIso03<0.3, no IP cut) used for the extra-lepton veto",
    "pfmet_pt": "PF MET (uncorrected)",
    "bestcvl_pt": "the selected jet with the highest UParT CvL among ALL selected jets (c-tagged or not)",
    "hf_tag_category": "0: no c-tagged and no b-tagged jet, 1: >= 1 c-tagged and no b-tagged jet, 2: >= 1 b-tagged jet (UParT medium)",
    "cjet_dr_genc": "dR(c-jet candidate, nearest outgoing hard-process gen c quark); overflow when there is none (bH, light/gluon cH)",
    "truecjet_pt_rank": "pT rank (0 = highest) of the highest-pT true c-jet (hadronFlavour 4) among the selected jets; overflow = none",
    "truecjet_cvl_rank": "CvL rank (0 = highest) of the same jet; overflow = no true c-jet",
    "tau_response": "reco tau_h pT (after TES) / gen visible tau_h pT, genuine tau_h only (genPartFlav 5, GenVisTau via Tau_genPartIdx)",
    "tau_response_raw": "as tau_response with the pT before the tau energy scale",
    "rho": "Rho_fixedGridRhoFastjetAll",
    "pu_ntrueint": "Pileup_nTrueInt (true mean number of interactions, MC)",
}


TEX_WORDS = {
    "rightarrow": "->", "leftarrow": "<-", "leq": "<=", "geq": ">=", "times": "x", "tau": "tau", "mu": "mu",
    "eta": "eta", "phi": "phi", "Delta": "Delta", "zeta": "zeta", "sum": "sum", "ell": "l", "nu": "nu", "pi": "pi",
    "text": "", "mathrm": "", "vec": "",
}


def strip_tex(label: str) -> str:
    """rough LaTeX -> text for the README glossary"""
    s = label.replace("$", "")
    s = re.sub(r"\\(?:mathrm|text)\{([^}]*)\}", r"\1", s)
    s = re.sub(r"\\([A-Za-z]+)", lambda m: TEX_WORDS.get(m.group(1), m.group(1)), s)
    s = s.replace("\\", "").replace("{", "").replace("}", "")
    return re.sub(r"\s+", " ", s).strip()


def md_to_html(text: str) -> list:
    """markdown subset (headings, pipe tables, paragraphs) -> html lines; every table is closed"""
    H, table = [], []

    def flush():
        nonlocal table
        if table:
            H.append("<table>" + "".join(table) + "</table>")
            table = []

    for line in text.splitlines():
        if line.startswith("|---"):
            continue
        if line.startswith("|"):
            cells = [c.strip() for c in line.strip("|").split("|")]
            tag = "th" if not table else "td"
            table.append("<tr>" + "".join(f"<{tag}>{html.escape(c)}</{tag}>" for c in cells) + "</tr>")
            continue
        flush()
        if line.startswith("## "):
            H.append(f"<h3>{html.escape(line[3:])}</h3>")
        elif line.startswith("# "):
            continue
        elif line.strip():
            H.append(f"<p>{html.escape(line)}</p>")
    flush()
    return H


def selection_text(config, category):
    cuts = config.event_selection["categories"][category]
    sels = config.event_selection["selections"]
    return [f"`{c}`: `{sels.get(c, '?')}`" for c in cuts]


def make_index(workflow, year):
    config = WorkflowConfigBuilder(workflow).build_workflow_config()
    root = plots_root(workflow, year)
    root.mkdir(parents=True, exist_ok=True)
    hc = config.histogram_config
    channel = channel_of(workflow)
    cat_folders = {c: category_folder(config, c) for c in hc.categories}
    hlt = config.event_selection["hlt_paths"]

    # ---------------------------------------------------------------- README.md
    L = [f"# {workflow} ({year}) — plots", "",
         f"Channel: {CHANNEL_TEXT[channel]}. Samples: cH -> tautau (H+c) and bH -> tautau (H+b), "
         f"lumi-scaled to 2024. Each `<variable>.png` has two panels: yields (left) and unit-area shapes (right).",
         "", "## Where is what", "",
         "| folder | content |", "|---|---|",
         "| `<category>/<group>/<variable>.png` | one plot per variable, grouped (see below) |",
         "| `<category>/tau_wp_study/` | variables split by the tighter tau_h working points (VSe / VSjet / VSmu Tight) |",
         "| `trigger_study/` | `trigger_eff.md`/`.json` + per category `<turn-on axis>_{all_paths,union,acceptance}.png` |",
         "| `cutflows/<category folder>.csv`, `..._raw.csv` | cutflow per process: genWeight-only Σw scaled to the luminosity (no scale factors; non-monotonic in small categories because of the negative weights) and unweighted event counts |",
         "| `index.html` | browsable thumbnails of everything above |", "",
         "## Categories", "",
         "| folder | category | selection |", "|---|---|---|"]
    for cat, folder in cat_folders.items():
        L.append(f"| `{folder}` | `{cat}` | " + "; ".join(selection_text(config, cat)) + " |")
    L += ["", f"HLT paths per primary dataset (`hlt_paths`, pure OR = `pass_union`): {hlt}", "",
          "`trigger` = OR of the gated trigger regions (HLT bit + offline plateau cuts on the selected pair + TrigObj matching): "
          + "; ".join(f"`{n}`: `{e}`" for n, e in config.event_selection["selections"].items() if n == "trigger" or n.startswith("trig_")),
          "", "## Groups", ""]
    for g in ordered_groups(config):
        axes = [a for a in hc.axes if group_of(config, a) == g and not a.startswith(("pass_", "in_", "l1_", "l2_", "turnon_")) and not a.endswith(("_wp", "_2d"))]
        if g == "tau_wp_study":
            axes = [l for l in hc.layout if l.endswith("_by_tauwp")]
        # two-axis layouts (2D maps / per-category curves), named after the layout
        axes += [f"{l} (2D)" for l, v in hc.layout.items() if is_2d_layout(l, v) and group_of(config, l) == g]
        if axes:
            L.append(f"- **{g}/**: " + ", ".join(axes))
    L += ["", "Conventions: *leg 1* (\"lepton\") = the muon / electron / leading tau_h of the selected pair, "
          "*tau* = the tau_h leg; `met_*` = PuppiMET after JEC and tau-energy-scale propagation; "
          "jets 1-3 = the three highest-pT selected jets; `cjet_*` = the c-jet candidate (highest UParT CvL).",
          "", "## Glossary", "", "| file | label | definition |", "|---|---|---|"]
    for name, axis in hc.axes.items():
        if name.startswith(("pass_", "in_", "l1_", "l2_", "turnon_")) or name.endswith(("_wp", "_2d")):
            continue
        L.append(f"| `{name}` | {strip_tex(axis.label)} | {DEFINITIONS.get(name, '')} |")
    L += ["", "Turn-on axes of the trigger study (`trigger_study/`): "
          + ", ".join(f"`{n}` ({strip_tex(hc.axes[n].label)})" for n in hc.axes if n.startswith("turnon_"))]
    (root / "README.md").write_text("\n".join(L) + "\n")

    # ---------------------------------------------------------------- index.html
    def rel(p: Path):
        return html.escape(str(p.relative_to(root)))

    H = [f"<!doctype html><html><head><meta charset='utf-8'><title>{workflow} {year}</title>",
         "<style>body{font-family:sans-serif;margin:20px} details{margin:6px 0} summary{cursor:pointer;font-weight:bold}"
         " .grid{display:flex;flex-wrap:wrap;gap:8px} .grid a{display:block;text-align:center;font-size:11px;width:310px}"
         " img{width:300px;border:1px solid #ddd} table{border-collapse:collapse;font-size:12px} td,th{border:1px solid #ccc;padding:2px 6px}"
         " h2{margin-top:30px}</style></head><body>",
         f"<h1>{workflow} ({year}) — {html.escape(CHANNEL_TEXT[channel])}</h1>",
         "<p><a href='README.md'>README (glossary, categories)</a> · <a href='trigger_study/trigger_eff.md'>trigger tables (md)</a></p>"]
    # trigger tables (markdown -> html tables)
    md = root / "trigger_study" / "trigger_eff.md"
    if md.exists():
        H.append("<h2>Trigger study</h2>")
        H += md_to_html(md.read_text())
        for folder in cat_folders.values():
            d = root / "trigger_study" / folder
            pngs = sorted(d.glob("*.png")) if d.exists() else []
            if pngs:
                H.append(f"<details><summary>trigger_study / {html.escape(folder)} ({len(pngs)})</summary><div class='grid'>")
                H += [f"<a href='{rel(p)}'><img src='{rel(p)}' loading='lazy'>{html.escape(p.stem)}</a>" for p in pngs]
                H.append("</div></details>")
    # per category / group
    H.append("<h2>Distributions</h2>")
    for cat, folder in cat_folders.items():
        d = root / folder
        if not d.exists():
            continue
        H.append(f"<details><summary>{html.escape(folder)}  (category <code>{html.escape(cat)}</code>)</summary>")
        for suffix, name in (("", "cutflow (Σw, lumi-scaled, genWeight only)"), ("_raw", "cutflow (unweighted events)")):
            cf = root / "cutflows" / f"{folder}{suffix}.csv"
            if cf.exists():
                H.append(f"<details><summary>{name}</summary><table>")
                with open(cf) as f:
                    for i, row in enumerate(csv.reader(f)):
                        tag = "th" if i == 0 else "td"
                        H.append("<tr>" + "".join(f"<{tag}>{html.escape(c[:10] if i else c)}</{tag}>" for c in row) + "</tr>")
                H.append("</table></details>")
        for g in ordered_groups(config):
            pngs = sorted((d / g).glob("*.png")) if (d / g).exists() else []
            if not pngs:
                continue
            H.append(f"<details><summary>{html.escape(g)} ({len(pngs)})</summary><div class='grid'>")
            H += [f"<a href='{rel(p)}'><img src='{rel(p)}' loading='lazy'>{html.escape(p.stem)}</a>" for p in pngs]
            H.append("</div></details>")
        H.append("</details>")
    H.append("</body></html>")
    (root / "index.html").write_text("\n".join(H))
    logging.info(f"index written: {root / 'README.md'}, {root / 'index.html'}")

    # ---------------------------------------------------------------- top-level index
    top = ["<!doctype html><html><head><meta charset='utf-8'><title>H->tautau plots</title></head><body>",
           "<h1>cH, H&rarr;&tau;&tau; plots</h1><ul>"]
    for idx in sorted(PLOTS_ROOT.glob("*/*/index.html")):
        wf, yr = idx.parent.parent.name, idx.parent.name
        top.append(f"<li><a href='{wf}/{yr}/index.html'>{wf} ({yr})</a> — <a href='{wf}/{yr}/README.md'>README</a></li>")
    top.append("</ul></body></html>")
    (PLOTS_ROOT / "index.html").write_text("\n".join(top))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    p = argparse.ArgumentParser()
    p.add_argument("-w", "--workflow", required=True)
    p.add_argument("-y", "--year", default="2024")
    a = p.parse_args()
    make_index(a.workflow, a.year)
