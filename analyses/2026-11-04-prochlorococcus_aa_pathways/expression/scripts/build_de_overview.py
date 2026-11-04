"""Overview: amino-acid pathway (and central carbon) genes x experiments -> significantly up / down / not significant / not measured.

Rows = genes on the pathway map (grounding/scripts/pathway_steps.py), in pathway order, grouped by family; genes absent
from every strain with expression data are left out. Columns = the 67 experiments, grouped by treatment type.
A cell aggregates every timepoint and every paralog of that gene in that experiment's strain.
Reads expression/data/map_de_rows.csv, experiments.csv, genome_wide_status.csv, grounding/data/ko_presence.csv.
Writes expression/data/de_overview_cells.csv, expression/reports/de_overview.html, expression/figures/de_overview.svg.
"""
from pathlib import Path
import html
import sys
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
G, D, R, F = HERE.parent / "grounding" / "data", HERE / "data", HERE / "reports", HERE / "figures"
F.mkdir(exist_ok=True)
sys.path.insert(0, str(HERE.parent / "grounding" / "scripts"))
from pathway_steps import ROUTES  # noqa: E402

RT = {r["id"]: r for r in ROUTES}
SECTIONS = [
    ("Calvin cycle, glycolysis, TCA", ["rubisco", "pgk", "tkt", "prs", "lowglyc", "pyk", "pdh", "ppc", "gltA", "acn", "icd",
                                       "ogdh", "scs", "sdh", "fum", "mdh"]),
    ("Glutamate family (Glu, Gln, Pro, Arg)", ["gogat", "gln", "pro", "arg"]),
    ("Aspartate family (Asp, Asn, Lys, Thr, Met, Ile)", ["asp", "asn", "asa", "lys", "hom", "thr", "met", "ile"]),
    ("Serine family (Ser, Gly, Cys)", ["serine", "glyA", "cys"]),
    ("Pyruvate family (Ala, Val, Leu)", ["ala", "kiv", "val", "leu"]),
    ("Aromatic family (Phe, Tyr, Trp)", ["shikimate", "phe", "tyr", "trp"]),
    ("Histidine", ["his"]),
    ("Degradation / exits", [r["id"] for r in ROUTES if r["kind"] == "deg"]),
]
short = lambda o: o.replace("Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)", "SS120").replace("Prochlorococcus ", "")

pres = pd.read_csv(G / "ko_presence.csv", encoding="utf-8")
de = pd.read_csv(D / "map_de_rows.csv", encoding="utf-8")
exp = pd.read_csv(D / "experiments.csv", encoding="utf-8")
gw = pd.read_csv(D / "genome_wide_status.csv", encoding="utf-8")
gw["frac"] = (gw.significant_up + gw.significant_down) / gw.genes_measured
gwmax = gw.groupby("experiment_id").frac.max()
expr_orgs = set(exp.organism_name)
tags_of = pres.groupby(["organism", "ko_id"]).locus_tag.apply(set).to_dict()
de_by = {k: g for k, g in de.groupby(["experiment_id", "locus_tag"])}

# rows: unique genes in section order
rows, seen = [], set()
for sec, rids in SECTIONS:
    for rid in rids:
        r = RT[rid]
        for label, alts in r["genes"]:
            key = (label, tuple(alts))
            if key in seen:
                continue
            seen.add(key)
            have = {o: set().union(*[tags_of.get((o, k), set()) for k in alts]) for o in expr_orgs}
            if not any(have.values()):
                continue
            rows.append(dict(section=sec, gene=label, step=f'{r["src"]} → {r["dst"]}', kind=r["kind"], tags=have,
                             kos=";".join(k.split(":")[-1] for k in alts)))

# columns
exp["tt"] = exp.treatment_type.str.split(";").str[0]
exp["org"] = exp.organism_name.map(short)
TT_ORDER = ["nitrogen", "phosphorus", "iron", "carbon", "light", "darkness", "coculture", "viral", "salt", "plastic", "chemical", "growth_phase"]
exp["tto"] = exp.tt.map({t: i for i, t in enumerate(TT_ORDER)}).fillna(99)
exp = exp.sort_values(["tto", "org", "omics_type", "experiment_name"]).reset_index(drop=True)
cols = exp.to_dict("records")

cells = []
for i, r in enumerate(rows):
    for j, e in enumerate(cols):
        tg = r["tags"].get(e["organism_name"], set())
        parts = [de_by[(e["experiment_id"], t)] for t in tg if (e["experiment_id"], t) in de_by]
        sub = pd.concat(parts) if parts else pd.DataFrame()
        if not tg:
            st = "no_gene"
        elif sub.empty:
            st = "not_measured"
        else:
            u, d = (sub.expression_status == "significant_up").any(), (sub.expression_status == "significant_down").any()
            st = "both" if u and d else "up" if u else "down" if d else "ns"
        det = "" if sub.empty else "; ".join(
            f'{x.locus_tag} {x.timepoint if isinstance(x.timepoint, str) else ""}: {x.expression_status.replace("significant_", "").replace("not_significant", "ns")} ({x.log2fc:+.2f})'
            for x in sub.sort_values(["locus_tag", "timepoint_order"]).itertuples())
        cells.append(dict(row=i, col=j, section=r["section"], gene=r["gene"], experiment_id=e["experiment_id"],
                          organism=e["org"], treatment_type=e["tt"], status=st, detail=det))
cdf = pd.DataFrame(cells)
cdf.to_csv(D / "de_overview_cells.csv", index=False, encoding="utf-8")

# ---------------------------------------------------------------- SVG
INK, INK2, MUTED, RULE, SURF = "#1f1f1d", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
UP, DN, NS = "#eb6834", "#2a78d6", "#e6e5e1"
cw, ch, left, top = 14, 13, 330, 330
sec_gap, tt_gap = 16, 6
xs, x, last = [], left, None
for e in cols:
    if last is not None and e["tt"] != last:
        x += tt_gap
    xs.append(x); x += cw; last = e["tt"]
right = x
ys, y, last = [], top, None
for r in rows:
    if last is not None and r["section"] != last:
        y += sec_gap
    ys.append(y); y += ch; last = r["section"]
W, H = right + 150, y + 200
esc = html.escape
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="system-ui,Segoe UI,sans-serif">',
     f'<rect width="{W}" height="{H}" fill="{SURF}"/>']
add = o.append
add(f'<text x="20" y="28" font-size="18" font-weight="700" fill="{INK}">Which amino-acid pathway genes change, and under which conditions</text>')
add(f'<text x="20" y="48" font-size="12" fill="{INK2}">Each cell = one gene in one experiment, any timepoint, any paralog. Orange = significantly up, blue = down, split = up at one timepoint and down at another,</text>')
add(f'<text x="20" y="64" font-size="12" fill="{INK2}">grey = measured, not significant; blank = not in that experiment\'s table (◆ partial tables list only significant genes); faint dot = strain lacks the gene. Significance = each paper\'s own call.</text>')
# column group labels + experiment labels
grp = exp.groupby("tt", sort=False).indices
for tt, idx in grp.items():
    x0, x1 = xs[idx[0]], xs[idx[-1]] + cw
    add(f'<line x1="{x0}" x2="{x1}" y1="{top - 8}" y2="{top - 8}" stroke="{INK}" stroke-width="2"/>')
    add(f'<text transform="translate({(x0 + x1) / 2},{top - 14}) rotate(-90)" font-size="12" font-weight="700" fill="{INK}" text-anchor="start" dy="4">{"" if len(idx) < 2 else ""}</text>')
for j, e in enumerate(cols):
    part = "◆ " if e["table_scope"] != "all_detected_genes" else ""
    lab = f'{part}{e["org"]} · {str(e["treatment"])[:30]}'
    add(f'<text class="h" data-t="{esc(e["experiment_name"] + " — " + str(e["treatment"]) + " vs " + str(e["control"]) + " · " + e["omics_type"] + " · table: " + e["table_scope"])}" '
        f'transform="translate({xs[j] + cw / 2 + 4},{top - 14}) rotate(-62)" font-size="10" fill="{INK2}">{esc(lab)}</text>')
for tt, idx in grp.items():
    x0, x1 = xs[idx[0]], xs[idx[-1]] + cw
    add(f'<text x="{(x0 + x1) / 2}" y="{y + 18}" text-anchor="middle" font-size="11" font-weight="700" fill="{INK}">{ {"growth_phase": "growth", "chemical": "chem.", "phosphorus": "phosph."}.get(tt, tt)}</text>')
    add(f'<text x="{(x0 + x1) / 2}" y="{y + 32}" text-anchor="middle" font-size="10" fill="{MUTED}">{len(idx)}</text>')
# genome-wide context row
add(f'<text x="{left - 8}" y="{y + 52}" text-anchor="end" font-size="10.5" fill="{INK2}">genome-wide % significant (full tables)</text>')
for j, e in enumerate(cols):
    v = gwmax.get(e["experiment_id"])
    if e["table_scope"] == "all_detected_genes" and pd.notna(v):
        h = 30 * v
        add(f'<rect x="{xs[j] + 2}" y="{y + 72 - h}" width="{cw - 4}" height="{h}" fill="{MUTED}"/>')
add(f'<line x1="{left}" x2="{right}" y1="{y + 72}" y2="{y + 72}" stroke="{RULE}"/>')
# rows
last = None
for i, r in enumerate(rows):
    yy = ys[i]
    if r["section"] != last:
        add(f'<text x="20" y="{yy - 4}" font-size="12" font-weight="700" fill="{INK}">{esc(r["section"])}</text>')
        last = r["section"]
    sty = ' font-style="italic"' if r["kind"] == "deg" else ""
    add(f'<text x="{left - 8}" y="{yy + 10}" text-anchor="end" font-size="10.5"{sty} fill="{INK}">{esc(r["gene"])}</text>')
    add(f'<text x="40" y="{yy + 10}" font-size="9.5" fill="{MUTED}">{esc(r["step"][:40])}</text>')
    sub = cdf[cdf.row == i]
    nu, nd = int(sub.status.isin(["up", "both"]).sum()), int(sub.status.isin(["down", "both"]).sum())
    add(f'<text x="{right + 8}" y="{yy + 10}" font-size="10" fill="{INK2}"><tspan fill="{UP}" font-weight="700">↑{nu}</tspan>  <tspan fill="{DN}" font-weight="700">↓{nd}</tspan></text>')
    for c in sub.itertuples():
        x0, y0 = xs[c.col] + 1, yy + 1
        t = esc(f"{c.gene} · {cols[c.col]['org']} · {cols[c.col]['experiment_name']}: {c.status}" + (f" — {c.detail}" if c.detail else ""))
        w, h = cw - 2, ch - 2
        if c.status == "up":
            add(f'<rect class="c" data-t="{t}" x="{x0}" y="{y0}" width="{w}" height="{h}" rx="2" fill="{UP}"/>')
        elif c.status == "down":
            add(f'<rect class="c" data-t="{t}" x="{x0}" y="{y0}" width="{w}" height="{h}" rx="2" fill="{DN}"/>')
        elif c.status == "both":
            add(f'<g class="c" data-t="{t}"><path d="M{x0},{y0} h{w} L{x0},{y0 + h} z" fill="{UP}"/><path d="M{x0 + w},{y0} v{h} h{-w} z" fill="{DN}"/></g>')
        elif c.status == "ns":
            add(f'<rect class="c" data-t="{t}" x="{x0}" y="{y0}" width="{w}" height="{h}" rx="2" fill="{NS}"/>')
        elif c.status == "not_measured":
            add(f'<rect class="c" data-t="{t}" x="{x0}" y="{y0}" width="{w}" height="{h}" rx="2" fill="none" stroke="{RULE}" stroke-width=".8"/>')
        else:
            add(f'<circle cx="{x0 + w / 2}" cy="{y0 + h / 2}" r="1.2" fill="{MUTED}"/>')
# legend
lx = 20
for k, (col, txt) in enumerate([(UP, "up"), (DN, "down"), (None, "up + down"), (NS, "not significant"), ("none", "not in table")]):
    yy = y + 95 + k * 18
    if col is None:
        add(f'<path d="M{lx},{yy} h12 L{lx},{yy + 11} z" fill="{UP}"/><path d="M{lx + 12},{yy} v11 h-12 z" fill="{DN}"/>')
    else:
        add(f'<rect x="{lx}" y="{yy}" width="12" height="11" rx="2" fill="{col}" stroke="{RULE if col == "none" else "none"}"/>')
    add(f'<text x="{lx + 18}" y="{yy + 10}" font-size="11" fill="{INK2}">{txt}</text>')
add('</svg>')
svg = "\n".join(o)
(F / "de_overview.svg").write_text(svg, encoding="utf-8")
page = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>Prochlorococcus amino-acid pathway genes — differential expression overview</title>
<style>body{{background:{SURF};margin:16px 24px;font:14px system-ui,Segoe UI,sans-serif}}.fig{{overflow:auto}}
#tip{{position:fixed;pointer-events:none;background:#fff;border:1px solid {RULE};border-radius:6px;padding:6px 9px;font-size:12px;max-width:560px;box-shadow:0 4px 14px rgba(0,0,0,.15);display:none}}
.c:hover{{outline:1.5px solid {INK}}}</style></head><body><div class="fig">{svg}</div><div id="tip"></div>
<script>const t=document.getElementById("tip");document.querySelectorAll(".c,.h").forEach(r=>{{
r.addEventListener("mousemove",e=>{{t.style.display="block";t.textContent=r.dataset.t;t.style.left=Math.min(e.clientX+12,innerWidth-580)+"px";t.style.top=Math.min(e.clientY+12,innerHeight-120)+"px"}});
r.addEventListener("mouseleave",()=>t.style.display="none")}});</script></body></html>"""
(R / "de_overview.html").write_text(page, encoding="utf-8")
print("rows", len(rows), "cols", len(cols))
print(cdf.status.value_counts().to_dict())
s = cdf.groupby(["section", "treatment_type"]).status.apply(lambda s: f"{s.isin(['up','both']).sum()}↑ {s.isin(['down','both']).sum()}↓").unstack()
pd.set_option("display.width", 250); print(s.to_string())
