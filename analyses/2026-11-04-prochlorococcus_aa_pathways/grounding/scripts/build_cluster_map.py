"""Cluster map: amino-acid pathway genes (rows, in path order per amino acid) x 45 Prochlorococcus genomes (columns,
hierarchically clustered on gene presence). Top panel: is each amino acid's biosynthesis route complete per genome.

Reads grounding/data/ko_presence.csv, genomes.csv and the step definitions in pathway_steps.py. No KG calls.
Writes grounding/data/cluster_map_matrix.csv, grounding/reports/cluster_map.html, grounding/figures/cluster_map.svg.
"""
from pathlib import Path
import html
import json
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.spatial.distance import pdist

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pathway_steps import ROUTES  # noqa: E402

D, R, F = HERE.parent / "data", HERE.parent / "reports", HERE.parent / "figures"
RT = {r["id"]: r for r in ROUTES}

# biosynthesis path (route ids in order) and degradation routes per amino acid [interpretation: textbook routes]
AA = [
    ("Glu", ["gogat"], ["gdh"]), ("Gln", ["gogat", "gln"], ["glsA"]), ("Pro", ["gogat", "pro"], ["putA"]),
    ("Arg", ["gogat", "arg"], ["spe", "rocF"]),
    ("Asp", ["asp"], ["aspA"]), ("Asn", ["asp", "asn"], ["ans"]), ("Lys", ["asp", "asa", "lys"], ["ldc"]),
    ("Thr", ["asp", "asa", "hom", "thr"], ["thrdeg"]), ("Met", ["asp", "asa", "hom", "met"], ["mdeA"]),
    ("Ile", ["asp", "asa", "hom", "thr", "ile"], ["bkd"]),
    ("Ser", ["serine"], ["sdaA"]), ("Gly", ["serine", "glyA"], ["gcv"]), ("Cys", ["serine", "cys"], []),
    ("Ala", ["ala"], ["dadA"]), ("Val", ["kiv", "val"], ["bkd"]), ("Leu", ["kiv", "leu"], ["bkd"]),
    ("Phe", ["shikimate", "phe"], ["phhA"]), ("Tyr", ["shikimate", "tyr"], ["hpd"]), ("Trp", ["shikimate", "trp"], ["tnaA"]),
    ("His", ["prs", "his"], ["hutH"]),
]

pres = pd.read_csv(D / "ko_presence.csv", encoding="utf-8")
gen = pd.read_csv(D / "genomes.csv", encoding="utf-8").set_index("organism_name")
genomes = sorted(pres.organism.unique())
N = len(genomes)
clade = {o: (gen.clade.get(o) if o in gen.index else None) for o in genomes}
clade = {o: (c if isinstance(c, str) else "?") for o, c in clade.items()}
short = lambda o: o.replace("Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)", "SS120").replace("Prochlorococcus ", "")
tags = pres.groupby(["organism", "ko_id"]).locus_tag.apply(set).to_dict()

def copies(o, alts):
    s = set()
    for k in alts:
        s |= tags.get((o, k), set())
    return s

rows = []  # one row per (amino acid, part, route, gene)
for aa, bio, deg in AA:
    for part, rids in (("biosynthesis", bio), ("degradation", deg)):
        for rid in rids:
            r = RT[rid]
            for label, alts in r["genes"]:
                lt = [copies(o, alts) for o in genomes]
                rows.append(dict(aa=aa, part=part, route=rid, step=f'{r["src"]} → {r["dst"]}', gene=label,
                                 kos=";".join(k.split(":")[1] for k in alts), n=[len(x) for x in lt],
                                 tags=[",".join(sorted(x)) for x in lt]))

# completeness: all biosynthesis genes of the amino acid present
comp = {aa: [all(r["n"][j] > 0 for r in rows if r["aa"] == aa and r["part"] == "biosynthesis") for j in range(N)] for aa, _, _ in AA}
missing_genes = {aa: [[r["gene"] for r in rows if r["aa"] == aa and r["part"] == "biosynthesis" and r["n"][j] == 0] for j in range(N)]
                 for aa, _, _ in AA}

# cluster genomes on the binary presence of distinct (route, gene) rows
uniq = {}
for r in rows:
    uniq.setdefault((r["route"], r["gene"]), r["n"])
X = (np.array(list(uniq.values())).T > 0).astype(float)          # genomes x genes
var = X.std(axis=0) > 0
Lk = linkage(pdist(X[:, var], "jaccard") if var.any() else np.zeros(N * (N - 1) // 2), "average")
order = list(leaves_list(Lk))
dn = dendrogram(Lk, no_plot=True)

# matrix CSV
mat = pd.DataFrame([{**{k: r[k] for k in ("aa", "part", "route", "step", "gene", "kos")},
                     **{short(genomes[j]): r["n"][j] for j in order}} for r in rows])
mat.to_csv(D / "cluster_map_matrix.csv", index=False, encoding="utf-8")
pd.DataFrame({aa: [int(c[j]) for j in order] for aa, c in comp.items()}, index=[short(genomes[j]) for j in order]).T \
    .to_csv(D / "cluster_map_completeness.csv", encoding="utf-8")

# ---------------------------------------------------------------- SVG
INK, INK2, MUTED, RULE, SURF = "#1f1f1d", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
P1, P2, ABS, MISS = "#7b74b5", "#470f84", "#ffffff", "#c9c6e3"
CL = {"HLI": "#2a9d8f", "HLII": "#8ab17d", "HLIII": "#e9c46a", "HLIV": "#c9a227", "LLI": "#6d597a", "LLII": "#b56576",
      "LLIII": "#355070", "LLIV": "#1d3557", "?": "#bbbbbb"}
cw, ch, left, top = 17, 13, 360, 150
dend_h = 60
nrow_c = len(AA)
y_comp = top + dend_h + 30
y_genes = y_comp + nrow_c * ch + 40
gap_aa = 10
# row y positions
ys, y, last = [], y_genes, None
for r in rows:
    if last is not None and (r["aa"], r["part"]) != last:
        y += gap_aa if r["aa"] != last[0] else 4
    ys.append(y); y += ch; last = (r["aa"], r["part"])
H = y + 40
W = left + N * cw + 220
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="system-ui,Segoe UI,sans-serif">',
     f'<rect width="{W}" height="{H}" fill="{SURF}"/>']
add = o.append
esc = html.escape
add(f'<text x="20" y="28" font-size="18" font-weight="700" fill="{INK}">Amino-acid pathway genes across {N} Prochlorococcus genomes</text>')
add(f'<text x="20" y="48" font-size="12" fill="{INK2}">Rows: genes in pathway order for each amino acid (biosynthesis, then degradation in italics). Columns: genomes, clustered by gene presence (average linkage, Jaccard).</text>')
add(f'<text x="20" y="64" font-size="12" fill="{INK2}">A gene shared by several amino acids (e.g. ilvE) appears in each of their paths. Gene = any of the listed KOs (PSTA by ortholog group). Source: KEGG KO annotations in the KG.</text>')
# dendrogram
icoord, dcoord = np.array(dn["icoord"]), np.array(dn["dcoord"])
dmax = dcoord.max() or 1
for xs, ds in zip(icoord, dcoord):
    pts = [(left + (x - 5) / 10 * cw + cw / 2, top + dend_h - d / dmax * dend_h) for x, d in zip(xs, ds)]
    add('<polyline points="' + " ".join(f"{a:.1f},{b:.1f}" for a, b in pts) + f'" fill="none" stroke="{MUTED}" stroke-width="1"/>')
# column labels + clade strip
for i, j in enumerate(order):
    x = left + i * cw
    add(f'<rect x="{x + 1}" y="{top + dend_h + 4}" width="{cw - 2}" height="10" fill="{CL.get(clade[genomes[j]], "#bbb")}"/>')
    add(f'<text transform="translate({x + cw / 2 + 4},{top - 4}) rotate(-60)" font-size="11" fill="{INK}">{esc(short(genomes[j]))}</text>')
add(f'<text x="{left - 8}" y="{top + dend_h + 13}" text-anchor="end" font-size="11" fill="{INK2}">clade</text>')
# completeness panel
add(f'<text x="20" y="{y_comp - 10}" font-size="13" font-weight="700" fill="{INK}">Biosynthesis route complete?</text>')
for a, (aa, _, _) in enumerate(AA):
    yy = y_comp + a * ch
    nc = sum(comp[aa])
    add(f'<text x="{left - 8}" y="{yy + 10}" text-anchor="end" font-size="11" font-weight="700" fill="{INK}">{aa}  <tspan fill="{INK2}" font-weight="400">{nc}/{N}</tspan></text>')
    for i, j in enumerate(order):
        ok = comp[aa][j]
        miss = ", ".join(missing_genes[aa][j])
        t = f"{aa} biosynthesis in {short(genomes[j])}: " + ("complete" if ok else f"missing {miss}")
        add(f'<rect class="c" data-t="{esc(t)}" x="{left + i * cw + 1}" y="{yy + 1}" width="{cw - 2}" height="{ch - 2}" rx="2" '
            f'fill="{P2 if ok else MISS}"/>')
        if not ok:
            add(f'<text x="{left + i * cw + cw / 2}" y="{yy + 10}" text-anchor="middle" font-size="10" font-weight="700" fill="{INK}" pointer-events="none">✕</text>')
# gene panel
add(f'<text x="20" y="{y_genes - 14}" font-size="13" font-weight="700" fill="{INK}">Genes (number = gene copies when more than one)</text>')
last = None
for r, yy in zip(rows, ys):
    if r["aa"] != (last or {}).get("aa"):
        add(f'<text x="20" y="{yy + 10}" font-size="13" font-weight="700" fill="{INK}">{r["aa"]}</text>')
        add(f'<line x1="20" x2="{left + N * cw}" y1="{yy - gap_aa / 2}" y2="{yy - gap_aa / 2}" stroke="{RULE}"/>')
    npres = sum(1 for v in r["n"] if v)
    partial = npres < N
    lab = f'{r["gene"]}'
    sty = ' font-style="italic"' if r["part"] == "degradation" else ""
    add(f'<text x="{left - 52}" y="{yy + 10}" text-anchor="end" font-size="10.5"{sty} fill="{INK if partial else INK2}" '
        f'font-weight="{700 if partial else 400}">{esc(lab)}</text>')
    add(f'<text x="{left - 8}" y="{yy + 10}" text-anchor="end" font-size="10.5" fill="{INK if partial else MUTED}">{npres}/{N}</text>')
    add(f'<text x="70" y="{yy + 10}" font-size="9.5" fill="{MUTED}">{esc(r["step"][:38])}</text>')
    for i, j in enumerate(order):
        v = r["n"][j]
        t = f'{r["aa"]} {r["part"]} · {r["gene"]} ({r["kos"]}) · {r["step"]} · {short(genomes[j])} ({clade[genomes[j]]}): ' + \
            (f'{v} copy: {r["tags"][j]}' if v == 1 else f'{v} copies: {r["tags"][j]}' if v else "absent")
        fill = P2 if v > 1 else P1 if v == 1 else ABS
        add(f'<rect class="c" data-t="{esc(t)}" x="{left + i * cw + 1}" y="{yy + 1}" width="{cw - 2}" height="{ch - 2}" rx="2" fill="{fill}" '
            f'stroke="{RULE if not v else "none"}" stroke-width="1"/>')
        if v > 1:
            add(f'<text x="{left + i * cw + cw / 2}" y="{yy + 10}" text-anchor="middle" font-size="9" fill="#fff" pointer-events="none">{v}</text>')
    last = r
# legend
lx = left + N * cw + 20
items = [(P1, "present (1 copy)"), (P2, "present, ≥2 copies (number)"), (ABS, "absent"), (MISS, "route incomplete (✕)")]
for k, (c, t) in enumerate(items):
    add(f'<rect x="{lx}" y="{y_comp + k * 20}" width="14" height="12" rx="2" fill="{c}" stroke="{RULE}"/>'
        f'<text x="{lx + 20}" y="{y_comp + k * 20 + 10}" font-size="11" fill="{INK2}">{t}</text>')
for k, (c, col) in enumerate([(c, CL[c]) for c in CL if c != "?" and c in set(clade.values())]):
    add(f'<rect x="{lx}" y="{y_comp + 100 + k * 18}" width="14" height="10" fill="{col}"/>'
        f'<text x="{lx + 20}" y="{y_comp + 109 + k * 18}" font-size="11" fill="{INK2}">{c}</text>')
add(f'<text x="{lx}" y="{y_comp + 92}" font-size="11" font-weight="700" fill="{INK}">clade</text>')
add('</svg>')
svg = "\n".join(o)
(F / "cluster_map.svg").write_text(svg, encoding="utf-8")
page = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>Prochlorococcus amino-acid pathway genes — cluster map</title>
<style>body{{background:{SURF};margin:16px 24px;font:14px system-ui,Segoe UI,sans-serif}}.fig{{overflow:auto}}
#tip{{position:fixed;pointer-events:none;background:#fff;border:1px solid {RULE};border-radius:6px;padding:6px 9px;font-size:12px;max-width:520px;box-shadow:0 4px 14px rgba(0,0,0,.15);display:none}}
rect.c:hover{{stroke:{INK};stroke-width:1.5}}</style></head><body><div class="fig">{svg}</div><div id="tip"></div>
<script>const t=document.getElementById("tip");document.querySelectorAll("rect.c").forEach(r=>{{
r.addEventListener("mousemove",e=>{{t.style.display="block";t.textContent=r.dataset.t;t.style.left=Math.min(e.clientX+12,innerWidth-540)+"px";t.style.top=(e.clientY+12)+"px"}});
r.addEventListener("mouseleave",()=>t.style.display="none")}});</script></body></html>"""
(R / "cluster_map.html").write_text(page, encoding="utf-8")
print("rows", len(rows), "genomes", N)
print({aa: sum(c) for aa, c in comp.items()})
inc = {aa: sorted({g for j in range(N) for g in missing_genes[aa][j]}) for aa in comp}
print({aa: v for aa, v in inc.items() if v})
