"""Pathway map: amino-acid biosynthesis/degradation routes meeting central carbon metabolism (TCA, glycolysis, PPP),
each arrow annotated with how many Prochlorococcus genomes carry all its genes.

Reads grounding/data/map_gene_presence.csv + map_route_presence.csv (counts) and pathway_steps.py (layout).
Writes grounding/reports/pathway_map.html (interactive) and grounding/figures/pathway_map.svg (static).
"""
from pathlib import Path
import html
import json
import sys
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pathway_steps import NODES, DEG_NODES, ROUTES, CX, CY  # noqa: E402

D, R, F = HERE.parent / "data", HERE.parent / "reports", HERE.parent / "figures"
R.mkdir(exist_ok=True); F.mkdir(exist_ok=True)
gp = pd.read_csv(D / "map_gene_presence.csv", encoding="utf-8")
rp = pd.read_csv(D / "map_route_presence.csv", encoding="utf-8").set_index("route")
N = int(rp["of"].iloc[0])
short = lambda o: o.replace("Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)", "SS120").replace("Prochlorococcus ", "")

INK, INK2, MUTED, RULE, SURF, CARD = "#1f1f1d", "#52514e", "#8a8984", "#d6d5d0", "#fcfcfb", "#f1f0ec"
TIERS = [(N, "#470f84", "#fff", f"all {N}"), (23, "#7b74b5", "#fff", f"23–{N - 1}"), (1, "#c9c6e3", INK, "1–22")]
def tier(n):
    for lo, fill, ink, _ in TIERS:
        if n >= lo:
            return fill, ink
    return None, MUTED

W, H = 1780, 1250
out = []
add = out.append
esc = html.escape

def box(label, kind):
    w = 7.2 * len(label) + (22 if kind == "aa" else 18)
    return w, (28 if kind == "aa" else 24)

def node_xy(name):
    if name in NODES:
        return NODES[name][0]
    return DEG_NODES[name]

def trim(p, q, name, pad=5):
    """Point where segment p->q leaves the box of node `name` centred at p."""
    if name in NODES:
        w, h = box(name, NODES[name][1])
    else:
        w, h = 7.0 * len(name) + 8, 18
    dx, dy = q[0] - p[0], q[1] - p[1]
    if dx == 0 and dy == 0:
        return p
    tx = (w / 2 + pad) / abs(dx) if dx else 1e9
    ty = (h / 2 + pad) / abs(dy) if dy else 1e9
    t = min(tx, ty, 1)
    return (p[0] + dx * t, p[1] + dy * t)

add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="system-ui,Segoe UI,sans-serif">')
add(f'<rect width="{W}" height="{H}" fill="{SURF}"/>')
add('<defs>')
for nm, col in (("a-c", MUTED), ("a-b", INK), ("a-d", INK2)):
    add(f'<marker id="{nm}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
        f'<path d="M0,0 L10,5 L0,10 z" fill="{col}"/></marker>')
add('</defs>')

# family backdrops (labels only; light panels)
FAM = [("Serine family", 60, 210, 560, 165), ("Aromatic family", 1120, 330, 600, 250), ("Histidine", 1240, 85, 360, 70),
       ("Pyruvate family", 1000, 565, 380, 205), ("Glutamate family", 1120, 975, 380, 235),
       ("Aspartate family", 190, 445, 360, 525)]
for lab, x, y, w, h in FAM:
    add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="14" fill="{CARD}"/>')
    add(f'<text x="{x + 12}" y="{y + 17}" font-size="12" fill="{MUTED}" font-weight="600">{lab}</text>')
add(f'<circle cx="{CX}" cy="{CY}" r="185" fill="none" stroke="{RULE}" stroke-width="18" stroke-opacity=".45"/>')
add(f'<text x="{CX}" y="{CY - 4}" text-anchor="middle" font-size="16" font-weight="700" fill="{INK2}">TCA cycle</text>')
add(f'<text x="{CX}" y="{CY + 16}" text-anchor="middle" font-size="11" fill="{MUTED}">(as annotated in the KG)</text>')

# edges first, nodes on top, badges last
badges = []
for r in ROUTES:
    p, q = node_xy(r["src"]), node_xy(r["dst"])
    c = r.get("via")
    a = trim(p, c or q, r["src"])
    b = trim(q, c or p, r["dst"], pad=7)
    col, mk = {"central": (MUTED, "a-c"), "bio": (INK, "a-b"), "deg": (INK2, "a-d")}[r["kind"]]
    dash = ' stroke-dasharray="6 5"' if r["kind"] == "deg" else ""
    d = f"M{a[0]:.0f},{a[1]:.0f} Q{c[0]},{c[1]} {b[0]:.0f},{b[1]:.0f}" if c else f"M{a[0]:.0f},{a[1]:.0f} L{b[0]:.0f},{b[1]:.0f}"
    add(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{1.6 if r["kind"] == "central" else 2}"{dash} marker-end="url(#{mk})"/>')
    mid = (0.25 * a[0] + 0.5 * c[0] + 0.25 * b[0], 0.25 * a[1] + 0.5 * c[1] + 0.25 * b[1]) if c else ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
    badges.append((r, mid))

for name, ((x, y), kind) in NODES.items():
    w, h = box(name, kind)
    if kind == "aa":
        add(f'<rect x="{x - w / 2}" y="{y - h / 2}" width="{w}" height="{h}" rx="{h / 2}" fill="{INK}"/>'
            f'<text x="{x}" y="{y + 5}" text-anchor="middle" font-size="14" font-weight="700" fill="#fff">{esc(name)}</text>')
    elif kind == "central":
        add(f'<rect x="{x - w / 2}" y="{y - h / 2}" width="{w}" height="{h}" rx="6" fill="#fff" stroke="{MUTED}" stroke-width="1.4"/>'
            f'<text x="{x}" y="{y + 4.5}" text-anchor="middle" font-size="12.5" font-weight="600" fill="{INK}">{esc(name)}</text>')
    else:
        add(f'<rect x="{x - w / 2}" y="{y - h / 2}" width="{w}" height="{h}" rx="6" fill="#fff" stroke="{RULE}" stroke-width="1.2"/>'
            f'<text x="{x}" y="{y + 4.5}" text-anchor="middle" font-size="12" fill="{INK2}">{esc(name)}</text>')
for name, (x, y) in DEG_NODES.items():
    add(f'<text x="{x}" y="{y + 4}" text-anchor="middle" font-size="11.5" font-style="italic" fill="{INK2}">{esc(name)}</text>')

tips = {}
for r, (mx, my) in badges:
    rid = r["id"]
    n = int(rp.loc[rid, "n_complete"])
    fill, ink = tier(n)
    g = gp[gp.route == rid]
    genes_txt = [f'{x.gene}' + ("" if x.n_strains == N else f" ({x.n_strains})") for x in g.itertuples()]
    note = r.get("note")
    lines, cur = [], ""
    for t in genes_txt:
        if len(cur) + len(t) > 34 and cur:
            lines.append(cur); cur = ""
        cur = (cur + " · " if cur else "") + t
    lines.append(cur)
    if note:
        lines.append(note)
    lab = f"{n}/{N}"
    bw = 8 * len(lab) + 12
    add(f'<g class="b" data-r="{rid}" style="cursor:help">')
    if fill:
        add(f'<rect x="{mx - bw / 2}" y="{my - 10}" width="{bw}" height="20" rx="10" fill="{fill}" stroke="{SURF}" stroke-width="2"/>')
    else:
        add(f'<rect x="{mx - bw / 2}" y="{my - 10}" width="{bw}" height="20" rx="10" fill="#fff" stroke="{MUTED}" stroke-width="1.2" stroke-dasharray="3 2"/>')
    add(f'<text x="{mx}" y="{my + 4.5}" text-anchor="middle" font-size="12" font-weight="700" fill="{ink}">{lab}</text>')
    for i, ln in enumerate(lines):
        add(f'<text x="{mx}" y="{my + 23 + 12 * i}" text-anchor="middle" font-size="10.5" fill="{INK2}" '
            f'paint-order="stroke" stroke="{SURF}" stroke-width="3">{esc(ln)}</text>')
    add('</g>')
    tips[rid] = dict(title=f'{r["src"]} → {r["dst"]}', kind=r["kind"], n=n, N=N,
                     incomplete=[short(o) for o in str(rp.loc[rid, "incomplete_in"]).split(";") if o and o != "nan"],
                     genes=[dict(g=x.gene, n=int(x.n_strains), ko=x.per_ko, names=x.kg_names,
                                 miss=[short(o) for o in str(x.missing_in).split(";") if o and o != "nan"]) for x in g.itertuples()])

# legend
lx, ly = 30, 1095
add(f'<text x="{lx}" y="{ly}" font-size="13" font-weight="700" fill="{INK}">Badge = genomes carrying every gene of the arrow (of {N} Prochlorococcus genomes)</text>')
for i, (lo, fill, ink, txt) in enumerate(TIERS + [(0, None, MUTED, "0 (no gene found)")]):
    x = lx + i * 150
    if fill:
        add(f'<rect x="{x}" y="{ly + 12}" width="34" height="18" rx="9" fill="{fill}"/>')
    else:
        add(f'<rect x="{x}" y="{ly + 12}" width="34" height="18" rx="9" fill="#fff" stroke="{MUTED}" stroke-dasharray="3 2"/>')
    add(f'<text x="{x + 42}" y="{ly + 25}" font-size="12" fill="{INK2}">{txt}</text>')
ly2 = ly + 52
for i, (txt, col, dash) in enumerate([("central carbon metabolism", MUTED, ""), ("amino-acid biosynthesis", INK, ""),
                                       ("degradation / exit", INK2, ' stroke-dasharray="6 5"')]):
    x = lx + i * 230
    add(f'<line x1="{x}" y1="{ly2}" x2="{x + 40}" y2="{ly2}" stroke="{col}" stroke-width="2"{dash}/>'
        f'<text x="{x + 48}" y="{ly2 + 4}" font-size="12" fill="{INK2}">{txt}</text>')
add(f'<text x="{lx}" y="{ly2 + 28}" font-size="11.5" fill="{INK2}">Gene names under a badge; "(n)" = carried by only n genomes.</text>')
add(f'<text x="{lx}" y="{ly2 + 44}" font-size="11.5" fill="{INK2}">Multi-step arrows collapse a linear route. Layout = textbook biochemistry;</text>')
add(f'<text x="{lx}" y="{ly2 + 60}" font-size="11.5" fill="{INK2}">counts = KEGG KO annotations in the KG (0 can be an annotation gap).</text>')
add('</svg>')
svg = "\n".join(out)
(F / "pathway_map.svg").write_text(svg, encoding="utf-8")

table = gp.assign(missing=gp.missing_in.fillna("").map(lambda s: ", ".join(short(o) for o in s.split(";") if o)))[
    ["route", "kind", "src", "dst", "gene", "n_strains", "of", "per_ko", "kg_names", "missing"]].values.tolist()
page = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>Prochlorococcus amino-acid pathway map — gene presence across genomes</title>
<style>body{{background:{SURF};color:{INK};font:14px/1.45 system-ui,Segoe UI,sans-serif;margin:24px 32px}}
h1{{font-size:21px;margin:0 0 4px}}h2{{font-size:16px;margin:28px 0 6px}}.meta{{color:{MUTED};font-size:12px;max-width:1150px}}
.fig{{overflow:auto;border:1px solid {RULE};border-radius:10px}}
#tip{{position:fixed;pointer-events:none;background:#fff;color:{INK};border:1px solid {RULE};border-radius:6px;padding:8px 10px;font-size:12px;max-width:560px;box-shadow:0 4px 14px rgba(0,0,0,.15);display:none;z-index:9}}
table{{border-collapse:collapse;font-size:12px;width:100%}}td,th{{border-bottom:1px solid {RULE};padding:3px 6px;text-align:left;vertical-align:top}}th{{color:{INK2};position:sticky;top:0;background:{SURF}}}
.tw{{max-height:520px;overflow:auto;border:1px solid {RULE};border-radius:8px}} g.b:hover rect{{stroke:{INK};stroke-width:2}}</style></head><body>
<h1>How <i>Prochlorococcus</i> amino-acid pathways meet central carbon metabolism</h1>
<div class="meta">Each arrow is one reaction or a linear run of reactions. The badge counts how many of the {N} <i>Prochlorococcus</i> genomes in the KG carry <b>every</b> gene on that arrow (KEGG KO annotation; one gene, PSTA, counted by its Cyanorak ortholog group because its KO in the KG is a hydrogenase). Where an enzyme has known alternative forms (e.g. mdh or mqo), any one counts. Hover a badge for per-gene counts, KOs and the genomes that lack it. Pathway layout is textbook biochemistry, not KG data. A 0 means no gene with that KO is annotated, which can also be an annotation gap.</div>
<div class="fig">{svg}</div>
<h2>Every gene on the map</h2>
<div class="tw"><table id="t"></table></div><div id="tip"></div>
<script>
const TIPS={json.dumps(tips)}, TAB={json.dumps(table)};
const tip=document.getElementById("tip");
document.querySelectorAll("g.b").forEach(g=>{{const d=TIPS[g.dataset.r];
 g.addEventListener("mousemove",e=>{{tip.style.display="block";tip.style.left=Math.min(e.clientX+14,innerWidth-580)+"px";tip.style.top=Math.min(e.clientY+14,innerHeight-260)+"px";
 tip.innerHTML=`<b>${{d.title}}</b> <span style="color:{MUTED}">(${{d.kind}})</span><br>complete in <b>${{d.n}}/${{d.N}}</b> genomes`+(d.incomplete.length&&d.incomplete.length<20?` — incomplete in: ${{d.incomplete.join(", ")}}`:"")+
 `<table style="margin-top:6px">`+d.genes.map(x=>`<tr><td><b>${{x.g}}</b></td><td>${{x.n}}/${{d.N}}</td><td style="color:{INK2}">${{x.ko}}</td><td style="color:{MUTED}">${{x.miss.length&&x.miss.length<12?"missing: "+x.miss.join(", "):""}}</td></tr>`).join("")+`</table>`;}});
 g.addEventListener("mouseleave",()=>tip.style.display="none");}});
const H=["route","type","from","to","gene","genomes","of","per-KO genome counts","KG names of the KOs","missing in"];
document.getElementById("t").innerHTML="<tr>"+H.map(h=>`<th>${{h}}</th>`).join("")+"</tr>"+TAB.map(r=>"<tr>"+r.map(v=>`<td>${{v}}</td>`).join("")+"</tr>").join("");
</script></body></html>"""
(R / "pathway_map.html").write_text(page, encoding="utf-8")
print("wrote", R / "pathway_map.html", "and", F / "pathway_map.svg")
