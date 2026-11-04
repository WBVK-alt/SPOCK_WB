"""Direction map: which arrows of the amino-acid pathway map can run in both directions.

Per-gene reversibility below is textbook biochemistry (physiological direction; ATP-coupled, decarboxylating and
other strongly exergonic steps = one-way) [interpretation] — NOT KG data; the KG's reaction records were not used.
An arrow is 'both' if every gene on it is reversible, 'forward' if any is one-way; degradation arrows are the
reverse route (amino acid -> central metabolism). Gene presence (n/45) from grounding/data/map_gene_presence.csv.
Writes grounding/data/route_direction.csv, grounding/reports/direction_map.html, grounding/figures/direction_map.svg.
"""
from pathlib import Path
import html
import sys
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pathway_steps import NODES, DEG_NODES, ROUTES, CX, CY  # noqa: E402

D, R, F = HERE.parent / "data", HERE.parent / "reports", HERE.parent / "figures"
REV = {  # True = reversible step
    "rbcL": 0, "rbcS": 0, "pgk": 1, "gap2": 1, "tktA": 1, "prs": 0, "gpm": 1, "eno": 1, "pyk": 0,
    "pdhA": 0, "pdhB": 0, "aceF": 0, "lpd": 0, "ppc": 0, "gltA": 0, "acnB": 1, "icd": 0, "sucA": 0, "sucB": 0,
    "sucC": 1, "sucD": 1, "sdhA": 1, "sdhB": 1, "fumC": 1, "mdh/mqo": 0,
    "glsF/gltB": 0, "glnA": 0, "proB": 0, "proA": 1, "proC": 0,
    "argJ": 1, "argB": 0, "argC": 1, "argD": 1, "carA": 0, "carB": 0, "argF": 0, "argG": 0, "argH": 1,
    "aspC": 1, "asnB": 0, "lysC": 0, "asd": 1, "dapA": 0, "dapB": 0, "dapL": 1, "dapF": 1, "lysA": 0,
    "hom": 1, "thrB": 0, "thrC": 0, "metX/metA": 0, "metY/metB": 0, "metE/metH": 0,
    "ilvA": 0, "ilvB": 0, "ilvN": 0, "ilvC": 0, "ilvD": 0, "ilvE": 1,
    "serA": 1, "serC/PSTA": 1, "serB/psp": 0, "glyA": 1, "cysE": 0, "cysK": 0,
    "ald/alaA/AGXT": 1, "leuA": 0, "leuC": 1, "leuD": 1, "leuB": 0,
    "aroG": 0, "aroB": 0, "aroQ": 1, "aroE": 1, "aroK": 0, "aroA": 1, "aroC": 0, "CM": 0, "pheA": 0, "tyrA": 0,
    "trpE": 0, "trpG": 0, "trpD": 0, "trpC": 0, "trpF": 0, "trpA": 0, "trpB": 0,
    "hisG": 0, "hisIE": 0, "hisA": 0, "hisH": 0, "hisF": 0, "hisB": 0, "hisC": 1, "hisN": 0, "hisD": 0,
    "gcvT": 1, "gcvP": 1, "gcvH": 1, "sdaA": 0, "dadA": 0, "ltaE/tdh": 0, "aspA": 1, "ansA/iaaA": 0, "glsA": 0,
    "gdhA": 1, "putA": 0, "speA": 0, "speB": 0, "rocF": 0, "hutH": 0, "tnaA": 1, "phhA": 0, "hpd": 0,
    "bkdA1": 0, "bkdA2": 0, "bkdB": 0, "ldcC": 0, "mdeA": 0,
}
gp = pd.read_csv(D / "map_gene_presence.csv", encoding="utf-8")
npres = {(r.route, r.gene): int(r.n_strains) for r in gp.itertuples()}
N = int(gp["of"].iloc[0])

rows = []
for r in ROUTES:
    genes = [g for g, _ in r["genes"]]
    missing = [g for g in genes if g not in REV]
    assert not missing, missing
    nrev = sum(REV[g] for g in genes)
    cls = "exit" if r["kind"] == "deg" else ("both" if nrev == len(genes) else "forward")
    rows.append(dict(route=r["id"], kind=r["kind"], src=r["src"], dst=r["dst"], direction=cls, n_steps=len(genes),
                     n_reversible=nrev, genes=" ".join(f'{g}{"⇄" if REV[g] else "→"}({npres[(r["id"], g)]})' for g in genes),
                     present_any=any(npres[(r["id"], g)] > 0 for g in genes)))
rd = pd.DataFrame(rows)
rd.to_csv(D / "route_direction.csv", index=False, encoding="utf-8")

INK, INK2, MUTED, RULE, SURF = "#1f1f1d", "#52514e", "#8a8984", "#d6d5d0", "#fcfcfb"
BOTH, FWD, EXIT = "#61409b", "#1f1f1d", "#8a8984"
W, H = 1780, 1290
esc = html.escape
box = lambda l, k: (7.2 * len(l) + (22 if k == "aa" else 18), 28 if k == "aa" else 24)
node_xy = lambda n: NODES[n][0] if n in NODES else DEG_NODES[n]


def trim(p, q, name, pad=6):
    w, h = box(name, NODES[name][1]) if name in NODES else (7.0 * len(name) + 8, 18)
    dx, dy = q[0] - p[0], q[1] - p[1]
    t = min((w / 2 + pad) / abs(dx) if dx else 1e9, (h / 2 + pad) / abs(dy) if dy else 1e9, 1)
    return (p[0] + dx * t, p[1] + dy * t)


o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="system-ui,Segoe UI,sans-serif">',
     f'<rect width="{W}" height="{H}" fill="{SURF}"/><defs>']
add = o.append
for nm, col in (("mb", BOTH), ("mf", FWD), ("me", EXIT), ("mx", RULE)):
    add(f'<marker id="{nm}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{col}"/></marker>')
add('</defs>')
add(f'<text x="30" y="34" font-size="20" font-weight="700" fill="{INK}">Which directions can each route run?</text>')
add(f'<text x="30" y="54" font-size="12" fill="{INK2}">Reversibility per enzyme step = textbook biochemistry (interpretation, not KG data). ⇄ = reversible step, → = one-way step; (n) = genomes carrying the gene, of {N}.</text>')
add(f'<circle cx="{CX}" cy="{CY}" r="185" fill="none" stroke="#e4e3df" stroke-width="18" stroke-opacity=".45"/>')
add(f'<text x="{CX}" y="{CY}" text-anchor="middle" font-size="16" font-weight="700" fill="{INK2}">TCA cycle</text>')
mids = []
for r, row in zip(ROUTES, rows):
    p, q, c = node_xy(r["src"]), node_xy(r["dst"]), r.get("via")
    a, b = trim(p, c or q, r["src"]), trim(q, c or p, r["dst"], pad=8)
    absent = not row["present_any"]
    col, mk = {"both": (BOTH, "mb"), "forward": (FWD, "mf"), "exit": (EXIT, "me")}[row["direction"]]
    if absent:
        col, mk = RULE, "mx"
    start = f' marker-start="url(#{mk})"' if row["direction"] == "both" else ""
    dash = ' stroke-dasharray="6 5"' if row["direction"] == "exit" else ""
    wdt = 3.2 if row["direction"] == "both" else 2
    d = f"M{a[0]:.0f},{a[1]:.0f} Q{c[0]},{c[1]} {b[0]:.0f},{b[1]:.0f}" if c else f"M{a[0]:.0f},{a[1]:.0f} L{b[0]:.0f},{b[1]:.0f}"
    add(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{wdt}"{dash}{start} marker-end="url(#{mk})"/>')
    mids.append((row, (0.25 * a[0] + 0.5 * c[0] + 0.25 * b[0], 0.25 * a[1] + 0.5 * c[1] + 0.25 * b[1]) if c else ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)))
for row, (mx, my) in mids:
    toks = row["genes"].split(" ")
    lines, cur = [], ""
    for t in toks:
        if cur and len(cur) + len(t) > 30:
            lines.append(cur); cur = ""
        cur = (cur + " " if cur else "") + t
    lines.append(cur)
    tag = {"both": "⇄ both", "forward": f'→ forward ({row["n_reversible"]}/{row["n_steps"]} steps ⇄)' if row["n_reversible"] else "→ forward",
           "exit": "exit (reverse route)"}[row["direction"]]
    colr = {"both": BOTH, "forward": FWD, "exit": EXIT}[row["direction"]] if row["present_any"] else MUTED
    add(f'<g class="c" data-t="{esc(row["src"] + " → " + row["dst"] + ": " + tag + " · " + row["genes"])}">')
    add(f'<text x="{mx}" y="{my - 3}" text-anchor="middle" font-size="11" font-weight="700" fill="{colr}" paint-order="stroke" stroke="{SURF}" stroke-width="4">{esc(tag)}</text>')
    for i, ln in enumerate(lines):
        add(f'<text x="{mx}" y="{my + 10 + 11 * i}" text-anchor="middle" font-size="9.5" fill="{INK2 if row["present_any"] else MUTED}" paint-order="stroke" stroke="{SURF}" stroke-width="3">{esc(ln)}</text>')
    add('</g>')
for name, ((x, y), kind) in NODES.items():
    w, h = box(name, kind)
    if kind == "aa":
        add(f'<rect x="{x - w / 2}" y="{y - h / 2}" width="{w}" height="{h}" rx="{h / 2}" fill="{INK}"/><text x="{x}" y="{y + 5}" text-anchor="middle" font-size="14" font-weight="700" fill="#fff">{esc(name)}</text>')
    else:
        add(f'<rect x="{x - w / 2}" y="{y - h / 2}" width="{w}" height="{h}" rx="6" fill="#fff" stroke="{MUTED if kind == "central" else RULE}" stroke-width="1.3"/>'
            f'<text x="{x}" y="{y + 4.5}" text-anchor="middle" font-size="12.5" font-weight="{600 if kind == "central" else 400}" fill="{INK if kind == "central" else INK2}">{esc(name)}</text>')
for name, (x, y) in DEG_NODES.items():
    add(f'<text x="{x}" y="{y + 4}" text-anchor="middle" font-size="11.5" font-style="italic" fill="{INK2}">{esc(name)}</text>')
lx, ly = 30, 1110
for i, (col, txt, kind) in enumerate([(BOTH, "both directions: every step reversible", "b"), (FWD, "forward only: ≥1 one-way step", "f"),
                                      (EXIT, "exit / reverse route (degradation)", "e"), (RULE, "no gene found in any genome", "x")]):
    y = ly + i * 20
    add(f'<line x1="{lx}" x2="{lx + 44}" y1="{y}" y2="{y}" stroke="{col}" stroke-width="{3.2 if kind == "b" else 2}"'
        f'{" stroke-dasharray=\"6 5\"" if kind == "e" else ""}{" marker-start=\"url(#mb)\"" if kind == "b" else ""} marker-end="url(#m{kind})"/>'
        f'<text x="{lx + 54}" y="{y + 4}" font-size="12" fill="{INK2}">{txt}</text>')
nb = rd[(rd.kind == "bio")]
add(f'<text x="{lx}" y="{ly + 92}" font-size="12" fill="{INK}">Biosynthesis arrows: {int((nb.direction == "both").sum())} both directions, '
    f'{int((nb.direction == "forward").sum())} forward only (of {len(nb)}).</text>')
add('</svg>')
svg = "\n".join(o)
(F / "direction_map.svg").write_text(svg, encoding="utf-8")
page = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>Prochlorococcus amino-acid pathways — direction</title>
<style>body{{background:{SURF};margin:16px 24px;font:14px system-ui,Segoe UI,sans-serif}}.fig{{overflow:auto}}
#tip{{position:fixed;pointer-events:none;background:#fff;border:1px solid {RULE};border-radius:6px;padding:6px 9px;font-size:12px;max-width:560px;box-shadow:0 4px 14px rgba(0,0,0,.15);display:none}}</style></head>
<body><div class="fig">{svg}</div><div id="tip"></div><script>const t=document.getElementById("tip");document.querySelectorAll(".c").forEach(r=>{{
r.addEventListener("mousemove",e=>{{t.style.display="block";t.textContent=r.dataset.t;t.style.left=Math.min(e.clientX+12,innerWidth-580)+"px";t.style.top=(e.clientY+12)+"px"}});
r.addEventListener("mouseleave",()=>t.style.display="none")}});</script></body></html>"""
(R / "direction_map.html").write_text(page, encoding="utf-8")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print(rd.groupby(["kind", "direction"]).size().to_string())
print(rd[rd.direction == "both"][["route", "src", "dst", "genes"]].to_string())
