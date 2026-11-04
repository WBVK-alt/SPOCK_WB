"""Cell illustration: how fixed carbon reaches 2-oxoglutarate through the open TCA cycle and meets nitrogen (GS-GOGAT).

Genome counts (n/45) from grounding/data/map_gene_presence.csv and map_route_presence.csv [KG]. The layout, the
idea that the open TCA cycle supplies 2-oxoglutarate for N assimilation, and the NtcA signal are [interpretation].
The NH4+ transporter is drawn for orientation only (not counted). Tolonen 2006 note from the KG literature index.
Writes grounding/figures/cell_carbon_nitrogen.svg and grounding/reports/cell_carbon_nitrogen.html.
"""
from pathlib import Path
import html
import pandas as pd

HERE = Path(__file__).resolve().parent
D, R, F = HERE.parent / "data", HERE.parent / "reports", HERE.parent / "figures"
gp = pd.read_csv(D / "map_gene_presence.csv", encoding="utf-8")
rp = pd.read_csv(D / "map_route_presence.csv", encoding="utf-8").set_index("route")
N = int(rp["of"].iloc[0])
G = {(r.route, r.gene): int(r.n_strains) for r in gp.itertuples()}
route_n = lambda rid: int(rp.loc[rid, "n_complete"])

INK, INK2, MUTED, RULE = "#1f1f1d", "#52514e", "#8a8984", "#d6d5d0"
C, NIT, MISS = "#2f7d4f", "#61409b", "#b5b3ad"      # carbon green, nitrogen violet, missing grey
CYTO, MEMB, THY = "#f4f7f1", "#9fb79a", "#cfe0c8"
W, H = 1500, 980
esc = html.escape
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="system-ui,Segoe UI,sans-serif">',
     '<rect width="100%" height="100%" fill="#fcfcfb"/><defs>']
add = o.append
for nm, col in (("ac", C), ("an", NIT), ("am", MISS), ("ai", INK2)):
    add(f'<marker id="{nm}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="{col}"/></marker>')
add('</defs>')

# --- cell -----------------------------------------------------------------------------------------
cx, cy, rx, ry = 760, 520, 650, 400
add(f'<ellipse cx="{cx}" cy="{cy}" rx="{rx + 14}" ry="{ry + 14}" fill="none" stroke="{MEMB}" stroke-width="3"/>')
add(f'<ellipse cx="{cx}" cy="{cy}" rx="{rx}" ry="{ry}" fill="{CYTO}" stroke="{MEMB}" stroke-width="3"/>')
for k in range(3):
    add(f'<ellipse cx="{cx}" cy="{cy}" rx="{rx - 26 - 12 * k}" ry="{ry - 24 - 12 * k}" fill="none" stroke="{THY}" stroke-width="2.5" stroke-dasharray="120 40" stroke-dashoffset="{30 * k}"/>')
add(f'<text x="{cx + rx - 40}" y="{cy - ry + 60}" text-anchor="end" font-size="11" fill="{MUTED}">thylakoid membranes (photosynthesis)</text>')
add(f'<text x="{cx}" y="{cy + ry + 40}" text-anchor="middle" font-size="12" fill="{MUTED}">Prochlorococcus cell (schematic)</text>')
add(f'<text x="30" y="34" font-size="20" font-weight="700" fill="{INK}">2-Oxoglutarate: where fixed carbon meets nitrogen</text>')
add(f'<text x="30" y="56" font-size="12" fill="{INK2}">Numbers = genomes carrying the gene, of {N} Prochlorococcus genomes (KEGG annotations in the KG). Layout and arrows = interpretation.</text>')

# sun + CO2
add('<g>' + "".join(f'<line x1="{110 + 40 * __import__("math").cos(a)}" y1="{130 + 40 * __import__("math").sin(a)}" x2="{110 + 56 * __import__("math").cos(a)}" y2="{130 + 56 * __import__("math").sin(a)}" stroke="#d4a017" stroke-width="3"/>' for a in [i * 0.785 for i in range(8)])
    + f'<circle cx="110" cy="130" r="30" fill="#f2c94c"/></g>')
add(f'<text x="110" y="205" text-anchor="middle" font-size="12" fill="{INK2}">light</text>')


def node(x, y, label, col=INK, w=None, fill="#fff", bold=True, size=13):
    w = w or 8 * len(label) + 22
    add(f'<rect x="{x - w / 2}" y="{y - 15}" width="{w}" height="30" rx="15" fill="{fill}" stroke="{col}" stroke-width="2"/>'
        f'<text x="{x}" y="{y + 5}" text-anchor="middle" font-size="{size}" font-weight="{700 if bold else 400}" fill="{col if fill == "#fff" else "#fff"}">{esc(label)}</text>')


def arrow(x1, y1, x2, y2, col, mk, w=2.6, dash="", curve=None):
    d = f"M{x1},{y1} Q{curve[0]},{curve[1]} {x2},{y2}" if curve else f"M{x1},{y1} L{x2},{y2}"
    add(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{w}"{dash} marker-end="url(#{mk})"/>')


def enz(x, y, text, n, col=INK2, anchor="middle"):
    miss = n == 0
    t = f"{text} {n}/{N}"
    add(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="11.5" font-weight="600" fill="{MISS if miss else col}" '
        f'paint-order="stroke" stroke="{CYTO}" stroke-width="4">{"✕ " if miss else ""}{esc(t)}</text>')


# CO2 -> carboxysome
add(f'<text x="200" y="300" font-size="16" font-weight="700" fill="{C}">CO₂</text>')
arrow(235, 300, 300, 300, C, "ac")
hx, hy = 360, 300
pts = " ".join(f"{hx + 52 * __import__('math').cos(a + 0.5236):.0f},{hy + 52 * __import__('math').sin(a + 0.5236):.0f}" for a in [i * 1.0472 for i in range(6)])
add(f'<polygon points="{pts}" fill="#e3eedf" stroke="{C}" stroke-width="2"/>')
add(f'<text x="{hx}" y="{hy - 4}" text-anchor="middle" font-size="11" font-weight="700" fill="{C}">carboxysome</text>')
add(f'<text x="{hx}" y="{hy + 12}" text-anchor="middle" font-size="10.5" fill="{INK2}">RuBisCO</text>')
enz(hx, hy + 72, "rbcL/rbcS", route_n("rubisco"), C)

# Calvin -> glycolysis
node(360, 420, "3-Phosphoglycerate", C)
arrow(360, 355, 360, 402, C, "ac")
node(360, 530, "PEP", C)
arrow(360, 437, 360, 512, C, "ac")
enz(380, 480, "gpm, eno", route_n("lowglyc"), C, "start")
node(360, 650, "Pyruvate", C)
arrow(360, 547, 360, 632, C, "ac")
enz(340, 595, "pyk", route_n("pyk"), C, "end")
node(560, 420, "Acetyl-CoA", C)
arrow(385, 640, 540, 437, C, "ac", curve=(470, 600))
enz(505, 500, "pdh", route_n("pdh"), C, "start")

# open TCA
node(610, 560, "Oxaloacetate", C)
arrow(385, 535, 555, 560, C, "ac", curve=(470, 575))
enz(470, 590, "ppc", route_n("ppc"), C)
node(760, 470, "Citrate", C)
arrow(605, 425, 722, 462, C, "ac")
enz(655, 420, "gltA", route_n("gltA"), C, "start")
arrow(650, 545, 725, 486, C, "ac")
node(910, 470, "Isocitrate", C)
arrow(795, 470, 862, 470, C, "ac")
enz(828, 455, "acnB", route_n("acn"), C)
node(1000, 590, "2-Oxoglutarate", INK, w=150, fill=INK)
arrow(940, 488, 985, 570, C, "ac")
enz(955, 545, "icd", route_n("icd"), C, "end")
node(560, 680, "Malate", C)
arrow(568, 663, 600, 578, C, "ac")
enz(545, 625, "mqo", G[("mdh", "mdh/mqo")], C, "end")
node(660, 770, "Fumarate", C)
arrow(640, 755, 580, 697, C, "ac")
enz(590, 745, "fumC", route_n("fum"), C, "end")
# missing half
arrow(1000, 607, 975, 712, MISS, "am", 2, ' stroke-dasharray="6 5"')
add(f'<text x="965" y="730" text-anchor="middle" font-size="12" fill="{MISS}">Succinyl-CoA</text>')
arrow(925, 742, 860, 790, MISS, "am", 2, ' stroke-dasharray="6 5"')
add(f'<text x="830" y="805" text-anchor="middle" font-size="12" fill="{MISS}">Succinate</text>')
arrow(790, 798, 715, 778, MISS, "am", 2, ' stroke-dasharray="6 5"')
enz(1005, 680, "sucA/sucB", G[("ogdh", "sucA")], C, "start")
enz(915, 790, "sucC/sucD", G[("scs", "sucC")], C, "start")
enz(755, 818, "sdhA/B complete in", route_n("sdh"), MISS)
add(f'<text x="790" y="610" text-anchor="middle" font-size="13" font-weight="700" fill="{INK2}">TCA cycle — open</text>')
add(f'<text x="790" y="627" text-anchor="middle" font-size="11" fill="{INK2}">no gene for 2-oxoglutarate → succinyl-CoA</text>')
add(f'<text x="790" y="642" text-anchor="middle" font-size="11" fill="{INK2}">in any genome: carbon stops at 2-oxoglutarate</text>')

# nitrogen
add(f'<rect x="{cx + rx - 14}" y="300" width="28" height="44" rx="6" fill="#efe9f7" stroke="{NIT}" stroke-width="2"/>')
add(f'<text x="1490" y="290" text-anchor="end" font-size="16" font-weight="700" fill="{NIT}">NH₄⁺</text>')
add(f'<text x="1490" y="362" text-anchor="end" font-size="10.5" fill="{MUTED}">transporter</text>')
add(f'<text x="1490" y="376" text-anchor="end" font-size="10.5" fill="{MUTED}">(drawn, not counted)</text>')
arrow(cx + rx + 40, 322, cx + rx - 20, 322, NIT, "an")
node(1220, 470, "Glutamine", NIT)
node(1220, 640, "Glutamate", NIT, w=120, fill=NIT)
arrow(cx + rx - 20, 330, 1240, 452, NIT, "an", curve=(1300, 360))
arrow(1260, 625, 1255, 487, NIT, "an", curve=(1300, 560))
enz(1300, 555, "glnA (GS)", route_n("gln"), NIT, "start")
add(f'<text x="1300" y="572" font-size="10.5" fill="{INK2}">+ NH₄⁺</text>')
arrow(1180, 487, 1165, 625, NIT, "an", curve=(1130, 540))
arrow(1076, 596, 1155, 632, NIT, "an", curve=(1100, 640))
enz(1128, 515, "glsF (GOGAT)", route_n("gogat"), NIT, "end")
add(f'<text x="1105" y="700" text-anchor="middle" font-size="10.5" fill="{INK2}">Gln + 2-OG → 2 Glu</text>')
add(f'<text x="1110" y="430" text-anchor="middle" font-size="12" font-weight="700" fill="{NIT}">GS–GOGAT cycle</text>')
# glutamate family + transaminases
node(1300, 760, "Pro", NIT)
node(1300, 815, "Arg", NIT)
arrow(1250, 655, 1290, 742, NIT, "an")
arrow(1235, 657, 1285, 797, NIT, "an")
enz(1325, 760, f"proB/A/C", route_n("pro"), NIT, "start")
enz(1325, 815, f"arg genes", route_n("arg"), NIT, "start")
add(f'<rect x="830" y="820" width="320" height="58" rx="10" fill="#fff" stroke="{NIT}" stroke-width="1.5" stroke-dasharray="4 3"/>')
add(f'<text x="990" y="842" text-anchor="middle" font-size="12" font-weight="700" fill="{NIT}">Asp, Ala, Val, Leu, Ile, Phe, Tyr, His …</text>')
add(f'<text x="990" y="860" text-anchor="middle" font-size="10.5" fill="{INK2}">get their –NH₂ from glutamate via transaminases</text>')
add(f'<text x="990" y="873" text-anchor="middle" font-size="10.5" fill="{INK2}">aspC {G[("asp", "aspC")]}/{N} · ilvE {G[("val", "ilvE")]}/{N} · hisC {G[("his", "hisC")]}/{N}</text>')
arrow(1185, 655, 1080, 815, NIT, "an", curve=(1150, 760))
add(f'<text x="1155" y="760" font-size="10.5" font-weight="700" fill="{NIT}">–NH₂</text>')

# signal (interpretation)
add(f'<rect x="900" y="250" width="250" height="62" rx="10" fill="#fff" stroke="{INK2}" stroke-width="1.2" stroke-dasharray="3 3"/>')
add(f'<text x="1025" y="272" text-anchor="middle" font-size="12" font-weight="700" fill="{INK}">NtcA (nitrogen regulator)</text>')
add(f'<text x="1025" y="289" text-anchor="middle" font-size="10.5" fill="{INK2}">2-OG builds up when N is short →</text>')
add(f'<text x="1025" y="303" text-anchor="middle" font-size="10.5" fill="{INK2}">switches on N-uptake genes [interpretation]</text>')
arrow(1000, 572, 1015, 316, INK2, "ai", 1.5, ' stroke-dasharray="4 4"', curve=(1040, 450))
# KG note
add(f'<rect x="700" y="365" width="260" height="46" rx="8" fill="#fff" stroke="{RULE}"/>')
add(f'<text x="830" y="384" text-anchor="middle" font-size="10.5" fill="{INK}">N stress: acnB and icd upregulated in MED4</text>')
add(f'<text x="830" y="399" text-anchor="middle" font-size="10" fill="{MUTED}">Tolonen et al. 2006 (KG literature index)</text>')

# legend
lx, ly = 30, 900
for i, (col, txt, dash) in enumerate([(C, "carbon path", ""), (NIT, "nitrogen path", ""), (MISS, f"no gene in any of the {N} genomes", ' stroke-dasharray="6 5"')]):
    x = lx + i * 230
    add(f'<line x1="{x}" x2="{x + 40}" y1="{ly}" y2="{ly}" stroke="{col}" stroke-width="2.6"{dash}/><text x="{x + 48}" y="{ly + 4}" font-size="12" fill="{INK2}">{txt}</text>')
add(f'<text x="{lx}" y="{ly + 26}" font-size="11" fill="{INK2}">Counts for multi-gene steps = genomes carrying every gene. sdh is complete in only some genomes; whether a bypass closes the cycle is not checked [gap].</text>')
add('</svg>')
svg = "\n".join(o)
(F / "cell_carbon_nitrogen.svg").write_text(svg, encoding="utf-8")
(R / "cell_carbon_nitrogen.html").write_text(f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>2-Oxoglutarate: where carbon meets nitrogen</title></head>'
                                             f'<body style="margin:16px;background:#fcfcfb">{svg}</body></html>', encoding="utf-8")
print("wrote", F / "cell_carbon_nitrogen.svg")
