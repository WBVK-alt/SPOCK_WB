"""One 16:9 slide (1920x1080 SVG): chronology of the 2026-11-04 session and the rationale behind each step.
Writes process/process_slide.svg.
"""
from pathlib import Path
import html
import textwrap

OUT = Path(__file__).resolve().parent / "process_slide.svg"
esc = html.escape
INK, INK2, MUTED, RULE, SURF = "#1f1f1d", "#52514e", "#8a8984", "#d6d5d0", "#fcfcfb"
C = {"scope": "#2f7d4f", "genome": "#61409b", "expr": "#7b74b5", "interp": "#355070", "check": "#c9a227"}

STEPS = [  # phase, title, what we did, why
    ("scope", "The question", "Amino-acid synthesis and breakdown genes in Prochlorococcus, and when they change",
     "Fix the scope before touching data"),
    ("scope", "Ground in the KG", "Checked KG version; listed 45 genomes and 140 experiments",
     "The KG is the only data source — know what it holds first"),
    ("genome", "Gene inventory", "Genes per pathway: 14 KEGG pathways + GO terms, every genome",
     "Two annotation systems: KEGG mixes synthesis and breakdown, GO separates them"),
    ("expr", "Expression pull", "All up/down calls for these genes: 67 experiments, 10 strains",
     "Answer 'under which conditions' — with table type and genome-wide change as context"),
    ("genome", "Pathway map", "Each step tied to its enzyme; genomes with every gene counted",
     "See how pathways meet the TCA cycle; a 0 was checked by name search before trusting it"),
    ("interp", "Breakdown & papers", "Looked behind the KEGG 'degradation' sets; searched 60 papers",
     "Answer your questions from the KG, and flag what it cannot answer"),
    ("genome", "Cluster map", "Gene presence, 20 amino acids × 45 strains, clustered",
     "Find which strains have incomplete pathways and whether gaps cluster"),
    ("expr", "Expression maps", "Overview (101 genes × 67) and three condition maps",
     "Put the calls on the pathways; one square per experiment so opposite treatments are not pooled"),
    ("interp", "Direction map", "Each route marked forward, both ways, or exit",
     "Your question; textbook biochemistry, labelled as not from the KG"),
    ("interp", "Carbon + nitrogen", "Summary and a cell diagram of 2-oxoglutarate",
     "Explain the pattern in plain terms before building on it"),
    ("check", "QC", "Recomputed every map square from fresh KG queries with new code",
     "Show the calculation; numbers must reproduce independently"),
    ("check", "Evidence check", "Each claim tested against KG reactions, expression and metabolite data",
     "Separate what the KG shows from interpretation and literature"),
]
RULES = ["The KG is the only data source", "Every number comes from a script and a saved file",
         "Locus tags, not gene names", "Check a zero before believing it", "Tag KG data vs interpretation",
         "Verify before claiming"]

W, H = 1920, 1080
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="system-ui,Segoe UI,sans-serif">',
     f'<rect width="{W}" height="{H}" fill="{SURF}"/>']
add = o.append
add(f'<text x="80" y="88" font-size="44" font-weight="700" fill="{INK}">How we worked: Prochlorococcus amino-acid pathways</text>')
add(f'<text x="80" y="128" font-size="22" fill="{INK2}">Twelve steps, in order — what we did, and why we did it next</text>')

# timeline: 6 cards above, 6 below a central line
n = len(STEPS)
cw, chh = 252, 262
lx0, lx1, ly = 80 + cw / 2, W - 80 - cw / 2, 560
add(f'<line x1="{lx0 - 30}" x2="{lx1 + 30}" y1="{ly}" y2="{ly}" stroke="{INK2}" stroke-width="4" stroke-linecap="round"/>')
add(f'<path d="M{lx1 + 30},{ly - 12} L{lx1 + 52},{ly} L{lx1 + 30},{ly + 12} z" fill="{INK2}"/>')
step_x = (lx1 - lx0) / (n - 1)
for i, (ph, title, what, why) in enumerate(STEPS):
    x = lx0 + i * step_x
    col = C[ph]
    up = i % 2 == 0
    add(f'<circle cx="{x}" cy="{ly}" r="22" fill="{col}" stroke="{SURF}" stroke-width="4"/>')
    add(f'<text x="{x}" y="{ly + 7}" text-anchor="middle" font-size="19" font-weight="700" fill="#fff">{i + 1}</text>')
    cx = x - cw / 2
    cy = ly - 44 - chh if up else ly + 44
    add(f'<line x1="{x}" x2="{x}" y1="{ly - 22 if up else ly + 22}" y2="{cy + chh if up else cy}" stroke="{col}" stroke-width="2"/>')
    add(f'<rect x="{cx}" y="{cy}" width="{cw}" height="{chh}" rx="14" fill="#fff" stroke="{RULE}" stroke-width="1.5"/>')
    add(f'<rect x="{cx}" y="{cy}" width="{cw}" height="8" rx="4" fill="{col}"/>')
    add(f'<text x="{cx + 18}" y="{cy + 42}" font-size="19.5" font-weight="700" fill="{INK}">{esc(title)}</text>')
    ty = cy + 72
    for s in textwrap.wrap(what, 27):
        add(f'<text x="{cx + 18}" y="{ty}" font-size="15.5" fill="{INK}">{esc(s)}</text>'); ty += 20
    ty += 8
    add(f'<text x="{cx + 18}" y="{ty}" font-size="13" font-weight="700" fill="{col}" letter-spacing=".5">WHY</text>'); ty += 21
    for s in textwrap.wrap(why, 28):
        add(f'<text x="{cx + 18}" y="{ty}" font-size="15" font-style="italic" fill="{INK2}">{esc(s)}</text>'); ty += 19

# phase legend
lx = 80
for nm, k in (("Scope", "scope"), ("Map the genomes", "genome"), ("Expression", "expr"), ("Interpret", "interp"), ("Check", "check")):
    add(f'<circle cx="{lx + 9}" cy="{168}" r="9" fill="{C[k]}"/><text x="{lx + 26}" y="{174}" font-size="18" fill="{INK2}">{nm}</text>')
    lx += 60 + 10 * len(nm)

# working rules
ry = 1036
add(f'<text x="80" y="{ry - 46}" font-size="19" font-weight="700" fill="{INK}">Rules applied at every step</text>')
x = 80
for r in RULES:
    w = 7.9 * len(r) + 28
    add(f'<rect x="{x}" y="{ry - 26}" width="{w}" height="38" rx="19" fill="#fff" stroke="{INK2}" stroke-width="1.2"/>'
        f'<text x="{x + w / 2}" y="{ry}" text-anchor="middle" font-size="16" fill="{INK}">{esc(r)}</text>')
    x += w + 10
add('</svg>')
OUT.write_text("\n".join(o), encoding="utf-8")
print("wrote", OUT, "rules end x", x)
