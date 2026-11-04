"""One 16:9 slide (1920x1080 SVG) summarising the 2026-11-04 session: process -> findings -> corrections & gaps.
Condensed from process_overview.svg (build_process_svg.py); numbers as recorded in the analysis files and QC log.
Writes process/process_slide.svg.
"""
from pathlib import Path
import html

OUT = Path(__file__).resolve().parent / "process_slide.svg"
esc = html.escape
INK, INK2, MUTED, RULE, SURF = "#1f1f1d", "#52514e", "#8a8984", "#d6d5d0", "#fcfcfb"
WARN = "#b5651d"
PHASES = [
    ("1 · Scope", "#2f7d4f", ["45 genomes, 140 experiments", "KG 0.1.0-alpha.9"]),
    ("2 · Map the genomes", "#61409b", ["14 KEGG pathways + GO terms", "Pathway map: 120 genes", "Cluster map: 45 strains"]),
    ("3 · Expression", "#7b74b5", ["67 experiments, 10 strains", "Overview: 101 genes × 67", "3 condition maps"]),
    ("4 · Interpret", "#355070", ["Degradation · literature", "Direction map", "2-OG cell diagram"]),
    ("5 · Check", "#c9a227", ["QC: 1,560/1,560 squares match", "Evidence check: 2-OG claim", "vs KG reactions"]),
]
FINDINGS = [
    ("Builds everything", "18 of 20 amino acids: complete biosynthesis in 41–45 of 45 genomes", "Asn 20/45 · His 0/45 (hisN unannotated)"),
    ("Breaks down almost nothing", "12 of 18 degradation routes: no gene in any genome", "Only 5 of 24 biosynthesis routes run both ways"),
    ("TCA cycle open at 2-OG", "No sucA/B or sucC/D in any genome; 2-OG feeds GS-GOGAT", "14 of 20 amino acids take N via Glu / Gln / NH₄⁺"),
    ("Nitrogen stress", "glnA and ntcA up in 6 of 7 N experiments", "acnB up 4/7 · icd mixed (3 up, 2 down)"),
]
CORR = ["Mitomycin: 39 up / 42 down (was double-counted as 41/47)",
        "Darkness lean-down is NATL2A only; extended light also leans down",
        "Constant-dark experiment: 94 of 97 genes not in its table — drop it",
        "PSTA has a wrong KO in the KG; PAC1 needs a direct query"]
GAPS = ["No carbon-flux data · no 2-OG measurement under N stress",
        "2-OG → NtcA signal: literature, not KG",
        "hisN gene and a possible TCA bypass: not searched yet",
        "No statistics vs genome-wide background yet"]

W, H = 1920, 1080
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="system-ui,Segoe UI,sans-serif">',
     f'<rect width="{W}" height="{H}" fill="{SURF}"/>']
add = o.append
add(f'<text x="80" y="92" font-size="44" font-weight="700" fill="{INK}">Prochlorococcus amino-acid pathways: how we got here</text>')
add(f'<text x="80" y="134" font-size="22" fill="{INK2}">From one question to a checked set of maps — 45 genomes, 67 experiments, all numbers from the multiomics KG</text>')

# process chevrons
x0, y0, pw, ph, gap = 80, 180, 336, 200, 20
for i, (t, c, items) in enumerate(PHASES):
    x = x0 + i * (pw + gap)
    tip = 26 if i < len(PHASES) - 1 else 0
    notch = 26 if i > 0 else 0
    add(f'<path d="M{x},{y0} h{pw - tip} l{tip},{ph / 2} l{-tip},{ph / 2} h{-(pw - tip)} l{notch},{-ph / 2} z" fill="{c}"/>')
    add(f'<text x="{x + notch + 22}" y="{y0 + 46}" font-size="26" font-weight="700" fill="#fff">{esc(t)}</text>')
    for k, s in enumerate(items):
        add(f'<text x="{x + notch + 22}" y="{y0 + 92 + 32 * k}" font-size="20" fill="#fff">{esc(s)}</text>')

# findings
fy = 430
add(f'<text x="80" y="{fy}" font-size="28" font-weight="700" fill="{INK}">What the KG shows</text>')
fw = (W - 160 - 3 * 24) / 4
for i, (t, a, b) in enumerate(FINDINGS):
    x = 80 + i * (fw + 24)
    add(f'<rect x="{x}" y="{fy + 22}" width="{fw}" height="228" rx="16" fill="#fff" stroke="{RULE}" stroke-width="1.5"/>')
    add(f'<rect x="{x}" y="{fy + 22}" width="{fw}" height="8" rx="4" fill="#61409b"/>')
    add(f'<text x="{x + 24}" y="{fy + 74}" font-size="26" font-weight="700" fill="{INK}">{esc(t)}</text>')
    for j, s in enumerate([a, b]):
        words, line, ln = s.split(), "", []
        for w_ in words:
            if len(line) + len(w_) > 34:
                ln.append(line); line = ""
            line = (line + " " + w_).strip()
        ln.append(line)
        yy = fy + 116 + j * 74
        for k, s2 in enumerate(ln):
            add(f'<text x="{x + 24}" y="{yy + 26 * k}" font-size="{20 if j == 0 else 18}" font-weight="{600 if j == 0 else 400}" fill="{INK if j == 0 else INK2}">{esc(s2)}</text>')

# corrections + gaps
cy = 740
cw = (W - 160 - 24) / 2
for i, (title, items, col, dash) in enumerate([("⚠ Corrected after QC", CORR, WARN, ""), ("Gaps — not in the KG / not done", GAPS, MUTED, ' stroke-dasharray="8 6"')]):
    x = 80 + i * (cw + 24)
    add(f'<rect x="{x}" y="{cy}" width="{cw}" height="228" rx="16" fill="#fff" stroke="{col}" stroke-width="2"{dash}/>')
    add(f'<text x="{x + 28}" y="{cy + 46}" font-size="26" font-weight="700" fill="{col}">{esc(title)}</text>')
    for k, s in enumerate(items):
        add(f'<text x="{x + 28}" y="{cy + 92 + 38 * k}" font-size="20" fill="{INK}">• {esc(s)}</text>')

add(f'<text x="80" y="1036" font-size="18" fill="{MUTED}">Status: Plan-phase data collection — no proposal, statistics or critic review yet. '
    f'Files: analyses/2026-11-04-prochlorococcus_aa_pathways · 10 local commits.</text>')
add('</svg>')
OUT.write_text("\n".join(o), encoding="utf-8")
print("wrote", OUT)
