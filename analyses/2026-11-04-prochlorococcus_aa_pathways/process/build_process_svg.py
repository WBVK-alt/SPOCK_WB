"""One-image overview of how the 2026-11-04 Prochlorococcus amino-acid analysis was built, stage by stage.

Content is a hand-written record of the session (questions, tools, outputs, results, corrections), with commit ids from
`git log`. Numbers quoted are the ones recorded in the analysis files and QC log. Writes process/process_overview.svg.
"""
from pathlib import Path
import html
import textwrap

OUT = Path(__file__).resolve().parent / "process_overview.svg"
esc = html.escape
INK, INK2, MUTED, RULE, SURF, CARD = "#1f1f1d", "#52514e", "#8a8984", "#d6d5d0", "#fcfcfb", "#ffffff"
PHASE = {"Scope": "#2f7d4f", "Map the genomes": "#61409b", "Expression": "#7b74b5", "Interpret": "#355070", "Check": "#c9a227"}
WARN = "#b5651d"

STAGES = [
    ("Scope", "Your question", "All KG data on Prochlorococcus amino-acid biosynthesis and degradation: genes per pathway, and under which conditions they go up or down.",
     "Loaded the research-methodology rules · kg_release_info", "—", "KG 0.1.0-alpha.9, compatibility ok", None, ""),
    ("Scope", "What does the KG hold?", "Which strains and experiments exist?",
     "list_organisms · list_experiments", "—", "45 Prochlorococcus genomes · 140 experiments (42 genome-only)", None, ""),
    ("Map the genomes", "Gene inventory", "How many genes per pathway?",
     "genes_by_ontology: 14 KEGG amino-acid pathways + 40 GO 'biosynthetic / catabolic process' terms",
     "grounding/data/counts_*.csv", "Median 2 (Trp metabolism) to 24 genes per pathway; GO catabolic terms nearly empty",
     "PAC1 skipped: its name also matches EQPAC1", "c929334"),
    ("Expression", "First expression pull", "Under which conditions are these genes up / down?",
     "differential_expression_by_gene × 10 strains · list_experiments", "expression/reports/summary.html",
     "8,241 rows, 67 experiments; every pathway changes under most conditions; whole-genome change 1–88%",
     "32 of 67 experiments are partial tables", "c929334"),
    ("Map the genomes", "Pathway master map", "Draw the pathways and how they meet the TCA cycle, with strain counts.",
     "Step list from textbook biochemistry (59 arrows, 120 genes as KEGG KOs) · KO presence via API (44) + Cypher (PAC1), cross-checked identical (4,541 rows) · genes_by_function name search in MED4",
     "grounding/reports/pathway_map.html", "TCA open: sucA/B 0/45, sucC/D 0/45 · His 0/45 only because hisN is unannotated · degradation mostly 0/45",
     "Alternative enzymes found (mqo, aroH, metA/B, psp); PSTA has a wrong KO → counted by ortholog group", "e815a24"),
    ("Interpret", "Where is degradation?", "Are these only biosynthetic pathways?",
     "Genes behind the KEGG 'degradation' pathways", "chat", "Mostly shared enzymes (ALDH, lpd, transaminases); real ones few: hdc 43/44, DAO 42/44",
     "My map missed hdc, DAO, leucine DH, mmsB", ""),
    ("Interpret", "Literature check", "Do papers discuss the open TCA cycle or the His gap?",
     "list_publications · discussed_by_publication over all 60 papers", "chat",
     "No paper discusses either; closest: Tolonen 2006 — acnB/icd 'synthesis of 2-oxoglutarate', up under N stress", None, ""),
    ("Map the genomes", "Cluster map", "Which strains have which genes, and which pathways are incomplete?",
     "151 gene rows (path order) × 45 strains, Jaccard + average-linkage clustering", "grounding/reports/cluster_map.html",
     "18 of 20 amino acids complete in 41–45 genomes; gaps sit in 6 strains (HNLC2 11, HNLC1 7, MIT9313 4 …); Asn 20/45",
     "HNLC genomes may be incomplete assemblies — not checked", "743d50d"),
    ("Expression", "Expression overview", "Which genes change, under which conditions?",
     "Fresh DE pull for every map locus tag · 101 genes × 67 experiments, any timepoint, any paralog",
     "expression/reports/de_overview.html", "Nitrogen and darkness columns densest; iron / viral / phosphorus near-blank (partial tables)", None, "ae28444"),
    ("Expression", "Expression on the pathway map", "Three maps: nitrogen, extended light/darkness, mitomycin.",
     "13 experiments chosen in chat · one square per gene per experiment", "expression/reports/de_maps.html",
     "N: protein leans up (56/19), coculture RNA down (14/35) · NATL2A darkness leans down", "Experiment lists are my choice — to confirm", "2cd527c"),
    ("Interpret", "Direction map", "Which routes can run forward, backward, or both?",
     "Per-step reversibility from textbook biochemistry (not KG)", "grounding/reports/direction_map.html",
     "5 of 24 biosynthesis arrows run both ways (Asp, Ala, Val, Gly, homoserine); 12 of 18 exits have no gene", "Reversibility not from the KG", "ef0998a"),
    ("Interpret", "Summary and the 2-oxoglutarate idea", "Patterns and gaps; explain 'carbon meets nitrogen at 2-OG'.",
     "Synthesis of all files · cell illustration", "grounding/reports/cell_carbon_nitrogen.html",
     "Build everything, break down almost nothing; open TCA ends at 2-OG, which GS-GOGAT uses to take up NH₄⁺", None, "d72b053"),
    ("Check", "QC of the expression maps", "Show the calculation, step by step.",
     "Standalone script, fresh KG queries, independent code; compares every square with the saved file",
     "expression/qc/qc_de_maps.py · qc_log.txt", "1,560 / 1,560 squares identical — QC passed",
     "Corrections: mitomycin 39↑/42↓ (not 41/47, double count) · darkness lean is NATL2A only; extended light also leans down · constant-dark has 94/97 genes not in table · 'both' was left out · icd mixed (3↑ 2↓)",
     "8ab1fb8"),
    ("Check", "Evidence check of the 2-OG explanation", "Is each claim actually in the KG?",
     "metabolites_by_gene (reactions) · genes_by_metabolite (L-glutamate) · assays_by_metabolite (2-OG) · DE incl. ntcA",
     "grounding/reports/evidence_2og.html", "Reactions verified (glnA R00253, GOGAT R00021, icd R00267); 14/20 amino acids get N via Glu/Gln/NH₄⁺; glnA and ntcA up in 6/7 N experiments",
     "Not in KG: carbon flux to 2-OG · 2-OG build-up under N stress · 2-OG → NtcA signal", "552472d"),
]

W, X0, CW, NX = 1500, 150, 860, 1040
lines = []
add = lines.append


def wrap(s, n):
    return textwrap.wrap(s, n) or [""]


# measure
y = 150
geom = []
for st in STAGES:
    ph, title, q, how, out, res, note, sha = st
    body = len(wrap(q, 110)) + len(wrap(how, 118)) + len(wrap(res, 112)) + 1
    h = 46 + 16 * body + 10
    nh = 30 + 15 * len(wrap(note, 46)) if note else 0
    geom.append((y, max(h, nh)))
    y += max(h, nh) + 18
H = y + 190

add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="system-ui,Segoe UI,sans-serif">')
add(f'<rect width="{W}" height="{H}" fill="{SURF}"/>')
add(f'<text x="40" y="50" font-size="24" font-weight="700" fill="{INK}">How we got here: Prochlorococcus amino-acid pathways, 2026-11-04</text>')
add(f'<text x="40" y="76" font-size="13" fill="{INK2}">Each card = one step of the session: the question that started it, how it was done (KG tools), what it produced, and what it showed.</text>')
add(f'<text x="40" y="96" font-size="13" fill="{INK2}">Side notes = gaps, choices and corrections found along the way. Commit ids = local git history of analyses/2026-11-04-prochlorococcus_aa_pathways.</text>')
lx = 40
for i, (ph, c) in enumerate(PHASE.items()):
    add(f'<rect x="{lx}" y="114" width="14" height="14" rx="3" fill="{c}"/><text x="{lx + 20}" y="126" font-size="12.5" fill="{INK}">{ph}</text>')
    lx += 30 + 8 * len(ph) + 20
add(f'<line x1="{X0 - 60}" x2="{X0 - 60}" y1="{geom[0][0]}" y2="{geom[-1][0] + geom[-1][1]}" stroke="{RULE}" stroke-width="3"/>')

for k, (st, (yy, h)) in enumerate(zip(STAGES, geom), 1):
    ph, title, q, how, out, res, note, sha = st
    c = PHASE[ph]
    add(f'<circle cx="{X0 - 60}" cy="{yy + 22}" r="17" fill="{c}"/><text x="{X0 - 60}" y="{yy + 27}" text-anchor="middle" font-size="14" font-weight="700" fill="#fff">{k}</text>')
    add(f'<rect x="{X0 - 20}" y="{yy}" width="{CW}" height="{h}" rx="10" fill="{CARD}" stroke="{RULE}"/>')
    add(f'<rect x="{X0 - 20}" y="{yy}" width="6" height="{h}" rx="3" fill="{c}"/>')
    add(f'<text x="{X0}" y="{yy + 24}" font-size="15.5" font-weight="700" fill="{INK}">{esc(title)}</text>')
    add(f'<text x="{X0 + CW - 40}" y="{yy + 24}" text-anchor="end" font-size="11" fill="{MUTED}">{esc(ph)}{"  ·  " + sha if sha else ""}</text>')
    ty = yy + 44
    for s in wrap(q, 110):
        add(f'<text x="{X0}" y="{ty}" font-size="12.5" font-style="italic" fill="{INK2}">{esc(s)}</text>'); ty += 16
    for i, s in enumerate(wrap(how, 118)):
        add(f'<text x="{X0}" y="{ty}" font-size="11.5" fill="{INK2}">{"How: " if i == 0 else ""}{esc(s)}</text>'); ty += 16
    for i, s in enumerate(wrap(res, 112)):
        add(f'<text x="{X0}" y="{ty}" font-size="12.5" font-weight="600" fill="{INK}">{"→ " if i == 0 else "   "}{esc(s)}</text>'); ty += 16
    if out != "—" and out != "chat":
        add(f'<text x="{X0}" y="{ty}" font-size="11" fill="{c}" font-family="Consolas,monospace">{esc(out)}</text>')
    elif out == "chat":
        add(f'<text x="{X0}" y="{ty}" font-size="11" fill="{MUTED}">(answered in chat)</text>')
    if note:
        nl = wrap(note, 46)
        nh = 26 + 15 * len(nl)
        warn = ph == "Check" and note.startswith("Corrections")
        col = WARN if warn or "wrong" in note or "missed" in note else MUTED
        add(f'<line x1="{X0 - 20 + CW}" x2="{NX}" y1="{yy + 22}" y2="{yy + 22}" stroke="{RULE}" stroke-dasharray="3 3"/>')
        add(f'<rect x="{NX}" y="{yy}" width="400" height="{nh}" rx="8" fill="#fff" stroke="{col}" stroke-width="1.4" stroke-dasharray="{"" if col == WARN else "4 3"}"/>')
        add(f'<text x="{NX + 12}" y="{yy + 18}" font-size="11" font-weight="700" fill="{col}">{"⚠ correction / issue" if col == WARN else "note · gap · choice"}</text>')
        for i, s in enumerate(nl):
            add(f'<text x="{NX + 12}" y="{yy + 34 + 15 * i}" font-size="11.5" fill="{INK}">{esc(s)}</text>')

fy = geom[-1][0] + geom[-1][1] + 40
add(f'<rect x="40" y="{fy}" width="{W - 80}" height="130" rx="12" fill="#fff" stroke="{INK}" stroke-width="1.5"/>')
add(f'<text x="60" y="{fy + 28}" font-size="15" font-weight="700" fill="{INK}">Where this leaves us</text>')
for i, s in enumerate([
    "Supported by the KG: biosynthesis genes for 18/20 amino acids in 41–45 genomes; open TCA as annotated; GS-GOGAT and transaminase reactions; glnA and ntcA induced under N stress.",
    "Not yet in hand: carbon flux, 2-OG under N stress, hisN identity, a TCA bypass, HNLC genome completeness, any test against the genome-wide background.",
    "Status: Plan-phase data collection only — no proposal.md, no statistics, no critic review yet. Next: confirm experiment sets, then enrichment per condition / bypass and hisN search.",
]):
    add(f'<text x="60" y="{fy + 54 + 22 * i}" font-size="12.5" fill="{INK2}">• {esc(s)}</text>')
add('</svg>')
OUT.write_text("\n".join(lines), encoding="utf-8")
print("wrote", OUT, "height", H)
