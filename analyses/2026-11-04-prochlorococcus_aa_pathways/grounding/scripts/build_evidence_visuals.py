"""Visuals for the evidence check of '2-oxoglutarate links carbon and nitrogen'.

A. Evidence map — the claim chain, each arrow shaded by KG support (supported / partly / not in KG).
B. Heatmap — gltA, acnB, icd, glsF (GOGAT), glnA, ntcA across the 7 nitrogen experiments of de_maps (KG DE rows).
C. Nitrogen donors — per amino acid, the MED4 gene + KEGG reaction that moves N from glutamate / glutamine / NH4+,
   each reaction verified against the KG (metabolites_by_gene) before drawing.
Run from repo root. Writes grounding/data/evidence_heatmap.csv, evidence_n_donors.csv,
grounding/reports/evidence_2og.html, grounding/figures/evidence_{map,heatmap,donors}.svg.
Evidence grades in panel A summarise the chat assessment of 2026-11-04 [interpretation of KG results].
"""
from pathlib import Path
import html
import pandas as pd
import multiomics_explorer as m

HERE = Path(__file__).resolve().parent
D, R, F = HERE.parent / "data", HERE.parent / "reports", HERE.parent / "figures"
esc = html.escape
INK, INK2, MUTED, RULE, SURF = "#1f1f1d", "#52514e", "#8a8984", "#d6d5d0", "#fcfcfb"
UP, DN, MID = "#eb6834", "#2a78d6", "#f0efec"
E_YES, E_PART, E_NO = "#470f84", "#9e9ac8", "#b5b3ad"


def svg_open(w, h):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" font-family="system-ui,Segoe UI,sans-serif">',
            f'<rect width="{w}" height="{h}" fill="{SURF}"/>']

# =============================================================== A. evidence map
o = svg_open(1300, 640)
add = o.append
add('<defs>' + "".join(f'<marker id="m{k}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="{c}"/></marker>'
                        for k, c in (("y", E_YES), ("p", E_PART), ("n", E_NO))) + '</defs>')
add(f'<text x="24" y="32" font-size="18" font-weight="700" fill="{INK}">A · How much of the 2-oxoglutarate story is in the KG?</text>')
N = {"CO2": (110, 200, "CO₂ (photosynthesis)"), "OAA": (330, 200, "Citrate / Isocitrate"), "OG": (560, 200, "2-Oxoglutarate"),
     "GLU": (820, 200, "Glutamate"), "GLN": (820, 80, "Glutamine"), "NH4": (1100, 80, "NH₄⁺"),
     "AA": (1100, 200, "Asp, Phe, Tyr, Val, Leu,\nIle, His, Lys, Arg, Gly"), "PRO": (1100, 300, "Pro, Arg (glutamate family)"),
     "NTCA": (560, 420, "NtcA (N regulator)"), "NUP": (880, 420, "N-uptake genes (glnA, amt1)"), "NSTR": (220, 420, "Nitrogen shortage")}
def nbox(k, bold=False, dark=False):
    x, y, t = N[k]
    lines = t.split("\n"); w = max(len(s) for s in lines) * 7.6 + 26; h = 26 + 15 * (len(lines) - 1)
    add(f'<rect x="{x - w / 2}" y="{y - h / 2}" width="{w}" height="{h}" rx="12" fill="{INK if dark else "#fff"}" stroke="{INK2}" stroke-width="1.4"/>')
    for i, s in enumerate(lines):
        add(f'<text x="{x}" y="{y + 5 - 7.5 * (len(lines) - 1) + 15 * i}" text-anchor="middle" font-size="12.5" font-weight="{700 if bold or dark else 500}" fill="{"#fff" if dark else INK}">{esc(s)}</text>')
for k in N:
    nbox(k, dark=(k in ("OG", "GLU")))
def link(a, b, grade, label, dy=-10, curve=None, lx=None, ly=None):
    (x1, y1, _), (x2, y2, _) = N[a], N[b]
    col, mk, dash, w = {"y": (E_YES, "y", "", 3.2), "p": (E_PART, "p", ' stroke-dasharray="8 5"', 3), "n": (E_NO, "n", ' stroke-dasharray="2 5"', 2.6)}[grade]
    # shorten
    import math
    dx, dy_ = x2 - x1, y2 - y1; L = math.hypot(dx, dy_)
    sx, sy, ex, ey = x1 + dx / L * 85, y1 + dy_ / L * 22, x2 - dx / L * 85, y2 - dy_ / L * 22
    d = f"M{sx:.0f},{sy:.0f} Q{curve[0]},{curve[1]} {ex:.0f},{ey:.0f}" if curve else f"M{sx:.0f},{sy:.0f} L{ex:.0f},{ey:.0f}"
    add(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{w}" stroke-linecap="round"{dash} marker-end="url(#m{mk})"/>')
    tx, ty = lx or (sx + ex) / 2, ly or (sy + ey) / 2 + dy
    for i, s in enumerate(label.split("\n")):
        add(f'<text x="{tx}" y="{ty + 13 * i}" text-anchor="middle" font-size="10.5" fill="{INK2}" paint-order="stroke" stroke="{SURF}" stroke-width="4">{esc(s)}</text>')
link("CO2", "OAA", "n", "carbon actually flows this way?\nno flux data in KG", dy=-28)
link("OAA", "OG", "y", "icd: Isocitrate → 2-OG + CO₂\n(R00267; 45/45 genomes)", dy=-28)
link("OG", "GLU", "y", "GOGAT (glsF): Gln + 2-OG → 2 Glu\n(R00021; 45/45)", dy=-28)
link("NH4", "GLN", "y", "", dy=-10)
link("GLU", "GLN", "y", "glnA: Glu + NH₄⁺ → Gln\n(R00253; 45/45)", curve=(760, 140), lx=700, ly=128)
link("GLN", "GLU", "y", "", curve=(880, 140))
link("GLU", "AA", "p", "transaminases (aspC, ilvE, hisC, dapL…)\nmost amino acids — not Ala, Ser", dy=-30)
link("GLU", "PRO", "y", "proB, argJ, argD", dy=14)
link("NSTR", "NTCA", "y", "ntcA induced ≤22× (MED4)\n≤12× (MIT9313)", dy=-26)
link("NTCA", "NUP", "y", "glnA up in 6/7 N experiments", dy=-14)
link("OG", "NTCA", "n", "2-OG builds up and\nswitches on NtcA?\nnot in KG", curve=(470, 330), lx=430, ly=300)
add(f'<text x="690" y="252" text-anchor="middle" font-size="11" fill="{INK2}">acnB up in 4/7 N experiments (never down)</text>')
add(f'<text x="690" y="267" text-anchor="middle" font-size="11" fill="{INK2}">icd up 3/7, down 2/7 · Tolonen 2006: both up in MED4</text>')
ly = 520
for i, (c, dash, w, t) in enumerate([(E_YES, "", 3.2, "supported by KG data (gene + reaction, or expression)"),
                                      (E_PART, ' stroke-dasharray="8 5"', 3, "partly supported"),
                                      (E_NO, ' stroke-dasharray="2 5"', 2.6, "not in the KG — literature / interpretation")]):
    add(f'<line x1="24" x2="74" y1="{ly + 22 * i}" y2="{ly + 22 * i}" stroke="{c}" stroke-width="{w}"{dash} stroke-linecap="round"/>'
        f'<text x="84" y="{ly + 22 * i + 4}" font-size="12" fill="{INK2}">{t}</text>')
add(f'<text x="24" y="{ly + 90}" font-size="11" fill="{MUTED}">Reactions = KEGG reaction records for MED4 genes in the KG; genome counts = KEGG KO presence across 45 Prochlorococcus genomes.</text>')
add('</svg>')
svg_a = "\n".join(o)

# =============================================================== B. heatmap
P = "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_"
EXPS = [(P + "rnaseq_axenic", "MED4 low-N\nRNA axenic"), (P + "rnaseq_coculture", "MED4 low-N\nRNA coculture"),
        (P + "proteomics_axenic", "MED4 low-N\nprotein axenic"), (P + "proteomics_coculture", "MED4 low-N\nprotein cocult."),
        ("10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq", "MED4 N-depleted\nRNA ◆"),
        ("10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray", "MED4 N-depriv.\narray (Tolonen)"),
        ("10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray", "MIT9313 N-depriv.\narray (Tolonen)")]
pres = pd.read_csv(D / "ko_presence.csv", encoding="utf-8")
GK = [("gltA", "K01647", "citrate synthase"), ("acnB", "K01682", "aconitase"), ("icd", "K00031", "isocitrate DH → 2-OG"),
      ("glsF", "K00284", "GOGAT"), ("glnA", "K01915", "glutamine synthetase")]
NTCA = {"Prochlorococcus MED4": "PMM0246", "Prochlorococcus MIT9313": "PMT1831"}   # ntcA, as in 2026-11-03-ntca_urea_response
rows = []
for org in ("Prochlorococcus MED4", "Prochlorococcus MIT9313"):
    tag_of = {g: sorted(pres[(pres.organism == org) & (pres.ko_id == f"kegg.orthology:{k}")].locus_tag) for g, k, _ in GK}
    tag_of["ntcA"] = [NTCA[org]]
    tags = sorted({t for v in tag_of.values() for t in v})
    ids = [e for e, _ in EXPS if (org.split()[-1].lower() in e)]
    d = m.differential_expression_by_gene(organism=org, locus_tags=tags, experiment_ids=ids, limit=None)
    assert d["returned"] == d["total_matching"]
    df = pd.DataFrame(d["results"])
    assert set(df[df.locus_tag == NTCA[org]].gene_name.dropna()) <= {"ntcA"}, "ntcA locus tag check"
    for g, ts in tag_of.items():
        for e in ids:
            sub = df[(df.experiment_id == e) & df.locus_tag.isin(ts)]
            if sub.empty:
                rows.append(dict(gene=g, experiment_id=e, organism=org, status="not_measured", log2fc=None, timepoint="", n_tp=0, n_sig=0)); continue
            sig = sub[sub.expression_status != "not_significant"]
            pick = (sig if len(sig) else sub).loc[lambda s: s.log2fc.abs().idxmax()]
            up, dn = (sub.expression_status == "significant_up").any(), (sub.expression_status == "significant_down").any()
            rows.append(dict(gene=g, experiment_id=e, organism=org, status="both" if up and dn else "up" if up else "down" if dn else "ns",
                             log2fc=round(float(pick.log2fc), 2), timepoint=str(pick.timepoint), n_tp=len(sub), n_sig=len(sig)))
hm = pd.DataFrame(rows)
hm.to_csv(D / "evidence_heatmap.csv", index=False, encoding="utf-8")
genes = ["gltA", "acnB", "icd", "glsF", "glnA", "ntcA"]
glab = dict([(g, f"{g} · {t}") for g, _, t in GK] + [("ntcA", "ntcA · N regulator")])
cw, ch, left, top = 118, 46, 230, 120
o = svg_open(left + cw * len(EXPS) + 40, top + ch * len(genes) + 120)
add = o.append
add(f'<text x="24" y="32" font-size="18" font-weight="700" fill="{INK}">B · Under nitrogen stress: the 2-OG route and its nitrogen partners</text>')
add(f'<text x="24" y="52" font-size="12" fill="{INK2}">Number = log2 fold change (stress vs control) at the timepoint with the largest significant change (or largest change if none significant).</text>')
add(f'<text x="24" y="68" font-size="12" fill="{INK2}">Colour = paper\'s call at any timepoint; "k/n" = significant timepoints of n. Grey = measured, never significant.</text>')
for j, (e, lab) in enumerate(EXPS):
    for i, s in enumerate(lab.split("\n")):
        add(f'<text x="{left + j * cw + cw / 2}" y="{top - 26 + 13 * i}" text-anchor="middle" font-size="11" fill="{INK}">{esc(s)}</text>')
for i, g in enumerate(genes):
    y = top + i * ch
    add(f'<text x="{left - 10}" y="{y + ch / 2 + 4}" text-anchor="end" font-size="12.5" font-weight="{700 if g in ("acnB", "icd") else 500}" fill="{INK}">{esc(glab[g])}</text>')
    for j, (e, _) in enumerate(EXPS):
        r = hm[(hm.gene == g) & (hm.experiment_id == e)].iloc[0]
        x = left + j * cw
        fill = {"up": UP, "down": DN, "both": None, "ns": "#e6e5e1", "not_measured": "#fff"}[r.status]
        if r.status == "both":
            add(f'<path d="M{x + 2},{y + 2} h{cw - 4} L{x + 2},{y + ch - 2} z" fill="{UP}"/><path d="M{x + cw - 2},{y + 2} v{ch - 4} h{-(cw - 4)} z" fill="{DN}"/>')
        else:
            add(f'<rect x="{x + 2}" y="{y + 2}" width="{cw - 4}" height="{ch - 4}" rx="4" fill="{fill}" stroke="{RULE if r.status == "not_measured" else "none"}"/>')
        ink = "#fff" if r.status in ("up", "down", "both") else INK2
        txt = "not in table" if r.status == "not_measured" else f"{r.log2fc:+.2f}"
        add(f'<text x="{x + cw / 2}" y="{y + ch / 2 + 1}" text-anchor="middle" font-size="13" font-weight="700" fill="{ink}">{txt}</text>')
        if r.status != "not_measured":
            add(f'<text x="{x + cw / 2}" y="{y + ch / 2 + 15}" text-anchor="middle" font-size="9.5" fill="{ink}">{r.n_sig}/{r.n_tp} tp · {esc(r.timepoint)}</text>')
yl = top + ch * len(genes) + 30
for k, (c, t) in enumerate([(UP, "up"), (DN, "down"), ("both", "up + down"), ("#e6e5e1", "not significant")]):
    x = 24 + k * 150
    if c == "both":
        add(f'<path d="M{x},{yl} h16 L{x},{yl + 12} z" fill="{UP}"/><path d="M{x + 16},{yl} v12 h-16 z" fill="{DN}"/>')
    else:
        add(f'<rect x="{x}" y="{yl}" width="16" height="12" rx="2" fill="{c}"/>')
    add(f'<text x="{x + 22}" y="{yl + 10}" font-size="12" fill="{INK2}">{t}</text>')
add(f'<text x="24" y="{yl + 34}" font-size="11" fill="{MUTED}">◆ = filtered source table. MIT9313 genes: PMT locus tags of the same KOs; ntcA = PMM0246 / PMT1831.</text>')
add('</svg>')
svg_b = "\n".join(o)

# =============================================================== C. nitrogen donors
DONORS = [  # amino acid, route class, MED4 gene, locus tag, KEGG reaction expected for that gene
    ("Glu", "NH₄⁺ + 2-OG (GS-GOGAT)", "glsF", "PMM1512", "kegg.reaction:R00021"),
    ("Gln", "NH₄⁺ on glutamate", "glnA", "PMM0920", "kegg.reaction:R00253"),
    ("Asp", "from glutamate (transaminase)", "aspC", "PMM0674", "kegg.reaction:R00355"),
    ("Phe", "from glutamate (transaminase)", "aspC", "PMM0674", "kegg.reaction:R00694"),
    ("Tyr", "from glutamate (transaminase)", "aspC", "PMM0674", "kegg.reaction:R00734"),
    ("Val", "from glutamate (transaminase)", "ilvE", "PMM0878", "kegg.reaction:R01214"),
    ("Leu", "from glutamate (transaminase)", "ilvE", "PMM0878", "kegg.reaction:R01090"),
    ("Ile", "from glutamate (transaminase)", "ilvE", "PMM0878", "kegg.reaction:R02199"),
    ("His", "from glutamate (transaminase)", "hisC", "PMM0186", "kegg.reaction:R03243"),
    ("Lys", "from glutamate (transaminase)", "dapL", "PMM1500", "kegg.reaction:R07613"),
    ("Arg", "from glutamate (transaminase)", "argD", "PMM1301", "kegg.reaction:R02283"),
    ("Gly", "from glutamate (transaminase)", "spt", "PMM0919", "kegg.reaction:R00372"),
    ("Pro", "glutamate is the backbone", "proB", "PMM0788", "kegg.reaction:R00239"),
    ("Trp", "from glutamine (amidotransferase)", "trpE", "PMM1577", "kegg.reaction:R00986"),
    ("Ala", "NH₄⁺ direct (ald) — not via Glu", "ald", None, None),
    ("Asn", "from Asp (+ Gln); asnB in 20/45", "asnB", None, None),
    ("Ser", "transaminase PSTA — KO wrong in KG", "PSTA", "PMM0035", None),
    ("Thr", "via Asp (N already from Glu)", "—", None, None),
    ("Met", "via Asp (N already from Glu)", "—", None, None),
    ("Cys", "via Ser", "—", None, None),
]
tags = sorted({t for *_, t, r in DONORS if t and r})
mb = m.metabolites_by_gene(organism="Prochlorococcus MED4", locus_tags=tags, evidence_sources=["metabolism"], limit=None)
have = {(x["locus_tag"], x["reaction_id"]) for x in mb["results"]}
dn = []
for aa, cls, g, t, rxn in DONORS:
    ok = (t, rxn) in have if rxn else None
    if rxn:
        assert ok, f"{aa}: {g} {t} lacks {rxn} in KG"
    dn.append(dict(amino_acid=aa, route=cls, gene=g, locus_tag=t or "", kegg_reaction=rxn or "", verified_in_kg=ok))
dn = pd.DataFrame(dn)
dn.to_csv(D / "evidence_n_donors.csv", index=False, encoding="utf-8")
CLS = [("NH₄⁺", ["NH₄⁺ + 2-OG (GS-GOGAT)", "NH₄⁺ on glutamate"], "#470f84"),
       ("from glutamate", ["from glutamate (transaminase)", "glutamate is the backbone"], "#61409b"),
       ("from glutamine", ["from glutamine (amidotransferase)"], "#7b74b5"),
       ("other / not shown in KG", None, "#c9c6e3")]
def cls_of(route):
    for k, (nm, rs, c) in enumerate(CLS):
        if rs is None or route in rs:
            return k
o = svg_open(1200, 760)
add = o.append
add(f'<text x="24" y="32" font-size="18" font-weight="700" fill="{INK}">C · Where each amino acid gets its nitrogen (MED4, KEGG reactions in the KG)</text>')
add(f'<text x="24" y="52" font-size="12" fill="{INK2}">Filled = gene + reaction verified in the KG by this script; outlined = route inferred or not represented in the KG.</text>')
colx = [60, 340, 620, 900]
for k, (nm, _, c) in enumerate(CLS):
    add(f'<text x="{colx[k]}" y="92" font-size="13" font-weight="700" fill="{INK}">{esc(nm)}</text>')
    sub = dn[dn.route.map(cls_of) == k]
    for i, r in enumerate(sub.itertuples()):
        y = 110 + i * 46
        filled = bool(r.verified_in_kg)
        add(f'<rect x="{colx[k]}" y="{y}" width="250" height="38" rx="8" fill="{c if filled else "#fff"}" stroke="{c}" stroke-width="1.6"/>')
        tc = "#fff" if filled else INK
        add(f'<text x="{colx[k] + 12}" y="{y + 17}" font-size="13" font-weight="700" fill="{tc}">{esc(r.amino_acid)}</text>')
        add(f'<text x="{colx[k] + 52}" y="{y + 17}" font-size="11" fill="{tc}">{esc(r.gene)}{" (" + r.locus_tag + ")" if r.locus_tag else ""}</text>')
        add(f'<text x="{colx[k] + 12}" y="{y + 31}" font-size="10" fill="{tc}">{esc((r.kegg_reaction.replace("kegg.reaction:", "") + " · ") if r.kegg_reaction else "")}{esc(r.route)}</text>')
n_glu = int(dn.route.map(cls_of).isin([0, 1, 2]).sum())
add(f'<text x="24" y="730" font-size="12.5" fill="{INK}">{n_glu} of 20 amino acids: nitrogen arrives via GS-GOGAT, glutamate or glutamine, verified as KEGG reactions of MED4 genes. '
    f'Thr, Met, Cys, Asn get it second-hand; Ala and Ser are the exceptions.</text>')
add('</svg>')
svg_c = "\n".join(o)

for nm, s in (("evidence_map", svg_a), ("evidence_heatmap", svg_b), ("evidence_donors", svg_c)):
    (F / f"{nm}.svg").write_text(s, encoding="utf-8")
(R / "evidence_2og.html").write_text(
    f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>2-Oxoglutarate links carbon and nitrogen — what the KG supports</title></head>'
    f'<body style="margin:16px 24px;background:{SURF};font-family:system-ui">{svg_a}<br>{svg_b}<br>{svg_c}</body></html>', encoding="utf-8")
pd.set_option("display.width", 200)
print(hm.pivot(index="gene", columns="experiment_id", values="status").rename(columns=lambda c: c[-30:]).to_string())
print(dn.to_string())
