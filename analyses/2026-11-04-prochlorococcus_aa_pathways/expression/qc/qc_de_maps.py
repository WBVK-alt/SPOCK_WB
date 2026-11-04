"""QC for expression/reports/de_maps.html — recompute every square on the three condition maps from scratch.

Run from the repo root (the explorer reads KG credentials from .env there):
    uv run python -X utf8 analyses/2026-11-04-prochlorococcus_aa_pathways/expression/qc/qc_de_maps.py

What it does (each step prints what it found):
  0. KG identity.
  1. The three experiment sets, checked to exist in the KG, with treatment / control / table scope.
  2. The gene list on the map (gene label -> KEGG KO alternatives), from pathway_steps.py — the only input shared
     with the original pipeline (it IS the definition of the map).
  3. Locus tags per gene in each strain, re-queried from the KG (genes_by_ontology on the KOs; PSTA via its
     Cyanorak ortholog group, because its KO in the KG is a hydrogenase).
  4. Differential-expression rows for those locus tags in those experiments, re-queried from the KG.
  5. One square per (gene, experiment), by this rule over all paralogs and all timepoints:
        no locus tag in the strain                     -> no_gene        (dot)
        locus tag(s), but no row in the experiment     -> not_measured   (white box)
        any row significant_up AND any significant_down-> both           (split)
        any row significant_up                         -> up             (orange)
        any row significant_down                       -> down           (blue)
        rows exist, none significant                   -> ns             (grey)
     'significant_up/down' is the paper's own call as stored in the KG (expression_status); no new statistics.
  6. Worked example: every raw row behind one gene in one condition set.
  7. Tallies: per experiment (each gene once) and per map (each gene slot on the map, so a gene drawn on several
     arrows — e.g. ilvE — counts once per arrow; this is what build_de_maps.py printed).
  8. Comparison with the saved file the maps were drawn from (expression/data/de_overview_cells.csv).
Outputs: expression/qc/qc_cells.csv, expression/qc/qc_mismatches.csv.
"""
from pathlib import Path
import sys
import pandas as pd
import multiomics_explorer as m

HERE = Path(__file__).resolve().parent
A = HERE.parents[1]                                   # analysis root
sys.path.insert(0, str(A / "grounding" / "scripts"))
from pathway_steps import ROUTES  # noqa: E402

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80); pd.set_option("display.max_rows", 400)
hr = lambda t: print(f"\n{'=' * 100}\n{t}\n{'=' * 100}")

# ---------------------------------------------------------------- 0
hr("STEP 0 · KG identity")
k = m.kg_release_info()
print("KG", k["kg"]["version"], "| explorer", k["explorer_version"], "| verdict", k["verdict"])

# ---------------------------------------------------------------- 1
hr("STEP 1 · experiment sets")
P = "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_"
SETS = {
    "nitrogen": [P + "rnaseq_axenic", P + "rnaseq_coculture", P + "proteomics_axenic", P + "proteomics_coculture",
                 "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq",
                 "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray",
                 "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray"],
    "lightdark": ["10.64898/2026.04.15.718746_extended_darkness_med4",
                  "10.1128/mSystems.00040-18_darkness_extended_darkness_natl2a_rnaseq_axenic",
                  "10.1128/mSystems.00040-18_darkness_extended_darkness_natl2a_rnaseq_coculture",
                  "10.1371/journal.pone.0165375_light_stress_constant_dark_med4_rnaseq_dark",
                  "10.64898/2026.04.15.718746_extended_light_med4"],
    "mitomycin": ["10.1016/j.cell.2022.12.006_mmc_dna_damage_mit0604"],
}
all_ids = [e for v in SETS.values() for e in v]
ex = m.list_experiments(experiment_ids=all_ids, verbose=True, limit=None)
assert not ex["not_found"], ex["not_found"]
E = {e["experiment_id"]: e for e in ex["results"]}
for s, ids in SETS.items():
    print(f"\n[{s}]")
    for i, e in enumerate(ids, 1):
        x = E[e]
        print(f"  [{i}] {x['organism_name'].replace('Prochlorococcus ', '')} | {x['treatment']} vs {x['control']} | "
              f"{x['omics_type']} | table={x['table_scope']} | timepoints={len(x.get('timepoints') or []) or 1}")

# ---------------------------------------------------------------- 2
hr("STEP 2 · genes on the map (label -> KO alternatives)")
GENES = {}
for r in ROUTES:
    for label, alts in r["genes"]:
        if label in GENES:
            assert GENES[label] == list(alts), f"label {label} has two different KO sets"
        GENES[label] = list(alts)
print(f"{len(GENES)} distinct genes on {len(ROUTES)} arrows; examples:")
for g in ["glnA", "glsF/gltB", "icd", "mdh/mqo", "serC/PSTA"]:
    print(f"  {g:10s} -> {GENES[g]}")

# ---------------------------------------------------------------- 3
hr("STEP 3 · locus tags per gene, per strain (fresh KG query)")
orgs = sorted({E[e]["organism_name"] for e in all_ids})
kos = sorted({k for v in GENES.values() for k in v if k.startswith("kegg.")})
ogs = sorted({k for v in GENES.values() for k in v if not k.startswith("kegg.")})
tags = {}                                  # (organism, KO or OG) -> set(locus_tag)
for o in orgs:
    r = m.genes_by_ontology(ontology="kegg", organism=o, term_ids=kos, min_gene_set_size=1, max_gene_set_size=None, limit=None)
    assert r["returned"] == r["total_matching"]
    for x in r["results"]:
        tags.setdefault((o, x["term_id"]), set()).add(x["locus_tag"])
    print(f"  {o}: {r['total_matching']} (gene, KO) rows")
for og in ogs:
    h = m.genes_by_homolog_group(group_ids=[og], limit=None)
    for x in h["results"]:
        if x["organism_name"] in orgs:
            tags.setdefault((x["organism_name"], og), set()).add(x["locus_tag"])
    print(f"  ortholog group {og}: {sum(1 for x in h['results'] if x['organism_name'] in orgs)} genes in these strains")
gene_tags = {(o, g): set().union(*[tags.get((o, k), set()) for k in alts]) for o in orgs for g, alts in GENES.items()}
print("\nexample, MED4:")
for g in ["glnA", "glsF/gltB", "icd", "mdh/mqo", "serC/PSTA", "sucA", "asnB"]:
    print(f"  {g:10s} -> {sorted(gene_tags[('Prochlorococcus MED4', g)]) or 'none (no gene)'}")

# ---------------------------------------------------------------- 4
hr("STEP 4 · DE rows (fresh KG query)")
rows = []
for o in orgs:
    ids = [e for e in all_ids if E[e]["organism_name"] == o]
    t = sorted(set().union(*[gene_tags[(o, g)] for g in GENES]))
    d = m.differential_expression_by_gene(organism=o, locus_tags=t, experiment_ids=ids, limit=None)
    assert d["returned"] == d["total_matching"]
    rows += d["results"]
    print(f"  {o}: {len(t)} locus tags, {len(ids)} experiments -> {d['total_matching']} rows")
de = pd.DataFrame(rows)
print("\nrows per experiment x status:")
print(de.groupby(["experiment_id", "expression_status"]).size().unstack(fill_value=0).rename(index=lambda s: s[-60:]).to_string())

# ---------------------------------------------------------------- 5
hr("STEP 5 · one square per (gene, experiment)")
def square(g, e):
    o = E[e]["organism_name"]
    t = gene_tags[(o, g)]
    if not t:
        return "no_gene", ""
    sub = de[(de.experiment_id == e) & de.locus_tag.isin(t)]
    if sub.empty:
        return "not_measured", ""
    up, dn = (sub.expression_status == "significant_up").any(), (sub.expression_status == "significant_down").any()
    st = "both" if up and dn else "up" if up else "down" if dn else "ns"
    return st, "; ".join(f"{x.locus_tag} {x.timepoint}: {x.expression_status} ({x.log2fc:+.2f})" for x in sub.itertuples())

cells = pd.DataFrame([dict(set=s, exp_no=i, experiment_id=e, gene=g, **dict(zip(("status", "rows"), square(g, e))))
                      for s, ids in SETS.items() for i, e in enumerate(ids, 1) for g in GENES])
cells.to_csv(HERE / "qc_cells.csv", index=False, encoding="utf-8")
print(f"{len(cells)} squares written to qc_cells.csv")

# ---------------------------------------------------------------- 6
hr("STEP 6 · worked example: glnA, nitrogen set")
for _, c in cells[(cells.gene == "glnA") & (cells.set == "nitrogen")].iterrows():
    print(f"  [{c.exp_no}] -> {c.status.upper():12s} from: {c.rows or '(no rows)'}")

# ---------------------------------------------------------------- 7
hr("STEP 7 · tallies")
print("per experiment, each gene counted once:")
per = cells.groupby(["set", "exp_no"]).status.value_counts().unstack(fill_value=0)
print(per.reindex(columns=["up", "down", "both", "ns", "not_measured", "no_gene"], fill_value=0).to_string())
slots = [g for r in ROUTES for g, _ in r["genes"]]
print(f"\nper map, each gene SLOT counted ({len(slots)} slots; genes on several arrows count once per arrow):")
for s in SETS:
    c = cells[cells.set == s].set_index(["experiment_id", "gene"]).status
    tally = pd.Series([c[(e, g)] for e in SETS[s] for g in slots]).value_counts()
    print(f"  {s:10s}", {k: int(tally.get(k, 0)) for k in ("up", "down", "both", "ns", "not_measured")})

# ---------------------------------------------------------------- 8
hr("STEP 8 · compare with the saved file the maps were drawn from")
saved = pd.read_csv(A / "expression" / "data" / "de_overview_cells.csv", encoding="utf-8")
saved = saved[saved.experiment_id.isin(all_ids)][["experiment_id", "gene", "status"]]
cmp_ = cells.merge(saved, on=["experiment_id", "gene"], how="left", suffixes=("_qc", "_saved"))
# genes absent from all 10 expression strains were left out of the saved file; the map draws them as dots
cmp_["status_saved"] = cmp_.status_saved.fillna("no_gene")
mis = cmp_[cmp_.status_qc != cmp_.status_saved]
mis.to_csv(HERE / "qc_mismatches.csv", index=False, encoding="utf-8")
print(f"squares compared: {len(cmp_)} | identical: {len(cmp_) - len(mis)} | mismatches: {len(mis)}")
if len(mis):
    print(mis[["set", "exp_no", "gene", "status_qc", "status_saved"]].to_string())
print("\nQC", "PASSED" if mis.empty else "FOUND DIFFERENCES — see qc_mismatches.csv")
