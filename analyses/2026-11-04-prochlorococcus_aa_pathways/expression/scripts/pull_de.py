"""Every KG differential-expression row for Prochlorococcus amino-acid pathway genes.

Gene set = grounding/data (KEGG amino-acid pathways + GO aa biosynthetic/catabolic terms), per genome.
Significance = each paper's own call as stored in the KG (expression_status). No new statistics.
Experiment context (treatment, control, table_scope, omics) from list_experiments.
"""
from pathlib import Path
import ast
import pandas as pd
import multiomics_explorer as m

HERE = Path(__file__).resolve().parents[1]
G = HERE.parent / "grounding" / "data"
OUT = HERE / "data"
calls = []

kegg = pd.read_csv(G / "kegg_pathway_genes.csv", encoding="utf-8")
go = pd.read_csv(G / "go_aa_genes.csv", encoding="utf-8")
memb = pd.concat([
    kegg.assign(set_kind="KEGG pathway", set_name=kegg.pathway)[["organism", "locus_tag", "gene_name", "product", "set_kind", "set_name"]],
    go.assign(set_kind="GO process", set_name=go.amino_acid + " " + go.direction)[["organism", "locus_tag", "gene_name", "product", "set_kind", "set_name"]],
]).drop_duplicates(["organism", "locus_tag", "set_kind", "set_name"])
memb.to_csv(OUT / "gene_set_membership.csv", index=False, encoding="utf-8")

gen = pd.read_csv(G / "genomes.csv", encoding="utf-8")
gen["omics"] = gen.omics_types.apply(ast.literal_eval)
with_expr = gen[gen.omics.apply(lambda o: any(x not in ("BIOINFORMATIC",) for x in o))].organism_name.tolist()

rows, nomeas = [], []
for g in with_expr:
    tags = sorted(memb[memb.organism == g].locus_tag.unique())
    if not tags:
        calls.append(dict(tool="differential_expression_by_gene", organism=g, note="no aa genes (genome skipped in grounding)")); continue
    de = m.differential_expression_by_gene(organism=g, locus_tags=tags, limit=None)
    assert de["returned"] == de["total_matching"], (g, de["returned"], de["total_matching"])
    calls.append(dict(tool="differential_expression_by_gene", organism=g, n_locus_tags=len(tags),
                      total_matching=de["total_matching"], matching_genes=de["matching_genes"],
                      experiment_count=de["experiment_count"], no_expression=len(de["no_expression"])))
    rows += [dict(organism=g, **r) for r in de["results"]]
    nomeas += [dict(organism=g, locus_tag=t) for t in de["no_expression"]]
rows = pd.DataFrame(rows)
rows["treatment_type"] = rows.treatment_type.apply(lambda v: ";".join(v) if isinstance(v, list) else v)
rows.to_csv(OUT / "de_rows.csv", index=False, encoding="utf-8")
pd.DataFrame(nomeas).to_csv(OUT / "not_measured.csv", index=False, encoding="utf-8")

# experiment context
ex = m.list_experiments(experiment_ids=sorted(rows.experiment_id.unique()), verbose=True, limit=None)
calls.append(dict(tool="list_experiments[verbose]", total_matching=ex["total_matching"], returned=ex["returned"]))
keep = ["experiment_id", "experiment_name", "organism_name", "publication_doi", "treatment_type", "treatment", "control",
        "omics_type", "table_scope", "control_timing", "coculture_partner", "light_condition", "experimental_context",
        "is_time_course", "gene_count", "design_notes"]
exp = pd.DataFrame([{k: (";".join(map(str, r[k])) if isinstance(r.get(k), list) else r.get(k)) for k in keep} for r in ex["results"]])
exp.to_csv(OUT / "experiments.csv", index=False, encoding="utf-8")
pd.DataFrame(calls).to_csv(OUT / "kg_calls.csv", index=False, encoding="utf-8")

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 70)
print(pd.DataFrame(calls).to_string())
print(rows.expression_status.value_counts())
print(exp[["experiment_id", "treatment_type", "table_scope", "omics_type"]].to_string())
print("missing keys:", [k for k in keep if k not in ex["results"][0]])
