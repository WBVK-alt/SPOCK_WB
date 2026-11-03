"""Q3: are Alteromonas amino-acid genes switched on / changed in coculture with Prochlorococcus?

Only KG experiment comparing Alteromonas in coculture with Prochlorococcus vs axenic:
  HOT1A3 + MED4 vs axenic HOT1A3, RNA-seq, days 11/18/31 (Weissberg et al. 2025 preprint).
Gene set = HOT1A3 genes in grounding/data (KEGG amino-acid pathways + GO aa biosynthetic/catabolic terms).
Significance = the paper's own call as stored in the KG (expression_status). No new statistics, except
the package's pathway_enrichment (ORA) over the 14 KEGG amino-acid pathways.
"""
from pathlib import Path
import pandas as pd
import multiomics_explorer as m

HERE = Path(__file__).resolve().parents[1]
G = HERE.parent / "grounding" / "data"
OUT = HERE / "data"
ORG = "Alteromonas macleodii HOT1A3"
EXP = "10.1101/2025.11.24.690089_coculture_prochlorococcus_med4_hot1a3_rnaseq"
calls = []

kegg = pd.read_csv(G / "kegg_pathway_genes.csv", encoding="utf-8"); kegg = kegg[kegg.organism == ORG]
go = pd.read_csv(G / "go_aa_genes.csv", encoding="utf-8"); go = go[go.organism == ORG]
memb = pd.concat([
    kegg.assign(set_kind="KEGG pathway", set_name=kegg.pathway)[["locus_tag", "gene_name", "product", "set_kind", "set_name"]],
    go.assign(set_kind="GO process", set_name=go.amino_acid + " " + go.direction,
              amino_acid=go.amino_acid, direction=go.direction)[["locus_tag", "gene_name", "product", "set_kind", "set_name", "amino_acid", "direction"]],
]).drop_duplicates(["locus_tag", "set_kind", "set_name"])
memb.to_csv(OUT / "gene_set_membership.csv", index=False, encoding="utf-8")
tags = sorted(memb.locus_tag.unique())

de = m.differential_expression_by_gene(organism=ORG, experiment_ids=[EXP], locus_tags=tags, limit=None)
assert de["returned"] == de["total_matching"], (de["returned"], de["total_matching"])
calls.append(dict(tool="differential_expression_by_gene", experiment=EXP, n_locus_tags=len(tags),
                  total_matching=de["total_matching"], matching_genes=de["matching_genes"],
                  no_expression=len(de["no_expression"]), not_found=len(de["not_found"])))
rows = pd.DataFrame(de["results"])
rows.to_csv(OUT / "de_rows.csv", index=False, encoding="utf-8")
pd.Series(de["no_expression"], name="locus_tag").to_csv(OUT / "not_measured.csv", index=False, encoding="utf-8")

# genome-wide status counts per timepoint, for reference
ex = m.list_experiments(experiment_ids=[EXP], verbose=True)["results"][0]
gw = pd.DataFrame([dict(timepoint=t["timepoint"], timepoint_order=t["timepoint_order"], **t["genes_by_status"],
                        genes_measured=t["gene_count"]) for t in ex["timepoints"]])
gw.to_csv(OUT / "genome_wide_status.csv", index=False, encoding="utf-8")
calls.append(dict(tool="list_experiments", experiment=EXP, total_matching=1))

# per set x timepoint summary (paper's call)
tps = gw.sort_values("timepoint_order").timepoint.tolist()
summ = []
for (kind, name), sub in memb.groupby(["set_kind", "set_name"]):
    st = set(sub.locus_tag)
    for tp in tps:
        r = rows[(rows.timepoint == tp) & rows.locus_tag.isin(st)]
        summ.append(dict(set_kind=kind, set_name=name, timepoint=tp, n_genes=len(st), n_measured=r.locus_tag.nunique(),
                         n_up=(r.expression_status == "significant_up").sum(),
                         n_down=(r.expression_status == "significant_down").sum(),
                         n_not_sig=(r.expression_status == "not_significant").sum(),
                         up_genes=";".join(sorted(r[r.expression_status == "significant_up"].locus_tag)),
                         down_genes=";".join(sorted(r[r.expression_status == "significant_down"].locus_tag))))
summ = pd.DataFrame(summ)
summ.to_csv(OUT / "set_by_timepoint_summary.csv", index=False, encoding="utf-8")

# package ORA over the 14 KEGG amino-acid pathways (background = genes the experiment measured)
pw = sorted(kegg.pathway_id.unique())
enr = m.pathway_enrichment(organism=ORG, experiment_ids=[EXP], ontology="kegg", level=2, term_ids=pw,
                           direction="both", min_gene_set_size=5, include_nonsignificant=True)
env = enr.to_envelope(limit=None)
calls.append(dict(tool="pathway_enrichment[kegg level2, 14 aa pathways]", experiment=EXP,
                  total_matching=env.get("total_matching"), n_significant=env.get("n_significant")))
ora = m.to_dataframe(env)
ora.to_csv(OUT / "ora_kegg_aa_pathways.csv", index=False, encoding="utf-8")
pd.DataFrame(calls).to_csv(OUT / "kg_calls.csv", index=False, encoding="utf-8")

pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 60)
print(pd.DataFrame(calls).to_string()); print(gw.to_string())
print(rows.expression_status.value_counts()); print(ora.columns.tolist()); print(ora.head(60).to_string())
