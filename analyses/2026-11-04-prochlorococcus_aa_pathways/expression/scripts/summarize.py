"""Count tables from expression/data/de_rows.csv + gene_set_membership.csv. Counts only, no statistics.

pathway_by_experiment.csv : set x experiment -> genes in set, genes with a row, genes up / down at >=1 timepoint
pathway_by_timepoint.csv  : same, per timepoint
pathway_condition_overview.csv : per set, which treatment types / experiments have >=1 gene up or down
"""
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
D = HERE / "data"
rows = pd.read_csv(D / "de_rows.csv", encoding="utf-8")
memb = pd.read_csv(D / "gene_set_membership.csv", encoding="utf-8")
exp = pd.read_csv(D / "experiments.csv", encoding="utf-8")

j = rows.merge(memb[["organism", "locus_tag", "set_kind", "set_name"]], on=["organism", "locus_tag"])
set_size = memb.groupby(["organism", "set_kind", "set_name"]).locus_tag.nunique().rename("n_genes_in_set")

def agg(g):
    up = sorted(g[g.expression_status == "significant_up"].locus_tag.unique())
    dn = sorted(g[g.expression_status == "significant_down"].locus_tag.unique())
    return pd.Series(dict(n_with_row=g.locus_tag.nunique(), n_up=len(up), n_down=len(dn),
                          up_genes=";".join(up), down_genes=";".join(dn)))

keys = ["organism", "set_kind", "set_name", "experiment_id"]
pe = j.groupby(keys).apply(agg, include_groups=False).reset_index()
pe = pe.merge(set_size.reset_index(), on=["organism", "set_kind", "set_name"])
pe = pe.merge(exp[["experiment_id", "treatment_type", "treatment", "control", "omics_type", "table_scope", "publication_doi"]], on="experiment_id", how="left")
pe.to_csv(D / "pathway_by_experiment.csv", index=False, encoding="utf-8")

pt = j.groupby(keys + ["timepoint", "timepoint_order"], dropna=False).apply(agg, include_groups=False).reset_index()
pt.to_csv(D / "pathway_by_timepoint.csv", index=False, encoding="utf-8")

ov = []
for (kind, name), g in pe.groupby(["set_kind", "set_name"]):
    u, d = g[g.n_up > 0], g[g.n_down > 0]
    ov.append(dict(set_kind=kind, set_name=name, n_experiments_with_rows=g.experiment_id.nunique(),
                   n_experiments_any_up=u.experiment_id.nunique(), n_experiments_any_down=d.experiment_id.nunique(),
                   treatment_types_up=";".join(sorted(set(";".join(u.treatment_type.astype(str)).split(";")) - {""})),
                   treatment_types_down=";".join(sorted(set(";".join(d.treatment_type.astype(str)).split(";")) - {""}))))
pd.DataFrame(ov).to_csv(D / "pathway_condition_overview.csv", index=False, encoding="utf-8")

# gene-level: per gene x experiment, status per timepoint
gl = rows.copy()
gl["call"] = gl.timepoint.fillna("").astype(str) + ":" + gl.expression_status.str.replace("significant_", "").str.replace("not_significant", "ns") \
    + "(" + gl.log2fc.round(2).astype(str) + ")"
gene = gl.groupby(["organism", "locus_tag", "gene_name", "experiment_id"], dropna=False).agg(
    n_up=("expression_status", lambda s: (s == "significant_up").sum()),
    n_down=("expression_status", lambda s: (s == "significant_down").sum()),
    calls=("call", " | ".join)).reset_index()
sets = memb.groupby(["organism", "locus_tag"]).set_name.apply(lambda s: "; ".join(sorted(set(s)))).rename("sets")
gene = gene.merge(sets.reset_index(), on=["organism", "locus_tag"], how="left")
gene.to_csv(D / "gene_by_experiment.csv", index=False, encoding="utf-8")

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 90)
print(len(pe), len(pt), len(gene))
print(pd.DataFrame(ov)[lambda x: x.set_kind == "KEGG pathway"].to_string())

# KEGG pathways, experiments with full source tables only: # experiments with >=1 gene up / down, by treatment type
full = pe[(pe.set_kind == "KEGG pathway") & (pe.table_scope == "all_detected_genes")]
tt = lambda d: "; ".join(f"{k} {v}" for k, v in d.treatment_type.value_counts().items())
fo = pd.DataFrame([dict(pathway=p, n_experiments=g.experiment_id.nunique(),
                        n_exp_any_up=int((g.n_up > 0).sum()), n_exp_any_down=int((g.n_down > 0).sum()),
                        up_by_treatment_type=tt(g[g.n_up > 0]), down_by_treatment_type=tt(g[g.n_down > 0]))
                   for p, g in full.groupby("set_name")])
fo.to_csv(D / "kegg_pathway_conditions_full_tables.csv", index=False, encoding="utf-8")
print(fo.to_string())
