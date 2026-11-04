"""Report: Prochlorococcus amino-acid pathway genes — counts per genome and up/down calls per experiment.

Reads only grounding/data and expression/data CSVs. Writes expression/reports/summary.html.
"""
from pathlib import Path
import json
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
G, D, R = HERE.parent / "grounding" / "data", HERE / "data", HERE / "reports"
R.mkdir(exist_ok=True)
rd = lambda p: pd.read_csv(p, encoding="utf-8")

gen = rd(G / "genomes.csv")
kc = rd(G / "counts_kegg_pathway_by_genome.csv").set_index("pathway")
gc = rd(G / "counts_go_aa_by_genome.csv").set_index("row")
skipped = rd(G / "skipped_genomes.csv")
pe = rd(D / "pathway_by_experiment.csv")
exp = rd(D / "experiments.csv")
gw = rd(D / "genome_wide_status.csv")
gene = rd(D / "gene_by_experiment.csv")
memb = rd(D / "gene_set_membership.csv")
calls = rd(D / "kg_calls.csv")

# genomes ordered by clade then strain
cl = gen.set_index("organism_name").clade
cols = sorted(kc.columns, key=lambda o: (str(cl.get(o)), o))
short = lambda o: o.replace("Prochlorococcus marinus subsp. marinus CCMP1375 (SS120)", "SS120").replace("Prochlorococcus ", "")
genomes = [dict(name=short(o), clade=str(cl.get(o)), expr=o in set(pe.organism)) for o in cols]
kegg_mat = dict(rows=kc.index.tolist(), vals=kc[cols].values.tolist())
n_go_rows = len(gc)
gc = gc[(gc[cols].sum(axis=1) > 0)]
go_mat = dict(rows=gc.index.tolist(), vals=gc[cols].values.tolist(), dropped=int(n_go_rows - len(gc)))

# experiments: genome-wide context (max fraction significant over timepoints)
gw["frac"] = (gw.significant_up + gw.significant_down) / gw.genes_measured
gwx = gw.groupby("experiment_id").agg(gw_max=("frac", "max"), n_tp=("timepoint", lambda s: s.notna().sum()))
exp = exp.merge(gwx, left_on="experiment_id", right_index=True, how="left")
exp["org"] = exp.organism_name.map(short)
exp = exp.sort_values(["org", "treatment_type", "omics_type", "experiment_name"])
names = memb.drop_duplicates(["organism", "locus_tag"]).set_index(["organism", "locus_tag"]).gene_name


def lab(o, t):
    n = names.get((o, t))
    return f"{n} ({t})" if isinstance(n, str) else t


def glist(o, s):
    return ", ".join(lab(o, t) for t in s.split(";")) if isinstance(s, str) and s else ""


cells = {}
for r in pe.itertuples():
    cells.setdefault(r.set_kind, {})[f"{r.experiment_id}||{r.set_name}"] = [
        int(r.n_up), int(r.n_down), int(r.n_with_row), int(r.n_genes_in_set),
        glist(r.organism, r.up_genes), glist(r.organism, r.down_genes)]
exps = [dict(id=e.experiment_id, org=e.org, tt=e.treatment_type, name=e.experiment_name, tr=str(e.treatment), ctl=str(e.control),
             om=e.omics_type, scope=e.table_scope, gw=None if (pd.isna(e.gw_max) or e.table_scope != "all_detected_genes") else round(100 * e.gw_max, 1),
             ntp=int(e.n_tp) if pd.notna(e.n_tp) else 0, doi=e.publication_doi) for e in exp.itertuples()]

gene["exp_name"] = gene.experiment_id.map(exp.set_index("experiment_id").experiment_name)
gene["tt"] = gene.experiment_id.map(exp.set_index("experiment_id").treatment_type)
gtab = [[short(r.organism), r.locus_tag, r.gene_name if isinstance(r.gene_name, str) else "", r.sets, r.tt, r.exp_name,
         int(r.n_up), int(r.n_down), r.calls] for r in gene.itertuples()]

de_calls = calls[calls.tool == "differential_expression_by_gene"]
tiles = dict(n_genomes=len(cols), n_pro=int(len(gen)), n_kegg=len(kc), n_expr_org=int(pe.organism.nunique()),
             n_exp=int(exp.experiment_id.nunique()), n_rows=int(de_calls.total_matching.sum()),
             skipped=", ".join(skipped.organism.map(short)) or "none",
             partial=int(exp.table_scope.ne("all_detected_genes").sum()))
payload = dict(genomes=genomes, kegg=kegg_mat, go=go_mat, exps=exps, cells=cells, go_sets=gc.index.tolist(), gtab=gtab, tiles=tiles)

HTML = (HERE / "scripts" / "report_template.html").read_text(encoding="utf-8")
(R / "summary.html").write_text(HTML.replace("__DATA__", json.dumps(payload, default=str)), encoding="utf-8")
print("wrote", R / "summary.html", (R / "summary.html").stat().st_size // 1024, "KB")
