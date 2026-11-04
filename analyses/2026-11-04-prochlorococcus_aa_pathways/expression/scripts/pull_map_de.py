"""DE rows for every locus tag on the pathway map (grounding/data/ko_presence.csv), in each strain with expression data.
Significance = each paper's own call as stored in the KG."""
from pathlib import Path
import pandas as pd
import multiomics_explorer as m

HERE = Path(__file__).resolve().parents[1]
G, D = HERE.parent / "grounding" / "data", HERE / "data"
pres = pd.read_csv(G / "ko_presence.csv", encoding="utf-8")
orgs = sorted(pd.read_csv(D / "experiments.csv", encoding="utf-8").organism_name.unique())
rows, calls = [], []
for o in orgs:
    tags = sorted(pres[pres.organism == o].locus_tag.unique())
    de = m.differential_expression_by_gene(organism=o, locus_tags=tags, limit=None)
    assert de["returned"] == de["total_matching"], o
    calls.append(dict(organism=o, n_tags=len(tags), rows=de["total_matching"], genes=de["matching_genes"], no_expression=len(de["no_expression"])))
    rows += [dict(organism=o, **r) for r in de["results"]]
df = pd.DataFrame(rows)
df["treatment_type"] = df.treatment_type.apply(lambda v: ";".join(v) if isinstance(v, list) else v)
df.to_csv(D / "map_de_rows.csv", index=False, encoding="utf-8")
pd.DataFrame(calls).to_csv(D / "kg_calls_map_de.csv", index=False, encoding="utf-8")
print(pd.DataFrame(calls).to_string()); print(df.expression_status.value_counts())
