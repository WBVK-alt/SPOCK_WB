"""Genome-wide significant up/down counts per experiment x timepoint (KG rollup), as context for the pathway counts."""
from pathlib import Path
import pandas as pd
import multiomics_explorer as m

D = Path(__file__).resolve().parents[1] / "data"
ids = sorted(pd.read_csv(D / "experiments.csv", encoding="utf-8").experiment_id)
ex = m.list_experiments(experiment_ids=ids, verbose=True, limit=None)
assert ex["returned"] == len(ids)
out = []
for e in ex["results"]:
    for t in (e.get("timepoints") or [{}]):
        s = t.get("genes_by_status") or e.get("genes_by_status") or {}
        out.append(dict(experiment_id=e["experiment_id"], timepoint=t.get("timepoint"), timepoint_order=t.get("timepoint_order"),
                        genes_measured=t.get("gene_count", e.get("gene_count")), significant_up=s.get("significant_up"),
                        significant_down=s.get("significant_down"), not_significant=s.get("not_significant")))
pd.DataFrame(out).to_csv(D / "genome_wide_status.csv", index=False, encoding="utf-8")
print(pd.DataFrame(out).head(10).to_string()); print(len(out))
