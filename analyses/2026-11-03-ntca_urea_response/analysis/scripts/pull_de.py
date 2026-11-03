"""Pull all-gene DE rows for the 4 experiments; add within-experiment/timepoint |log2fc| percentile."""
from pathlib import Path
import pandas as pd
import multiomics_explorer as m

OUT = Path(__file__).resolve().parents[1] / "data"
EXP = {
    ("MED4", "urea"): "10.1038/msb4100087_growth_medium_growth_on_urea_as_med4_microarray",
    ("MIT9313", "urea"): "10.1038/msb4100087_growth_medium_growth_on_urea_as_mit9313_microarray",
    ("MED4", "ndep"): "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray",
    ("MIT9313", "ndep"): "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray",
}
log = []
for (strain, arm), eid in EXP.items():
    res, off = [], 0
    while True:
        r = m.differential_expression_by_gene(organism=strain, experiment_ids=[eid], limit=5000, offset=off)
        res += r["results"]; off += r["returned"]
        if r["returned"] == 0 or off >= r["total_matching"]:
            break
    df = pd.DataFrame(res)
    assert len(df) == r["total_matching"], (len(df), r["total_matching"])
    df["abs_log2fc"] = df["log2fc"].abs()
    df["abs_pct"] = df.groupby("timepoint", dropna=False)["abs_log2fc"].rank(pct=True) * 100
    df["strain"] = strain
    df.to_csv(OUT / f"de_all_{strain}_{arm}.csv", index=False, encoding="utf-8")
    log.append(dict(tool="differential_expression_by_gene", experiment_id=eid, organism=strain,
                    total_matching=r["total_matching"], saved=len(df),
                    n_genes=df["locus_tag"].nunique(), timepoints=";".join(map(str, df["timepoint"].unique())),
                    n_significant=int(df["expression_status"].str.startswith("significant").sum())))
L = pd.DataFrame(log)
L.to_csv(OUT / "pull_de_kg_calls.csv", index=False, encoding="utf-8")
print(L.to_string())
