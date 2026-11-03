"""Resolve ntcA + target panel and ribosomal-protein set per strain."""
from pathlib import Path
import json
import pandas as pd
import multiomics_explorer as m

OUT = Path(__file__).resolve().parents[1] / "data"
STRAINS = ["MED4", "MIT9313"]
NAMES = ["ntcA", "glnA", "glnB", "amt1", "urtA", "ureC", "pipX"]
rows, log = [], []
for s in STRAINS:
    for n in NAMES:
        r = m.resolve_gene(n, organism=s)
        log.append(dict(tool="resolve_gene", identifier=n, organism=s,
                        total_matching=r["total_matching"], returned=r["returned"]))
        hits = r["results"]
        if r["total_matching"] == 1:
            h = hits[0]
            rows.append(dict(strain=s, locus_tag=h["locus_tag"], gene_name=h["gene_name"],
                             product=h["product"], role="ntcA" if n == "ntcA" else "target",
                             query_name=n, resolution_note="resolve_gene: unique hit"))
        else:
            tags = ";".join(h["locus_tag"] for h in hits)
            rows.append(dict(strain=s, locus_tag="", gene_name="", product="",
                             role="ntcA" if n == "ntcA" else "target", query_name=n,
                             resolution_note=f"resolve_gene: {r['total_matching']} hits [{tags}]; not assigned"))
    # ribosomal set: genes_by_function('ribosomal protein'), keep gene_name prefix rpl/rps/rpm
    res, off = [], 0
    while True:
        g = m.genes_by_function("ribosomal protein", organism=s, limit=1000, offset=off)
        res += g["results"]; off += g["returned"]
        if g["returned"] == 0 or off >= g["total_matching"]:
            break
    assert len(res) == g["total_matching"], (len(res), g["total_matching"])
    log.append(dict(tool="genes_by_function", identifier="ribosomal protein", organism=s,
                    total_matching=g["total_matching"], returned=len(res)))
    seen = set()
    for h in res:
        gn = (h.get("gene_name") or "")
        if gn.lower().startswith(("rpl", "rps", "rpm")) and h["locus_tag"] not in seen:
            seen.add(h["locus_tag"])
            rows.append(dict(strain=s, locus_tag=h["locus_tag"], gene_name=gn, product=h.get("product", ""),
                             role="ribosomal", query_name="",
                             resolution_note="genes_by_function('ribosomal protein') + gene_name prefix rpl/rps/rpm"))
df = pd.DataFrame(rows)
df.to_csv(OUT / "panel.csv", index=False, encoding="utf-8")
pd.DataFrame(log).to_csv(OUT / "panel_kg_calls.csv", index=False, encoding="utf-8")
print(pd.DataFrame(log).to_string())
print(df[df.role != "ribosomal"].to_string())
print(df.groupby(["strain", "role"]).size())
