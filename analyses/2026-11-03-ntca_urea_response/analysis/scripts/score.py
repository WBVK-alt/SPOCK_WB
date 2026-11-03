"""Join panel to DE, classify calls, run method check / falsification / ribosomal expected-negative."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "methods"))
from classify import classify_call, ribosomal_fail

D = ROOT / "analysis" / "data"
panel = pd.read_csv(ROOT / "methods" / "data" / "panel.csv", encoding="utf-8")
panel = panel[panel.locus_tag.notna() & (panel.locus_tag != "")]
parts = []
for strain in ["MED4", "MIT9313"]:
    for arm in ["urea", "ndep"]:
        de = pd.read_csv(D / f"de_all_{strain}_{arm}.csv", encoding="utf-8")
        de["experiment"] = arm
        p = panel[panel.strain == strain][["strain", "locus_tag", "gene_name", "role"]]
        parts.append(p.merge(de[["locus_tag", "experiment", "timepoint", "timepoint_hours", "log2fc",
                                 "expression_status", "abs_pct"]], on="locus_tag", how="left"))
pde = pd.concat(parts, ignore_index=True)
pde["significant"] = pde.expression_status.fillna("").str.startswith("significant")
pde["call_class"] = [classify_call(l, s, a) if pd.notna(l) else "no_data"
                     for l, s, a in zip(pde.log2fc, pde.significant, pde.abs_pct)]
pde = pde.sort_values(["strain", "role", "gene_name", "experiment", "timepoint_hours"])
cols = ["strain", "experiment", "timepoint", "locus_tag", "gene_name", "role", "log2fc",
        "expression_status", "abs_pct", "call_class"]
pde[cols].to_csv(D / "panel_de.csv", index=False, encoding="utf-8")

checks = []
for strain in ["MED4", "MIT9313"]:
    s = pde[pde.strain == strain]
    nd = s[s.experiment == "ndep"]
    for g in ["ntcA", "glnA", "amt1", "glnB"]:
        r = nd[nd.gene_name == g]
        up = r[r.call_class == "called_up"]
        checks.append(dict(check="method_check_ndep_called_up", strain=strain, gene=g,
                           locus_tag=";".join(r.locus_tag.unique()),
                           numbers=f"called_up at {len(up)}/{len(r)} timepoints: {';'.join(up.timepoint.astype(str))}; max log2fc {r.log2fc.max():.3f}",
                           result="pass" if len(up) >= 1 else "fail"))
    u = s[s.experiment == "urea"]
    ntca_up = (u[u.gene_name == "ntcA"].call_class == "called_up").any()
    tg = u[u.gene_name.isin(["glnA", "amt1", "glnB"])]
    n_tg_up = int((tg.call_class == "called_up").sum())
    falsified = bool(ntca_up and n_tg_up >= 2)
    checks.append(dict(check="falsification_urea_ntcA_and_2of3_up", strain=strain, gene="ntcA+glnA/amt1/glnB",
                       locus_tag=";".join(u[u.gene_name.isin(["ntcA", "glnA", "amt1", "glnB"])].locus_tag),
                       numbers=f"ntcA called_up={ntca_up}; targets called_up={n_tg_up}/3",
                       result="falsified" if falsified else "not_falsified"))
    rb = u[u.role == "ribosomal"].dropna(subset=["log2fc"])
    n = len(rb)
    fp, fn = (rb.log2fc > 0).mean(), (rb.log2fc < 0).mean()
    med = rb.abs_pct.median()
    n_called = int(rb.significant.sum())
    fail = ribosomal_fail(rb.abs_pct, rb.log2fc)
    checks.append(dict(check="ribosomal_expected_negative_urea", strain=strain, gene="rpl/rps/rpm",
                       locus_tag=f"{n} genes",
                       numbers=f"n={n}; median_abs_pct={med:.1f}; frac_pos={fp:.3f}; frac_neg={fn:.3f}; n_called={n_called}",
                       result="fail" if fail else "pass"))
C = pd.DataFrame(checks)
C.to_csv(D / "checks.csv", index=False, encoding="utf-8")
print(C.to_string())
print(pde.groupby(["strain", "role", "experiment"]).size())
