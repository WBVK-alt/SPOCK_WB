"""Heatmap of panel log2FC (urea + N-dep timepoints) per strain; called cells outlined."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[1]
pde = pd.read_csv(ROOT / "data" / "panel_de.csv", encoding="utf-8")
pde = pde[pde.role != "ribosomal"].copy()
pde["col"] = np.where(pde.experiment == "urea", "urea", "N-dep " + pde.timepoint.astype(str))
COLS = ["urea"] + [f"N-dep {t}" for t in ["0h", "3h", "6h", "12h", "24h", "48h"]]
GENES = ["ntcA", "glnA", "glnB", "amt1", "pipX", "urtA", "ureC"]
cmap = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#f0efec", "#eb6834"])
vmax = np.ceil(pde.log2fc.abs().max())
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#888", "axes.linewidth": 0.6})
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
for ax, strain in zip(axes, ["MED4", "MIT9313"]):
    s = pde[pde.strain == strain]
    lab = {g: f"{g} ({s[s.gene_name == g].locus_tag.iloc[0]})" for g in GENES}
    M = s.pivot_table(index="gene_name", columns="col", values="log2fc").reindex(index=GENES, columns=COLS)
    S = s.assign(sig=s.call_class.isin(["called_up", "called_down"])).pivot_table(
        index="gene_name", columns="col", values="sig", aggfunc="max").reindex(index=GENES, columns=COLS)
    im = ax.imshow(M.values, cmap=cmap, norm=norm, aspect="auto")
    for i in range(len(GENES)):
        for j in range(len(COLS)):
            if S.values[i, j] == 1:
                ax.add_patch(plt.Rectangle((j - 0.44, i - 0.44), 0.88, 0.88, fill=False, ec="black", lw=1.6))
    ax.set_xticks(range(len(COLS)), COLS, rotation=40, ha="right")
    ax.set_yticks(range(len(GENES)), [lab[g] for g in GENES])
    ax.axvline(0.5, color="white", lw=3)
    ax.set_title(f"Prochlorococcus {strain}", loc="left", fontsize=10)
    ax.tick_params(length=0)
cb = fig.colorbar(im, ax=axes, shrink=0.8, label="log2FC (vs NH4+)")
cb.outline.set_visible(False)
fig.text(0.01, -0.04, "Black outline = called significant in the KG (paper's call). Orange = up, blue = down.",
         fontsize=8, color="#444")
for ext in ["png", "svg"]:
    fig.savefig(ROOT / "figures" / f"panel_log2fc_heatmap.{ext}", dpi=200, bbox_inches="tight")
print("saved")
