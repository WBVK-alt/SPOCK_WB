"""Build analysis/reports/summary.html (static, figure inline as SVG)."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
svg = (ROOT / "figures" / "panel_log2fc_heatmap.svg").read_text(encoding="utf-8")
svg = svg[svg.find("<svg"):]
pde = pd.read_csv(ROOT / "data" / "panel_de.csv", encoding="utf-8")
urea = pde[(pde.experiment == "urea") & (pde.role != "ribosomal")]
checks = pd.read_csv(ROOT / "data" / "checks.csv", encoding="utf-8")
css = """body{font-family:system-ui,sans-serif;max-width:1100px;margin:2em auto;color:#222}
table{border-collapse:collapse;font-size:13px;margin:1em 0}td,th{border-bottom:1px solid #ddd;padding:4px 8px;text-align:left}
th{background:#f5f5f3}.fig svg{width:100%;height:auto}"""
html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>Analysis summary: ntcA on urea</title>
<style>{css}</style></head><body>
<h1>Analysis milestone: ntcA and target panel on urea vs NH4+</h1>
<p>Source: KG, Tolonen et al. 2006 (10.1038/msb4100087), four microarray experiments. Data: analysis/data/.</p>
<h2>Figure: panel log2FC</h2><div class="fig">{svg}</div>
<h2>Panel rows, urea contrast</h2>{urea.to_html(index=False, float_format=lambda x: f'{x:.3f}')}
<h2>Checks (checks.csv)</h2>{checks.to_html(index=False)}
</body></html>"""
(ROOT / "reports" / "summary.html").write_text(html, encoding="utf-8")
print("written", len(urea), len(checks))
