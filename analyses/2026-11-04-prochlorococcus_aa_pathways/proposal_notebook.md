# Proposal notebook — Prochlorococcus amino-acid biosynthesis and degradation pathways

## Dialogue so far (summary)
- Researcher's request (2026-11-04): all available KG data on Prochlorococcus amino-acid biosynthesis and degradation pathways; total genes per pathway; conditions under which these genes are up- or down-regulated.
- Treated as Plan-phase data collection (counts only, no interpretation), mirroring the grounding step of 2026-11-03-alteromonas_aa_pathways.

## KG grounding (2026-11-04)
- `kg_release_info`: KG 0.1.0-alpha.9, explorer 0.1.0-alpha.6, verdict ok.
- `list_organisms(['Prochlorococcus'])`: 45 genome strains; `list_experiments(organism='Prochlorococcus')`: 140 experiments (42 bioinformatic).
- Gene sets: 14 KEGG pathways under kegg.subcategory:09105; GO-BP exact-name terms '<aa> biosynthetic/catabolic process' (20 amino acids).
- `grounding/scripts/pull_inventory.py` → `grounding/data/` (44 genomes; PAC1 skipped, see gaps_and_friction.md). KEGG genes per pathway, median across genomes: from 2 (tryptophan metabolism) to 24 (cysteine and methionine; glycine, serine and threonine).
- `expression/scripts/pull_de.py` → 8,241 DE rows, 67 experiments, 10 strains (MED4 33 experiments, NATL2A 9, MIT9313 8, SS120 6, MIT9312 4, MIT0604 3, AS9601, NATL1A, MIT9301, MIT9303 1 each). Status: 1,003 significant_up, 1,408 significant_down, 5,830 not_significant rows.
- `expression/scripts/pull_genome_wide.py` → genome-wide up/down per experiment × timepoint (context).
- `expression/scripts/summarize.py` → per pathway × experiment / timepoint counts; `kegg_pathway_conditions_full_tables.csv` (35 experiments with full source tables, 8 strains).
- Report: `expression/reports/summary.html` (built by `build_report.py` + `report_template.html`).
- Observation [KG]: in experiments with full tables, every KEGG amino-acid pathway has ≥1 gene up in 5–19 and ≥1 gene down in 8–19 of 26–34 experiments; nitrogen experiments contribute the most. Genome-wide change in the same experiments ranges from ~1% to ~88% of measured genes, so raw counts are not yet comparable to background.
