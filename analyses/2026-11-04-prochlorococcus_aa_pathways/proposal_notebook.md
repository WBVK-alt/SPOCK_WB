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

## Pathway map (2026-11-04, researcher asked for a drawing of how the pathways meet the TCA cycle, with strain counts)
- `grounding/scripts/pathway_steps.py` (layout + KO per step, textbook biochemistry [interpretation]) → `pull_ko_presence.py` (45 genomes; PAC1 via Cypher, cross-checked) → `build_pathway_map.py` → `grounding/reports/pathway_map.html`, `grounding/figures/pathway_map.svg`.
- Alternative enzymes added after MED4 name search: mqo (malate → OAA), PSTA (by ortholog group) and psp (serine), aroH (chorismate mutase), metA / metB (methionine), AGXT (alanine).

## Expression overview (2026-11-04, researcher asked for an overview of DE genes across conditions)
- `expression/scripts/pull_map_de.py` → `map_de_rows.csv` (DE rows for every map locus tag in the 10 strains with expression data; 949 up, 1,298 down, 5,018 not significant rows).
- `expression/scripts/build_de_overview.py` → `expression/reports/de_overview.html`, `expression/figures/de_overview.svg`, `data/de_overview_cells.csv` (101 genes × 67 experiments; any timepoint, any paralog).

## Expression pathway maps by condition (2026-11-04, researcher asked for three maps)
- `expression/scripts/build_de_maps.py` → `expression/reports/de_maps.html` (tabs), `expression/figures/de_map_{nitrogen,lightdark,mitomycin}.svg`.
- Experiment sets (chosen in chat, researcher to confirm): nitrogen = MED4 low-N starvation ×4 (RNA/protein, axenic/coculture), MED4 N-depleted (ismej.2017.88), N deprivation MED4 + MIT9313 (msb4100087); excluded urea/cyanate/nitrite as N source and azaserine. Light/dark = MED4 extended darkness, NATL2A extended darkness axenic + coculture, MED4 constant dark vs constant light (uninfected), MED4 extended light (kept as its own square, opposite treatment). Mitomycin = MIT0604 MMC vs untreated (MMC-vs-UV contrast excluded: no untreated control).

## Direction map (2026-11-04)
- `grounding/scripts/build_direction_map.py` → `grounding/reports/direction_map.html`, `grounding/figures/direction_map.svg`, `data/route_direction.csv`. Per-step reversibility is textbook biochemistry [interpretation], not KG data (KG reaction records not checked — the probe was interrupted).
