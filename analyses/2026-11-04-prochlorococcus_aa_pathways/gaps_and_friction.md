# Gaps and friction

- `[tool]` **PAC1 cannot be addressed on its own.** `organism="Prochlorococcus PAC1"` raises "matches multiple organisms" (also hits EQPAC1). PAC1 was left out of the gene inventory (44 of 45 genomes). Seen in grounding/scripts/pull_inventory.py. Same failure as AltDE / AltDE1 in 2026-11-03-alteromonas_aa_pathways (second occurrence).
- `[pitfall]` **Run scripts from the repo root.** The explorer reads `.env` from the working directory; running a script from inside the analysis folder fails with a connection refused to localhost:7687.
- `[kg-gap]` **No KEGG modules.** Pathway completeness cannot be read from the KG directly (same as the Alteromonas analysis).
- `[pitfall]` **Genome-wide "% significant" is meaningless for partial tables.** For significant_only / filtered_subset experiments the rollup reads ~100%; the report blanks it for those experiments.
- `[kg-gap]` **MIT0801 has metabolomics but no expression rows** for these genes; MIT9303 has 1 row and MIT9301 12 rows (significant-only tables).
