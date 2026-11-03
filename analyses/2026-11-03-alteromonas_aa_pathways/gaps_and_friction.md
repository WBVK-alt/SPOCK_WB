# Gaps and friction

- `[tool]` **AltDE cannot be addressed on its own.** Every `organism=` string for *Alteromonas mediterranea* AltDE ("AltDE", the full preferred name, "mediterranea DE") raises "matches multiple organisms" because the word match also hits AltDE1. AltDE was excluded from the grounding inventory (9 of 10 Alteromonas genome strains used). Seen in grounding/scripts/pull_inventory.py.
- `[kg-gap]` **No KEGG modules.** The KEGG ontology carries pathway (level 2) and KO (level 3) terms but no module definitions, so per-step pathway completeness cannot be read from the KG directly; it would need a step list built in the analysis.
- `[kg-gap]` **No Alteromonas metabolite measurements.** `list_organisms` reports `has_metabolomics: 0` for all 12 Alteromonas entries — amino-acid release by Alteromonas cannot be checked in the KG.
