# Proposal notebook — Alteromonas amino-acid pathways vs Prochlorococcus partner needs

## Dialogue so far (summary)
- Researcher's purpose: check whether Alteromonas can make an amino acid its partner needs.
- Scope: all 20 amino acids; biosynthesis and degradation.
- Partner: Prochlorococcus (coculture partners in the KG: MED4 with HOT1A3, NATL2A with MIT1002, MIT9312 with EZ55).
- Researcher instruction: collect and show the data first, no interpretation or inference.

## KG grounding (2026-11-03)
- `kg_release_info`: KG 0.1.0-alpha.9, explorer 0.1.0-alpha.6, verdict ok.
- `list_organisms(['Alteromonas'])`: 12 entries — 10 genome strains, 1 reference proteome (MarRef v6), 1 treatment label; 0 with metabolomics.
- `list_publications(organism='Alteromonas')`: 13 publications (RNA-seq 9, proteomics 4, bioinformatic 2, vesicle proteomics 1, exoproteomics 1).
- `ontology_term_details('kegg.subcategory:09105')`: 14 KEGG amino-acid metabolism pathways. No KEGG module terms seen.
- Raw inventory pulled by `grounding/scripts/pull_inventory.py` → `grounding/data/` (every call and its row count in `grounding/data/kg_calls.csv`); rendered by `grounding/scripts/build_report.py` → `grounding/reports/inventory.html`.
  - Genomes: 9 Alteromonas (AltDE excluded, see gaps_and_friction.md) + Prochlorococcus MED4, NATL2A, MIT9312.
  - GO-BP: 40 searches (20 amino acids × biosynthetic/catabolic), exact term-name match; 'aspartate catabolic process' and 'isoleucine catabolic process' have no GO term.
  - Not yet checked: whether `genes_by_ontology` on a GO term includes genes annotated only to child terms.
