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

## Exploration — Q3 (researcher picked it from three offered questions; asked for an infographic and an answer)
- Question: are Alteromonas amino-acid genes switched on, and do they change, in coculture with Prochlorococcus vs axenic?
- `list_experiments(organism='Alteromonas')`: 63 experiments. Only one compares Alteromonas in coculture with Prochlorococcus against axenic Alteromonas: `10.1101/2025.11.24.690089_coculture_prochlorococcus_med4_hot1a3_rnaseq` (HOT1A3 + MED4, RNA-seq, days 11/18/31, all detected genes, matched sampling day).
  - MIT1002 + NATL2A: only 24/48 h vs 12 h coculture (UP-only), darkness vs diel, dark-tolerant vs parental — no axenic arm. [gap]
  - EZ55 + MIT9312: partner vs mean-of-three-partners, and CO2 within coculture — the paper's coculture-vs-axenic table is not among KG experiments. [gap]
- Script `exploration_q3_coculture_expression/scripts/pull_de.py` → `data/`; report `reports/q3_infographic.html`.
- Counts [KG]: 230 HOT1A3 amino-acid genes, all measured on all 3 days. Significant (paper's call): day 11 10 up / 7 down (7.4%; genome 7.1%); day 18 22 up / 18 down (17.4%; genome 13.5%); day 31 8 up / 3 down (4.8%; genome 11.6%).
- `pathway_enrichment` over the 14 KEGG amino-acid pathways × 3 days × up/down = 84 tests: 0 with adjusted p < 0.05 (lowest: arginine and proline metabolism, down, day 18, adj. p 0.071).
