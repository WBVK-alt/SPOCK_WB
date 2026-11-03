# Proposal notebook — brainstorming record

## Grounding queries (2026-11-03, KG 0.1.0-alpha.9, explorer 0.1.0-alpha.6, verdict ok)

1. `list_experiments(search_text="urea")` → 4 matches. Two are urea-vs-NH4+ microarrays
   (Tolonen 2006; MED4 and MIT9313). Two are Pandhal light proteomics, which match "urea"
   only through the text and are irrelevant.
2. `resolve_gene("ntcA")` → 102 matches across strains. MED4 = PMM0246, MIT9313 = PMT1831.
3. `differential_expression_by_gene(MED4 urea experiment, significant_only)` → 19 significant
   genes. ntcA is not among them. Mostly hli, desC4 and unknowns.
4. `differential_expression_by_gene(PMM0246 / PMT1831, treatment_type=nitrogen)`:
   - PMM0246 urea: log2FC −2.27, padj 1, rank 6 (not significant)
   - PMT1831 urea: log2FC +0.86, padj 1, rank 189 (not significant)
   - N deprivation (same paper): MED4 +4.45 (6 h), +3.02 (12 h), +1.68 (24 h) significant;
     MIT9313 +3.03/+3.56/+3.54/+2.91 (6/12/24/48 h) significant
   - Also seen: PMT1831 significantly up on nitrite (+1.43), PMM0246 flat on cyanate (+0.03)

## Alternatives considered
- Using the later N-starvation studies (2025 preprint, ISME 2017) as the positive control:
  rejected. They use a different platform/baseline design; the same-paper control is cleaner.
- Regulon inference from binding motifs: out of scope (not in KG tooling for this exercise).

## Open item noted
- The MED4 urea ntcA value (−2.27, rank 6) is large but not significant. That suggests
  high replicate variance, or a significance call that needs more than fold change.
  Look at this in the analysis milestone; do not interpret it now.
