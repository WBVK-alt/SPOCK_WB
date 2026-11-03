# Critical review: methods + analysis (data-integrity + interpretation), 2026-11-03

No Blockers. The data checks were clean: no truncation, call_class matches expression_status in all 840 rows, checks.csv matches panel_de.csv, and the figure's outlines are drawn from the call column.

| # | Severity | Finding | Disposition |
|---|---|---|---|
| 1 | Concern | The MIT9313 glnB disposition quietly sets aside the pre-registered rule ("if not up, extraction is wrong"). glnB is called *down* in both arms (N-dep 48h −1.754; urea −1.619), not just "not up". Row counts rule out truncation only. | **Fixed.** The notebook and paper now state this as a deviation, report glnB as down in both experiments, and say a probe or mapping artefact is not excluded. The falsification outcome is unchanged (0/3 targets called up either way). |
| 2 | Concern | "Noisy ntcA on MED4" rests on two values and can't be tested with binary padj. Calling −2.27 a "drop" contradicts the proposal's non-call wording. | **Fixed.** "Noisy" is now labelled a hypothesis. The urea value is reported as "large change, measured but not called significant". |
| 3 | Note | MIT9313 mid-percentile non-calls (ntcA, amt1, glnA) were left out. | **Fixed.** Added to the paper's Results with no direction claim. |
| 4 | Note | The MIT9313 ribosomal spread is elevated (median pct 62.2 vs 27.2), and the rpm prefix goes beyond the proposal's rpl/rps. | **Fixed (caveat added)** for the spread. rpm is accepted: those genes are ribosomal proteins (every product in the set contains "ribosomal"). |
| 5 | Note | 138 MIT9313 rows use the `PMT_` locus format. | **Deferred.** No effect here: every panel locus matched. |
