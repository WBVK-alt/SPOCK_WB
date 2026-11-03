# Analysis milestone: notebook

**Co-define:** a single combined sign-off for all Run milestones (see methods/notebook.md).

## Run-manifest (from the coding subagent, verbatim excerpt)
- Scripts run in this order: `analysis/scripts/pull_de.py`, `score.py`, `figure.py`, `build_reports.py` (each as `uv run python <path>` from the repo root), then `scripts/build_report_index.py`.
- In each DE pull the saved rows equal total_matching: MED4 urea 1696, MIT9313 urea 2241, MED4 N-dep 10176 (6 timepoints), MIT9313 N-dep 13446 (6 timepoints).
- `analysis/data/panel_de.csv` has 840 rows. `analysis/data/checks.csv` has 12 rows. Figure: `analysis/figures/panel_log2fc_heatmap.{png,svg}`. Report: `analysis/reports/summary.html`.
- Anomalies the subagent listed (facts): MIT9313 glnB (PMT1481) is called up at 0/6 N-dep timepoints (max +0.30) and is called down on urea (−1.62). MED4 ntcA is −2.17 at N-dep 0h (not called), and no other panel gene exceeds |0.80| at 0h. MED4 urtA is called up at N-dep 3–48h but is −0.16 on urea; MIT9313 urtA is +3.20 on urea (called). padj takes only two values, {0.01, 1.0}.

## Results (from checks.csv / panel_de.csv) [KG]
| Strain | ntcA on urea | Targets on urea | Falsification rule | Ribosomal expected-negative |
|---|---|---|---|---|
| MED4 | −2.27, not called, 99.7th pct → *large, not called* | all 6 small (pct 8–67), none called | not met (0/3 of glnA/amt1/glnB up) | pass: median pct 27.2, 62% negative, 0 called |
| MIT9313 | +0.86, not called, 91.6th pct → *mid* | urtA **called up** +3.20; glnB **called down** −1.62; amt1 +1.02, glnA −0.76 (mid, not called) | not met (0/3 up) | pass: median pct 62.2, 64% negative, 1 called |

Method check (N-deprivation, called up at ≥1 timepoint): MED4 ntcA, glnA, glnB and amt1 all pass. MIT9313 ntcA, glnA and amt1 pass; **glnB fails** (0/6).

## Surprises
- **Deviation from a pre-registered rule.** MIT9313 glnB fails the method check. The proposal (method-check paragraph) said a failure means "extraction or locus mapping is wrong". Against that: in the same strain and experiment, ntcA, glnA and amt1 pass; every pull's row count equals total_matching; and PMT1481 resolves uniquely to glnB. glnB is not merely flat: it is called **down** at N-dep 48h (−1.75) and on urea (−1.62), so it moves the same way in both experiments. My reading is a gene-specific result, not a pipeline fault [interpretation]. Matching row counts rule out truncation only; a probe or mapping problem is not excluded (unverified). Consequence: glnB is not informative in MIT9313, so that strain's falsification rule effectively becomes ntcA + glnA + amt1. Neither the original rule nor the reduced one is met.
- MED4 ntcA shows a similar uncalled drop at N-dep 0h (−2.17), a point where cells had only just been resuspended and other panel genes barely moved. One *hypothesis* is that the MED4 ntcA measurement is noisy; another is a shared or offset reference sample. This rests on two values, and binary padj cannot test it [interpretation][gap]. Per the proposal, the urea value is reported only as a large change that was measured but not called significant, not as a drop.
- MIT9313 urtA, the urea transporter, is strongly up on urea. That fits a urea-specific uptake response that does not involve broad N stress [interpretation]. MED4 urtA does not change on urea.

## Decisions
- glnB/MIT9313 is treated as uninformative rather than as an extraction failure, for the reasons above. The critic is asked to challenge this.

## Decide-gate checklist
- Outputs: listed above. Report: analysis/reports/summary.html; index.html regenerated.
- Results presented: the tables above, plus the figure.
- Verification: the main thread re-read panel_de.csv for PMM0246/PMT1831/PMT1481/PMT2229 and the ribosomal summaries, and they match checks.csv. In MED4 urea de_all, 52.2% of log2fc values are negative, so the sign is intact.
- Critical review (data-integrity + interpretation): 2 Concerns, 3 Notes, no Blockers. See critical_review.md; all fixed in wording. The headline is unchanged.
