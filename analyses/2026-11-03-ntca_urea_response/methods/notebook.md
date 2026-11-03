# Methods milestone: notebook

**Co-define:** For time, the researcher agreed to a single sign-off at the end for all Run milestones ("approved, but again in 10 mins", 2026-11-03). Scope is the approved proposal, with nothing added.

## Run-manifest (from the coding subagent, verbatim excerpt)
- `uv run python analyses/2026-11-03-ntca_urea_response/methods/scripts/build_panel.py`. It ran twice: the first run read only 2000 of 2383 rows for MIT9313 genes_by_function. The fix pages through every row, and that rerun produced the output.
- `uv run python -m pytest -q analyses/2026-11-03-ntca_urea_response/methods/test_classify.py` gave **4 passed**. The tests cover boundaries at abs_pct 74.99/75/94.99/95, a median of exactly 75 (counts as no fail), and a sign share of exactly 80% (no fail) vs 90% and 100% (fail).
- resolve_gene was called 14 times ({ntcA, glnA, glnB, amt1, urtA, ureC, pipX} × {MED4, MIT9313}). Every call returned total_matching=1.
- The ribosomal set is genes_by_function("ribosomal protein"), paged fully (MED4 1557/1557, MIT9313 2383/2383), then filtered to gene_name prefix rpl/rps/rpm. That leaves 53 genes per strain.
- Output: `methods/data/panel.csv`, 120 rows (per strain: 1 ntcA, 6 targets, 53 ribosomal).

## Decisions
- The panel is fixed by name, so there was no discovery step. Each name resolved to exactly one locus tag per strain.
- **Deviation:** *gifA/B*, which the proposal listed "if present", were left out of the subagent brief. They are recorded as not checked; the analysis does not depend on them.

## Decide-gate checklist
- Outputs: methods/scripts/build_panel.py, methods/classify.py, methods/test_classify.py, methods/data/panel.csv (+ panel_kg_calls.csv)
- Verification: pytest passed 4/4. The main thread spot-checked PMM0246/PMT1831/PMT2229/PMT1481 in panel_de.csv against the grounding MCP values, and all matched.
- Critical review: combined with the analysis milestone pass (see analysis/critical_review.md, or the clean note in analysis/notebook.md).
