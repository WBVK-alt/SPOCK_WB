# Proposal critical review (interpretation only), 2026-11-03

Fresh-context critic pass over proposal.md (draft d42d295). No Blockers.

| # | Severity | Finding (summary of critic's text) | Disposition |
|---|---|---|---|
| 1 | Concern | "Bounded negative / no detectable change" misdescribes MED4 ntcA (log2FC −2.27, rank 6 of 1,696, padj 1): this is a large change the paper did not call significant, not the absence of a change. | **Fixed.** Falsifiability section now separates "not called, small \|log2FC\|" (below 75th percentile, a bounded negative) from "not called, large \|log2FC\|" (above 95th percentile, reported as a large measured change not called significant). It names the significance call, not platform sensitivity, as the limit. |
| 2 | Concern | No pre-stated reading for a decrease, or for the two strains disagreeing in sign. | **Fixed.** Added pre-stated readings for "down" (with and without the targets) and for the strains disagreeing in sign. |
| 3 | Concern | The ribosomal expected-negative can hardly fail: the base call rate is 1.1% / 3.1%. | **Fixed.** It is now judged on the \|log2FC\| percentile distribution and sign consistency, with a stated fail condition (median percentile >75 or >80% sharing one sign). |
| 4 | Note | The N-deprivation positive control is a method check, not a falsifiability check. | **Fixed.** Relabelled "Method check". Falsifiability is now the urea prediction (ntcA plus at least two targets called UP). |
| 5 | Note | Whether the NH4+ reference arrays are shared is unverified, and the urea arms are only "assumed" matched. | **Fixed (wording) / deferred (check).** The table now says "assumed". The methods milestone checks sample sharing if the KG holds sample IDs; otherwise it is logged as a gap. |
