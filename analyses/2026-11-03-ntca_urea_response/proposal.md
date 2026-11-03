# Proposal — Does NtcA respond when *Prochlorococcus* grows on urea instead of ammonium?

Status: approved 2026-11-03 (design review: all four recommendations accepted as written)

## Question

**Locked question:** In *Prochlorococcus* cultures acclimated to urea as the sole
nitrogen source, is the *ntcA* transcript (global nitrogen regulator) changed relative
to matched cultures on ammonium — and do NtcA's known target genes move with it?

**Purpose:** a quick standalone exercise of the research process (researcher's
request, ~10 min budget). Not part of a paper; no one's hypothesis is anchored to
the answer.

## KG entries (enumerated from the KG, release 0.1.0-alpha.9)

| Role | Experiment id | Organism | Platform | Design |
|---|---|---|---|---|
| Test | `10.1038/msb4100087_growth_medium_growth_on_urea_as_med4_microarray` | MED4 | microarray | urea vs NH4+, steady state, exponential; 1,696 genes, 19 significant |
| Test | `10.1038/msb4100087_growth_medium_growth_on_urea_as_mit9313_microarray` | MIT9313 | microarray | urea vs NH4+, steady state, exponential; 2,241 genes, 69 significant |
| Positive control | `10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray` | MED4 | microarray | −N vs +NH4+, 0–48 h time course |
| Positive control | `10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray` | MIT9313 | microarray | −N vs +NH4+, 0–48 h time course |

Publication: Tolonen et al. 2006, *Mol Syst Biol* (10.1038/msb4100087). All four
experiments come from one paper, one lab and one array platform. All have
table_scope `all_detected_genes`.

Genes: *ntcA* = **PMM0246** (MED4), **PMT1831** (MIT9313) [KG, `resolve_gene`].

Grounding observations [KG]:
- PMM0246 on urea: log2FC −2.27, padj 1 (not significant), rank 6 of 1,696 by |log2FC|.
- PMT1831 on urea: log2FC +0.86, padj 1 (not significant), rank 189.
- In the same paper's N-deprivation series, *ntcA* is significantly up in both strains
  (MED4 up to +4.45 at 6 h; MIT9313 up to +3.56 at 12 h).

## Framing

**Hypothesis.** Urea is a reduced N source that cells convert to ammonium
inside the cell, so steady-state growth on urea should not cause N stress. Expected
result: *ntcA* and its N-stress targets are **not** strongly induced on urea, unlike
under N deprivation. The alternative is that urea acts as a mild N limitation (slower
assimilation), which would show up as *ntcA* induction together with its targets.

Pre-stated readings for the other patterns:
- *ntcA* **down** on urea together with lower *glnA*/*amt1* would mean urea-grown cells are
  *less* N-stressed than cells on ammonium [interpretation]. *ntcA* down with targets flat
  would have no N-status reading; it would be reported as a gene-level observation.
- **Strains disagree in sign** (one up, one down): reported per strain as strain-specific.
  This would not count as evidence for or against N stress.

**Approach.**
1. Pull DE rows (log2FC, the paper's significance call, rank) for *ntcA* plus a fixed
   NtcA-target panel in the two urea experiments and the two N-deprivation experiments.
2. Compare the two contrasts gene by gene: urea vs ammonium against N deprivation vs ammonium.
3. Put each gene's urea |log2FC| in context as a percentile rank within its
   experiment's full gene distribution.

**NtcA-target panel (validation set).** Names are resolved to locus tags at the
methods milestone with `resolve_gene`; genes absent in a strain are reported as absent:
*glnA* (glutamine synthetase), *glnB* (PII), *amt1* (ammonium transporter),
*urtA* (urea ABC transporter substrate-binding), *ureC* (urease α), *pipX*,
*gifA/B* (if present).
- Expected in N deprivation (method check): *ntcA*, *glnA*, *amt1*, *glnB* significantly UP.
  If they are not, the extraction is wrong.
- Expected on urea [interpretation]: *urtA*/*ureC* may rise as substrate-specific
  responses. Stress targets (*glnA*, *amt1*) are expected small.

**Statistics plan.** No new statistical test. We use the paper's own significance
calls as stored in the KG. In this dataset padj looks coded (0.01 / 1), so we treat
it as a binary call, not a continuous p-value. This analysis is descriptive: n = 1
gene of interest × 2 strains, so a formal test would add nothing. The within-experiment
percentile of |log2FC| is reported as context only.

**Method check (not falsifiability).** In N deprivation, *ntcA*/*glnA*/*amt1*
should be significantly up. The *ntcA* part already passes in grounding. If the
targets fail, the extraction or the locus mapping is wrong, and we fix that before
reading the urea arm.

**Falsifiability check.**
- The urea prediction is falsified if *ntcA* **and** at least two of *glnA*/*amt1*/*glnB*
  are called significantly UP on urea in the same strain.
- Pre-registered expected-negative: ribosomal-protein genes (*rpl*/*rps*,
  all members in each strain) should not shift together on urea. This is judged on
  their |log2FC| distribution, not just hit count, because the base call rate is only
  1.1% (MED4) and 3.1% (MIT9313). **Fail condition:** a median within-experiment
  |log2FC| percentile above 75, **or** more than 80% of members sharing one sign. Either
  would mean the urea contrast carries a broad growth-rate shift, so *ntcA*/target
  changes cannot be read as N-specific.
- **Two kinds of non-call, reported separately.** "Not called, small |log2FC|" (below
  the experiment's 75th percentile) is a bounded negative: no detectable change.
  "Not called, large |log2FC|" (above the 95th percentile) is reported as **measured
  large change, not called significant by the paper's criterion**. That is a
  replicate-variance / calling-rule limit, not evidence of no change. Grounding
  already puts MED4 *ntcA* (−2.27, rank 6 of 1,696) in this second class.

**Design comparability.**

| Experiment | Role | Treatment arm | Reference arm and its state | Shared samples / model | Can confirm / cannot |
|---|---|---|---|---|---|
| Urea MED4 / MIT9313 | Test | log-phase, acclimated to urea | log-phase on NH4+ (Pro99), exponential; **assumed** matched by growth state (KG design notes: harvest day not stated) | Unverified whether the NH4+ reference arrays are shared with the cyanate/nitrite experiments or the N-deprivation +NH4+ arm. The methods milestone checks this if the KG holds sample IDs; otherwise it is a **gap** | Can show steady-state differences; cannot show transient induction during the switch |
| N-deprivation MED4 | Positive control | −N, 0–48 h | +NH4+ in parallel, matched | One time-course series | Shows the platform/extraction detects *ntcA* induction |
| N-deprivation MIT9313 | Positive control | −N, 0–48 h | +NH4+; **24 h point uses the 12 h control** (control_timing "mixed") | One series | Same as above; avoid over-reading the 24 h point |

Key caveat: the test contrasts are at acclimated steady state and the control contrasts
are acute stress. Comparing them shows whether the urea state *looks* N-stressed. It
does not show whether NtcA responded transiently during the switch.

**Methods milestone type.** Not discovery. The entity set is fixed (one gene plus a
named panel), so one spot-run against the MCP values above is enough verification.

## Out of scope

Other N sources (cyanate, nitrite), coculture experiments, proteomics, and
NtcA-binding-site (regulon) inference.
