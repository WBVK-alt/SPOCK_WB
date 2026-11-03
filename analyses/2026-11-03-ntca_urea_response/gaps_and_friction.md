# Gaps and friction

- [pitfall] genes_by_function silently stops at the page limit: the first MIT9313 run read 2000 of 2383 rows. Always compare total_matching with the number of rows returned.
- [kg-gap] The KG holds only the paper's final significance call (padj takes two values, 0.01 and 1.0). It has no replicate values, so a large change the paper did not call can't be investigated further.
- [kg-gap] The KG holds no sample IDs that would show whether the NH4+ reference arrays are shared between the urea and N-deprivation experiments.
- [tool] The artifact viewer did not render a JS-drawn SVG chart; static inline SVG worked.
- [bkm] Using the same paper and platform for the N-deprivation series as the positive control kept the method check free of platform confounds.
