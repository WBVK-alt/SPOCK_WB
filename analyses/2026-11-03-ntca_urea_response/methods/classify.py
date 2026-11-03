"""Call classification and ribosomal expected-negative rule (proposal.md)."""
import numpy as np


def classify_call(log2fc, significant, abs_pct):
    """abs_pct: within-experiment percentile (0-100) of |log2fc|."""
    if significant:
        return "called_up" if log2fc > 0 else "called_down"
    if abs_pct < 75:
        return "not_called_small"
    if abs_pct >= 95:
        return "not_called_large"
    return "not_called_mid"


def ribosomal_fail(abs_pcts, signs):
    """True if median abs percentile > 75 OR > 80% of members share one sign."""
    abs_pcts = np.asarray(abs_pcts, dtype=float)
    signs = np.sign(np.asarray(signs, dtype=float))
    n = len(signs)
    frac_pos = (signs > 0).sum() / n
    frac_neg = (signs < 0).sum() / n
    return bool(np.median(abs_pcts) > 75 or frac_pos > 0.8 or frac_neg > 0.8)
