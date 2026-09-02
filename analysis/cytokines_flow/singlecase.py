"""
Single-case vs small-control-sample inference for the VEXAS/MAS study.

Implements the estimators named in the manuscript Methods:

  * Crawford-Howell (1998) modified t-test, df = n_controls - 1
  * standardised difference  delta = (x - M) / SD   (control-SD units)
  * Crawford-Garthwaite (2007) Bayesian percentile: posterior mean and
    95% credible interval for the proportion of the control population
    falling below the patient
  * leave-one-control-out jackknife (worst-case two-sided p)
  * model-free no-overlap statistic (lowest patient value / highest control)
  * Benjamini-Hochberg FDR across analytes

Reference implementations for cross-checking: R package `singcar`
(functions TD and BTD).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
from scipy import stats

__all__ = [
    "CaseResult",
    "crawford_howell",
    "bayesian_percentile",
    "jackknife_worst_p",
    "no_overlap",
    "benjamini_hochberg",
]


@dataclass
class CaseResult:
    """Result of one single-case vs control-sample comparison."""

    case: float
    n_controls: int
    control_mean: float
    control_sd: float
    t: float
    df: int
    p_two_sided: float
    delta: float
    pct_below: float
    pct_below_lo: float
    pct_below_hi: float

    def as_dict(self) -> dict:
        return asdict(self)


def crawford_howell(case: float, controls) -> CaseResult:
    """Crawford-Howell modified t-test for a single case vs a control sample.

    The variance term carries the extra ``(n + 1) / n`` factor that accounts
    for the control mean and SD being estimated from a finite (here: very
    small) sample, so the reference distribution is t with ``n - 1`` df
    rather than the normal.

    Parameters
    ----------
    case : float
        The single patient observation.
    controls : array-like
        Control observations. Must contain at least 2 finite values.

    Returns
    -------
    CaseResult
    """
    c = np.asarray(controls, dtype=float)
    c = c[np.isfinite(c)]
    n = c.size
    if n < 2:
        raise ValueError(f"need >= 2 finite control values, got {n}")
    if not np.isfinite(case):
        raise ValueError("case value must be finite")

    mean = c.mean()
    sd = c.std(ddof=1)
    if sd == 0:
        raise ValueError("control SD is zero; t statistic undefined")

    df = n - 1
    t = (case - mean) / (sd * np.sqrt((n + 1) / n))
    p_two = float(2 * stats.t.sf(abs(t), df))
    delta = (case - mean) / sd
    pct, lo, hi = bayesian_percentile(case, c)

    return CaseResult(
        case=float(case),
        n_controls=int(n),
        control_mean=float(mean),
        control_sd=float(sd),
        t=float(t),
        df=int(df),
        p_two_sided=p_two,
        delta=float(delta),
        pct_below=pct,
        pct_below_lo=lo,
        pct_below_hi=hi,
    )


def bayesian_percentile(
    case: float,
    controls,
    n_draws: int = 200_000,
    seed: int = 20260123,
) -> tuple[float, float, float]:
    """Crawford-Garthwaite Bayesian percentile with a 95% credible interval.

    Monte-Carlo draws from the standard non-informative posterior for a
    normal control population: ``sigma^2 | data`` scaled-inverse-chi-square,
    ``mu | sigma^2, data`` normal. Each draw yields a percentile
    ``Phi((x - mu) / sigma)``; the posterior mean of that quantity equals
    ``1 - p_one_sided`` from the Crawford-Howell test, which is asserted in
    the test suite.

    Returns
    -------
    (posterior_mean, ci_lower, ci_upper) as percentages.
    """
    c = np.asarray(controls, dtype=float)
    c = c[np.isfinite(c)]
    n = c.size
    mean = c.mean()
    sd = c.std(ddof=1)
    df = n - 1

    rng = np.random.default_rng(seed)
    # sigma^2 | data  ~  (n-1) s^2 / chi^2_{n-1}
    chi2 = rng.chisquare(df, size=n_draws)
    sigma = np.sqrt(df * sd**2 / chi2)
    # mu | sigma^2, data  ~  N(xbar, sigma^2 / n)
    mu = rng.normal(mean, sigma / np.sqrt(n))

    z = (case - mu) / sigma
    pct = stats.norm.cdf(z) * 100.0
    lo, hi = np.percentile(pct, [2.5, 97.5])
    return float(pct.mean()), float(lo), float(hi)


def jackknife_worst_p(case: float, controls) -> tuple[float, int]:
    """Leave-one-control-out jackknife.

    Refits the Crawford-Howell test n times, each time dropping one control,
    and returns the worst-case (largest) two-sided p together with the index
    of the dropped control that produced it. A conclusion that survives this
    is not resting on a single control observation.
    """
    c = np.asarray(controls, dtype=float)
    c = c[np.isfinite(c)]
    worst_p, worst_i = -1.0, -1
    for i in range(c.size):
        kept = np.delete(c, i)
        p = crawford_howell(case, kept).p_two_sided
        if p > worst_p:
            worst_p, worst_i = p, i
    return float(worst_p), int(worst_i)


def no_overlap(case_values, controls) -> float:
    """Model-free separation: lowest patient value / highest control value.

    Assumption-free companion to the parametric test. Values > 1 mean the
    patient's *worst* observation still exceeds every control.
    """
    cv = np.asarray(case_values, dtype=float)
    cv = cv[np.isfinite(cv)]
    c = np.asarray(controls, dtype=float)
    c = c[np.isfinite(c)]
    return float(cv.min() / c.max())


def benjamini_hochberg(pvals) -> np.ndarray:
    """Benjamini-Hochberg step-up FDR adjustment.

    Self-contained so the repository has no statsmodels dependency; the test
    suite cross-checks it against ``statsmodels.stats.multitest.multipletests``
    when that package is available.
    """
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / np.arange(1, n + 1)
    # enforce monotonicity from the largest p downwards
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adj, 0, 1)
    return out
