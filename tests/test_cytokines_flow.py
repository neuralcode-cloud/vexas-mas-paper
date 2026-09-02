"""
Regression tests pinning the values published in Tables S3-S4.

Run:  pytest -q

These tests fail if a code change alters any number in the manuscript: a run
that reproduces the published result passes, and one that does not fails naming
the value that differs.
"""

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "analysis" / "cytokines_flow"))

from singlecase import (  # noqa: E402
    benjamini_hochberg,
    crawford_howell,
    jackknife_worst_p,
    no_overlap,
)
from run_analysis import (  # noqa: E402
    build_primary,
    build_trajectories,
    load,
)


# ---------------------------------------------------------------- estimator

def test_crawford_howell_matches_closed_form():
    """t = (x - M) / (SD * sqrt((n+1)/n)) on df = n-1."""
    controls = [1.0, 2.0, 3.0, 4.0]
    r = crawford_howell(10.0, controls)
    c = np.asarray(controls)
    expected_t = (10.0 - c.mean()) / (c.std(ddof=1) * np.sqrt(5 / 4))
    assert r.df == 3
    assert r.t == pytest.approx(expected_t)
    assert r.p_two_sided == pytest.approx(2 * stats.t.sf(abs(expected_t), 3))
    assert r.delta == pytest.approx((10.0 - c.mean()) / c.std(ddof=1))


def test_bayesian_posterior_mean_equals_one_minus_one_sided_p():
    """Crawford-Garthwaite posterior mean == 1 - one-sided p (to MC error)."""
    controls = [0.5, 1.2, 0.9, 1.5]
    r = crawford_howell(4.0, controls)
    one_sided = stats.t.sf(r.t, r.df)
    assert r.pct_below / 100 == pytest.approx(1 - one_sided, abs=2e-3)
    assert r.pct_below_lo < r.pct_below < r.pct_below_hi


def test_crawford_howell_rejects_degenerate_input():
    with pytest.raises(ValueError):
        crawford_howell(1.0, [2.0])            # n < 2
    with pytest.raises(ValueError):
        crawford_howell(1.0, [2.0, 2.0, 2.0])  # zero SD


def test_benjamini_hochberg_against_statsmodels():
    p = [0.0002, 0.0321, 0.0042, 0.0053, 0.0475, 0.0073, 0.0165, 0.0149, 0.6688]
    mine = benjamini_hochberg(p)
    sm = pytest.importorskip("statsmodels.stats.multitest")
    theirs = sm.multipletests(p, method="fdr_bh")[1]
    assert np.allclose(mine, theirs)


def test_benjamini_hochberg_is_monotone_and_bounded():
    rng = np.random.default_rng(0)
    p = np.sort(rng.uniform(0, 1, 40))
    adj = benjamini_hochberg(p)
    assert np.all(np.diff(adj) >= -1e-12)
    assert adj.min() >= 0 and adj.max() <= 1
    assert np.all(adj >= p - 1e-12)


def test_no_overlap_uses_worst_patient_value():
    assert no_overlap([100, 50, 200], [10, 5, 2]) == pytest.approx(5.0)


def test_jackknife_returns_worst_case():
    controls = [1.0, 1.1, 1.2, 5.0]
    worst_p, idx = jackknife_worst_p(20.0, controls)
    per_p = [
        crawford_howell(20.0, np.delete(controls, i)).p_two_sided
        for i in range(4)
    ]
    assert worst_p == pytest.approx(max(per_p))
    assert idx == int(np.argmax(per_p))


# ------------------------------------------------------- published values

# Table S4: analyte -> (delta, p_BH, pct_below, no_overlap) as printed.
PUBLISHED_PRIMARY = {
    "CXCL9 (pg/mL)":            (24.5, 0.002, 99.9, 19),
    "CXCL10 (pg/mL)":           (4.2,  0.041, 98.4, 4.1),
    "IFN-g (pg/mL)":            (8.8,  0.016, 99.8, 463),
    "IL-6 (pg/mL)":             (8.1,  0.016, 99.7, 23),
    "IL-10 (pg/mL)":            (3.6,  0.053, 97.6, 7.7),
    "TNF-a (pg/mL)":            (7.3,  0.016, 99.6, 8.2),
    "HLH-like population (%)":  (5.5,  0.025, 99.2, 15),
    "HLH-like (cells/uL)":      (5.7,  0.025, 99.3, 19),
    "CD8+ T cells (fraction)":  (-0.5, 0.669, 33.4, None),
}

# Table S3: analyte -> (delta_T1, delta_T2, delta_T3) as printed.
PUBLISHED_TRAJECTORY = {
    "CXCL9 (pg/mL)":           (24.5, 17.4, 13.7),
    "CXCL10 (pg/mL)":          (4.2,  3.1,  2.9),
    "IL-6 (pg/mL)":            (8.1,  8.2,  8.5),
    "IL-10 (pg/mL)":           (3.6,  4.8,  5.4),
    "TNF-a (pg/mL)":           (7.3,  7.5,  7.5),
    "HLH-like population (%)": (5.5,  5.3,  5.0),
    "HLH-like (cells/uL)":     (5.7,  5.8,  5.4),
    "CD8+ T cells (fraction)": (-0.5, 0.7,  0.0),
}


@pytest.fixture(scope="module")
def primary():
    return build_primary(load()).set_index("analyte")


@pytest.fixture(scope="module")
def trajectories():
    return build_trajectories(load()).set_index("analyte")


@pytest.mark.parametrize("analyte", list(PUBLISHED_PRIMARY))
def test_primary_table_matches_publication(primary, analyte):
    delta, p_bh, pct, overlap = PUBLISHED_PRIMARY[analyte]
    row = primary.loc[analyte]
    assert row["delta"] == pytest.approx(delta, abs=0.05)
    assert row["p_BH"] == pytest.approx(p_bh, abs=0.001)
    assert row["pct_controls_below"] == pytest.approx(pct, abs=0.1)
    if overlap is not None:
        rel = 0.05 if overlap < 100 else 0.01
        assert row["no_overlap_ratio"] == pytest.approx(overlap, rel=rel)


@pytest.mark.parametrize("analyte", list(PUBLISHED_TRAJECTORY))
def test_trajectory_table_matches_publication(trajectories, analyte):
    d1, d2, d3 = PUBLISHED_TRAJECTORY[analyte]
    row = trajectories.loc[analyte]
    assert row["delta_T1"] == pytest.approx(d1, abs=0.05)
    assert row["delta_T2"] == pytest.approx(d2, abs=0.05)
    assert row["delta_T3"] == pytest.approx(d3, abs=0.05)


def test_ifng_has_only_one_valid_timepoint(trajectories):
    """IFN-g: T1 is analysed; T2 and T3 are excluded."""
    row = trajectories.loc["IFN-g (pg/mL)"]
    assert row["n_valid_timepoints"] == 1
    assert np.isfinite(row["delta_T1"])
    assert np.isnan(row["delta_T2"]) and np.isnan(row["delta_T3"])


def test_il10_is_not_significant_after_fdr(primary):
    """IL-10: raw p 0.047, BH-adjusted p 0.053."""
    row = primary.loc["IL-10 (pg/mL)"]
    assert row["p_two_sided"] < 0.05
    assert row["p_BH"] > 0.05
    assert not row["significant_BH_005"]


def test_cd8_fraction_is_null_result(primary):
    row = primary.loc["CD8+ T cells (fraction)"]
    assert abs(row["delta"]) < 1.0
    assert row["p_BH"] > 0.5
    assert not row["significant_BH_005"]


def test_excluded_rows_never_enter_any_statistic():
    """Excluded rows do not enter the control mean or SD."""
    df = load()
    excluded = df[df.included == 0]
    assert len(excluded) == 8  # 2 IFN-g artefacts + 6 flare analytes
    assert set(excluded.exclusion_reason.dropna()) == {
        "anti-IFNg antibody interference",
        "flare sample - descriptive only",
    }
    ifng_ctrl = df[
        (df.analyte == "IFNg")
        & (df.group == "VEXAS_control")
        & (df.included == 1)
    ]["value"]
    assert len(ifng_ctrl) == 4
    assert 22960 not in set(df[df.included == 1]["value"])


def test_flare_il6_is_descriptive_only():
    """The 100-fold flare IL-6 claim, and its exclusion from inference."""
    df = load()
    flare = df[(df.group == "VEXAS_flare") & (df.analyte == "IL6")]["value"].iloc[0]
    mas_t1 = df[
        (df.group == "VEXAS_MAS") & (df.analyte == "IL6") & (df.timepoint == "T1")
    ]["value"].iloc[0]
    assert flare == 22960
    assert flare / mas_t1 == pytest.approx(101.6, abs=0.5)
    assert df[(df.group == "VEXAS_flare")]["included"].eq(0).all()


def test_control_subject_linkage_is_consistent():
    """Each CTRL label carries one individual's cytokine and flow profile.

    The two source worksheets list the four controls in different column orders
    (cytokines: C1, C2, C3, C4; flow cytometry: C1, C3, C2, C4), so the two
    sheets must be joined by subject identity rather than by column position.
    These are the resulting per-subject pairings.
    """
    df = load()
    expected = {
        "CTRL1": dict(CXCL10=321.0,  HLH_population=0.0299, HLH_like_cells=0.005863, CD8_fraction=0.1961),
        "CTRL2": dict(CXCL10=188.0,  HLH_population=0.0164, HLH_like_cells=0.005338, CD8_fraction=0.3255),
        "CTRL3": dict(CXCL10=54.7,   HLH_population=0.0056, HLH_like_cells=0.001454, CD8_fraction=0.2597),
        "CTRL4": dict(CXCL10=149.0,  HLH_population=0.0112, HLH_like_cells=0.002542, CD8_fraction=0.227),
    }
    ctrl = df[df.group == "VEXAS_control"]
    for subject, analytes in expected.items():
        for analyte, value in analytes.items():
            got = ctrl[(ctrl.subject == subject) & (ctrl.analyte == analyte)]["value"]
            assert len(got) == 1, f"{subject}/{analyte}: expected exactly one row"
            assert float(got.iloc[0]) == pytest.approx(value), (
                f"{subject}/{analyte}: {got.iloc[0]} != {value} "
                "(subject-linkage swap between source worksheets?)"
            )


def test_each_control_has_a_complete_profile():
    """No control may be missing an analyte, which would shift a group moment."""
    df = load()
    ctrl = df[(df.group == "VEXAS_control") & (df.included == 1)]
    counts = ctrl.groupby("subject")["analyte"].nunique()
    assert set(counts.index) == {"CTRL1", "CTRL2", "CTRL3", "CTRL4"}
    assert counts.eq(9).all(), f"incomplete control profiles: {counts.to_dict()}"


def test_mann_whitney_floor():
    """At 3 case observations vs 4 controls the smallest attainable two-sided
    rank-test p is 2/C(7,3) = 0.057."""
    from scipy.stats import mannwhitneyu

    df = load()
    a = df[(df.group == "VEXAS_MAS") & (df.analyte == "CXCL9") & (df.included == 1)]
    b = df[(df.group == "VEXAS_control") & (df.analyte == "CXCL9")]
    p = mannwhitneyu(a["value"], b["value"], alternative="two-sided").pvalue
    assert p == pytest.approx(2 / 35, abs=1e-4)
    assert p > 0.05
