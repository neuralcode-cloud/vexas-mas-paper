"""
Reproduce every statistic reported in Tables S3-S4 and Figure 2B-C.

Usage
-----
    python src/run_analysis.py

Writes to outputs/:
    table_s3_trajectories.csv   per-time-point standardised differences
    table_s4_primary.csv        T1-primary comparison vs VEXAS controls
    robustness.csv              jackknife, scale sensitivity, no-overlap
    figure_2c.png               cytokine and HLH panels, case vs controls
    figure_s4_flare.png         descriptive 3-condition comparison
    session_info.txt            package versions for the record
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from singlecase import (  # noqa: E402
    benjamini_hochberg,
    crawford_howell,
    jackknife_worst_p,
    no_overlap,
)

ROOT = Path(__file__).resolve().parents[2]      # repository root
DATA = ROOT / "data" / "cytokines_flow.csv"
OUT = ROOT / "outputs" / "cytokines_flow"

# Analyte display names and the order used in the manuscript tables.
ANALYTES = [
    ("CXCL9", "CXCL9 (pg/mL)"),
    ("CXCL10", "CXCL10 (pg/mL)"),
    ("IFNg", "IFN-g (pg/mL)"),
    ("IL6", "IL-6 (pg/mL)"),
    ("IL10", "IL-10 (pg/mL)"),
    ("TNFa", "TNF-a (pg/mL)"),
    ("HLH_population", "HLH-like population (%)"),
    ("HLH_like_cells", "HLH-like (cells/uL)"),
    ("CD8_fraction", "CD8+ T cells (fraction)"),
]

PRIMARY_TP = "T1"  # pre-specified primary comparison (pre-antibody, treatment-naive)


def load() -> pd.DataFrame:
    df = pd.read_csv(DATA)
    df["value"] = pd.to_numeric(df["value"], errors="raise")
    return df


def transform(values, scale: str) -> np.ndarray:
    """Apply the analysis scale declared per-analyte in the data file."""
    v = np.asarray(values, dtype=float)
    if scale == "log10":
        if np.any(v <= 0):
            raise ValueError("non-positive value on a log10 analyte")
        return np.log10(v)
    if scale == "linear":
        return v
    raise ValueError(f"unknown scale {scale!r}")


def slices(df: pd.DataFrame, analyte: str):
    """Return (case rows, control values, scale) for one analyte."""
    a = df[(df.analyte == analyte) & (df.included == 1)]
    case = a[a.group == "VEXAS_MAS"].sort_values("timepoint")
    ctrl = a[a.group == "VEXAS_control"]
    scale = a["scale"].iloc[0]
    return case, ctrl["value"].to_numpy(dtype=float), scale


def build_primary(df: pd.DataFrame) -> pd.DataFrame:
    """Table S4: T1-primary comparison, BH-adjusted across the nine analytes."""
    rows = []
    for key, label in ANALYTES:
        case, ctrl, scale = slices(df, key)
        prim = case[case.timepoint == PRIMARY_TP]
        if prim.empty:
            raise ValueError(f"{key}: no included {PRIMARY_TP} observation")
        x = float(prim["value"].iloc[0])

        r = crawford_howell(transform([x], scale)[0], transform(ctrl, scale))
        overlap = no_overlap(case["value"].to_numpy(dtype=float), ctrl)
        rows.append(
            dict(
                analyte=label,
                controls_min=ctrl.min(),
                controls_max=ctrl.max(),
                patient_T1=x,
                delta=r.delta,
                p_two_sided=r.p_two_sided,
                pct_controls_below=r.pct_below,
                crI_low=r.pct_below_lo,
                crI_high=r.pct_below_hi,
                no_overlap_ratio=overlap,
                scale=scale,
            )
        )
    out = pd.DataFrame(rows)
    out["p_BH"] = benjamini_hochberg(out["p_two_sided"].to_numpy())
    out["significant_BH_005"] = out["p_BH"] < 0.05
    return out


def build_trajectories(df: pd.DataFrame) -> pd.DataFrame:
    """Table S3: per-time-point delta as a within-patient robustness check.

    Reported per time point, NOT treated as independent replicates. Analytes
    with a single valid time point (IFN-g) yield a value at T1 only; the
    remaining cells are genuinely undefined and emitted as NaN.
    """
    rows = []
    for key, label in ANALYTES:
        case, ctrl, scale = slices(df, key)
        rec = {"analyte": label, "scale": scale}
        for tp in ("T1", "T2", "T3"):
            sel = case[case.timepoint == tp]
            if sel.empty:
                rec[f"delta_{tp}"] = np.nan
                rec[f"p_{tp}"] = np.nan
                continue
            x = transform([float(sel["value"].iloc[0])], scale)[0]
            r = crawford_howell(x, transform(ctrl, scale))
            rec[f"delta_{tp}"] = r.delta
            rec[f"p_{tp}"] = r.p_two_sided
        rec["n_valid_timepoints"] = int(case.shape[0])
        rows.append(rec)
    return pd.DataFrame(rows)


def build_robustness(df: pd.DataFrame) -> pd.DataFrame:
    """Three robustness checks named in the Methods."""
    rows = []
    for key, label in ANALYTES:
        case, ctrl, scale = slices(df, key)
        prim = case[case.timepoint == PRIMARY_TP]
        x = float(prim["value"].iloc[0])

        # 1. leave-one-control-out jackknife (worst-case p)
        wp, wi = jackknife_worst_p(transform([x], scale)[0], transform(ctrl, scale))

        # 2. raw- vs log-scale sensitivity
        p_log = crawford_howell(
            transform([x], "log10")[0], transform(ctrl, "log10")
        ).p_two_sided if (x > 0 and np.all(ctrl > 0)) else np.nan
        p_lin = crawford_howell(x, ctrl).p_two_sided

        # 3. model-free no-overlap
        rows.append(
            dict(
                analyte=label,
                primary_scale=scale,
                p_primary=crawford_howell(
                    transform([x], scale)[0], transform(ctrl, scale)
                ).p_two_sided,
                jackknife_worst_p=wp,
                jackknife_dropped_control=f"CTRL{wi + 1}",
                p_log10_scale=p_log,
                p_linear_scale=p_lin,
                no_overlap_ratio=no_overlap(
                    case["value"].to_numpy(dtype=float), ctrl
                ),
            )
        )
    return pd.DataFrame(rows)


def make_figures(df: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.titlesize": 8,
            "axes.labelsize": 8,
            "legend.fontsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.dpi": 300,
        }
    )
    CASE_C, CTRL_C = "#c0392b", "#7f8c8d"

    # ---- Figure 2C: case (T1->T3) vs controls, per analyte ----
    keys = [k for k, _ in ANALYTES]
    fig, axes = plt.subplots(3, 3, figsize=(7.2, 6.6))
    for ax, key in zip(axes.ravel(), keys):
        case, ctrl, scale = slices(df, key)
        cv = case["value"].to_numpy(dtype=float)
        ax.scatter(
            np.full(ctrl.size, 0) + np.linspace(-0.08, 0.08, ctrl.size),
            ctrl,
            s=18, facecolors="none", edgecolors=CTRL_C, linewidths=1.0, zorder=2,
        )
        ax.plot(
            np.full(cv.size, 1) + np.linspace(-0.08, 0.08, cv.size),
            cv, "-o", color=CASE_C, ms=4, lw=1.0, zorder=3,
        )
        if scale == "log10":
            ax.set_yscale("log")
        ax.set_xlim(-0.45, 1.45)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["VEXAS\nctrl (n=4)", "VEXAS+MAS\nT1\u2192T3"])
        ax.set_title(dict(ANALYTES)[key].split(" (")[0], loc="left")
        ax.set_ylabel(dict(ANALYTES)[key].split("(")[-1].rstrip(")"))
    handles = [
        Line2D([], [], marker="o", ls="none", mfc="none", mec=CTRL_C,
               label="VEXAS controls (n=4)"),
        Line2D([], [], marker="o", ls="-", color=CASE_C,
               label="VEXAS/MAS case (T1\u2192T3)"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.005))
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(OUT / "figure_2c.png", bbox_inches="tight")
    plt.close(fig)

    # ---- Figure S4: descriptive MAS vs flare vs no-flare ----
    cyto = ["CXCL9", "CXCL10", "IFNg", "IL6", "IL10", "TNFa"]
    fig2, ax = plt.subplots(figsize=(7.2, 3.6))
    for j, key in enumerate(cyto):
        a = df[df.analyte == key]
        mas = a[(a.group == "VEXAS_MAS") & (a.included == 1)]["value"].to_numpy(float)
        nofl = a[(a.group == "VEXAS_control") & (a.subject == "CTRL4")]["value"]
        flare = a[a.group == "VEXAS_flare"]["value"]
        ax.scatter(np.full(mas.size, j) + np.linspace(-0.1, 0.1, mas.size),
                   mas, s=32, color=CASE_C, zorder=3)
        if not nofl.empty:
            ax.scatter(j - 0.28, float(nofl.iloc[0]), s=46, marker="s",
                       color="#2980b9", zorder=3)
        if not flare.empty:
            ax.scatter(j + 0.28, float(flare.iloc[0]), s=46, marker="^",
                       color="#e67e22", zorder=3)
    ax.set_yscale("log")
    ax.set_ylabel("Serum concentration (pg/mL)")
    ax.set_xticks(range(len(cyto)))
    ax.set_xticklabels(["CXCL9", "CXCL10", "IFN-\u03b3", "IL-6", "IL-10", "TNF-\u03b1"],
                       rotation=20, ha="right")
    ax.set_title("Descriptive comparison: VEXAS/MAS vs VEXAS flare vs VEXAS no-flare",
                 loc="left")
    ax.margins(x=0.06)
    h2 = [
        Line2D([], [], marker="o", ls="none", color=CASE_C, label="VEXAS/MAS case (T1\u2192T3)"),
        Line2D([], [], marker="^", ls="none", color="#e67e22", label="VEXAS in flare"),
        Line2D([], [], marker="s", ls="none", color="#2980b9", label="VEXAS without flare"),
    ]
    ax.legend(handles=h2, frameon=False, loc="upper right")
    fig2.tight_layout()
    fig2.savefig(OUT / "figure_s4_flare.png", bbox_inches="tight")
    plt.close(fig2)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = load()

    primary = build_primary(df)
    traj = build_trajectories(df)
    robust = build_robustness(df)

    primary.to_csv(OUT / "table_s4_primary.csv", index=False)
    traj.to_csv(OUT / "table_s3_trajectories.csv", index=False)
    robust.to_csv(OUT / "robustness.csv", index=False)
    make_figures(df)

    import platform
    import scipy

    with open(OUT / "session_info.txt", "w") as fh:
        fh.write(
            f"python  {platform.python_version()}\n"
            f"numpy   {np.__version__}\n"
            f"pandas  {pd.__version__}\n"
            f"scipy   {scipy.__version__}\n"
        )

    show = primary[["analyte", "delta", "p_two_sided", "p_BH",
                    "pct_controls_below", "no_overlap_ratio"]]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
