"""
Pins the proteomics numbers reported in the manuscript, and documents the one
figure that could not be reproduced from the supplied files.

Run:  pytest tests/test_proteomics.py -q
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
NPX = ROOT / "data" / "proteomics" / "olink_long.csv"
META = ROOT / "data" / "proteomics" / "sample_metadata.csv"

pytestmark = pytest.mark.skipif(not NPX.exists(), reason="proteomics data absent")


@pytest.fixture(scope="module")
def npx() -> pd.DataFrame:
    df = pd.read_csv(NPX)
    return df.loc[:, [c for c in df.columns if not c.startswith("Unnamed")]]


@pytest.fixture(scope="module")
def meta() -> pd.DataFrame:
    return pd.read_csv(META)


# ------------------------------------------------------------------ cohort

def test_sample_and_assay_counts(npx):
    """24 plasma samples, 1033 assays."""
    assert npx.SampleID.nunique() == 24
    assert npx.Assay.nunique() == 1033
    assert npx.OlinkID.nunique() == 1033


def test_assay_count_in_the_deposited_table(npx):
    """The deposited NPX table carries 1033 distinct assays."""
    assert npx.Assay.nunique() == 1033


def test_cohort_composition(meta):
    """11 patients, 24 samples: onset N=8, post-DMARDs N=5, post-DMT N=5,
    HR-MDS N=3, HR-MDS post-azacitidine N=3."""
    assert meta.subject.nunique() == 11
    assert len(meta) == 24
    counts = meta.group.value_counts().to_dict()
    assert counts["VEXAS onset"] == 8
    assert counts["VEXAS post-DMARDs"] == 5
    assert counts["VEXAS post-DMT"] == 5
    assert counts["HR- MDS"] == 3
    assert counts["HR- MDS post AZA"] == 3


def test_every_npx_sample_has_metadata(npx, meta):
    assert set(npx.SampleID.unique()) == set(meta.sample_id)


def test_index_case_is_flagged_and_has_three_timepoints(meta):
    """The index case is UPN10 in the proteomics numbering and contributes
    three longitudinal samples."""
    idx = meta[meta.is_index_case == 1]
    assert set(idx.subject) == {"UPN10"}
    assert len(idx) == 3, f"index case has {len(idx)} samples, expected 3"


def test_all_qc_passing(npx):
    """The deposited table contains only QC-passing assays and samples."""
    assert set(npx.AssayQC.unique()) == {"PASS"}
    assert set(npx.SampleQC.unique()) == {"PASS"}
    assert set(npx.SampleType.unique()) == {"SAMPLE"}


# --------------------------------------------------------------- LOD filter

def _below_lod_counts(npx: pd.DataFrame) -> pd.Series:
    return npx[npx.NPX < npx.LOD].groupby("Assay").size()


def test_lod_100_percent_rule_gives_15_and_1018(npx):
    """Assays below LOD in 100% of samples: 15 discarded, 1018 retained.
    This is the count reported in the manuscript."""
    freq = _below_lod_counts(npx)
    n = npx.SampleID.nunique()
    n_all = int((freq == n).sum())
    assert n_all == 15
    assert npx.Assay.nunique() - n_all == 1018


def test_lod_95_percent_rule_gives_24_and_1009(npx):
    """A ">95% of samples" rule discards 24 assays, retaining 1009.

    With 24 samples one sample is 4.17%, so 23/24 = 95.83% also exceeds 95%;
    nine assays fall between this threshold and the 100% one. Both counts are
    pinned so either can be recomputed. See analysis/proteomics/README.md.
    """
    freq = _below_lod_counts(npx)
    n = npx.SampleID.nunique()
    pct = (freq / n * 100).round(2)
    discarded = set(pct[pct > 95].index)
    assert len(discarded) == 24
    assert npx.Assay.nunique() - len(discarded) == 1009

    borderline = sorted(discarded - set(freq[freq == n].index))
    assert len(borderline) == 9
    assert borderline == [
        "CENPS", "CSF2", "DDX60", "HDAC9", "HRAS",
        "IL20RA", "LTA4H", "NUB1", "SIKE1",
    ]


def test_lod_values_present_for_every_assay(npx):
    """Every assay in the deposited table has an LOD value."""
    assert npx.LOD.notna().all()


# ------------------------------------------------------------ de-identification

def test_no_lab_accession_numbers(meta):
    """No laboratory specimen accession numbers in the sample metadata.

    Cells are coerced with str() individually rather than via
    DataFrame.astype(str), whose handling of missing values varies by pandas
    version.
    """
    assert "ExpSampleID" not in meta.columns
    joined = " ".join(str(v) for v in meta.to_numpy().ravel())
    # Specimen accessions have the shape <2-3 letters> <digits>/<year>.
    accession = re.compile(r"\b[A-Z]{2,3}\s?\d{3,5}\s?/\s?\d{2,4}\b")
    hits = accession.findall(joined)
    assert not hits, f"laboratory accession-like values present: {hits[:5]}"


def test_subjects_are_pseudonymous(meta):
    assert meta.subject.str.fullmatch(r"UPN\d+").all()


def test_timepoints_are_relative(meta):
    assert meta.timepoint.str.fullmatch(r"T\d").all()
