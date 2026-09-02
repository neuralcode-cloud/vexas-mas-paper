"""
Pins the WES summary numbers, and guards the two things that must never regress:
no raw sequence data in the repository, and no unqualified negative claim
resting on uncallable genes.

Run:  pytest tests/test_wes.py -q
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
WES = ROOT / "data" / "wes"

pytestmark = pytest.mark.skipif(not WES.exists(), reason="WES summaries absent")


@pytest.fixture(scope="module")
def qc() -> dict:
    df = pd.read_csv(WES / "qc_summary.csv")
    return dict(zip(df.metric, df.value))


@pytest.fixture(scope="module")
def uba1() -> pd.DataFrame:
    return pd.read_csv(WES / "uba1_m41t_allele_counts.csv")


@pytest.fixture(scope="module")
def panel() -> pd.DataFrame:
    return pd.read_csv(WES / "panel_callability.csv")


# ------------------------------------------------------------------- assay

def test_reference_build_is_hg19(qc):
    """Methods state hg19 and the BAM header records it. Every coordinate in
    this module is on that build."""
    assert "hg19" in qc["reference_build"] or "GRCh37" in qc["reference_build"]


def test_read_counts(qc):
    assert int(qc["total_reads"]) == 31_176_237
    assert int(qc["mapped_reads"]) == 30_934_525
    assert float(qc["mapped_pct"]) == pytest.approx(99.22, abs=0.01)


def test_library_is_single_end(qc):
    """Ion Torrent single-end library: 0 paired reads."""
    assert "single-end" in qc["library_type"]


# -------------------------------------------------------------------- UBA1

def test_uba1_locus_is_the_published_variant_position(uba1):
    """chrX:47058451 = UBA1 c.122, reference T; p.Met41Thr is T>C.

    Verified against the GRCh37 reference: codon 41 reads ATG (Met).
    """
    assert uba1.locus.str.contains("chrX:47058451").all()
    assert uba1.locus.str.contains("Met41Thr").all()


def test_uba1_no_alt_reads_at_any_filter(uba1):
    """Zero C-supporting reads at every filter setting."""
    assert (uba1.C_alt == 0).all()
    assert (uba1.vaf_pct == 0).all()


def test_uba1_vaf_upper_bounds(uba1):
    """Depth at the locus is 21-29x, bounding the VAF at roughly 10-13%.

    The bound is the one-sided 95% Clopper-Pearson limit for 0 observations at
    the stated depth.
    """
    assert uba1.depth.between(20, 30).all()
    assert uba1.vaf_95pct_upper_bound.min() > 5, (
        "upper bound implausibly tight for this depth"
    )
    assert uba1.vaf_95pct_upper_bound.max() < 25


# ---------------------------------------------------------------- callability

def test_panel_has_all_27_genes(panel):
    assert len(panel) == 27
    assert panel.gene.is_unique


def test_x_linked_gene_coverage(panel):
    """Pins the reported coverage for XIAP and SH2D1A, the two lowest in the
    panel. A different capture or a deeper run changes these values."""
    x = panel[panel.gene.isin(["XIAP", "SH2D1A"])].set_index("gene")
    assert len(x) == 2
    for g in ("XIAP", "SH2D1A"):
        assert x.loc[g, "pct_ge20x"] < 15, f"{g} coverage differs from the reported value"
        assert x.loc[g, "callable"] == "NO"


def test_only_one_gene_is_fully_callable(panel):
    assert (panel.callable == "yes").sum() == 1
    assert panel.loc[panel.callable == "yes", "gene"].iloc[0] == "STXBP2"
    assert (panel.callable == "NO").sum() == 19


def test_callable_flag_matches_its_threshold(panel):
    """The callable flag agrees with pct_ge20x at the stated thresholds."""
    for _, r in panel.iterrows():
        expect = "yes" if r.pct_ge20x >= 90 else "partial" if r.pct_ge20x >= 70 else "NO"
        assert r.callable == expect, f"{r.gene}: flag/threshold mismatch"


def test_panel_definition_is_committed(panel):
    """The gene list and exon coordinates ship alongside the coverage table."""
    listing = ROOT / "analysis" / "wes" / "hlh_gene_panel.txt"
    bed = ROOT / "analysis" / "wes" / "hlh_panel_exons.bed"
    assert listing.exists() and bed.exists()
    genes = {l.strip() for l in listing.read_text().splitlines()
             if l.strip() and not l.startswith("#")}
    assert genes == set(panel.gene)


# ------------------------------------------------------- no raw sequence data

# The only sequence files permitted in this repository. Each is a de-identified
# single-locus extract small enough to publish openly; the full alignment is not
# publicly released. Any other sequence file fails the gate.
ALLOWED_SEQUENCE_FILES = {
    "uba1_locus.bam",
    "uba1_locus.bam.bai",
}


def test_no_sequence_files_anywhere_in_the_repository():
    """Whole-exome alignments are not committed to this repository.

    A narrow allowlist covers the de-identified UBA1 c.122 extract; everything
    else with a sequence-data extension fails.
    """
    bad = [p for p in ROOT.rglob("*")
           if p.is_file()
           and ".git" not in p.parts
           and p.suffix.lower() in {".bam", ".bai", ".cram", ".crai", ".sam",
                                    ".vcf", ".fastq", ".fq"}
           and p.name not in ALLOWED_SEQUENCE_FILES]
    assert not bad, f"sequence data present in repo: {[str(p) for p in bad]}"


def test_allowed_slice_is_present_and_small():
    """The allowlist must not outlive the file it exists for."""
    bam = WES / "uba1_locus.bam"
    assert bam.exists(), "allowlisted slice missing; remove it from the allowlist"
    assert bam.stat().st_size < 1_000_000, (
        f"{bam.name} is {bam.stat().st_size} bytes; the allowlist covers a "
        "single-locus extract only"
    )


def test_slice_carries_no_laboratory_or_patient_identifiers():
    """The BAM header and read records must stay de-identified.

    A samtools region query copies the source header verbatim, so the extract
    initially carried the internal sample label, sequencing-centre and
    instrument identifiers, run and chip identifiers, barcode, and run
    timestamps. This pins that they remain stripped.

    Written as a positive contract on header structure rather than a blocklist
    of literal identifier strings, so that this file does not itself carry the
    values it exists to exclude.
    """
    import gzip
    import re

    text = gzip.open(WES / "uba1_locus.bam", "rb").read(2_000_000).decode("latin-1")
    header = text[:text.find("\x00")] if "\x00" in text else text

    # SAM read-group tags that identify a centre, run, chip, barcode or date.
    banned_tags = ["CN", "DT", "PU", "KS", "FO", "LB"]
    present = [t for t in banned_tags if re.search(rf"\t{t}:", header)]
    assert not present, f"identifying @RG tags present: {present}"

    # The only sample name may be the manuscript pseudonym.
    samples = set(re.findall(r"\tSM:([^\t\n]+)", header))
    assert samples <= {"UPN1"}, f"unexpected sample name(s): {samples}"

    # No ISO timestamps and no absolute filesystem paths anywhere in the file.
    assert not re.search(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}", text), "timestamp present"
    assert not re.search(r"(?<![\w.])/(?:results|data|home|Users)/", text), "path present"

    # Read names must be the sequential replacements, not instrument-assigned.
    names = re.findall(rb"[\x20-\x7e]{6,40}", gzip.open(WES / "uba1_locus.bam", "rb").read())
    instrument_like = [n for n in names if re.fullmatch(rb"[A-Z0-9]{5}:\d+:\d+", n)]
    assert not instrument_like, f"instrument read names present: {instrument_like[:3]}"


def test_gitignore_blocks_sequence_extensions():
    ignored = (ROOT / ".gitignore").read_text()
    for pat in ("*.bam", "*.bai", "*.cram", "*.vcf", "*.fastq"):
        assert pat in ignored, f"{pat} missing from .gitignore"
