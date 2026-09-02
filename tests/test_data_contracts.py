"""
Repository-wide data contracts.

These tests apply to every module's data file, including files added later.
They pass when a module's data is absent and fail when a file present in the
repository breaks one of the conventions below.

Run:  pytest tests/ -q
"""

from __future__ import annotations

import csv
import datetime as dt
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

# Subject labels permitted anywhere in the repository's data.
SUBJECT_OK = re.compile(r"^(UPN\d+|CTRL\d+|HRMDS\d+|POOL|NA)$")

# Anything that looks like a real-world date: these are identifying under HIPAA
# safe-harbour; timepoints are relative labels instead.
DATE_LIKE = re.compile(
    r"\b(?:\d{1,2}[./-]\d{1,2}[./-]\d{2,4}"      # 23/01/2026, 1-2-26
    r"|\d{4}[./-]\d{1,2}[./-]\d{1,2})\b"          # 2026-01-23
)

# A month name adjacent to a day or year number. The adjacent number is part
# of the pattern: bare month abbreviations occur in HGNC gene symbols
# (DECR1, SEPTIN8, MARCHF5, SEPT9, JUN).
_MON = (
    r"(?:jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)"
    r"(?:uary|ruary|ch|il|e|y|ust|tember|ober|ember)?"
)
MONTHS = re.compile(
    rf"\b(?:\d{{1,2}}[\s,.-]+{_MON}\b|{_MON}\.?[\s,-]+\d{{1,4}}\b)",
    re.I,
)


def data_files() -> list[Path]:
    return sorted(DATA.rglob("*.csv")) if DATA.exists() else []


def read_rows(path: Path) -> tuple[list[str], list[dict]]:
    with path.open(newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh)
        return (r.fieldnames or []), list(r)


def test_at_least_one_data_file_exists():
    assert data_files(), "no CSV data files found under data/"


@pytest.mark.parametrize("path", data_files(), ids=lambda p: p.name)
def test_no_date_like_values(path: Path):
    """Timepoints must be relative labels (T1, T2), never calendar dates."""
    text = path.read_text(encoding="utf-8")
    offenders = set(DATE_LIKE.findall(text)) | set(MONTHS.findall(text))
    assert not offenders, (
        f"{path.name}: date-like values present ({sorted(offenders)[:5]}). "
        "Use relative timepoint labels; dates are identifying."
    )


@pytest.mark.parametrize("path", data_files(), ids=lambda p: p.name)
def test_subject_labels_are_pseudonymous(path: Path):
    """Subject columns may only contain UPN*/CTRL*-style labels."""
    header, rows = read_rows(path)
    col = next((c for c in header if c.lower() in ("subject", "subject_id", "patient")), None)
    if col is None:
        pytest.skip(f"{path.name} has no subject column")
    bad = sorted({r[col] for r in rows if r[col] and not SUBJECT_OK.match(r[col])})
    assert not bad, (
        f"{path.name}: non-pseudonymous subject labels {bad}. "
        "Expected UPN<n> / CTRL<n> / HRMDS<n>."
    )


@pytest.mark.parametrize("path", data_files(), ids=lambda p: p.name)
def test_exclusions_carry_a_reason(path: Path):
    """An excluded row carries a non-empty exclusion_reason."""
    header, rows = read_rows(path)
    if "included" not in header:
        pytest.skip(f"{path.name} has no included column")
    assert "exclusion_reason" in header, (
        f"{path.name}: has 'included' but no 'exclusion_reason' column"
    )
    missing = [
        i + 2
        for i, r in enumerate(rows)
        if str(r["included"]).strip() == "0" and not (r["exclusion_reason"] or "").strip()
    ]
    assert not missing, f"{path.name}: excluded rows without a reason at lines {missing}"


@pytest.mark.parametrize("path", data_files(), ids=lambda p: p.name)
def test_included_flag_is_binary(path: Path):
    header, rows = read_rows(path)
    if "included" not in header:
        pytest.skip(f"{path.name} has no included column")
    vals = {str(r["included"]).strip() for r in rows}
    assert vals <= {"0", "1"}, f"{path.name}: 'included' must be 0/1, found {sorted(vals)}"


@pytest.mark.parametrize("path", data_files(), ids=lambda p: p.name)
def test_no_free_text_clinical_narrative(path: Path):
    """No data field exceeds the length cap."""
    _, rows = read_rows(path)
    for i, row in enumerate(rows, start=2):
        for col, val in row.items():
            if col == "exclusion_reason" or val is None:
                continue
            assert len(str(val)) <= 120, (
                f"{path.name} line {i}, column '{col}': value is "
                f"{len(str(val))} chars. Long free text risks carrying "
                "identifying clinical narrative."
            )


def test_placeholder_modules_declare_their_status():
    """Each analysis module has a README; one without code states that."""
    mods = [p for p in (ROOT / "analysis").iterdir() if p.is_dir() and not p.name.startswith((".", "_"))]
    assert mods, "no analysis modules found"
    for m in mods:
        readme = m / "README.md"
        code = [p for p in m.rglob("*") if p.suffix in {".py", ".R", ".sh", ".Rmd"}]
        if code:
            continue  # module has code; nothing to assert about placeholder wording
        assert readme.exists(), f"analysis/{m.name}/ has neither code nor a README"
        text = readme.read_text(encoding="utf-8").lower()
        declared = any(
            phrase in text
            for phrase in ("placeholder", "awaiting data", "status: awaiting",
                           "not yet added", "code needs to be added",
                           "to be added")
        )
        assert declared, (
            f"analysis/{m.name}/README.md has no code and does not declare an "
            "incomplete status -- a reader cannot tell it is unfinished. State "
            "'placeholder' or 'awaiting data' near the top."
        )
