"""Public repository text: no internal references, README grouped by submission type."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Internal curation-spec rule IDs and versions, internal paths. Sequence data (examples/, fixtures) is not scanned.
# S-rules need context: "ribosomal protein S12" is a real product name.
PATTERN = re.compile(r"\b[BCDFM][0-9]{2}\b|(?<!protein )\bS[0-9]{2}\b|v1\.[0-9]+\b|査定|nsub-validation|work-submit|/lustre|/home/w3const")
SUFFIXES = {".py", ".md", ".json", ".tsv", ".toml"}


def iter_files():
    for root in [ROOT / "src", ROOT / "tests", ROOT / "docs"]:
        yield from (p for p in root.rglob("*") if p.suffix in SUFFIXES and "fixtures" not in p.parts)
    yield ROOT / "README.md"
    yield ROOT / "pyproject.toml"


def test_no_internal_references():
    hits = [f"{p.relative_to(ROOT)}:{n}: {m.group()}"
            for p in iter_files() if p.name != Path(__file__).name
            for n, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1)
            for m in [PATTERN.search(line)] if m]
    assert hits == []


def test_readme_groups_tools_by_submission_type():
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "塩基配列（nuc）" in text
    assert "Nucleotide sequences (nuc)" in text
