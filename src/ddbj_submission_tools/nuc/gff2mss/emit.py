from __future__ import annotations

from .model import MssFeature, MssQualifier


def feature_rows(feat: MssFeature) -> list[list[str]]:
    quals = feat.qualifiers or [MssQualifier("", "")]
    rows: list[list[str]] = []
    for i, q in enumerate(quals):
        col2 = feat.key if i == 0 else ""
        col3 = feat.location if i == 0 else ""
        rows.append(["", col2, col3, q.key, q.value])
    return rows


def emit_fasta(seqs: dict[str, object]) -> str:
    out: list[str] = []
    for name, seq in seqs.items():
        out.append(f">{name}")
        s = str(seq)
        for i in range(0, len(s), 60):
            out.append(s[i : i + 60])
        out.append("//")
    return "\n".join(out) + "\n"
