from __future__ import annotations

from Bio import SeqIO

from ..io import open_text
from .context import RepairContext
from .registry import get_operation, list_operations, REGISTRY
from .report import Candidate

# internal-utr-reconcile-exons runs first: internal-stop-to-misc retypes a CDS to
# misc_feature, after which the internal-UTR anomaly is no longer detectable.
DEFAULT_ORDER = ["internal-utr-reconcile-exons", "internal-stop-to-misc",
                 "utr-absent-to-partial-mrna", "missing-start-stop-to-partial-cds"]


def apply_all_names() -> list[str]:
    """Operation names for `--apply all`: DEFAULT_ORDER minus destructive operations.

    An operation is destructive when it requires an explicit opt-in because it deletes
    annotation or asserts a submitter judgment; such operations must be named explicitly.
    """
    return [n for n in DEFAULT_ORDER if n in REGISTRY and not REGISTRY[n].destructive]


def load_sequences(path: str) -> dict:
    with open_text(path) as fh:
        return {rec.id: rec.seq for rec in SeqIO.parse(fh, "fasta")}


def run_detect(doc, ctx: RepairContext, names=None) -> list[Candidate]:
    if names:
        ops = [get_operation(n) for n in names]
    else:
        ops = [op for op in list_operations() if op.detect_by_default]
    out: list[Candidate] = []
    for op in ops:
        out.extend(op.detect(doc, ctx))
    return out


def run_apply(doc, ctx: RepairContext, names) -> list:
    changes: list = []
    for n in names:
        changes.extend(get_operation(n).apply(doc, ctx, None))
    return changes
