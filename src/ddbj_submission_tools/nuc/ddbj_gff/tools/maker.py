"""MAKER: gene models plus evidence alignments (match/match_part/...) and contig rows; drop the evidence."""
import sys

from ._base import ToolSpec, run_tool

SPEC = ToolSpec(
    "maker",
    config={"drop_types": ("match", "match_part", "protein_match", "expressed_sequence_match", "contig")},
    passes=("pass_drop_types",),
)

if __name__ == "__main__":
    sys.exit(run_tool(SPEC))
