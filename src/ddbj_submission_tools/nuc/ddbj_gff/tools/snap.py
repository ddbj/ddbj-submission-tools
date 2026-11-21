"""SNAP (zff2gff3.pl output): 8-column CDS rows grouped by Name; build gene/mRNA/CDS and compute phases."""
import sys

from ._base import ToolSpec, run_tool

SPEC = ToolSpec("snap", reader="snap", config={"compute_phase": True}, passes=("pass_compute_phase",))

if __name__ == "__main__":
    sys.exit(run_tool(SPEC))
