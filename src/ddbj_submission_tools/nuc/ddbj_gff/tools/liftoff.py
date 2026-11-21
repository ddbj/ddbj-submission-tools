"""Liftoff: gene/mRNA/exon/CDS with CDS phase '.'; compute phases."""
import sys

from ._base import ToolSpec, run_tool

SPEC = ToolSpec("liftoff", config={"compute_phase": True}, passes=("pass_compute_phase",))

if __name__ == "__main__":
    sys.exit(run_tool(SPEC))
