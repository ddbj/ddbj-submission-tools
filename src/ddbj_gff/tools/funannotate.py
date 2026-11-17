"""funannotate: already gene -> mRNA -> exon/CDS; pass-through entry point for a uniform pipeline."""
import sys

from ._base import ToolSpec, run_tool

SPEC = ToolSpec("funannotate")

if __name__ == "__main__":
    sys.exit(run_tool(SPEC))
