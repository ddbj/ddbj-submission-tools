"""ginger: parentless mRNA -> CDS models; wrap each mRNA in a gene."""
import sys

from ._base import ToolSpec, run_tool

SPEC = ToolSpec("ginger", config={"wrap_mrna_in_gene": True}, passes=("pass_wrap_mrna_in_gene",))

if __name__ == "__main__":
    sys.exit(run_tool(SPEC))
