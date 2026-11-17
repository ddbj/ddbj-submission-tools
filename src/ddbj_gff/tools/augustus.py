"""augustus: gene/transcript/CDS with start_codon/stop_codon/intron rows (GFF3 or GTF);
stop codon may sit outside the CDS; per-chromosome runs may reuse IDs."""
import sys

from ._base import ToolSpec, run_tool

SPEC = ToolSpec(
    "augustus",
    dedupe_ids=True,
    config={"drop_types": ("intron", "start_codon", "stop_codon"),
            "extend_cds_by_stop_codon": True, "coerce_transcript_to_mrna": True},
    passes=("pass_extend_cds_by_stop_codon", "pass_drop_types", "pass_coerce_transcript_to_mrna"),
)

if __name__ == "__main__":
    sys.exit(run_tool(SPEC))
