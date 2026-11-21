"""gff_from_tool <tool> <input> [--out] [--report]: convert one gene predictor's GFF/GTF to standard GFF3."""
import importlib
import pkgutil
import sys

from . import __path__ as _tools_path
from ._base import run_tool

TOOLS = sorted(m.name for m in pkgutil.iter_modules(_tools_path)
               if not m.name.startswith("_") and hasattr(importlib.import_module(f"{__package__}.{m.name}"), "SPEC"))


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] not in TOOLS:
        sys.stderr.write("usage: gff_from_tool <tool> <input> [--out OUT] [--report REPORT]\n"
                         f"tools: {', '.join(TOOLS)}\n")
        return 0 if argv[:1] in (["-h"], ["--help"]) else 2
    mod = importlib.import_module(f"{__package__}.{argv[0]}")
    return run_tool(mod.SPEC, argv[1:])
