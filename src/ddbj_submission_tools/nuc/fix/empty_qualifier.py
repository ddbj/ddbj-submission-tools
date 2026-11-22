"""Put a temporary qualifier on .ann feature lines that have none.

A feature line such as

    <tab>5'UTR<tab>1..211<tab><tab>

carries a Location but no Qualifier at all. This is a correct DDBJ
annotation: 5'UTR, 3'UTR, D-loop, centromere, telomere and other
region-only features are often written without any qualifier, and the
ANN3170 warning of the DDBJ validator ('No qualifier found after the
location column') is not something to fix.

The mssann parser, however, still rejects such a line ('empty Qualifier
column'), so nuc.fix.edit cannot be run on the file. Until the parser is
fixed, this script works around that on the raw text: put a temporary
qualifier in with --set, run nuc.fix.edit, then take it out again with
--clear so that the file goes back to its submitted form.

  --set FEATURE:QUALIFIER=VALUE     fill the empty Qualifier and Value
                                    columns (FEATURE may be '*')
  --clear FEATURE:QUALIFIER=VALUE   empty those columns again, but only
                                    where the pair matches exactly
  --delete-unmatched                drop the whole line instead

Every other line is written out unchanged. Line endings become LF.

Usage:
    python -m ddbj_submission_tools.nuc.fix.empty_qualifier FILE.ann [...] --dry-run
    python -m ddbj_submission_tools.nuc.fix.empty_qualifier FILE.ann \\
        --set "5'UTR:note=TMP" --set "3'UTR:note=TMP" --in-place
    python -m ddbj_submission_tools.nuc.fix.empty_qualifier FILE.ann \\
        --clear "5'UTR:note=TMP" --clear "3'UTR:note=TMP" --in-place

Exit code: 0 if no problem remains, 1 if some line could not be
repaired, 2 on a usage or file error.
"""

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path

NUM_COLUMNS = 5
BACKUP_SUFFIX = ".orig"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Put a temporary qualifier on .ann feature lines "
                    "that have none, or take it out again.")
    parser.add_argument("ann", type=Path, nargs="+",
                        help="annotation file(s) (.ann)")
    parser.add_argument("--set", dest="sets", action="append", default=[],
                        metavar="FEATURE:QUALIFIER=VALUE",
                        help="qualifier to put on a feature with no "
                             "qualifier; FEATURE may be '*' (repeatable)")
    parser.add_argument("--clear", dest="clears", action="append",
                        default=[], metavar="FEATURE:QUALIFIER=VALUE",
                        help="empty the Qualifier and Value columns of a "
                             "feature line carrying exactly this pair, "
                             "undoing --set (repeatable)")
    parser.add_argument("--delete-unmatched", action="store_true",
                        help="delete lines that no --set rule matched")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true",
                      help="report only; write nothing")
    mode.add_argument("--outdir", type=Path, default=None,
                      help="write the repaired copies to this directory")
    mode.add_argument("--in-place", action="store_true",
                      help="overwrite the file, keeping FILE.ann.orig")
    return parser.parse_args(argv)


def parse_sets(specs):
    """Turn --set strings into {feature_key: (qualifier, value)}."""
    table = {}
    for spec in specs:
        if ":" not in spec or "=" not in spec.split(":", 1)[1]:
            raise ValueError(f"--set must be FEATURE:QUALIFIER=VALUE: {spec}")
        feature, rest = spec.split(":", 1)
        qualifier, value = rest.split("=", 1)
        if not feature or not qualifier:
            raise ValueError(f"empty FEATURE or QUALIFIER in --set: {spec}")
        if feature in table:
            raise ValueError(f"--set given twice for feature '{feature}'")
        table[feature] = (qualifier, value)
    return table


def split_row(line):
    """Split a line into exactly 5 columns (pad short rows, keep extras)."""
    cols = line.rstrip("\r\n").split("\t")
    if len(cols) < NUM_COLUMNS:
        cols += [""] * (NUM_COLUMNS - len(cols))
    return cols


def repair_lines(path, lines, table, clear_table, delete_unmatched):
    """Return (new_lines, changes, unfixed).

    changes is a list of (line_no, action, text) for the report; unfixed
    counts the lines that stayed broken.
    """
    new_lines = []
    changes = []
    unfixed = 0
    for line_no, line in enumerate(lines, start=1):
        if line.strip() == "":
            new_lines.append(line.rstrip("\r\n") + "\n")
            continue
        cols = split_row(line)
        feature, qualifier = cols[1], cols[3]

        if feature != "" and qualifier != "":
            rule = clear_table.get(feature) or clear_table.get("*")
            if rule is not None and (qualifier, cols[4]) == rule:
                cols[3] = cols[4] = ""
                changes.append((line_no, "clear",
                                f"{feature} {cols[2]} -> (no qualifier)"))
            new_lines.append("\t".join(cols[:NUM_COLUMNS]) + "\n")
            continue

        if feature == "" or qualifier != "":
            new_lines.append("\t".join(cols[:NUM_COLUMNS]) + "\n")
            continue

        rule = table.get(feature) or table.get("*")
        if rule is not None:
            cols[3], cols[4] = rule
            new_lines.append("\t".join(cols[:NUM_COLUMNS]) + "\n")
            changes.append((line_no, "set",
                            f"{feature} {cols[2]} -> {rule[0]} {rule[1]}"))
        elif delete_unmatched:
            changes.append((line_no, "delete", f"{feature} {cols[2]}"))
        else:
            new_lines.append("\t".join(cols[:NUM_COLUMNS]) + "\n")
            if not table:
                continue       # a --clear-only run leaves these alone
            unfixed += 1
            print(f"Warning: {path}:{line_no}: feature '{feature}' has no "
                  "qualifier and no --set rule matches", file=sys.stderr)
    return new_lines, changes, unfixed


def write_out(path, new_lines, args):
    """Write the repaired file, via a temporary file in the same directory."""
    if args.outdir is not None:
        args.outdir.mkdir(parents=True, exist_ok=True)
        target = args.outdir / path.name
        if target.resolve() == path.resolve():
            print(f"Error: --outdir is the input directory: {args.outdir}",
                  file=sys.stderr)
            return False
    else:
        target = path
        shutil.copy2(path, str(path) + BACKUP_SUFFIX)

    fd, tmp = tempfile.mkstemp(dir=str(target.parent), suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            fh.writelines(new_lines)
        os.replace(tmp, target)
    except Exception:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    return True


def main(argv=None):
    args = parse_args(argv)
    try:
        table = parse_sets(args.sets)
        clear_table = parse_sets(args.clears)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    total_changes = 0
    total_unfixed = 0
    for path in args.ann:
        if not path.is_file():
            print(f"Error: file not found: {path}", file=sys.stderr)
            return 2
        with open(path, encoding="utf-8") as fh:
            lines = fh.readlines()

        new_lines, changes, unfixed = repair_lines(
            path, lines, table, clear_table, args.delete_unmatched)
        for line_no, action, text in changes:
            print(f"{path}:{line_no}: {action}: {text}", file=sys.stderr)
        if unfixed:
            print(f"Warning: {path}: {unfixed} line(s) still have no "
                  "qualifier", file=sys.stderr)

        total_changes += len(changes)
        total_unfixed += unfixed

        if not args.dry_run and changes:
            if not write_out(path, new_lines, args):
                return 2

    verb = "would change" if args.dry_run else "changed"
    print(f"# {verb} {total_changes} line(s) in {len(args.ann)} file(s), "
          f"{total_unfixed} left unfixed", file=sys.stderr)
    return 1 if total_unfixed else 0


if __name__ == "__main__":
    sys.exit(main())
