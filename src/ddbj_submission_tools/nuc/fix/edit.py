"""Edit MSS annotation files (.ann) in batch with a JSON rule file.

Rules can change anything in the file: the COMMON block, entry names,
feature keys, locations, qualifier names and qualifier values (see
mssann/edit.py for the rule syntax and rules/edit_rules_example.json).

Original files are never changed unless --in-place is given. Choose one:
    --dry-run       show what would change; write nothing
    --outdir DIR    write edited files to DIR (same file names)
    --in-place      replace the files, keeping FILE.orig backups

If the rules rename entries, the matching .fasta (same base name) gets
the same renames, so .ann and .fasta stay consistent.

Usage:
    python -m ddbj_submission_tools.nuc.fix.edit RULES.json FILE.ann [FILE.ann ...] --dry-run
    python -m ddbj_submission_tools.nuc.fix.edit RULES.json DIR/*.ann --outdir fixed/
    python -m ddbj_submission_tools.nuc.fix.edit RULES.json FILE.ann --in-place
"""

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path

from ddbj_submission_tools.nuc import mssann
from ddbj_submission_tools.nuc.mssann import edit, fasta
from ddbj_submission_tools.nuc.mssann.writer import WriteError, write_entries


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Edit .ann files in batch with a JSON rule file.")
    parser.add_argument("rules", type=Path, help="JSON rule file")
    parser.add_argument("ann", type=Path, nargs="+", help="annotation file(s)")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true",
                      help="report changes only")
    mode.add_argument("--outdir", type=Path,
                      help="write edited files to this directory")
    mode.add_argument("--in-place", action="store_true",
                      help="overwrite files, keeping FILE.orig backups")
    parser.add_argument("--keep-quotes", action="store_true",
                        help='write spreadsheet quoting "..." back '
                             "(default: quotes are removed)")
    parser.add_argument("--show", type=int, default=5,
                        help="example changes shown per rule (default 5)")
    return parser.parse_args(argv)


def edit_file(path, editor, out_path, keep_quotes):
    """Parse, edit and write one file. Returns (n_quoted, dropped)."""
    seen_names = set()
    n_quoted = 0

    def edited_entries():
        nonlocal n_quoted
        for entry in mssann.parse(path):
            for feat in entry.features:
                n_quoted += feat.location_quoted + sum(q.quoted for q in feat.qualifiers)
            editor.apply(entry)
            if entry.name in seen_names:
                raise WriteError(f"duplicate entry name '{entry.name}' after "
                                 "editing")
            seen_names.add(entry.name)
            yield entry

    if out_path is None:                       # dry run: consume, write nothing
        dropped = []
        for entry in edited_entries():
            dropped += [(entry, f) for f in entry.features if not f.qualifiers]
        return n_quoted, dropped
    with open(out_path, "w", encoding="utf-8", newline="\n") as fh:
        dropped = write_entries(edited_entries(), fh, keep_quotes=keep_quotes)
    return n_quoted, dropped


def main(argv=None):
    args = parse_args(argv)
    try:
        rules = edit.load_rules(args.rules)
    except (OSError, edit.RuleError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    if not rules:
        print("Error: no enabled rules", file=sys.stderr)
        return 2
    if args.outdir:
        args.outdir.mkdir(parents=True, exist_ok=True)
        if any(p.resolve().parent == args.outdir.resolve() for p in args.ann):
            print("Error: --outdir must differ from the input directory",
                  file=sys.stderr)
            return 2

    editor = edit.Editor(rules, max_examples=args.show)
    n_failed = 0
    for path in args.ann:
        before = editor.total_changes
        fasta_in = path.with_suffix(".fasta")
        if editor.renames_entries and not fasta_in.is_file():
            print(f"Error: {path}: rules rename entries but {fasta_in} is "
                  "missing; skipped", file=sys.stderr)
            n_failed += 1
            continue

        # Write to a temporary file first; move it only on success
        tmp = None
        if not args.dry_run:
            target_dir = args.outdir if args.outdir else path.parent
            fd, tmp = tempfile.mkstemp(suffix=".ann.tmp", dir=target_dir)
            os.close(fd)
        try:
            n_quoted, dropped = edit_file(path, editor, tmp, args.keep_quotes)
        except (mssann.AnnParseError, WriteError) as exc:
            print(f"Error: {path}: {exc}; file not written", file=sys.stderr)
            if tmp:
                os.unlink(tmp)
            n_failed += 1
            continue

        for entry, feat in dropped:
            print(f"Warning: {path}: entry {entry.name} line {feat.line_no}: "
                  f"{feat.key} has no qualifier left and was removed",
                  file=sys.stderr)

        if not args.dry_run:
            if args.in_place:
                shutil.copy2(path, path.with_suffix(".ann.orig"))
                os.replace(tmp, path)
            else:
                os.replace(tmp, args.outdir / path.name)
            if editor.renames_entries:
                fasta_out = (fasta_in if args.in_place
                             else args.outdir / fasta_in.name)
                if args.in_place:
                    shutil.copy2(fasta_in, fasta_in.with_suffix(".fasta.orig"))
                    tmp_f = fasta_in.with_suffix(".fasta.tmp")
                    fasta.rename_entries(fasta_in, tmp_f,
                                         editor.rename_entry_name)
                    os.replace(tmp_f, fasta_in)
                else:
                    fasta.rename_entries(fasta_in, fasta_out,
                                         editor.rename_entry_name)

        n_changes = editor.total_changes - before
        quoted = (f", {n_quoted} quoted field(s) "
                  f"{'kept' if args.keep_quotes else 'unquoted'}"
                  if n_quoted else "")
        print(f"{path}: {n_changes} change(s){quoted}", file=sys.stderr)

    # Summary per rule with examples
    print("\n# Summary per rule", file=sys.stderr)
    for rule in editor.rules:
        print(f"[{rule.index}] {rule.name}: {rule.count} change(s)",
              file=sys.stderr)
        for c in rule.examples:
            print(f"    {c.entry} line {c.line_no} {c.what}: "
                  f"'{c.before}' -> '{c.after}'", file=sys.stderr)
        if rule.count > len(rule.examples):
            print(f"    ... {rule.count - len(rule.examples)} more",
                  file=sys.stderr)
    for message in editor.errors:
        print(f"Warning: {message}", file=sys.stderr)

    mode = ("dry run (nothing written)" if args.dry_run else
            "in place (backups *.orig)" if args.in_place else
            f"written to {args.outdir}")
    print(f"# {len(args.ann) - n_failed}/{len(args.ann)} file(s) processed, "
          f"{editor.total_changes} change(s), {mode}", file=sys.stderr)
    return 1 if n_failed else 0


if __name__ == "__main__":
    sys.exit(main())
