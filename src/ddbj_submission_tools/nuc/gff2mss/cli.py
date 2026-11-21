from __future__ import annotations
import argparse
import os


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="gff2mss",
                                 description="Convert canonical INSDC GFF3 + FASTA to DDBJ MSS (.ann/.fasta)")
    ap.add_argument("--gff", required=True, help="canonical INSDC GFF3 (run gff_normalize first)")
    ap.add_argument("--fasta", required=True, help="genome FASTA; every sequence becomes an entry")
    ap.add_argument("--mss-config", required=True,
                    help="TOML with annotation settings: [locus_tag] prefix (required), [cds], [product], [transcript]")
    ap.add_argument("--common", required=True,
                    help="common JSON (DBLINK, SUBMITTER, REFERENCE, SOURCE, ASSEMBLY_GAP, ...), same as mss_builder")
    ap.add_argument("--sequence-roles", "--sequence_roles", dest="sequence_roles",
                    help="sequence role TSV (chromosomes, organelles, plasmids); omit for WGS. Same as mss_builder")
    ap.add_argument("--submission-category", "--submission_category", dest="submission_category", default="",
                    help="submission category (e.g. WGS, GNM); overrides _submission_category in the common JSON")
    ap.add_argument("--locus-tag-start", type=int, default=None,
                    help="override [locus_tag].start (e.g. continue organelle numbering after nuclear)")
    ap.add_argument("--out", required=True, help="output prefix; writes <out>.ann and <out>.fasta")
    args = ap.parse_args(argv)

    from ddbj_submission_tools.nuc.gff2mss.assemble import build_ann_text  # after parse_args: keeps --help fast
    from ddbj_submission_tools.nuc.gff2mss.emit import emit_fasta

    ann_text, seqs = build_ann_text(args.gff, args.fasta, args.mss_config, args.common,
                                    args.sequence_roles, args.submission_category,
                                    locus_tag_start=args.locus_tag_start)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(f"{args.out}.ann", "w", encoding="utf-8") as fh:
        fh.write(ann_text)
    with open(f"{args.out}.fasta", "w", encoding="utf-8") as fh:
        fh.write(emit_fasta(seqs))
    print(f"[gff2mss] -> {args.out}.ann / {args.out}.fasta")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
