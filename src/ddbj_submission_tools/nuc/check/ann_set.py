"""Summarize a set of MSS annotation files that belong to one submission.

A submission may arrive as many .ann/.fasta pairs (one entry per pair).
This script reads them all and reports what differs between the files and
what has to be checked before the entries are merged into one annotation
file:

  * COMMON blocks: which qualifiers are identical in every file and which
    ones differ (DBLINK, ST_COMMENT ...). Only the identical ones may stay
    in COMMON after a merge; the others have to be written per entry.
  * source features: qualifier sets, organism / isolate / entry-name
    agreement, sequence length against the source location.
  * features: counts per key, partial locations (< or >), features without
    locus_tag.
  * products: every distinct product value with its count, and the ones
    whose wording is questionable for INSDC (ambiguous "or", a database
    name, an uncertainty word ...).
  * tRNA: length, product / gene wording, pseudogene qualifiers.

Usage:
    python -m ddbj_submission_tools.nuc.check.ann_set DIR [--products] [--trna]
    python -m ddbj_submission_tools.nuc.check.ann_set DIR/*.ann
"""

import argparse
import collections
import re
import sys
from pathlib import Path

LOC_NUM = re.compile(r"\d+")
# Product wording that DDBJ asks submitters to avoid.
PRODUCT_SUSPECT = [
    (r"\bor\b", "ambiguous 'or'"),
    (r"\b(VFDB|CARD|PHROG|phrog|Pfam|UniProt)\b", "database name"),
    (r"\bputative\b|\bprobable\b|\bpossible\b", "uncertainty word"),
    (r"\bunknown\b|\buncharacteri[sz]ed\b", "unknown"),
    (r"\bdomain\b|\bprotein protein\b", "wording"),
    (r"^\s|\s$|\s{2,}", "extra whitespace"),
    (r"\.$", "ends with a period"),
]


def parse(path):
    """Parse a five-column .ann file into (common, entries).

    common  -- list of (feature key, qualifier, value) of the COMMON entry
    entries -- OrderedDict: entry name -> list of features, each feature a
               dict with key, loc and q (list of (qualifier, value)).
    """
    common = []
    entries = collections.OrderedDict()
    entry = feat = None
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            cols = (raw.rstrip("\n").rstrip("\r").split("\t") + [""] * 5)[:5]
            if cols[0]:
                entry = cols[0]
                if entry != "COMMON":
                    entries.setdefault(entry, [])
            if cols[1]:
                feat = {"key": cols[1], "loc": cols[2], "q": []}
                if entry != "COMMON":
                    entries[entry].append(feat)
            if cols[3] or cols[4]:
                if entry == "COMMON":
                    common.append((feat["key"] if feat else "", cols[3],
                                   cols[4]))
                elif feat is not None:
                    feat["q"].append((cols[3], cols[4]))
    return common, entries


def seq_length(ann_path):
    """Length of the sequence in the .fasta file next to the .ann file."""
    fasta = ann_path.with_suffix(".fasta")
    if not fasta.is_file():
        return None
    n = 0
    with open(fasta, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">") or line == "//":
                continue
            n += len(line)
    return n


def qval(feat, name):
    return [v for k, v in feat["q"] if k == name]


def span(loc):
    nums = [int(n) for n in LOC_NUM.findall(loc)]
    return (min(nums), max(nums)) if nums else (None, None)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("path", nargs="+", type=Path,
                    help="a directory of .ann files, or the files themselves")
    ap.add_argument("--products", action="store_true",
                    help="list every distinct product value")
    ap.add_argument("--trna", action="store_true",
                    help="list every tRNA feature")
    args = ap.parse_args(argv)

    files = []
    for p in args.path:
        files += sorted(p.glob("*.ann")) if p.is_dir() else [p]
    if not files:
        print("No .ann file found", file=sys.stderr)
        return 2

    common_values = collections.OrderedDict()      # (feature, qual) -> values
    source_qsets = collections.Counter()
    products = collections.Counter()
    feature_keys = collections.Counter()
    trnas = []
    notes = collections.Counter()
    functions = collections.Counter()
    problems = []
    n_entries = 0

    for path in files:
        common, entries = parse(path)
        for feature, qual, value in common:
            common_values.setdefault((feature, qual), set()).add(value)
        length = seq_length(path)
        for name, feats in entries.items():
            n_entries += 1
            if name != path.stem:
                problems.append(f"{path.name}: entry name '{name}' differs "
                                f"from the file name")
            src = next((f for f in feats if f["key"] == "source"), None)
            if src is None:
                problems.append(f"{path.name}: no source feature")
                continue
            source_qsets[tuple(k for k, _ in src["q"])] += 1
            if length is not None and span(src["loc"])[1] != length:
                problems.append(f"{path.name}: source {src['loc']} but the "
                                f"sequence is {length} nt")
            organism = (qval(src, "organism") or [""])[0]
            isolate = (qval(src, "isolate") or qval(src, "strain") or [""])[0]
            if isolate and not organism.endswith(isolate):
                problems.append(f"{path.name}: organism '{organism}' does not "
                                f"end with '{isolate}'")
            for f in feats:
                feature_keys[f["key"]] += 1
                if f["key"] in ("CDS", "tRNA", "rRNA", "ncRNA", "tmRNA") \
                        and not qval(f, "locus_tag"):
                    problems.append(f"{path.name}: {f['key']} {f['loc']} "
                                    "without locus_tag")
                for v in qval(f, "product"):
                    products[v] += 1
                for v in qval(f, "note"):
                    notes[v] += 1
                for v in qval(f, "function"):
                    functions[v] += 1
                if "<" in f["loc"] or ">" in f["loc"]:
                    problems.append(f"{path.name}: partial location "
                                    f"{f['key']} {f['loc']}")
                if f["key"] == "tRNA":
                    lo, hi = span(f["loc"])
                    trnas.append((path.stem, f["loc"], hi - lo + 1,
                                  (qval(f, "product") or [""])[0],
                                  (qval(f, "gene") or [""])[0],
                                  (qval(f, "pseudogene") or [""])[0]))

    print(f"== {len(files)} files, {n_entries} entries")

    print("\n== COMMON")
    for (feature, qual), values in common_values.items():
        if len(values) == 1:
            print(f"   same   {feature:12s} {qual:24s} {next(iter(values))}")
        else:
            shown = ", ".join(sorted(values)[:4])
            more = " ..." if len(values) > 4 else ""
            print(f"   DIFFER {feature:12s} {qual:24s} "
                  f"{len(values)} values: {shown}{more}")
    print("   (only the 'same' lines may stay in COMMON after a merge)")

    print("\n== source qualifier sets")
    for qset, n in source_qsets.most_common():
        print(f"   {n:5d}  {', '.join(qset)}")

    print("\n== feature keys")
    for key, n in feature_keys.most_common():
        print(f"   {n:7d}  {key}")

    print("\n== products: "
          f"{len(products)} distinct, {sum(products.values())} total")
    print(f"   hypothetical protein: {products.get('hypothetical protein', 0)}")
    for value, n in sorted(products.items()):
        why = [w for pat, w in PRODUCT_SUSPECT
               if re.search(pat, value, re.IGNORECASE)]
        if why:
            print(f"   [{', '.join(why)}] {n:5d}  {value}")
        elif args.products:
            print(f"   {n:5d}  {value}")

    print("\n== function values")
    for value, n in functions.most_common():
        print(f"   {n:7d}  {value}")

    print("\n== note values (top 10)")
    for value, n in notes.most_common(10):
        print(f"   {n:7d}  {value}")

    print(f"\n== tRNA: {len(trnas)} features")
    for stem, loc, length, product, gene, pseudo in trnas:
        odd = length < 50 or length > 150 or pseudo or "OTHER" in product
        if odd or args.trna:
            print(f"   {stem:20s} {loc:28s} {length:4d} nt  "
                  f"product={product!r} gene={gene!r} pseudogene={pseudo!r}")

    print(f"\n== Other points to check: {len(problems)}")
    kinds = collections.Counter(re.sub(r"^.*?: ", "", p).split(" ")[0]
                                for p in problems)
    for kind, n in kinds.most_common():
        print(f"   {n:7d}  {kind} ...")
    for p in problems[:200]:
        print(f"   - {p}")
    if len(problems) > 200:
        print(f"   ... {len(problems) - 200} more")
    return 0


if __name__ == "__main__":
    sys.exit(main())
