"""Review an MSS annotation file (.ann) and report items for curation.

Usage: python -m ddbj_submission_tools.nuc.check.review FILE.ann

The checks complement the DDBJ validator. They focus on content that the
validator does not inspect: qualifier names that contain spaces, locus_tag
order, ff_definition / chromosome on source, transl_except positions,
product name wording, RNA counts and overlapping features.
"""
import argparse
import collections
import re
import sys

# Qualifier names that can legally contain a space.
SPACED_QUALIFIERS = {"sequence read archive", "Genome Coverage",
                     "Assembly Method", "Sequencing Technology",
                     "Annotation Pipeline", "Annotation Software Version",
                     "Annotation Method", "Coverage", "Assembly Name"}

# Qualifiers that are written without a value.
VALUELESS_QUALIFIERS = {"pseudo", "environmental_sample", "proviral", "focus",
                        "germline", "macronuclear", "transgenic",
                        "ribosomal_slippage", "trans_splicing", "circular_RNA"}

# Wording that INSDC / NCBI product naming guidelines discourage.
PRODUCT_SUSPECT = [
    (r"\bputative\b|\bprobable\b|\bpossible\b", "uncertainty word"),
    (r"\bhomolog(ue)?\b", "homolog"),
    (r"\bfragment\b|\bpartial\b", "fragment/partial"),
    (r"\buncharacteri[sz]ed\b|\bunknown\b", "unknown/uncharacterized"),
    (r"\.$", "ends with a period"),
    (r"\bprotein protein\b", "duplicated word"),
    (r"\bEC[ :]?\d", "EC number in product"),
    (r"^\s|\s$|\s{2,}", "extra whitespace"),
    (r"[\[\]{}]", "brackets"),
    (r"\b[A-Z]{2,}\d*_\d+\b|\bWP_\d+", "accession or locus_tag-like token"),
]

LOC_NUM = re.compile(r"\d+")


def parse(path):
    """Parse the five-column .ann file into entries and features."""
    entries = collections.OrderedDict()
    entry = feat = None
    with open(path, encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n").rstrip("\r")
            cols = (line.split("\t") + [""] * 5)[:5]
            if len(line.split("\t")) > 5:
                # More than five tab-separated columns: likely a stray tab
                entries.setdefault("_warnings", []).append(
                    ("extra columns", lineno, line))
            if cols[0]:
                entry = cols[0]
                entries.setdefault(entry, [])
            if cols[1]:
                feat = {"key": cols[1], "loc": cols[2], "q": [],
                        "line": lineno, "entry": entry}
                entries[entry].append(feat)
            if cols[3] or cols[4]:
                feat["q"].append((cols[3], cols[4], lineno))
    return entries


def span(loc):
    """Return (min, max) coordinates of a location string."""
    nums = [int(n) for n in LOC_NUM.findall(loc)]
    return (min(nums), max(nums)) if nums else (None, None)


def qvals(feat, name):
    return [v for k, v, _ in feat["q"] if k == name]


def check_transl_except(feat, out):
    """Check that each transl_except pos is an in-frame codon of the CDS."""
    loc = feat["loc"]
    comp = loc.startswith("complement(")
    if "join" in loc or "order" in loc:
        return  # multi-interval CDS: not handled here
    s, e = span(loc)
    cs = int((qvals(feat, "codon_start") or ["1"])[0])
    for v in qvals(feat, "transl_except"):
        m = re.match(r"\(pos:(complement\()?(\d+)\.\.(\d+)\)?,aa:(\w+)\)", v)
        if not m:
            out.append(("transl_except syntax", feat, v))
            continue
        p1, p2 = int(m.group(2)), int(m.group(3))
        if p2 - p1 != 2 or bool(m.group(1)) != comp:
            out.append(("transl_except not a codon on the CDS strand", feat, v))
            continue
        # Distance from the CDS 5' end, corrected by codon_start
        offset = (e - p2) if comp else (p1 - s)
        offset -= cs - 1
        if offset < 0 or offset % 3 or p1 < s or p2 > e:
            out.append(("transl_except out of frame/range", feat, v))
            continue
        codon_no = offset // 3 + 1
        for note in qvals(feat, "note"):
            n = re.search(r"position (\d+)", note)
            if n and int(n.group(1)) != codon_no:
                out.append(("transl_except/note codon number mismatch",
                            feat, f"{v} = codon {codon_no}; note: {note}"))
        out.append(("transl_except OK (info)", feat,
                    f"{v} = codon {codon_no}"))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("ann")
    args = ap.parse_args(argv)
    entries = parse(args.ann)
    findings = []
    for w in entries.pop("_warnings", []):
        print(f"[WARN] {w[0]} at line {w[1]}")

    seq_entries = [e for e in entries if e != "COMMON"]
    print("== Entries")
    chrom_like = []
    for name in seq_entries:
        feats = entries[name]
        keys = collections.Counter(f["key"] for f in feats)
        src = [f for f in feats if f["key"] == "source"]
        srcq = {k: v for k, v, _ in src[0]["q"]} if src else {}
        print(f"{name}: {dict(keys)}  source={src[0]['loc'] if src else None}")
        print(f"   source qualifiers: {sorted(srcq)}")
        if "plasmid" not in srcq and "organelle" not in srcq:
            chrom_like.append((name, srcq))
        if src and "ff_definition" not in srcq:
            findings.append(("source without ff_definition", src[0], name))
        # rRNA product tally per entry
        rrna = collections.Counter(v for f in feats if f["key"] == "rRNA"
                                   for v in qvals(f, "product"))
        print(f"   rRNA: {dict(rrna)}")
        # Features outside the source range (TOPOLOGY has no location)
        if src:
            smax = span(src[0]["loc"])[1]
            for f in feats:
                if not f["loc"]:
                    continue
                lo, hi = span(f["loc"])
                if hi and hi > smax:
                    findings.append(("location beyond source end", f, f["loc"]))
    # Several chromosome-like entries need a chromosome qualifier
    if len(chrom_like) > 1:
        for name, srcq in chrom_like:
            if "chromosome" not in srcq:
                src = [f for f in entries[name] if f["key"] == "source"]
                if src:
                    findings.append(("multiple chromosomes without "
                                     "chromosome", src[0], name))

    print("== Qualifier names")
    qnames = collections.Counter(k for e in entries.values() for f in e
                                 for k, _, _ in f["q"])
    for k, n in sorted(qnames.items()):
        print(f"   {n:6d}  {k}")
    for e in entries.values():
        for f in e:
            for k, v, ln in f["q"]:
                if " " in k and k not in SPACED_QUALIFIERS:
                    # Typically a space typed instead of a tab
                    findings.append(("qualifier name with space", f,
                                     f"line {ln}: '{k}' value='{v}'"))
                elif v == "" and f["key"] != "TOPOLOGY" \
                        and k not in VALUELESS_QUALIFIERS:
                    findings.append(("empty value", f, f"line {ln}: {k}"))

    print("== locus_tag")
    tags = [(f, v) for e in seq_entries for f in entries[e]
            for v in qvals(f, "locus_tag")]
    dup = [t for t, n in collections.Counter(v for _, v in tags).items()
           if n > 1]
    prefixes = collections.Counter(v.split("_")[0] for _, v in tags)
    print(f"   total={len(tags)} prefixes={dict(prefixes)} "
          f"duplicates={dup[:10]}")
    for key in ("CDS", "rRNA", "tRNA"):
        nums = [(f, int(re.sub(r"\D", "", v.split("_", 1)[-1]) or 0))
                for f, v in tags if f["key"] == key]
        breaks = [(a[0]["line"], a[1], b[1]) for a, b in zip(nums, nums[1:])
                  if b[1] <= a[1]]
        steps = collections.Counter(b[1] - a[1]
                                    for a, b in zip(nums, nums[1:]))
        print(f"   {key}: n={len(nums)} order breaks={breaks[:5]} "
              f"steps={steps.most_common(4)}")
    for e in seq_entries:
        for f in entries[e]:
            if f["key"] in ("CDS", "rRNA", "tRNA", "ncRNA", "tmRNA") \
                    and not qvals(f, "locus_tag"):
                findings.append(("feature without locus_tag", f, f["loc"]))

    print("== CDS")
    cds = [f for e in seq_entries for f in entries[e]
           if f["key"] == "CDS" and f["loc"]]
    qsets = collections.Counter(tuple(sorted({k for k, _, _ in f["q"]}))
                                for f in cds)
    for s, n in qsets.most_common():
        print(f"   {n:6d}  {s}")
    for f in cds:
        if not qvals(f, "product"):
            findings.append(("CDS without product", f, f["loc"]))
        if qvals(f, "transl_table") != ["11"]:
            findings.append(("transl_table not 11", f,
                             qvals(f, "transl_table")))
        lo, hi = span(f["loc"])
        if not re.search(r"join|order|[<>]", f["loc"]) and (hi - lo + 1) % 3:
            findings.append(("CDS length not a multiple of 3", f, f["loc"]))
        if qvals(f, "transl_except"):
            check_transl_except(f, findings)
    if cds:
        hyp = sum(1 for f in cds
                  if qvals(f, "product") == ["hypothetical protein"])
        print(f"   CDS={len(cds)} hypothetical={hyp} "
              f"({100 * hyp / len(cds):.1f}%)")
    infer = collections.Counter(
        re.sub(r"RefSeq:WP_\d+\.\d+", "RefSeq:WP_*", v)
        for f in cds for v in qvals(f, "inference"))
    print("   inference forms:", dict(infer))
    notes = collections.Counter(v for f in cds for v in qvals(f, "note"))
    print("   notes:", dict(notes))
    genes = [(qvals(f, "locus_tag"), qvals(f, "gene"), qvals(f, "product"))
             for f in cds if qvals(f, "gene")]
    print("   gene:", genes)

    print("== Products (suspect wording)")
    prods = collections.Counter(v for f in cds for v in qvals(f, "product"))
    for p in sorted(prods):
        for pat, why in PRODUCT_SUSPECT:
            if re.search(pat, p, re.I if "putative" in pat else 0):
                print(f"   [{why}] {prods[p]:4d}  {p}")

    print("== Duplicate / overlapping features (same entry, same key)")
    for e in seq_entries:
        bykey = collections.defaultdict(list)
        for f in entries[e]:
            if f["key"] not in ("source", "TOPOLOGY") and f["loc"]:
                bykey[f["key"]].append(f)
        for key, fs in bykey.items():
            locs = collections.Counter(f["loc"] for f in fs)
            for loc, n in locs.items():
                if n > 1:
                    findings.append(("identical location", fs[0],
                                     f"{e} {key} {loc} x{n}"))
            iv = sorted((span(f["loc"]) + (f,) for f in fs),
                        key=lambda x: x[0])
            for a, b in zip(iv, iv[1:]):
                ov = a[1] - b[0] + 1
                if ov > 150:
                    findings.append((f"{key} overlap > 150 bp", b[2],
                                     f"{e} {a[2]['loc']} / {b[2]['loc']} "
                                     f"({ov} bp)"))

    print("== Findings")
    tally = collections.Counter(k for k, _, _ in findings)
    for k, n in tally.items():
        print(f"   {n:5d}  {k}")
    for k, f, detail in findings:
        tag = (qvals(f, "locus_tag") or [""])[0]
        print(f"   - {k}: line {f['line']} {tag} {detail}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
