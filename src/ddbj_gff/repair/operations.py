"""Concrete repair operations, registered on import.

Each operation is a `detect(doc, ctx) -> list[Candidate]` (non-destructive scan)
paired with an `apply(doc, ctx, selection) -> list[Change]` (mutating; `selection`
is either a caller-provided subset of candidates or `None` to re-detect).
"""
from __future__ import annotations

from collections import defaultdict

from ..normalize.report import Change
from ..model import Feature, Span, cds_exon_intron_conflict
from .cds import coding_sequence, protein_of, has_internal_stop
from .context import RepairContext
from .registry import Operation, register
from .report import Candidate
from .partial import is_partial, partial_attrs


def _child_spans(mrna, ftype):
    return [s for c in mrna.children if c.type == ftype for s in c.spans]


# --- utr-absent-to-partial-mrna -------------------------------------------

def _detect_utr(doc, ctx: RepairContext) -> list[Candidate]:
    out: list[Candidate] = []
    for f in doc.features:
        if f.type != "mRNA" or is_partial(f):
            continue
        exon = _child_spans(f, "exon")
        cds = _child_spans(f, "CDS")
        if not exon or not cds:
            continue
        exon_lo, exon_hi = min(s.start for s in exon), max(s.end for s in exon)
        cds_lo, cds_hi = min(s.start for s in cds), max(s.end for s in cds)
        strand = exon[0].strand
        left_partial = exon_lo == cds_lo
        right_partial = exon_hi == cds_hi
        five, three = ((right_partial, left_partial) if strand == "-"
                       else (left_partial, right_partial))
        if not (five or three):
            continue
        m_lo = min(s.start for s in f.spans) if f.spans else exon_lo
        m_hi = max(s.end for s in f.spans) if f.spans else exon_hi
        ends = ", ".join(e for e, on in (("5'", five), ("3'", three)) if on)
        out.append(Candidate(
            "utr-absent-to-partial-mrna", f.id, exon[0].seqid,
            detail=f"mRNA {f.id!r} missing UTR on {ends} end -> mark partial",
            payload={"five": five, "three": three, "strand": strand,
                     "start": m_lo, "end": m_hi}))
    return out


def _apply_utr(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_utr(doc, ctx)
    changes: list[Change] = []
    for c in cands:
        f = doc.feature_index.get(c.feature_id)
        if f is None or f.type != "mRNA" or is_partial(f):
            continue
        p = c.payload
        attrs = partial_attrs(p["five"], p["three"], p["strand"], p["start"], p["end"])
        f.attributes.update(attrs)
        changes.append(Change("mark-partial", f.id or "?",
                              f"mRNA marked partial ({', '.join(sorted(attrs))})"))
    return changes


register(Operation("utr-absent-to-partial-mrna",
                   "Mark an mRNA partial on ends where a UTR is absent (structural).",
                   requires_sequence=False, detect=_detect_utr, apply=_apply_utr))


# --- missing-start-stop-to-partial-cds ------------------------------------

def _startstop_partial(f, ctx: RepairContext) -> tuple:
    """(five_prime_partial, three_prime_partial) for a CDS, judged from its sequence.

    Caller must ensure ctx.sequences is not None and f's seqid is present in it.
    Shared with the short-cds-drop-mrna exemption test so the two never drift apart.
    """
    coding, table = coding_sequence(f, ctx)
    codon_start = f.codon_start or 1
    coding_full = coding[: len(coding) - len(coding) % 3]  # whole in-frame codons only
    first = coding_full[:3]
    last = coding_full[-3:] if len(coding_full) >= 3 else ""
    five = codon_start > 1 or len(first) < 3 or first not in table.start_codons
    three = len(last) < 3 or last not in table.stop_codons
    return five, three


def _terminal_cds(cds_group, *, five_prime: bool) -> Feature:
    """The group member holding the 5' (or 3') end -- by strand, not by list order."""
    minus = cds_group[0].spans[0].strand == "-"
    want_rightmost = minus if five_prime else not minus
    if want_rightmost:
        return max(cds_group, key=lambda c: max(s.end for s in c.spans))
    return min(cds_group, key=lambda c: min(s.start for s in c.spans))


def _joined_cds(cds_group) -> Feature:
    """One synthetic CDS holding every span of *cds_group*, for a whole-CDS verdict.

    In canonical GFF a single biological CDS is often several features sharing an mRNA
    parent (AGAT gives every row its own ID), and completeness exists only across all of
    them -- it is also the unit gff2mss emits and the validator measures.

    Feature.ordered_spans() puts minus-strand spans in descending order, so extraction
    follows the biological reading direction, and codon_start reads the 5'-most *span's*
    phase -- both fall out of handing over every span. Only transl_table needs the 5'-most
    member's attributes.
    """
    spans = [s for c in cds_group for s in c.spans]
    head = _terminal_cds(cds_group, five_prime=True)
    return Feature(head.id, head.source, "CDS", spans, dict(head.attributes), [])


def _cds_groups(doc) -> list:
    """Each biological CDS as a list of features: an mRNA's CDS children, or a lone CDS.

    A CDS with no mRNA parent forms a group of one, so a document with no transcript layer
    is judged exactly as it was before grouping existed.
    """
    out, grouped = [], set()
    for m in doc.features:
        if m.type != "mRNA":
            continue
        grp = [ch for ch in m.children if ch.type == "CDS"]
        if not grp:
            continue
        grouped.update(id(c) for c in grp)
        out.append(grp)
    out.extend([f] for f in doc.features if f.type == "CDS" and id(f) not in grouped)
    return out


def _group_homogeneous(cds_group, ctx: RepairContext) -> "tuple[str, str] | None":
    """(seqid, strand) shared by every span of every member, or None if not judgeable.

    coding_sequence resolves the parent sequence from the *joined* feature's first span
    only, so a group whose members disagree on seqid or strand would have every other
    member extracted from the wrong sequence (or read backwards) with no warning. A member
    holding no span at all cannot be joined either -- there is nothing to extract. None of
    this is guessable, so the caller must skip the whole group rather than judge it, or fall
    back to judging each member alone.
    """
    if ctx.sequences is None or any(not c.spans for c in cds_group):
        return None
    seqids = {s.seqid for c in cds_group for s in c.spans}
    strands = {s.strand for c in cds_group for s in c.spans}
    if len(seqids) != 1 or len(strands) != 1:
        return None
    seqid = next(iter(seqids))
    if seqid not in ctx.sequences:
        return None
    return seqid, next(iter(strands))


def _detect_startstop(doc, ctx: RepairContext) -> list[Candidate]:
    out: list[Candidate] = []
    if ctx.sequences is None:
        return out
    # keyed by id(f): a CDS shared by two mRNAs (Parent=mA,mB) is reached from two groups,
    # each judged against a different partner set; marks accumulate across ALL groups via
    # OR on the five/three flags, so a feature that is a partial terminal in any transcript
    # ends up marked once, regardless of which mRNA is visited first.
    marks: dict = {}
    for grp in _cds_groups(doc):
        if any(is_partial(c) or c.is_trans_spliced for c in grp):
            continue
        if _group_homogeneous(grp, ctx) is None:
            continue
        five, three = _startstop_partial(_joined_cds(grp), ctx)
        if not (five or three):
            continue
        # each end is declared on the member that actually holds it; both land on the same
        # member when the group has one. seqid/strand come from that member itself, not from
        # the group or from grp[0] -- the payload must name that feature's own ends.
        for on, key in ((five, "five"), (three, "three")):
            if not on:
                continue
            f = _terminal_cds(grp, five_prime=(key == "five"))
            if not f.id:
                continue   # feature_index can't resolve an ID-less feature; nothing to mark
            entry = marks.setdefault(id(f), {"f": f, "seqid": f.spans[0].seqid,
                                             "strand": f.spans[0].strand,
                                             "five": False, "three": False})
            entry[key] = True
    for m in marks.values():
        f = m["f"]
        ends = ", ".join(e for e, on in (("5' (no start codon)", m["five"]),
                                         ("3' (no stop codon)", m["three"])) if on)
        out.append(Candidate(
            "missing-start-stop-to-partial-cds", f.id, m["seqid"],
            detail=f"CDS {f.id!r} partial: {ends}",
            payload={"five": m["five"], "three": m["three"], "strand": m["strand"],
                     "start": min(s.start for s in f.spans),
                     "end": max(s.end for s in f.spans)}))
    return out


def _apply_startstop(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_startstop(doc, ctx)
    changes: list[Change] = []
    for c in cands:
        f = doc.feature_index.get(c.feature_id)
        if f is None or f.type != "CDS" or is_partial(f):
            continue
        p = c.payload
        attrs = partial_attrs(p["five"], p["three"], p["strand"], p["start"], p["end"])
        f.attributes.update(attrs)
        changes.append(Change("mark-partial", f.id or "?",
                              f"CDS marked partial ({', '.join(sorted(attrs))})"))
    return changes


register(Operation("missing-start-stop-to-partial-cds",
                   "Mark a CDS partial when its sequence lacks a start or stop codon.",
                   requires_sequence=True, detect=_detect_startstop, apply=_apply_startstop))


# --- internal-stop-to-misc -------------------------------------------------

def _detect_internal_stop(doc, ctx: RepairContext) -> list[Candidate]:
    out: list[Candidate] = []
    if ctx.sequences is None:
        return out
    for f in doc.features:
        if f.type != "CDS" or not f.spans or f.is_trans_spliced:
            continue
        if f.spans[0].seqid not in ctx.sequences:
            continue
        if has_internal_stop(protein_of(f, ctx)):
            out.append(Candidate(
                "internal-stop-to-misc", f.id, f.spans[0].seqid,
                detail=f"CDS {f.id!r} has an internal stop codon -> misc_feature",
                payload={}))
    return out


def _apply_internal_stop(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_internal_stop(doc, ctx)
    changes: list[Change] = []
    for c in cands:
        f = doc.feature_index.get(c.feature_id)
        if f is None or f.type != "CDS":
            continue
        f.type = "misc_feature"
        note = (f"nonfunctional CDS: internal stop codon(s) detected in {f.id}; "
                f"not translated")
        f.attributes.setdefault("Note", []).append(note)
        changes.append(Change("rename-type", f.id or "?", f"CDS -> misc_feature ({note})"))
    return changes


register(Operation("internal-stop-to-misc",
                   "Retype a CDS with an internal stop codon to misc_feature (+Note).",
                   requires_sequence=True, detect=_detect_internal_stop,
                   apply=_apply_internal_stop))


# --- internal-utr-* (shared detection) ------------------------------------

def _detect_internal_utr(doc, op_name: str, remedy: str) -> list[Candidate]:
    """Candidates for mRNAs with a CDS intron that is not an exon intron.

    A CDS intron spanning an exon means a non-coding (UTR) exon sits *between*
    coding exons — a UTR interior to the CDS, which is malformed.
    """
    out: list[Candidate] = []
    for f in doc.features:
        if f.type != "mRNA":
            continue
        if f.is_trans_spliced:
            continue
        conflicts = cds_exon_intron_conflict(f)
        if not conflicts:
            continue
        lo, hi = conflicts[0]
        out.append(Candidate(
            op_name, f.id, f.spans[0].seqid if f.spans else "",
            detail=(f"mRNA {f.id!r}: {len(conflicts)} CDS intron(s) span an exon "
                    f"(internal UTR); e.g. {lo}..{hi} -> {remedy}"),
            payload={"conflicts": [[a, b] for a, b in conflicts]}))
    return out


def _reconcile_exon_spans(mrna) -> list:
    """New exon (start, end) list: the CDS splice structure with terminal UTR retained."""
    exon_segs = sorted((s.start, s.end) for s in _child_spans(mrna, "exon"))
    cds_segs = sorted((s.start, s.end) for s in _child_spans(mrna, "CDS"))
    cds_lo, cds_hi = cds_segs[0][0], cds_segs[-1][1]
    low = [(s, e) for (s, e) in exon_segs if e < cds_lo]
    high = [(s, e) for (s, e) in exon_segs if s > cds_hi]
    straddle_lo = next(((s, e) for (s, e) in exon_segs if s <= cds_lo <= e), None)
    straddle_hi = next(((s, e) for (s, e) in exon_segs if s <= cds_hi <= e), None)
    coding = []
    n = len(cds_segs)
    for i, (cs, ce) in enumerate(cds_segs):
        s, e = cs, ce
        if i == 0 and straddle_lo is not None and straddle_lo[0] < cs:
            s = straddle_lo[0]
        if i == n - 1 and straddle_hi is not None and straddle_hi[1] > ce:
            e = straddle_hi[1]
        coding.append((s, e))
    return sorted(low + coding + high)


def _remove_features(doc, feats) -> set:
    """Drop *feats* from doc.features/feature_index/roots; return their id() set."""
    ids = {id(f) for f in feats}
    for f in feats:
        if f.id and doc.feature_index.get(f.id) is f:
            del doc.feature_index[f.id]
    doc.features = [x for x in doc.features if id(x) not in ids]
    doc.roots = [x for x in doc.roots if id(x) not in ids]
    return ids


def _descendants(feat) -> list:
    """All features below *feat* (children, their children, …), depth-first, each once.

    A feature with two parents inside the subtree is reachable by more than one path. It is
    returned a single time, so a caller counting what it deleted or kept counts features
    rather than paths.
    """
    out, seen, stack = [], set(), list(feat.children)
    while stack:
        ch = stack.pop()
        if id(ch) in seen:
            continue
        seen.add(id(ch))
        out.append(ch)
        stack.extend(ch.children)
    return out


def _drop_subtree(doc, feat, *,
                  remove_childless_gene: bool = True) -> tuple[list, list, "Feature | None"]:
    """Delete *feat* and the descendants exclusive to it, keeping those shared elsewhere.

    A descendant survives when a parent of it lies outside *feat*'s subtree, or is itself a
    survivor. Every survivor is detached from *feat* **and from each deleted descendant** on
    all three links (`parent_ids`, `attributes["Parent"]`, `parents`), so nothing is left
    naming a removed feature. *feat* is detached from its gene parent, and a gene left with
    no children is removed too — unless *remove_childless_gene* is False, which the
    duplicate-isoform drop needs: it always leaves a surviving sibling, so its gene must
    never go.

    Returns `(exclusive, shared, removed_gene)` — the features deleted below *feat*, those
    kept and detached, and the gene removed as a consequence (None if none was). Callers own
    the wording of their Change messages, which differs per operation.
    """
    gene = doc.feature_index.get(feat.parent_ids[0]) if feat.parent_ids else None
    subtree = _descendants(feat)
    doomed = {id(x) for x in subtree} | {id(feat)}
    # A descendant survives when a parent of it lies outside the subtree, or is itself a
    # survivor. Grown to a fixed point so an entire branch hanging off a kept feature is
    # rescued, not only the kept feature at its top.
    keep: set[int] = set()
    while True:
        grown = {id(ch) for ch in subtree if id(ch) not in keep
                 and any(id(p) not in doomed or id(p) in keep for p in ch.parents)}
        if not grown:
            break
        keep |= grown
    shared = [ch for ch in subtree if id(ch) in keep]
    exclusive = [ch for ch in subtree if id(ch) not in keep]
    # A survivor may still name a feature that is going away -- *feat* itself, or a deleted
    # node it hung off alongside an outside parent. Drop every such link so nothing is left
    # pointing at a removed feature.
    going = [feat] + exclusive
    going_ids = {id(f) for f in going}
    # Match names only for features the index actually resolves to: the parser tolerates two
    # features sharing an ID when their types differ (it flags id-type-mismatch and keeps the
    # shadow out of feature_index), and a shadow's name must not strip its live twin's link.
    going_names = {f.id for f in going if f.id and doc.feature_index.get(f.id) is f}
    for ch in shared:
        ch.parent_ids = [p for p in ch.parent_ids if p not in going_names]
        if "Parent" in ch.attributes:
            ch.attributes["Parent"] = list(ch.parent_ids)
        ch.parents = [p for p in ch.parents if id(p) not in going_ids]
    _remove_features(doc, exclusive + [feat])
    if gene is None:
        return exclusive, shared, None
    gene.children = [ch for ch in gene.children if id(ch) != id(feat)]
    if remove_childless_gene and not gene.children:
        _remove_features(doc, [gene])
        return exclusive, shared, gene
    return exclusive, shared, None


def _detect_utr_reconcile(doc, ctx: RepairContext) -> list[Candidate]:
    return _detect_internal_utr(doc, "internal-utr-reconcile-exons",
                                "rebuild exons to CDS structure (keep terminal UTR)")


def _apply_utr_reconcile(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_utr_reconcile(doc, ctx)
    changes: list[Change] = []
    for c in cands:
        m = doc.feature_index.get(c.feature_id)
        if m is None or m.type != "mRNA" or not cds_exon_intron_conflict(m):
            continue
        new_spans = _reconcile_exon_spans(m)
        old = [ch for ch in m.children if ch.type in ("exon", "intron")]
        exclusive, shared = [], []
        for ch in old:
            (shared if [p for p in ch.parent_ids if p != m.id] else exclusive).append(ch)
        for ch in shared:                      # keep it; just detach this mRNA
            ch.parent_ids = [p for p in ch.parent_ids if p != m.id]
            if "Parent" in ch.attributes:
                ch.attributes["Parent"] = list(ch.parent_ids)
            ch.parents = [p for p in ch.parents if id(p) != id(m)]
        n_intron = sum(1 for ch in exclusive if ch.type == "intron")
        _remove_features(doc, exclusive)
        detached = {id(ch) for ch in old}
        m.children = [ch for ch in m.children if id(ch) not in detached]
        seqid, strand = m.spans[0].seqid, m.spans[0].strand
        for i, (s, e) in enumerate(new_spans, 1):
            eid = f"{m.id}.exon{i}"
            ex = Feature(eid, m.source, "exon", [Span(seqid, s, e, strand)],
                         {"ID": [eid], "Parent": [m.id]}, [m.id])
            ex.parents = [m]
            m.children.append(ex)
            doc.features.append(ex)
            doc.feature_index[eid] = ex
        m.spans = [Span(seqid, min(s for s, _ in new_spans),
                        max(e for _, e in new_spans), strand)]
        msg = (f"rebuilt {len(new_spans)} exon(s) of mRNA {m.id!r} to CDS structure "
               f"(terminal UTR kept)")
        if n_intron:
            msg += f"; dropped {n_intron} intron feature(s) (not re-derived)"
        if shared:
            msg += f"; kept {len(shared)} child feature(s) shared with another parent"
        changes.append(Change("reconcile-internal-utr", m.id or "?", msg))
    return changes


register(Operation("internal-utr-reconcile-exons",
                   "Rebuild an mRNA's exons to the CDS splice structure when a UTR sits "
                   "interior to the CDS (terminal UTR kept; CDS unchanged). "
                   "Existing intron features are dropped (not re-derived).",
                   requires_sequence=False, detect=_detect_utr_reconcile,
                   apply=_apply_utr_reconcile))


# --- internal-utr-drop-mrna (destructive) ---------------------------------

def _detect_utr_drop(doc, ctx: RepairContext) -> list[Candidate]:
    return _detect_internal_utr(doc, "internal-utr-drop-mrna",
                                "drop the mRNA and its children")


def _apply_utr_drop(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_utr_drop(doc, ctx)
    changes: list[Change] = []
    for c in cands:
        m = doc.feature_index.get(c.feature_id)
        if m is None or m.type != "mRNA" or not cds_exon_intron_conflict(m):
            continue
        exclusive, shared, dropped_gene = _drop_subtree(doc, m)
        msg = (f"dropped mRNA {m.id!r} (internal UTR) and "
               f"{len(exclusive)} child feature(s)")
        if shared:
            msg += f"; kept {len(shared)} child feature(s) shared with another parent"
        changes.append(Change("drop-internal-utr", m.id or "?", msg))
        if dropped_gene is not None:
            changes.append(Change("drop-internal-utr", dropped_gene.id or "?",
                                  f"removed gene {dropped_gene.id!r} "
                                  f"(no remaining transcripts)"))
    return changes


register(Operation("internal-utr-drop-mrna",
                   "Delete an mRNA (and its children, and a gene left childless) when a UTR "
                   "sits interior to the CDS. Destructive: excluded from --apply all.",
                   requires_sequence=False, detect=_detect_utr_drop,
                   apply=_apply_utr_drop, destructive=True))


# --- duplicate-mrna-keep-longest-cds (destructive) -------------------------

def _cds_segs(mrna) -> list:
    return sorted((s.start, s.end) for s in _child_spans(mrna, "CDS"))


def _exon_key(mrna) -> tuple:
    return tuple(sorted((s.start, s.end) for s in _child_spans(mrna, "exon")))


def _seg_introns(segs) -> set:
    return {(a[1] + 1, b[0] - 1) for a, b in zip(segs, segs[1:]) if b[0] > a[1] + 1}


def _is_five_prime_truncation(s_segs, l_segs, strand) -> bool:
    """True when *s_segs* is *l_segs* truncated at its 5' end only (a later start codon).

    Identical CDSs also satisfy this, so the byte-identical case needs no separate branch.
    A start codon that moves into a later exon drops whole segments, which is why segment
    counts are not compared. The minus strand compares the genomic start and every other
    strand value (`+`, `.`, `?`) compares the genomic end.
    """
    if not s_segs or not l_segs:
        return False
    if not all(any(ls <= s and e <= le for ls, le in l_segs) for s, e in s_segs):
        return False                                  # every coding base of S is coding in L
    if not _seg_introns(s_segs) <= _seg_introns(l_segs):
        return False                                  # S adds no splice junction
    if strand == "-":
        return min(s for s, _ in s_segs) == min(s for s, _ in l_segs)
    return max(e for _, e in s_segs) == max(e for _, e in l_segs)


def _duplicate_mrna_losers(doc) -> dict:
    """{doomed mRNA id: (keep_id, cds_len, keep_cds_len)} for conforming duplicate groups.

    A group is the mRNAs of one gene that share an exon structure. It is conforming when
    every other member's CDS is a 5' truncation of the group's longest CDS; a divergent
    group is left alone. Trans-spliced mRNAs and groups containing a CDS-less member are
    skipped.
    """
    groups = defaultdict(list)
    for f in doc.features:
        if f.type != "mRNA" or not f.id or not f.parent_ids or f.is_trans_spliced:
            continue
        key = _exon_key(f)
        if key:
            groups[(f.parent_ids[0], key)].append(f)
    out: dict = {}
    for members in groups.values():
        if len(members) < 2:
            continue
        segs = {m.id: _cds_segs(m) for m in members}
        if not all(segs[m.id] for m in members):
            continue                                  # a CDS-less member -> leave the group alone
        lengths = {m.id: sum(e - s + 1 for s, e in segs[m.id]) for m in members}
        keep = sorted(members, key=lambda m: (-lengths[m.id], m.id))[0]
        strand = keep.spans[0].strand if keep.spans else "+"
        others = [m for m in members if id(m) != id(keep)]
        if not all(_is_five_prime_truncation(segs[m.id], segs[keep.id], strand) for m in others):
            continue                                  # divergent -> not handled
        for m in others:
            out[m.id] = (keep.id, lengths[m.id], lengths[keep.id])
    return out


def _detect_duplicate_mrna(doc, ctx: RepairContext) -> list[Candidate]:
    out: list[Candidate] = []
    for mid, (keep, ln, keep_ln) in sorted(_duplicate_mrna_losers(doc).items()):
        m = doc.feature_index.get(mid)
        out.append(Candidate(
            "duplicate-mrna-keep-longest-cds", mid,
            m.spans[0].seqid if (m is not None and m.spans) else "",
            detail=(f"mRNA {mid!r} duplicates {keep!r} (identical exon structure); "
                    f"CDS {ln} bp vs {keep_ln} bp -> drop the shorter"),
            payload={"keep": keep, "cds_len": ln, "keep_cds_len": keep_ln}))
    return out


def _apply_duplicate_mrna(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_duplicate_mrna(doc, ctx)
    losers = _duplicate_mrna_losers(doc)          # judge against the document as it is now
    changes: list[Change] = []
    for c in cands:
        if c.feature_id not in losers:
            continue
        m = doc.feature_index.get(c.feature_id)
        if m is None or m.type != "mRNA":
            continue
        keep_id, ln, keep_ln = losers[c.feature_id]
        # the kept isoform is always a surviving sibling, so the gene never goes
        exclusive, shared, _ = _drop_subtree(doc, m, remove_childless_gene=False)
        msg = (f"dropped duplicate mRNA {m.id!r} (CDS {ln} bp); kept {keep_id!r} "
               f"(CDS {keep_ln} bp), {len(exclusive)} child feature(s) removed")
        if shared:
            msg += f"; kept {len(shared)} child feature(s) shared with another parent"
        changes.append(Change("drop-duplicate-mrna", m.id or "?", msg))
    return changes


register(Operation("duplicate-mrna-keep-longest-cds",
                   "Collapse mRNAs of one gene that share an exon structure down to the one "
                   "with the longest CDS (the others differ only by a later start codon). "
                   "Run after --apply all. Destructive: excluded from --apply all.",
                   requires_sequence=False, detect=_detect_duplicate_mrna,
                   apply=_apply_duplicate_mrna, destructive=True))


# --- short-intron-artificial-location (destructive) ------------------------

_SHORT_INTRON_MIN = 10
_ARTIFICIAL_LOCATION_VALUE = "low-quality sequence region"


def _has_location_exception(feat) -> bool:
    """True when the feature already declares artificial_location or ribosomal_slippage."""
    return bool(feat.attributes.get("artificial_location")
                or feat.attributes.get("ribosomal_slippage"))


def _shortest_short_intron(segs) -> "int | None":
    """Length of the shortest sub-threshold intron implied by *segs*, else None."""
    short = [e - s + 1 for s, e in _seg_introns(sorted(segs)) if e - s + 1 < _SHORT_INTRON_MIN]
    return min(short) if short else None


def _short_intron_targets(doc) -> list:
    """[(feature, shortest_intron_len)] for features needing artificial_location.

    An exon-derived short intron, or an explicit `intron` child under the threshold,
    belongs to the mRNA; a CDS-derived one belongs to every CDS child of that mRNA
    (gff2mss collapses them into one .ann CDS feature). The sources are independent.
    """
    out: list = []
    for m in doc.features:
        if m.type != "mRNA" or not m.id:
            continue
        if m.is_trans_spliced or any(c.is_trans_spliced for c in m.children):
            continue
        lens = []
        # gff2mss builds the mRNA location from exons, falling back to CDS spans when the
        # mRNA has no exon children, so mirror that fallback here.
        mrna_spans = _child_spans(m, "exon") or _child_spans(m, "CDS")
        exon_len = _shortest_short_intron([(s.start, s.end) for s in mrna_spans])
        if exon_len is not None:
            lens.append(exon_len)
        lens.extend(s.end - s.start + 1 for s in _child_spans(m, "intron")
                    if s.end - s.start + 1 < _SHORT_INTRON_MIN)
        if lens and not _has_location_exception(m):
            out.append((m, min(lens)))
        cds_len = _shortest_short_intron([(s.start, s.end) for s in _child_spans(m, "CDS")])
        if cds_len is not None:
            for c in m.children:
                if c.type == "CDS" and c.id and not _has_location_exception(c):
                    out.append((c, cds_len))
    return out


def _detect_short_intron(doc, ctx: RepairContext) -> list[Candidate]:
    out: list[Candidate] = []
    for feat, n in _short_intron_targets(doc):
        out.append(Candidate(
            "short-intron-artificial-location", feat.id,
            feat.spans[0].seqid if feat.spans else "",
            detail=(f"{feat.type} {feat.id!r} has an intron of {n} bp "
                    f"(< {_SHORT_INTRON_MIN} bp) -> "
                    f"artificial_location={_ARTIFICIAL_LOCATION_VALUE}"),
            payload={"shortest_intron": n, "value": _ARTIFICIAL_LOCATION_VALUE}))
    return out


def _apply_short_intron(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_short_intron(doc, ctx)
    targets = {id(f) for f, _ in _short_intron_targets(doc)}   # judge the document as it is now
    changes: list[Change] = []
    for c in cands:
        f = doc.feature_index.get(c.feature_id)
        if f is None or id(f) not in targets:
            continue
        f.attributes["artificial_location"] = [_ARTIFICIAL_LOCATION_VALUE]
        changes.append(Change("add-artificial-location", f.id or "?",
                              f"{f.type} {f.id!r}: artificial_location="
                              f"{_ARTIFICIAL_LOCATION_VALUE} "
                              f"(intron < {_SHORT_INTRON_MIN} bp)"))
    return changes


register(Operation("short-intron-artificial-location",
                   "Mark a CDS or mRNA containing an intron shorter than 10 bp with "
                   "artificial_location=low-quality sequence region (a submitter judgment that "
                   "the location is intentional). Run after --apply all. Destructive: excluded "
                   "from --apply all.",
                   requires_sequence=False, detect=_detect_short_intron,
                   apply=_apply_short_intron, destructive=True))


# --- partial rRNA (destructive; two mutually exclusive remedies) ------------

# Mirrors DDBJ validator rule ANN2580 (apps/ddbj/rules/annotation.py), which flags an rRNA
# when this literal substring appears in ANY of its qualifiers while the location still
# claims to be complete. Written by pybarrnap / DFAST, e.g.
# "aligned only 58.21 percent of the 18S ribosomal RNA".
_PARTIAL_RRNA_MARKER = "aligned only"


def _partial_rrna_targets(doc) -> list:
    """[(rRNA feature, the attribute value that matched)] for partial-annotated rRNAs.

    A target is an rRNA that records a partial alignment but is not yet declared partial;
    one already carrying partial=true is what the validator asks for and is left alone.
    """
    out: list = []
    for f in doc.features:
        if f.type != "rRNA" or not f.id or is_partial(f):
            continue
        hit = next((v for vals in f.attributes.values() for v in vals
                    if _PARTIAL_RRNA_MARKER in v), None)
        if hit is not None:
            out.append((f, hit))
    return out


def _detect_partial_rrna(doc, op_name: str, remedy: str) -> list[Candidate]:
    """Candidates for rRNAs recording a partial alignment while claiming a complete location."""
    out: list[Candidate] = []
    for f, hit in _partial_rrna_targets(doc):
        product = " ".join(f.attributes.get("product", [])) or f.type
        out.append(Candidate(
            op_name, f.id, f.spans[0].seqid if f.spans else "",
            detail=(f"rRNA {f.id!r} ({product}) records a partial alignment "
                    f"({hit!r}) but its location claims to be complete -> {remedy}"),
            payload={"marker": hit}))
    return out


def _detect_rrna_mark(doc, ctx: RepairContext) -> list[Candidate]:
    return _detect_partial_rrna(doc, "rrna-partial-mark-location",
                                "declare it partial (partial=true + start_range/end_range)")


def _apply_rrna_mark(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_rrna_mark(doc, ctx)
    targets = {id(f) for f, _ in _partial_rrna_targets(doc)}
    changes: list[Change] = []
    for c in cands:
        f = doc.feature_index.get(c.feature_id)
        if f is None or id(f) not in targets:
            continue
        lo = min(s.start for s in f.spans)
        hi = max(s.end for s in f.spans)
        strand = f.spans[0].strand
        # a partial alignment does not say which end is missing -> declare both
        f.attributes.update(partial_attrs(True, True, strand, lo, hi))
        changes.append(Change("mark-rrna-partial", f.id or "?",
                              f"rRNA {f.id!r} declared partial at both ends "
                              f"(partial=true, start_range, end_range)"))
    return changes


register(Operation("rrna-partial-mark-location",
                   "Declare an rRNA partial (partial=true + start_range/end_range) when it "
                   "records a partial alignment but its location claims to be complete. "
                   "Alternative to rrna-partial-drop. Run after --apply all. Destructive: "
                   "excluded from --apply all.",
                   requires_sequence=False, detect=_detect_rrna_mark,
                   apply=_apply_rrna_mark, destructive=True))


def _detect_rrna_drop(doc, ctx: RepairContext) -> list[Candidate]:
    return _detect_partial_rrna(doc, "rrna-partial-drop",
                                "drop the rRNA and its children")


def _apply_rrna_drop(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_rrna_drop(doc, ctx)
    targets = {id(f) for f, _ in _partial_rrna_targets(doc)}
    changes: list[Change] = []
    for c in cands:
        f = doc.feature_index.get(c.feature_id)
        if f is None or id(f) not in targets:
            continue
        exclusive, shared, dropped_gene = _drop_subtree(doc, f)
        msg = (f"dropped rRNA {f.id!r} (partial alignment) and "
               f"{len(exclusive)} child feature(s)")
        if shared:
            msg += f"; kept {len(shared)} child feature(s) shared with another parent"
        changes.append(Change("drop-rrna-partial", f.id or "?", msg))
        if dropped_gene is not None:
            changes.append(Change("drop-rrna-partial", dropped_gene.id or "?",
                                  f"removed gene {dropped_gene.id!r} "
                                  f"(no remaining features)"))
    return changes


register(Operation("rrna-partial-drop",
                   "Delete an rRNA (and its children, and a gene left childless) when it "
                   "records a partial alignment. Alternative to rrna-partial-mark-location. "
                   "Run after --apply all. Destructive: excluded from --apply all.",
                   requires_sequence=False, detect=_detect_rrna_drop,
                   apply=_apply_rrna_drop, destructive=True))


# --- add-inference (destructive: asserts an evidence basis) -----------------

def _split_filter(value: str) -> list[str]:
    """Comma-separated filter list -> entries, whitespace stripped, blanks dropped."""
    return [p.strip() for p in value.split(",") if p.strip()]


def _inference_targets(doc, ctx: RepairContext) -> list:
    """Features matching the type/source filters that carry no inference yet."""
    types = _split_filter(ctx.inference_types)
    sources = _split_filter(ctx.inference_sources)
    any_type = "*" in types
    out: list = []
    for f in doc.features:
        if not f.id or f.attributes.get("inference"):
            continue
        if not any_type and f.type not in types:
            continue
        if sources and f.source not in sources:
            continue
        out.append(f)
    return out


def _detect_inference(doc, ctx: RepairContext) -> list[Candidate]:
    out: list[Candidate] = []
    for f in _inference_targets(doc, ctx):
        out.append(Candidate(
            "add-inference", f.id, f.spans[0].seqid if f.spans else "",
            detail=(f"{f.type} {f.id!r} (source {f.source!r}) has no inference "
                    f"-> record {ctx.inference_value!r}"),
            payload={"value": ctx.inference_value, "source": f.source}))
    return out


def _apply_inference(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_inference(doc, ctx)
    targets = {id(f) for f in _inference_targets(doc, ctx)}
    changes: list[Change] = []
    for c in cands:
        f = doc.feature_index.get(c.feature_id)
        if f is None or id(f) not in targets:
            continue
        value = c.payload.get("value", ctx.inference_value)
        f.attributes["inference"] = [value]
        changes.append(Change("add-inference", f.id or "?",
                              f"{f.type} {f.id!r}: inference={value!r}"))
    return changes


register(Operation("add-inference",
                   "Record an evidence basis (/inference) on features selected by type and GFF "
                   "source, when they carry none. Defaults to type mRNA and the value "
                   "'ab initio prediction'; --inference-types '*' matches any type. Clears "
                   "ANN4300 for short CDSs. Run after --apply all. Destructive: excluded from "
                   "--apply all.",
                   requires_sequence=False, detect=_detect_inference,
                   apply=_apply_inference, destructive=True, detect_by_default=False))


# --- short-cds-drop-mrna (destructive) --------------------------------------

def _cds_counts_as_partial(cds_group, ctx: RepairContext) -> bool:
    """Whether a transcript's CDS is exempt from short-cds-drop-mrna as a partial CDS.

    Judged over the whole joined CDS -- the unit the converter emits and the validator
    measures -- with the same _startstop_partial the partiality operation uses, so the two
    cannot disagree. Uses the same _group_homogeneous guard as missing-start-stop-to-partial-cds,
    so a group with a spans-less member, or one that disagrees on seqid or strand, is treated
    as unjudgeable here too rather than crashing or guessing. Falls back to the attribute when
    there is no sequence for this seqid. An explicit partial=true on any member always exempts,
    even when the sequence disagrees.
    """
    if any(is_partial(c) for c in cds_group):
        return True
    if _group_homogeneous(cds_group, ctx) is None:
        return False
    return any(_startstop_partial(_joined_cds(cds_group), ctx))


def _short_cds_targets(doc, ctx: RepairContext) -> list:
    """[(mRNA, total CDS length)] for transcripts whose complete CDS is too short.

    The length the validator sees exists only at the mRNA level: in canonical GFF one CDS is
    several rows sharing an mRNA parent, which gff2mss joins into a single .ann feature.
    ANN4300 exempts partial CDSs, so an mRNA with one is left alone.
    """
    out: list = []
    for m in doc.features:
        if m.type != "mRNA" or not m.id or m.is_trans_spliced:
            continue
        cds = [ch for ch in m.children if ch.type == "CDS"]
        if not cds or any(ch.is_trans_spliced for ch in cds):
            continue
        if _cds_counts_as_partial(cds, ctx):
            continue
        total = sum(s.end - s.start + 1 for s in _child_spans(m, "CDS"))
        if total >= ctx.min_cds_length:
            continue
        out.append((m, total))
    return out


def _detect_short_cds(doc, ctx: RepairContext) -> list[Candidate]:
    out: list[Candidate] = []
    for m, total in _short_cds_targets(doc, ctx):
        out.append(Candidate(
            "short-cds-drop-mrna", m.id, m.spans[0].seqid if m.spans else "",
            detail=(f"mRNA {m.id!r}: complete CDS is {total} bp, below "
                    f"{ctx.min_cds_length} bp -> drop the transcript"),
            payload={"cds_len": total, "threshold": ctx.min_cds_length}))
    return out


def _apply_short_cds(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_short_cds(doc, ctx)
    targets = {id(m) for m, _ in _short_cds_targets(doc, ctx)}
    changes: list[Change] = []
    for c in cands:
        m = doc.feature_index.get(c.feature_id)
        if m is None or id(m) not in targets:
            continue
        exclusive, shared, dropped_gene = _drop_subtree(doc, m)
        msg = (f"dropped mRNA {m.id!r} (CDS {c.payload.get('cds_len', '?')} bp, below "
               f"{c.payload.get('threshold', ctx.min_cds_length)} bp) and "
               f"{len(exclusive)} child feature(s)")
        if shared:
            msg += f"; kept {len(shared)} child feature(s) shared with another parent"
        changes.append(Change("drop-short-cds", m.id or "?", msg))
        if dropped_gene is not None:
            changes.append(Change("drop-short-cds", dropped_gene.id or "?",
                                  f"removed gene {dropped_gene.id!r} "
                                  f"(no remaining transcripts)"))
    return changes


register(Operation("short-cds-drop-mrna",
                   "Delete a transcript whose complete CDS totals fewer than --min-cds-length "
                   "bases (default 90, per ANN4300), together with its children and a gene left "
                   "childless. Partial and trans-spliced CDSs are exempt. Alternative to recording "
                   "evidence with add-inference. Run after --apply all. Destructive: excluded from "
                   "--apply all.",
                   requires_sequence=False, detect=_detect_short_cds,
                   apply=_apply_short_cds, destructive=True))
