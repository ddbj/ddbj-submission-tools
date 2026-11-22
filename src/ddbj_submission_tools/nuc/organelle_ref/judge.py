"""Verdicts: identification against the references, boundaries against the reference annotation, proposals.

The numbers are guides; the person running the check decides. Proposals are fix.edit rules,
all disabled, and only from RefSeq references.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from ddbj_submission_tools.nuc.mssann.annotation import Feature
from ddbj_submission_tools.nuc.mssann.location import LocationError, parse_location
from ddbj_submission_tools.nuc.organelle_ref.align import Hsp, q2s, s2q
from ddbj_submission_tools.nuc.organelle_ref.refs import RefSummary

FINDING_RESULTS = frozenset({"PARTIAL_MAP", "SPECIES_MISMATCH", "OTHER_SPECIES_MATCH", "LOW_IDENTITY",
                             "BOUNDARY", "CODON_START"})
KEYS = ("CDS", "rRNA", "tRNA")


@dataclass
class Finding:
    entry: str
    check: str          # identify | boundary
    result: str
    detail: str
    ref: str = ""
    proposal: str = ""


@dataclass
class Thresholds:
    min_coverage: float = 0.95
    delta: float = 1.0
    other_identity: float = 99.0
    min_identity: float = 97.0


def hit_table(hits: list[tuple[RefSummary, Hsp]], n: int = 5) -> str:
    return "; ".join(f"{r.species} {r.accession} {h.identity:.1f}% {h.qend - h.qstart + 1}/{h.qlen}"
                     for r, h in hits[:n])


def identify(entry: str, claimed: str | None, hits: list[tuple[RefSummary, Hsp]], th: Thresholds) -> list[Finding]:
    """hits sorted by identity, highest first."""
    out: list[Finding] = []
    add = lambda result, detail, ref="": out.append(Finding(entry, "identify", result, detail, ref))
    if not hits:
        add("PARTIAL_MAP", "no hit to any reference genome of the genus")
        return out
    best_cov = max(h.coverage for _, h in hits)
    if best_cov < th.min_coverage:
        add("PARTIAL_MAP", f"at most {best_cov:.0%} of the sequence aligns to a reference: {hit_table(hits)}")
    full = [(r, h) for r, h in hits if h.coverage >= th.min_coverage] or hits
    best_ref, best = full[0]
    if claimed is None:
        add("BEST_MATCH", f"species not resolved; closest: {hit_table(full)}", best_ref.accession)
        return out
    mine = [(r, h) for r, h in hits if r.species == claimed]     # even below the coverage: it is loaded
    if not mine:
        add("CLAIMED_ABSENT", f"no complete genome of {claimed}; closest: {hit_table(full)}")
        if best.identity >= th.other_identity:
            add("OTHER_SPECIES_MATCH", f"{best.identity:.1f}% to {best_ref.species} ({best_ref.accession})",
                best_ref.accession)
        return out
    my_ref, my = mine[0]
    if best_ref.species != claimed and best.identity - my.identity >= th.delta:
        add("SPECIES_MISMATCH", f"{best_ref.species} {best.identity:.1f}% vs {claimed} {my.identity:.1f}%: "
                                f"{hit_table(full)}", best_ref.accession)
    if my.identity < th.min_identity:
        add("LOW_IDENTITY", f"{my.identity:.1f}% to {claimed} ({my_ref.accession}): {hit_table(full)}",
            my_ref.accession)
    return out


def _ref_features(record, key):
    """(parts, strand, crosses_origin, mixed_strands) per reference feature of *key*; parts are 1-based."""
    for f in record.features:
        if f.type != key:
            continue
        parts = f.location.parts
        strands = {p.strand or 1 for p in parts}
        strand = strands.pop() if len(strands) == 1 else 0
        starts = [int(p.start) for p in parts]
        crosses = bool(strand) and len(parts) > 1 and starts != sorted(starts, reverse=(strand == -1))
        yield [(int(p.start) + 1, int(p.end)) for p in parts], strand, crosses, strand == 0


def _overlap(parts, a, b):
    return sum(max(0, min(e, b) - max(s, a) + 1) for s, e in parts)


def _escape_loc(loc: str) -> str:
    return "^" + re.escape(loc) + "$"


def _int(value, default=1):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


def _incomplete_stop(seq: str, loc, codon_start: int) -> bool:
    """True if a one-span CDS ends with a partial codon T or TA (completed by polyadenylation)."""
    if len(loc.spans) != 1:
        return False
    nt = seq[loc.start - 1:loc.end]
    nt = (_revcomp(nt) if loc.strand == -1 else nt)[codon_start - 1:].upper()
    rem = len(nt) % 3
    return rem in (1, 2) and nt[-rem:] in ("T", "TA")


def boundary(entry: str, features: list[Feature], seq: str, h: Hsp, ref: tuple) -> tuple[list[Finding], list[dict]]:
    """Compare CDS / rRNA / tRNA of an entry (sequence *seq*) with the reference annotation through HSP *h*.

    For a circular reference the CLI aligns to the genome plus its first part again, so HSP positions
    beyond the genome length are taken modulo the length.
    """
    summary, record = ref
    qlen, rlen = len(seq), len(record.seq)
    wrap = lambda p: (p - 1) % rlen + 1

    def from_ref(rpos):
        q = s2q(h, rpos)
        return q if q is not None else s2q(h, rpos + rlen)

    findings: list[Finding] = []
    rules: list[dict] = []
    for feat in features:
        if feat.key not in KEYS or not feat.location:
            continue
        loc_text = feat.location
        add = lambda result, detail, proposal="": findings.append(
            Finding(entry, "boundary", result, f"{feat.key} {loc_text}: {detail}", summary.accession, proposal))
        try:
            loc = parse_location(loc_text, qlen)
        except LocationError as e:
            add("NOT_COMPARED", f"location not parsed ({e})")
            continue
        if loc.strand == 0 or not loc.spans:
            add("NOT_COMPARED", "mixed strands")
            continue
        low_span = min(loc.spans, key=lambda s: s.start)
        high_span = max(loc.spans, key=lambda s: s.end)
        qlow, qhigh = loc.start, loc.end
        plow, phigh = low_span.partial_start, high_span.partial_end
        raw_low, raw_high = q2s(h, qlow), q2s(h, qhigh)
        if raw_low is None or raw_high is None:
            add("NOT_COMPARED", "feature end outside the aligned region")
            continue
        if (raw_low - 1) // rlen != (raw_high - 1) // rlen:
            add("NOT_COMPARED", "the feature spans the origin of the reference")
            continue
        slow, shigh = wrap(raw_low), wrap(raw_high)
        a, b = sorted((slow, shigh))
        cands = [(_overlap(parts, a, b), parts, strand, crosses, mixed)
                 for parts, strand, crosses, mixed in _ref_features(record, feat.key)]
        cands = [c for c in cands if c[0] > 0]
        if not cands:
            add("NOT_COMPARED", f"no {feat.key} of the reference overlaps it")
            continue
        _, parts, rstrand, crosses, mixed = max(cands, key=lambda c: c[0])
        if mixed:
            add("NOT_COMPARED", "the reference feature is on both strands (trans-spliced); compare by hand")
            continue
        if crosses:
            add("NOT_COMPARED", "the reference feature crosses the origin")
            continue
        multi = len(parts) > 1
        if multi and len(loc.spans) == 1:
            add("NOT_COMPARED", f"the reference {feat.key} has more than one part (intron); compare by hand")
            continue
        fstrand = loc.strand * h.strand
        if fstrand != rstrand:
            add("BOUNDARY", "on the other strand than the reference feature")
            continue
        rlow, rhigh = min(s for s, _ in parts), max(e for _, e in parts)
        r5, r3 = (rlow, rhigh) if rstrand == 1 else (rhigh, rlow)
        # the submitted ends in gene orientation: (query position, partial, reference position)
        if loc.strand == 1:
            e5, e3 = (qlow, plow, slow), (qhigh, phigh, shigh)
        else:
            e5, e3 = (qhigh, phigh, shigh), (qlow, plow, slow)
        have = _int(feat.get("codon_start"))
        problems = []
        for name, (qpos, partial, spos), rpos in (("5'", e5, r5), ("3'", e3, r3)):
            if partial:
                inside = from_ref(rpos)
                if inside is not None:
                    problems.append(("BOUNDARY", f"{name} end is partial, but the reference {feat.key} ends at "
                                                 f"{inside} inside the sequence"))
            elif spos != rpos:
                short = (rpos - spos) * rstrand            # > 0: the submitted CDS stops before the reference
                if feat.key == "CDS" and name == "3'" and short in (1, 2) and _incomplete_stop(seq, loc, have):
                    problems.append(("INCOMPLETE_STOP", f"3' end {short} base(s) short of the reference stop "
                                                        "(incomplete stop codon T/TA)"))
                else:
                    problems.append(("BOUNDARY", f"{name} end maps to {spos}, reference {feat.key} at {rpos}"))

        # the location to propose (RefSeq only, one span on both sides)
        new = None
        can_propose = summary.refseq and len(loc.spans) == 1 and not multi
        if can_propose and any(r == "BOUNDARY" for r, _ in problems):
            q5, q3 = from_ref(r5), from_ref(r3)
            if loc.strand == 1:
                (nl, pl), (nh, ph) = ((q5, False) if q5 else (1, True)), ((q3, False) if q3 else (qlen, True))
                new5 = (nl, pl)
            else:
                (nl, pl), (nh, ph) = ((q3, False) if q3 else (1, True)), ((q5, False) if q5 else (qlen, True))
                new5 = (nh, ph)
            new = f"{'<' if pl else ''}{nl}..{'>' if ph else ''}{nh}"
            new = f"complement({new})" if loc.strand == -1 else new
            if new == loc_text:
                new = None
        # codon_start for the final location: 1 for a complete 5' end, else from the reference frame
        codon_start = None
        if feat.key == "CDS" and not multi:
            q5pos, p5 = new5 if new else (e5[0], e5[1])
            expected = 1
            if p5:
                raw = q2s(h, q5pos)
                expected = None if raw is None else (3 - abs(wrap(raw) - r5) % 3) % 3 + 1
            if expected is not None and expected != have:
                codon_start = expected
                problems.append(("CODON_START", f"codon_start {have}, the reference frame gives {expected}"
                                                + (" for the proposed location" if new else "")))
        proposal = ""
        if can_propose and (new or codon_start is not None):
            note = f"{summary.accession} {feat.key} {rlow}..{rhigh}{'' if rstrand == 1 else ' (complement)'}"
            common = {"name": f"{entry} {feat.key} {loc_text}", "enabled": False, "entry": re.escape(entry),
                      "feature": feat.key}
            # set_value first: it selects the feature by its current location, which replace_location changes
            if codon_start is not None:
                rules.append({"action": "set_value", **common, "qualifier": "codon_start",
                              "where_location": _escape_loc(loc_text), "value": str(codon_start),
                              "add_if_missing": True, "note": note})
            if new:
                rules.append({"action": "replace_location", **common, "where_location": _escape_loc(loc_text),
                              "pattern": loc_text, "replacement": new, "note": note})
            proposal = " ".join(([new] if new else []) + ([f"codon_start {codon_start}"] if codon_start else []))
        for result, detail in problems:
            add(result, detail, proposal if result in ("BOUNDARY", "CODON_START") else "")
    return findings, rules
