"""mssann: modules for DDBJ MSS submission files.

Usage (Biopython-like):
    from ddbj_submission_tools.nuc import mssann
    common = mssann.read_common("file.ann")
    for entry in mssann.parse("file.ann"):
        ...
    ann = mssann.read("file.ann")

    from ddbj_submission_tools.nuc.mssann import fasta
    for rec in fasta.parse("file.fasta"):
        print(rec.name, rec.length)

    loc = mssann.parse_location("complement(join(1..10,20..>30))")
    print(loc.strand, loc.start, loc.end, loc.spans)

    site = mssann.parse_location("4618897^1", seq_length=4618897, circular=True)

    seq = fasta.read_sequence("file.fasta", "chr01")
    cds_seq = loc.extract(seq)
"""

from . import edit, fasta, location, rrna, seqfeature, strain, trna, writer
from .location import (Location, LocationError, Span, parse_location,
                       reverse_complement)
from .annotation import (
    COMMON,
    DDBJ_FEATURES,
    PSEUDOGENE_VALUES,
    AnnParseError,
    AnnotationFile,
    Entry,
    Feature,
    Qualifier,
    effective_source,
    find_entry,
    is_circular,
    parse,
    read,
    read_common,
    split_row,
    unquote,
)


__all__ = [
    "COMMON", "DDBJ_FEATURES", "PSEUDOGENE_VALUES", "AnnParseError", "AnnotationFile", "Entry",
    "Feature", "Location", "LocationError", "Qualifier", "Span",
    "effective_source", "fasta", "find_entry", "is_circular", "location",
    "parse_location", "reverse_complement", "edit", "rrna", "seqfeature", "strain", "trna", "writer",
    "parse", "read", "read_common", "split_row", "unquote",
]
