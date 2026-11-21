"""Parse ddbj-validator reports: validation_report.json (--json) or validation_report_details.txt lines.

    CODE:LEVEL:ENTRY:FEATURE[:QUALIFIER]:Line[N]:message
    CODE:LEVEL:ALL:ALL[:QUALIFIER]:[FILE]:message          (file-level, no line)

A report for several files has one section per file set, headed "N. NAME" (NAME = the .ann base name).
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass

_LINE = re.compile(
    r"^(?P<code>[A-Z]+\d+):(?P<level>ERR|WAR|INFO?):(?P<entry>[^:]*):(?P<feature>[^:]*):"
    r"(?:(?P<qual>[^:\[]*):)?(?:Line\[(?P<line>\d+)\]|\[FILE\]|[^:]*):(?P<message>.*)$")
_SECTION = re.compile(r"^\d+\. (?P<name>\S.*)$")


@dataclass(frozen=True)
class Finding:
    code: str
    level: str          # "ERR" | "WAR" | "INF"
    entry: str
    feature: str
    line: int | None    # 1-based line of the feature's first row in the .ann
    message: str
    file_set: str | None = None   # section of a multi-file report
    location: str | None = None   # location of the reported feature (JSON reports only)


def parse_details(text: str, file_set: str | None = None) -> list[Finding]:
    """All findings, or with `file_set` only those in that file's section.

    A report without section headers (one file set) is taken as a whole. Raises ValueError when the report
    has sections but none for `file_set`.
    """
    out: list[Finding] = []
    sections: list[str] = []
    for raw in text.splitlines():
        s = _SECTION.match(raw.strip())
        if s:
            sections.append(s.group("name").strip())
            continue
        m = _LINE.match(raw.strip())
        if not m:
            continue
        line = m.group("line")
        out.append(Finding(m.group("code"), m.group("level")[:3], m.group("entry"), m.group("feature"),
                           int(line) if line else None, m.group("message").strip(),
                           sections[-1] if sections else None))
    if file_set is None or not sections:
        return out
    if file_set not in sections:
        raise ValueError(f"report has no section for {file_set} (sections: {', '.join(sections)})")
    return [f for f in out if f.file_set == file_set]


def parse_json(text: str, file_set: str | None = None) -> list[Finding]:
    """Findings of validation_report.json, or with `file_set` only those of NAME.ann."""
    out: list[Finding] = []
    for m in json.loads(text).get("messages", []):
        name = os.path.basename(m.get("file") or "")      # "in/NAME.ann" and "NAME.ann" both occur
        name = name[:-4] if name.endswith(".ann") else name
        if file_set is not None and name != file_set:
            continue
        out.append(Finding(m["id"], m["level"][:3].upper(), m.get("entry") or "",
                           m.get("feature_type") or m.get("target") or "", m.get("line_number"),
                           m.get("message", ""), name or None, m.get("location") or None))
    return out


def parse_report(text: str, file_set: str | None = None) -> list[Finding]:
    """JSON or text report, by its first character."""
    parse = parse_json if text.lstrip().startswith("{") else parse_details
    return parse(text, file_set)


def group_by_code(findings: list[Finding]) -> dict[str, list[Finding]]:
    groups: dict[str, list[Finding]] = {}
    for f in findings:
        groups.setdefault(f.code, []).append(f)
    return groups
