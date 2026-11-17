from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class MssQualifier:
    key: str
    value: str


@dataclass
class MssFeature:
    key: str
    location: str
    qualifiers: list[MssQualifier] = field(default_factory=list)
