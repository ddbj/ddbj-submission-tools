from __future__ import annotations


class LocusTagAssigner:
    def __init__(self, prefix: str, width: int = 6, start: int = 10, step: int = 10):
        self._prefix = prefix
        self._width = width
        self._next = start
        self._step = step
        self._given: dict[int, tuple] = {}   # id(feature) -> (feature, tag); holding feature keeps its id unique

    @classmethod
    def from_config(cls, cfg) -> "LocusTagAssigner":
        return cls(cfg.locus_tag_prefix, cfg.locus_tag_width,
                   cfg.locus_tag_start, cfg.locus_tag_step)

    def assign(self, feature) -> str:
        existing = feature.locus_tag
        if existing:
            return existing
        if id(feature) not in self._given:
            self._given[id(feature)] = (feature, f"{self._prefix}_{self._next:0{self._width}d}")
            self._next += self._step
        return self._given[id(feature)][1]

    def tag_of(self, feature) -> str | None:
        """The tag already given to (or carried by) feature; None if it has none yet."""
        return feature.locus_tag or self._given.get(id(feature), (None, None))[1]
