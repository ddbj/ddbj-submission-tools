"""Rule-based editing of MSS annotation entries.


Rules are read from a JSON file:

    {"rules": [ {"action": "...", ...filters..., ...parameters...}, ... ]}

Every rule may have "name" (shown in reports), "note" (free text) and
"enabled" (default true). Rules are applied in file order to every entry,
including COMMON (use "entry": "COMMON" to target COMMON only).

Filters (all optional; regular expressions, matched against the WHOLE
string unless noted):
    entry           entry name ("COMMON" for the COMMON entry)
    feature         feature key (e.g. "CDS", "CDS|tRNA", "SUBMITTER")
    qualifier       qualifier name
    where_value     regex SEARCHED in the qualifier value
    where_not_value regex SEARCHED in the qualifier value; matching values
                    are skipped (e.g. "^Holliday" to protect a proper name)
    where_location  regex SEARCHED in the feature location
    where_has       {"qualifier": NAME, "value": REGEX (searched, optional)}
                    the feature must have such a qualifier
    where_entry_has {"feature": KEY, "qualifier": NAME, "value": REGEX}
                    the ENTRY must contain a feature KEY (full match) with
                    such a qualifier (value searched, optional); e.g. set the
                    source ff_definition according to the CDS gene

Text replacement (replace_value, replace_location, rename_entry):
    pattern, replacement   literal text by default
    regex: true            pattern is a regex; replacement may use \\1

Actions:
    Qualifier level (filters: entry, feature, qualifier, where_value,
    where_not_value, where_location, where_has)
      replace_value     pattern -> replacement in matching values
      transform_value   function: lower_first_word | upper_first | lower |
                        upper | strip
      set_value         value: new value; add_if_missing: true appends the
                        qualifier to matching features that lack it
      add_qualifier     qualifier (exact name), value; if_missing: true
                        adds only when the feature has no such qualifier
      delete_qualifier  remove matching qualifiers
      rename_qualifier  new_name
    Feature level (filters: entry, feature, where_location, where_has)
      replace_location  pattern -> replacement; the result is checked
                        with parse_location and reverted if invalid
      rename_feature    new_key
      delete_feature    remove matching features
    Entry level (filter: entry)
      add_feature       key, location (optional), qualifiers [[name, value],
                        ...], position: "end" (default) | "start" |
                        "before:KEY" | "after:KEY"; if_missing: true skips
                        entries that already have the feature key
      rename_entry      pattern -> replacement on the entry name (never
                        COMMON); the .fasta entry names must be renamed too
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Pattern, Union

from .annotation import COMMON, Entry, Feature, Qualifier
from .location import LocationError, parse_location

PathLike = Union[str, Path]

QUALIFIER_ACTIONS = {"replace_value", "transform_value", "set_value",
                     "add_qualifier", "delete_qualifier", "rename_qualifier"}
FEATURE_ACTIONS = {"replace_location", "rename_feature", "delete_feature"}
ENTRY_ACTIONS = {"add_feature", "rename_entry"}
ACTIONS = QUALIFIER_ACTIONS | FEATURE_ACTIONS | ENTRY_ACTIONS

COMMON_KEYS = {"action", "name", "note", "enabled"}
FILTER_KEYS = {"entry", "feature", "qualifier", "where_value",
               "where_not_value", "where_location", "where_has",
               "where_entry_has"}
FEATURE_FILTER_KEYS = FILTER_KEYS - {"qualifier", "where_value",
                                    "where_not_value"}
TEXT_KEYS = {"pattern", "replacement", "regex"}
ALLOWED_KEYS = {
    "replace_value": FILTER_KEYS | TEXT_KEYS,
    "transform_value": FILTER_KEYS | {"function"},
    "set_value": FILTER_KEYS | {"value", "add_if_missing"},
    "add_qualifier": FILTER_KEYS | {"value", "if_missing"},
    "delete_qualifier": FILTER_KEYS,
    "rename_qualifier": FILTER_KEYS | {"new_name"},
    "replace_location": FEATURE_FILTER_KEYS | TEXT_KEYS,
    "rename_feature": FEATURE_FILTER_KEYS | {"new_key"},
    "delete_feature": FEATURE_FILTER_KEYS,
    "add_feature": {"entry", "where_entry_has", "key", "location",
                    "qualifiers", "position", "if_missing"},
    "rename_entry": {"entry"} | TEXT_KEYS,
}
REQUIRED_KEYS = {
    "replace_value": {"qualifier", "pattern", "replacement"},
    "transform_value": {"qualifier", "function"},
    "set_value": {"qualifier", "value"},
    "add_qualifier": {"qualifier", "value"},
    "delete_qualifier": {"qualifier"},
    "rename_qualifier": {"qualifier", "new_name"},
    "replace_location": {"pattern", "replacement"},
    "rename_feature": {"feature", "new_key"},
    "delete_feature": {"feature"},
    "add_feature": {"key", "qualifiers"},
    "rename_entry": {"pattern", "replacement"},
}


class RuleError(ValueError):
    """Raised for an invalid rule file."""


# ---------------------------------------------------------------------------
# Value transforms
# ---------------------------------------------------------------------------

def _lower_first_word(text: str) -> str:
    """Lower-case the first letter if the first word looks like an ordinary
    capitalised word (e.g. "Transmembrane" -> "transmembrane"). Words with
    other capitals or digits (DNA, ATP-binding, IS256, EpsG) are kept."""
    if not text:
        return text
    first = text.split(" ", 1)[0]
    if (len(first) > 1 and first[0].isupper()
            and not any(c.isupper() or c.isdigit() for c in first[1:])):
        return text[0].lower() + text[1:]
    return text


TRANSFORMS: Dict[str, Callable[[str], str]] = {
    "lower_first_word": _lower_first_word,
    "upper_first": lambda s: s[:1].upper() + s[1:],
    "lower": str.lower,
    "upper": str.upper,
    "strip": str.strip,
}


# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

@dataclass
class Change:
    """One applied change, for reports."""
    rule: int
    entry: str
    line_no: int
    what: str
    before: str
    after: str


@dataclass
class Rule:
    index: int                  # 1-based position in the rule file
    action: str
    spec: dict
    name: str
    entry_re: Optional[Pattern] = None
    feature_re: Optional[Pattern] = None
    qualifier_re: Optional[Pattern] = None
    where_value_re: Optional[Pattern] = None
    where_not_value_re: Optional[Pattern] = None
    where_location_re: Optional[Pattern] = None
    where_has: Optional[tuple] = None
    where_entry_has: Optional[tuple] = None
    count: int = 0
    examples: List["Change"] = field(default_factory=list)

    # -- matching ------------------------------------------------------------
    def match_entry(self, entry: Entry) -> bool:
        if self.entry_re is not None and not self.entry_re.fullmatch(entry.name):
            return False
        if self.where_entry_has:
            key_re, name, value_re = self.where_entry_has
            if not any(key_re.fullmatch(f.key) and any(
                    q.key == name and (value_re is None or value_re.search(q.value))
                    for q in f.qualifiers) for f in entry.features):
                return False
        return True

    def match_feature(self, feat: Feature) -> bool:
        if self.feature_re and not self.feature_re.fullmatch(feat.key):
            return False
        if self.where_location_re and not (
                feat.location and self.where_location_re.search(feat.location)):
            return False
        if self.where_has:
            name, value_re = self.where_has
            if not any(q.key == name and (value_re is None
                                          or value_re.search(q.value))
                       for q in feat.qualifiers):
                return False
        return True

    def match_qualifier(self, q: Qualifier) -> bool:
        if self.qualifier_re and not self.qualifier_re.fullmatch(q.key):
            return False
        if self.where_value_re and not self.where_value_re.search(q.value):
            return False
        if self.where_not_value_re and self.where_not_value_re.search(q.value):
            return False
        return True

    # -- text replacement ----------------------------------------------------
    def substitute(self, text: str) -> str:
        if self.spec.get("regex"):
            return re.sub(self.spec["pattern"], self.spec["replacement"], text)
        return text.replace(self.spec["pattern"], self.spec["replacement"])


def _compile(rule_no: int, key: str, pattern: Optional[str]) -> Optional[Pattern]:
    if pattern is None:
        return None
    try:
        return re.compile(pattern)
    except re.error as exc:
        raise RuleError(f"rule {rule_no}: bad regex in '{key}': {exc}") from exc


def parse_rules(data: dict) -> List[Rule]:
    """Validate a decoded rule file and return enabled rules."""
    if not isinstance(data, dict) or not isinstance(data.get("rules"), list):
        raise RuleError('rule file must be {"rules": [ ... ]}')
    rules: List[Rule] = []
    for i, spec in enumerate(data["rules"], start=1):
        if not isinstance(spec, dict):
            raise RuleError(f"rule {i}: must be an object")
        action = spec.get("action")
        if action not in ACTIONS:
            raise RuleError(f"rule {i}: unknown action {action!r} "
                            f"(choose from {', '.join(sorted(ACTIONS))})")
        unknown = set(spec) - COMMON_KEYS - ALLOWED_KEYS[action]
        if unknown:
            raise RuleError(f"rule {i} ({action}): unknown key(s) "
                            f"{', '.join(sorted(unknown))}")
        missing = REQUIRED_KEYS[action] - set(spec)
        if missing:
            raise RuleError(f"rule {i} ({action}): missing key(s) "
                            f"{', '.join(sorted(missing))}")
        if action == "transform_value" and spec["function"] not in TRANSFORMS:
            raise RuleError(f"rule {i}: unknown function {spec['function']!r} "
                            f"(choose from {', '.join(TRANSFORMS)})")
        if action == "add_feature":
            quals = spec["qualifiers"]
            if (not isinstance(quals, list) or not quals or not all(
                    isinstance(q, list) and len(q) == 2 for q in quals)):
                raise RuleError(f"rule {i}: qualifiers must be a non-empty "
                                "list of [name, value]")
        if spec.get("enabled", True) is False:
            continue

        rule = Rule(index=i, action=action, spec=spec,
                    name=spec.get("name") or f"rule {i} ({action})")
        rule.entry_re = _compile(i, "entry", spec.get("entry"))
        rule.feature_re = _compile(i, "feature", spec.get("feature"))
        if action in QUALIFIER_ACTIONS:
            rule.qualifier_re = _compile(i, "qualifier", spec.get("qualifier"))
            rule.where_value_re = _compile(i, "where_value",
                                           spec.get("where_value"))
            rule.where_not_value_re = _compile(i, "where_not_value",
                                               spec.get("where_not_value"))
        rule.where_location_re = _compile(i, "where_location",
                                          spec.get("where_location"))
        if "where_has" in spec:
            wh = spec["where_has"]
            if not isinstance(wh, dict) or "qualifier" not in wh:
                raise RuleError(f"rule {i}: where_has needs 'qualifier'")
            rule.where_has = (wh["qualifier"],
                              _compile(i, "where_has.value", wh.get("value")))
        if "where_entry_has" in spec:
            weh = spec["where_entry_has"]
            if (not isinstance(weh, dict) or "feature" not in weh
                    or "qualifier" not in weh):
                raise RuleError(f"rule {i}: where_entry_has needs 'feature' "
                                "and 'qualifier'")
            rule.where_entry_has = (
                _compile(i, "where_entry_has.feature", weh["feature"]),
                weh["qualifier"],
                _compile(i, "where_entry_has.value", weh.get("value")))
        if spec.get("regex"):
            _compile(i, "pattern", spec["pattern"])
        rules.append(rule)
    return rules


def load_rules(path: PathLike) -> List[Rule]:
    """Read and validate a JSON rule file."""
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except json.JSONDecodeError as exc:
        raise RuleError(f"{path}: invalid JSON: {exc}") from exc
    return parse_rules(data)


# ---------------------------------------------------------------------------
# Applying rules
# ---------------------------------------------------------------------------

@dataclass
class Editor:
    """Apply rules to entries.

    Counts every change per rule (Rule.count) and in total
    (total_changes), but keeps only the first max_examples changes of each
    rule (Rule.examples), so memory stays small for large batches.
    """
    rules: List[Rule]
    max_examples: int = 5
    total_changes: int = 0
    errors: List[str] = field(default_factory=list)

    def _log(self, rule: Rule, entry: Entry, line_no: int, what: str,
             before: str, after: str) -> None:
        rule.count += 1
        self.total_changes += 1
        if len(rule.examples) < self.max_examples:
            rule.examples.append(Change(rule.index, entry.name, line_no, what,
                                        before, after))

    def apply(self, entry: Entry) -> Entry:
        """Apply all rules to one entry (in place) and return it."""
        for rule in self.rules:
            if not rule.match_entry(entry):
                continue
            if rule.action in QUALIFIER_ACTIONS:
                self._apply_qualifier_rule(rule, entry)
            elif rule.action in FEATURE_ACTIONS:
                self._apply_feature_rule(rule, entry)
            else:
                self._apply_entry_rule(rule, entry)
        return entry

    # -- qualifier level -------------------------------------------------------
    def _apply_qualifier_rule(self, rule: Rule, entry: Entry) -> None:
        spec, action = rule.spec, rule.action
        for feat in entry.features:
            if not rule.match_feature(feat):
                continue
            targets = [q for q in feat.qualifiers if rule.match_qualifier(q)]

            if action == "replace_value":
                for q in targets:
                    new = rule.substitute(q.value)
                    if new != q.value:
                        self._log(rule, entry, q.line_no, f"{feat.key}/{q.key}",
                                  q.value, new)
                        q.value = new
            elif action == "transform_value":
                func = TRANSFORMS[spec["function"]]
                for q in targets:
                    new = func(q.value)
                    if new != q.value:
                        self._log(rule, entry, q.line_no, f"{feat.key}/{q.key}",
                                  q.value, new)
                        q.value = new
            elif action == "set_value":
                for q in targets:
                    if q.value != spec["value"]:
                        self._log(rule, entry, q.line_no, f"{feat.key}/{q.key}",
                                  q.value, spec["value"])
                        q.value = spec["value"]
                if not targets and spec.get("add_if_missing"):
                    self._add(rule, entry, feat, spec["qualifier"], spec["value"])
            elif action == "add_qualifier":
                exists = any(q.key == spec["qualifier"] for q in feat.qualifiers)
                if not (spec.get("if_missing") and exists):
                    self._add(rule, entry, feat, spec["qualifier"], spec["value"])
            elif action == "delete_qualifier":
                for q in targets:
                    self._log(rule, entry, q.line_no, f"{feat.key}/{q.key}",
                              q.value, "(deleted)")
                feat.qualifiers = [q for q in feat.qualifiers if q not in targets]
            elif action == "rename_qualifier":
                for q in targets:
                    self._log(rule, entry, q.line_no, f"{feat.key}/qualifier",
                              q.key, spec["new_name"])
                    q.key = spec["new_name"]

    def _add(self, rule: Rule, entry: Entry, feat: Feature, key: str,
             value: str) -> None:
        line_no = feat.qualifiers[-1].line_no if feat.qualifiers else feat.line_no
        feat.qualifiers.append(Qualifier(key=key, value=value, line_no=line_no))
        self._log(rule, entry, feat.line_no, f"{feat.key}/{key}", "(added)", value)

    # -- feature level ---------------------------------------------------------
    def _apply_feature_rule(self, rule: Rule, entry: Entry) -> None:
        spec, action = rule.spec, rule.action
        kept: List[Feature] = []
        for feat in entry.features:
            if not rule.match_feature(feat):
                kept.append(feat)
                continue
            if action == "replace_location" and feat.location is not None:
                new = rule.substitute(feat.location)
                if new != feat.location:
                    try:
                        parse_location(new)
                    except LocationError as exc:
                        self.errors.append(
                            f"{rule.name}: entry {entry.name} line "
                            f"{feat.line_no}: new location rejected ({exc}); "
                            f"kept '{feat.location}'")
                    else:
                        self._log(rule, entry, feat.line_no,
                                  f"{feat.key}/location", feat.location, new)
                        feat.location = new
            elif action == "rename_feature":
                self._log(rule, entry, feat.line_no, "feature key", feat.key,
                          spec["new_key"])
                feat.key = spec["new_key"]
            elif action == "delete_feature":
                self._log(rule, entry, feat.line_no, f"{feat.key}",
                          feat.location or "", "(deleted)")
                continue
            kept.append(feat)
        entry.features = kept

    # -- entry level -----------------------------------------------------------
    def _apply_entry_rule(self, rule: Rule, entry: Entry) -> None:
        spec, action = rule.spec, rule.action
        if action == "rename_entry":
            if entry.name == COMMON:
                return
            new = rule.substitute(entry.name)
            if new != entry.name:
                self._log(rule, entry, entry.line_no, "entry name",
                          entry.name, new)
                entry.name = new
            return

        # add_feature
        key = spec["key"]
        if spec.get("if_missing") and entry.get_feature(key) is not None:
            return
        line_no = entry.line_no
        feat = Feature(key=key, location=spec.get("location") or None,
                       line_no=line_no,
                       qualifiers=[Qualifier(k, v, line_no)
                                   for k, v in spec["qualifiers"]])
        position = spec.get("position", "end")
        index = len(entry.features)
        if position == "start":
            index = 0
        elif position.startswith(("before:", "after:")):
            where, ref = position.split(":", 1)
            found = [i for i, f in enumerate(entry.features) if f.key == ref]
            if found:
                index = found[0] if where == "before" else found[-1] + 1
        entry.features.insert(index, feat)
        self._log(rule, entry, line_no, f"{key}", "(added)",
                  "; ".join(f"{k}={v}" for k, v in spec["qualifiers"]))

    # -- helpers ---------------------------------------------------------------
    @property
    def renames_entries(self) -> bool:
        return any(r.action == "rename_entry" for r in self.rules)

    def rename_entry_name(self, name: str) -> str:
        """Apply only the rename_entry rules to a name (for .fasta headers)."""
        for rule in self.rules:
            if rule.action != "rename_entry" or name == COMMON:
                continue
            if rule.entry_re is None or rule.entry_re.fullmatch(name):
                name = rule.substitute(name)
        return name
