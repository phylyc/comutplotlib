"""Config-driven classification of free-text SIF fields into category codes.

The rule tables live in ``comutplotlib/config/sif_classification.json`` (or a
user override, see :mod:`comutplotlib.config`). This module only interprets
them; the institution/study vocabulary itself is data, not code.
"""

from comutplotlib.config import load_config

_RULES = load_config("sif_classification.json", default={})


def _matches(rule: dict, fields: dict) -> bool:
    op = rule["op"]
    if op == "default":
        return True
    if op == "all":
        return all(_matches(r, fields) for r in rule["conditions"])
    if op == "any":
        return any(_matches(r, fields) for r in rule["conditions"])

    value = fields.get(rule["field"])
    if op == "contains":
        return isinstance(value, str) and any(p in value for p in rule["values"])
    if op == "equals":
        if "values" in rule:
            return value in rule["values"]
        return value == rule["value"]
    if op == "na_or_in":
        # mirrors the historical ``isinstance(x, float) or x in [...]`` guard
        return isinstance(value, float) or value in rule.get("values", [])
    raise ValueError(f"Unknown classification rule op: {op!r}")


def _resolve(return_spec, fields):
    if isinstance(return_spec, str) and return_spec.startswith("@"):
        return fields.get(return_spec[1:], "")
    return return_spec


def _apply_transform(value, transform):
    if transform == "replace_space_dash" and isinstance(value, str):
        return value.replace(" ", "-")
    return value


def classify(kind: str, fields: dict):
    """Classify ``fields`` using the ordered rule table named ``kind``.

    Returns the first matching rule's resolved code, or ``None`` if no table
    exists for ``kind`` (e.g. config missing).
    """
    table = _RULES.get(kind)
    if not table:
        return None
    for rule in table["rules"]:
        if _matches(rule, fields):
            return _apply_transform(_resolve(rule["return"], fields), rule.get("transform"))
    return None

