"""Conservative in-memory normalization of LXCat BOLSIG+ process records.

The contributor file is never rewritten. Each change is reported to the
caller so a private, hash-pinned complete set can be reviewed before use.
"""

from __future__ import annotations

import re
from typing import Any


_SIMPLE_FORMULA = re.compile(r"(?:[A-Z][a-z]?\d*)+(?:\^?[+-])?")


def split_lxcat_products(label: str) -> list[str]:
    """Split actual product separators, preserving ionic and state-label pluses."""

    products: list[str] = []
    start = 0
    depth = 0
    for position, character in enumerate(label):
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
            if depth < 0:
                return [label]
        elif character == "+" and depth == 0:
            remaining = label[position + 1:].lstrip()
            if (
                not remaining
                or (position > 0 and label[position - 1] == "^")
                or remaining.startswith("+")
            ):
                continue
            products.append(label[start:position].strip())
            start = position + 1
    if depth != 0:
        return [label]
    products.append(label[start:].strip())
    return products


def normalize_lxcat_bolsig_text(
    text: str, target: str
) -> tuple[str, list[dict[str, Any]]]:
    """Remove a duplicate arrow contained only in a trailing product note.

    BOLOS expects exactly one arrow on the process target/product line. Some
    contributor exports repeat the reaction inside a parenthetical note. This
    function leaves all numeric tables and source comments byte-for-byte
    unchanged and returns an auditable in-memory edit log.
    """

    lines = text.splitlines(keepends=True)
    changes: list[dict[str, Any]] = []
    prefix = f"{target} -> "
    for index, original in enumerate(lines):
        line = original.rstrip("\r\n")
        if not line.startswith(prefix) or line.count("->") <= 1:
            continue
        note_start = line.find(" (", len(prefix))
        if (
            note_start < 0
            or not line.endswith(")")
            or "->" not in line[note_start:]
            or line[:note_start].count("->") != 1
        ):
            raise ValueError(
                f"Ambiguous {target} collision label on line {index + 1}; "
                "cannot safely remove a duplicate reaction annotation"
            )
        cleaned = line[:note_start]
        lines[index] = cleaned + original[len(line):]
        changes.append(
            {
                "line_number": index + 1,
                "kind": "trailing_reaction_annotation_removed",
                "before": line,
                "after": cleaned,
            }
        )
    return "".join(lines), changes


def _charge_token(token: str) -> str:
    if token.endswith(("^+", "^-")):
        return token
    if token.endswith(("+", "-")):
        return f"{token[:-1]}^{token[-1]}"
    return token


def normalize_lxcat_processes(
    processes: list[dict[str, Any]], target: str
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Normalize explicit space-delimited products and identical duplicate rows.

    An unequal pair of cross sections at the same energy is ambiguous and is
    rejected rather than averaged. A normalized label describes the source's
    *declared* products; it does not validate a branching ratio or yield.
    """

    normalized: list[dict[str, Any]] = []
    changes: list[dict[str, Any]] = []
    for index, original in enumerate(processes):
        process = dict(original)
        if process.get("target") != target:
            normalized.append(process)
            continue

        raw_label = process.get("product")
        if raw_label is not None:
            raw = str(raw_label)
            # A standalone plus is a product separator; a suffix plus is a
            # charge. Do not split internal excited-state descriptions such as
            # "H2 (v=0 - v=1)" into invented products.
            spaced_products = re.split(r"\s+\+\s+", raw)
            tokens = spaced_products if len(spaced_products) > 1 else raw.split()
            label = (
                "+".join(_charge_token(token) for token in tokens)
                if all(_SIMPLE_FORMULA.fullmatch(token) for token in tokens)
                else raw
            )
            if label != raw_label:
                process["product"] = label
                changes.append(
                    {
                        "process_index": index,
                        "kind": "product_label",
                        "before": raw_label,
                        "after": label,
                    }
                )

        data = process.get("data", [])
        unique_rows: list[list[float]] = []
        duplicate_count = 0
        for row in data:
            energy, cross_section = (float(value) for value in row)
            if unique_rows and energy == unique_rows[-1][0]:
                if cross_section != unique_rows[-1][1]:
                    raise ValueError(
                        f"{target} process {index} has different cross sections "
                        f"at the same energy {energy} eV"
                    )
                duplicate_count += 1
                continue
            unique_rows.append([energy, cross_section])
        if duplicate_count:
            process["data"] = unique_rows
            changes.append(
                {
                    "process_index": index,
                    "kind": "identical_duplicate_energy_rows_removed",
                    "count": duplicate_count,
                }
            )
        normalized.append(process)
    return normalized, changes
