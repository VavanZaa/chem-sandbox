"""Audit an externally downloaded broad CHON network without running it.

This reads the unmodified KIDA V23 archive extracted under .tools/research.
It is an applicability/integrity check only: no reaction is imported into the
blind solver and no Miller-Urey product or yield is used for network selection.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

from app.prebiotic.precursor_benchmark import stribling_miller_measurement_profile


DEFAULT_DIRECTORY = Path(__file__).resolve().parents[2] / ".tools" / "research" / "V23_KIDA"
ELEMENTS = ("C", "H", "O", "N")


def compositions(directory: Path) -> dict[str, tuple[int, ...]]:
    records: dict[str, tuple[int, ...]] = {}
    for line in (directory / "composes.dat").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 7:
            continue
        name = fields[1]
        amount = tuple(int(value) for value in fields[3:7])
        if name in records:
            raise ValueError(f"Duplicate composition for {name}")
        records[name] = amount
    return records


def reaction_slots(line: str) -> tuple[list[str], list[str]]:
    if len(line) < 110:
        raise ValueError(f"Reaction line has fewer than ten fixed-width slots: {line!r}")
    slots = [line[offset * 11 + 1 : offset * 11 + 11].strip() for offset in range(10)]
    return [name for name in slots[:5] if name], [name for name in slots[5:] if name]


def audit(directory: Path) -> dict:
    species = compositions(directory)
    readme = (directory / "readme.txt").read_text(encoding="utf-8", errors="replace")
    documented = re.search(
        r"contains (\d+) chemical reactions \((\d+) reversible and (\d+) irreversible\)",
        readme,
    )
    archive = directory.parent / "V23_KIDA.7z"
    counts: Counter[str] = Counter()
    participation: Counter[str] = Counter()
    missing_species: Counter[str] = Counter()
    unbalanced: list[dict] = []
    malformed: list[dict] = []
    for path in sorted(directory.glob("reactions_*.dat")):
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if not line.strip():
                continue
            try:
                reactants, products = reaction_slots(line)
            except ValueError as error:
                malformed.append({"file": path.name, "line": line_number, "error": str(error)})
                continue
            counts[path.name] += 1
            for name in reactants + products:
                if name != "M":
                    participation[name] += 1
                    if name not in species:
                        missing_species[name] += 1
            if any(name not in species and name != "M" for name in reactants + products):
                continue
            left = tuple(sum(species[name][position] for name in reactants if name != "M")
                         for position in range(len(ELEMENTS)))
            right = tuple(sum(species[name][position] for name in products if name != "M")
                          for position in range(len(ELEMENTS)))
            if left != right:
                unbalanced.append({
                    "file": path.name, "line": line_number,
                    "reactants": reactants, "products": products,
                    "left_CHON": left, "right_CHON": right,
                })
    declared_feed = ("CH4", "N2", "H2O", "H2")
    irreversible_files = {"reactions_3.dat", "reactions_4.dat", "reactions_10.dat", "reactions_11.dat"}
    irreversible_count = sum(counts[name] for name in irreversible_files)
    reversible_count = sum(counts.values()) - irreversible_count
    readme_counts = (
        {
            "total": int(documented[1]),
            "reversible": int(documented[2]),
            "irreversible": int(documented[3]),
        }
        if documented else None
    )
    discharge_temperature_K = float(
        stribling_miller_measurement_profile()["apparatus"]["discharge_temperature"]["value"]
    )
    return {
        "network": "Veillet et al. V23 KIDA CHON, neutral chemistry only",
        "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest() if archive.exists() else None,
        "source_url": "https://kida.astrochem-tools.org/uploads/models/V23_KIDA.7z",
        "not_imported_into_blind_solver": True,
        "published_network_validity_K": [500, 2500],
        "current_declared_discharge_temperature_K": discharge_temperature_K,
        "in_published_temperature_domain": 500 <= discharge_temperature_K <= 2500,
        "neutral_species_count": len(species),
        "reaction_count": sum(counts.values()),
        "reversible_reaction_count": reversible_count,
        "irreversible_reaction_count": irreversible_count,
        "readme_reaction_counts": readme_counts,
        "archive_matches_readme_counts": (
            readme_counts == {
                "total": sum(counts.values()),
                "reversible": reversible_count,
                "irreversible": irreversible_count,
            }
            if readme_counts else None
        ),
        "reaction_counts_by_file": dict(sorted(counts.items())),
        "declared_feed_covered": all(name in species for name in declared_feed),
        "declared_feed_reaction_participation": {
            name: participation[name] for name in declared_feed
        },
        "missing_species_in_compositions": dict(missing_species),
        "malformed_lines": malformed[:20],
        "malformed_line_count": len(malformed),
        "unbalanced_reactions": unbalanced[:20],
        "unbalanced_reaction_count": len(unbalanced),
        "limitations": [
            "This archive provides kinetics but not the thermodynamic data required to reverse 1237 reactions.",
            "The published validation range starts at 500 K, above the current 400 K scenario.",
            "It is neutral chemistry; it does not replace the charge-balanced plasma branch or gas-liquid kinetics.",
            "Reaction identity and atom balance are necessary checks, not kinetic validation.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", type=Path, default=DEFAULT_DIRECTORY)
    arguments = parser.parse_args()
    print(json.dumps(audit(arguments.directory), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
