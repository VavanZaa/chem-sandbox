"""Audit an external VULCAN low-T NCHO mechanism without importing it.

The VULCAN repository does not publish a license. Keep its downloaded text in
private research storage; this script reports structural evidence only and
must never turn target-specific reactions into the blind apparatus mechanism.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import re

from app.prebiotic.precursor_benchmark import stribling_miller_measurement_profile


NETWORK = (
    Path(__file__).resolve().parents[2]
    / ".tools"
    / "research"
    / "VULCAN_NCHO_photo_network_lowT.txt"
)
EXPECTED_SHA256 = "f0964ef5ed0cbe4493fcf4ebe1dd209fc8d62864dfc670ea357a62441497704e"
ROW = re.compile(r"^\s*(\d+)\s*\[\s*(.*?)\s*\]\s*(.*)$")
ATOM = re.compile(r"([A-Z][a-z]?)(\d*)")
RANGE = re.compile(r"(?<!\d)(\d{2,4})\s*-\s*(\d{2,4})(?!\d)")
NUMBER = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?$")


def composition(species: str) -> Counter[str]:
    if species == "M":
        return Counter()
    bare = species.split("_", 1)[0]
    atoms = ATOM.findall(bare)
    if not atoms or "".join(element + count for element, count in atoms) != bare:
        raise ValueError(f"Unparsed chemical formula: {species}")
    return Counter({element: sum(int(count or 1) for atom, count in atoms if atom == element)
                    for element, _ in atoms})


def parse_reaction(line: str) -> dict | None:
    match = ROW.match(line)
    if match is None:
        return None
    identifier = int(match.group(1))
    equation = match.group(2).strip()
    if "->" not in equation:
        raise ValueError(f"Reaction {identifier} lacks an arrow")
    left, right = equation.split("->", 1)
    reactants = [item.strip() for item in left.split("+")]
    products = [item.strip() for item in right.split("+")]
    left_atoms: Counter[str] = Counter()
    right_atoms: Counter[str] = Counter()
    for species in reactants:
        left_atoms.update(composition(species))
    for species in products:
        right_atoms.update(composition(species))
    tail = match.group(3).strip()
    first_three = tail.split()[:3]
    has_thermal_rate = len(first_three) == 3 and all(NUMBER.fullmatch(item) for item in first_three)
    ranges = [(int(low), int(high)) for low, high in RANGE.findall(tail)]
    return {
        "id": identifier,
        "equation": equation,
        "reactants": reactants,
        "products": products,
        "atom_balanced": left_atoms == right_atoms,
        "thermal_rate_triplet_present": has_thermal_rate,
        "explicit_temperature_ranges_K": ranges,
        "unknown_or_estimated_provenance": bool(re.search(r"\?\?|\bEstimate\b|\bnot\b", tail)),
    }


def audit_network(path: Path = NETWORK) -> dict:
    if not path.is_file():
        raise FileNotFoundError(path)
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != EXPECTED_SHA256:
        raise ValueError(f"VULCAN research copy SHA-256 differs: {actual}")
    rows = [item for line in path.read_text(encoding="utf-8").splitlines()
            if (item := parse_reaction(line)) is not None]
    temperature = float(
        stribling_miller_measurement_profile()["apparatus"]["discharge_temperature"]["value"]
    )
    thermal = [item for item in rows if item["thermal_rate_triplet_present"]]
    temperature_covered = [
        item for item in thermal
        if any(low <= temperature <= high for low, high in item["explicit_temperature_ranges_K"])
    ]
    reactant_species = {species for item in thermal for species in item["reactants"]}
    product_species = {species for item in thermal for species in item["products"]}
    return {
        "status": "research_only_not_qualified_for_blind_solver",
        "source": "https://github.com/shami-EEG/VULCAN/blob/master/thermo/NCHO_photo_network_lowT.txt",
        "sha256": actual,
        "declared_discharge_temperature_K": temperature,
        "parsed_reaction_rows": len(rows),
        "thermal_rate_rows": len(thermal),
        "nonthermal_or_special_rows": len(rows) - len(thermal),
        "atom_unbalanced_ids": [item["id"] for item in rows if not item["atom_balanced"]],
        "thermal_rows_with_explicit_range_covering_discharge_temperature": len(temperature_covered),
        "thermal_rows_without_explicit_covering_range": len(thermal) - len(temperature_covered),
        "thermal_rows_with_unknown_or_estimated_provenance": sum(
            item["unknown_or_estimated_provenance"] for item in thermal
        ),
        "feed_species_available_as_reactants": {
            name: name in reactant_species for name in ("CH4", "N2", "H2O", "H2")
        },
        "distinct_thermal_reactant_species": len(reactant_species - {"M"}),
        "distinct_thermal_product_species": len(product_species - {"M"}),
        "requires_photolysis_not_present_in_spark": any(
            item["id"] >= 627 for item in rows
        ),
        "mechanism_data_license_identified": False,
        "not_imported_into_apparatus": True,
        "not_a_validated_400K_wet_spark_mechanism": True,
    }


if __name__ == "__main__":
    print(json.dumps(audit_network(), indent=2, sort_keys=True))
