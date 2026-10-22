"""Element-only HCN ceilings, independent of kinetics or discharge assumptions."""
from __future__ import annotations

from collections import defaultdict
import math
from typing import Any

from .schema import PrebioticApparatusRequest


def hcn_atom_ceiling_from_atoms(
    atoms: dict[str, float], initial_cyanide_mol: float = 0.0,
    *, scheduled_additions_included: bool = False,
    atom_inventory_scope: str = "supplied_atom_inventory",
) -> dict[str, Any]:
    """Count C/N limits for HCN+CN- and the additional H limit for HCN."""
    if any(not math.isfinite(value) or value < 0 for value in atoms.values()):
        raise ValueError("Atom inventories must be finite and nonnegative")
    if not math.isfinite(initial_cyanide_mol) or initial_cyanide_mol < 0:
        raise ValueError("Initial cyanide-family inventory must be finite and nonnegative")
    carbon, nitrogen, hydrogen = (atoms.get(key, 0.0) for key in ("C", "N", "H"))
    cyanide_ceiling = min(carbon, nitrogen)
    molecular_ceiling = min(carbon, nitrogen, hydrogen)
    if initial_cyanide_mol > cyanide_ceiling + 1e-12:
        raise ValueError("Initial cyanide-family amount exceeds supplied C/N atoms")
    return {
        "status": "element_only_upper_bound",
        "available_atom_mol": {"C": carbon, "N": nitrogen, "H": hydrogen},
        "cyanide_family_total_ceiling_mol": cyanide_ceiling,
        "molecular_hcn_total_ceiling_mol": molecular_ceiling,
        "maximum_net_new_cyanide_family_mol": max(cyanide_ceiling - initial_cyanide_mol, 0.0),
        "initial_and_scheduled_cyanide_family_mol": initial_cyanide_mol,
        "limiting_atoms_for_cyanide_family": [
            key for key, amount in (("C", carbon), ("N", nitrogen))
            if abs(amount - cyanide_ceiling) <= 1e-12 * max(1.0, cyanide_ceiling)
        ],
        "scheduled_additions_included": scheduled_additions_included,
        "atom_inventory_scope": atom_inventory_scope,
        "reaction_path_or_energy_constraint_included": False,
        "predicted_yield": False,
        "interpretation": (
            "One HCN or CN- contains one C and one N; molecular HCN also needs H. "
            "This is only an element-count ceiling. It ignores competing products, "
            "incomplete conversion, sampling loss, energy barriers and transport, "
            "and is not an achievable yield."
        ),
    }


def hcn_atom_ceiling(request: PrebioticApparatusRequest) -> dict[str, Any]:
    """Return a generous bound from initial feed and every scheduled addition."""
    species = {item.id: item for item in request.species}
    atoms: dict[str, float] = defaultdict(float)
    initial_cyanide = 0.0
    for item in request.initial_inventory:
        amount = item.amount.as_unit("mol")
        if amount < 0:
            raise ValueError("Initial amounts cannot be negative")
        for element, count in species[item.species_id].composition.items():
            atoms[element] += amount * count
        if item.species_id in {"HCN", "CN-"}:
            initial_cyanide += amount
    for intervention in request.interventions:
        for item in intervention.additions:
            amount = item.amount.as_unit("mol")
            if amount < 0:
                raise ValueError("Scheduled additions cannot be negative")
            for element, count in species[item.species_id].composition.items():
                atoms[element] += amount * count
            if item.species_id in {"HCN", "CN-"}:
                initial_cyanide += amount
    return hcn_atom_ceiling_from_atoms(
        atoms, initial_cyanide, scheduled_additions_included=True,
        atom_inventory_scope="all_initial_compartments_and_scheduled_additions",
    )
