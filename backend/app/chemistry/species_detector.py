from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import networkx as nx
from ase import Atoms

from .bond_detector import BondDetector
from .stereochemistry import identity_from_geometry


SUBSCRIPT = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


@dataclass(frozen=True, slots=True)
class DetectedSpecies:
    formula: str
    atom_indices: tuple[int, ...]
    connectivity_hash: str
    display_formula: str
    name: str
    inferred_charge: int | None
    resolved: bool
    note: str | None = None
    sort_order: int = 100
    canonical_smiles: str | None = None
    isomeric_smiles: str | None = None
    stereochemistry: str | None = None
    stereo_assignments: tuple[str, ...] = ()
    stereo_status: str = "unknown"
    aromatic: bool | None = None
    identity_basis: str = "connectivity_only"
    identity_key: str | None = None


class SpeciesDetector:
    def __init__(self, bonds: BondDetector | None = None) -> None:
        self.bonds = bonds or BondDetector()

    @staticmethod
    def _formula(symbols: list[str], unicode: bool = False) -> str:
        counts = Counter(symbols)
        order: list[str]
        if "C" in counts:
            order = ["C"] + (["H"] if "H" in counts else []) + sorted(
                symbol for symbol in counts if symbol not in {"C", "H"}
            )
        else:
            order = sorted(counts)
        output = ""
        for symbol in order:
            output += symbol
            if counts[symbol] > 1:
                count = str(counts[symbol])
                output += count.translate(SUBSCRIPT) if unicode else count
        return output

    def identify(self, atoms: Atoms) -> list[DetectedSpecies]:
        graph = self.bonds.graph(atoms)
        species: list[DetectedSpecies] = []
        for component in nx.connected_components(graph):
            indices = tuple(sorted(component))
            subgraph = graph.subgraph(indices).copy()
            graph_hash = nx.weisfeiler_lehman_graph_hash(subgraph, node_attr="element")
            symbols = [atoms[index].symbol for index in indices]
            formula = self._formula(symbols)
            identity = self._classify(formula, indices, graph, atoms)
            if atoms.info.get("chemistry_mode") == "general":
                known = atoms.info.get("species_catalog", {}).get(f"{formula}:{graph_hash}")
                instance_key = ",".join(map(str, indices))
                instance = atoms.info.get("species_instances", {}).get(instance_key)
                if instance and instance.get("connectivity_hash") != graph_hash:
                    instance = None
                charge = instance["charge"] if instance else known["charge"] if known else None
                suffix = ""
                if charge:
                    suffix = (str(abs(charge)).translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")) if abs(charge) != 1 else "") + ("⁺" if charge > 0 else "⁻")
                if instance:
                    geometry_identity = identity_from_geometry(atoms, instance)
                    stereo_status = geometry_identity["stereo_status"]
                    stereo = geometry_identity.get("stereochemistry")
                    stereo_suffix = ""
                    if stereo_status == "specified" and stereo:
                        stereo_suffix = f" [{stereo}]"
                    elif stereo_status == "unassigned":
                        stereo_suffix = " [stereo ?]"
                    elif stereo_status == "indeterminate":
                        stereo_suffix = " [stereo indeterminate]"
                    aromatic = bool(instance.get("aromatic"))
                    identity = {
                        "display_formula": self._formula(symbols, unicode=True) + suffix + stereo_suffix,
                        "name": instance["name"],
                        "inferred_charge": charge,
                        "resolved": stereo_status not in {"unassigned", "indeterminate"},
                        "sort_order": 100,
                        "note": (
                            "Stereochemistry was not specified in the input SMILES; this sampled conformer is not assigned a stereoisomer identity."
                            if stereo_status == "unassigned"
                            else "The input bond-order template is intact; stereochemistry is assigned from mapped 3D coordinates."
                            if stereo_status == "specified"
                            else "The input aromatic/bond-order template is intact."
                        ),
                        "canonical_smiles": instance.get("canonical_smiles"),
                        "isomeric_smiles": geometry_identity.get("isomeric_smiles"),
                        "stereochemistry": stereo,
                        "stereo_assignments": tuple(geometry_identity.get("stereo_assignments", ())),
                        "stereo_status": stereo_status,
                        "aromatic": aromatic,
                        "identity_basis": "mapped_input_template_and_3d_geometry",
                        "identity_key": geometry_identity["identity_key"],
                    }
                else:
                    identity = {"display_formula": self._formula(symbols, unicode=True) + suffix,
                                "name": known["name"] if known else "Unassigned sampled fragment",
                                "inferred_charge": charge,
                                "resolved": bool(known), "sort_order": 100,
                                "identity_basis": "connectivity_only",
                                "identity_key": None,
                                "note": "Legacy input identity matched by element/connectivity only; bond order and stereochemistry are not established." if known else "Geometry-derived connectivity; identity, charge and bond orders need validation."}
            species.append(
                DetectedSpecies(
                    formula=formula,
                    atom_indices=indices,
                    connectivity_hash=graph_hash,
                    **identity,
                )
            )
        return sorted(species, key=lambda item: (item.sort_order, item.formula, item.connectivity_hash, item.identity_key or ""))

    def _classify(
        self,
        formula: str,
        indices: tuple[int, ...],
        graph: nx.Graph,
        atoms: Atoms,
    ) -> dict[str, object]:
        """Assign bounded v0.1 display identities from measured connectivity.

        Charges are chemically inferred presentation metadata. They are not output
        by MACE and are kept separate from the connectivity signature and energy.
        """

        component = set(indices)
        symbols = atoms.get_chemical_symbols()

        def neighbors(index: int, symbol: str | None = None) -> list[int]:
            found = [item for item in graph.neighbors(index) if item in component]
            return found if symbol is None else [item for item in found if symbols[item] == symbol]

        def carbon_oxo_center(hydroxyl_count: int) -> bool:
            carbons = [index for index in indices if symbols[index] == "C"]
            if len(carbons) != 1:
                return False
            carbon = carbons[0]
            oxygens = neighbors(carbon, "O")
            if len(oxygens) != 3 or neighbors(carbon, "H"):
                return False
            hydroxyls = 0
            for oxygen in oxygens:
                oxygen_h = neighbors(oxygen, "H")
                if len(oxygen_h) > 1:
                    return False
                if oxygen_h and all(len(neighbors(hydrogen)) == 1 for hydrogen in oxygen_h):
                    hydroxyls += 1
            return hydroxyls == hydroxyl_count

        def acetate_backbone(hydroxyl_count: int) -> bool:
            carbons = [index for index in indices if symbols[index] == "C"]
            if len(carbons) != 2 or not graph.has_edge(*carbons):
                return False
            methyl = next((carbon for carbon in carbons if len(neighbors(carbon, "H")) == 3), None)
            if methyl is None:
                return False
            carboxyl = carbons[1] if carbons[0] == methyl else carbons[0]
            oxygens = neighbors(carboxyl, "O")
            if len(oxygens) != 2:
                return False
            return sum(1 for oxygen in oxygens if len(neighbors(oxygen, "H")) == 1) == hydroxyl_count

        resolved: dict[str, object]
        if formula == "C2H4O2" and acetate_backbone(1):
            resolved = {
                "display_formula": "CH₃COOH",
                "name": "Acetic acid",
                "inferred_charge": 0,
                "resolved": True,
                "sort_order": 10,
            }
        elif formula == "C2H3O2" and acetate_backbone(0):
            resolved = {
                "display_formula": "CH₃COO⁻",
                "name": "Acetate ion",
                "inferred_charge": -1,
                "resolved": True,
                "sort_order": 11,
            }
        elif formula == "CHO3" and carbon_oxo_center(1):
            resolved = {
                "display_formula": "HCO₃⁻",
                "name": "Bicarbonate ion",
                "inferred_charge": -1,
                "resolved": True,
                "sort_order": 20,
            }
        elif formula == "CH2O3" and carbon_oxo_center(2):
            resolved = {
                "display_formula": "H₂CO₃",
                "name": "Carbonic acid",
                "inferred_charge": 0,
                "resolved": True,
                "sort_order": 21,
            }
        elif formula == "CO2":
            carbon = next((index for index in indices if symbols[index] == "C"), -1)
            is_co2 = carbon >= 0 and len(neighbors(carbon, "O")) == 2 and len(indices) == 3
            resolved = {
                "display_formula": "CO₂" if is_co2 else "CO₂ (unresolved)",
                "name": "Carbon dioxide" if is_co2 else "Unresolved CO₂ connectivity",
                "inferred_charge": 0 if is_co2 else None,
                "resolved": is_co2,
                "sort_order": 30 if is_co2 else 90,
            }
        elif formula == "H2O":
            oxygen = next((index for index in indices if symbols[index] == "O"), -1)
            is_water = oxygen >= 0 and len(neighbors(oxygen, "H")) == 2 and len(indices) == 3
            resolved = {
                "display_formula": "H₂O" if is_water else "H₂O (unresolved)",
                "name": "Water" if is_water else "Unresolved H₂O connectivity",
                "inferred_charge": 0 if is_water else None,
                "resolved": is_water,
                "sort_order": 40 if is_water else 90,
            }
        elif formula == "Na" and len(indices) == 1:
            resolved = {
                "display_formula": "Na⁺",
                "name": "Sodium ion",
                "inferred_charge": 1,
                "resolved": True,
                "sort_order": 50,
            }
        elif formula == "H" and len(indices) == 1:
            resolved = {
                "display_formula": "H (unbound)",
                "name": "Unbound hydrogen — unresolved",
                "inferred_charge": None,
                "resolved": False,
                "note": "An isolated H atom is not accepted as a solvated aqueous proton.",
                "sort_order": 99,
            }
        else:
            display = self._formula([symbols[index] for index in indices], unicode=True)
            resolved = {
                "display_formula": f"{display} (unresolved)",
                "name": "Unresolved connectivity",
                "inferred_charge": None,
                "resolved": False,
                "note": "No bounded v0.1 species rule matches this connectivity.",
                "sort_order": 90,
            }
        resolved.setdefault("note", None)
        return resolved

    def state_signature(self, atoms: Atoms) -> str:
        return "|".join(
            f"{item.formula}:{item.connectivity_hash}:{item.identity_key or 'connectivity-only'}"
            for item in self.identify(atoms)
        )

    def state_formula(self, atoms: Atoms, unicode: bool = True) -> str:
        pieces: list[str] = []
        for item in self.identify(atoms):
            symbols = [atoms[index].symbol for index in item.atom_indices]
            pieces.append(self._formula(symbols, unicode=unicode))
        return " + ".join(pieces)

    def state_label(self, atoms: Atoms) -> str:
        counts: Counter[str] = Counter(item.display_formula for item in self.identify(atoms))
        ordered = list(dict.fromkeys(item.display_formula for item in self.identify(atoms)))
        return " + ".join(
            f"{counts[label]} {label}" if counts[label] > 1 else label for label in ordered
        )

    def state_species(self, atoms: Atoms) -> list[dict[str, object]]:
        return [
            {
                "formula": item.formula,
                "display_formula": item.display_formula,
                "name": item.name,
                "inferred_charge": item.inferred_charge,
                "resolved": item.resolved,
                "note": item.note,
                "atom_indices": list(item.atom_indices),
                "canonical_smiles": item.canonical_smiles,
                "isomeric_smiles": item.isomeric_smiles,
                "stereochemistry": item.stereochemistry,
                "stereo_assignments": list(item.stereo_assignments),
                "stereo_status": item.stereo_status,
                "aromatic": item.aromatic,
                "identity_basis": item.identity_basis,
            }
            for item in self.identify(atoms)
        ]

    def has_isolated_hydrogen(self, atoms: Atoms) -> bool:
        return any(item.formula == "H" and len(item.atom_indices) == 1 for item in self.identify(atoms))
