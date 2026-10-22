"""Product-agnostic neutral CHNO candidate enumeration from aqueous inventory.

The generator deliberately contains no product names, product SMILES, reaction
templates, or historical benchmark yields.  It combines the elemental ledgers
of species that are actually present in water, enumerates neutral closed-shell
acyclic valence graphs, and only then applies generic functional-group labels.

This is a candidate-discovery gate, not a kinetic model: a graph appearing in
the output does not establish a pathway, rate, concentration, or yield.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement, product
from typing import Any

import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions


SUPPORTED_ELEMENTS = {"C", "H", "N", "O"}
HEAVY_ELEMENTS = {"C", "N", "O"}
MAX_HEAVY_ATOMS = 6
MAX_FORMULAS = 80
MAX_STRUCTURES_PER_FORMULA = 2_000
MAX_VALENCE = {"C": 4, "N": 3, "O": 2}

POST_DISCOVERY_MOTIFS = {
    "carboxylic_acid": Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),
    "primary_amine": Chem.MolFromSmarts("[NX3;H2]"),
    "alpha_amino_carboxylic_acid": Chem.MolFromSmarts(
        "[NX3;H2,H1][CX4H1,H2][CX3](=O)[OX2H1]"
    ),
    "nitrile": Chem.MolFromSmarts("[C]#[N]"),
    "aldehyde": Chem.MolFromSmarts("[CX3H1](=O)"),
}


def _composition_key(composition: Counter[str]) -> tuple[tuple[str, int], ...]:
    return tuple(sorted((element, int(count)) for element, count in composition.items() if count))


def _formula(composition: Counter[str]) -> str:
    order = ["C", "H", "N", "O"]
    return "".join(
        f"{element}{composition[element] if composition[element] != 1 else ''}"
        for element in order
        if composition[element]
    )


def _molecule_composition(molecule: Chem.Mol) -> Counter[str]:
    explicit = Chem.AddHs(molecule)
    return Counter(atom.GetSymbol() for atom in explicit.GetAtoms())


def enumerate_neutral_tree_isomers(
    composition: dict[str, int] | Counter[str],
    *,
    maximum_structures: int = MAX_STRUCTURES_PER_FORMULA,
) -> list[dict[str, Any]]:
    """Enumerate neutral, closed-shell, acyclic CHNO graphs for one formula."""

    target = Counter({key: int(value) for key, value in composition.items() if value})
    if set(target) - SUPPORTED_ELEMENTS:
        return []
    heavy = [element for element in ("C", "N", "O") for _ in range(target[element])]
    if not 1 <= len(heavy) <= MAX_HEAVY_ATOMS:
        return []
    if target["H"] < 0:
        return []

    if len(heavy) == 1:
        trees = [nx.Graph([(0, 0)])]
    elif len(heavy) == 2:
        tree = nx.Graph()
        tree.add_nodes_from(range(2))
        tree.add_edge(0, 1)
        trees = [tree]
    else:
        trees = (
            nx.from_prufer_sequence(sequence)
            for sequence in product(range(len(heavy)), repeat=len(heavy) - 2)
        )

    unique: dict[str, dict[str, Any]] = {}
    for tree in trees:
        if len(heavy) == 1:
            edge_list: list[tuple[int, int]] = []
        else:
            edge_list = sorted(tree.edges())
        for orders in product((1, 2, 3), repeat=len(edge_list)):
            bond_valence = [0] * len(heavy)
            for (left, right), order in zip(edge_list, orders, strict=True):
                bond_valence[int(left)] += order
                bond_valence[int(right)] += order
            if any(
                bond_valence[position] > MAX_VALENCE[element]
                for position, element in enumerate(heavy)
            ):
                continue
            if sum(
                MAX_VALENCE[element] - bond_valence[position]
                for position, element in enumerate(heavy)
            ) != target["H"]:
                continue
            molecule = Chem.RWMol()
            for element in heavy:
                molecule.AddAtom(Chem.Atom(element))
            for (left, right), order in zip(edge_list, orders, strict=True):
                molecule.AddBond(
                    int(left),
                    int(right),
                    {
                        1: Chem.BondType.SINGLE,
                        2: Chem.BondType.DOUBLE,
                        3: Chem.BondType.TRIPLE,
                    }[order],
                )
            candidate = molecule.GetMol()
            try:
                Chem.SanitizeMol(candidate)
            except (Chem.KekulizeException, Chem.AtomValenceException):
                continue
            if Chem.GetFormalCharge(candidate) != 0 or Descriptors.NumRadicalElectrons(candidate):
                continue
            if _molecule_composition(candidate) != target:
                continue
            stereoisomers = tuple(
                EnumerateStereoisomers(
                    candidate,
                    options=StereoEnumerationOptions(
                        onlyUnassigned=True,
                        unique=True,
                        maxIsomers=16,
                    ),
                )
            ) or (candidate,)
            for stereoisomer in stereoisomers:
                Chem.AssignStereochemistry(stereoisomer, cleanIt=True, force=True)
                smiles = Chem.MolToSmiles(
                    stereoisomer, canonical=True, isomericSmiles=True
                )
                if "." in smiles or smiles in unique:
                    continue
                motifs = [
                    label
                    for label, pattern in POST_DISCOVERY_MOTIFS.items()
                    if pattern is not None and stereoisomer.HasSubstructMatch(pattern)
                ]
                chiral_centers = [
                    {"atom_index": int(position), "descriptor": descriptor}
                    for position, descriptor in Chem.FindMolChiralCenters(
                        stereoisomer,
                        includeUnassigned=True,
                        useLegacyImplementation=False,
                    )
                ]
                unique[smiles] = {
                    "smiles": smiles,
                    "formula": rdMolDescriptors.CalcMolFormula(stereoisomer),
                    "motifs": motifs,
                    "heavy_atoms": len(heavy),
                    "ring_count": int(rdMolDescriptors.CalcNumRings(stereoisomer)),
                    "stereochemistry": chiral_centers,
                }
                if len(unique) >= maximum_structures:
                    return sorted(unique.values(), key=lambda item: item["smiles"])
    return sorted(unique.values(), key=lambda item: item["smiles"])


def discover_aqueous_candidates(
    liquid_mol: np.ndarray,
    prepared: dict[str, Any],
    *,
    minimum_amount_mol: float = 1e-14,
) -> dict[str, Any]:
    """Enumerate structures from formulas supplied by captured liquid species."""

    available: list[dict[str, Any]] = []
    for position, amount in enumerate(liquid_mol):
        if amount <= minimum_amount_mol:
            continue
        identity = prepared["identities"][position]
        if identity.get("charge", 0) != 0:
            # Acid-base ions are carried by the material ledger but cannot be
            # treated as neutral molecular feedstocks for tree enumeration.
            continue
        composition = Counter(identity["composition"])
        if set(composition) - SUPPORTED_ELEMENTS or not set(composition) & HEAVY_ELEMENTS:
            continue
        if sum(composition[element] for element in HEAVY_ELEMENTS) > 3:
            continue
        available.append(
            {
                "id": identity["id"],
                "amount_mol": float(amount),
                "composition": composition,
            }
        )

    water_position = prepared["index"].get("H2O")
    water_available = bool(
        water_position is not None and liquid_mol[water_position] > minimum_amount_mol
    )
    formula_sources: dict[tuple[tuple[str, int], ...], set[tuple[str, ...]]] = defaultdict(set)
    for size in (2, 3):
        for group in combinations(available, size):
            composition = sum((item["composition"] for item in group), Counter())
            heavy_atoms = sum(composition[element] for element in HEAVY_ELEMENTS)
            if composition["C"] and 1 <= heavy_atoms <= MAX_HEAVY_ATOMS:
                formula_sources[_composition_key(composition)].add(
                    tuple(sorted(item["id"] for item in group))
                )
    if water_available:
        water = Counter({"H": 2, "O": 1})
        for size in (1, 2):
            for group in combinations_with_replacement(available, size):
                composition = sum((item["composition"] for item in group), water.copy())
                heavy_atoms = sum(composition[element] for element in HEAVY_ELEMENTS)
                if composition["C"] and 1 <= heavy_atoms <= MAX_HEAVY_ATOMS:
                    formula_sources[_composition_key(composition)].add(
                        tuple(sorted([*(item["id"] for item in group), "H2O"]))
                    )

    records = []
    for key in sorted(formula_sources, key=lambda item: (sum(value for name, value in item if name != "H"), item))[:MAX_FORMULAS]:
        composition = Counter(dict(key))
        structures = enumerate_neutral_tree_isomers(composition)
        if not structures:
            continue
        records.append(
            {
                "formula": _formula(composition),
                "composition": dict(composition),
                "source_combinations": [list(group) for group in sorted(formula_sources[key])],
                "candidate_count": len(structures),
                "enumeration_truncated": len(structures) >= MAX_STRUCTURES_PER_FORMULA,
                "candidates": structures,
            }
        )

    motif_counts = Counter(
        motif
        for record in records
        for candidate in record["candidates"]
        for motif in candidate["motifs"]
    )
    return {
        "status": "product_agnostic_candidates_enumerated" if records else "no_supported_candidates",
        "aqueous_source_species": [
            {"id": item["id"], "amount_mol": item["amount_mol"]} for item in available
        ],
        "formula_count": len(records),
        "candidate_structure_count": sum(record["candidate_count"] for record in records),
        "post_discovery_motif_counts": dict(motif_counts),
        "formula_records": records,
        "limitations": [
            "Only neutral, closed-shell, acyclic CHNO structures with at most six heavy atoms are enumerated; unassigned stereocenters are expanded into explicit stereoisomers.",
            "Candidate presence is not pathway feasibility, thermodynamic stability, concentration, or yield.",
            "No reaction template, known product label, benchmark yield, or target structure is used for generation.",
            "Energy ranking, transition-state searches, aqueous activities, pH and kinetics are pending.",
        ],
        "scientific_firewall": {
            "known_products_accessible_to_generator": False,
            "reaction_templates_accessible_to_generator": False,
            "historical_yields_accessible_to_generator": False,
            "motifs_applied_after_generation": True,
            "quantitative_yield_claim_allowed": False,
        },
    }
