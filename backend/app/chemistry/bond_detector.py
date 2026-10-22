from __future__ import annotations

from dataclasses import dataclass

import networkx as nx
import numpy as np
from ase import Atoms
from ase.data import covalent_radii


@dataclass(frozen=True, slots=True)
class BondChange:
    created: tuple[tuple[int, int], ...]
    removed: tuple[tuple[int, int], ...]


class BondDetector:
    """Distance-based connectivity using a documented covalent-radius rule."""

    def __init__(self, scale: float = 1.20, minimum_distance: float = 0.55) -> None:
        self.scale = scale
        self.minimum_distance = minimum_distance

    def detect(self, atoms: Atoms) -> set[tuple[int, int]]:
        positions = np.asarray(atoms.positions)
        numbers = np.asarray(atoms.numbers)
        symbols = atoms.get_chemical_symbols()
        bonds: set[tuple[int, int]] = set()
        hydrogen_candidates: dict[int, list[tuple[float, int]]] = {}
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                distance = float(np.linalg.norm(positions[i] - positions[j]))
                cutoff = self.scale * float(covalent_radii[numbers[i]] + covalent_radii[numbers[j]])
                if self.minimum_distance <= distance <= cutoff:
                    # Sodium is an ionic counterion in the bounded v0.1 chemistry.
                    # A short Na···O contact must not merge it into a covalent species.
                    if "Na" in {symbols[i], symbols[j]} and atoms.info.get("chemistry_mode") != "general":
                        continue
                    # Hydrogen is monovalent here. Hydrogen-bond proximity can put H
                    # within the radius cutoff of two oxygens, so keep only its nearest
                    # plausible covalent partner instead of creating a bridge.
                    if symbols[i] == "H" and symbols[j] == "H":
                        if atoms.info.get("chemistry_mode") == "general":
                            hydrogen_candidates.setdefault(i, []).append((distance, j))
                            hydrogen_candidates.setdefault(j, []).append((distance, i))
                        continue
                    if symbols[i] == "H":
                        hydrogen_candidates.setdefault(i, []).append((distance, j))
                    elif symbols[j] == "H":
                        hydrogen_candidates.setdefault(j, []).append((distance, i))
                    else:
                        bonds.add((i, j))
        for hydrogen, candidates in hydrogen_candidates.items():
            _, partner = min(candidates, key=lambda item: item[0])
            if symbols[partner] == "H" and min(hydrogen_candidates[partner])[1] != hydrogen:
                continue
            bonds.add(tuple(sorted((hydrogen, partner))))
        return bonds

    def minimum_pair_distance(self, atoms: Atoms) -> float:
        if len(atoms) < 2:
            return float("inf")
        positions = np.asarray(atoms.positions)
        return min(
            float(np.linalg.norm(positions[i] - positions[j]))
            for i in range(len(atoms))
            for j in range(i + 1, len(atoms))
        )

    def graph(self, atoms: Atoms) -> nx.Graph:
        graph = nx.Graph()
        for index, symbol in enumerate(atoms.get_chemical_symbols()):
            graph.add_node(index, element=symbol)
        graph.add_edges_from(self.detect(atoms))
        return graph

    def changes(self, before: Atoms, after: Atoms) -> BondChange:
        before_bonds = self.detect(before)
        after_bonds = self.detect(after)
        return BondChange(
            created=tuple(sorted(after_bonds - before_bonds)),
            removed=tuple(sorted(before_bonds - after_bonds)),
        )
