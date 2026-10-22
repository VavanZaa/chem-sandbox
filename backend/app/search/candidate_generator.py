from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import networkx as nx
from ase import Atoms
from ase.data import covalent_radii

from app.chemistry.bond_detector import BondDetector


@dataclass(slots=True)
class CandidateProposal:
    atoms: Atoms
    method: str
    description: str
    assumptions: list[str]


class CandidateGenerator:
    """Chemically guided, product-agnostic candidate proposals.

    Rules operate on local element/connectivity motifs. They do not contain a
    product formula or a reactant-to-product lookup. The interatomic calculator
    supplies all reported energies and decides whether relaxation retains a state.
    """

    def __init__(self, bonds: BondDetector | None = None) -> None:
        self.bonds = bonds or BondDetector()

    def generate(self, atoms: Atoms, rng: np.random.Generator, limit: int) -> list[CandidateProposal]:
        """Interleave applicable operator families; never stop at the first exit.

        Randomize within the nearest-motif window, reproducibly by seed. This
        keeps a populous proton family from starving water-loss candidates.
        """
        families = [self.proton_transfers(atoms, limit=limit), self.carbon_oxo_dehydrations(atoms, limit=limit)]
        if atoms.info.get("chemistry_mode") == "general":
            families.append(self.bond_edits(atoms, limit))
        for family in families:
            rng.shuffle(family)
        rng.shuffle(families)
        proposals = []
        while any(families) and len(proposals) < limit:
            for family in families:
                if family and len(proposals) < limit:
                    proposals.append(family.pop())
        return proposals

    def bond_edits(self, atoms: Atoms, limit: int) -> list[CandidateProposal]:
        """Geometric association/dissociation seeds, not predicted reactions."""
        graph = self.bonds.graph(atoms)
        proposals = []
        components = [set(part) for part in nx.connected_components(graph)]
        for left_index, left in enumerate(components):
            for right in components[left_index + 1:]:
                pairs = sorted((float(np.linalg.norm(atoms.positions[i] - atoms.positions[j])), i, j)
                               for i in left for j in right)
                for distance, i, j in pairs[:2]:
                    if distance > 6 or distance < .55:
                        continue
                    candidate = atoms.copy()
                    vector = atoms.positions[j] - atoms.positions[i]
                    target = 1.05 * (covalent_radii[atoms.numbers[i]] + covalent_radii[atoms.numbers[j]])
                    candidate.positions[list(right)] -= vector / distance * (distance - target)
                    proposals.append(CandidateProposal(candidate, "association_scan", f"Approach fragments at atoms {i}–{j}",
                        ["Geometric encounter proposal; valence, charge and reaction barrier are not established."]))
                    if len(proposals) >= limit:
                        return proposals
        for i, j in sorted(graph.edges):
            cut = graph.copy()
            cut.remove_edge(i, j)
            if nx.has_path(cut, i, j):
                continue
            moving = list(nx.node_connected_component(cut, j))
            vector = atoms.positions[j] - atoms.positions[i]
            vector /= max(np.linalg.norm(vector), 1e-10)
            candidate = atoms.copy()
            candidate.positions[moving] += 1.7 * vector
            proposals.append(CandidateProposal(candidate, "dissociation_scan", f"Stretch bond {i}–{j}",
                ["Exploratory bond cleavage may produce radicals or ions outside a neural potential's training domain."]))
            if len(proposals) >= limit:
                break
        return proposals

    def perturb(self, atoms: Atoms, rng: np.random.Generator, strength: float = 0.10) -> Atoms:
        perturbed = atoms.copy()
        tags = np.asarray(perturbed.get_tags())
        for tag in np.unique(tags):
            indices = np.flatnonzero(tags == tag)
            center = perturbed.positions[indices].mean(axis=0)
            axis = rng.normal(size=3)
            axis /= max(float(np.linalg.norm(axis)), 1e-12)
            angle = float(rng.normal(0.0, strength))
            cross = np.array(
                [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
            )
            rotation = np.eye(3) + np.sin(angle) * cross + (1 - np.cos(angle)) * (cross @ cross)
            perturbed.positions[indices] = (perturbed.positions[indices] - center) @ rotation.T + center
            perturbed.positions[indices] += rng.normal(0.0, strength * 0.3, size=3)
        return perturbed

    def proton_transfers(self, atoms: Atoms, limit: int = 8) -> list[CandidateProposal]:
        bonds = self.bonds.detect(atoms)
        neighbors: dict[int, set[int]] = {i: set() for i in range(len(atoms))}
        for i, j in bonds:
            neighbors[i].add(j)
            neighbors[j].add(i)
        symbols = atoms.get_chemical_symbols()
        donors: list[tuple[int, int]] = []
        acceptors = [i for i, symbol in enumerate(symbols) if symbol in {"O", "N"}]
        for hydrogen, symbol in enumerate(symbols):
            if symbol != "H":
                continue
            for neighbor in neighbors[hydrogen]:
                if symbols[neighbor] in {"O", "N"}:
                    donors.append((neighbor, hydrogen))

        ranked: list[tuple[float, CandidateProposal]] = []
        for donor, hydrogen in donors:
            for acceptor in acceptors:
                if acceptor == donor or acceptor in neighbors[hydrogen]:
                    continue
                donor_acceptor = float(np.linalg.norm(atoms.positions[donor] - atoms.positions[acceptor]))
                if donor_acceptor > 5.0:
                    continue
                candidate = atoms.copy()
                direction = atoms.positions[donor] - atoms.positions[acceptor]
                direction /= max(float(np.linalg.norm(direction)), 1e-12)
                candidate.positions[hydrogen] = atoms.positions[acceptor] + direction * 0.98
                score = donor_acceptor + 0.35 * len(neighbors[acceptor])
                ranked.append(
                    (
                        score,
                        CandidateProposal(
                            atoms=candidate,
                            method="proton_transfer_scan",
                            description=f"Moved H{hydrogen} from {symbols[donor]}{donor} to {symbols[acceptor]}{acceptor}",
                            assumptions=[
                                "Heteroatom-bound hydrogens may transfer to nearby O/N acceptors.",
                                "A 0.98 Å trial acceptor-H distance seeds optimization; it is not a transition state.",
                            ],
                        ),
                    )
                )
        ranked.sort(key=lambda item: item[0])
        return [proposal for _, proposal in ranked[:limit]]

    def carbon_oxo_dehydrations(self, atoms: Atoms, limit: int = 6) -> list[CandidateProposal]:
        """Enumerate water-loss proposals from C(O)(O)(O) motifs.

        This local rule names neither carbon dioxide nor a target product. It is a
        selected bond rearrangement of the kind explicitly permitted by v0.1.
        """

        bonds = self.bonds.detect(atoms)
        neighbors: dict[int, set[int]] = {i: set() for i in range(len(atoms))}
        for i, j in bonds:
            neighbors[i].add(j)
            neighbors[j].add(i)
        symbols = atoms.get_chemical_symbols()
        proposals: list[CandidateProposal] = []
        for carbon, symbol in enumerate(symbols):
            if symbol != "C":
                continue
            oxygen_neighbors = [i for i in neighbors[carbon] if symbols[i] == "O"]
            if len(oxygen_neighbors) < 3:
                continue
            hydroxyls: list[tuple[int, int]] = []
            for oxygen in oxygen_neighbors:
                hydrogens = [i for i in neighbors[oxygen] if symbols[i] == "H"]
                hydroxyls.extend((oxygen, hydrogen) for hydrogen in hydrogens)
            for donor_oxygen, transfer_h in hydroxyls:
                for leaving_oxygen, leaving_h in hydroxyls:
                    if donor_oxygen == leaving_oxygen:
                        continue
                    candidate = atoms.copy()
                    radial = atoms.positions[leaving_oxygen] - atoms.positions[carbon]
                    radial /= max(float(np.linalg.norm(radial)), 1e-12)
                    target_oxygen = atoms.positions[carbon] + radial * 2.70
                    shift = target_oxygen - atoms.positions[leaving_oxygen]
                    candidate.positions[leaving_oxygen] += shift
                    candidate.positions[leaving_h] += shift
                    existing_vector = candidate.positions[leaving_h] - candidate.positions[leaving_oxygen]
                    existing_vector /= max(float(np.linalg.norm(existing_vector)), 1e-12)
                    trial_vector = -existing_vector + np.array([0.0, 0.35, 0.15])
                    trial_vector /= max(float(np.linalg.norm(trial_vector)), 1e-12)
                    candidate.positions[transfer_h] = candidate.positions[leaving_oxygen] + 0.97 * trial_vector
                    proposals.append(
                        CandidateProposal(
                            atoms=candidate,
                            method="selected_bond_rearrangement",
                            description=(
                                f"Scanned water loss around C{carbon}: broke C{carbon}-O{leaving_oxygen} "
                                f"and moved H{transfer_h}"
                            ),
                            assumptions=[
                                "C bonded to at least three O atoms with two O-H groups is eligible for a water-loss proposal.",
                                "The broken C-O distance is seeded at 2.70 Å; optimized connectivity is measured afterward.",
                            ],
                        )
                    )
                    if len(proposals) >= limit:
                        return proposals
        return proposals
