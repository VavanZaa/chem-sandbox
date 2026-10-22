from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from ase import Atoms

from app.experiment.schema import ExperimentDefinition


@dataclass(slots=True)
class BuiltSystem:
    atoms: Atoms
    fragment_names: list[str]
    notes: list[str]


class MoleculeBuilder:
    """Build the bounded v0.1 atomistic starting configuration.

    Amounts describe the macroscopic experiment. To keep the search bounded, the
    atomistic proxy uses one formula unit of each reactant plus one explicit water.
    This is a documented sampling assumption, not a claim of bulk-solvent realism.
    """

    @staticmethod
    def _bicarbonate() -> Atoms:
        return Atoms(
            symbols=["C", "O", "O", "O", "H"],
            positions=[
                [0.00, 0.00, 0.00],
                [1.23, 0.00, 0.00],
                [-0.62, 1.07, 0.00],
                [-0.70, -1.20, 0.00],
                [-1.66, -1.25, 0.00],
            ],
        )

    @staticmethod
    def _sodium() -> Atoms:
        return Atoms("Na", positions=[[0.0, 0.0, 0.0]])

    @staticmethod
    def _acetic_acid() -> Atoms:
        return Atoms(
            symbols=["C", "C", "O", "O", "H", "H", "H", "H"],
            positions=[
                [0.00, 0.00, 0.00],
                [1.50, 0.00, 0.00],
                [2.18, 1.02, 0.00],
                [2.16, -1.12, 0.00],
                [3.10, -1.14, 0.00],
                [-0.39, 0.98, 0.00],
                [-0.39, -0.49, 0.85],
                [-0.39, -0.49, -0.85],
            ],
        )

    @staticmethod
    def _water() -> Atoms:
        return Atoms(
            symbols=["O", "H", "H"],
            positions=[[0.0, 0.0, 0.0], [0.76, 0.59, 0.0], [-0.76, 0.59, 0.0]],
        )

    def build(self, experiment: ExperimentDefinition, seed: int | None = None) -> BuiltSystem:
        if experiment.mode == "general":
            from .general_builder import build_general
            return build_general(experiment)
        rng = np.random.default_rng(experiment.search.seed if seed is None else seed)
        fragments = [
            ("bicarbonate", self._bicarbonate(), np.array([0.0, 0.0, 0.0])),
            ("sodium", self._sodium(), np.array([0.0, 0.0, 3.2])),
            # Acid O-H points toward the exposed bicarbonate oxygen.
            ("acetic_acid", self._acetic_acid(), np.array([-5.05, 2.20, 0.15])),
            ("water", self._water(), np.array([1.4, 3.7, 0.5])),
        ]
        combined = Atoms()
        fragment_ids: list[int] = []
        fragment_names: list[str] = []
        for fragment_id, (name, fragment, offset) in enumerate(fragments):
            fragment = fragment.copy()
            fragment.positions += offset
            combined += fragment
            fragment_ids.extend([fragment_id] * len(fragment))
            fragment_names.append(name)
        combined.set_tags(fragment_ids)
        combined.positions += rng.normal(0.0, 0.015, size=combined.positions.shape)
        combined.set_cell([18.0, 18.0, 18.0])
        combined.center()
        combined.pbc = False
        return BuiltSystem(
            atoms=combined,
            fragment_names=fragment_names,
            notes=[
                "One NaHCO3, one CH3COOH, and one explicit H2O form the bounded atomistic proxy.",
                "Macroscopic amounts are retained as experiment metadata but are not converted into bulk molecule counts.",
                "Formal ionic charge is not represented by MACE-MP; sodium is modeled as a nearby counterion.",
            ],
        )
