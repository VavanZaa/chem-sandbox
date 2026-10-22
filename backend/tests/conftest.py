from __future__ import annotations

import numpy as np
import pytest
from ase import Atoms

from app.models.base_calculator import CalculationResult, InteratomicCalculator


class DeterministicCalculator(InteratomicCalculator):
    @property
    def metadata(self) -> dict[str, str]:
        return {"engine": "test", "model": "quadratic", "device": "cpu", "dtype": "float64"}

    def energy(self, atoms: Atoms) -> float:
        return float(np.square(atoms.positions - atoms.positions.mean(axis=0)).sum() * 0.001)

    def forces(self, atoms: Atoms) -> np.ndarray:
        return np.zeros_like(atoms.positions)

    def optimize(self, atoms: Atoms, *, fmax: float = 0.35, steps: int = 12) -> CalculationResult:
        copied = atoms.copy()
        return CalculationResult(copied, self.energy(copied), self.forces(copied), True, 1, 0.0, self.metadata)


@pytest.fixture()
def deterministic_calculator() -> DeterministicCalculator:
    return DeterministicCalculator()

