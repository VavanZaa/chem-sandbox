"""Evidence-gated molecular thermochemistry and free-energy differences."""

from .calculator import evaluate_activation_free_energy, evaluate_ensemble, evaluate_reaction_free_energy
from .harmonic import harmonic_vibrational_corrections

__all__ = [
    "evaluate_activation_free_energy",
    "evaluate_ensemble",
    "evaluate_reaction_free_energy",
    "harmonic_vibrational_corrections",
]
