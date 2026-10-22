from __future__ import annotations

from dataclasses import dataclass

from ase import Atoms

from .bond_detector import BondChange, BondDetector
from .species_detector import SpeciesDetector


@dataclass(frozen=True, slots=True)
class ReactionEvent:
    before_signature: str
    after_signature: str
    before_formula: str
    after_formula: str
    changes: BondChange

    @property
    def connectivity_changed(self) -> bool:
        return self.before_signature != self.after_signature


class ReactionDetector:
    def __init__(self, bond_detector: BondDetector | None = None) -> None:
        self.bonds = bond_detector or BondDetector()
        self.species = SpeciesDetector(self.bonds)

    def compare(self, before: Atoms, after: Atoms) -> ReactionEvent:
        return ReactionEvent(
            before_signature=self.species.state_signature(before),
            after_signature=self.species.state_signature(after),
            before_formula=self.species.state_label(before),
            after_formula=self.species.state_label(after),
            changes=self.bonds.changes(before, after),
        )
