import numpy as np
from ase import Atoms

from app.chemistry.bond_detector import BondDetector
from app.chemistry.molecule_builder import MoleculeBuilder
from app.chemistry.species_detector import SpeciesDetector
from app.experiment.schema import default_experiment
from app.search.candidate_generator import CandidateGenerator


def test_bond_detection_and_clash_floor() -> None:
    detector = BondDetector()
    water = Atoms("OH2", positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]])
    assert detector.detect(water) == {(0, 1), (0, 2)}
    assert detector.minimum_pair_distance(water) > detector.minimum_distance


def test_bond_detection_keeps_sodium_ionic_and_hydrogen_monovalent() -> None:
    detector = BondDetector()
    contact = Atoms(
        ["Na", "O", "H", "O"],
        positions=[[0, 0, 0], [2.2, 0, 0], [3.16, 0, 0], [4.18, 0, 0]],
    )
    bonds = detector.detect(contact)
    assert (0, 1) not in bonds
    assert sum(1 for bond in bonds if 2 in bond) == 1
    assert (1, 2) in bonds


def test_species_signature_is_translation_invariant() -> None:
    detector = SpeciesDetector()
    water = Atoms("OH2", positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]])
    moved = water.copy()
    moved.positions += np.array([10.0, -2.0, 4.0])
    assert detector.state_signature(water) == detector.state_signature(moved)
    assert detector.state_formula(water) == "H₂O"


def test_guided_search_can_propose_second_level_rearrangement_without_product_formula() -> None:
    atoms = MoleculeBuilder().build(default_experiment()).atoms
    generator = CandidateGenerator()
    proton_candidates = generator.proton_transfers(atoms, limit=24)
    assert proton_candidates
    second_level = [
        candidate
        for proton in proton_candidates
        for candidate in generator.carbon_oxo_dehydrations(proton.atoms, limit=6)
    ]
    assert second_level
    assert all("CO2" not in candidate.description for candidate in second_level)


def test_v01_species_labels_use_condensed_formulas_and_inferred_ions() -> None:
    atoms = MoleculeBuilder().build(default_experiment()).atoms
    detector = SpeciesDetector()
    assert detector.state_label(atoms) == "CH₃COOH + HCO₃⁻ + H₂O + Na⁺"
    identified = {item.name: item for item in detector.identify(atoms)}
    assert identified["Acetic acid"].formula == "C2H4O2"
    assert identified["Acetic acid"].display_formula == "CH₃COOH"
    assert identified["Bicarbonate ion"].inferred_charge == -1
    assert identified["Sodium ion"].inferred_charge == 1


def test_isolated_hydrogen_is_marked_unresolved() -> None:
    detector = SpeciesDetector()
    hydrogen = Atoms("H", positions=[[0, 0, 0]])
    assert detector.has_isolated_hydrogen(hydrogen)
    species = detector.identify(hydrogen)[0]
    assert species.resolved is False
    assert species.inferred_charge is None
