import numpy as np
import pytest
from ase import Atoms
from ase.io import read, write
from pydantic import ValidationError
from app.chemistry.molecule_builder import MoleculeBuilder
from app.chemistry.species_detector import SpeciesDetector
from app.experiment.schema import ExperimentDefinition
from app.search.search_runner import SearchRunner
from app.search.refinement import PathKineticsRequest, predict_path_populations, safe_structure
from app.search.consensus import build_discovery_comparison
from app.graph.reaction_graph import ReactionGraphBuilder


def definition(smiles="CCO", engine="mace_mp"):
    return ExperimentDefinition(mode="general", molecules=[{"name": "Input", "smiles": smiles, "amount_mol": .05}], compute={"engine": engine})


def test_general_builder_is_reproducible_and_names_roundtrip(tmp_path):
    first = MoleculeBuilder().build(definition())
    second = MoleculeBuilder().build(definition())
    assert len(first.atoms) == 9
    np.testing.assert_allclose(first.atoms.positions, second.atoms.positions)
    assert SpeciesDetector().identify(first.atoms)[0].name == "Input"
    path = tmp_path / "input.xyz"
    write(path, first.atoms)
    restored = read(path)
    assert restored.info["chemistry_mode"] == "general"
    assert SpeciesDetector().identify(restored)[0].name == "Input"


@pytest.mark.parametrize("smiles", ["not a molecule", "*", "C" * 40])
def test_invalid_or_oversized_molecules_rejected(smiles):
    with pytest.raises(ValidationError):
        definition(smiles)


@pytest.mark.parametrize("smiles", ["[Na+].[Cl-]", "[CH3]", "[NH4+]", "[Fe]"])
def test_ani_domain_rejected(smiles):
    with pytest.raises(ValidationError, match="ANI"):
        definition(smiles, "ani2x")


@pytest.mark.parametrize("smiles", ["[Na+]", "[CH3]"])
def test_aimnet_general_domain_rejected(smiles):
    with pytest.raises(ValidationError, match="AIMNet2"):
        definition(smiles, "aimnet2")


@pytest.mark.parametrize("smiles", ["Cl", "[NH4+]"])
def test_aimnet_rxn_domain_rejected(smiles):
    with pytest.raises(ValidationError, match="AIMNet2-RXN"):
        definition(smiles, "aimnet2_rxn")


def test_all_compatible_defers_to_explicit_domain_selection():
    experiment = definition("[Na+].[Cl-]", "all_compatible")
    assert experiment.compute.engine == "all_compatible"


def test_general_unknown_fragments_remain_visible(tmp_path, deterministic_calculator):
    atoms = Atoms("H", positions=[[0, 0, 0]])
    atoms.info["chemistry_mode"] = "general"
    runner = SearchRunner(deterministic_calculator, output_dir=tmp_path / "runs" / "test")
    runner.output_dir.mkdir(parents=True)
    result = runner._evaluate(atoms, run_index=0, name="test", steps=1, fmax=.35, method="test")
    assert result is not None
    assert not SpeciesDetector().identify(result[0].atoms)[0].resolved
    assert result[0].metadata["structure_reference"]


def test_general_search_preserves_paired_structures(tmp_path, deterministic_calculator):
    experiment = definition()
    experiment.search.runs = 1
    experiment.search.max_evaluations_per_run = 4
    result = SearchRunner(deterministic_calculator, output_dir=tmp_path / "runs" / "test").run(experiment)
    assert result["status"] == "completed"
    assert result["stoichiometry"] is None
    assert result["graph"]["edges"]
    assert all(node["comparison_key"] for node in result["graph"]["nodes"])
    assert all(edge["comparison_key"] for edge in result["graph"]["edges"])
    for edge in result["graph"]["edges"]:
        before = safe_structure(edge["source_reference"], tmp_path / "runs")
        after = safe_structure(edge["target_reference"], tmp_path / "runs")
        assert before.get_chemical_symbols() == after.get_chemical_symbols()
    with pytest.raises(ValueError):
        safe_structure("../../outside.xyz", tmp_path / "runs")


def test_unverified_or_missing_engine_barriers_never_become_zero_probability():
    graph = {"nodes": [{"id": "A"}, {"id": "B"}], "edges": [{"id": "ab", "source": "A", "target": "B"}]}
    refinement = {"graph": graph, "settings": {"engines": ["mace_mp", "ani2x"]},
                  "records": [{"engine": "mace_mp", "edge_id": "ab", "eligible_for_exploratory_kinetics": False}]}
    req = PathKineticsRequest(initial_state_id="A", target_state_id="B", acknowledge_electronic_barrier_proxy=True)
    result = predict_path_populations(refinement, req)
    assert result["status"] == "insufficient_barrier_evidence"
    assert len(result["unavailable"]) == 2
    assert not result["predictions"]
    assert result["actual_yield"] is None


def test_cuda_requested_does_not_silently_fallback(monkeypatch):
    import torch
    from app.models.mace_adapter import MACECalculatorAdapter
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    with pytest.raises(RuntimeError, match="CPU fallback was not performed"):
        MACECalculatorAdapter(device="cuda")._resolve_device()
    assert MACECalculatorAdapter(device="auto")._resolve_device() == "cpu"


def test_input_charge_display_and_neutral_hydrogen():
    ions = MoleculeBuilder().build(definition('[Na+].[Cl-]')).atoms
    assert {s.display_formula for s in SpeciesDetector().identify(ions)} == {'Na⁺', 'Cl⁻'}
    hydrogen = Atoms('H2', positions=[[0, 0, 0], [.74, 0, 0]])
    hydrogen.info['chemistry_mode'] = 'general'
    assert len(SpeciesDetector().identify(hydrogen)) == 1


def test_verified_paths_only_enable_explicit_proxy_not_bulk_yield():
    from app.search.refinement import RefinementRequest
    with pytest.raises(ValidationError):
        RefinementRequest(engines=['mace_mp', 'mace_mp'])
    refinement = {'graph': {'nodes': [{'id': 'A'}, {'id': 'B'}], 'edges': [{'id': 'ab', 'source': 'A', 'target': 'B'}]},
                  'settings': {'engines': ['mace_mp', 'ani2x']},
                  'records': [{'engine': 'mace_mp', 'edge_id': 'ab', 'eligible_for_exploratory_kinetics': True, 'electronic_barrier_ev': .8}]}
    req = PathKineticsRequest(initial_state_id='A', target_state_id='B', acknowledge_electronic_barrier_proxy=True)
    result = predict_path_populations(refinement, req)
    assert len(result['predictions']) == 1
    assert result['predictions'][0]['engine'] == 'mace_mp'
    assert 0 < result['predictions'][0]['target_state_population'] < 1
    assert result['unavailable'][0]['engine'] == 'ani2x'
    assert result['calibrated_product_probabilities'] is None


def test_discovery_agreement_is_not_confidence():
    def one(engine, states, paths):
        nodes = [{"id": key, "comparison_key": key, "label": key, "formula": key} for key in states]
        edges = [{"id": f"{a}_{b}", "source": a, "target": b, "comparison_key": f"{a}=>{b}",
                  "delta_energy_ev": delta} for a, b, delta in paths]
        return {"status": "completed", "model": {"engine": engine}, "graph": {"nodes": nodes, "edges": edges},
                "runtime_seconds": 1, "statistics": {"failed_optimizations": 0}}
    comparison = build_discovery_comparison([
        ("mace_mp", one("MACE", ["A", "B"], [("A", "B", 1.0)]), "one"),
        ("ani2x", one("ANI", ["A", "C"], [],), "two"),
    ], [])
    shared = next(row for row in comparison["states"] if row["comparison_key"] == "A")
    assert shared["agreement_fraction"] == 1
    assert comparison["calibrated_confidence"] is None
    assert "probability" in comparison["warnings"][0]


@pytest.mark.parametrize(
    ("smiles", "name"),
    [
        ("c1ccccc1", "Benzene"),
        ("c1ccc2ccccc2c1", "Naphthalene"),
        ("n1ccccc1", "Pyridine"),
        ("COc1ccccc1", "Anisole"),
    ],
)
def test_aromatic_input_identity_is_retained(smiles, name):
    experiment = definition(smiles)
    experiment.molecules[0].name = name
    detected = SpeciesDetector().identify(MoleculeBuilder().build(experiment).atoms)
    assert len(detected) == 1
    assert detected[0].name == name
    assert detected[0].aromatic is True
    assert detected[0].canonical_smiles
    assert detected[0].identity_basis == "mapped_input_template_and_3d_geometry"
    assert detected[0].resolved


def test_enantiomers_have_distinct_identity_and_graph_nodes():
    left = MoleculeBuilder().build(definition("C[C@H](O)C(=O)O")).atoms
    right = MoleculeBuilder().build(definition("C[C@@H](O)C(=O)O")).atoms
    detector = SpeciesDetector()
    left_species = detector.identify(left)[0]
    right_species = detector.identify(right)[0]
    assert {left_species.stereochemistry, right_species.stereochemistry} == {"R", "S"}
    assert left_species.isomeric_smiles != right_species.isomeric_smiles
    assert detector.state_signature(left) != detector.state_signature(right)
    graph = ReactionGraphBuilder()
    left_id = graph.add_state(left, 0.0, 0, converged=True, max_force=0.0, structure_reference=None)
    right_id = graph.add_state(right, 0.0, 0, converged=True, max_force=0.0, structure_reference=None)
    assert left_id != right_id


def test_reflecting_3d_coordinates_flips_enantiomer_identity():
    original = MoleculeBuilder().build(definition("C[C@H](O)C(=O)O")).atoms
    reflected = original.copy()
    reflected.positions[:, 0] *= -1
    detector = SpeciesDetector()
    before = detector.identify(original)[0]
    after = detector.identify(reflected)[0]
    assert {before.stereochemistry, after.stereochemistry} == {"R", "S"}
    assert detector.state_signature(original) != detector.state_signature(reflected)


def test_unspecified_stereochemistry_stays_unresolved():
    atoms = MoleculeBuilder().build(definition("CC(O)C(=O)O")).atoms
    species = SpeciesDetector().identify(atoms)[0]
    assert species.stereo_status == "unassigned"
    assert species.stereochemistry == "unassigned"
    assert not species.resolved
    assert "stereo ?" in species.display_formula


def test_e_and_z_alkenes_have_distinct_identity():
    trans = MoleculeBuilder().build(definition("C/C=C/C")).atoms
    cis = MoleculeBuilder().build(definition("C/C=C\\C")).atoms
    detector = SpeciesDetector()
    trans_species = detector.identify(trans)[0]
    cis_species = detector.identify(cis)[0]
    assert {trans_species.stereochemistry, cis_species.stereochemistry} == {"E", "Z"}
    assert detector.state_signature(trans) != detector.state_signature(cis)


@pytest.mark.parametrize(
    ("first_smiles", "mirror_smiles"),
    [
        ("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"),
        ("CC1=CCC[C@H](C1)C(=C)C", "CC1=CCC[C@@H](C1)C(=C)C"),
        (
            "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
            "OC[C@@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O",
        ),
    ],
)
def test_chiral_amino_acid_terpene_and_sugar_pairs_are_distinct(first_smiles, mirror_smiles):
    detector = SpeciesDetector()
    first = MoleculeBuilder().build(definition(first_smiles)).atoms
    mirror = MoleculeBuilder().build(definition(mirror_smiles)).atoms
    first_species = detector.identify(first)[0]
    mirror_species = detector.identify(mirror)[0]
    assert first_species.stereo_status == mirror_species.stereo_status == "specified"
    assert first_species.stereochemistry != mirror_species.stereochemistry
    assert detector.state_signature(first) != detector.state_signature(mirror)


def test_stereochemical_identity_survives_structure_roundtrip(tmp_path):
    atoms = MoleculeBuilder().build(definition("N[C@@H](C)C(=O)O")).atoms
    path = tmp_path / "alanine.xyz"
    write(path, atoms)
    restored = read(path)
    detector = SpeciesDetector()
    assert detector.state_signature(restored) == detector.state_signature(atoms)
    assert detector.identify(restored)[0].stereochemistry in {"R", "S"}
