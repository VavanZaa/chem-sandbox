from copy import deepcopy

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.thermochemistry.calculator import (
    evaluate_activation_free_energy,
    evaluate_ensemble,
    evaluate_reaction_free_energy,
)
from app.thermochemistry.harmonic import harmonic_vibrational_corrections
from app.thermochemistry.schema import (
    ActivationFreeEnergyRequest,
    ReactionFreeEnergyRequest,
    ThermochemicalEnsembleRequest,
    HarmonicCorrectionRequest,
)


def evidence(value: float, unit: str) -> dict:
    return {
        "value": value,
        "unit": unit,
        "source": "calculation artifact sha256:example",
        "method": "documented test method",
        "basis": "computed",
    }


def conformer_payload(
    conformer_id: str = "conf-1", *, electronic_kj_mol: float = -100.0, kind: str = "minimum"
) -> dict:
    transition = kind == "transition_state"
    return {
        "id": conformer_id,
        "degeneracy": 1,
        "identity": {
            "formula": "H2O",
            "elements": {"H": 2, "O": 1},
            "charge": 0,
            "multiplicity": 1,
            "phase": "gas",
            "standard_state": "1 atm ideal gas",
            "electronic_method": "test-level electronic method",
        },
        "stationary_point": {
            "kind": kind,
            "atom_count": 3,
            "linear": False,
            "geometry_optimized": True,
            "optimizer_converged": True,
            "max_force": evidence(0.01, "eV/angstrom"),
            "frequencies": {
                "values_cm_inv": [-500.0, 1200.0, 3600.0] if transition else [1600.0, 3700.0, 3800.0],
                "source": "frequency artifact sha256:example",
                "method": "analytic Hessian test",
                "basis": "computed",
            },
            "reaction_coordinate_verified": transition,
            "coordinate_verification_source": "IRC artifact sha256:example" if transition else None,
        },
        "electronic_energy": {
            **evidence(electronic_kj_mol, "kJ/mol"),
            "method": "test-level electronic method",
        },
        "zero_point_energy": evidence(50.0, "kJ/mol"),
        "thermal_enthalpy_correction": evidence(10.0, "kJ/mol"),
        "entropy": evidence(190.0, "J/mol/K"),
        "standard_state_correction": evidence(0.0, "kJ/mol"),
    }


def ensemble_payload(
    ensemble_id: str = "water-minima", *, electronic_kj_mol: float = -100.0, kind: str = "minimum"
) -> dict:
    return {
        "id": ensemble_id,
        "temperature": {"value": 298.15, "unit": "K"},
        "conformers": [conformer_payload(electronic_kj_mol=electronic_kj_mol, kind=kind)],
        "conformer_search_complete": True,
        "conformer_search_protocol": "exhaustive test fixture search with declared window",
    }


def test_certified_ensemble_and_population_weighting() -> None:
    payload = ensemble_payload()
    second = conformer_payload("conf-2", electronic_kj_mol=-95.0)
    second["degeneracy"] = 2
    payload["conformers"].append(second)
    result = evaluate_ensemble(ThermochemicalEnsembleRequest.model_validate(payload))
    assert result["certified"] is True
    assert result["ensemble_gibbs_free_energy_kj_mol"] is not None
    assert sum(item["population"] for item in result["conformer_populations"]) == pytest.approx(1.0)
    assert result["conformer_populations"][0]["population"] > result["conformer_populations"][1]["population"]


def test_missing_correction_or_imaginary_minimum_refuses_gibbs_claim() -> None:
    payload = ensemble_payload()
    payload["conformers"][0]["entropy"] = None
    payload["conformers"][0]["stationary_point"]["frequencies"]["values_cm_inv"][0] = -200
    result = evaluate_ensemble(ThermochemicalEnsembleRequest.model_validate(payload))
    assert result["certified"] is False
    assert result["ensemble_gibbs_free_energy_kj_mol"] is None
    problems = result["conformers"][0]["problems"]
    assert "missing entropy" in problems
    assert any("imaginary" in item for item in problems)


def test_incomplete_conformer_search_is_conditional_not_certified() -> None:
    payload = ensemble_payload()
    payload["conformer_search_complete"] = False
    result = evaluate_ensemble(ThermochemicalEnsembleRequest.model_validate(payload))
    assert result["certified"] is False
    assert result["ensemble_gibbs_free_energy_kj_mol"] is None
    assert result["conditional_submitted_ensemble_gibbs_kj_mol"] is not None


def test_activation_free_energy_requires_verified_transition_state() -> None:
    request = ActivationFreeEnergyRequest.model_validate(
        {
            "reactant": ensemble_payload(electronic_kj_mol=-100.0),
            "transition_state": ensemble_payload("water-ts", electronic_kj_mol=-50.0, kind="transition_state"),
        }
    )
    result = evaluate_activation_free_energy(request)
    assert result["certified"] is True
    assert result["activation_gibbs_free_energy_kj_mol"] == pytest.approx(50.0)

    invalid = request.model_dump(mode="json")
    invalid["transition_state"]["conformers"][0]["stationary_point"]["reaction_coordinate_verified"] = False
    invalid["transition_state"]["conformers"][0]["stationary_point"]["coordinate_verification_source"] = None
    result = evaluate_activation_free_energy(ActivationFreeEnergyRequest.model_validate(invalid))
    assert result["activation_gibbs_free_energy_kj_mol"] is None
    assert result["claim"] == "No ΔG‡ claim is made."


def test_balanced_reaction_free_energy_and_unbalanced_rejection() -> None:
    reactant = ensemble_payload("reactant", electronic_kj_mol=-100.0)
    product = ensemble_payload("product", electronic_kj_mol=-90.0)
    request = ReactionFreeEnergyRequest.model_validate(
        {
            "reactants": [{"coefficient": 1, "ensemble": reactant}],
            "products": [{"coefficient": 1, "ensemble": product}],
        }
    )
    result = evaluate_reaction_free_energy(request)
    assert result["certified"] is True
    assert result["reaction_gibbs_free_energy_kj_mol"] == pytest.approx(10.0)

    unbalanced = deepcopy(product)
    unbalanced["conformers"][0]["identity"]["elements"] = {"H": 2}
    result = evaluate_reaction_free_energy(
        ReactionFreeEnergyRequest.model_validate(
            {
                "reactants": [{"coefficient": 1, "ensemble": reactant}],
                "products": [{"coefficient": 1, "ensemble": unbalanced}],
            }
        )
    )
    assert result["certified"] is False
    assert result["reaction_gibbs_free_energy_kj_mol"] is None


def test_thermochemistry_api() -> None:
    client = TestClient(app)
    assert client.post("/api/thermochemistry/ensemble", json=ensemble_payload()).status_code == 200
    activation = {
        "reactant": ensemble_payload(electronic_kj_mol=-100.0),
        "transition_state": ensemble_payload("water-ts", electronic_kj_mol=-50.0, kind="transition_state"),
    }
    assert client.post("/api/thermochemistry/activation", json=activation).json()["certified"] is True


def test_harmonic_vibrational_corrections_are_component_only() -> None:
    request = HarmonicCorrectionRequest.model_validate(
        {
            "temperature": {"value": 298.15, "unit": "K"},
            "frequencies": {
                "values_cm_inv": [25.0, 1600.0, 3700.0],
                "source": "frequency artifact sha256:example",
                "method": "analytic Hessian test",
                "basis": "computed",
            },
            "stationary_point_kind": "minimum",
            "low_frequency_treatment": "raise_to_cutoff",
            "low_frequency_cutoff_cm_inv": 50.0,
        }
    )
    result = harmonic_vibrational_corrections(request)
    assert result["zero_point_energy_kj_mol"] > 0
    assert result["vibrational_entropy_J_mol_K"] > 0
    assert result["low_frequency_replacements"] == [{"original_cm_inv": 25.0, "used_cm_inv": 50.0}]
    assert any("vibrational" in item for item in result["limitations"])
