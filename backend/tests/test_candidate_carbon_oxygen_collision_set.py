"""Private CO2/O2 complete-set intake; no HCN-yield calibration occurs here."""

import hashlib
from pathlib import Path

import pytest

from app.prebiotic.boltzmann import (
    IST_CO2_O2_DATASET_ID,
    IST_CO2_O2_SHA256,
    _load_ist_co2_o2_processes,
    collision_dataset_status,
    solve_self_consistent_eedf,
)
from app.prebiotic.schema import BoltzmannEEDFRequest


PRIVATE_EXPORT = (
    Path(__file__).resolve().parents[1]
    / "data" / "collision" / "private" / "IST_Lisbon_CO2_O2_2026-09-26"
    / "Cross section.txt"
)


def _request(energy_cells: int = 400) -> BoltzmannEEDFRequest:
    return BoltzmannEEDFRequest(
        collision_set=IST_CO2_O2_DATASET_ID,
        mole_fractions={
            "H2": 0.01, "CO2": 0.05, "O2": 0.001,
            "H2O": 0.02, "N2": 0.919,
        },
        gas_temperature_K=298.15,
        reduced_electric_field_Td=200,
        energy_cells=energy_cells,
        acknowledge_two_term_limit=True,
    )


def test_private_complete_sets_are_hash_pinned_and_not_mixed_with_partial_o2() -> None:
    assert hashlib.sha256(PRIVATE_EXPORT.read_bytes()).hexdigest() == IST_CO2_O2_SHA256
    status = collision_dataset_status(IST_CO2_O2_DATASET_ID)
    assert status["available"] is True
    coverage = status["coverage"]
    assert coverage["all_required_targets_complete"] is True
    assert coverage["process_count"] == 108
    assert coverage["targets"]["O2"]["process_count"] == 14
    assert coverage["targets"]["CO2"]["process_count"] == 17
    assert coverage["targets"]["O2"]["explicit_neutral_fragment_labels"] == [
        "O(3P)+O(3P)", "O(3P)+O(1D)",
    ]
    assert coverage["targets"]["O2"]["resolved_ground_state_fragment_labels"] == [
        "O(3P)+O(3P)",
    ]
    assert coverage["targets"]["CO2"]["explicit_neutral_fragment_labels"] == [
        "CO(X) + O(1D)", "CO(a3P) + O(3P)",
    ]
    assert coverage["targets"]["CO2"]["resolved_ground_state_fragment_labels"] == []
    assert "no LXCat tables redistributed" in status["redistribution"]
    assert "not automatically grounded" in status["scientific_scope"].lower()
    co2_attachment = next(
        item for item in _load_ist_co2_o2_processes()
        if item["target"] == "CO2" and item["kind"] == "ATTACHMENT"
    )
    assert co2_attachment["threshold"] == 3.9


def test_co2_o2_eedf_solves_with_grid_sensitivity_but_is_not_a_yield() -> None:
    coarse = solve_self_consistent_eedf(
        _request(400), include_same_set_cantera_reference=False
    )
    fine = solve_self_consistent_eedf(
        _request(800), include_same_set_cantera_reference=False
    )
    assert coarse["collision_set"]["selected_process_count"] == 89
    assert len(coarse["rate_coefficients"]) == 84
    assert coarse["eedf"]["mean_energy_eV"] == pytest.approx(
        fine["eedf"]["mean_energy_eV"], rel=0.01
    )
    def co2_fragment_rate(result: dict) -> float:
        return next(
            item["rate_coefficient_m3_s"]
            for item in result["rate_coefficients"]
            if item["target"] == "CO2"
            and item["product_label"] == "CO(X) + O(1D)"
        )
    assert co2_fragment_rate(coarse) == pytest.approx(co2_fragment_rate(fine), rel=0.01)
    assert coarse["scientific_firewall"]["quantitative_yield_claim_allowed"] is False
    co2_attachment = next(
        item for item in coarse["rate_coefficients"]
        if item["target"] == "CO2" and item["kind"] == "attachment"
    )
    assert co2_attachment["threshold_eV"] == 3.9


def test_co2_o2_same_curves_cross_solver_is_numerical_only() -> None:
    result = solve_self_consistent_eedf(_request())
    reference = result["same_collision_set_reference"]
    assert reference["available_for_this_mixture"] is True
    assert reference["private_collision_data_persisted"] is False
    assert reference["rate_comparison"]["material_within_factor_ten_fraction"] == 1.0
    assert "not an external swarm experiment" in reference["comparison_scope"]
    assert any("solver-sensitive" in warning for warning in result["numerical_warnings"])


@pytest.mark.parametrize("field_Td", [1, 9.99, 200.01])
def test_co2_o2_set_rejects_unreviewed_field_overlap(field_Td: float) -> None:
    with pytest.raises(ValueError, match="jointly documented only for 10-200 Td"):
        BoltzmannEEDFRequest(
            **(_request().model_dump() | {"reduced_electric_field_Td": field_Td})
        )


def test_old_collision_set_still_rejects_co2_feed() -> None:
    with pytest.raises(ValueError, match="unsupported targets"):
        BoltzmannEEDFRequest(
            **(_request().model_dump() | {"collision_set": "lxcat_morgan_2026-09-22"})
        )
