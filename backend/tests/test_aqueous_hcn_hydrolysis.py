"""Held-out HCN hydrolysis screening from independent aqueous kinetics."""

import copy
import math

import pytest

from app.prebiotic.aqueous_hcn_hydrolysis import (
    evaluate_aqueous_hcn_hydrolysis_screen,
    hcn_hydrolysis_rate_screen,
    load_hcn_hydrolysis_rate_card,
)


def _exposure(*, cyanide_mol: float = 1e-6, pH_status: str = "ideal_limited_species_pH_estimated") -> dict:
    place = {
        "molecular_hcn_mol": cyanide_mol,
        "cyanide_ion_mol": 0.0,
        "water_mol": 2.775,
        "temperature_K": 298.15,
        "pH_ideal_aqueous_model": 4.53,
        "pKa_HCN": 9.21,
        "pKw": 13.994,
        "pH_status": pH_status,
    }
    zero = place | {"molecular_hcn_mol": 0.0}
    return {"step_snapshots": [{
        "cycle": 1, "substep": 1, "duration_s": 925 * 60,
        "compartments": {
            "boiler_liquid": place,
            "condensate_trap": zero,
            "collected_condensate": zero,
        },
    }]}


def test_primary_table_rate_and_parameter_range_are_unit_explicit() -> None:
    card = load_hcn_hydrolysis_rate_card()
    assert "2608.00025" in card["source"]
    assert card["acid"]["activation_energy_J_mol_from_table_1"] == 60000
    assert card["abstract_acid_activation_energy_J_mol"] == 67100
    rate = hcn_hydrolysis_rate_screen(298.15, 4.53, 9.21, 13.994)
    assert rate["nominal_s_inv"] == pytest.approx(8.263e-10, rel=0.002)
    assert rate["parameter_low_s_inv"] < rate["nominal_s_inv"] < rate["parameter_high_s_inv"]
    assert rate["abstract_acid_energy_sensitivity_s_inv"] > rate["nominal_s_inv"]
    assert 0 < rate["molecular_hcn_fraction"] <= 1


def test_hydrolysis_screen_is_time_dependent_but_does_not_mutate_prediction() -> None:
    exposure = _exposure()
    frozen = copy.deepcopy(exposure)
    result = evaluate_aqueous_hcn_hydrolysis_screen(exposure)
    assert exposure == frozen
    assert result["status"] == "read_only_step_endpoint_hydrolysis_screen"
    rate = hcn_hydrolysis_rate_screen(298.15, 4.53, 9.21, 13.994)
    expected = 1e-6 * -math.expm1(-rate["nominal_s_inv"] * 925 * 60)
    assert result["combined_screen_mol"]["nominal_mol"] == pytest.approx(expected)
    assert result["combined_screen_mol"]["parameter_low_mol"] < expected
    assert result["combined_screen_mol"]["parameter_high_mol"] > expected
    assert result["prediction_inventory_mutated"] is False
    assert result["kinetic_products_injected"] is False
    assert result["thermal_ledger_closed_by_this_screen"] is False


def test_missing_ph_or_outside_temperature_fails_closed_for_nonzero_cyanide() -> None:
    missing = evaluate_aqueous_hcn_hydrolysis_screen(
        _exposure(pH_status="speciation_unavailable_for_inventory")
    )
    assert missing["status"] == "unavailable_for_one_or_more_cyanide_exposures"
    assert missing["unavailable_exposure_count"] == 1
    assert len(missing["unavailable_exposures"]) == 1
    too_hot = _exposure()
    too_hot["step_snapshots"][0]["compartments"]["boiler_liquid"]["temperature_K"] = 373.15
    result = evaluate_aqueous_hcn_hydrolysis_screen(too_hot)
    assert result["status"] == "unavailable_for_one_or_more_cyanide_exposures"
    assert "outside the screening window" in result["unavailable_exposures"][0]["reason"]
    assert result["combined_screen_mol"]["nominal_mol"] == 0


def test_long_unsupported_exposure_keeps_count_but_limits_details() -> None:
    exposure = _exposure(pH_status="speciation_unavailable_for_inventory")
    exposure["step_snapshots"] *= 100
    result = evaluate_aqueous_hcn_hydrolysis_screen(exposure)
    assert result["unavailable_exposure_count"] == 100
    assert len(result["unavailable_exposures"]) == 20


def test_faster_alkaline_loss_and_zero_cyanide_negative_control() -> None:
    acidic = hcn_hydrolysis_rate_screen(298.15, 4.53, 9.21, 13.994)
    alkaline = hcn_hydrolysis_rate_screen(298.15, 12.0, 9.21, 13.994)
    assert alkaline["nominal_s_inv"] > acidic["nominal_s_inv"]
    zero = evaluate_aqueous_hcn_hydrolysis_screen(_exposure(cyanide_mol=0))
    assert zero["status"] == "read_only_step_endpoint_hydrolysis_screen"
    assert zero["combined_screen_mol"]["nominal_mol"] == 0
