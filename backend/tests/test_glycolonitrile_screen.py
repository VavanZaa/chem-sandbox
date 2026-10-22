"""Glycolonitrile screen is conditional and cannot change the reactor."""

from copy import deepcopy

import pytest

from app.prebiotic.glycolonitrile_screen import (
    endpoint_glycolonitrile_equilibrium,
    endpoint_glycolonitrile_association_timescale,
    evaluate_glycolonitrile_step_exposure,
    evaluate_glycolonitrile_endpoint_screen,
)


def _endpoint(**overrides):
    values = dict(
        aqueous_cyanide_mol=1e-6,
        aqueous_formaldehyde_pool_mol=1e-6,
        liquid_volume_L=1.0,
        temperature_K=298.15,
        pH=7.0,
        pKa_HCN=9.21,
    )
    values.update(overrides)
    return endpoint_glycolonitrile_equilibrium(**values)


def test_fixed_pH_equilibrium_closes_stoichiometry_and_mass_action() -> None:
    report = _endpoint()
    bound = report["conditional_glycolonitrile_mol"]
    assert 0 < bound < 1e-6
    assert bound + report["conditional_unbound_cyanide_mol"] == pytest.approx(1e-6)
    assert bound + report["conditional_unbound_formaldehyde_pool_mol"] == pytest.approx(1e-6)
    hcn_M = (
        report["hcn_fraction_at_frozen_pH"]
        * report["conditional_unbound_cyanide_mol"]
    )
    formaldehyde_M = report["conditional_unbound_formaldehyde_pool_mol"]
    assert hcn_M * formaldehyde_M / bound == pytest.approx(
        report["dissociation_K_M"], rel=1e-12
    )
    assert report["dissociation_K_M"] == pytest.approx(2.1e-6, rel=0.05)


def test_absent_precursor_and_source_temperature_guard() -> None:
    assert _endpoint(aqueous_cyanide_mol=0)["conditional_glycolonitrile_mol"] == 0
    assert _endpoint(aqueous_formaldehyde_pool_mol=0)["conditional_glycolonitrile_mol"] == 0
    with pytest.raises(ValueError, match="measured range"):
        _endpoint(temperature_K=280.0)
    with pytest.raises(ValueError, match="nonnegative"):
        _endpoint(aqueous_cyanide_mol=-1)


def test_forward_timescale_uses_only_cn_minus_and_free_formaldehyde() -> None:
    report = endpoint_glycolonitrile_association_timescale(
        aqueous_free_formaldehyde_mol=1e-10,
        liquid_volume_L=0.05,
        temperature_K=298.15,
        pH=4.7,
        pKa_HCN=9.21,
        duration_s=925 * 60,
    )
    expected_fraction = 10 ** (4.7 - 9.21) / (1 + 10 ** (4.7 - 9.21))
    expected_rate = 3.5e5 * expected_fraction * (1e-10 / 0.05)
    assert report["status"] == "source_pH_bounded_constant_endpoint_association_timescale"
    assert report["pseudo_first_order_forward_rate_per_s"] == pytest.approx(expected_rate)
    assert report["constant_endpoint_one_way_conversion_fraction_over_duration"] < 0.01
    assert report["is_reactor_yield_or_trajectory_bound"] is False
    assert endpoint_glycolonitrile_association_timescale(
        aqueous_free_formaldehyde_mol=0,
        liquid_volume_L=0.05,
        temperature_K=298.15,
        pH=4.7,
        pKa_HCN=9.21,
        duration_s=925 * 60,
    )["constant_endpoint_one_way_conversion_fraction_over_duration"] == 0
    with pytest.raises(ValueError, match="298.15 K"):
        endpoint_glycolonitrile_association_timescale(
            aqueous_free_formaldehyde_mol=1e-9,
            liquid_volume_L=0.05,
            temperature_K=310.0,
            pH=4.7,
            pKa_HCN=9.21,
            duration_s=1,
        )


def test_hydrated_and_free_formaldehyde_are_counted_once_without_mutation() -> None:
    result = {
        "final_compartments_mol": {
            "boiler_liquid": {"HCN": 1e-6, "CH2O": 1e-9, "CH4O2_AQ": 9e-9},
            "condensate_trap": {},
            "collected_condensate": {},
        },
        "aqueous_cyanide_speciation": {"compartments": {
            "boiler_liquid": {
                "status": "liquid_ideal_acid_base_equilibrium_estimated",
                "cyanide_total_mol": 1e-6,
                "liquid_volume_L": 0.05,
                "temperature_K": 298.15,
                "pH_ideal_aqueous_model": 7.0,
                "pKa_HCN": 9.21,
            },
            "condensate_trap": {"status": "no_liquid_water_inventory"},
            "collected_condensate": {"status": "no_liquid_water_inventory"},
        }},
    }
    before = deepcopy(result)
    report = evaluate_glycolonitrile_endpoint_screen(result, duration_s=925 * 60)
    assert result == before
    assert report["compartments"]["boiler_liquid"]["aqueous_formaldehyde_pool_input_mol"] == pytest.approx(1e-8)
    assert report["compartments"]["boiler_liquid"]["conditional_glycolonitrile_mol"] <= 1e-8
    assert report["compartments"]["boiler_liquid"]["association_timescale"]["within_measured_pH"] is False
    assert report["compartments"]["boiler_liquid"]["association_timescale"]["aqueous_free_formaldehyde_mol"] == 1e-9
    assert report["compartments"]["condensate_trap"]["status"] == "aqueous_speciation_unavailable"
    assert report["known_product_or_yield_used_by_simulator"] is False


def test_step_opportunity_separates_source_domain_from_extrapolation() -> None:
    def step(pH, cyanide, free_formaldehyde):
        return {
            "duration_s": 60.0,
            "compartments": {
                "boiler_liquid": {
                    "molecular_hcn_mol": cyanide,
                    "cyanide_ion_mol": 0.0,
                    "aqueous_free_formaldehyde_mol": free_formaldehyde,
                    "liquid_volume_L": 0.05,
                    "temperature_K": 298.15,
                    "pH_ideal_aqueous_model": pH,
                    "pKa_HCN": 9.21,
                    "pH_status": "ideal_limited_species_pH_estimated",
                },
                "condensate_trap": {
                    "molecular_hcn_mol": 0.0, "cyanide_ion_mol": 0.0,
                    "aqueous_free_formaldehyde_mol": 0.0,
                },
                "collected_condensate": {
                    "molecular_hcn_mol": 0.0, "cyanide_ion_mol": 0.0,
                    "aqueous_free_formaldehyde_mol": 0.0,
                },
            },
        }
    exposure = {"step_snapshots": [
        step(4.7, 1e-6, 1e-9), step(5.2, 2e-6, 1e-9),
    ]}
    before = deepcopy(exposure)
    report = evaluate_glycolonitrile_step_exposure(exposure)
    assert exposure == before
    boiler = report["by_compartment"]["boiler_liquid"]
    assert boiler["source_pH_bounded_gross_forward_opportunity_mol"] > 0
    assert boiler["out_of_pH_domain_conditional_gross_forward_opportunity_mol"] > 0
    assert boiler["outside_measured_pH_steps"] == 1
    assert report["all_cyanide_bearing_steps_within_measured_pH"] is False
    assert report["product_injected_into_solver"] is False
    exposure["step_snapshots"] = [step(4.7, 0.0, 1e-9)]
    empty = evaluate_glycolonitrile_step_exposure(exposure)
    assert empty["status"] == "no_aqueous_cyanide_to_screen"
    assert empty["all_cyanide_bearing_steps_within_measured_pH"] is None
