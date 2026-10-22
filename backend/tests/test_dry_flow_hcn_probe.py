"""Dry-flow gas-packet checks; these are not experimental yield validations."""

import pytest
import cantera as ct

from app.prebiotic.heavy_gas import (
    gri_ideal_gas_enthalpy_inventory,
    gri_ideal_gas_internal_energy_inventory,
)
from scripts.run_dry_flow_hcn_probe import _gri_inventory, run_dry_flow_probe


def test_dry_flow_probe_fails_closed_without_resolved_n2_product_states() -> None:
    result = run_dry_flow_probe(substeps=4)
    assert 0 <= result["outlet_hcn_mol"] < 1e-12
    assert 0 <= result["outlet_hcn_per_feed_carbon"] < 1
    assert result["material_ledger"]["closed"] is True
    assert result["energy_ledger"]["closed_upper_bound"] is True
    assert result["energy_ledger"]["complete_thermal_energy_balance_closed"] is False
    assert result["energy_ledger"]["gas_internal_energy_inventory_coverage_complete"] is True
    assert result["energy_ledger"]["ideal_gas_internal_energy_change_J"] == pytest.approx(
        result["energy_ledger"]["final_ideal_gas_internal_energy_J"]
        - result["energy_ledger"]["initial_ideal_gas_internal_energy_J"]
    )
    assert result["energy_ledger"]["conditional_unmeasured_non_electrical_energy_into_gas_J"] == pytest.approx(
        result["energy_ledger"]["ideal_gas_internal_energy_change_J"]
        - result["energy_ledger"]["declared_discharge_J"]
    )
    energy = result["energy_ledger"]
    assert energy["steady_flow_enthalpy_coverage_complete"] is True
    assert energy["stream_enthalpy_change_J_per_packet"] == pytest.approx(
        energy["outlet_stream_enthalpy_J_per_packet"]
        - energy["inlet_stream_enthalpy_J_per_packet"]
    )
    assert energy["stream_enthalpy_change_rate_W"] == pytest.approx(
        energy["stream_enthalpy_change_J_per_packet"]
        / result["conditions"]["simulated_hot_packet_residence_s"]
    )
    assert energy["conditional_electrical_power_minus_stream_enthalpy_rate_W"] == pytest.approx(
        result["conditions"]["average_discharge_power_W"]
        - energy["stream_enthalpy_change_rate_W"]
    )
    assert all(
        abs(residual) <= result["material_ledger"]["tolerance_mol_by_conserved_quantity"][key]
        for key, residual in result["material_ledger"]["residual_mol_by_conserved_quantity"].items()
    )
    assert result["plasma"]["density_was_power_capped"] is True
    assert result["source"]["measured_hcn_used_as_solver_input"] is False
    assert result["scientific_firewall"]["published_hcn_mechanism_used"] is False
    assert result["scientific_firewall"]["experimental_magnitude_validated"] is False
    assert result["scientific_firewall"]["impulse_substep_result_is_timestep_converged"] is False
    flow_basis = result["flow_temperature_basis_diagnostic"]
    assert flow_basis["nominal_volume_flow_implied_by_reported_residence_sccm"] == pytest.approx(200.4, rel=1e-3)
    assert flow_basis["ideal_gas_residence_at_uniform_scenario_temperature_s"] < result["conditions"]["nominal_residence_s"]
    assert flow_basis["internally_consistent_hot_packet_flow_and_power_basis"] is True
    assert result["conditions"]["simulated_hot_packet_residence_s"] == pytest.approx(
        result["inlet_total_mol"] / flow_basis["inferred_feed_molar_flow_mol_s"]
    )
    assert result["conditions"]["power_times_residence_J_per_gas_packet"] == pytest.approx(
        result["conditions"]["average_discharge_power_W"]
        * result["conditions"]["simulated_hot_packet_residence_s"]
    )
    assert result["conditions"]["nominal_ac_frequency_Hz"] == 22_000
    assert result["conditions"]["nominal_ac_cycles_per_hot_packet"] == pytest.approx(
        result["conditions"]["simulated_hot_packet_residence_s"] * 22_000
    )
    assert result["conditions"]["nominal_ac_cycles_per_chemistry_substep"] == pytest.approx(
        result["conditions"]["nominal_ac_cycles_per_hot_packet"]
        / result["conditions"]["discharge_chemistry_substeps"]
    )
    assert result["conditions"]["individual_microdischarge_pulses_resolved"] is False
    assert flow_basis["comparison_basis_reconciled"] is False
    assert "plug-flow packet" in result["conditions"]["reactor_idealization"]


def test_inventory_limited_source_has_stable_methane_products_but_trace_hcn() -> None:
    coarse = run_dry_flow_probe(
        source_forcing="inventory_limited", integration_max_step_fraction=0.05
    )
    refined = run_dry_flow_probe(
        source_forcing="inventory_limited", integration_max_step_fraction=0.0125
    )
    assert coarse["outlet_hcn_mol"] < 1e-12
    assert coarse["outlet_packet_diagnostics"]["methane_conversion_fraction"] == pytest.approx(
        refined["outlet_packet_diagnostics"]["methane_conversion_fraction"], rel=1e-5
    )
    assert coarse["outlet_packet_diagnostics"]["species_mol"]["C2H6"] == pytest.approx(
        refined["outlet_packet_diagnostics"]["species_mol"]["C2H6"], rel=1e-5
    )
    assert coarse["material_ledger"]["closed"] is True
    assert coarse["energy_ledger"]["closed_upper_bound"] is True
    assert coarse["energy_ledger"]["resolved_fragment_J"] <= coarse["energy_ledger"]["allocated_collision_J"]
    assert coarse["heavy_chemistry"]["solver"]["method"].endswith("inventory-limited plasma channels")
    assert coarse["scientific_firewall"]["experimental_magnitude_validated"] is False


def test_evolving_eedf_source_converges_without_opening_ledgers() -> None:
    coarse = run_dry_flow_probe(source_forcing="evolving_eedf", substeps=32)
    refined = run_dry_flow_probe(source_forcing="evolving_eedf", substeps=64)
    assert coarse["outlet_hcn_mol"] < 1e-12
    assert refined["outlet_hcn_mol"] < 1e-12
    assert coarse["outlet_packet_diagnostics"]["methane_conversion_fraction"] == pytest.approx(
        refined["outlet_packet_diagnostics"]["methane_conversion_fraction"], rel=2e-4
    )
    assert coarse["outlet_packet_diagnostics"]["species_mol"]["C2H6"] == pytest.approx(
        refined["outlet_packet_diagnostics"]["species_mol"]["C2H6"], rel=2e-4
    )
    assert refined["plasma"]["snapshot_count"] == 64
    covered_range = refined["plasma"]["covered_gas_fraction_range"]
    assert covered_range is not None and covered_range[0] < covered_range[1]
    assert refined["material_ledger"]["closed"] is True
    assert refined["energy_ledger"]["closed_upper_bound"] is True
    assert refined["source"]["measured_hcn_used_as_solver_input"] is False
    assert refined["scientific_firewall"]["experimental_magnitude_validated"] is False


def test_dry_outlet_diagnostics_preserve_all_species_and_carbon() -> None:
    result = run_dry_flow_probe(source_forcing="evolving_eedf", substeps=16)
    outlet = result["outlet_packet_diagnostics"]
    excitation = result["plasma"]["n2_excitation_opportunity_mol_by_product_label"]
    assert result["plasma"]["n2_excitation_opportunity_is_propagated"] is False
    assert excitation["N2(A3,v0-4)"] > 0
    assert excitation["N2(A3,v5-9)"] > 0
    assert excitation["N2(A3,v10-)"] > 0
    assert excitation["N2(a'1)"] > 0
    assert result["heavy_chemistry"]["methane_neutral_pathways"] is None
    assert outlet["measured_product_values_used_as_inputs"] is False
    assert sum(outlet["species_mole_fraction"].values()) == pytest.approx(1)
    assert outlet["outlet_pressure_Pa_from_fixed_volume_inventory"] > 0
    assert outlet["pressure_drift_fraction_from_inlet"] == pytest.approx(
        outlet["outlet_pressure_Pa_from_fixed_volume_inventory"] / 101_325 - 1
    )
    assert sum(outlet["carbon_yield_fraction_by_species"].values()) == pytest.approx(1, abs=1e-8)
    assert outlet["species_mol"]["HCN"] == pytest.approx(result["outlet_hcn_mol"])
    assert outlet["carbon_yield_fraction_by_species"]["CH4"] == pytest.approx(
        1 - outlet["methane_conversion_fraction"]
    )
    product_yields = outlet["product_carbon_yield_fraction_by_species"]
    assert outlet["methane_conversion_fraction"] == pytest.approx(
        outlet["methane_electron_fragment_source_fraction_of_feed"]
        + outlet["methane_net_neutral_chemistry_consumption_fraction_of_feed"],
        abs=1e-10,
    )
    assert outlet["methane_electron_fragment_source_consumption_mol"] == pytest.approx(
        result["heavy_chemistry"]["source_consumption_mol_by_target"]["CH4"]
    )
    product_selectivities = outlet["product_carbon_selectivity_fraction_by_species"]
    assert "CH4" not in product_yields
    assert sum(product_yields.values()) == pytest.approx(outlet["methane_conversion_fraction"], abs=1e-8)
    assert sum(product_selectivities.values()) == pytest.approx(1, abs=1e-8)
    hydrogen_selectivities = outlet["product_hydrogen_selectivity_fraction_by_species"]
    assert sum(hydrogen_selectivities.values()) == pytest.approx(1, abs=1e-8)
    assert hydrogen_selectivities["H2"] == pytest.approx(
        2 * outlet["species_mol"]["H2"]
        / (4 * result["inlet_total_mol"] * 0.01 * outlet["methane_conversion_fraction"])
    )
    assert product_yields["HCN"] == pytest.approx(result["outlet_hcn_per_feed_carbon"])
    assert product_yields["C2H6"] == pytest.approx(
        2 * outlet["species_mol"]["C2H6"] / (result["inlet_total_mol"] * 0.01)
    )
    assert product_selectivities["HCN"] < product_selectivities["C2H6"]
    assert outlet["species_mol"]["C2H6"] > 0
    assert result["plasma"]["unresolved_n2_total_outcome_opportunity_mol"] > 0
    assert result["plasma"]["unresolved_n2_total_outcome_threshold_energy_J"] > 0
    assert 0 < result["plasma"]["conditional_sub14_526eV_ground_n_opportunity_mol"] < (
        result["plasma"]["unresolved_n2_total_outcome_opportunity_mol"]
    )
    assert result["plasma"]["unresolved_n2_opportunity_is_propagated"] is False
    assert result["plasma"]["conditional_n2d_fate_status"] == "unavailable"
    assert result["plasma"]["conditional_n2d_fate_unavailable_reason"] == (
        "outside reviewed near-298 K diagnostic domain"
    )
    channels = result["heavy_chemistry"]
    assert channels["source_extents_mol_by_channel"]
    assert set(channels["source_extents_mol_by_channel"]) == set(channels["active_fragment_channel_ids"])
    assert all(channel.startswith("CH4_") for channel in channels["active_fragment_channel_ids"])
    assert channels["excluded_fragment_channels_by_id"]["N2_dissociation_outcome_000"]["reason"] == (
        "unresolved_electronic_product_state"
    )
    no_carbon = run_dry_flow_probe(discharge_on=False, methane_fraction=0)
    assert no_carbon["outlet_packet_diagnostics"]["methane_conversion_fraction"] is None
    assert no_carbon["outlet_packet_diagnostics"][
        "methane_electron_fragment_source_fraction_of_feed"
    ] is None
    assert no_carbon["outlet_packet_diagnostics"]["carbon_yield_fraction_by_species"] == {}
    assert no_carbon["outlet_packet_diagnostics"]["product_carbon_selectivity_fraction_by_species"] == {}
    assert no_carbon["outlet_packet_diagnostics"]["product_hydrogen_selectivity_fraction_by_species"] == {}
    assert no_carbon["plasma"]["unresolved_n2_total_outcome_opportunity_mol"] == 0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"discharge_on": False},
        {"methane_fraction": 0.0},
        {"methane_fraction": 1.0},
    ],
)
@pytest.mark.parametrize(
    "source_forcing", ["impulse_substeps", "inventory_limited", "evolving_eedf"]
)
def test_dry_flow_probe_negative_controls(kwargs: dict, source_forcing: str) -> None:
    result = run_dry_flow_probe(substeps=4, source_forcing=source_forcing, **kwargs)
    assert result["outlet_hcn_mol"] < 1e-12
    assert result["material_ledger"]["closed"] is True
    assert result["energy_ledger"]["closed_upper_bound"] is True


def test_dry_flow_probe_rejects_invalid_conditions() -> None:
    with pytest.raises(ValueError, match="fraction"):
        run_dry_flow_probe(methane_fraction=1.1)
    with pytest.raises(ValueError, match="substeps"):
        run_dry_flow_probe(substeps=0)
    with pytest.raises(ValueError, match="maximum-step"):
        run_dry_flow_probe(integration_max_step_fraction=0)
    with pytest.raises(ValueError, match="Flow reference temperature"):
        run_dry_flow_probe(flow_reference_temperature_K=0)
    with pytest.raises(ValueError, match="Uniform gas scenario temperature"):
        run_dry_flow_probe(uniform_gas_temperature_K=300)


def test_flow_reference_changes_hot_packet_time_consistently() -> None:
    cold_reference = run_dry_flow_probe(discharge_on=False, flow_reference_temperature_K=273.15)
    warm_reference = run_dry_flow_probe(discharge_on=False, flow_reference_temperature_K=298.15)
    cold_time = cold_reference["conditions"]["simulated_hot_packet_residence_s"]
    warm_time = warm_reference["conditions"]["simulated_hot_packet_residence_s"]
    assert cold_time / warm_time == pytest.approx(273.15 / 298.15)
    assert cold_reference["inlet_total_mol"] == pytest.approx(warm_reference["inlet_total_mol"])
    assert cold_reference["outlet_hcn_mol"] == pytest.approx(0)
    assert warm_reference["outlet_hcn_mol"] == pytest.approx(0)


def test_gas_internal_energy_uses_cantera_species_molar_units() -> None:
    prepared, methane = _gri_inventory(1.0, 1.0)
    report = gri_ideal_gas_internal_energy_inventory(methane, prepared, 400.0)
    enthalpy = gri_ideal_gas_enthalpy_inventory(methane, prepared, 400.0)
    phase = ct.Solution("gri30.yaml")
    phase.TP = 400.0, ct.one_atm
    expected_J = phase.partial_molar_int_energies[phase.species_index("CH4")] / 1000.0
    assert report["covered_internal_energy_J"] == pytest.approx(expected_J)
    assert report["all_nonzero_species_covered"] is True
    assert enthalpy["covered_enthalpy_J"] == pytest.approx(
        phase.partial_molar_enthalpies[phase.species_index("CH4")] / 1000.0
    )
    assert enthalpy["covered_enthalpy_J"] - report["covered_internal_energy_J"] == pytest.approx(
        ct.gas_constant * 400.0 / 1000.0
    )
