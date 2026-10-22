"""O(1D) first-collision audit, independent of historical HCN observations."""

import math

import pytest

from app.prebiotic.state_fate import (
    o1d_fate_diagnostic,
    project_thermalized_co2_o1d_events,
)


_ATOMS = {
    "O": {"O": 1}, "N2": {"N": 2}, "CO2": {"C": 1, "O": 2},
    "O2": {"O": 2}, "H2O": {"H": 2, "O": 1}, "H2": {"H": 2},
    "OH": {"O": 1, "H": 1}, "H": {"H": 1},
    "N2O": {"N": 2, "O": 1}, "NO": {"N": 1, "O": 1},
    "NH3": {"N": 1, "H": 3}, "NH2": {"N": 1, "H": 2},
}


def _bath() -> dict[str, float]:
    return {
        "N2": 0.919, "CO2": 0.05, "O2": 0.001,
        "H2O": 0.02, "H2": 0.01,
    }


def test_competing_state_fates_sum_to_one_and_have_subpulse_lifetime() -> None:
    result = o1d_fate_diagnostic(
        _bath(), gas_volume_L=25.0, gas_temperature_K=298.15,
        active_discharge_duration_s=0.001,
    )
    assert result["status"] == "conditional_fate_resolved"
    hazards = {
        row["bath_species"]: row["pseudo_first_order_hazard_s_inv"]
        for row in result["channels"]
    }
    physical = (
        hazards["N2"] + hazards["CO2"] + hazards["O2"]
        + 0.010 * hazards["H2O"] + 0.007 * hazards["H2"]
    )
    total = sum(hazards.values())
    assert result["physical_ground_O_fraction"] == pytest.approx(physical / total)
    assert sum(row["first_collision_fraction"] for row in result["channels"]) == pytest.approx(1)
    assert result["physical_ground_O_fraction"] + result["reactive_OH_H_fraction"] == pytest.approx(1)
    assert result["mean_first_collision_lifetime_s"] < 1e-5
    assert result["lifetime_to_active_discharge_ratio"] < 0.01
    assert result["products_injected_into_solver"] is False


def test_every_reviewed_fate_channel_balances_atoms_and_has_provenance() -> None:
    result = o1d_fate_diagnostic(_bath(), gas_volume_L=25, gas_temperature_K=298.15)
    assert "NASA/JPL Evaluation 17" in result["source"]
    for channel in result["channels"]:
        reactant_atoms: dict[str, int] = {"O": 1}
        for element, count in _ATOMS[channel["bath_species"]].items():
            reactant_atoms[element] = reactant_atoms.get(element, 0) + count
        product_atoms: dict[str, int] = {}
        for product, coefficient in channel["products"].items():
            for element, count in _ATOMS[product].items():
                product_atoms[element] = product_atoms.get(element, 0) + coefficient * count
        assert product_atoms == reactant_atoms
        if channel.get("minor_quench_products"):
            quenched_atoms: dict[str, int] = {}
            for product, coefficient in channel["minor_quench_products"].items():
                for element, count in _ATOMS[product].items():
                    quenched_atoms[element] = quenched_atoms.get(element, 0) + coefficient * count
            assert quenched_atoms == reactant_atoms
        if channel.get("alternate_products"):
            alternate_atoms: dict[str, int] = {}
            for product, coefficient in channel["alternate_products"].items():
                for element, count in _ATOMS[product].items():
                    alternate_atoms[element] = (
                        alternate_atoms.get(element, 0) + coefficient * count
                    )
            assert alternate_atoms == reactant_atoms


def test_reviewed_n2o_and_ammonia_branches_project_with_atom_closure() -> None:
    bath = _bath() | {"N2O": 0.002, "NH3": 0.002}
    fate = o1d_fate_diagnostic(bath, gas_volume_L=25, gas_temperature_K=298.15)
    assert fate["status"] == "conditional_fate_resolved"
    n2o = next(item for item in fate["channels"] if item["bath_species"] == "N2O")
    ammonia = next(item for item in fate["channels"] if item["bath_species"] == "NH3")
    assert n2o["rate_coefficient_cm3_molecule_s"] == pytest.approx(1.27e-10)
    assert n2o["branch_fraction_to_listed_products"] == pytest.approx(4.95 / 12.7)
    assert n2o["alternate_branch_fraction"] == pytest.approx(7.75 / 12.7)
    assert ammonia["products"] == {"OH": 1, "NH2": 1}
    assert fate["reactive_N2O_NH3_fraction"] == pytest.approx(
        n2o["first_collision_fraction"] + ammonia["first_collision_fraction"]
    )
    assert (
        fate["physical_ground_O_fraction"]
        + fate["reactive_OH_H_fraction"]
        + fate["reactive_N2O_NH3_fraction"]
    ) == pytest.approx(1)
    projection = project_thermalized_co2_o1d_events(
        events_mol=1e-5, neutral_gas_mol=bath,
        gas_volume_L=25, gas_temperature_K=298.15,
    )
    delta = projection["signed_gas_change_mol"]
    assert delta["N2O"] < 0 and delta["NH3"] < 0
    assert delta["N2"] > 0 and delta["O2"] > 0 and delta["NO"] > 0
    assert delta["NH2"] > 0
    assert all(abs(value) < 1e-18 for value in projection["atom_residual_mol"].values())
    assert projection["charge_ledger_closed"] is True


def test_reviewed_methane_total_rate_does_not_invent_isomer_branch() -> None:
    bath = _bath() | {"CH4": 0.002}
    fate = o1d_fate_diagnostic(bath, gas_volume_L=25, gas_temperature_K=298.15)
    assert fate["status"] == "conditional_fate_partially_resolved"
    methane = next(item for item in fate["channels"] if item["bath_species"] == "CH4")
    assert methane["rate_coefficient_cm3_molecule_s"] == pytest.approx(1.75e-10)
    assert methane["products"] is None
    assert fate["unresolved_CH4_fate_fraction"] == pytest.approx(
        methane["first_collision_fraction"]
    )
    assert "physical_ground_O_fraction" not in fate
    with pytest.raises(ValueError, match="fully resolved"):
        project_thermalized_co2_o1d_events(
            events_mol=1e-5, neutral_gas_mol=bath,
            gas_volume_L=25, gas_temperature_K=298.15,
        )


def test_measured_water_and_hydrogen_minor_quench_is_not_counted_as_radicals() -> None:
    result = o1d_fate_diagnostic(_bath(), gas_volume_L=1, gas_temperature_K=298.15)
    water = next(item for item in result["channels"] if item["bath_species"] == "H2O")
    hydrogen = next(item for item in result["channels"] if item["bath_species"] == "H2")
    assert water["rate_coefficient_cm3_molecule_s"] == pytest.approx(
        1.70e-10 * math.exp(36 / 298.15)
    )
    assert hydrogen["rate_coefficient_cm3_molecule_s"] == pytest.approx(1.35e-10)
    assert water["ground_O_yield_per_O1D_removed"] == pytest.approx(0.010)
    assert hydrogen["ground_O_yield_per_O1D_removed"] == pytest.approx(0.007)
    assert water["ground_O_yield_reported_plus_minus"] == pytest.approx(0.003)
    assert hydrogen["ground_O_yield_reported_plus_minus"] == pytest.approx(0.007)
    assert "Vranckx" in water["source"]
    assert result["reactive_OH_H_fraction"] == pytest.approx(
        water["first_collision_fraction"] * 0.990
        + hydrogen["first_collision_fraction"] * 0.993
    )


def test_dilution_changes_lifetime_but_not_first_collision_branches() -> None:
    dense = o1d_fate_diagnostic(_bath(), gas_volume_L=1, gas_temperature_K=298.15)
    sparse = o1d_fate_diagnostic(_bath(), gas_volume_L=2, gas_temperature_K=298.15)
    assert sparse["mean_first_collision_lifetime_s"] == pytest.approx(
        2 * dense["mean_first_collision_lifetime_s"]
    )
    assert sparse["physical_ground_O_fraction"] == pytest.approx(
        dense["physical_ground_O_fraction"]
    )


def test_measured_co_removal_shortens_lifetime_without_inventing_products() -> None:
    baseline = o1d_fate_diagnostic(
        _bath(), gas_volume_L=1, gas_temperature_K=298.15
    )
    result = o1d_fate_diagnostic(
        _bath() | {"CO": 0.001}, gas_volume_L=1, gas_temperature_K=298.15
    )
    assert result["status"] == "conditional_fate_partially_resolved"
    assert result["mean_first_collision_lifetime_s"] < baseline["mean_first_collision_lifetime_s"]
    co = next(row for row in result["channels"] if row["bath_species"] == "CO")
    assert co["kind"] == "unresolved_products"
    assert co["products"] is None
    assert "Davidson" in co["source"]
    assert co["rate_coefficient_cm3_molecule_s"] == pytest.approx(5.81e-11, rel=0.01)
    assert result["unresolved_CO_fate_fraction"] == pytest.approx(co["first_collision_fraction"])
    rate_low, rate_high = result["co_rate_parameter_range_cm3_molecule_s"]
    lifetime_low, lifetime_high = result["lifetime_range_from_CO_rate_only_s"]
    assert rate_low < co["rate_coefficient_cm3_molecule_s"] < rate_high
    assert lifetime_low < result["mean_first_collision_lifetime_s"] < lifetime_high
    assert "not a confidence interval" in result["co_rate_range_scope"]
    assert (
        result["known_physical_ground_O_fraction"]
        + result["known_reactive_OH_H_fraction"]
        + result["unresolved_CO_fate_fraction"]
    ) == pytest.approx(1)
    assert "physical_ground_O_fraction" not in result
    assert result["products_injected_into_solver"] is False


def test_measured_hcn_removal_resolves_quench_not_reactive_products() -> None:
    baseline = o1d_fate_diagnostic(
        _bath(), gas_volume_L=1, gas_temperature_K=298.15
    )
    result = o1d_fate_diagnostic(
        _bath() | {"HCN": 0.01}, gas_volume_L=1, gas_temperature_K=298.15
    )
    assert result["status"] == "conditional_fate_partially_resolved"
    assert result["mean_first_collision_lifetime_s"] < baseline["mean_first_collision_lifetime_s"]
    hcn = next(row for row in result["channels"] if row["bath_species"] == "HCN")
    assert hcn["products"] is None
    assert "Strekowski" in hcn["source"]
    assert hcn["rate_coefficient_cm3_molecule_s"] == pytest.approx(
        1.08e-10 * math.exp(105 / 298.15)
    )
    ground_yield = 0.150 * math.exp(200 / 298.15)
    assert hcn["ground_O_yield_per_O1D_removed"] == pytest.approx(ground_yield)
    assert result["hcn_ground_O_quench_fraction_of_all_collisions"] == pytest.approx(
        hcn["first_collision_fraction"] * ground_yield
    )
    assert result["unresolved_HCN_reactive_fraction"] == pytest.approx(
        hcn["first_collision_fraction"] * (1 - ground_yield)
    )
    assert (
        result["known_physical_ground_O_fraction"]
        + result["known_reactive_OH_H_fraction"]
        + result["unresolved_HCN_reactive_fraction"]
    ) == pytest.approx(1)
    assert "physical_ground_O_fraction" not in result
    assert result["products_injected_into_solver"] is False


def test_co_and_hcn_unresolved_fractions_remain_separate() -> None:
    result = o1d_fate_diagnostic(
        _bath() | {"CO": 0.01, "HCN": 0.01},
        gas_volume_L=1, gas_temperature_K=298.15,
    )
    assert result["status"] == "conditional_fate_partially_resolved"
    assert 0 < result["unresolved_CO_fate_fraction"] < 1
    assert 0 < result["unresolved_HCN_reactive_fraction"] < 1
    assert (
        result["known_physical_ground_O_fraction"]
        + result["known_reactive_OH_H_fraction"]
        + result["unresolved_CO_fate_fraction"]
        + result["unresolved_HCN_reactive_fraction"]
    ) == pytest.approx(1)


def test_unreviewed_bath_partner_fails_closed_without_renormalizing() -> None:
    result = o1d_fate_diagnostic(
        _bath() | {"CO": 0.001, "H2O2": 0.0001},
        gas_volume_L=1, gas_temperature_K=298.15
    )
    assert result["status"] == "unavailable"
    assert result["unreviewed_bath_species"] == ["H2O2"]
    assert result["unreviewed_bath_mole_fraction"] == pytest.approx(0.0001 / 1.0011)
    assert "physical_ground_O_fraction" not in result
    assert "mean_first_collision_lifetime_s" not in result
    stress = result["unreviewed_rate_stress"]
    assert [item["hypothetical_rate_cm3_molecule_s"] for item in stress] == [1e-9, 1e-8]
    assert 0 < stress[0]["hypothetical_unknown_hazard_fraction"] < stress[1][
        "hypothetical_unknown_hazard_fraction"
    ] < 1
    assert "not measured" in result["unreviewed_rate_stress_scope"]
    assert result["products_injected_into_solver"] is False


def test_unreviewed_rate_stress_is_not_a_volume_or_product_artifact() -> None:
    mixture = _bath() | {"H2O2": 0.0001}
    first = o1d_fate_diagnostic(
        mixture, gas_volume_L=1, gas_temperature_K=298.15,
    )
    second = o1d_fate_diagnostic(
        mixture, gas_volume_L=2, gas_temperature_K=298.15,
    )
    assert second["reviewed_total_removal_hazard_s_inv"] == pytest.approx(
        first["reviewed_total_removal_hazard_s_inv"] / 2
    )
    for left, right in zip(
        first["unreviewed_rate_stress"], second["unreviewed_rate_stress"], strict=True,
    ):
        assert left["hypothetical_unknown_hazard_fraction"] == pytest.approx(
            right["hypothetical_unknown_hazard_fraction"]
        )
    assert "channels" not in first


def test_out_of_domain_temperature_does_not_reuse_298k_rates() -> None:
    result = o1d_fate_diagnostic(
        _bath(), gas_volume_L=1, gas_temperature_K=373.15
    )
    assert result["status"] == "unavailable"
    assert "physical_ground_O_fraction" not in result


def test_isolated_co2_o1d_projection_closes_atoms_without_injecting_yield() -> None:
    amount = 1e-5
    projected = project_thermalized_co2_o1d_events(
        events_mol=amount, neutral_gas_mol=_bath(),
        gas_volume_L=25, gas_temperature_K=298.15,
    )
    delta = projected["signed_gas_change_mol"]
    assert delta["CO2"] == pytest.approx(-amount)
    assert delta["CO"] == pytest.approx(amount)
    assert delta["O"] > 0
    assert delta["H2O"] < 0
    assert delta["H2"] < 0
    assert delta["OH"] > 0
    assert all(abs(value) < 1e-18 for value in projected["atom_residual_mol"].values())
    assert projected["atom_ledger_closed"] is True
    assert projected["charge_ledger_closed"] is True
    assert projected["electronic_excitation_energy_initially_carried_J"] == pytest.approx(
        amount * 1.967 * 96485.332, rel=0.001,
    )
    assert projected["co2_to_co_plus_o1d_standard_enthalpy_requirement_eV_per_event"] == pytest.approx(
        7.48, rel=0.001,
    )
    assert projected["electron_collision_energy_allocation_checked"] is False
    assert projected["excitation_energy_booked_as_thermal_heat"] is False
    assert projected["complete_thermal_energy_balance_closed"] is False
    assert projected["products_injected_into_apparatus"] is False


def test_isolated_co2_o1d_projection_refuses_incomplete_or_depleted_bath() -> None:
    with pytest.raises(ValueError, match="fully resolved"):
        project_thermalized_co2_o1d_events(
            events_mol=1e-5, neutral_gas_mol=_bath() | {"H2O2": 1e-4},
            gas_volume_L=25, gas_temperature_K=298.15,
        )
    with pytest.raises(ValueError, match="1% frozen-bath CO2"):
        project_thermalized_co2_o1d_events(
            events_mol=1e-3, neutral_gas_mol=_bath(),
            gas_volume_L=25, gas_temperature_K=298.15,
        )


@pytest.mark.parametrize("volume,temperature", [(0, 298.15), (1, float("nan"))])
def test_invalid_conditions_are_rejected(volume: float, temperature: float) -> None:
    with pytest.raises(ValueError):
        o1d_fate_diagnostic(_bath(), gas_volume_L=volume, gas_temperature_K=temperature)
