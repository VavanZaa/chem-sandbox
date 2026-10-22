"""Conditional localization must preserve ledgers and refuse yield validation."""

import pytest

from scripts.audit_dry_localized_exposure import (
    audit_dry_localization_bracket,
    run_dry_localized_exposure_probe,
)
from scripts.run_dry_flow_hcn_probe import run_dry_flow_probe


def test_unit_active_fraction_reduces_exactly_to_uniform_packet() -> None:
    localized = run_dry_localized_exposure_probe(
        active_volume_fraction=1.0, local_electron_density_m3=7.3e15, substeps=4
    )
    uniform = run_dry_flow_probe(
        reduced_electric_field_Td=435.0,
        electron_density_ceiling_m3=7.3e15,
        source_forcing="evolving_eedf",
        substeps=4,
    )
    assert localized["outlet"]["ch4_conversion_fraction"] == pytest.approx(
        uniform["outlet_packet_diagnostics"]["methane_conversion_fraction"], rel=1e-10
    )
    assert localized["outlet"]["c2h6_mol"] == pytest.approx(
        uniform["outlet_packet_diagnostics"]["species_mol"]["C2H6"], rel=1e-10
    )
    assert localized["ledger"]["material_closed"] is True
    assert localized["ledger"]["charge_closed"] is True
    assert localized["ledger"]["collision_energy_upper_bound_closed"] is True


def test_electrical_localization_bracket_is_conditional_not_validated() -> None:
    bracket = audit_dry_localization_bracket()
    assert bracket["hcn_assay_loaded"] is False
    assert bracket["product_fitted_streamer_fraction_used"] is False
    assert bracket["independent_occupancy_measured"] is False
    assert bracket["mixing_encounter_count_independently_measured"] is False
    assert [row["prediction"]["outlet"]["ch4_conversion_fraction"] for row in bracket["rows"]] == pytest.approx(
        [0.06204119249, 0.14688390673], rel=1e-5
    )
    for row in bracket["rows"]:
        prediction = row["prediction"]
        assert row["fraction_used_as_active_volume_sensitivity"] == pytest.approx(
            row["conditional_power_compatible_joint_fraction_upper"]
        )
        assert prediction["scenario"]["mixing_encounter_count_independently_measured"] is False
        assert prediction["scenario"]["varying_substeps_changes_the_physical_mixing_scenario"] is True
        assert prediction["quantitative_hcn_or_conversion_validation_allowed"] is False
        assert prediction["published_hcn_measurements_used_as_inputs"] is False
        assert prediction["outlet"]["hcn_below_1e-12_mol_reporting_floor"] is True
        assert prediction["n2_unresolved_product_state_source_excluded"] is True
        assert prediction["ledger"]["material_closed"] is True
        assert prediction["ledger"]["charge_closed"] is True
        assert prediction["ledger"]["collision_energy_upper_bound_closed"] is True


def test_localized_negative_controls_preserve_conservation() -> None:
    controls = (
        {"discharge_on": False},
        {"methane_fraction": 0.0},
        {"methane_fraction": 1.0},
    )
    for override in controls:
        result = run_dry_localized_exposure_probe(
            active_volume_fraction=0.03,
            local_electron_density_m3=1.5e15,
            substeps=2,
            **override,
        )
        assert result["outlet"]["hcn_below_1e-12_mol_reporting_floor"] is True
        assert result["n2_unresolved_product_state_source_excluded"] is True
        assert result["ledger"]["material_closed"] is True
        assert result["ledger"]["charge_closed"] is True
        assert result["ledger"]["collision_energy_upper_bound_closed"] is True


def test_unmeasured_encounter_count_changes_result() -> None:
    coarse = run_dry_localized_exposure_probe(
        active_volume_fraction=0.007876522468522865,
        local_electron_density_m3=7.3e15,
        substeps=4,
    )
    finer = run_dry_localized_exposure_probe(
        active_volume_fraction=0.007876522468522865,
        local_electron_density_m3=7.3e15,
        substeps=8,
    )
    assert abs(
        coarse["outlet"]["ch4_conversion_fraction"]
        - finer["outlet"]["ch4_conversion_fraction"]
    ) > 0.001
    assert coarse["scenario"]["mixing_encounter_count"] == 4
    assert finer["scenario"]["mixing_encounter_count"] == 8
    assert coarse["quantitative_hcn_or_conversion_validation_allowed"] is False
    assert finer["quantitative_hcn_or_conversion_validation_allowed"] is False
