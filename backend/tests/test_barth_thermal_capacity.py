"""The water thermal screen cannot be mistaken for a closed heat balance."""

import pytest

from scripts.audit_barth_thermal_capacity import audit_barth_water_thermal_capacity


def test_barth_water_screen_exposes_isothermal_energy_requirement_without_assay() -> None:
    report = audit_barth_water_thermal_capacity()
    assert report["feed_card_contains_hcn_assay"] is False
    assert report["complete_thermal_energy_balance_closed"] is False
    assert report["actual_heat_rejection_inferred_J"] is None
    assert "not measured" in report["nominal_discharge_energy_basis"]
    assert report["declared_spark_duty_fraction"] == pytest.approx(0.1)
    assert report["declared_average_deposited_power_W"] == pytest.approx(1.47)
    assert report["declared_total_deposited_energy_J"] == pytest.approx(81_585)
    assert 205 < report["initial_water_heat_capacity_J_K"] < 212
    assert 15_000 < report["liquid_sensible_enthalpy_to_saturation_J"] < 16_000
    assert report["saturation_temperature_at_initial_pressure_K"] > 372
    assert report["declared_input_to_liquid_sensible_capacity_ratio"] > 5


def test_barth_water_screen_rejects_duration_outside_assayed_run() -> None:
    with pytest.raises(ValueError, match="within the reported run"):
        audit_barth_water_thermal_capacity(duration_min=926)
    with pytest.raises(ValueError, match="positive duration"):
        audit_barth_water_thermal_capacity(duration_min=0)
    minute = audit_barth_water_thermal_capacity(duration_min=1)
    assert minute["declared_duration_s"] == 60
    assert minute["declared_total_deposited_energy_J"] == pytest.approx(88.2)
