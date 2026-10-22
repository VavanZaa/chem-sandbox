"""Assay-free dry N2 support diagnostics; never a fragment-source test."""

import pytest

from scripts.audit_dry_n2_source_energy_support import audit_dry_n2_source_energy_support


def test_dry_energy_window_is_partial_and_cannot_inject_atoms() -> None:
    audit = audit_dry_n2_source_energy_support()
    assert audit["hcn_assay_read"] is False
    assert audit["atom_source_injected"] is False
    assert audit["source_state_branch_inferred"] is False
    assert audit["hcn_yield_predicted"] is False
    assert audit["feed_mole_fractions"] == {"CH4": 0.01, "N2": 0.99}
    fields = audit["field_sweep"]
    assert [item["reduced_electric_field_Td"] for item in fields] == [189.0, 200.0, 300.0, 435.0]
    assert all(item["eedf_converged"] for item in fields)
    assert all(item["product_state_branching_resolved"] is False for item in fields)
    assert all(
        0 < item["product_state_energy_support"]["fraction_of_computed_rate_below_measured_window"] < 1
        for item in fields
    )
    assert fields[1]["product_state_energy_support"]["fraction_of_computed_rate_below_measured_window"] == pytest.approx(
        0.558758, abs=0.005
    )
    grid = audit["grid_sweep_at_200_Td"]
    assert all(item["eedf_converged"] for item in grid)
    rates = [item["total_n2_neutral_dissociation_rate_coefficient_m3_s"] for item in grid]
    assert max(rates) / min(rates) < 1.01
