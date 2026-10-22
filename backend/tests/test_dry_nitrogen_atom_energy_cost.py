"""Read-only N2 energy-cost screen does not invent atom branching or HCN."""

import pytest

from scripts.audit_dry_nitrogen_atom_energy_cost import (
    _intrinsic_dissociation_cost,
    audit_dry_nitrogen_atom_energy_cost,
)


def test_pure_n2_energy_cost_grid_and_independent_scope() -> None:
    audit = audit_dry_nitrogen_atom_energy_cost()
    rows = audit["pure_n2_reference_condition_grid_check"]
    costs = [row["electron_drift_eV_per_total_n2_outcome"] for row in rows]
    assert all(row["eedf_converged"] for row in rows)
    assert max(costs) / min(costs) < 1.02
    assert costs[-1] / 2 == pytest.approx(215.0555, rel=0.015)
    assert all(
        row["conditional_eV_per_atom_if_one_ground_atom_per_outcome"] == pytest.approx(
            2 * row["conditional_eV_per_atom_if_two_ground_atoms_per_outcome"]
        ) for row in rows
    )
    assert audit["independent_comparison"]["like_for_like_with_dry_ch4_n2_reactor"] is False
    assert audit["independent_comparison"]["shares_cosby_total_cross_section_ancestry"] is True
    assert audit["nitrogen_atom_source_injected"] is False
    assert audit["hcn_assay_loaded"] is False
    assert audit["hcn_yield_predicted"] is False


def test_dry_scenario_is_separate_and_branch_unresolved() -> None:
    audit = audit_dry_nitrogen_atom_energy_cost()
    dry = audit["dry_one_percent_methane_average_field_scenario"]
    assert dry["methane_mole_fraction"] == pytest.approx(0.01)
    assert dry["reduced_electric_field_Td"] == pytest.approx(189.0)
    assert dry["eedf_converged"] is True
    assert dry["electron_drift_eV_per_total_n2_outcome"] > 0
    assert dry["product_state_branch_known"] is False
    assert dry["actual_ground_atom_energy_cost_determined"] is False


def test_energy_cost_requires_nitrogen_target() -> None:
    with pytest.raises(ValueError, match="N2 target fraction"):
        _intrinsic_dissociation_cost(
            methane_fraction=1.0, temperature_K=300.0, field_Td=140.0,
            energy_max_eV=60.0, energy_cells=400,
        )
