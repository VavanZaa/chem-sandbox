"""The field sweep remains independent of held-out product observations."""

from scripts.audit_dry_field_sweep import FIELDS_TD, audit_dry_field_sweep


def test_dry_field_sweep_keeps_n2_unresolved_and_ledgers_closed() -> None:
    report = audit_dry_field_sweep()
    assert report["measured_hcn_assay_loaded"] is False
    assert [row["reduced_electric_field_Td"] for row in report["rows"]] == list(FIELDS_TD)
    for row in report["rows"]:
        assert row["material_ledger_closed"] is True
        assert row["collision_energy_upper_bound_closed"] is True
        assert row["n2_product_state_source_excluded"] is True
        assert row["hcn_amount_below_1e-12_mol_reporting_floor"] is True
        assert row["hcn_carbon_yield_fraction_below_reporting_floor"] < row["c2h6_carbon_yield_fraction"]
        assert row["unresolved_n2_total_outcome_opportunity_mol"] > 0
        assert 0 < row["conditional_sub14_526eV_ground_n_opportunity_mol"] < (
            row["unresolved_n2_total_outcome_opportunity_mol"]
        )
