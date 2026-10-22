"""Temperature sensitivity stays assay-free and conserves material/energy."""

import pytest

from scripts.audit_dry_temperature_sweep import TEMPERATURES_K, audit_dry_temperature_sweep


def test_dry_temperature_sweep_preserves_reactor_accounting() -> None:
    report = audit_dry_temperature_sweep()
    assert report["measured_hcn_assay_loaded"] is False
    assert [row["uniform_gas_temperature_K"] for row in report["rows"]] == list(TEMPERATURES_K)
    for row in report["rows"]:
        assert row["material_ledger_closed"] is True
        assert row["collision_energy_upper_bound_closed"] is True
        assert row["n2_product_state_source_excluded"] is True
        assert row["hcn_amount_below_1e-12_mol_reporting_floor"] is True
        assert row["hcn_carbon_yield_fraction_below_reporting_floor"] < row["c2h6_carbon_yield_fraction"]
    assert report["rows"][1]["equals_downstream_probe_value"] is True
    assert report["rows"][0]["equals_downstream_probe_value"] is False
    assert report["rows"][2]["equals_downstream_probe_value"] is False
    assert report["rows"][0]["ideal_gas_packet_residence_s"] / report["rows"][2][
        "ideal_gas_packet_residence_s"
    ] == pytest.approx(700.0 / 450.0)
