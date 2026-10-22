"""Flow-reference scenarios are independent of the HCN assay."""

import pytest

from scripts.audit_dry_flow_reference import (
    GF40_LISTED_REFERENCE_TEMPERATURES_K,
    PROJECT_LEGACY_ASSUMPTION_K,
    audit_dry_flow_reference,
)


def test_gf40_reference_sensitivity_is_physically_consistent() -> None:
    report = audit_dry_flow_reference()
    assert report["measured_hcn_assay_loaded"] is False
    assert report["actual_experiment_reference_temperature_identified"] is False
    rows = report["rows"]
    assert [row["reference_temperature_K"] for row in rows] == [
        *GF40_LISTED_REFERENCE_TEMPERATURES_K,
        PROJECT_LEGACY_ASSUMPTION_K,
    ]
    assert all(row["listed_GF40_option"] for row in rows[:-1])
    assert rows[-1]["listed_GF40_option"] is False
    for row in rows:
        assert row["material_ledger_closed"] is True
        assert row["collision_energy_upper_bound_closed"] is True
        assert row["n2_product_state_source_excluded"] is True
        assert row["hcn_amount_below_1e-12_mol_reporting_floor"] is True
        assert row["declared_energy_J_per_packet"] == pytest.approx(
            28.1 * row["ideal_gas_packet_residence_s"]
        )
    assert [row["inferred_feed_molar_flow_mol_s"] for row in rows] == sorted(
        (row["inferred_feed_molar_flow_mol_s"] for row in rows), reverse=True
    )
    assert [row["ideal_gas_packet_residence_s"] for row in rows] == sorted(
        row["ideal_gas_packet_residence_s"] for row in rows
    )
    assert rows[0]["ideal_gas_packet_residence_s"] / rows[-1][
        "ideal_gas_packet_residence_s"
    ] == pytest.approx(273.15 / 298.15)
