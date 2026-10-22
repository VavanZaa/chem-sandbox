"""Thermal CH+N2 comparison must not become an HCN-yield calibration."""

from scripts.audit_dry_ch_n2_routes import audit_dry_ch_n2_routes


def test_dry_thermal_ch_n2_routes_are_assay_free_and_not_a_product_split() -> None:
    audit = audit_dry_ch_n2_routes()
    assert audit["uniform_gas_temperature_K_scenario"] == 555.15
    assert audit["gri_direct_rate_m3_kmol_s"] > 0
    assert audit["gri_association_effective_rate_m3_kmol_s"] > 0
    assert 0 < audit["direct_over_association_forward_rate_ratio"] < 1e-6
    assert audit["outlet_hcn_below_1e_12_mol_reporting_floor"] is True
    assert audit["n2_total_dissociation_source_excluded"] is True
    assert audit["direct_over_association_ratio_is_complete_product_branch"] is False
    assert audit["hcnn_ncn_hncn_fate_resolved"] is False
    assert audit["hcn_assay_loaded"] is False
    assert audit["solver_or_exposure_modified"] is False
