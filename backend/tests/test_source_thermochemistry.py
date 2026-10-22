"""Check necessary fragment-energy safeguards without loading HCN assay data."""

from copy import deepcopy

from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.schema import PrebioticApparatusRequest
from scripts.audit_source_thermochemistry import (
    audit_barth_feed_source_thermochemistry,
    audit_source_thermochemistry,
)
from scripts.probe_barth_feed_pathway import make_probe_payload


def test_accepted_fragment_enthalpy_floor_and_n2_parent_event_capacity() -> None:
    report = audit_barth_feed_source_thermochemistry()
    assert report["status"] == "accepted_source_thermochemical_floor_passed"
    assert report["assay_loaded"] is False
    assert report["deficits"] == []
    assert report["n2_parent_coverage_deficits"] == []
    assert len(report["n2_outcome_parent_manifold_coverage"]) == 2
    assert report["n2_parent_minimum_capacity_to_outcome_ratio"] is None
    assert report["material_ledger_closed"] is True
    assert report["complete_thermal_energy_balance_closed"] is False

    n2_outcomes = [
        row for row in report["channels"]
        if row["channel_id"] == "N2_dissociation_outcome_000"
    ]
    assert n2_outcomes == []


def test_source_audit_fails_for_energy_and_parent_event_deficits() -> None:
    payload = make_probe_payload(
        duration_s=60, cycles=1, substeps=2,
        gas_water_kla_per_s=1e-4, gas_water_solver="coupled_bdf",
        primary_source_forcing="inventory_limited",
        formaldehyde_hydration_model="winkelman_2002_pH4_to_5_extrapolation_diagnostic",
    )
    baseline = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    temperature_K = float(payload["apparatus"]["discharge_temperature"]["value"])

    underfunded = deepcopy(baseline)
    h2o = next(
        source
        for source in underfunded["plasma_model"]["heavy_species_cycle_summaries"][0]["source_extents"]
        if source["target"] == "H2O"
    )
    h2o["energy_cost_eV_per_event"] = 0.1
    energy_report = audit_source_thermochemistry(underfunded, temperature_K)
    assert energy_report["status"] == "accepted_source_thermochemical_floor_failed"
    assert any(item["target"] == "H2O" for item in energy_report["deficits"])

    missing_parents = deepcopy(baseline)
    # Synthetic accepted event tests the audit's parent-capacity detector;
    # the actual state-unresolved Tabata outcome is never injected.
    for step in missing_parents["plasma_model"]["heavy_species_cycle_summaries"]:
        step["source_extents"].append({
            "channel_id": "synthetic_N2_outcome_test_only",
            "target": "N2", "products": ["N", "N"],
            "energy_cost_eV_per_event": 0.0,
            "energy_loss_accounted_by_parent_manifold": True,
            "accepted_extent_mol": 1e-6,
        })
    for step in missing_parents["plasma_model"]["cycle_summaries"]:
        for source in step["source_terms"]:
            if source["target"] == "N2" and source["role"] == "excitation":
                source["uncoupled_active_window_opportunity_mol"] = 0.0
    parent_report = audit_source_thermochemistry(missing_parents, temperature_K)
    assert parent_report["status"] == "accepted_source_thermochemical_floor_failed"
    assert len(parent_report["n2_parent_coverage_deficits"]) == 2
