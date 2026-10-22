"""Private IST-Lisbon H2 intake and opt-in EEDF/chemistry boundaries."""

from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from app.prebiotic import boltzmann
from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.boltzmann import (
    COMMUNITY_CH4_DATASET_ID,
    IST_H2_COMMUNITY_CH4_DATASET_ID,
    collision_dataset_status,
    solve_self_consistent_eedf,
)
from app.prebiotic.collision_intake import normalize_lxcat_bolsig_text
from app.prebiotic.coupled_plasma import _channel_identity
from app.prebiotic.heavy_gas import _neutral_fragment_extents, _pavan_source_product_name
from app.prebiotic.schema import BoltzmannEEDFRequest, PrebioticApparatusRequest
from scripts.audit_candidate_methane_collision_set import audit_candidate_collision_set


H2_CANDIDATE = (
    Path(__file__).resolve().parents[1]
    / "data" / "collision" / "private" / "IST_Lisbon_H2_2026-09-23"
    / "Cross section.txt"
)
H2_SHA256 = "2721e8cb5b0683b614c9094c279c104f38e987fcb29455fb168cb39e0e511468"


def test_hydrogen_candidate_is_structurally_parseable_without_activation() -> None:
    assert hashlib.sha256(H2_CANDIDATE.read_bytes()).hexdigest() == H2_SHA256
    report = audit_candidate_collision_set(
        H2_CANDIDATE, "H2", normalize_lxcat_export=True
    )
    assert report["candidate_sha256"] == H2_SHA256
    assert report["process_count"] == 22
    assert report["kind_counts"] == {
        "ELASTIC": 1, "EXCITATION": 19, "IONIZATION": 2
    }
    assert report["passes_preliminary_structure_check"] is True
    assert len(report["normalization_log"]) == 3
    assert report["normalization_log"][0]["kind"] == "trailing_reaction_annotation_removed"
    assert "H(1s)+H(2p)" in report["explicit_neutral_fragment_labels"]
    assert report["active_solver_changed"] is False


def test_hydrogen_annotation_normalizer_fails_closed_on_ambiguous_arrows() -> None:
    with pytest.raises(ValueError, match="Ambiguous H2 collision label"):
        normalize_lxcat_bolsig_text("H2 -> H -> H\n", "H2")


def test_hydrogen_annotation_normalizer_preserves_numeric_tables() -> None:
    raw = "H2 -> H + H+ (H2(X) -> H(1S) + H+)\r\n1.0 2e-21\r\n"
    normalized, changes = normalize_lxcat_bolsig_text(raw, "H2")
    assert normalized == "H2 -> H + H+\r\n1.0 2e-21\r\n"
    assert changes == [
        {
            "line_number": 1,
            "kind": "trailing_reaction_annotation_removed",
            "before": "H2 -> H + H+ (H2(X) -> H(1S) + H+)",
            "after": "H2 -> H + H+",
        }
    ]


def test_excited_hydrogen_fragments_are_not_silently_grounded() -> None:
    identity = _channel_identity({"H": 2}, "excitation", "H(1s)+H(2p)")
    assert identity["atoms_and_product_charge_balanced"] is True
    assert identity["reactive"] is True
    assert identity["products"][0]["atomic_h_ground_state"] is True
    assert identity["products"][0]["excited"] is False
    assert identity["products"][1]["excited"] is True

    source = {
        "id": "synthetic_h2_excited_fragment",
        "target": "H2",
        "product_label": "H(1s)+H(2p)",
        "uncoupled_active_window_opportunity_mol": 0.1,
        "independent_threshold_energy_lower_bound_J": 1000.0,
        "threshold_eV": 10.0,
        **identity,
    }
    prepared = {"index": {"H2": 0, "H": 1}}
    initial = np.array([1.0, 0.0])
    injected, accepted, _, allocation, warnings = _neutral_fragment_extents(
        initial, prepared, {"source_terms": [source]}, 1e6
    )
    assert injected == pytest.approx(initial)
    assert accepted == []
    assert allocation["excluded_fragment_channels"][0]["channel_id"] == source["id"]
    assert allocation["excluded_fragment_channels"][0]["reason"] == "unresolved_excited_state"
    assert any("unresolved excited" in warning for warning in warnings)
    assert _pavan_source_product_name(source, identity["products"][0], {"H", "H2"}) == "H"
    assert _pavan_source_product_name(source, identity["products"][1], {"H", "H2"}) is None

    ground_identity = _channel_identity({"H": 2}, "excitation", "H(1s)+H(1s)")
    ground_source = source | {
        "id": "synthetic_h2_ground_fragments",
        "product_label": "H(1s)+H(1s)",
        **ground_identity,
    }
    ground_injected, ground_accepted, _, _, _ = _neutral_fragment_extents(
        initial, prepared, {"source_terms": [ground_source]}, 1e6
    )
    assert ground_injected == pytest.approx([0.9, 0.2])
    assert len(ground_accepted) == 1


def test_morgan_water_ev_note_is_not_mistaken_for_excited_oh() -> None:
    identity = _channel_identity({"H": 2, "O": 1}, "excitation", "H+OH(7eV)")
    assert identity["atoms_and_product_charge_balanced"] is True
    assert identity["products"][1]["excited"] is False
    assert identity["products"][1]["energy_annotation_eV"] == 7.0
    source = {
        "id": "synthetic_morgan_water_ground_fragment",
        "target": "H2O",
        "product_label": "H+OH(7eV)",
        "uncoupled_active_window_opportunity_mol": 0.1,
        "independent_threshold_energy_lower_bound_J": 1000.0,
        "threshold_eV": 7.0,
        **identity,
    }
    prepared = {"index": {"H2O": 0, "H": 1, "OH": 2}}
    initial = np.array([1.0, 0.0, 0.0])
    injected, accepted, _, _, _ = _neutral_fragment_extents(
        initial, prepared, {"source_terms": [source]}, 1e6
    )
    assert injected == pytest.approx([0.9, 0.1, 0.1])
    assert len(accepted) == 1

    mismatched = source | {"threshold_eV": 9.0}
    rejected, accepted, _, allocation, _ = _neutral_fragment_extents(
        initial, prepared, {"source_terms": [mismatched]}, 1e6
    )
    assert rejected == pytest.approx(initial)
    assert accepted == []
    assert allocation["excluded_fragment_channels"][0]["reason"] == "energy_annotation_threshold_mismatch"


def test_charged_hydrogen_fragments_are_not_neutralized_in_primary_branch() -> None:
    identity = _channel_identity({"H": 2}, "ionization", "H+H^+")
    assert identity["atoms_and_product_charge_balanced"] is True
    source = {
        "id": "synthetic_h2_dissociative_ionization",
        "target": "H2",
        "product_label": "H+H^+",
        "uncoupled_active_window_opportunity_mol": 0.1,
        "independent_threshold_energy_lower_bound_J": 1000.0,
        "threshold_eV": 16.0,
        **identity,
    }
    initial = np.array([1.0, 0.0])
    injected, accepted, _, allocation, _ = _neutral_fragment_extents(
        initial, {"index": {"H2": 0, "H": 1}}, {"source_terms": [source]}, 1e6
    )
    assert injected == pytest.approx(initial)
    assert accepted == []
    assert allocation["excluded_fragment_channels"][0]["reason"] == (
        "charged_products_require_explicit_ion_chemistry"
    )


def test_opt_in_h2_whole_set_replacement_and_field_domain() -> None:
    assert hashlib.sha256(H2_CANDIDATE.read_bytes()).hexdigest() == H2_SHA256
    assert boltzmann.IST_H2_COMMUNITY_CH4_DATASET_SHA256 == hashlib.sha256(
        (
            boltzmann.DATASET_SHA256
            + H2_SHA256
            + boltzmann.COMMUNITY_CH4_SHA256
            + boltzmann.PHELPS_N2_SHA256
            + boltzmann.TABATA_N2_SHA256
        ).encode()
    ).hexdigest()
    status = collision_dataset_status(IST_H2_COMMUNITY_CH4_DATASET_ID)
    assert status["available"] is True
    assert status["coverage"]["process_count"] == 82
    assert status["coverage"]["targets"]["H2"]["process_count"] == 22
    assert status["coverage"]["targets"]["H2"]["has_explicit_neutral_fragment_source"] is True
    assert status["coverage"]["targets"]["H2"]["has_resolved_ground_state_fragment_source"] is False
    assert status["coverage"]["targets"]["H2O"]["has_resolved_ground_state_fragment_source"] is True
    processes = boltzmann._load_processes(IST_H2_COMMUNITY_CH4_DATASET_ID)
    assert Counter(item["target"] for item in processes)["H2"] == 22
    assert Counter(item["target"] for item in processes)["CH4"] == 16

    request = {
        "collision_set": IST_H2_COMMUNITY_CH4_DATASET_ID,
        "mole_fractions": {"H2": 1.0},
        "gas_temperature_K": 400,
        "reduced_electric_field_Td": 200,
        "energy_max_eV": 60,
        "energy_cells": 400,
        "initial_electron_temperature_eV": 2,
        "maximum_iterations": 200,
        "relative_tolerance": 1e-6,
        "acknowledge_two_term_limit": True,
    }
    new = solve_self_consistent_eedf(BoltzmannEEDFRequest.model_validate(request))
    old = solve_self_consistent_eedf(
        BoltzmannEEDFRequest.model_validate(request | {"collision_set": COMMUNITY_CH4_DATASET_ID})
    )
    assert new["solver"]["converged"] is True
    assert new["eedf"]["mean_energy_eV"] != pytest.approx(old["eedf"]["mean_energy_eV"], rel=1e-4)
    same_set = new["same_collision_set_reference"]
    assert same_set["available_for_this_mixture"] is True
    assert same_set["collision_set_sha256"] == boltzmann.IST_H2_COMMUNITY_CH4_DATASET_SHA256
    assert same_set["selected_eedf_process_count"] == 22
    assert abs(same_set["relative_mean_energy_difference"]) < 0.02
    assert same_set["rate_comparison"]["matched_process_count"] == 21
    assert same_set["rate_comparison"]["material_within_factor_ten_fraction"] == 1.0
    assert same_set["private_collision_data_persisted"] is False
    assert any(rate["target"] == "H2" and rate["product_label"] == "H(1s)+H(2p)"
               for rate in new["rate_coefficients"])
    for field in (0.5, 201):
        with pytest.raises(ValueError, match="1-200 Td"):
            BoltzmannEEDFRequest.model_validate(request | {"reduced_electric_field_Td": field})


def test_five_gas_same_collision_set_two_term_engine_comparison() -> None:
    request = BoltzmannEEDFRequest.model_validate({
        "collision_set": IST_H2_COMMUNITY_CH4_DATASET_ID,
        "mole_fractions": {"H2": 0.35, "H2O": 0.25, "CH4": 0.2, "NH3": 0.1, "N2": 0.1},
        "gas_temperature_K": 400,
        "reduced_electric_field_Td": 200,
        "energy_max_eV": 60,
        "energy_cells": 400,
        "initial_electron_temperature_eV": 2,
        "maximum_iterations": 200,
        "relative_tolerance": 1e-6,
        "acknowledge_two_term_limit": True,
    })
    result = solve_self_consistent_eedf(request)
    same_set = result["same_collision_set_reference"]
    assert same_set["available_for_this_mixture"] is True
    assert same_set["selected_eedf_process_count"] == 81
    assert same_set["collision_set_dataset_id"] == IST_H2_COMMUNITY_CH4_DATASET_ID
    assert abs(same_set["relative_mean_energy_difference"]) < 0.05
    rate_check = same_set["rate_comparison"]
    assert rate_check["matched_process_count"] == 77
    assert rate_check["material_process_count"] > 50
    assert rate_check["material_within_factor_ten_fraction"] == 1.0
    assert {row["target"] for row in rate_check["largest_reactive_channel_by_target"]} == {
        "H2", "H2O", "CH4", "NH3", "N2"
    }
    assert "not an external swarm experiment" in same_set["comparison_scope"]


def test_low_field_hydrogen_nonfinite_rate_is_rejected() -> None:
    request = BoltzmannEEDFRequest.model_validate({
        "collision_set": IST_H2_COMMUNITY_CH4_DATASET_ID,
        "mole_fractions": {"H2": 1.0},
        "gas_temperature_K": 400,
        "reduced_electric_field_Td": 1,
        "energy_max_eV": 60,
        "energy_cells": 400,
        "initial_electron_temperature_eV": 2,
        "maximum_iterations": 200,
        "relative_tolerance": 1e-6,
        "acknowledge_two_term_limit": True,
    })
    with pytest.raises(ValueError, match="nonfinite or negative collision rate"):
        solve_self_consistent_eedf(request)


def test_opt_in_h2_excited_sources_do_not_enter_ground_state_apparatus() -> None:
    profile = Path(__file__).resolve().parents[1] / "data" / "benchmarks" / "stribling_miller_measurement_profile.json"
    payload = json.loads(profile.read_text(encoding="utf-8"))
    payload["apparatus"]["cycles"] = 1
    payload["spark"]["boltzmann"]["collision_set"] = IST_H2_COMMUNITY_CH4_DATASET_ID
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    snapshot = result["cycles"][0]["plasma_snapshot"]
    assert any(source["target"] == "H2" and source["product_label"] == "H(1s)+H(2p)"
               for source in snapshot["source_terms"])
    heavy = result["cycles"][0]["heavy_species_snapshot"]
    assert any(item["reason"] == "unresolved_excited_state" and item["product_label"] == "H(1s)+H(2p)"
               for item in heavy["collision_energy_allocation"]["excluded_fragment_channels"])
    assert not any(item["target"] == "H2" for item in heavy["source_extents"])
    assert result["material_ledger"]["closed"] is True
    assert result["plasma_model"]["propagated_charged_branch"]["material_ledger"]["closed"] is True
    assert result["scientific_firewall"]["quantitative_yield_claim_allowed"] is False
