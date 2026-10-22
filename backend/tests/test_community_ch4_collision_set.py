"""Private complete-set intake checks; never copy contributor tables into tests."""

from collections import Counter
import hashlib
from pathlib import Path

import pytest

from app.prebiotic import boltzmann
from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.boltzmann import (
    COMMUNITY_CH4_DATASET_ID,
    COMMUNITY_CH4_DATASET_SHA256,
    COMMUNITY_CH4_PATH,
    COMMUNITY_CH4_SHA256,
    REACTIVE_N2_DATASET_ID,
    collision_dataset_status,
    solve_self_consistent_eedf,
)
from app.prebiotic.collision_intake import normalize_lxcat_processes
from app.prebiotic.coupled_plasma import _channel_identity
from app.prebiotic.schema import BoltzmannEEDFRequest, PrebioticApparatusRequest
from scripts.audit_candidate_methane_collision_set import audit_candidate_collision_set


def _request(dataset_id: str, cells: int = 400) -> BoltzmannEEDFRequest:
    return BoltzmannEEDFRequest.model_validate(
        {
            "collision_set": dataset_id,
            "mole_fractions": {"CH4": 0.2, "H2": 0.4, "H2O": 0.3, "NH3": 0.1},
            "gas_temperature_K": 400,
            "reduced_electric_field_Td": 300,
            "energy_max_eV": 60,
            "energy_cells": cells,
            "initial_electron_temperature_eV": 2,
            "maximum_iterations": 200,
            "relative_tolerance": 1e-6,
            "acknowledge_two_term_limit": True,
        }
    )


def test_private_community_ch4_is_complete_and_replaces_morgan_ch4() -> None:
    assert COMMUNITY_CH4_PATH.is_file()
    assert hashlib.sha256(COMMUNITY_CH4_PATH.read_bytes()).hexdigest() == COMMUNITY_CH4_SHA256
    assert COMMUNITY_CH4_DATASET_SHA256 == hashlib.sha256(
        (
            boltzmann.DATASET_SHA256
            + COMMUNITY_CH4_SHA256
            + boltzmann.PHELPS_N2_SHA256
            + boltzmann.TABATA_N2_SHA256
        ).encode()
    ).hexdigest()
    status = collision_dataset_status(COMMUNITY_CH4_DATASET_ID)
    assert status["available"] is True
    assert status["coverage"]["all_required_targets_complete"] is True
    methane = status["coverage"]["targets"]["CH4"]
    assert methane["process_count"] == 16
    assert methane["channel_counts"] == {
        "elastic": 1, "effective": 0, "momentum": 0,
        "excitation": 6, "dissociation": 0, "ionization": 7,
        "attachment": 2,
    }
    assert methane["neutral_dissociation_present"] is True
    assert methane["unresolved_dissociation_channel_count"] == 0
    assert set(methane["explicit_neutral_fragment_labels"]) == {
        "C+H2+H2", "CH+H2+H", "CH3+H", "CH2+H2"
    }
    all_processes = boltzmann._load_processes(COMMUNITY_CH4_DATASET_ID)
    assert Counter(item["target"] for item in all_processes)["CH4"] == 16
    assert len(all_processes) == 77  # 11 Morgan CH4 curves replaced by 16 Community curves
    for process in all_processes:
        if process["target"] != "CH4" or process["kind"] == "ELASTIC":
            continue
        identity = _channel_identity(
            {"C": 1, "H": 4}, process["kind"].lower(), process["product"]
        )
        assert identity["atoms_and_product_charge_balanced"], process["product"]


def test_audit_records_only_in_memory_normalization() -> None:
    report = audit_candidate_collision_set(
        COMMUNITY_CH4_PATH, "CH4", normalize_lxcat_export=True
    )
    assert report["candidate_sha256"] == COMMUNITY_CH4_SHA256
    assert report["passes_preliminary_structure_check"] is True
    assert len(report["normalization_log"]) == 15
    assert report["candidate_data_copied_or_installed"] is False
    assert report["active_solver_changed"] is False


def test_normalizer_rejects_conflicting_duplicate_points() -> None:
    with pytest.raises(ValueError, match="different cross sections"):
        normalize_lxcat_processes(
            [{"target": "CH4", "kind": "EXCITATION", "product": "CH3 H",
              "data": [(2.0, 1e-21), (2.0, 2e-21)]}],
            "CH4",
        )


def test_private_ch4_missing_or_tampered_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(boltzmann, "COMMUNITY_CH4_PATH", tmp_path / "missing.txt")
    assert collision_dataset_status(COMMUNITY_CH4_DATASET_ID)["available"] is False
    (tmp_path / "tampered.txt").write_text("not a collision set", encoding="utf-8")
    monkeypatch.setattr(boltzmann, "COMMUNITY_CH4_PATH", tmp_path / "tampered.txt")
    assert collision_dataset_status(COMMUNITY_CH4_DATASET_ID)["available"] is False
    with pytest.raises(ValueError, match="reviewed SHA-256"):
        boltzmann._load_processes(COMMUNITY_CH4_DATASET_ID)


def test_new_complete_set_changes_eedf_and_produces_balanced_methane_fragments() -> None:
    new = solve_self_consistent_eedf(_request(COMMUNITY_CH4_DATASET_ID))
    old = solve_self_consistent_eedf(_request(REACTIVE_N2_DATASET_ID))
    assert new["solver"]["converged"] is True
    assert new["collision_set"]["eedf_process_count"] == 51
    assert new["collision_set"]["post_solve_outcome_process_count"] == 0
    assert new["eedf"]["mean_energy_eV"] != pytest.approx(
        old["eedf"]["mean_energy_eV"], rel=1e-4
    )
    fragments = [
        item for item in new["rate_coefficients"]
        if item["target"] == "CH4" and item["product_label"] == "CH3+H"
    ]
    assert len(fragments) == 1
    assert fragments[0]["rate_coefficient_m3_s"] > 0


def test_complete_methane_set_propagates_through_one_apparatus_cycle() -> None:
    profile = (
        Path(__file__).resolve().parents[1]
        / "data" / "benchmarks" / "stribling_miller_measurement_profile.json"
    )
    import json

    payload = json.loads(profile.read_text(encoding="utf-8"))
    payload["apparatus"]["cycles"] = 1
    payload["spark"]["boltzmann"]["collision_set"] = COMMUNITY_CH4_DATASET_ID
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    snapshot = result["cycles"][0]["plasma_snapshot"]
    assert "Community database" in snapshot["collision_set"]["source"]
    methane_source = next(
        source for source in snapshot["source_terms"]
        if source["target"] == "CH4" and source["product_label"] == "CH3+H"
    )
    assert methane_source["atoms_and_product_charge_balanced"] is True
    assert methane_source["conditional_source_rate_mol_s"] > 0
    accepted = result["cycles"][0]["heavy_species_snapshot"]["source_extents"]
    assert any(
        item["target"] == "CH4" and item["products"] == ["CH3", "H"]
        and item["accepted_extent_mol"] > 0
        for item in accepted
    )
    assert result["material_ledger"]["closed"] is True
    assert result["plasma_model"]["propagated_charged_branch"]["material_ledger"]["closed"] is True
    assert result["scientific_firewall"]["quantitative_yield_claim_allowed"] is False
