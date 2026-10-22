from pathlib import Path

import pytest

from app.prebiotic.boltzmann import DATASET_PATH, DATASET_SHA256, _coverage
from scripts import audit_candidate_methane_collision_set as candidate_audit
from scripts.audit_candidate_methane_collision_set import (
    audit_candidate_collision_set,
    audit_candidate_methane_collision_set,
)


def test_existing_morgan_set_fails_fragment_intake_without_changing_solver() -> None:
    report = audit_candidate_methane_collision_set(DATASET_PATH)
    assert report["candidate_sha256"] == DATASET_SHA256
    assert report["eedf_input_coverage"] is True
    assert report["explicit_neutral_fragment_labels"] == []
    assert len(report["unresolved_dissociation_labels"]) == 4
    assert report["passes_preliminary_structure_check"] is False
    assert report["candidate_data_copied_or_installed"] is False
    assert report["active_solver_changed"] is False


def test_missing_candidate_is_not_silently_treated_as_coverage(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        audit_candidate_methane_collision_set(tmp_path / "missing.txt")


def test_generic_audit_does_not_mistake_suggested_ammonia_products_for_resolved_channels(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        candidate_audit.parser,
        "parse",
        lambda _: [
            {
                "target": "NH3",
                "kind": "EXCITATION",
                "product": "NH3(E4)",
                "comment": "Suggested dissociation channel NH2 + H",
                "data": [(7.48, 0.0), (8.0, 1e-21)],
            }
        ],
    )
    report = audit_candidate_collision_set(DATASET_PATH, "NH3")
    assert report["target"] == "NH3"
    assert report["explicit_neutral_fragment_labels"] == []
    assert report["suggested_dissociation_channel_labels_not_verified_as_products"] == ["NH3(E4)"]
    assert report["passes_preliminary_structure_check"] is False


def test_momentum_keyword_counts_as_eedf_input_and_neutral_labels_exclude_ions() -> None:
    processes = [
        {"target": "CH4", "kind": "MOMENTUM", "product": "CH4"},
        {"target": "CH4", "kind": "EXCITATION", "product": "CH3+H"},
        {"target": "CH4", "kind": "EXCITATION", "product": "CH3+H^-"},
        {"target": "CH4", "kind": "EXCITATION", "product": "CH4+"},
        {
            "target": "CH4",
            "kind": "EXCITATION",
            "product": "CH3+H",
            "comment": "Suggested dissociation channel; unverified branching",
        },
        {"target": "CH4", "kind": "IONIZATION", "product": "CH4^+"},
    ]
    coverage = _coverage(processes, ("CH4",))["targets"]["CH4"]
    assert coverage["complete"] is True
    assert coverage["momentum_transfer_present"] is True
    assert coverage["channel_counts"]["momentum"] == 1
    assert coverage["explicit_neutral_fragment_labels"] == ["CH3+H"]


def test_intake_rejects_incompatible_momentum_definitions(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    processes = [
        {"target": "NH3", "kind": kind, "data": [(0.0, 1e-20), (20.0, 1e-20)]}
        for kind in ("ELASTIC", "MOMENTUM")
    ]
    processes.extend(
        [
            {
                "target": "NH3",
                "kind": "EXCITATION",
                "product": "NH2+H",
                "data": [(0.0, 0.0), (20.0, 1e-21)],
            },
            {
                "target": "NH3",
                "kind": "IONIZATION",
                "product": "NH3^+",
                "data": [(0.0, 0.0), (20.0, 1e-21)],
            },
        ]
    )
    monkeypatch.setattr(candidate_audit.parser, "parse", lambda _: processes)
    report = audit_candidate_collision_set(DATASET_PATH, "NH3")
    assert report["eedf_input_coverage"] is True
    assert report["explicit_neutral_fragment_labels"] == ["NH2+H"]
    assert report["passes_preliminary_structure_check"] is False
    assert any("mixes ELASTIC" in issue for issue in report["structural_issues"])
