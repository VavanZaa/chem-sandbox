"""Conditional source-driven HCN tests; none are experimental-yield tests."""
from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.prebiotic.hcn_ceiling import hcn_atom_ceiling, hcn_atom_ceiling_from_atoms
from app.prebiotic.reduced_hcn import ReducedHCNRequest, audit_reduced_hcn, simulate_reduced_hcn
from app.prebiotic.schema import PrebioticApparatusRequest
from scripts.run_hcn_formation_gate import make_formation_scenario


def evidence(value: float, unit: str, *, basis: str = "hypothetical") -> dict:
    return {
        "value": value, "unit": unit, "source": "synthetic regression fixture",
        "method": "declared source-driven control", "basis": basis,
    }


def request_data() -> dict:
    return {
        "initial_gas_mol": {"CH4": 0.01, "N2": 0.02, "H2": 0.02, "H2O": 0.01},
        "gas_volume": evidence(1, "L"),
        "water_volume": evidence(0.05, "L"),
        "temperature": evidence(298.15, "K"),
        "duration": evidence(60, "s"),
        "fixed_pH": {
            "value": 7, "source": "synthetic regression fixture",
            "method": "declared fixed pH", "basis": "hypothetical",
        },
        "hcn_henry_hcp": evidence(0.12, "mol/(m^3*Pa)"),
        "hcn_liquid_side_kla": evidence(1e-3, "1/s"),
        "deposited_energy": evidence(2000, "J"),
        "sources": [
            {
                "id": "methane_cleavage", "parent": "CH4", "products": {"CH3": 1, "H": 1},
                "parent_conversion_hazard": evidence(1e-3, "1/s"),
                "minimum_event_energy": evidence(400_000, "J/mol"),
                "product_state_source": "synthetic ground-state assignment only",
                "product_state_basis": "hypothetical",
            },
            {
                "id": "nitrogen_cleavage", "parent": "N2", "products": {"N": 2},
                "parent_conversion_hazard": evidence(1e-3, "1/s"),
                "minimum_event_energy": evidence(900_000, "J/mol"),
                "product_state_source": "synthetic N(4S) assignment, not the N2 total collision table",
                "product_state_basis": "hypothetical",
            },
        ],
        "sample_count": 5,
        "hcn_assay_used_to_set_inputs": False,
    }


def test_source_driven_run_is_conditional_and_conserves_atoms() -> None:
    result = simulate_reduced_hcn(ReducedHCNRequest.model_validate(request_data()))
    assert result["status"] == "conditional_reduced_mechanistic_run"
    assert result["atom_ledger_closed"] is True
    assert result["experimental_yield_validated"] is False
    assert "source:nitrogen_cleavage:product_state" in result["hypothetical_input_fields"]
    assert result["numerical_convergence_assessed_in_this_run"] is False
    assert result["hcn_assay_used_to_set_inputs"] is False
    assert result["atom_ceiling"]["cyanide_family_total_ceiling_mol"] == pytest.approx(.01)
    assert result["final_total_hcn_plus_cn_mol"] <= .01
    assert result["minimum_fragment_energy_within_deposit"] is True
    assert result["final_aqueous_hcn_plus_cn_mol"] > 0
    assert set(result["source_extents_mol"]) == {"methane_cleavage", "nitrogen_cleavage"}


def test_no_carbon_or_no_nitrogen_has_zero_hcn_and_zero_ceiling() -> None:
    for absent in ("CH4", "N2"):
        payload = request_data()
        payload["initial_gas_mol"].pop(absent)
        result = simulate_reduced_hcn(ReducedHCNRequest.model_validate(payload))
        assert result["atom_ceiling"]["cyanide_family_total_ceiling_mol"] == 0
        assert result["final_total_hcn_plus_cn_mol"] < 1e-12
        assert result["atom_ledger_closed"] is True


def test_source_cannot_inject_hcn_or_unbalanced_fragments() -> None:
    payload = request_data()
    payload["sources"][0]["products"] = {"HCN": 1, "H": 1}
    with pytest.raises(ValueError, match="C-N product"):
        simulate_reduced_hcn(ReducedHCNRequest.model_validate(payload))
    payload = request_data()
    payload["sources"][0]["products"] = {"CH3": 1, "H": 2}
    with pytest.raises(ValueError, match="conserve atoms"):
        simulate_reduced_hcn(ReducedHCNRequest.model_validate(payload))


def test_source_energy_cannot_exceed_declared_deposit() -> None:
    payload = request_data()
    payload["deposited_energy"] = evidence(1e-9, "J")
    with pytest.raises(ValueError, match="exceeds deposited"):
        simulate_reduced_hcn(ReducedHCNRequest.model_validate(payload))


def test_hcn_assay_cannot_be_declared_as_source_input() -> None:
    payload = request_data()
    payload["hcn_assay_used_to_set_inputs"] = True
    with pytest.raises(ValueError):
        ReducedHCNRequest.model_validate(payload)


def test_step_refinement_and_zero_source_control() -> None:
    payload = request_data()
    payload["maximum_solver_step_s"] = 5
    coarse = simulate_reduced_hcn(ReducedHCNRequest.model_validate(payload))
    payload["maximum_solver_step_s"] = 2.5
    fine = simulate_reduced_hcn(ReducedHCNRequest.model_validate(payload))
    assert fine["final_total_hcn_plus_cn_mol"] == pytest.approx(
        coarse["final_total_hcn_plus_cn_mol"], rel=0.01
    )
    payload["sources"] = []
    off = simulate_reduced_hcn(ReducedHCNRequest.model_validate(payload))
    assert off["final_total_hcn_plus_cn_mol"] < 1e-12
    assert off["active_effective_source_count"] == 0


def test_audit_refines_step_and_runs_assay_free_controls() -> None:
    audit = audit_reduced_hcn(ReducedHCNRequest.model_validate(request_data()))
    assert audit["numerical_checks_passed"] is True
    assert audit["solver_step_consistent"] is True
    assert audit["experimental_yield_validated"] is False
    assert audit["hcn_assay_used_to_set_inputs"] is False
    assert audit["refined_maximum_solver_step_s"] < audit["coarse_maximum_solver_step_s"]
    assert audit["source_off_control"]["below_reporting_floor"] is True
    assert all(
        item["status"] == "below_reporting_floor"
        for item in audit["atom_removal_controls"].values()
    )
    assert audit["refined_result"]["atom_ledger_closed"] is True


def test_audit_does_not_call_unrunnable_atom_controls_a_pass() -> None:
    payload = request_data()
    payload["initial_gas_mol"] = {"HCN": 1e-6}
    payload["sources"] = []
    audit = audit_reduced_hcn(ReducedHCNRequest.model_validate(payload))
    assert audit["audit_outcome"] == "incomplete"
    assert audit["numerical_checks_passed"] is False
    assert audit["negative_controls_complete"] is False


def test_seeded_hcn_contact_moves_material_without_forming_it() -> None:
    payload = request_data()
    payload["initial_gas_mol"] = {"N2": .02, "HCN": 1e-6}
    payload["sources"] = []
    captured = simulate_reduced_hcn(ReducedHCNRequest.model_validate(payload))
    assert captured["final_aqueous_hcn_plus_cn_mol"] > 0
    assert captured["final_total_hcn_plus_cn_mol"] == pytest.approx(1e-6, rel=1e-5)
    assert captured["atom_ceiling"]["maximum_net_new_cyanide_family_mol"] == pytest.approx(0)
    payload["hcn_liquid_side_kla"] = evidence(0, "1/s")
    dry = simulate_reduced_hcn(ReducedHCNRequest.model_validate(payload))
    assert dry["final_aqueous_hcn_plus_cn_mol"] == pytest.approx(0)


def test_atom_ceiling_separates_cn_from_molecular_hcn_and_additions() -> None:
    no_hydrogen = hcn_atom_ceiling_from_atoms({"C": 2, "N": 1, "H": 0})
    assert no_hydrogen["cyanide_family_total_ceiling_mol"] == 1
    assert no_hydrogen["molecular_hcn_total_ceiling_mol"] == 0
    base = make_formation_scenario()
    before = hcn_atom_ceiling(PrebioticApparatusRequest.model_validate(base))
    base["interventions"] = [{
        "cycle": 1, "additions": [{"species_id": "CH4", "compartment": "gas_loop",
                                "amount": {"value": 0.01, "unit": "mol"}}],
    }]
    after = hcn_atom_ceiling(PrebioticApparatusRequest.model_validate(base))
    assert after["available_atom_mol"]["C"] == pytest.approx(before["available_atom_mol"]["C"] + .01)
    assert after["cyanide_family_total_ceiling_mol"] >= before["cyanide_family_total_ceiling_mol"]
    assert after["predicted_yield"] is False


def test_reduced_and_ceiling_api() -> None:
    client = TestClient(app)
    reduced = client.post("/api/prebiotic/hcn/reduced-mechanistic", json=request_data())
    assert reduced.status_code == 200
    assert reduced.json()["experimental_yield_validated"] is False
    audit = client.post("/api/prebiotic/hcn/reduced-mechanistic/audit", json=request_data())
    assert audit.status_code == 200
    assert audit.json()["numerical_checks_passed"] is True
    ceiling = client.post("/api/prebiotic/hcn/atom-ceiling", json=make_formation_scenario())
    assert ceiling.status_code == 200
    assert ceiling.json()["status"] == "element_only_upper_bound"
