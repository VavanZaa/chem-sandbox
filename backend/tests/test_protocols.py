import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError

from app.main import app
from app.protocols.runner import simulate_protocol
from app.protocols.schema import LabProtocolRequest
from app.protocols.sbml_adapter import inspect_sbml


def protocol_payload():
    return {
        "id": "staged_hexose_cycle_execution_test",
        "species": [
            {"id": "glucose", "name": "Glucose stereoisomer", "smiles": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"},
            {"id": "fructose", "name": "Open-chain fructose", "smiles": "OCC(=O)C(O)C(O)C(O)CO"},
            {"id": "mannitol", "name": "Mannitol", "smiles": "OCC(O)C(O)C(O)C(O)CO"},
            {"id": "sorbose", "name": "Sorbose identity proxy", "smiles": "OCC(=O)C(O)C(O)C(O)CO"},
            {"id": "h2", "name": "Hydrogen", "smiles": "[H][H]"},
        ],
        "reactions": [
            {"id": "isomerize_1", "name": "Declared glucose to fructose isomerization",
             "stoichiometry": {"glucose": -1, "fructose": 1},
             "rate_law": {"kind": "mass_action", "rate_constant": {"value": .08, "unit": "1/s"},
                          "source": "hypothetical execution fixture", "basis": "hypothetical"}},
            {"id": "hydrogenate", "name": "Declared fructose hydrogenation",
             "stoichiometry": {"fructose": -1, "h2": -1, "mannitol": 1},
             "rate_law": {"kind": "mass_action", "rate_constant": {"value": 2, "unit": "L/mol/s"},
                          "source": "hypothetical execution fixture", "basis": "hypothetical"}},
            {"id": "dehydrogenate", "name": "Declared mannitol dehydrogenation",
             "stoichiometry": {"mannitol": -1, "sorbose": 1, "h2": 1},
             "rate_law": {"kind": "mass_action", "rate_constant": {"value": .04, "unit": "1/s"},
                          "source": "hypothetical execution fixture", "basis": "hypothetical"}},
            {"id": "isomerize_2", "name": "Declared sorbose to glucose isomerization",
             "stoichiometry": {"sorbose": -1, "glucose": 1},
             "rate_law": {"kind": "mass_action", "rate_constant": {"value": .02, "unit": "1/s"},
                          "source": "hypothetical execution fixture", "basis": "hypothetical"}},
        ],
        "initial_amounts": {"glucose": {"value": .01, "unit": "mol"}},
        "initial_volume": {"value": .1, "unit": "L"},
        "initial_temperature": {"value": 298.15, "unit": "K"},
        "stages": [
            {"id": "convert_1", "name": "First conversion", "operation": "reaction",
             "duration": {"value": 8, "unit": "s"}, "active_reactions": ["isomerize_1"]},
            {"id": "mid_addition", "name": "Add hydrogen after stage one", "operation": "addition",
             "additions": [{"species_id": "h2", "amount": {"value": .02, "unit": "mol"}}]},
            {"id": "convert_2", "name": "Hydrogenation stage", "operation": "reaction",
             "duration": {"value": 12, "unit": "s"}, "active_reactions": ["hydrogenate"]},
            {"id": "sample_1", "name": "Remove analytical sample", "operation": "sample",
             "withdraw_fraction": .05},
            {"id": "convert_3", "name": "Dehydrogenation stage", "operation": "reaction",
             "duration": {"value": 6, "unit": "s"}, "active_reactions": ["dehydrogenate"]},
            {"id": "hold_1", "name": "Timed hold", "operation": "hold",
             "duration": {"value": 3, "unit": "s"}},
            {"id": "convert_4", "name": "Close declared cycle", "operation": "reaction",
             "duration": {"value": 5, "unit": "s"}, "active_reactions": ["isomerize_2"]},
        ],
        "points_per_reaction_stage": 11,
        "include_sbml": True,
        "cross_check_roadrunner": True,
        "acknowledge_evidence_limits": True,
    }


def test_long_staged_protocol_mid_chain_addition_and_cross_check():
    result = simulate_protocol(LabProtocolRequest.model_validate(protocol_payload()))
    assert result["status"] == "conditional_protocol_simulation"
    assert len(result["stages"]) == 7
    assert result["stages"][1]["additions_mol"]["h2"] == pytest.approx(.02)
    assert result["stages"][3]["removed_mol"]["glucose"] > 0
    assert result["material_ledger"]["closed"]
    assert max(abs(value) for value in result["material_ledger"]["residual"].values()) < 1e-9
    assert result["scientific_boundary"]["recipe_used_for_blind_discovery"] is False
    assert result["roadrunner_cross_check"]["status"] == "passed"
    assert result["sbml"]["validation"]["status"] == "valid"
    assert result["sbml"]["document"].startswith("<?xml")


def test_protocol_species_stereochemistry_is_preserved():
    result = simulate_protocol(LabProtocolRequest.model_validate(protocol_payload()))
    glucose = next(item for item in result["species_identity"] if item["id"] == "glucose")
    assert glucose["stereo_status"] == "specified"
    assert len(glucose["stereochemistry"]) >= 4


def test_protocol_rejects_unbalanced_reaction():
    payload = protocol_payload()
    payload["reactions"][0]["stoichiometry"] = {"glucose": -1, "mannitol": 1}
    with pytest.raises(ValueError, match="does not conserve"):
        simulate_protocol(LabProtocolRequest.model_validate(payload))


def test_protocol_rejects_unknown_mid_stage_addition():
    payload = protocol_payload()
    payload["stages"][1]["additions"][0]["species_id"] = "unknown"
    with pytest.raises(ValidationError, match="unknown species"):
        LabProtocolRequest.model_validate(payload)


def test_protocol_api_and_sbml_inspection():
    response = TestClient(app).post("/api/protocols/simulate", json=protocol_payload())
    assert response.status_code == 200, response.text
    result = response.json()
    inspected = inspect_sbml(result["sbml"]["document"])
    assert inspected["validation"]["status"] == "valid"
    assert len(inspected["model"]["species"]) == 5
    assert len(inspected["model"]["reactions"]) == 4
    api_inspection = TestClient(app).post("/api/protocols/sbml/inspect", json={"sbml": result["sbml"]["document"]})
    assert api_inspection.status_code == 200
    assert api_inspection.json()["sha256"] == result["sbml"]["sha256"]


def test_protocol_stage_operations_are_strict():
    payload = protocol_payload()
    payload["stages"][0]["duration"] = {"value": 0, "unit": "s"}
    with pytest.raises(ValidationError, match="positive duration"):
        LabProtocolRequest.model_validate(payload)


def test_mid_protocol_catalyst_addition_enables_michaelis_menten_stage():
    payload = protocol_payload()
    payload["species"].append({"id": "enzyme", "name": "Declared enzyme pool",
                               "composition": {"C": 10, "H": 20, "N": 2, "O": 5}, "charge": 0,
                               "identity_note": "Composition-only execution fixture, not a real enzyme identity"})
    payload["reactions"] = [{
        "id": "enzyme_isomerization", "name": "Declared enzyme-catalyzed isomerization",
        "stoichiometry": {"glucose": -1, "fructose": 1},
        "rate_law": {"kind": "michaelis_menten", "substrate_id": "glucose", "catalyst_id": "enzyme",
                     "kcat": {"value": 2, "unit": "1/s"}, "km": {"value": .02, "unit": "mol/L"},
                     "source": "hypothetical execution fixture", "basis": "hypothetical"},
    }]
    payload["stages"] = [
        {"id": "add_enzyme", "name": "Add catalyst after initialization", "operation": "addition",
         "additions": [{"species_id": "enzyme", "amount": {"value": .0002, "unit": "mol"}}]},
        {"id": "enzyme_stage", "name": "Enzyme stage", "operation": "reaction",
         "duration": {"value": 5, "unit": "s"}, "active_reactions": ["enzyme_isomerization"]},
    ]
    result = simulate_protocol(LabProtocolRequest.model_validate(payload))
    assert result["stages"][0]["additions_mol"]["enzyme"] == pytest.approx(.0002)
    assert result["final_amounts_mol"]["fructose"] > 0
    assert result["final_amounts_mol"]["enzyme"] == pytest.approx(.0002)
    assert result["roadrunner_cross_check"]["status"] == "passed"


def test_eyring_stage_and_temperature_change_cross_check():
    payload = protocol_payload()
    payload["reactions"] = [{
        "id": "eyring_isomerization", "name": "Declared Eyring stage",
        "stoichiometry": {"glucose": -1, "fructose": 1},
        "rate_law": {"kind": "eyring_mass_action",
                     "activation_free_energy": {"value": 82, "unit": "kJ/mol"},
                     "reference_temperature": {"value": 298.15, "unit": "K"},
                     "activation_entropy": {"value": 0, "unit": "J/mol/K"},
                     "source": "hypothetical execution fixture", "basis": "hypothetical"},
    }]
    payload["stages"] = [{
        "id": "heated_stage", "name": "Change temperature then react", "operation": "reaction",
        "temperature": {"value": 320, "unit": "K"}, "duration": {"value": 2, "unit": "s"},
        "active_reactions": ["eyring_isomerization"],
    }]
    result = simulate_protocol(LabProtocolRequest.model_validate(payload))
    assert result["stages"][0]["temperature_K"] == pytest.approx(320)
    assert result["final_amounts_mol"]["fructose"] > 0
    assert result["roadrunner_cross_check"]["status"] == "passed"


def test_mass_action_rate_constant_units_follow_molecularity():
    payload = protocol_payload()
    payload["reactions"][1]["rate_law"]["rate_constant"] = {"value": 2, "unit": "1/s"}
    with pytest.raises(ValueError, match="not compatible"):
        simulate_protocol(LabProtocolRequest.model_validate(payload))


def test_twelve_step_cyclic_network_with_mid_cycle_feed():
    species = [
        {"id": f"pool_{index}", "name": f"Abstract conserved pool {index}",
         "composition": {"C": 4, "H": 6, "O": 5}, "charge": 0,
         "identity_note": "Topology stress fixture; not assigned to a real metabolite"}
        for index in range(12)
    ]
    reactions = []
    stages = []
    for index in range(12):
        next_index = (index + 1) % 12
        reaction_id = f"step_{index + 1}"
        reactions.append({
            "id": reaction_id, "name": f"Declared cyclic step {index + 1}",
            "stoichiometry": {f"pool_{index}": -1, f"pool_{next_index}": 1},
            "rate_law": {"kind": "mass_action", "rate_constant": {"value": .15, "unit": "1/s"},
                         "source": "hypothetical twelve-step topology stress fixture", "basis": "hypothetical"},
        })
        stages.append({
            "id": f"stage_{index + 1}", "name": f"Sequential stage {index + 1}", "operation": "reaction",
            "duration": {"value": 2, "unit": "s"}, "active_reactions": [reaction_id],
        })
        if index == 5:
            stages.append({
                "id": "mid_cycle_feed", "name": "Add pool material after step six", "operation": "addition",
                "additions": [{"species_id": "pool_6", "amount": {"value": .002, "unit": "mol"}}],
            })
    request = LabProtocolRequest.model_validate({
        "id": "twelve_step_cycle_stress", "species": species, "reactions": reactions,
        "initial_amounts": {"pool_0": {"value": .01, "unit": "mol"}},
        "initial_volume": {"value": .1, "unit": "L"},
        "initial_temperature": {"value": 310, "unit": "K"}, "stages": stages,
        "points_per_reaction_stage": 5, "cross_check_roadrunner": True,
        "acknowledge_evidence_limits": True,
    })
    result = simulate_protocol(request)
    assert len(result["stages"]) == 13
    assert result["material_ledger"]["added_mol"]["pool_6"] == pytest.approx(.002)
    assert result["material_ledger"]["closed"]
    assert result["roadrunner_cross_check"]["status"] == "passed"
