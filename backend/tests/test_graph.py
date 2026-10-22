from ase import Atoms

from app.chemistry.reaction_detector import ReactionDetector
from app.graph.reaction_graph import ReactionGraphBuilder


def test_duplicate_states_and_frequency_serialization() -> None:
    reactant = Atoms("OH", positions=[[0, 0, 0], [0.97, 0, 0]])
    product = Atoms("OH", positions=[[0, 0, 0], [3.0, 0, 0]])
    graph = ReactionGraphBuilder()
    source = graph.add_state(reactant, -2.0, 0, converged=True, max_force=0.0, structure_reference=None)
    duplicate = graph.add_state(reactant, -2.1, 1, converged=True, max_force=0.0, structure_reference=None)
    target = graph.add_state(product, -1.0, 0, converged=True, max_force=0.0, structure_reference=None)
    assert source == duplicate
    event = ReactionDetector().compare(reactant, product)
    graph.add_transition(source, target, event, delta_energy_ev=1.0, method="candidate_search", run_index=0)
    graph.add_transition(source, target, event, delta_energy_ev=1.2, method="candidate_search", run_index=1)
    graph.mark_starting(source, 0)
    graph.mark_starting(source, 1)
    graph.mark_terminal(target, 0)
    graph.mark_terminal(target, 1)
    payload = graph.serialize(total_runs=2)
    assert len(payload["nodes"]) == 2
    assert payload["edges"][0]["transition_count"] == 2
    assert payload["edges"][0]["run_count"] == 2
    assert payload["edges"][0]["search_frequency"] == 1.0
    assert payload["edges"][0]["run_frequency"] == 1.0
    assert payload["edges"][0]["barrier_ev"] is None
    by_id = {node["id"]: node for node in payload["nodes"]}
    assert by_id[source]["starting_runs"] == 2
    assert by_id[source]["role"] == "reactant"
    assert by_id[target]["terminal_runs"] == 2
    assert by_id[target]["terminal_run_frequency"] == 1.0
    assert by_id[target]["role"] == "endpoint"
    assert by_id[target]["endpoint_reasons"] == {"no_proposal": 2}
