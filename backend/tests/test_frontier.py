import pytest
from app.search.frontier import FrontierLimits, explore


def run_graph(graph, **limits):
    transitions = []
    def evaluate(source, payload, target, attempt):
        transitions.append((source, target))
        if target.startswith("!"):
            return None, None, target[1:]
        return target, target, None
    result = explore("A", "A", propose=lambda state: graph.get(state, []), evaluate=evaluate,
                     limits=FrontierLimits(**({"depth": 4, "evaluations": 30, "states": 20, "width": 10} | limits)))
    return result, transitions


def test_all_branches_expanded_and_cycles_retained_without_looping():
    result, transitions = run_graph({"A": ["B", "C"], "B": ["D"], "C": ["E"], "D": ["A"]})
    assert set(transitions) == {("A", "B"), ("A", "C"), ("B", "D"), ("C", "E"), ("D", "A")}
    assert result["expanded_states"] == 5
    assert {(item["state_id"], item["reason"]) for item in result["endpoints"]} == {("D", "cycle_closed"), ("E", "no_proposal")}


def test_merging_branches_are_not_mislabeled_as_cycles():
    result, transitions = run_graph({"A": ["B", "C"], "B": ["D"], "C": ["D"]})
    assert ("C", "D") in transitions
    assert {item["reason"] for item in result["endpoints"]} == {"branch_merged", "no_proposal"}
    assert result["expanded_states"] == 4


@pytest.mark.parametrize("limits, reason", [({"depth": 1}, "depth_budget"), ({"evaluations": 1}, "evaluation_budget"),
                                             ({"states": 1}, "state_budget"), ({"width": 1}, "frontier_budget")])
def test_every_budget_reports_unexplored_branches(limits, reason):
    result, _ = run_graph({"A": ["B", "C"], "B": ["D"]}, **limits)
    assert result["status"] == "budget_exhausted"
    assert reason in {item["reason"] for item in result["endpoints"]}
    assert result["evaluations"] <= limits.get("evaluations", 30)


@pytest.mark.parametrize("targets, reason", [(["A", "A"], "proposal_reverted"), (["!model_failure"], "model_failure"),
                                               (["!unsupported_connectivity"], "unsupported_connectivity"),
                                               (["A", "!unsupported_connectivity"], "no_valid_exit")])
def test_failed_and_reverted_proposals_accounted(targets, reason):
    result, _ = run_graph({"A": targets})
    assert result["endpoints"][0]["reason"] == reason
    assert sum(result["attempt_counts"].values()) == len(targets)
