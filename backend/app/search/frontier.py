"""Bounded breadth-first exploration, independent of chemistry and energy ranking.

Every discovered connectivity is expanded once per seed, unless explicitly
stopped by a budget. Back edges are retained by the evaluator but not requeued.
Search leaves are computational endpoints, never claims of chemical stability.
"""
from collections import Counter, deque
from dataclasses import dataclass
from typing import Any, Callable


@dataclass(frozen=True)
class FrontierLimits:
    depth: int
    evaluations: int
    states: int
    width: int


def explore(
    initial_id: str,
    initial: Any,
    *,
    propose: Callable[[Any], list],
    evaluate: Callable[[str, Any, Any, int], tuple[str | None, Any, str | None]],
    limits: FrontierLimits,
) -> dict:
    queue = deque([(initial_id, initial, 0)])
    seen = {initial_id}
    links: dict[str, set[str]] = {}
    admitted_by_depth = Counter({0: 1})
    endpoints, expansions = [], []
    attempts = Counter()
    evaluations = 0

    def stop(state_id: str, depth: int, reason: str) -> None:
        endpoints.append({"state_id": state_id, "depth": depth, "reason": reason})

    def closes_cycle(source: str, target: str) -> bool:
        pending, visited = [target], set()
        while pending:
            node = pending.pop()
            if node == source:
                return True
            if node not in visited:
                visited.add(node)
                pending.extend(links.get(node, set()) - visited)
        return False

    while queue:
        state_id, payload, depth = queue.popleft()
        if depth >= limits.depth or evaluations >= limits.evaluations:
            stop(state_id, depth, "depth_budget" if depth >= limits.depth else "evaluation_budget")
            continue
        proposals = propose(payload)
        counts = Counter()
        new_children = 0
        for proposal in proposals:
            if evaluations >= limits.evaluations:
                counts["not_evaluated_budget"] += len(proposals) - sum(counts.values())
                break
            evaluations += 1
            target, candidate, failure = evaluate(state_id, payload, proposal, evaluations)
            if target is None:
                counts[failure or "model_failure"] += 1
                continue
            if target == state_id:
                counts["proposal_reverted"] += 1
                continue
            cyclic = closes_cycle(state_id, target)
            links.setdefault(state_id, set()).add(target)
            if target in seen:
                counts["cycle_revisit" if cyclic else "branch_merge"] += 1
                continue
            seen.add(target)
            counts["new_connectivity"] += 1
            new_children += 1
            if len(seen) > limits.states:
                stop(target, depth + 1, "state_budget")
            elif admitted_by_depth[depth + 1] >= limits.width:
                stop(target, depth + 1, "frontier_budget")
            else:
                admitted_by_depth[depth + 1] += 1
                queue.append((target, candidate, depth + 1))
        attempts.update(counts)
        expansions.append({"state_id": state_id, "depth": depth, "proposals": len(proposals), "counts": dict(counts)})
        if counts["not_evaluated_budget"]:
            stop(state_id, depth, "evaluation_budget")
        elif not new_children:
            if not proposals:
                reason = "no_proposal"
            elif counts["cycle_revisit"]:
                reason = "cycle_closed"
            elif counts["branch_merge"]:
                reason = "branch_merged"
            elif counts["proposal_reverted"] == len(proposals):
                reason = "proposal_reverted"
            elif counts["model_failure"] == len(proposals):
                reason = "model_failure"
            elif counts["unsupported_connectivity"] == len(proposals):
                reason = "unsupported_connectivity"
            else:
                reason = "no_valid_exit"
            stop(state_id, depth, reason)
    return {
        "status": "budget_exhausted" if any(item["reason"].endswith("budget") for item in endpoints) else "explored",
        "evaluations": evaluations,
        "discovered_states": len(seen),
        "expanded_states": len(expansions),
        "endpoints": endpoints,
        "expansions": expansions,
        "attempt_counts": dict(attempts),
    }
