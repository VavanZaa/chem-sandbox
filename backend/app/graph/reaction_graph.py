from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
from typing import Any

import networkx as nx
from ase import Atoms

from app.chemistry.reaction_detector import ReactionEvent
from app.chemistry.species_detector import SpeciesDetector


@dataclass(slots=True)
class StateRecord:
    id: str
    signature: str
    label: str
    formula: str
    species: list[dict[str, Any]]
    interpretation_quality: str
    energy_ev: float
    relative_energy_ev: float = 0.0
    occurrences: int = 0
    runs_seen: set[int] = field(default_factory=set)
    starting_runs: set[int] = field(default_factory=set)
    terminal_runs: set[int] = field(default_factory=set)
    endpoint_reasons: dict[str, int] = field(default_factory=dict)
    converged: bool = False
    max_force: float | None = None
    structure_reference: str | None = None


class ReactionGraphBuilder:
    def __init__(self) -> None:
        self.graph = nx.DiGraph()
        self.species = SpeciesDetector()
        self._signature_to_id: dict[str, str] = {}
        self._next_id = 1

    def add_state(
        self,
        atoms: Atoms,
        energy_ev: float,
        run_index: int,
        *,
        converged: bool,
        max_force: float | None,
        structure_reference: str | None,
    ) -> str:
        signature = self.species.state_signature(atoms)
        label = self.species.state_label(atoms)
        species = self.species.state_species(atoms)
        state_id = self._signature_to_id.get(signature)
        if state_id is None:
            state_id = f"state_{self._next_id:03d}"
            self._next_id += 1
            record = StateRecord(
                id=state_id,
                signature=signature,
                label=label,
                formula=self.species.state_formula(atoms, unicode=False),
                species=species,
                interpretation_quality=(
                    "resolved" if all(item["resolved"] for item in species) else "uncertain"
                ),
                energy_ev=energy_ev,
                occurrences=1,
                runs_seen={run_index},
                converged=converged,
                max_force=max_force,
                structure_reference=structure_reference,
            )
            self._signature_to_id[signature] = state_id
            self.graph.add_node(state_id, record=record)
        else:
            record = self.graph.nodes[state_id]["record"]
            record.occurrences += 1
            record.runs_seen.add(run_index)
            if energy_ev < record.energy_ev:
                record.energy_ev = energy_ev
                record.converged = converged
                record.max_force = max_force
                record.structure_reference = structure_reference
        return state_id

    def add_transition(
        self,
        source: str,
        target: str,
        event: ReactionEvent,
        *,
        delta_energy_ev: float,
        method: str,
        run_index: int,
        source_reference: str | None = None,
        target_reference: str | None = None,
    ) -> None:
        if source == target or not event.connectivity_changed:
            return
        if self.graph.has_edge(source, target):
            edge = self.graph.edges[source, target]
            edge["transition_count"] += 1
            edge["runs_seen"].add(run_index)
            edge["delta_energy_samples"].append(delta_energy_ev)
        else:
            self.graph.add_edge(
                source,
                target,
                transition_count=1,
                runs_seen={run_index},
                delta_energy_samples=[delta_energy_ev],
                bonds_created=[list(pair) for pair in event.changes.created],
                bonds_removed=[list(pair) for pair in event.changes.removed],
                method=method,
                barrier_ev=None,
                source_reference=source_reference,
                target_reference=target_reference,
            )

    def mark_terminal(self, state_id: str, run_index: int, reason: str = "no_proposal") -> None:
        """Record a branch stop; several states can be endpoints in one seed."""

        if state_id in self.graph:
            self.graph.nodes[state_id]["record"].terminal_runs.add(run_index)
            reasons = self.graph.nodes[state_id]["record"].endpoint_reasons
            reasons[reason] = reasons.get(reason, 0) + 1

    def mark_starting(self, state_id: str, run_index: int) -> None:
        """Record the measured starting connectivity for one search run."""

        if state_id in self.graph:
            self.graph.nodes[state_id]["record"].starting_runs.add(run_index)

    def serialize(self, total_runs: int) -> dict[str, Any]:
        records = [data["record"] for _, data in self.graph.nodes(data=True)]
        if not records:
            return {"nodes": [], "edges": []}
        minimum_energy = min(record.energy_ev for record in records)
        nodes: list[dict[str, Any]] = []
        for state_id in nx.topological_sort(self.graph) if nx.is_directed_acyclic_graph(self.graph) else self.graph.nodes:
            record: StateRecord = self.graph.nodes[state_id]["record"]
            record.relative_energy_ev = record.energy_ev - minimum_energy
            payload = asdict(record)
            payload["comparison_key"] = hashlib.sha256(record.signature.encode("utf-8")).hexdigest()[:20]
            payload["runs_seen"] = len(record.runs_seen)
            payload["run_indices"] = sorted(record.runs_seen)
            payload["endpoint_run_indices"] = sorted(record.terminal_runs)
            payload["run_frequency"] = len(record.runs_seen) / max(total_runs, 1)
            payload["starting_runs"] = len(record.starting_runs)
            payload["terminal_runs"] = len(record.terminal_runs)
            payload["terminal_run_frequency"] = len(record.terminal_runs) / max(total_runs, 1)
            payload.pop("signature", None)
            if record.starting_runs:
                payload["role"] = "reactant"
            elif record.terminal_runs:
                payload["role"] = "endpoint"
            else:
                payload["role"] = "intermediate"
            nodes.append(payload)

        edges: list[dict[str, Any]] = []
        for source, target, data in self.graph.edges(data=True):
            outgoing_total = sum(
                self.graph.edges[source, neighbor]["transition_count"] for neighbor in self.graph.successors(source)
            )
            samples = data["delta_energy_samples"]
            edges.append(
                {
                    "id": f"{source}__{target}",
                    "source": source,
                    "target": target,
                    "comparison_key": (
                        hashlib.sha256(self.graph.nodes[source]["record"].signature.encode("utf-8")).hexdigest()[:20]
                        + "=>" +
                        hashlib.sha256(self.graph.nodes[target]["record"].signature.encode("utf-8")).hexdigest()[:20]
                    ),
                    "transition_count": data["transition_count"],
                    "run_count": len(data["runs_seen"]),
                    "search_frequency": data["transition_count"] / max(outgoing_total, 1),
                    "run_frequency": len(data["runs_seen"]) / max(total_runs, 1),
                    "delta_energy_ev": sum(samples) / len(samples),
                    "energy_assessment": (
                        "favorable"
                        if sum(samples) / len(samples) < -0.05
                        else "unfavorable"
                        if sum(samples) / len(samples) > 0.05
                        else "approximately_neutral"
                    ),
                    "barrier_ev": data["barrier_ev"],
                    "source_reference": data.get("source_reference"),
                    "target_reference": data.get("target_reference"),
                    "method": data["method"],
                    "bonds_created": data["bonds_created"],
                    "bonds_removed": data["bonds_removed"],
                }
            )
        return {"nodes": nodes, "edges": edges}
