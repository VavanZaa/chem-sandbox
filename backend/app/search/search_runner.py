from __future__ import annotations

import json
import math
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
from ase import Atoms
from ase.io import write

from app.chemistry.bond_detector import BondDetector
from app.chemistry.molecule_builder import MoleculeBuilder
from app.chemistry.reaction_detector import ReactionDetector
from app.chemistry.species_detector import SpeciesDetector
from app.experiment.compiler import ExperimentCompiler
from app.experiment.schema import ExperimentDefinition
from app.experiment.stoichiometry import theoretical_yield
from app.graph.graph_statistics import build_statistics
from app.graph.reaction_graph import ReactionGraphBuilder
from app.models.base_calculator import CalculationResult, InteratomicCalculator
from .candidate_generator import CandidateGenerator, CandidateProposal
from .frontier import FrontierLimits, explore

ProgressCallback = Callable[[str, int, str], None]


class ScientificValidationError(ValueError):
    """A finite calculation outside the supported species/connectivity domain."""


class SearchRunner:
    def __init__(self, calculator: InteratomicCalculator, *, output_dir: Path,
                 progress: ProgressCallback | None = None) -> None:
        self.calculator = calculator
        self.output_dir = output_dir
        self.progress = progress or (lambda stage, percent, message: None)
        self.builder = MoleculeBuilder()
        self.compiler = ExperimentCompiler()
        self.bonds = BondDetector()
        self.reactions = ReactionDetector(self.bonds)
        self.species = SpeciesDetector(self.bonds)
        self.candidates = CandidateGenerator(self.bonds)
        self.network = ReactionGraphBuilder()
        self.valid_states = 0
        self.failed_optimizations = 0
        self.rejected_candidate_states = 0
        self.errors: list[dict[str, Any]] = []
        self.last_failure = "model_failure"

    def _log(self, event: dict[str, Any]) -> None:
        with (self.output_dir / "events.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(json.dumps({"timestamp": time.time(), **event}, ensure_ascii=False) + "\n")

    def _store_structure(self, atoms: Atoms, run_index: int, name: str) -> str:
        structures = self.output_dir / "structures"
        structures.mkdir(parents=True, exist_ok=True)
        path = structures / f"run_{run_index:03d}_{name}.xyz"
        write(path, atoms)
        return str(path.relative_to(self.output_dir.parent.parent)).replace("\\", "/")

    def _evaluate(self, atoms: Atoms, *, run_index: int, name: str, steps: int,
                  fmax: float, method: str) -> tuple[CalculationResult, str] | None:
        reference = None
        try:
            result = self.calculator.optimize(atoms, fmax=fmax, steps=steps)
            if not math.isfinite(result.energy_ev) or not np.isfinite(result.forces_ev_per_angstrom).all():
                raise ValueError("model returned non-finite energy or forces")
            if not np.isfinite(result.atoms.positions).all():
                raise ValueError("model returned non-finite coordinates")
            if result.atoms.get_chemical_symbols() != atoms.get_chemical_symbols():
                raise ValueError("model changed atom inventory/order")
            reference = self._store_structure(result.atoms, run_index, name)
            minimum_distance = self.bonds.minimum_pair_distance(result.atoms)
            if minimum_distance < self.bonds.minimum_distance:
                raise ScientificValidationError(f"atom clash: {minimum_distance:.3f} Å")
            unresolved = [item for item in self.species.identify(result.atoms) if not item.resolved]
            if unresolved and atoms.info.get("chemistry_mode") != "general":
                reasons = ", ".join(f"{item.formula}: {item.name}" for item in unresolved)
                raise ScientificValidationError(f"unsupported connectivity ({reasons})")
            if steps > 0:
                self.valid_states += 1
            result.metadata["structure_reference"] = reference
            self._log({"kind": "calculation", "run": run_index, "name": name, "method": method,
                       "energy_ev": result.energy_ev, "max_force": result.max_force,
                       "optimizer_converged": result.converged, "optimization_steps": result.steps,
                       "minimum_distance_angstrom": minimum_distance, "structure_reference": reference})
            return result, reference
        except Exception as error:
            rejected = isinstance(error, ScientificValidationError)
            self.last_failure = "unsupported_connectivity" if rejected else "model_failure"
            self.rejected_candidate_states += int(rejected)
            self.failed_optimizations += int(not rejected)
            failure = {"kind": "scientific_rejection" if rejected else "calculation_error",
                       "run": run_index, "name": name, "method": method,
                       "error": f"{type(error).__name__}: {error}", "structure_reference": reference}
            self.errors.append(failure)
            self._log(failure)
            return None

    def _add_state(self, payload: tuple[CalculationResult, str], run_index: int) -> str:
        state, reference = payload
        return self.network.add_state(state.atoms, state.energy_ev, run_index,
                                      converged=state.converged, max_force=state.max_force,
                                      structure_reference=reference)

    def _record_transition(self, before: CalculationResult, before_id: str,
                           after: CalculationResult, after_id: str,
                           proposal: CandidateProposal, run_index: int) -> None:
        event = self.reactions.compare(before.atoms, after.atoms)
        self.network.add_transition(before_id, after_id, event,
                                    delta_energy_ev=after.energy_ev - before.energy_ev,
                                    method=proposal.method, run_index=run_index,
                                    source_reference=before.metadata.get("structure_reference"),
                                    target_reference=after.metadata.get("structure_reference"))
        self._log({"kind": "candidate_transition", "run": run_index, "method": proposal.method,
                   "source": before_id, "target": after_id, "description": proposal.description,
                   "assumptions": proposal.assumptions, "connectivity_changed": event.connectivity_changed,
                   "before": event.before_formula, "after": event.after_formula,
                   "bonds_created": event.changes.created, "bonds_removed": event.changes.removed,
                   "delta_energy_ev": after.energy_ev - before.energy_ev})

    def run(self, experiment: ExperimentDefinition) -> dict[str, Any]:
        started = time.perf_counter()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "experiment.json").write_text(experiment.model_dump_json(indent=2), encoding="utf-8")
        plan = self.compiler.compile(experiment)
        search = experiment.search
        built = self.builder.build(experiment)
        write(self.output_dir / "starting_structure.xyz", built.atoms)
        rng = np.random.default_rng(search.seed)
        self.progress("preparing", 8, "Preparing molecular structures…")
        self._log({"kind": "run_start", "model": self.calculator.metadata, "seed": search.seed,
                   "builder_notes": built.notes, "compiler_assumptions": plan.assumptions})
        outcomes = []

        for run_index in range(plan.runs):
            percent = 12 + int(75 * run_index / plan.runs)
            self.progress("searching", percent, f"Exploring seed {run_index + 1}/{plan.runs} · branching frontier")
            start_atoms = self.candidates.perturb(built.atoms, rng, strength=plan.perturbation_strength)
            raw = self._evaluate(start_atoms, run_index=run_index, name="starting_single_point",
                                 steps=0, fmax=plan.fmax, method="single_point")
            if raw is None:
                outcome = {"run": run_index, "status": self.last_failure, "evaluations": 0,
                           "discovered_states": 0, "expanded_states": 0, "expansions": [], "attempt_counts": {},
                           "endpoints": [{"state_id": None, "depth": 0, "reason": self.last_failure}],
                           "initialization_evaluations": 1}
            else:
                raw_id = self._add_state(raw, run_index)
                self.network.mark_starting(raw_id, run_index)
                initial = self._evaluate(start_atoms, run_index=run_index, name="initial_relaxation",
                                         steps=plan.optimization_steps, fmax=plan.fmax, method="initial_relaxation")
                if initial is None:
                    outcome = {"run": run_index, "status": self.last_failure, "evaluations": 0,
                               "discovered_states": 1, "expanded_states": 0, "expansions": [], "attempt_counts": {},
                               "endpoints": [{"state_id": raw_id, "depth": 0, "reason": self.last_failure}],
                               "initialization_evaluations": 2}
                else:
                    initial_id = self._add_state(initial, run_index)
                    self._record_transition(raw[0], raw_id, initial[0], initial_id,
                                            CandidateProposal(start_atoms, "initial_relaxation", "Relaxed initial geometry", []), run_index)

                    def evaluate(source_id, source, proposal, attempt):
                        self.progress("searching", percent, f"Seed {run_index + 1}/{plan.runs} · candidate {attempt}/{search.max_evaluations_per_run}")
                        candidate = self._evaluate(proposal.atoms, run_index=run_index,
                                                   name=f"candidate_{attempt:04d}", steps=plan.candidate_optimization_steps,
                                                   fmax=plan.fmax, method=proposal.method)
                        if candidate is None:
                            return None, None, self.last_failure
                        target_id = self._add_state(candidate, run_index)
                        self._record_transition(source, source_id, candidate[0], target_id, proposal, run_index)
                        return target_id, candidate[0], None

                    outcome = explore(initial_id, initial[0],
                                      propose=lambda state: self.candidates.generate(state.atoms, rng, plan.candidates_per_state),
                                      evaluate=evaluate, limits=FrontierLimits(search.max_depth,
                                      search.max_evaluations_per_run, search.max_states_per_run, search.frontier_width))
                    outcome.update(run=run_index, initialization_evaluations=2)
            for endpoint in outcome["endpoints"]:
                if endpoint["state_id"]:
                    self.network.mark_terminal(endpoint["state_id"], run_index, endpoint["reason"])
            outcomes.append(outcome)
            self._log({"kind": "run_outcome", **outcome})

        self.progress("building", 93, "Accounting for branches and calculating stoichiometric ceilings…")
        graph = self.network.serialize(plan.runs)
        statistics = build_statistics(graph, total_runs=plan.runs, valid_states=self.valid_states,
                                      failed_optimizations=self.failed_optimizations,
                                      rejected_candidate_states=self.rejected_candidate_states, outcomes=outcomes)
        result = {"schema_version": "0.4", "experiment_id": experiment.experiment_id,
                  "status": "completed" if graph["nodes"] else "model_failure", "model": self.calculator.metadata,
                  "graph": graph, "statistics": statistics, "run_outcomes": outcomes,
                  "search_settings": search.model_dump(), "experiment": experiment.model_dump(mode="json"),
                  "stoichiometry": theoretical_yield(experiment) if experiment.mode == "bicarbonate" else None,
                  "scientific_notes": [*built.notes, *plan.assumptions], "errors": self.errors,
                  "runtime_seconds": time.perf_counter() - started, "barriers_calculated": False}
        (self.output_dir / "result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
        self.progress("complete", 100, "Bounded reaction network ready")
        return result
