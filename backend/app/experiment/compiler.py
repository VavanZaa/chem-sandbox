from __future__ import annotations

from dataclasses import dataclass

from .schema import ExperimentDefinition


@dataclass(frozen=True, slots=True)
class SearchPlan:
    runs: int
    optimization_steps: int
    fmax: float
    perturbation_strength: float
    candidates_per_state: int
    candidate_optimization_steps: int
    assumptions: tuple[str, ...]


class ExperimentCompiler:
    """Translate UI controls into bounded search parameters."""

    _effort_steps = {"low": 6, "medium": 10, "high": 10}
    _mixing_scale = {"gentle": 0.75, "standard": 1.0, "vigorous": 1.35}

    def compile(self, experiment: ExperimentDefinition) -> SearchPlan:
        search = experiment.search
        temperature_scale = (experiment.environment.temperature_K / 298.15) ** 0.5
        strength = 0.10 * temperature_scale * self._mixing_scale[experiment.environment.mixing]
        default_steps = self._effort_steps[search.effort]
        requested_steps = search.optimization_steps if search.optimization_steps > 0 else default_steps
        steps = requested_steps
        return SearchPlan(
            runs=search.runs,
            optimization_steps=steps,
            fmax=0.35,
            perturbation_strength=strength,
            candidates_per_state=search.max_candidates_per_state,
            candidate_optimization_steps=search.candidate_optimization_steps,
            assumptions=(
                "Temperature scales proposal perturbations; it is not a molecular-dynamics thermostat.",
                "Pressure does not change the finite-cluster energy model; it is used only in the ideal-gas yield ceiling.",
                "Mixing scales orientation/translation perturbations only.",
                "Relaxation step limits are user-configurable computational budgets, not evidence of stability.",
                "Breadth-first search expands multiple connectivities with cycle deduplication and explicit depth, width, state and evaluation budgets.",
                "Each seed explores a graph, not a physical trajectory. Branch endpoint frequencies overlap and are not yields.",
                ("General mode adds association/dissociation geometric scans; unknown fragments are retained but their identity, charge and bond orders remain unvalidated. Operators are not exhaustive."
                 if experiment.mode == "general" else "Local proton transfer and carbon-oxo water-loss proposals are incomplete; unsupported species are saved but excluded from the interpreted graph."),
            ),
        )
