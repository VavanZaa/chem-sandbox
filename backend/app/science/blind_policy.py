from __future__ import annotations

import hashlib
import json
from typing import Any

from app.experiment.schema import ExperimentDefinition


POLICY_VERSION = "blind-discovery-v1"

# These concepts belong only in the post-run evaluator. Pydantic's strict
# schemas reject them at the API boundary; this list is also exposed so the
# contract is machine-readable and testable.
FORBIDDEN_DISCOVERY_FIELDS = frozenset(
    {
        "answer",
        "expected_pathway",
        "expected_product",
        "known_product",
        "known_yield",
        "literature_product",
        "literature_yield",
        "product_recipe",
        "reaction_recipe",
        "reference_energy",
        "reference_product",
        "reference_yield",
        "target_product",
    }
)


def _canonical_payload(experiment: ExperimentDefinition) -> tuple[dict[str, Any], str]:
    payload = experiment.model_dump(mode="json", exclude_none=True)
    encoded = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return payload, hashlib.sha256(encoded).hexdigest()


def build_blind_manifest(experiment: ExperimentDefinition) -> dict[str, Any]:
    """Create the immutable audit record attached to every discovery result.

    This deliberately records policy and provenance, not benchmark answers.
    The hash lets an evaluator prove which inputs were frozen before references
    became visible.
    """

    payload, digest = _canonical_payload(experiment)
    return {
        "policy_version": POLICY_VERSION,
        "mode": experiment.knowledge_policy.mode,
        "input_sha256": digest,
        "input_fields": sorted(payload),
        "target_specific_knowledge_used": False,
        "reaction_recipe_used": False,
        "reference_access": experiment.knowledge_policy.reference_access,
        "reaction_family_priors": experiment.knowledge_policy.reaction_family_priors,
        "forbidden_discovery_fields": sorted(FORBIDDEN_DISCOVERY_FIELDS),
        "statement": (
            "The discovery search received structures, amounts, physical conditions, compute settings, "
            "and search budgets only. Target products, recipes, literature yields, and benchmark answers "
            "are reserved for immutable post-run evaluation."
        ),
    }


def scientific_policy() -> dict[str, Any]:
    return {
        "policy_version": POLICY_VERSION,
        "discovery_mode": "blind_physics_only",
        "allowed_discovery_inputs": [
            "reactant identity and starting geometry",
            "amount, concentration, phase, and charge/spin",
            "temperature, pressure, solvent, mixing, and reactor conditions",
            "model selection and finite compute/search budgets",
            "general physical laws, non-targeted local operators, and model parameters learned before benchmark selection",
        ],
        "post_run_only": [
            "expected products and named reaction recipes",
            "literature yields, selectivities, pathways, and reference energies",
            "benchmark pass/fail labels",
        ],
        "forbidden_discovery_fields": sorted(FORBIDDEN_DISCOVERY_FIELDS),
        "knowledge_assisted_mode_available": False,
    }
