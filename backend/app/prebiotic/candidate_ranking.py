"""Background, formula-local energy ranking of a frozen blind candidate set."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
import hashlib
import json
from pathlib import Path
import threading
import traceback
import uuid
from typing import Any, Callable

from ase import Atoms
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from app.models.registry import make_engine

from .post_discovery import freeze_discovery_artifact
from .schema import FrozenCandidateRankingRequest


EV_TO_KJ_MOL = 96.48533212331002
CARBOXYLIC_ACID = Chem.MolFromSmarts("C(=O)[O;H1]")
PRIMARY_OR_SECONDARY_AMINE = Chem.MolFromSmarts("[NX3;H1,H2]")
HISTORICAL_BENCHMARK_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "benchmarks"
    / "miller_urey_1953_v1.json"
)


def select_amino_carboxylic_candidates(discovery: dict[str, Any]) -> list[dict[str, Any]]:
    """Apply a generic functional-class filter only after blind enumeration."""

    selected: list[dict[str, Any]] = []
    for record in discovery.get("formula_records", []):
        for candidate in record.get("candidates", []):
            molecule = Chem.MolFromSmiles(candidate.get("smiles", ""))
            if (
                molecule is not None
                and CARBOXYLIC_ACID is not None
                and PRIMARY_OR_SECONDARY_AMINE is not None
                and molecule.HasSubstructMatch(CARBOXYLIC_ACID)
                and molecule.HasSubstructMatch(PRIMARY_OR_SECONDARY_AMINE)
            ):
                selected.append(
                    {
                        "formula": record["formula"],
                        "smiles": candidate["smiles"],
                        "stereochemistry": candidate.get("stereochemistry", []),
                    }
                )
    return sorted(selected, key=lambda item: (item["formula"], item["smiles"]))


def _starting_geometry(smiles: str, seed: int) -> Atoms:
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed % (2**31 - 1)
    if AllChem.EmbedMolecule(molecule, parameters) < 0:
        raise ValueError(f"RDKit could not embed candidate {smiles}")
    if AllChem.MMFFHasAllMoleculeParams(molecule):
        AllChem.MMFFOptimizeMolecule(molecule, maxIters=300)
    atoms = Atoms(
        [atom.GetSymbol() for atom in molecule.GetAtoms()],
        positions=np.asarray(molecule.GetConformer().GetPositions(), dtype=float),
    )
    atoms.info.update(total_charge=0, radical_electrons=0, has_ions=False)
    return atoms


def _achiral_connectivity_key(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid frozen candidate SMILES: {smiles}")
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)


def build_paired_starting_geometries(
    candidates: list[dict[str, Any]], seed_base: int
) -> tuple[dict[str, Atoms], list[dict[str, Any]]]:
    """Use reflected coordinates for simple R/S pairs in an achiral model."""

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for candidate in candidates:
        grouped.setdefault(
            (candidate["formula"], _achiral_connectivity_key(candidate["smiles"])), []
        ).append(candidate)
    starting: dict[str, Atoms] = {}
    mirror_pairs: list[dict[str, Any]] = []
    for group_index, ((formula, key), group) in enumerate(sorted(grouped.items())):
        descriptors = {
            center["descriptor"]
            for candidate in group
            for center in candidate.get("stereochemistry", [])
        }
        if len(group) == 2 and descriptors == {"R", "S"}:
            ordered = sorted(group, key=lambda item: item["smiles"])
            base = _starting_geometry(ordered[0]["smiles"], seed_base + group_index)
            reflected = base.copy()
            positions = reflected.get_positions()
            positions[:, 0] *= -1
            reflected.set_positions(positions)
            starting[ordered[0]["smiles"]] = base
            starting[ordered[1]["smiles"]] = reflected
            mirror_pairs.append(
                {
                    "formula": formula,
                    "achiral_connectivity_key": key,
                    "smiles": [item["smiles"] for item in ordered],
                }
            )
        else:
            for candidate_index, candidate in enumerate(group):
                starting[candidate["smiles"]] = _starting_geometry(
                    candidate["smiles"],
                    seed_base + group_index * 31 + candidate_index,
                )
    return starting, mirror_pairs


def rank_frozen_candidates(
    request: FrozenCandidateRankingRequest,
    progress: Callable[[str, int, str], None] | None = None,
) -> dict[str, Any]:
    digest, _ = freeze_discovery_artifact(request.discovery)
    if digest != request.discovery_artifact_sha256:
        raise ValueError("Frozen discovery hash does not match the submitted artifact")
    candidates = select_amino_carboxylic_candidates(request.discovery)
    if not candidates:
        raise ValueError("The frozen artifact contains no neutral amino-carboxylic candidates")
    seed_base = int(digest[:8], 16)
    starting, mirror_pairs = build_paired_starting_geometries(candidates, seed_base)
    engine_results: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    total_work = len(request.engines) * len(candidates)
    completed_work = 0
    for engine_id in request.engines:
        if progress:
            progress(engine_id, int(100 * completed_work / total_work), f"Loading {engine_id}")
        engine = None
        try:
            engine = make_engine(engine_id, request.device)
            records = []
            for item in candidates:
                if progress:
                    progress(
                        engine_id,
                        int(100 * completed_work / total_work),
                        f"Optimizing {item['smiles']}",
                    )
                result = engine.optimize(
                    starting[item["smiles"]],
                    fmax=request.force_threshold_eV_A,
                    steps=request.optimization_steps,
                )
                records.append(
                    {
                        **item,
                        "energy_eV": result.energy_ev,
                        "optimizer_converged": result.converged,
                        "optimization_steps": result.steps,
                        "max_force_eV_A": result.max_force,
                    }
                )
                completed_work += 1
            mirror_diagnostics = []
            for pair in mirror_pairs:
                pair_records = [
                    item for item in records if item["smiles"] in pair["smiles"]
                ]
                mean_energy = float(np.mean([item["energy_eV"] for item in pair_records]))
                difference_kJ_mol = (
                    abs(pair_records[0]["energy_eV"] - pair_records[1]["energy_eV"])
                    * EV_TO_KJ_MOL
                )
                for item in pair_records:
                    item["symmetry_corrected_energy_eV"] = mean_energy
                    item["achiral_mirror_degeneracy_enforced"] = True
                mirror_diagnostics.append(
                    {
                        **pair,
                        "raw_optimized_energy_difference_kJ_mol": difference_kJ_mol,
                        "degenerate_rank_enforced": True,
                    }
                )
            for item in records:
                item.setdefault("symmetry_corrected_energy_eV", item["energy_eV"])
                item.setdefault("achiral_mirror_degeneracy_enforced", False)
            for formula in sorted({item["formula"] for item in records}):
                group = [item for item in records if item["formula"] == formula]
                reference = min(item["symmetry_corrected_energy_eV"] for item in group)
                ordered = sorted(group, key=lambda item: item["symmetry_corrected_energy_eV"])
                previous_energy = None
                previous_rank = 0.0
                for position, item in enumerate(ordered, start=1):
                    energy = item["symmetry_corrected_energy_eV"]
                    if previous_energy is not None and abs(energy - previous_energy) <= 1e-8:
                        rank = previous_rank
                    else:
                        tied = [
                            candidate
                            for candidate in ordered
                            if abs(candidate["symmetry_corrected_energy_eV"] - energy) <= 1e-8
                        ]
                        rank = position + (len(tied) - 1) / 2
                    item["relative_energy_kJ_mol"] = (
                        energy - reference
                    ) * EV_TO_KJ_MOL
                    item["formula_local_rank"] = rank
                    previous_energy = energy
                    previous_rank = rank
            engine_results.append(
                {
                    "engine": engine_id,
                    "metadata": engine.metadata,
                    "candidates": records,
                    "all_optimizations_converged": all(
                        item["optimizer_converged"] for item in records
                    ),
                    "mirror_pair_diagnostics": mirror_diagnostics,
                }
            )
        except Exception as error:
            failures.append(
                {
                    "engine": engine_id,
                    "error": f"{type(error).__name__}: {error}",
                }
            )
            completed_work += len(candidates)
        finally:
            del engine
            try:
                import torch

                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            except ImportError:
                pass

    if not engine_results:
        raise RuntimeError("All requested candidate-ranking engines failed")
    consensus = []
    for candidate in candidates:
        ranks = []
        relative_energies = []
        for result in engine_results:
            record = next(
                item
                for item in result["candidates"]
                if item["smiles"] == candidate["smiles"]
            )
            ranks.append(record["formula_local_rank"])
            relative_energies.append(record["relative_energy_kJ_mol"])
        consensus.append(
            {
                **candidate,
                "engine_count": len(ranks),
                "mean_formula_local_rank": float(np.mean(ranks)),
                "rank_range": [min(ranks), max(ranks)],
                "relative_energy_range_kJ_mol": [
                    float(min(relative_energies)),
                    float(max(relative_energies)),
                ],
            }
        )
    consensus.sort(
        key=lambda item: (item["formula"], item["mean_formula_local_rank"], item["smiles"])
    )
    if progress:
        progress("complete", 100, "Formula-local ranking complete")
    return {
        "status": "frozen_candidate_energy_ranking_complete",
        "discovery_artifact_sha256": digest,
        "candidate_filter": "neutral primary/secondary amine plus carboxylic acid",
        "candidate_count": len(candidates),
        "engines_requested": request.engines,
        "engines_succeeded": [item["engine"] for item in engine_results],
        "engine_failures": failures,
        "engine_results": engine_results,
        "achiral_mirror_pairs": mirror_pairs,
        "formula_local_consensus": consensus,
        "scientific_firewall": {
            "discovery_hash_verified_before_ranking": True,
            "candidate_set_changed_by_ranking": False,
            "historical_product_labels_available_to_ranker": False,
            "energy_scales_mixed_between_engines": False,
            "achiral_enantiomer_degeneracy_enforced": True,
            "solution_free_energy_claimed": False,
            "pathway_or_yield_claimed": False,
        },
        "limitations": [
            "Vacuum electronic potential energies are compared only among candidates with the same formula and within one engine.",
            "A consensus rank is model agreement, not a probability, equilibrium population, reaction pathway, or yield.",
            "Solvation, protonation, conformer completeness, entropy, transition states, and analytical detectability are not included.",
        ],
    }


def evaluate_frozen_candidate_ranking(
    ranking_result: dict[str, Any], top_k: int = 3
) -> dict[str, Any]:
    """Reveal held-out identities only after a compact prediction is frozen."""

    if ranking_result.get("status") != "frozen_candidate_energy_ranking_complete":
        raise ValueError("A completed frozen candidate ranking is required")
    consensus = ranking_result.get("formula_local_consensus", [])
    if not consensus:
        raise ValueError("The ranking has no formula-local consensus candidates")
    prediction_artifact = {
        "scope": "post_discovery_class_conditional_vacuum_energy_ranking",
        "discovery_artifact_sha256": ranking_result["discovery_artifact_sha256"],
        "candidate_filter": ranking_result["candidate_filter"],
        "engines_succeeded": ranking_result["engines_succeeded"],
        "formula_local_consensus": [
            {
                "formula": item["formula"],
                "smiles": item["smiles"],
                "mean_formula_local_rank": item["mean_formula_local_rank"],
                "rank_range": item["rank_range"],
            }
            for item in consensus
        ],
    }
    frozen_bytes = json.dumps(
        prediction_artifact,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
    prediction_sha256 = hashlib.sha256(frozen_bytes).hexdigest()

    # This is the first read of historical identities in this evaluation path.
    # The prediction artifact and digest above are immutable inputs to scoring.
    benchmark = json.loads(HISTORICAL_BENCHMARK_PATH.read_text(encoding="utf-8"))
    by_smiles = {item["smiles"]: item for item in consensus}
    observations = []
    for observation in benchmark["evaluation_only_observations"]:
        hits = [
            by_smiles[smiles]
            for smiles in observation["acceptable_smiles"]
            if smiles in by_smiles
        ]
        best_rank = min(
            (float(item["mean_formula_local_rank"]) for item in hits),
            default=None,
        )
        observations.append(
            {
                "name": observation["name"],
                "candidate_present_after_class_filter": bool(hits),
                "candidate_smiles": [item["smiles"] for item in hits],
                "formula": hits[0]["formula"] if hits else None,
                "best_mean_formula_local_rank": best_rank,
                "rank_ranges": [item["rank_range"] for item in hits],
                "top_k_recovered": best_rank is not None and best_rank <= top_k,
                "source": observation["source"],
            }
        )
    denominator = len(observations)
    class_coverage = sum(
        item["candidate_present_after_class_filter"] for item in observations
    ) / denominator
    top_k_recall = sum(item["top_k_recovered"] for item in observations) / denominator
    return {
        "status": "class_conditional_held_out_evaluation_complete_not_reproduction",
        "benchmark_id": benchmark["benchmark_id"],
        "scope": "post_discovery_class_conditional",
        "prediction_artifact_sha256": prediction_sha256,
        "prediction_artifact_bytes": len(frozen_bytes),
        "discovery_artifact_sha256": ranking_result["discovery_artifact_sha256"],
        "top_k": top_k,
        "class_filter_coverage": class_coverage,
        "top_k_formula_local_recall": top_k_recall,
        "observation_comparison": observations,
        "validated_reproduction": False,
        "scientific_firewall": {
            "prediction_frozen_before_benchmark_read": True,
            "historical_product_labels_available_to_ranker": False,
            "benchmark_feedback_to_prediction": False,
            "class_filter_selected_post_discovery": True,
            "global_blind_product_prediction_claimed": False,
            "pathway_probability_or_yield_claimed": False,
        },
        "limitations": [
            "This evaluates an amino-carboxylic functional class selected after blind enumeration; it is not a global blind product prediction.",
            "Formula-local vacuum stability ranks do not measure formation kinetics, destruction, solubility, analytical response, abundance, or yield.",
            "A product can rank poorly in thermodynamic stability and still be formed by a kinetically favored pathway.",
        ],
    }


@dataclass(slots=True)
class CandidateRankingJob:
    id: str
    status: str = "queued"
    stage: str = "queued"
    progress: int = 0
    message: str = "Waiting for the model-ranking compute slot"
    created_at: str = field(default_factory=lambda: datetime.now(UTC).isoformat())
    updated_at: str = field(default_factory=lambda: datetime.now(UTC).isoformat())
    result: dict[str, Any] | None = None
    evaluation: dict[str, Any] | None = None
    error: str | None = None

    def payload(self) -> dict[str, Any]:
        return asdict(self)


class CandidateRankingJobManager:
    def __init__(self, runs_dir: Path | None = None) -> None:
        self.jobs: dict[str, CandidateRankingJob] = {}
        self.lock = threading.RLock()
        self.compute_lock = threading.Lock()
        self.runs_dir = runs_dir
        if self.runs_dir is not None:
            self.runs_dir.mkdir(parents=True, exist_ok=True)

    def _save(self, job: CandidateRankingJob) -> None:
        if self.runs_dir is None:
            return
        directory = self.runs_dir / job.id
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "job.json").write_text(
            json.dumps(job.payload(), indent=2, ensure_ascii=False), encoding="utf-8"
        )

    def create(self) -> CandidateRankingJob:
        job = CandidateRankingJob(id=uuid.uuid4().hex[:12])
        with self.lock:
            self.jobs[job.id] = job
        self._save(job)
        return job

    def get(self, job_id: str) -> CandidateRankingJob | None:
        with self.lock:
            job = self.jobs.get(job_id)
        if job is not None or self.runs_dir is None:
            return job
        saved = self.runs_dir / job_id / "job.json"
        if not saved.is_file():
            return None
        payload = json.loads(saved.read_text(encoding="utf-8"))
        job = CandidateRankingJob(**payload)
        with self.lock:
            self.jobs[job_id] = job
        return job

    def _progress(self, job: CandidateRankingJob, stage: str, percent: int, message: str) -> None:
        with self.lock:
            job.stage = stage
            job.progress = percent
            job.message = message
            job.updated_at = datetime.now(UTC).isoformat()

    def run(self, job_id: str, request: FrozenCandidateRankingRequest) -> None:
        job = self.get(job_id)
        if job is None:
            return
        with self.compute_lock:
            try:
                with self.lock:
                    job.status = "running"
                result = rank_frozen_candidates(
                    request,
                    lambda stage, percent, message: self._progress(
                        job, stage, percent, message
                    ),
                )
                with self.lock:
                    job.status = "completed"
                    job.result = result
                self._progress(job, "complete", 100, "Candidate energy ranking ready")
                self._save(job)
            except Exception as error:
                with self.lock:
                    job.status = "failed"
                    job.error = f"{type(error).__name__}: {error}"
                self._progress(job, "failed", 100, "Candidate energy ranking failed")
                self._save(job)
                traceback.print_exc()

    def evaluate(self, job_id: str) -> dict[str, Any]:
        job = self.get(job_id)
        if job is None:
            raise KeyError(job_id)
        if job.status != "completed" or job.result is None:
            raise ValueError("Candidate-ranking job must complete before held-out evaluation")
        evaluation = evaluate_frozen_candidate_ranking(job.result)
        with self.lock:
            job.evaluation = evaluation
            job.updated_at = datetime.now(UTC).isoformat()
        self._save(job)
        if self.runs_dir is not None:
            directory = self.runs_dir / job.id
            (directory / "held-out-evaluation.json").write_text(
                json.dumps(evaluation, indent=2, ensure_ascii=False), encoding="utf-8"
            )
        return evaluation
