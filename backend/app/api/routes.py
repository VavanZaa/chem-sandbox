from __future__ import annotations

from pathlib import Path
import re
import json

import numpy as np
import cantera as ct
import httpx
from ase.io import read
from fastapi import APIRouter, BackgroundTasks, HTTPException

from app.chemistry.bond_detector import BondDetector
from app.chemistry.species_detector import SpeciesDetector
from app.catalysis.microkinetics import export_catmap_setup, simulate_catalysis
from app.catalysis.schema import CatalysisSimulationRequest
from app.experiment.schema import SUBSTANCE_LIBRARY, ExperimentDefinition, default_experiment
from app.experiment.stoichiometry import theoretical_yield
from app.experiment.kinetics import KineticsRequest, simulate_kinetics
from app.experiment.efficiency import EfficiencyRequest, simulate_efficiency
from app.models.registry import hardware_status
from app.multiphase.gas_liquid import simulate_gas_liquid_transfer
from app.multiphase.reaktoro_bridge import (
    ReaktoroUnavailableError,
    equilibrate_with_reaktoro,
    reaktoro_status,
)
from app.multiphase.schema import ElectrolyteEquilibriumRequest, GasLiquidTransferRequest
from app.operations.capabilities import operation_capabilities
from app.operations.cantera_runner import simulate_cantera
from app.operations.cantera_schema import CanteraSimulationRequest
from app.operations.schema import OperationValidationRequest
from app.operations.validation import validate_operation
from app.pathways.rates import build_rate_artifact
from app.pathways.schema import RateArtifactRequest
from app.process.flowsheet import simulate_process_flowsheet
from app.process.idaes_bridge import idaes_capabilities
from app.process.schema import ProcessFlowsheetRequest
from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.hcn_ceiling import hcn_atom_ceiling
from app.prebiotic.reduced_hcn import ReducedHCNRequest, audit_reduced_hcn, simulate_reduced_hcn
from app.prebiotic.candidate_ranking import CandidateRankingJobManager
from app.prebiotic.boltzmann import solve_self_consistent_eedf
from app.prebiotic.capabilities import prebiotic_capabilities
from app.prebiotic.historical_benchmark import historical_material_time_profile
from app.prebiotic.precursor_benchmark import stribling_miller_measurement_profile
from app.prebiotic.recorded_validation import recorded_prebiotic_validation
from app.prebiotic.plasma import (
    calculate_electron_collision_rate,
    validate_electron_collision_dataset,
)
from app.prebiotic.schema import (
    BoltzmannEEDFRequest,
    ElectronCollisionDatasetRequest,
    ElectronCollisionRateRequest,
    PrebioticApparatusRequest,
    FrozenCandidateRankingRequest,
)
from app.protocols.runner import simulate_protocol
from app.protocols.schema import LabProtocolRequest, SBMLInspectionRequest
from app.protocols.sbml_adapter import inspect_sbml
from app.science.blind_policy import scientific_policy
from app.thermochemistry.calculator import (
    evaluate_activation_free_energy,
    evaluate_ensemble,
    evaluate_reaction_free_energy,
)
from app.thermochemistry.harmonic import harmonic_vibrational_corrections
from app.thermochemistry.schema import (
    ActivationFreeEnergyRequest,
    HarmonicCorrectionRequest,
    ReactionFreeEnergyRequest,
    ThermochemicalEnsembleRequest,
)
from app.search.refinement import RefinementRequest, PathKineticsRequest, predict_path_populations
from app.validation.calibration import evaluate_held_out, fit_calibration, freeze_predictions
from app.validation.schema import CalibrationRequest, HeldOutEvaluationRequest, PredictionSetRequest

from .jobs import JobManager


router = APIRouter(prefix="/api")
RUNS_DIR = Path(__file__).resolve().parents[2] / "data" / "runs"
jobs = JobManager(RUNS_DIR)
candidate_ranking_jobs = CandidateRankingJobManager(RUNS_DIR / "prebiotic-ranking")


@router.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "service": "chem-sandbox-api", "version": "0.26.0"}


@router.get("/compute")
def compute_status():
    return hardware_status()


@router.get("/scientific-policy")
def get_scientific_policy():
    return scientific_policy()


@router.get("/operations/capabilities")
def get_operation_capabilities():
    return operation_capabilities()


@router.post("/operations/validate")
def validate_operation_definition(request: OperationValidationRequest):
    try:
        return validate_operation(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/operations/cantera/simulate")
def run_cantera_operation(request: CanteraSimulationRequest):
    try:
        return simulate_cantera(request)
    except (ValueError, ct.CanteraError) as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/thermochemistry/ensemble")
def thermochemical_ensemble(request: ThermochemicalEnsembleRequest):
    try:
        return evaluate_ensemble(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/thermochemistry/activation")
def thermochemical_activation(request: ActivationFreeEnergyRequest):
    try:
        return evaluate_activation_free_energy(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/thermochemistry/reaction")
def thermochemical_reaction(request: ReactionFreeEnergyRequest):
    try:
        return evaluate_reaction_free_energy(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/thermochemistry/harmonic-corrections")
def thermochemical_harmonic_corrections(request: HarmonicCorrectionRequest):
    try:
        return harmonic_vibrational_corrections(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/pathways/rate-artifact")
def pathway_rate_artifact(request: RateArtifactRequest):
    try:
        return build_rate_artifact(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/catalysis/simulate")
def run_catalysis(request: CatalysisSimulationRequest):
    try:
        return simulate_catalysis(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/catalysis/catmap/export")
def export_catalysis_to_catmap(request: CatalysisSimulationRequest):
    try:
        return export_catmap_setup(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.get("/process/capabilities")
def process_capabilities():
    return idaes_capabilities()


@router.get("/prebiotic/capabilities")
def get_prebiotic_capabilities():
    return prebiotic_capabilities()


@router.get("/prebiotic/benchmarks/miller-urey-1953/profile")
def get_miller_urey_1953_profile():
    return historical_material_time_profile()


@router.get("/prebiotic/benchmarks/stribling-miller-1986/profile")
def get_stribling_miller_1986_profile():
    return stribling_miller_measurement_profile()


@router.get("/prebiotic/validation/recorded")
def get_recorded_prebiotic_validation():
    return recorded_prebiotic_validation()


@router.post("/prebiotic/post-discovery/rank")
def start_frozen_candidate_ranking(
    request: FrozenCandidateRankingRequest,
    background_tasks: BackgroundTasks,
):
    job = candidate_ranking_jobs.create()
    background_tasks.add_task(candidate_ranking_jobs.run, job.id, request)
    return job.payload()


@router.get("/prebiotic/post-discovery/rank/{job_id}")
def get_frozen_candidate_ranking(job_id: str):
    job = candidate_ranking_jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Candidate-ranking job not found")
    return job.payload()


@router.post("/prebiotic/post-discovery/rank/{job_id}/held-out-evaluation")
def evaluate_frozen_candidate_ranking_job(job_id: str):
    try:
        return candidate_ranking_jobs.evaluate(job_id)
    except KeyError as error:
        raise HTTPException(status_code=404, detail="Candidate-ranking job not found") from error
    except ValueError as error:
        raise HTTPException(status_code=409, detail=str(error)) from error


@router.post("/prebiotic/apparatus/simulate")
def run_prebiotic_apparatus(request: PrebioticApparatusRequest):
    try:
        return simulate_prebiotic_apparatus(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/prebiotic/hcn/atom-ceiling")
def calculate_hcn_atom_ceiling(request: PrebioticApparatusRequest):
    try:
        return hcn_atom_ceiling(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/prebiotic/hcn/reduced-mechanistic")
def run_reduced_hcn(request: ReducedHCNRequest):
    try:
        return simulate_reduced_hcn(request)
    except (ValueError, ct.CanteraError) as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/prebiotic/hcn/reduced-mechanistic/audit")
def run_reduced_hcn_audit(request: ReducedHCNRequest):
    try:
        return audit_reduced_hcn(request)
    except (ValueError, ct.CanteraError) as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/prebiotic/plasma/cross-sections/validate")
def validate_prebiotic_cross_sections(request: ElectronCollisionDatasetRequest):
    try:
        return validate_electron_collision_dataset(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/prebiotic/plasma/rate-coefficient")
def calculate_prebiotic_collision_rate(request: ElectronCollisionRateRequest):
    try:
        return calculate_electron_collision_rate(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/prebiotic/plasma/eedf/solve")
def solve_prebiotic_eedf(request: BoltzmannEEDFRequest):
    try:
        return solve_self_consistent_eedf(request)
    except (ValueError, RuntimeError) as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/process/simulate")
def run_process_flowsheet(request: ProcessFlowsheetRequest):
    try:
        return simulate_process_flowsheet(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/protocols/simulate")
def run_lab_protocol(request: LabProtocolRequest):
    try:
        return simulate_protocol(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/protocols/sbml/inspect")
def inspect_sbml_model(request: SBMLInspectionRequest):
    try:
        return inspect_sbml(request.sbml)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/validation/freeze")
def freeze_validation_predictions(request: PredictionSetRequest):
    return freeze_predictions(request)


@router.post("/validation/calibrate")
def calibrate_predictions(request: CalibrationRequest):
    try:
        return fit_calibration(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/validation/evaluate")
def evaluate_validation_predictions(request: HeldOutEvaluationRequest):
    try:
        return evaluate_held_out(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.get("/multiphase/capabilities")
def multiphase_capabilities():
    return {
        "native_gas_liquid_transfer": True,
        "electrolyte_equilibrium": reaktoro_status(),
    }


@router.post("/multiphase/gas-liquid")
def gas_liquid_transfer(request: GasLiquidTransferRequest):
    try:
        return simulate_gas_liquid_transfer(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/multiphase/electrolyte/equilibrate")
def electrolyte_equilibrium(request: ElectrolyteEquilibriumRequest):
    try:
        return equilibrate_with_reaktoro(request)
    except ReaktoroUnavailableError as error:
        raise HTTPException(status_code=503, detail=str(error)) from error
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error
    except httpx.HTTPError as error:
        raise HTTPException(status_code=502, detail=str(error)) from error


@router.post("/efficiency")
def efficiency(request: EfficiencyRequest):
    try:
        return simulate_efficiency(request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/experiments/{job_id}/refinement", status_code=202)
def start_refinement(job_id: str, settings: RefinementRequest, background_tasks: BackgroundTasks):
    result = stored_result(job_id)
    if not result.get("graph", {}).get("edges"):
        raise HTTPException(status_code=422, detail="No candidate pathways to refine")
    job = jobs.create()
    background_tasks.add_task(jobs.refine, job.id, result, settings)
    return job.payload()


@router.post("/refinements/{job_id}/populations")
def path_populations(job_id: str, request: PathKineticsRequest):
    result = stored_result(job_id)
    if "records" not in result:
        raise HTTPException(status_code=422, detail="This is not a path-refinement result")
    try:
        return predict_path_populations(result, request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.post("/estimates/stoichiometry")
def estimate_stoichiometry(experiment: ExperimentDefinition) -> dict:
    if experiment.mode != "bicarbonate":
        raise HTTPException(status_code=422, detail="General stoichiometric ceilings require a balanced kinetic network in the Efficiency model")
    return theoretical_yield(experiment)


def stored_result(job_id: str) -> dict:
    if not re.fullmatch(r"[A-Za-z0-9_-]{1,64}", job_id):
        raise HTTPException(status_code=400, detail="Invalid experiment job id")
    job = jobs.get(job_id)
    if job and job.result:
        return job.result
    path = RUNS_DIR / job_id / "result.json"
    if not path.is_file():
        raise HTTPException(status_code=404, detail="Completed experiment not found")
    return json.loads(path.read_text(encoding="utf-8"))


@router.get("/experiments/{job_id}/result")
def experiment_result(job_id: str) -> dict:
    return stored_result(job_id)


@router.post("/experiments/{job_id}/kinetics")
def experiment_kinetics(job_id: str, request: KineticsRequest) -> dict:
    result = stored_result(job_id)
    try:
        return simulate_kinetics(result["graph"], request)
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error


@router.get("/substances")
def substances() -> list[dict[str, object]]:
    return [{"id": key.value, **value} for key, value in SUBSTANCE_LIBRARY.items()]


@router.get("/experiments/default")
def default_experiment_definition() -> dict[str, object]:
    return default_experiment().model_dump(mode="json")


@router.post("/experiments", status_code=202)
def start_experiment(
    experiment: ExperimentDefinition, background_tasks: BackgroundTasks
) -> dict[str, object]:
    job = jobs.create()
    background_tasks.add_task(jobs.run, job.id, experiment)
    return job.payload()


@router.get("/experiments/{job_id}")
def experiment_status(job_id: str) -> dict[str, object]:
    job = jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Experiment job not found")
    return job.payload()


@router.get("/experiments/{job_id}/structures/{filename}")
def experiment_structure(job_id: str, filename: str) -> dict[str, object]:
    """Return one saved XYZ as browser-safe atom and bond data."""

    if not re.fullmatch(r"[A-Za-z0-9_-]{1,64}", job_id):
        raise HTTPException(status_code=400, detail="Invalid experiment job id")
    if Path(filename).name != filename or not filename.endswith(".xyz"):
        raise HTTPException(status_code=400, detail="Invalid structure filename")
    structures_dir = (RUNS_DIR / job_id / "structures").resolve()
    path = (structures_dir / filename).resolve()
    if path.parent != structures_dir or not path.is_file():
        raise HTTPException(status_code=404, detail="Saved structure not found")

    atoms = read(path)
    positions = np.asarray(atoms.positions, dtype=float)
    center = positions.mean(axis=0) if len(positions) else np.zeros(3)
    centered = positions - center
    bonds = BondDetector().detect(atoms)
    species = SpeciesDetector().state_species(atoms)
    return {
        "filename": filename,
        "atoms": [
            {
                "index": index,
                "element": atom.symbol,
                "position": [float(value) for value in centered[index]],
            }
            for index, atom in enumerate(atoms)
        ],
        "bonds": [{"source": source, "target": target} for source, target in sorted(bonds)],
        "species": species,
        "note": "Saved sampled geometry; bonds are inferred from covalent radii, not measured bond orders or proof of chemical identity.",
    }
