"""Compare matched paths using independent potentials and bounded CI-NEB.

Electronic barriers are not free-energy barriers. Convergence is exposed; no
transition state is certified without vibrational/solvation validation.
"""
from pathlib import Path
from typing import Literal
import json
import time

import numpy as np
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import FIRE
from ase.vibrations import Vibrations
from pydantic import BaseModel, ConfigDict, Field, model_validator
from app.chemistry.bond_detector import BondDetector
from app.chemistry.species_detector import SpeciesDetector
from app.experiment.kinetics import KineticsRequest, RateEvidence, simulate_kinetics


class RefinementRequest(BaseModel):
    model_config = ConfigDict(allow_inf_nan=False, extra="forbid")
    engines: list[Literal["mace_mp", "ani2x", "aimnet2", "aimnet2_rxn"]] = Field(
        default_factory=lambda: ["mace_mp", "ani2x", "aimnet2", "aimnet2_rxn"], min_length=1, max_length=4
    )
    device: Literal["auto", "cuda", "cpu"] = "auto"
    images: int = Field(default=5, ge=3, le=9)
    steps: int = Field(default=80, ge=1, le=300)
    endpoint_steps: int = Field(default=120, ge=1, le=500)
    max_edges: int = Field(default=8, ge=1, le=24)
    transition_state_certification: bool = False
    frequency_delta_angstrom: float = Field(default=0.01, ge=0.002, le=0.05)
    frequency_atom_limit: int = Field(default=30, ge=2, le=80)
    coordinate_displacement_angstrom: float = Field(default=0.12, ge=0.03, le=0.30)
    coordinate_follow_steps: int = Field(default=120, ge=10, le=500)

    @model_validator(mode="after")
    def distinct_engines(self):
        if len(set(self.engines)) != len(self.engines):
            raise ValueError("Choose each engine once")
        return self


def _frequency_value_cm_inv(value) -> float:
    value = complex(value)
    return -abs(float(value.imag)) if abs(value.imag) > 1e-8 else float(value.real)


def certify_neb_peak(peak_atoms, left, right, engine, settings: RefinementRequest,
                     output_dir: Path, prefix: str) -> dict:
    """Numerical Hessian plus imaginary-mode endpoint following.

    This is intentionally opt-in because central finite differences require
    many force calls. Endpoint following verifies connectivity on both sides;
    it is not relabelled as a continuous intrinsic reaction coordinate.
    """
    if len(peak_atoms) > settings.frequency_atom_limit:
        return {"status": "skipped_atom_limit", "certified": False,
                "note": f"{len(peak_atoms)} atoms exceeds the configured limit {settings.frequency_atom_limit}."}
    vibration_root = output_dir / "vibrations"
    vibration_root.mkdir(parents=True, exist_ok=True)
    candidate = peak_atoms.copy()
    candidate.calc = engine.ase_calculator(candidate)
    vibration = Vibrations(candidate, name=str(vibration_root / prefix),
                           delta=settings.frequency_delta_angstrom, nfree=2)
    try:
        vibration.run()
        raw_frequencies = vibration.get_frequencies()
        raw_modes = vibration.get_modes()
    except Exception as error:
        return {"status": "frequency_failed", "certified": False,
                "error": f"{type(error).__name__}: {error}"}
    finally:
        vibration.clean()
    frequencies = [_frequency_value_cm_inv(value) for value in raw_frequencies]
    moments = np.asarray(candidate.get_moments_of_inertia(), dtype=float)
    maximum_moment = float(moments.max()) if moments.size else 0.0
    minimum_moment = float(moments.min()) if moments.size else 0.0
    linear = bool(maximum_moment > 0 and minimum_moment / maximum_moment < 1e-6)
    remove_count = 5 if linear else 6
    removed = set(sorted(range(len(frequencies)), key=lambda index: abs(frequencies[index]))[:remove_count])
    kept_indices = [index for index in range(len(frequencies)) if index not in removed]
    vibrational_frequencies = [frequencies[index] for index in kept_indices]
    imaginary = [index for index in kept_indices if frequencies[index] < -20.0]
    evidence = {
        "status": "hessian_complete",
        "certified": False,
        "linear": linear,
        "frequencies_cm_inv": vibrational_frequencies,
        "significant_imaginary_modes_cm_inv": [frequencies[index] for index in imaginary],
        "frequency_method": "ASE central finite-difference Hessian of the selected neural-potential energy surface",
        "frequency_delta_angstrom": settings.frequency_delta_angstrom,
        "coordinate_following_verified": False,
        "continuous_irc_calculated": False,
    }
    if len(imaginary) != 1:
        evidence["status"] = "wrong_imaginary_mode_count"
        evidence["note"] = f"Expected one significant imaginary mode, found {len(imaginary)}."
        return evidence
    mode = np.asarray(raw_modes[imaginary[0]], dtype=float)
    maximum = float(np.linalg.norm(mode, axis=1).max(initial=0.0))
    if maximum <= 0 or not np.isfinite(mode).all():
        evidence["status"] = "invalid_imaginary_mode"
        return evidence
    displacement = mode * (settings.coordinate_displacement_angstrom / maximum)
    followed = []
    try:
        for sign in (-1.0, 1.0):
            displaced = peak_atoms.copy()
            displaced.positions += sign * displacement
            followed.append(engine.optimize(displaced, fmax=.08, steps=settings.coordinate_follow_steps))
    except Exception as error:
        evidence["status"] = "coordinate_following_failed"
        evidence["error"] = f"{type(error).__name__}: {error}"
        return evidence
    species = SpeciesDetector()
    endpoint_signatures = {species.state_signature(left.atoms), species.state_signature(right.atoms)}
    followed_signatures = {species.state_signature(item.atoms) for item in followed}
    followed_converged = all(item.converged for item in followed)
    verified = followed_converged and len(endpoint_signatures) == 2 and followed_signatures == endpoint_signatures
    evidence.update(
        status="certified_connectivity_following" if verified else "coordinate_following_endpoint_mismatch",
        certified=verified,
        coordinate_following_verified=verified,
        coordinate_following_converged=followed_converged,
        note=(
            "The unique imaginary mode relaxes toward both mapped endpoint connectivities. This is a connectivity-following check, not a continuous IRC."
            if verified
            else "Imaginary-mode displacements did not recover both mapped endpoint connectivities."
        ),
    )
    return evidence


def safe_structure(reference: str | None, runs_dir: Path):
    if not reference:
        raise ValueError("No paired geometry was saved for this edge; rerun the v0.3 search first")
    path = (runs_dir.parent / reference).resolve()
    if not path.is_relative_to(runs_dir.resolve()) or path.suffix != ".xyz" or not path.is_file():
        raise ValueError("Invalid saved structure reference")
    return read(path)


def refine_path(before, after, engine, settings: RefinementRequest, output_dir: Path, prefix: str):
    if before.get_chemical_symbols() != after.get_chemical_symbols():
        raise ValueError("NEB endpoints must have the same mapped atom inventory")
    species = SpeciesDetector()
    left = engine.optimize(before, steps=settings.endpoint_steps, fmax=.08)
    right = engine.optimize(after, steps=settings.endpoint_steps, fmax=.08)
    base = {"delta_energy_ev": right.energy_ev - left.energy_ev,
            "endpoint_converged": left.converged and right.converged,
            "source_max_force": left.max_force, "target_max_force": right.max_force,
            "eligible_for_exploratory_kinetics": False}
    if species.state_signature(left.atoms) != species.state_signature(before) or species.state_signature(right.atoms) != species.state_signature(after):
        return {**base, "status": "endpoint_connectivity_changed", "note": "Relaxation left a searched state; no barrier assigned to that original edge."}
    images = [left.atoms.copy()] + [left.atoms.copy() for _ in range(settings.images - 2)] + [right.atoms.copy()]
    calc = engine.ase_calculator(left.atoms)
    for atoms in images:
        atoms.calc = calc
    neb = NEB(images, climb=True, method="improvedtangent", allow_shared_calculator=True,
              remove_rotation_and_translation=True)
    neb.interpolate(method="idpp")
    optimizer = FIRE(neb, logfile=None, maxstep=.1)
    converged = bool(optimizer.run(fmax=.12, steps=settings.steps))
    energies = np.array([atoms.get_potential_energy() for atoms in images], dtype=float)
    forces = np.asarray(neb.get_forces())
    if not np.isfinite(energies).all() or not np.isfinite(forces).all():
        raise ValueError("Non-finite path calculation")
    if any(BondDetector().minimum_pair_distance(atoms) < .45 for atoms in images):
        raise ValueError("A path image contains an atom clash")
    peak = int(energies.argmax())
    transition_state_evidence = None
    if settings.transition_state_certification and converged and base["endpoint_converged"] and 0 < peak < len(images) - 1:
        transition_state_evidence = certify_neb_peak(
            images[peak], left, right, engine, settings, output_dir, prefix
        )
    references = []
    structures = output_dir / "structures"
    structures.mkdir(parents=True, exist_ok=True)
    for i, atoms in enumerate(images):
        atoms.calc = None
        path = structures / f"{prefix}_image_{i:02d}.xyz"
        write(path, atoms)
        references.append(str(path.relative_to(output_dir.parent.parent)).replace("\\", "/"))
    interior = 0 < peak < len(images) - 1
    if transition_state_evidence is not None:
        transition_state_evidence["structure_reference"] = references[peak]
    return {**base, "status": "converged_electronic_path" if converged and base["endpoint_converged"] and interior else "unverified_path",
            "neb_converged": converged, "neb_steps": optimizer.nsteps,
            "max_neb_force": float(np.linalg.norm(forces, axis=1).max()),
            "electronic_barrier_ev": float(energies.max() - energies[0]),
            "reverse_electronic_barrier_ev": float(energies.max() - energies[-1]),
            "relative_image_energies_ev": (energies - energies[0]).tolist(), "image_references": references,
            "transition_state_evidence": transition_state_evidence,
            "eligible_for_exploratory_kinetics": bool(converged and base["endpoint_converged"] and interior),
            "note": "CI-NEB potential-energy path. No Hessian, zero-point, entropy, solvent or standard-state correction; NOT a verified Gibbs activation barrier."}


def compare_engines(result, settings, engine_factory, runs_dir, output_dir, progress):
    started = time.perf_counter()
    output_dir.mkdir(parents=True, exist_ok=True)
    edges = result["graph"]["edges"][:settings.max_edges]
    records, engines = [], list(dict.fromkeys(settings.engines))
    for engine_id in engines:
        engine = engine_factory(engine_id, settings.device)
        for edge in edges:
            progress("refining", int(95 * len(records) / max(1, len(edges) * len(engines))), f"{engine_id} · refining {edge['id']}")
            try:
                before = safe_structure(edge.get("source_reference"), runs_dir)
                after = safe_structure(edge.get("target_reference"), runs_dir)
                payload = refine_path(before, after, engine, settings, output_dir, f"{engine_id}_{edge['id']}")
                record = {"engine": engine_id, "edge_id": edge["id"], "metadata": engine.metadata, **payload}
            except Exception as error:
                record = {"engine": engine_id, "edge_id": edge["id"], "status": "failed_or_out_of_domain",
                          "error": f"{type(error).__name__}: {error}", "eligible_for_exploratory_kinetics": False}
            records.append(record)
            with (output_dir / "events.jsonl").open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(record) + "\n")
    comparisons = []
    for edge in edges:
        rows = [row for row in records if row["edge_id"] == edge["id"] and "electronic_barrier_ev" in row]
        barriers = [row["electronic_barrier_ev"] for row in rows]
        comparisons.append({"edge_id": edge["id"], "models_with_path": len(rows),
                            "barrier_range_ev": [min(barriers), max(barriers)] if barriers else None,
                            "all_paths_eligible": len(rows) == len(engines) and all(row["eligible_for_exploratory_kinetics"] for row in rows)})
    payload = {"status": "completed", "graph": result["graph"], "settings": settings.model_dump(), "records": records,
               "comparisons": comparisons, "unrefined_edges": len(result["graph"]["edges"]) - len(edges),
               "runtime_seconds": time.perf_counter() - started, "calibrated_product_probabilities": None,
               "warnings": ["Models evaluate the same discovery network and paired atom mappings; discovery is not independently repeated with every model.",
                            "Never compare absolute energies between engines. Compare within-engine energy differences for the same mapped endpoints.",
                            "Inter-model ranges measure disagreement, not confidence intervals or calibrated product probabilities.",
                            "Failed or unconverged paths remain unknown, not zero-rate reactions. Partial networks cannot supply complete kinetics.",
                            "Paths use vacuum neural potentials. Temperature, solvent, pH, entropy and pressure free-energy corrections are not calculated."]}
    (output_dir / "result.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


class PathKineticsRequest(BaseModel):
    model_config = ConfigDict(allow_inf_nan=False, extra="forbid")
    temperature_K: float = Field(default=298.15, ge=250, le=450)
    duration_s: float = Field(default=1, gt=0, le=1e8)
    initial_state_id: str
    target_state_id: str
    acknowledge_electronic_barrier_proxy: Literal[True]


def predict_path_populations(refinement, request):
    graph = refinement["graph"]
    if request.target_state_id not in {node["id"] for node in graph["nodes"]}:
        raise ValueError("Unknown target state")
    predictions, unavailable = [], []
    for engine in refinement["settings"]["engines"]:
        records = [row for row in refinement["records"] if row["engine"] == engine and row["eligible_for_exploratory_kinetics"]]
        if {row["edge_id"] for row in records} != {edge["id"] for edge in graph["edges"]}:
            unavailable.append({"engine": engine, "reason": "A converged electronic path is required for every directed edge; unknown barriers cannot be set to zero."})
            continue
        try:
            rates = [RateEvidence(edge_id=row["edge_id"], activation_free_energy_ev=row["electronic_barrier_ev"],
                                  temperature_K=request.temperature_K, basis="hypothetical",
                                  source=f"{engine} CI-NEB electronic energy used as unvalidated ΔG‡ proxy; no entropy/solvent correction") for row in records]
            kinetics = simulate_kinetics(graph, KineticsRequest(initial_state_id=request.initial_state_id,
                temperature_K=request.temperature_K, duration_s=request.duration_s, rates=rates, acknowledge_conditional_model=True))
            final = kinetics["samples"][-1]["fractions"]
            predictions.append({"engine": engine, "state_conversion": 1 - final[request.initial_state_id],
                                "target_state_population": final[request.target_state_id], "kinetics": kinetics})
        except ValueError as error:
            unavailable.append({"engine": engine, "reason": str(error)})
    return {"status": "exploratory_only" if predictions else "insufficient_barrier_evidence", "predictions": predictions,
            "unavailable": unavailable, "actual_yield": None, "calibrated_product_probabilities": None,
            "warning": "These are model-conditioned cluster-state populations with ΔG‡ ≈ ΔE‡, not bulk efficiency or validated probabilities. Use the bulk efficiency model with condition-appropriate kinetic evidence for molar yields."}
