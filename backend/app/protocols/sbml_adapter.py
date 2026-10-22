"""SBML network export, validation, and libRoadRunner cross-checks."""
from __future__ import annotations

from hashlib import sha256
from typing import Any

import libsbml
import roadrunner

from .schema import LabProtocolRequest


def _set(result: int, description: str) -> None:
    if result != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise ValueError(f"libSBML failed to {description}: code {result}")


def _validation(document: libsbml.SBMLDocument) -> dict[str, Any]:
    document.checkConsistency()
    messages = []
    error_count = 0
    for index in range(document.getNumErrors()):
        error = document.getError(index)
        severity = int(error.getSeverity())
        if severity >= libsbml.LIBSBML_SEV_ERROR:
            error_count += 1
        if len(messages) < 50:
            messages.append({
                "severity": error.getSeverityAsString(),
                "category": error.getCategoryAsString(),
                "message": error.getMessage(),
                "line": error.getLine(),
            })
    return {
        "status": "valid" if error_count == 0 else "invalid",
        "errors": error_count,
        "diagnostics": messages,
        "libsbml_version": libsbml.getLibSBMLDottedVersion(),
    }


def build_sbml(request: LabProtocolRequest, prepared: dict[str, Any]) -> dict[str, Any]:
    document = libsbml.SBMLDocument(3, 2)
    model = document.createModel()
    _set(model.setId("chem_sandbox_protocol"), "set model id")
    _set(model.setName(request.id), "set model name")
    _set(model.setTimeUnits("second"), "set time units")
    _set(model.setExtentUnits("mole"), "set extent units")
    _set(model.setSubstanceUnits("mole"), "set substance units")
    _set(model.setVolumeUnits("litre"), "set volume units")

    compartment = model.createCompartment()
    _set(compartment.setId("vessel"), "set compartment id")
    _set(compartment.setConstant(False), "make vessel volume variable")
    _set(compartment.setSpatialDimensions(3), "set compartment dimensions")
    _set(compartment.setUnits("litre"), "set compartment units")
    _set(compartment.setSize(request.initial_volume.as_unit("L")), "set initial volume")

    sbml_species: dict[str, str] = {}
    initial_amounts = {key: value.as_unit("mol") for key, value in request.initial_amounts.items()}
    for number, item in enumerate(request.species):
        sbml_id = f"s_{number}"
        sbml_species[item.id] = sbml_id
        species = model.createSpecies()
        _set(species.setId(sbml_id), "set species id")
        _set(species.setName(f"{item.id}: {item.name}"), "set species name")
        _set(species.setCompartment("vessel"), "set species compartment")
        _set(species.setInitialAmount(initial_amounts.get(item.id, 0.0)), "set species initial amount")
        _set(species.setHasOnlySubstanceUnits(True), "use amount units")
        _set(species.setBoundaryCondition(False), "set floating species")
        _set(species.setConstant(False), "set variable species")
        _set(species.setSubstanceUnits("mole"), "set species units")

    temperature = model.createParameter()
    _set(temperature.setId("temperature"), "set temperature id")
    _set(temperature.setValue(request.initial_temperature.as_unit("K")), "set initial temperature")
    _set(temperature.setConstant(False), "make temperature variable")
    _set(temperature.setUnits("kelvin"), "set temperature units")

    parameter_ids: dict[tuple[str, str], str] = {}
    active_ids: dict[str, str] = {}

    def parameter(reaction_number: int, reaction_id: str, name: str, value: float, units: str | None = None) -> str:
        identifier = f"p_{reaction_number}_{name}"
        item = model.createParameter()
        _set(item.setId(identifier), f"set {name} id")
        _set(item.setName(f"{reaction_id}: {name}"), f"set {name} name")
        _set(item.setValue(value), f"set {name} value")
        _set(item.setConstant(True), f"set {name} constant")
        if units:
            _set(item.setUnits(units), f"set {name} units")
        parameter_ids[(reaction_id, name)] = identifier
        return identifier

    for reaction_number, (definition, metadata) in enumerate(zip(request.reactions, prepared["reaction_metadata"])):
        reaction = model.createReaction()
        reaction_id = f"r_{reaction_number}"
        _set(reaction.setId(reaction_id), "set reaction id")
        _set(reaction.setName(f"{definition.id}: {definition.name}"), "set reaction name")
        _set(reaction.setReversible(False), "set irreversible reaction")
        for species_id, coefficient in definition.stoichiometry.items():
            reference = reaction.createReactant() if coefficient < 0 else reaction.createProduct()
            _set(reference.setSpecies(sbml_species[species_id]), "set reaction species")
            _set(reference.setStoichiometry(abs(coefficient)), "set stoichiometry")
            _set(reference.setConstant(True), "set constant stoichiometry")
        active = parameter(reaction_number, definition.id, "active", 0.0, "dimensionless")
        active_ids[definition.id] = active
        orders = metadata["reactant_orders"]
        if metadata["rate_kind"] == "mass_action":
            k = parameter(reaction_number, definition.id, "k", metadata["rate_constant"])
            activity = " * ".join(
                f"({sbml_species[prepared['ids'][species_index]]}/vessel)^{order}"
                for species_index, order in orders.items()
            )
            formula = f"{active} * {k} * vessel * {activity}"
        elif metadata["rate_kind"] == "eyring_mass_action":
            dg = parameter(reaction_number, definition.id, "dg", metadata["activation_free_energy_kj_mol"])
            tref = parameter(reaction_number, definition.id, "tref", metadata["reference_temperature_K"], "kelvin")
            ds = parameter(reaction_number, definition.id, "ds", metadata["activation_entropy_j_mol_K"])
            activity = " * ".join(
                f"({sbml_species[prepared['ids'][species_index]]}/vessel)^{order}"
                for species_index, order in orders.items()
            )
            formula = (
                f"{active} * (1.380649e-23*temperature/6.62607015e-34) * "
                f"exp(-(({dg}-(temperature-{tref})*{ds}/1000)*1000)/(8.31446261815324*temperature)) "
                f"* vessel * {activity}"
            )
        else:
            kcat = parameter(reaction_number, definition.id, "kcat", metadata["kcat_s_inv"], "per_second")
            km = parameter(reaction_number, definition.id, "km", metadata["km_mol_L"])
            substrate = sbml_species[prepared["ids"][metadata["substrate_index"]]]
            catalyst = sbml_species[prepared["ids"][metadata["catalyst_index"]]]
            formula = f"{active} * {kcat} * {catalyst} * ({substrate}/vessel)/({km}+{substrate}/vessel)"
        law = reaction.createKineticLaw()
        math_ast = libsbml.parseL3Formula(formula)
        if math_ast is None:
            raise ValueError(f"Could not translate rate law for {definition.id} to SBML")
        _set(law.setMath(math_ast), "set kinetic law")

    document_xml = libsbml.writeSBMLToString(document)
    return {
        "document": document_xml,
        "sha256": sha256(document_xml.encode("utf-8")).hexdigest(),
        "validation": _validation(document),
        "species_ids": sbml_species,
        "active_parameter_ids": active_ids,
    }


def cross_check_roadrunner(
    request: LabProtocolRequest,
    prepared: dict[str, Any],
    sbml: str,
    native_final: dict[str, float],
) -> dict[str, Any]:
    ids = prepared["ids"]
    sbml_species = {species_id: f"s_{number}" for number, species_id in enumerate(ids)}
    active = {reaction.id: f"p_{number}_active" for number, reaction in enumerate(request.reactions)}
    try:
        runner = roadrunner.RoadRunner(sbml)
        runner.integrator = "cvode"
        runner.integrator.absolute_tolerance = 1e-13
        runner.integrator.relative_tolerance = 1e-9
        current_time = 0.0
        volume = request.initial_volume.as_unit("L")
        for stage in request.stages:
            if stage.temperature is not None:
                runner["temperature"] = stage.temperature.as_unit("K")
            if stage.volume is not None:
                volume = stage.volume.as_unit("L")
                runner["vessel"] = volume
            for addition in stage.additions:
                key = sbml_species[addition.species_id]
                runner[key] = float(runner[key]) + addition.amount.as_unit("mol")
            if stage.operation == "sample":
                fraction = 1 - stage.withdraw_fraction
                for key in sbml_species.values():
                    runner[key] = float(runner[key]) * fraction
                volume *= fraction
                runner["vessel"] = volume
            elif stage.operation == "transfer":
                for key in sbml_species.values():
                    runner[key] = float(runner[key]) * stage.recovery_fraction
                volume *= stage.recovery_fraction
                runner["vessel"] = volume
            enabled = set(stage.active_reactions) if stage.operation == "reaction" else set()
            for reaction_id, parameter_id in active.items():
                runner[parameter_id] = 1.0 if reaction_id in enabled else 0.0
            duration = stage.duration.as_unit("s")
            if duration > 0:
                runner.simulate(current_time, current_time + duration, max(2, request.points_per_reaction_stage))
                current_time += duration
        cross_final = {key: float(runner[sbml_id]) for key, sbml_id in sbml_species.items()}
        differences = {key: cross_final[key] - native_final[key] for key in ids}
        scale = max(max(native_final.values(), default=0), 1e-12)
        maximum = max((abs(value) for value in differences.values()), default=0.0)
        tolerance = max(1e-8 * scale, 1e-11)
        return {
            "status": "passed" if maximum <= tolerance else "disagreed",
            "engine": "libRoadRunner/CVODE",
            "version": roadrunner.__version__,
            "final_amounts_mol": cross_final,
            "difference_from_native_mol": differences,
            "max_abs_difference_mol": maximum,
            "tolerance_mol": tolerance,
        }
    except Exception as error:
        return {
            "status": "failed",
            "engine": "libRoadRunner/CVODE",
            "version": getattr(roadrunner, "__version__", "unknown"),
            "error": f"{type(error).__name__}: {error}",
        }


def inspect_sbml(sbml: str) -> dict[str, Any]:
    document = libsbml.readSBMLFromString(sbml)
    validation = _validation(document)
    model = document.getModel()
    if model is None:
        return {"validation": validation, "model": None}
    return {
        "validation": validation,
        "model": {
            "id": model.getId(),
            "name": model.getName(),
            "compartments": model.getNumCompartments(),
            "species": [{"id": item.getId(), "name": item.getName()} for item in model.getListOfSpecies()],
            "reactions": [{"id": item.getId(), "name": item.getName()} for item in model.getListOfReactions()],
            "events": [{"id": item.getId(), "name": item.getName()} for item in model.getListOfEvents()],
        },
        "sha256": sha256(sbml.encode("utf-8")).hexdigest(),
    }
