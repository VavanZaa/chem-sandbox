"""Reviewed cold-water cyanide equilibrium is diagnostic, not a yield model."""

from copy import deepcopy

import pytest
import numpy as np

from app.prebiotic import aqueous_speciation as speciation
from app.prebiotic.phase_partition import equilibrate_gas_water


def test_iapws_water_formula_matches_published_verification_point() -> None:
    data = speciation.load_acid_base_subset()
    assert speciation._water_pkw(300.0, 1.0, data) == pytest.approx(
        data["water"]["verification_pKw_300K_1g_cm3"], abs=1e-6
    )
    assert speciation._hcn_pka(298.15, data) == pytest.approx(9.21)
    assert speciation._ammonium_pka(298.15, data) == pytest.approx(9.245)
    assert speciation._carbonate_pkas(298.15, data) == pytest.approx((6.351, 10.329))


def test_cold_equilibrium_closes_charge_and_conserves_cyanide() -> None:
    report = speciation.speciate_hcn_ammonia_water(
        hcn_mol=1e-6, ammonia_mol=1e-6, water_mol=55.5, temperature_K=298.15
    )
    species = report["equilibrated_species_mol"]
    assert species["HCN"] + species["CN-"] == pytest.approx(1e-6)
    assert species["NH3"] + species["NH4+"] == pytest.approx(1e-6)
    assert abs(report["charge_balance_residual_mol"]) < 1e-15
    assert abs(report["hydrogen_balance_residual_mol"]) < 1e-15
    assert 7 < report["pH_ideal_aqueous_model"] < 9
    assert report["prediction_inventory_mutated"] is False


def test_ammonia_changes_cyanide_ion_fraction_without_external_pH() -> None:
    without = speciation.speciate_hcn_ammonia_water(
        hcn_mol=1e-6, ammonia_mol=0.0, water_mol=55.5, temperature_K=280.0
    )
    with_ammonia = speciation.speciate_hcn_ammonia_water(
        hcn_mol=1e-6, ammonia_mol=1e-3, water_mol=55.5, temperature_K=280.0
    )
    assert without["cyanide_ion_fraction"] < with_ammonia["cyanide_ion_fraction"]
    assert with_ammonia["pH_ideal_aqueous_model"] > without["pH_ideal_aqueous_model"]


def test_dissolved_carbon_dioxide_buffers_ammonia_and_changes_cyanide_fate() -> None:
    conditions = {
        "hcn_mol": 1e-6, "ammonia_mol": 1e-3,
        "water_mol": 55.5, "temperature_K": 298.15,
    }
    without_carbonate = speciation.speciate_hcn_ammonia_water(**conditions)
    with_carbonate = speciation.speciate_hcn_ammonia_water(
        **conditions, carbon_dioxide_mol=1e-3
    )
    amounts = with_carbonate["equilibrated_species_mol"]
    assert with_carbonate["pH_ideal_aqueous_model"] < without_carbonate["pH_ideal_aqueous_model"]
    assert with_carbonate["cyanide_ion_fraction"] < without_carbonate["cyanide_ion_fraction"]
    assert amounts["CO2"] + amounts["HCO3-"] + amounts["CO3--"] == pytest.approx(1e-3)
    assert amounts["HCO3-"] > amounts["CO3--"] > 0
    assert amounts["HCN"] + amounts["CN-"] == pytest.approx(1e-6)
    assert abs(with_carbonate["charge_balance_residual_mol"]) < 1e-14


def test_postprocess_does_not_mutate_completed_inventory() -> None:
    completed = {
        "operating_conditions": {"boiler_temperature_K": 373.15, "condenser_temperature_K": 280.0},
        "final_compartments_mol": {
            "boiler_liquid": {"HCN": 3e-7, "NH3": 0.0, "H2O": 55.5},
            "condensate_trap": {"HCN": 1e-7, "NH3": 1e-6, "H2O": 10.0},
            "collected_condensate": {"HCN": 2e-7, "NH3": 0.0, "H2O": 20.0},
        },
    }
    original = deepcopy(completed)
    report = speciation.evaluate_aqueous_cyanide_speciation(completed)
    assert completed == original
    assert report["prediction_inventory_mutated"] is False
    assert report["boiler_liquid_speciated"] is True
    assert report["compartments"]["boiler_liquid"]["cyanide_total_mol"] == 3e-7
    assert report["compartments"]["condensate_trap"]["cyanide_total_mol"] == 1e-7
    assert report["compartments"]["collected_condensate"]["cyanide_total_mol"] == 2e-7


def test_speciation_rejects_unsupported_or_invalid_conditions() -> None:
    with pytest.raises(ValueError, match="evidence range"):
        speciation.speciate_hcn_ammonia_water(
            hcn_mol=1e-6, ammonia_mol=0, water_mol=55.5, temperature_K=423.15
        )
    with pytest.raises(ValueError, match="nonnegative"):
        speciation.speciate_hcn_ammonia_water(
            hcn_mol=-1e-6, ammonia_mol=0, water_mol=55.5, temperature_K=298.15
        )


def test_unusable_diagnostic_inventory_does_not_break_primary_result() -> None:
    completed = {
        "operating_conditions": {"boiler_temperature_K": 373.15, "condenser_temperature_K": 280.0},
        "final_compartments_mol": {
            "boiler_liquid": {"HCN": -1e-3, "H2O": 1.0},
            "condensate_trap": {},
            "collected_condensate": {},
        },
    }
    report = speciation.evaluate_aqueous_cyanide_speciation(completed)
    assert report["compartments"]["boiler_liquid"]["status"] == "speciation_unavailable_for_inventory"
    assert report["prediction_inventory_mutated"] is False


def test_explicit_ion_update_is_idempotent_and_closes_atom_charge_ledger() -> None:
    species = (
        ("HCN", {"H": 1, "C": 1, "N": 1}),
        ("NH3", {"N": 1, "H": 3}),
        ("H2O", {"H": 2, "O": 1}),
    )
    prepared = {
        "ids": [name for name, _ in species],
        "index": {name: position for position, (name, _) in enumerate(species)},
        "identities": [
            {"id": name, "composition": composition, "charge": 0}
            for name, composition in species
        ],
        "transfer": [{} for _ in species],
        "stoichiometry": np.zeros((len(species), 0)),
    }
    speciation.augment_prepared_with_aqueous_ions(prepared)
    initial = np.zeros(len(prepared["ids"]))
    initial[prepared["index"]["HCN"]] = 1e-6
    initial[prepared["index"]["NH3"]] = 1e-3
    initial[prepared["index"]["H2O"]] = 55.5
    equilibrated, report = speciation.equilibrate_aqueous_ion_inventory(
        initial, prepared, 298.15
    )
    repeated, _ = speciation.equilibrate_aqueous_ion_inventory(
        equilibrated, prepared, 298.15
    )
    assert report["ledger_maximum_residual_mol"] < 1e-12
    assert np.max(np.abs(prepared["atomic"] @ (equilibrated - initial))) < 1e-12
    assert np.max(np.abs(repeated - equilibrated)) < 1e-10
    assert equilibrated[prepared["index"]["CN-"]] > 0
    assert equilibrated[prepared["index"]["NH4+"]] > 0


def test_carbonate_ion_update_conserves_carbon_oxygen_hydrogen_and_charge() -> None:
    species = (
        ("HCN", {"H": 1, "C": 1, "N": 1}),
        ("NH3", {"N": 1, "H": 3}),
        ("H2O", {"H": 2, "O": 1}),
        ("CO2", {"C": 1, "O": 2}),
    )
    prepared = {
        "ids": [name for name, _ in species],
        "index": {name: position for position, (name, _) in enumerate(species)},
        "identities": [
            {"id": name, "composition": composition, "charge": 0}
            for name, composition in species
        ],
        "transfer": [{} for _ in species],
        "stoichiometry": np.zeros((len(species), 0)),
    }
    speciation.augment_prepared_with_aqueous_ions(prepared)
    initial = np.zeros(len(prepared["ids"]))
    for name, amount in (("HCN", 1e-6), ("NH3", 1e-3), ("H2O", 55.5), ("CO2", 1e-3)):
        initial[prepared["index"][name]] = amount
    equilibrated, report = speciation.equilibrate_aqueous_ion_inventory(
        initial, prepared, 298.15
    )
    repeated, _ = speciation.equilibrate_aqueous_ion_inventory(
        equilibrated, prepared, 298.15
    )
    residuals = prepared["atomic"] @ (equilibrated - initial)
    assert dict(zip(prepared["conserved"], residuals))["charge"] == pytest.approx(0, abs=1e-12)
    assert np.max(np.abs(residuals)) < 1e-12
    assert report["carbonate_total_mol"] == pytest.approx(1e-3)
    assert equilibrated[prepared["index"]["HCO3-"]] > 0
    assert equilibrated[prepared["index"]["CO3--"]] > 0
    assert np.max(np.abs(repeated - equilibrated)) < 1e-10


def test_co2_nh3_gas_water_transfer_feeds_ions_without_making_atoms_or_charge() -> None:
    species = (
        ("H2O", {"H": 2, "O": 1}),
        ("CO2", {"C": 1, "O": 2}),
        ("NH3", {"N": 1, "H": 3}),
    )
    prepared = {
        "ids": [name for name, _ in species],
        "index": {name: position for position, (name, _) in enumerate(species)},
        "identities": [
            {"id": name, "composition": composition, "charge": 0}
            for name, composition in species
        ],
        "transfer": [{} for _ in species],
        "stoichiometry": np.zeros((len(species), 0)),
    }
    speciation.augment_prepared_with_aqueous_ions(prepared)
    gas = np.zeros(len(prepared["ids"]))
    liquid = np.zeros(len(prepared["ids"]))
    gas[prepared["index"]["CO2"]] = 1e-3
    gas[prepared["index"]["NH3"]] = 1e-3
    liquid[prepared["index"]["H2O"]] = 55.5
    starting_ledger = prepared["atomic"] @ (gas + liquid)
    gas, liquid, transfer = equilibrate_gas_water(
        gas, liquid, prepared, 298.15, 101_325.0, 5.0
    )
    liquid, equilibrium = speciation.equilibrate_aqueous_ion_inventory(
        liquid, prepared, 298.15
    )
    assert transfer["dataset_id"] == "sander_nist_henry_prebiotic_subset_v3"
    assert equilibrium["carbonate_total_mol"] > 0
    assert liquid[prepared["index"]["HCO3-"]] > 0
    assert liquid[prepared["index"]["NH4+"]] > 0
    assert np.max(np.abs(prepared["atomic"] @ (gas + liquid) - starting_ledger)) < 1e-12
