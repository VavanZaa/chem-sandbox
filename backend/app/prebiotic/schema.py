from __future__ import annotations

import re
from typing import Any, Literal

from pydantic import Field, model_validator

from app.operations.schema import QuantityInput, StrictModel
from app.thermochemistry.schema import EvidenceQuantity


SPECIES_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.+-]{0,63}$")
COMPARTMENTS = ("boiler_liquid", "gas_loop", "condensate_trap", "collected_condensate", "wall")
Compartment = Literal["boiler_liquid", "gas_loop", "condensate_trap", "collected_condensate", "wall"]
EvidenceBasis = Literal["computed", "experimental", "literature", "hypothetical"]


class PrebioticSpecies(StrictModel):
    id: str = Field(pattern=SPECIES_ID.pattern)
    name: str = Field(min_length=1, max_length=160)
    composition: dict[str, int] = Field(default_factory=dict, max_length=30)
    charge: int = Field(default=0, ge=-20, le=20)
    smiles: str | None = Field(default=None, min_length=1, max_length=2000)
    phase_behavior: Literal[
        "permanent_gas", "volatile", "condensable", "nonvolatile", "surface", "electron"
    ]
    evaporation_rate: EvidenceQuantity | None = None
    condensation_rate: EvidenceQuantity | None = None
    wall_loss_rate: EvidenceQuantity | None = None
    enantiomer_of: str | None = Field(default=None, pattern=SPECIES_ID.pattern)
    identity_note: str | None = Field(default=None, max_length=500)

    @model_validator(mode="after")
    def validate_species(self) -> "PrebioticSpecies":
        if self.phase_behavior == "electron":
            if self.composition or self.charge != -1 or self.smiles is not None:
                raise ValueError("Electron species require empty composition, charge -1 and no SMILES")
            if any(
                value is not None
                for value in (self.evaporation_rate, self.condensation_rate, self.wall_loss_rate)
            ):
                raise ValueError("Electron species cannot use bulk phase-transfer rates")
            return self
        if not self.composition:
            raise ValueError("Material species require a nonempty elemental composition")
        if any(
            not re.fullmatch(r"[A-Z][a-z]?", element) or count <= 0
            for element, count in self.composition.items()
        ):
            raise ValueError("Composition keys must be element symbols with positive counts")
        for label, quantity in (
            ("evaporation", self.evaporation_rate),
            ("condensation", self.condensation_rate),
            ("wall loss", self.wall_loss_rate),
        ):
            if quantity is not None and quantity.as_unit("1/s") < 0:
                raise ValueError(f"{label} rate cannot be negative")
        if self.phase_behavior == "surface" and (
            self.evaporation_rate is not None or self.condensation_rate is not None
        ):
            raise ValueError("Surface species cannot have bulk evaporation or condensation rates")
        if self.enantiomer_of == self.id:
            raise ValueError("A species cannot be its own enantiomer")
        return self


class ElectronCollisionDatasetRequest(StrictModel):
    id: str = Field(pattern=SPECIES_ID.pattern)
    target_species: str = Field(min_length=1, max_length=120)
    product_species: str | None = Field(default=None, max_length=120)
    kind: Literal["effective", "excitation", "ionization", "attachment"]
    energy_levels_eV: list[float] = Field(min_length=2, max_length=100_000)
    cross_sections_m2: list[float] = Field(min_length=2, max_length=100_000)
    source: str = Field(min_length=3, max_length=1000)
    method: str = Field(min_length=3, max_length=500)
    basis: Literal["experimental", "literature", "computed"]
    coverage: Literal["single_channel", "complete_target_set"] = "single_channel"
    mechanism_origin: Literal["user_declared", "blind_discovery_artifact"] = "user_declared"
    discovery_manifest_sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")

    @model_validator(mode="after")
    def validate_dataset(self) -> "ElectronCollisionDatasetRequest":
        if len(self.energy_levels_eV) != len(self.cross_sections_m2):
            raise ValueError("Energy levels and cross sections must have the same length")
        if self.energy_levels_eV[0] < 0 or any(
            right <= left for left, right in zip(self.energy_levels_eV, self.energy_levels_eV[1:])
        ):
            raise ValueError("Electron energy levels must be nonnegative and strictly increasing")
        if any(value < 0 for value in self.cross_sections_m2):
            raise ValueError("Electron collision cross sections cannot be negative")
        if max(self.cross_sections_m2) <= 0:
            raise ValueError("At least one electron collision cross section must be positive")
        if self.mechanism_origin == "blind_discovery_artifact" and self.discovery_manifest_sha256 is None:
            raise ValueError("Blind-discovery collision data require a frozen manifest hash")
        if self.mechanism_origin == "user_declared" and self.discovery_manifest_sha256 is not None:
            raise ValueError("User-declared collision data cannot claim a blind-discovery manifest")
        return self


class ElectronEnergyDistribution(StrictModel):
    mode: Literal["maxwellian", "tabulated"]
    mean_energy: QuantityInput | None = None
    energy_levels_eV: list[float] | None = Field(default=None, min_length=2, max_length=100_000)
    probability_density_eV_inv: list[float] | None = Field(
        default=None, min_length=2, max_length=100_000
    )
    source: str = Field(min_length=3, max_length=1000)
    method: str = Field(min_length=3, max_length=500)
    basis: EvidenceBasis

    @model_validator(mode="after")
    def validate_distribution(self) -> "ElectronEnergyDistribution":
        if self.mode == "maxwellian":
            if self.mean_energy is None:
                raise ValueError("Maxwellian EEDF requires a mean electron energy")
            if self.mean_energy.as_unit("eV") <= 0:
                raise ValueError("Mean electron energy must be positive")
            if self.energy_levels_eV is not None or self.probability_density_eV_inv is not None:
                raise ValueError("Maxwellian EEDF cannot also contain a tabulated distribution")
            return self
        if self.mean_energy is not None:
            raise ValueError("Tabulated EEDF derives its mean from the supplied distribution")
        if self.energy_levels_eV is None or self.probability_density_eV_inv is None:
            raise ValueError("Tabulated EEDF requires energy and probability-density arrays")
        if len(self.energy_levels_eV) != len(self.probability_density_eV_inv):
            raise ValueError("Tabulated EEDF energy and density arrays must have the same length")
        if self.energy_levels_eV[0] < 0 or any(
            right <= left for left, right in zip(self.energy_levels_eV, self.energy_levels_eV[1:])
        ):
            raise ValueError("EEDF energies must be nonnegative and strictly increasing")
        if any(value < 0 for value in self.probability_density_eV_inv):
            raise ValueError("EEDF probability density cannot be negative")
        if max(self.probability_density_eV_inv) <= 0:
            raise ValueError("EEDF probability density must contain a positive value")
        return self


class ElectronCollisionRateRequest(StrictModel):
    dataset: ElectronCollisionDatasetRequest
    eedf: ElectronEnergyDistribution
    electron_number_density: EvidenceQuantity

    @model_validator(mode="after")
    def validate_density(self) -> "ElectronCollisionRateRequest":
        if self.electron_number_density.as_unit("1/m^3") <= 0:
            raise ValueError("Electron number density must be positive")
        return self


class BoltzmannEEDFRequest(StrictModel):
    """Inputs for a composition-dependent, two-term Boltzmann EEDF solve."""

    solver: Literal["bolos_two_term"] = "bolos_two_term"
    collision_set: Literal[
        "lxcat_morgan_2026-09-22",
        "hybrid_morgan_phelps_2026-09-23",
        "hybrid_morgan_phelps_tabata_2026-09-23",
        "hybrid_morgan_community_phelps_tabata_2026-09-23",
        "hybrid_ist_h2_community_ch4_phelps_tabata_2026-09-23",
        "hybrid_ist_h2_co2_o2_phelps_tabata_2026-09-26",
    ] = "lxcat_morgan_2026-09-22"
    mole_fractions: dict[str, float] = Field(min_length=1, max_length=7)
    gas_temperature_K: float = Field(default=400.0, ge=50.0, le=5000.0)
    reduced_electric_field_Td: float = Field(gt=0.0, le=2000.0)
    energy_max_eV: float = Field(default=60.0, ge=20.0, le=500.0)
    energy_cells: int = Field(default=400, ge=100, le=4000)
    initial_electron_temperature_eV: float = Field(default=2.0, gt=0.0, le=100.0)
    maximum_iterations: int = Field(default=200, ge=10, le=2000)
    relative_tolerance: float = Field(default=1e-6, ge=1e-10, le=1e-2)
    acknowledge_two_term_limit: Literal[True]

    @model_validator(mode="after")
    def validate_mixture(self) -> "BoltzmannEEDFRequest":
        supported = {"H2", "H2O", "CH4", "NH3"}
        if self.collision_set in {
            "hybrid_morgan_phelps_2026-09-23",
            "hybrid_morgan_phelps_tabata_2026-09-23",
            "hybrid_morgan_community_phelps_tabata_2026-09-23",
            "hybrid_ist_h2_community_ch4_phelps_tabata_2026-09-23",
            "hybrid_ist_h2_co2_o2_phelps_tabata_2026-09-26",
        }:
            supported.add("N2")
        if self.collision_set == "hybrid_ist_h2_co2_o2_phelps_tabata_2026-09-26":
            supported.update({"CO2", "O2"})
        unknown = set(self.mole_fractions) - supported
        if unknown:
            raise ValueError(
                f"The selected collision set only covers {', '.join(sorted(supported))}; "
                f"unsupported targets: {sorted(unknown)}"
            )
        if any(value < 0 for value in self.mole_fractions.values()):
            raise ValueError("Mole fractions cannot be negative")
        if sum(self.mole_fractions.values()) <= 0:
            raise ValueError("At least one mole fraction must be positive")
        if (
            self.collision_set in {
                "hybrid_morgan_phelps_tabata_2026-09-23",
                "hybrid_morgan_community_phelps_tabata_2026-09-23",
                "hybrid_ist_h2_community_ch4_phelps_tabata_2026-09-23",
                "hybrid_ist_h2_co2_o2_phelps_tabata_2026-09-26",
            }
            and self.energy_cells < 200
        ):
            raise ValueError(
                "The reviewed N2 dissociation outcome requires at least 200 energy cells"
            )
        if (
            self.collision_set == "hybrid_ist_h2_community_ch4_phelps_tabata_2026-09-23"
            and not 1 <= self.reduced_electric_field_Td <= 200
        ):
            raise ValueError(
                "The IST-Lisbon H2 set without rotational channels is documented for 1-200 Td"
            )
        if (
            self.collision_set == "hybrid_ist_h2_co2_o2_phelps_tabata_2026-09-26"
            and not 10 <= self.reduced_electric_field_Td <= 200
        ):
            raise ValueError(
                "The selected IST-Lisbon H2/O2 sets without rotational channels "
                "are jointly documented only for 10-200 Td"
            )
        return self


class BoltzmannDischargeProgram(StrictModel):
    """Self-consistent EEDF controls whose mixture is supplied by the apparatus."""

    solver: Literal["bolos_two_term"] = "bolos_two_term"
    collision_set: Literal[
        "lxcat_morgan_2026-09-22",
        "hybrid_morgan_phelps_2026-09-23",
        "hybrid_morgan_phelps_tabata_2026-09-23",
        "hybrid_morgan_community_phelps_tabata_2026-09-23",
        "hybrid_ist_h2_community_ch4_phelps_tabata_2026-09-23",
        "hybrid_ist_h2_co2_o2_phelps_tabata_2026-09-26",
    ] = "lxcat_morgan_2026-09-22"
    reduced_electric_field_Td: float = Field(gt=0.0, le=2000.0)
    electron_number_density: EvidenceQuantity
    energy_max_eV: float = Field(default=60.0, ge=20.0, le=500.0)
    energy_cells: int = Field(default=400, ge=100, le=4000)
    initial_electron_temperature_eV: float = Field(default=2.0, gt=0.0, le=100.0)
    maximum_iterations: int = Field(default=200, ge=10, le=2000)
    relative_tolerance: float = Field(default=1e-6, ge=1e-10, le=1e-2)
    acknowledge_two_term_limit: Literal[True]

    @model_validator(mode="after")
    def validate_density(self) -> "BoltzmannDischargeProgram":
        if self.electron_number_density.as_unit("1/m^3") <= 0:
            raise ValueError("Electron number density must be positive")
        if (
            self.collision_set in {
                "hybrid_morgan_phelps_tabata_2026-09-23",
                "hybrid_morgan_community_phelps_tabata_2026-09-23",
                "hybrid_ist_h2_community_ch4_phelps_tabata_2026-09-23",
                "hybrid_ist_h2_co2_o2_phelps_tabata_2026-09-26",
            }
            and self.energy_cells < 200
        ):
            raise ValueError(
                "The reviewed N2 dissociation outcome requires at least 200 energy cells"
            )
        if (
            self.collision_set == "hybrid_ist_h2_community_ch4_phelps_tabata_2026-09-23"
            and not 1 <= self.reduced_electric_field_Td <= 200
        ):
            raise ValueError(
                "The IST-Lisbon H2 set without rotational channels is documented for 1-200 Td"
            )
        if (
            self.collision_set == "hybrid_ist_h2_co2_o2_phelps_tabata_2026-09-26"
            and not 10 <= self.reduced_electric_field_Td <= 200
        ):
            raise ValueError(
                "The selected IST-Lisbon H2/O2 sets without rotational channels "
                "are jointly documented only for 10-200 Td"
            )
        return self


class PrebioticRateLaw(StrictModel):
    kind: Literal["mass_action", "electron_impact_reduced", "electron_collision"]
    rate_constant: EvidenceQuantity | None = None
    collision_dataset: ElectronCollisionDatasetRequest | None = None
    target_species_id: str | None = Field(default=None, pattern=SPECIES_ID.pattern)
    electron_species_id: str | None = Field(default=None, pattern=SPECIES_ID.pattern)
    energy_cost: EvidenceQuantity | None = None

    @model_validator(mode="after")
    def validate_energy(self) -> "PrebioticRateLaw":
        collision_fields = (
            self.collision_dataset,
            self.target_species_id,
            self.electron_species_id,
        )
        if self.kind == "electron_collision":
            if self.rate_constant is not None or any(value is None for value in collision_fields):
                raise ValueError(
                    "Electron-collision channels require a dataset, target and electron species, not a fitted rate constant"
                )
        elif self.rate_constant is None or any(value is not None for value in collision_fields):
            raise ValueError("Declared-rate channels require only a rate_constant")
        if self.kind == "mass_action" and self.energy_cost is not None:
            raise ValueError("energy_cost is restricted to electron-impact channels")
        if self.energy_cost is not None and self.energy_cost.as_unit("J/mol") < 0:
            raise ValueError("Electron-impact energy cost cannot be negative")
        return self


class PrebioticReaction(StrictModel):
    id: str = Field(pattern=SPECIES_ID.pattern)
    name: str = Field(min_length=1, max_length=200)
    phase: Literal["gas", "liquid"]
    stoichiometry: dict[str, int] = Field(min_length=2, max_length=50)
    rate_law: PrebioticRateLaw

    @model_validator(mode="after")
    def validate_reaction(self) -> "PrebioticReaction":
        if any(value == 0 or abs(value) > 20 for value in self.stoichiometry.values()):
            raise ValueError("Reaction coefficients must be nonzero integers with magnitude at most 20")
        if not any(value < 0 for value in self.stoichiometry.values()) or not any(
            value > 0 for value in self.stoichiometry.values()
        ):
            raise ValueError("A reaction requires reactants and products")
        if self.rate_law.kind in {"electron_impact_reduced", "electron_collision"} and self.phase != "gas":
            raise ValueError("Electron-impact channels are restricted to the gas loop")
        return self


class CompartmentAmount(StrictModel):
    species_id: str = Field(pattern=SPECIES_ID.pattern)
    compartment: Compartment
    amount: QuantityInput


class ApparatusAddition(CompartmentAmount):
    pass


class CycleIntervention(StrictModel):
    cycle: int = Field(ge=1, le=100_000)
    additions: list[ApparatusAddition] = Field(default_factory=list, max_length=100)
    sample_compartment: Compartment | None = None
    sample_fraction: float = Field(default=0, ge=0, lt=1)
    note: str | None = Field(default=None, max_length=1000)

    @model_validator(mode="after")
    def validate_sampling(self) -> "CycleIntervention":
        if (self.sample_compartment is None) != (self.sample_fraction == 0):
            raise ValueError("Sampling requires both sample_compartment and a positive sample_fraction")
        if not self.additions and self.sample_fraction == 0:
            raise ValueError("An intervention requires an addition or sample")
        return self


class SparkProgram(StrictModel):
    mode: Literal["off", "reduced_order", "cross_section_eedf", "boltzmann_eedf"] = "off"
    voltage: QuantityInput | None = None
    current: QuantityInput | None = None
    pulse_width: QuantityInput | None = None
    frequency: QuantityInput | None = None
    coupling_efficiency: float | None = Field(default=None, gt=0, le=1)
    mean_electron_energy: QuantityInput | None = None
    eedf: ElectronEnergyDistribution | None = None
    boltzmann: BoltzmannDischargeProgram | None = None
    calorimetric_deposited_energy_per_cycle: EvidenceQuantity | None = None
    source: str | None = Field(default=None, max_length=1000)
    basis: EvidenceBasis | None = None

    @model_validator(mode="after")
    def validate_program(self) -> "SparkProgram":
        electrical = (
            self.voltage,
            self.current,
            self.pulse_width,
            self.frequency,
            self.coupling_efficiency,
            self.mean_electron_energy,
            self.eedf,
            self.boltzmann,
            self.calorimetric_deposited_energy_per_cycle,
            self.source,
            self.basis,
        )
        if self.mode == "off":
            if any(value is not None for value in electrical):
                raise ValueError("Spark parameters cannot be supplied when the spark is off")
            return self
        common = (
            self.voltage,
            self.current,
            self.pulse_width,
            self.frequency,
            self.coupling_efficiency,
            self.source,
            self.basis,
        )
        if any(value is None for value in common):
            raise ValueError("Active spark modes require complete electrical conditions and provenance")
        if self.mode == "reduced_order":
            if self.mean_electron_energy is None or self.eedf is not None or self.boltzmann is not None:
                raise ValueError(
                    "Reduced-order spark mode requires mean_electron_energy and no EEDF solver"
                )
            electron_energy = self.mean_electron_energy.as_unit("eV")
        elif self.mode == "cross_section_eedf":
            if self.eedf is None or self.mean_electron_energy is not None or self.boltzmann is not None:
                raise ValueError("Cross-section spark mode requires an EEDF and no standalone mean energy")
            electron_energy = 1.0
        else:
            if self.boltzmann is None or self.mean_electron_energy is not None or self.eedf is not None:
                raise ValueError(
                    "Boltzmann spark mode requires solver controls and derives the mixture from the gas loop"
                )
            electron_energy = self.boltzmann.initial_electron_temperature_eV
        voltage = self.voltage.as_unit("V")
        current = self.current.as_unit("A")
        width = self.pulse_width.as_unit("s")
        frequency = self.frequency.as_unit("Hz")
        if min(voltage, current, width, frequency, electron_energy) <= 0:
            raise ValueError("Spark voltage, current, pulse width, frequency and electron energy must be positive")
        if width * frequency > 1 + 1e-12:
            raise ValueError("Spark pulse width times frequency cannot exceed a duty fraction of one")
        if self.calorimetric_deposited_energy_per_cycle is not None:
            if self.calorimetric_deposited_energy_per_cycle.as_unit("J") <= 0:
                raise ValueError("Calorimetric deposited energy per cycle must be positive")
            if self.calorimetric_deposited_energy_per_cycle.basis == "computed":
                raise ValueError(
                    "Calorimetric deposited energy must be experimental, literature, or hypothetical evidence"
                )
        return self


class RecirculatingApparatus(StrictModel):
    cycles: int = Field(ge=1, le=100_000)
    operation_mode: Literal["staged", "continuous_overlap", "closed_stirred_flask"] = "staged"
    integration_substeps_per_cycle: int = Field(default=1, ge=1, le=48)
    primary_source_forcing: Literal["impulse", "time_averaged", "inventory_limited"] = "impulse"
    aqueous_acid_base_model: Literal["off", "hcn_nh3_water_equilibrium"] = "off"
    boiler_liquid_volume: QuantityInput
    gas_loop_volume: QuantityInput
    discharge_zone_volume: EvidenceQuantity | None = None
    condensate_trap_volume: QuantityInput
    gas_circulation_rate: QuantityInput
    pressure: QuantityInput
    boiler_temperature: QuantityInput
    discharge_temperature: QuantityInput
    condenser_temperature: QuantityInput
    boil_duration: QuantityInput
    discharge_duration: QuantityInput
    condense_duration: QuantityInput
    condensate_return_fraction: float = Field(default=1, ge=0, le=1)
    condensate_nonwater_return_fraction: float = Field(default=1, ge=0, le=1)
    condenser_phase_transfer_model: Literal[
        "step_equilibrium", "flow_limited_capture", "coupled_flow_transport"
    ] = "step_equilibrium"
    condenser_contact_water_holdup: EvidenceQuantity | None = None
    closed_flask_gas_water_transfer_model: Literal[
        "step_equilibrium", "finite_rate_contact"
    ] = "step_equilibrium"
    closed_flask_liquid_side_kla: EvidenceQuantity | None = None
    closed_flask_finite_contact_solver: Literal[
        "post_step_analytic", "coupled_bdf"
    ] = "post_step_analytic"
    closed_flask_formaldehyde_hydration: Literal[
        "off", "winkelman_2002_source_bounded",
        "winkelman_2002_pH4_to_5_extrapolation_diagnostic",
    ] = "off"

    @model_validator(mode="after")
    def validate_apparatus(self) -> "RecirculatingApparatus":
        positive = {
            "boiler liquid volume": self.boiler_liquid_volume.as_unit("L"),
            "gas loop volume": self.gas_loop_volume.as_unit("L"),
            "condensate trap volume": self.condensate_trap_volume.as_unit("L"),
            "pressure": self.pressure.as_unit("Pa"),
            "boiler temperature": self.boiler_temperature.as_unit("K"),
            "discharge temperature": self.discharge_temperature.as_unit("K"),
            "condenser temperature": self.condenser_temperature.as_unit("K"),
        }
        invalid = [name for name, value in positive.items() if value <= 0]
        if invalid:
            raise ValueError(f"Apparatus quantities must be positive: {', '.join(invalid)}")
        if self.gas_circulation_rate.as_unit("L/s") < 0:
            raise ValueError("Gas circulation rate cannot be negative")
        if self.discharge_zone_volume is not None:
            discharge_zone_volume_L = self.discharge_zone_volume.as_unit("L")
            gas_loop_volume_L = self.gas_loop_volume.as_unit("L")
            if discharge_zone_volume_L <= 0:
                raise ValueError("Discharge zone volume must be positive")
            if discharge_zone_volume_L > gas_loop_volume_L * (1 + 1e-12):
                raise ValueError("Discharge zone volume cannot exceed the gas loop volume")
        durations = (
            self.boil_duration.as_unit("s"),
            self.discharge_duration.as_unit("s"),
            self.condense_duration.as_unit("s"),
        )
        if any(value < 0 for value in durations) or sum(durations) <= 0:
            raise ValueError("Stage durations cannot be negative and at least one must be positive")
        if self.operation_mode == "continuous_overlap" and self.integration_substeps_per_cycle < 2:
            raise ValueError(
                "Continuous-overlap operation requires at least two integration substeps per cycle"
            )
        if self.operation_mode == "closed_stirred_flask":
            if self.integration_substeps_per_cycle < 2:
                raise ValueError("Closed stirred flask requires at least two integration substeps per cycle")
            if any(value != 0 for value in (durations[0], durations[2], self.gas_circulation_rate.as_unit("L/s"))):
                raise ValueError("Closed stirred flask cannot boil, condense, or circulate gas")
            temperatures = (
                self.boiler_temperature.as_unit("K"),
                self.discharge_temperature.as_unit("K"),
                self.condenser_temperature.as_unit("K"),
            )
            if max(temperatures) - min(temperatures) > 1e-6:
                raise ValueError("Closed stirred flask requires one isothermal gas/water temperature")
            if self.condenser_phase_transfer_model != "step_equilibrium":
                raise ValueError("Closed stirred flask uses only step-equilibrium gas/water contact")
            if self.closed_flask_gas_water_transfer_model == "finite_rate_contact":
                if self.closed_flask_liquid_side_kla is None:
                    raise ValueError("Finite-rate closed-flask contact requires declared liquid-side kLa")
                if self.closed_flask_liquid_side_kla.as_unit("1/s") < 0:
                    raise ValueError("Closed-flask liquid-side kLa cannot be negative")
            elif self.closed_flask_liquid_side_kla is not None:
                raise ValueError("Closed-flask liquid-side kLa requires finite-rate contact")
            if (
                self.closed_flask_gas_water_transfer_model != "finite_rate_contact"
                and self.closed_flask_finite_contact_solver != "post_step_analytic"
            ):
                raise ValueError("Coupled closed-flask contact requires finite-rate contact")
            if (
                self.closed_flask_formaldehyde_hydration != "off"
                and (
                    self.closed_flask_finite_contact_solver != "coupled_bdf"
                    or self.aqueous_acid_base_model != "hcn_nh3_water_equilibrium"
                )
            ):
                raise ValueError(
                    "Source-bounded formaldehyde hydration requires coupled BDF contact "
                    "and aqueous acid-base pH calculation"
                )
        elif (
            self.closed_flask_gas_water_transfer_model != "step_equilibrium"
            or self.closed_flask_liquid_side_kla is not None
            or self.closed_flask_finite_contact_solver != "post_step_analytic"
            or self.closed_flask_formaldehyde_hydration != "off"
        ):
            raise ValueError("Closed-flask gas/water transfer settings require closed-stirred-flask mode")
        if self.condenser_phase_transfer_model in {"flow_limited_capture", "coupled_flow_transport"}:
            if self.condenser_contact_water_holdup is None:
                raise ValueError("Flow-limited condenser capture requires a declared contact-water holdup")
            holdup_L = self.condenser_contact_water_holdup.as_unit("L")
            if holdup_L <= 0 or holdup_L > self.condensate_trap_volume.as_unit("L") * (1 + 1e-12):
                raise ValueError("Condenser contact-water holdup must be positive and fit inside the trap")
        return self


class PrebioticApparatusRequest(StrictModel):
    id: str = Field(min_length=1, max_length=100)
    species: list[PrebioticSpecies] = Field(min_length=1, max_length=500)
    reactions: list[PrebioticReaction] = Field(default_factory=list, max_length=2000)
    initial_inventory: list[CompartmentAmount] = Field(min_length=1, max_length=2000)
    apparatus: RecirculatingApparatus
    spark: SparkProgram = Field(default_factory=SparkProgram)
    interventions: list[CycleIntervention] = Field(default_factory=list, max_length=1000)
    mechanism_origin: Literal["user_declared", "blind_discovery_artifact"] = "user_declared"
    discovery_manifest_sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")
    knowledge_policy: Literal["blind_physics_only"]
    acknowledge_reduced_order_limit: Literal[True]

    @model_validator(mode="after")
    def validate_references(self) -> "PrebioticApparatusRequest":
        species_ids = [item.id for item in self.species]
        reaction_ids = [item.id for item in self.reactions]
        if len(species_ids) != len(set(species_ids)):
            raise ValueError("Prebiotic species ids must be unique")
        if len(reaction_ids) != len(set(reaction_ids)):
            raise ValueError("Prebiotic reaction ids must be unique")
        known = set(species_ids)
        species_by_id = {item.id: item for item in self.species}
        for species in self.species:
            if species.enantiomer_of is not None and species.enantiomer_of not in known:
                raise ValueError(f"Species {species.id} references an unknown enantiomer")
            if species.enantiomer_of is not None:
                partner = species_by_id[species.enantiomer_of]
                if partner.enantiomer_of != species.id:
                    raise ValueError("Enantiomer declarations must be reciprocal")
                if partner.composition != species.composition or partner.charge != species.charge:
                    raise ValueError("Declared enantiomers must have identical composition and charge")
        for reaction in self.reactions:
            unknown = set(reaction.stoichiometry) - known
            if unknown:
                raise ValueError(f"Reaction {reaction.id} references unknown species: {sorted(unknown)}")
            if reaction.rate_law.kind == "electron_collision":
                target = reaction.rate_law.target_species_id
                electron = reaction.rate_law.electron_species_id
                if target not in known or electron not in known:
                    raise ValueError("Electron-collision rate law references unknown target or electron species")
                if target == electron:
                    raise ValueError("Electron-collision target and electron species must be distinct")
                if species_by_id[electron].phase_behavior != "electron":
                    raise ValueError("electron_species_id must reference an electron species")
                if reaction.stoichiometry.get(target, 0) != -1:
                    raise ValueError("Electron-collision channels must consume exactly one target particle")
                dataset = reaction.rate_law.collision_dataset
                if dataset.kind == "effective":
                    raise ValueError("Effective/elastic collision curves cannot define a chemical channel")
                if dataset.mechanism_origin != self.mechanism_origin:
                    raise ValueError("Collision dataset origin must match the apparatus mechanism origin")
                if dataset.discovery_manifest_sha256 != self.discovery_manifest_sha256:
                    raise ValueError("Collision dataset manifest must match the apparatus manifest")
        for amount in self.initial_inventory:
            if amount.species_id not in known:
                raise ValueError("Initial inventory references an unknown species")
            if species_by_id[amount.species_id].phase_behavior == "electron" and amount.compartment != "gas_loop":
                raise ValueError("Electron inventory is restricted to the gas loop")
        for intervention in self.interventions:
            if intervention.cycle > self.apparatus.cycles:
                raise ValueError("Intervention cycle exceeds the apparatus cycle count")
            if any(addition.species_id not in known for addition in intervention.additions):
                raise ValueError("Intervention references an unknown species")
            if any(
                species_by_id[addition.species_id].phase_behavior == "electron"
                and addition.compartment != "gas_loop"
                for addition in intervention.additions
            ):
                raise ValueError("Electron additions are restricted to the gas loop")
        if self.mechanism_origin == "blind_discovery_artifact" and self.discovery_manifest_sha256 is None:
            raise ValueError("Blind-discovery mechanisms require a frozen discovery manifest hash")
        if self.mechanism_origin == "user_declared" and self.discovery_manifest_sha256 is not None:
            raise ValueError("User-declared mechanisms cannot claim a blind-discovery manifest")
        if (
            self.apparatus.operation_mode == "closed_stirred_flask"
            and self.spark.mode != "off"
            and self.apparatus.discharge_zone_volume is None
        ):
            raise ValueError("An active closed-flask spark requires a declared discharge-zone volume")
        reduced = any(
            reaction.rate_law.kind == "electron_impact_reduced" for reaction in self.reactions
        )
        collision = any(
            reaction.rate_law.kind == "electron_collision" for reaction in self.reactions
        )
        if reduced and self.spark.mode in {"cross_section_eedf", "boltzmann_eedf"}:
            raise ValueError("Reduced pseudo-rate channels cannot run under an EEDF solver mode")
        if collision and self.spark.mode not in {"cross_section_eedf", "off"}:
            raise ValueError(
                "Explicit cross-section reactions require cross_section_eedf mode or an off control; "
                "automatic Boltzmann source terms are reported separately"
            )
        if (
            self.apparatus.condenser_phase_transfer_model == "coupled_flow_transport"
            and self.spark.mode == "boltzmann_eedf"
            and self.apparatus.primary_source_forcing != "time_averaged"
        ):
            raise ValueError("Coupled flow transport requires time-averaged plasma source forcing")
        if self.apparatus.closed_flask_finite_contact_solver == "coupled_bdf":
            if self.spark.mode != "boltzmann_eedf":
                raise ValueError("Coupled closed-flask contact requires Boltzmann EEDF mode")
            if self.apparatus.primary_source_forcing not in {"time_averaged", "inventory_limited"}:
                raise ValueError("Coupled closed-flask contact requires distributed plasma source forcing")
        return self


class FrozenCandidateRankingRequest(StrictModel):
    discovery: dict[str, Any]
    discovery_artifact_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    engines: list[Literal["mace_mp", "ani2x", "aimnet2", "aimnet2_rxn"]] = Field(
        default_factory=lambda: ["mace_mp", "ani2x", "aimnet2", "aimnet2_rxn"],
        min_length=1,
        max_length=4,
    )
    device: Literal["auto", "cuda", "cpu"] = "auto"
    optimization_steps: int = Field(default=60, ge=5, le=500)
    force_threshold_eV_A: float = Field(default=0.15, gt=0, le=1)
    acknowledge_vacuum_energy_limit: Literal[True]

    @model_validator(mode="after")
    def unique_engines(self) -> "FrozenCandidateRankingRequest":
        if len(self.engines) != len(set(self.engines)):
            raise ValueError("Candidate-ranking engines must be unique")
        return self
