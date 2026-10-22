"""Independent balanced-reaction ceiling; NEVER used to steer discovery."""
from .schema import ExperimentDefinition, SUBSTANCE_LIBRARY


def theoretical_yield(experiment: ExperimentDefinition) -> dict:
    amounts = {}
    for item in experiment.substances:
        material = SUBSTANCE_LIBRARY[item.id]
        mass = item.amount * material.get("density_g_mL", 1) if item.unit == "mL" else item.amount
        amounts[item.id.value] = item.amount if item.unit == "mol" else mass / material["molar_mass_g_mol"]
    bicarbonate, acid = amounts["sodium_bicarbonate"], amounts["acetic_acid"]
    extent = min(bicarbonate, acid)
    limiting = "stoichiometric mixture" if abs(bicarbonate - acid) <= 1e-9 * max(bicarbonate, acid) else (
        "Acetic acid" if acid < bicarbonate else "Baking soda")
    return {
        "basis": "Balanced-reaction theoretical ceiling; independent of atomistic search",
        "equation": "NaHCO₃ + CH₃COOH → CH₃COONa + CO₂ + H₂O",
        "limiting_reagent": limiting, "input_moles": amounts,
        "extent_mol": extent, "co2_mol": extent, "co2_mass_g": extent * 44.0095,
        "co2_ideal_gas_L": extent * 0.082057366080960 * experiment.environment.temperature_K / experiment.environment.pressure_atm,
        "sodium_acetate_mol": extent, "water_formed_mol": extent,
        "excess_bicarbonate_mol": max(0, bicarbonate - extent), "excess_acid_mol": max(0, acid - extent),
        "temperature_K": experiment.environment.temperature_K, "pressure_atm": experiment.environment.pressure_atm,
        "actual_yield": None,
        "assumptions": ["Amounts refer to pure reagents, not the total mass of vinegar or another solution.",
                        "Assumes complete conversion of the limiting reagent via this one balanced reaction; competing products would share that inventory.",
                        "Gas volume assumes ideal dry CO₂ at the entered temperature and pressure, with no dissolved gas or collection loss.",
                        "Actual yield needs validated reaction rates/equilibria, concentrations, time, gas exchange and experimental calibration."]}
