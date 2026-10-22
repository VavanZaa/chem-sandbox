# Phase 0–2 backend API

## Phase 0: blind discovery and evaluation

- `GET /api/scientific-policy` returns the accepted blind inputs and post-run-only fields.
- Experiment definitions reject undeclared fields at every nested level.
- Completed discovery results include `blind_manifest` and `provenance`.
- `backend/evaluation/reference_only` contains benchmark answers and known endpoint recipes. No runtime package imports it.
- Offline evaluation creates `prediction.frozen.json`, verifies its SHA-256, then creates `evaluation.unlocked.json` with a separate reference hash.

The known-endpoint thermochemistry benchmark is not a blind product-discovery benchmark. Its reaction endpoints are evaluator inputs and never enter candidate generation.

## Phase 1: operation validation

`POST /api/operations/validate` accepts an `operation` and balanced `mechanism`.

Supported operation data:

- batch, semibatch, CSTR, and PFR;
- gas, ideal-solution, aqueous-electrolyte, liquid, solid, and surface phase declarations;
- initial moles, timed molar feeds, volume, temperature, pressure, duration/residence time;
- adiabatic, isothermal, or overall-UA heat transfer;
- ideal/stirred/axial-dispersion mixing declarations;
- equilibrium or kLa gas transfer;
- homogeneous catalyst loading or surface active-site amount.

Pint validates dimensions and normalizes accepted values to SI. Every elementary reaction receives element, charge, and surface-site residuals. `ready_for_solver` remains false for an unbalanced mechanism. This endpoint does not run kinetics.

## Phase 2: thermochemistry

- `POST /api/thermochemistry/harmonic-corrections` calculates vibrational harmonic ZPE, thermal energy, and entropy components from sourced frequencies. It clearly excludes translational, rotational, electronic, solvation, and standard-state terms.
- `POST /api/thermochemistry/ensemble` certifies stationary points and calculates conformer populations/free energy.
- `POST /api/thermochemistry/reaction` calculates a balanced reaction Gibbs energy.
- `POST /api/thermochemistry/activation` calculates activation Gibbs energy only for a verified transition-state ensemble relative to a matching reactant reference ensemble.

Certification requires converged optimized geometries, force tolerance, the correct number of modes, no significant imaginary mode for a minimum, exactly one for a transition state, reaction-coordinate evidence for a transition state, non-hypothetical sourced corrections, consistent composition/charge/spin/method/phase/standard state, and a declared complete conformer search. Solution states additionally require an explicit solvation model and solvation free energy.

When evidence is incomplete, the response contains `problems` and leaves the certified Gibbs field null. Conditional submitted-ensemble values are kept separately and never relabelled as validated reaction yields.
