# Phase 3–5 API guide

These endpoints are post-discovery tools. They consume explicit mechanisms or physical evidence and cannot write target products, recipes, or reference outcomes back into blind candidate generation.

## Phase 3: rate artifacts and transition-state evidence

`POST /api/pathways/rate-artifact` converts sourced forward/reverse activation Gibbs energies into Eyring rate constants. It requires reaction orders, temperature, both barriers, reaction Gibbs energy, evidence provenance, and a detailed-balance tolerance. Optional tabulated pressure multipliers are interpolated in log-pressure only inside their sourced range. Barrier uncertainty is propagated into multiplicative rate ranges and an active-learning priority score.

The result is certified only when evidence is non-hypothetical and

`ΔG‡forward − ΔG‡reverse = ΔGreaction`

within the declared tolerance. The endpoint refuses extrapolated pressure rates or thermodynamically inconsistent forward/reverse pairs.

`POST /api/experiments/{id}/refinement` now accepts `transition_state_certification: true`. For a converged NEB peak it can calculate a finite-difference Hessian, require exactly one significant imaginary mode, displace in both directions, optimize, and compare endpoint connectivity. This is connectivity-following evidence—not a continuous mass-weighted IRC.

## Phase 4: Cantera operation simulation

`POST /api/operations/cantera/simulate` accepts a validated Phase-1 operation definition plus either an allowlisted bundled mechanism (`gri30.yaml` or `h2o2.yaml`) or explicit Cantera YAML text. Supported paths are:

- ideal-gas closed batch, constant volume or constant pressure;
- ideal-gas steady CSTR with one feed;
- ideal-gas constant-pressure Lagrangian PFR with one feed.

Responses include mechanism SHA-256, samples, conversion, optional target formation/stoichiometric ceiling/yield, material and energy balances, and Cantera/SUNDIALS statistics. A response is accepted only when balances close. Semibatch operation, surface catalysis, gas transfer, overall-UA heat loss, mismatched feed conditioning, and non-ideal phases are rejected instead of approximated.

The declared Phase-1 mechanism may be an acknowledged conservation subset; the submitted Cantera mechanism is authoritative for kinetics. This separation is reported in provenance.

## Phase 5: multiphase and electrolytes

`POST /api/multiphase/gas-liquid` integrates one neutral solute with supplied Henry `Hcp` and `kLa` evidence. It supports a closed headspace or an open reservoir with fixed gas partial pressure and returns gas/dissolved histories, mass balance, and uncertainty scenarios. It does not calculate acid/base speciation, ionic activities, salting-out, chemical enhancement, or precipitation.

`GET /api/multiphase/capabilities` reports whether the isolated Reaktoro service is healthy. `POST /api/multiphase/electrolyte/equilibrate` forwards a hashed, strict equilibrium request only when `CHEM_SANDBOX_REAKTORO_URL` is configured. Otherwise it returns HTTP 503; no idealized fallback is manufactured. The Windows launcher configures the verified project-local service automatically.

The pinned service is in `runtimes/reaktoro/`, with Reaktoro 2.13.0 installed under `.tools/reaktoro-env`. It supports closed Gibbs equilibrium, optional fixed pH and gas fugacity, selected aqueous activity models, gas/mineral phases, requested species amounts, pH, and ionic strength. A real Pitzer closed-equilibrium smoke case passes; broader official-example reproduction is still required before discovered mechanisms may be coupled to it.

## Scientific boundary

Cantera predicts the submitted kinetic mechanism; it does not discover missing channels. Reaktoro predicts equilibrium for its selected database/activity model; it does not supply kinetics. Henry transfer predicts physical partitioning for the supplied coefficients; it does not derive those coefficients. Agreement or solver convergence is therefore conditional evidence, not calibrated confidence in real plant yield.
