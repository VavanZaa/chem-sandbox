# Real chemical operation implementation plan

Status: phases 0–8 backend foundations implemented with the explicit scope limits below. Last updated 2026-09-22.

## Goal and claim boundary

Chem Sandbox should progress from bounded molecular exploration to a traceable prediction workflow for laboratory and process conditions. It should estimate possible products, rates, conversion, selectivity, and uncertainty without turning search frequency or endpoint energy into a fake yield.

The target is a decision-support research tool, not an autonomous safety authority. A value is displayed as a quantitative prediction only when the required thermodynamic, kinetic, transport, and operating evidence exists. Missing evidence remains missing.

## Non-negotiable blind-discovery boundary

Target-specific knowledge must not influence internal discovery. The discovery service may receive reactant identities and geometries, amounts, charges/spins, physical conditions, model choices, compute budgets, and generic non-targeted local chemistry operators. It must not receive expected products, named recipes, literature pathways, reference energies, literature yields, or benchmark labels. Generic operators are declared in the manifest and must be selected before any held-out case is opened.

Every run therefore uses this lifecycle:

```text
frozen reactants + physical conditions + compute budget
                    |
                    v
       blind candidate generation and scoring
                    |
                    v
      immutable result + input SHA-256 manifest
                    |
         +----------+-----------+
         |                      |
         v                      v
 post-run reference       mechanism refinement
 evaluation only          and operation modeling
```

Rules:

1. API schemas reject unknown fields, including target products, recipes, known yields, and reference energies.
2. Each result carries a `blind_manifest` containing the frozen input hash and policy declaration.
3. Runtime discovery packages must never import benchmark/reference packages. CI enforces this statically.
4. Benchmark references are opened only after the prediction result has been written and hashed.
5. Held-out benchmark cases cannot be used to tune operators, budgets, thresholds, ranks, or model weights.
6. A future knowledge-assisted mode, if added, must be a separate workflow and score. It must never be reported as blind performance.
7. General pretrained model knowledge is permitted only when the checkpoint and training process predate benchmark selection. Provenance and possible training overlap must be disclosed.

The phase-0 implementation enforces these rules. Reference tables now live under `backend/evaluation/reference_only`, outside the installed application package. The offline evaluator writes a reference-free prediction artifact and SHA-256 before it can write a separately hashed unlocked evaluation artifact.

## Architecture

### A. Blind molecular discovery

- Generate conformers, orientations, association/dissociation proposals, bond rearrangements, proton/electron transfers, and surface events without a product target.
- Use every compatible engine independently. Preserve raw per-engine results rather than averaging incompatible absolute energies.
- Expand states using uncertainty-aware exploration, not expected-product similarity.
- Mark bounded-search absence as unknown, never impossible.

### B. Electronic-structure refinement

- Reoptimize important minima and candidate transition paths using a method appropriate to charge, spin, elements, phase, and catalyst.
- Verify minima with zero imaginary modes and transition states with exactly one relevant imaginary mode.
- Add zero-point, thermal, entropy, standard-state, and solvation corrections to obtain Gibbs energies under stated conditions.
- Enforce atom, charge, spin, and mapped-connectivity conservation. Preserve nonconvergence and domain warnings.

### C. Kinetic mechanism construction

- Convert verified elementary steps into forward/reverse rate laws with units and evidence provenance.
- Enforce microscopic reversibility when thermodynamic data support it.
- Represent pressure-dependent and falloff kinetics where necessary.
- Track rate uncertainty and correlations. Missing channels are unknown, not zero.

### D. Reactor and transport simulation

- Cantera: homogeneous batch, semibatch, constant-volume/pressure, CSTR, PFR, reacting surfaces, heat transfer, and gas flows.
- Reaktoro: electrolyte activities, pH, aqueous speciation, precipitation, and reactive transport in an isolated supported runtime.
- CatMAP or equivalent: catalytic surface microkinetics with explicit sites, coverages, and elementary surface events.
- IDAES: process flowsheets, recycle, separation, heat exchange, optimization, and equipment-level constraints.

Operation predictions consume a frozen, provenance-labelled mechanism. They cannot write products or references back into blind candidate generation.

### E. Calibration and uncertainty

- Calibrate kinetics only against a declared training/calibration set, never held-out evaluation cases.
- Propagate model, barrier, thermochemistry, operating, and numerical uncertainty separately.
- Report ensembles and sensitivity indices rather than presenting cross-model agreement as confidence.
- Add applicability-domain checks and out-of-distribution warnings.

## Implementation phases and acceptance gates

### Phase 0 — scientific firewall and reproducibility

Implementation status: complete for the backend foundation.

Deliverables:

- strict input schemas and target-knowledge denylist;
- immutable blind-run manifest and input digest in every result;
- reference-import boundary test;
- evaluator-only reference store and benchmark unlock-after-run logic;
- environment and model checkpoint lock file, hashes, seeds, hardware, and software provenance.

Gate: deliberately submitted target/reference fields receive HTTP 422; a completed result can be reproduced from its manifest; runtime discovery has no reference imports.

Implemented: strict schemas, blind manifests, full runtime import-boundary test, offline reference package, prediction-freeze/reference-unlock protocol, software/hardware/seed records, and SHA-256 hashes for executed MACE, ANI, and AIMNet checkpoints.

### Phase 1 — operation data model and unit system

Implementation status: validation foundation complete; reactor integration begins in phase 4.

Deliverables:

- reactor types: batch, semibatch, CSTR, and PFR;
- phase definitions, reactor volume, residence time, feed schedules, headspace, heat-transfer mode, agitation, gas transfer, and catalyst loading;
- Pint-backed dimensional validation and explicit standard states;
- material, charge, element, site, and energy balance diagnostics.

Gate: invalid dimensions fail before simulation; every accepted simulation closes balances to a declared tolerance.

Implemented: strict batch/semibatch/CSTR/PFR inputs, phases, inventories, timed feeds, heat transfer, mixing, gas transfer, catalysts, Pint conversion to SI, cross-reference validation, and elementary-reaction element/charge/surface-site balance diagnostics. The validator explicitly leaves energy balance pending until thermochemistry and a converged reactor trajectory exist.

### Phase 2 — free-energy quality

Implementation status: evidence and calculation foundation complete.

Deliverables:

- conformer ensembles and population weighting;
- frequency analysis, zero-point and thermal corrections;
- solution/gas standard-state corrections and configurable solvation evidence;
- charge, multiplicity, electronic-state, and method provenance;
- minimum/transition-state certification.

Gate: no quantity is labelled `delta_G` or `delta_G_dagger` without all required corrections and a provenance record.

Implemented: conformer ensembles and Boltzmann populations; frequency count and imaginary-mode checks; optimizer/force certification; charge, multiplicity, phase, method, solvation, and standard-state provenance; harmonic vibrational ZPE/thermal/entropy components; balanced reaction Gibbs energies; and verified transition-state activation Gibbs energies. Missing, hypothetical, inconsistent, or incomplete evidence returns diagnostics while certified ΔG/ΔG‡ fields remain null. Phase 3 adds optional automatic finite-difference Hessians and connectivity following; continuous IRC remains pending.

### Phase 3 — pathway and rate refinement

Implementation status: foundation complete; continuous IRC and dimer search remain future refinements.

Deliverables:

- NEB/dimer/transition-state search cascade;
- intrinsic reaction coordinate or connectivity-following verification;
- Eyring/Arrhenius rate generation with transmission factors and uncertainty;
- reverse-rate consistency and pressure dependence;
- active-learning queue for high-impact uncertain barriers.

Gate: kinetic networks reject unbalanced, duplicate, unit-inconsistent, or thermodynamically inconsistent channels.

Implemented: opt-in finite-difference Hessian certification of NEB peaks, exactly-one-imaginary-mode checks, and bidirectional imaginary-mode displacement followed by endpoint optimization/connectivity comparison; Eyring forward/reverse rate artifacts with reaction-order units, transmission coefficient, barrier uncertainty, microscopic-reversibility checks, sourced tabulated pressure interpolation, and an impact/uncertainty queue score. Connectivity following is not mislabeled as a continuous IRC, and hypothetical or inconsistent evidence cannot produce a certified rate artifact.

### Phase 4 — Cantera reactor integration

Implementation status: ideal-gas batch, CSTR, and Lagrangian PFR foundation complete; surface reactors, semibatch schedules, heat-loss trajectories, sensitivity, and YAML export remain.

Deliverables:

- exporter from a verified Chem Sandbox mechanism to a versioned Cantera YAML artifact;
- batch/CSTR/PFR adapters with temperature/pressure programs, heat loss, flow, and residence time;
- conversion, selectivity, yield, heat release, species histories, and solver diagnostics;
- sensitivity analysis and mass/element closure checks.

Gate: canonical Cantera examples reproduce documented profiles within numerical tolerance, followed by blind Chem Sandbox mechanism tests.

Implemented: built-in or supplied versioned Cantera YAML loading with SHA-256 provenance; constant-volume/constant-pressure batch, steady CSTR, and constant-pressure Lagrangian PFR execution; species histories, conversion, stoichiometric target ceilings, material/element and energy closure, solver statistics, and strict refusal of currently unsupported transfer/catalyst/UA/semibatch physics. Feed and reactor states must match until explicit conditioning equipment is modeled.

### Phase 5 — electrolyte and multiphase chemistry

Implementation status: native gas/liquid transfer and isolated Reaktoro 2.13 service active in a project-local Windows environment; published-example coverage is still incomplete.

Deliverables:

- isolated Reaktoro service with pinned environment;
- activities, acid/base speciation, ionic strength, precipitation, gas dissolution, and mass transfer;
- explicit CO2 headspace/dissolution model for the bicarbonate validation case.

Gate: published equilibrium/speciation examples are reproduced before coupling to discovered mechanisms.

Implemented: sourced Henry-Hcp/kLa gas-liquid transfer for closed headspace or open fixed-partial-pressure operation, BDF trajectories, conservation checks, and parameter-uncertainty scenarios. An isolated Reaktoro 2.13 service implements closed equilibrium, selected aqueous activity models, optional gas/mineral phases, fixed pH/fugacity constraints, requested species, pH and ionic-strength output, plus request hashing. Reaktoro 2.13.0 is installed under `.tools/reaktoro-env`; its Pitzer closed-equilibrium smoke case passes and `start.cmd` launches the service automatically. The main API still returns HTTP 503 when that service is stopped and never substitutes fabricated electrolyte speciation. Broader published-example reproduction remains the coupling gate.

### Phase 6 — catalysts

Implementation status: explicit single-site mean-field microkinetic foundation complete.

Deliverables:

- explicit homogeneous catalyst states or periodic catalyst surfaces;
- adsorption, desorption, diffusion, coverage, and site balances;
- CatMAP-compatible microkinetic export or an equivalent maintained solver;
- catalyst deactivation and transport hooks.

Gate: catalyst-free endpoint cancellation is never presented as catalytic validation; a catalyst result requires a catalyst-bound pathway and site balance.

Implemented: strict gas/surface species, one explicit vacant-site pseudo-species, adsorption/desorption/surface-reaction/diffusion steps, element/charge/site conservation, forward/reverse Eyring rates, microscopic reversibility, BDF coverage integration, target turnover frequencies, first-order deactivation and sourced transport-effectiveness hooks. Each elementary step must change a non-vacant bound intermediate, so an overall reaction with a catalyst label that cancels at both endpoints is rejected. A CatMAP setup fragment exporter emits reaction expressions, site definitions, gas pressures, and multisite occupancy while refusing to invent its separate energetic input table.

### Phase 7 — process flowsheets

Implementation status: one-reactor/separator/purge/recycle steady-state foundation complete with native and IDAES cross-check engines.

Deliverables:

- IDAES unit models, recycle, separators, utilities, and optimization;
- uncertainty-aware objectives and constraints;
- safety/property limits and infeasibility diagnostics.

Gate: every optimized flowsheet has converged unit balances and shows assumptions, solver status, and sensitivity to uncertain kinetics.

Implemented: sourced single-pass conversion and per-component separator recoveries, mixer/reactor/separator/purge/recycle balances, fixed-point recycle convergence, elemental closure, target production, overall conversion, utility duty, bounded recycle optimization, internal-flow feasibility limits, and low/high parameter scenarios. The official Windows IDAES 2.13/IPOPT extension bundle is installed project-locally; an independent IDAES `FlowsheetBlock` algebraic model is solved and compared with the native result. Rigorous property packages, phase equilibrium, equipment sizing, pressure drop, relief design, and arbitrary unit graphs remain outside this foundation.

### Phase 8 — calibration and validation program

Implementation status: prediction-freeze, calibration, and held-out evaluation foundation complete; validation tiers remain an expanding program rather than a one-time claim.

Validation tiers:

1. numerical invariants and analytic micro-networks;
2. library reference examples;
3. held-out gas-phase elementary reactions;
4. solution acid/base and esterification cases;
5. pressure/temperature industrial cases;
6. catalytic surface cases;
7. complex molecules such as glucose hydrogenation;
8. prospective experiments unknown during model development.

For every tier, freeze inputs and model versions first, then unlock references. Report coverage, failures, calibration, error distributions, uncertainty coverage, and compute cost—not only successful examples.

Implemented: SHA-256 prediction freezing, affine or offset-only weighted calibration, calibrated uncertainty scaling, immutable calibration digests, strict calibration/held-out ID separation, tamper detection, MAE/RMSE/bias/Gaussian-NLL metrics, and 50/80/95% interval coverage. The Phase 3–5 validation runner freezes raw results before importing reference-only checks; its current six checks cover Eyring detailed balance, analytic closed-headspace equilibrium, reactive Cantera batch/CSTR/PFR closure, and Reaktoro fixed-pH equilibrium.

## Initial validation cases

- Sodium bicarbonate + acetic acid in water: acid/base transfer, carbonic-acid decomposition, CO2 gas/liquid partition and escape. This is a multiphase validation, not suitable for gas-phase endpoint energies alone.
- Fischer esterification: reversible solution kinetics, water activity, acid catalysis, and equilibrium limitation.
- Water-gas shift: gas-phase equilibrium plus explicit catalytic surface mechanism for kinetic prediction.
- Glucose hydrogenation: conformers, stereochemistry, solvent, hydrogen pressure, catalyst surface, and competing products.

These cases cover increasing physical complexity; their recipes and reference outcomes remain evaluator-only.

## Installed foundation and runtime decisions

- Installed: Cantera 3.2.0 and IDAES 2.13.0 in the project virtual environment, plus the official project-local IDAES/IPOPT Windows extension bundle.
- Existing: ASE, RDKit, MACE-MP, ANI-2x, AIMNet2, AIMNet2-RXN, SciPy, NetworkX, PyTorch CUDA.
- Reaktoro 2.13.0 installed in a project-local win-64 Micromamba environment with an isolated service and verified Pitzer equilibrium smoke case. CatMAP text export is available, while execution of CatMAP itself and RMG-Py remain deferred to isolated runtimes.
- WSL is not installed and is not needed for the current Reaktoro win-64 service. Any future Linux-only RMG/CatMAP deployment remains a separate operating-system/runtime decision.

## Near-term work order

1. Couple verified surface mechanisms to Cantera surface reactors and packed-bed transport without weakening the Phase-6 site-balance gate.
2. Expand the Phase-7 graph from one recycle loop to unit-specific IDAES property packages, heat exchangers, pressure changers, phase separators, costing, and safety constraints.
3. Complete continuous IRC/dimer search, Cantera YAML export, heat-transfer programs, and global sensitivity reporting.
4. Expand prospective Phase-8 datasets and Reaktoro published-example coverage while keeping calibration and held-out IDs disjoint.

## Definition of a defensible prediction

A result may be called a predicted operation yield only when it contains: a balanced mechanism; method/domain provenance; verified or sourced forward and reverse kinetics; phase/activity model; reactor type and operating schedule; heat and relevant mass-transfer assumptions; solver convergence; balance closure; uncertainty/sensitivity; and a statement of validation domain. Otherwise the UI must call it a conditional scenario, exploratory population, or theoretical ceiling.
