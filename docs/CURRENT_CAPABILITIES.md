# Chem Sandbox capability boundary (2026-09-27)

This inventory separates a **calculation the software can run** from a
**prediction established by independent experiment**. It covers the current
local build, not the full ambitions in `REAL_OPERATION_IMPLEMENTATION_PLAN.md`.
The backend has 398 collected tests. A full 397-test run passed immediately
before the last reduced-HCN contact regression was added; the final focused
HCN/API regression run passes separately. The frontend production build
passes, and the installed PyTorch build sees the local RTX 4080. Those checks
confirm software operation, not scientific accuracy.

## In the visible app

| Workspace | What it can do now | What its output does not prove |
| --- | --- | --- |
| Network and Outcomes | Enter supported molecules/SMILES, edit amounts and units, perform a bounded product-candidate search, inspect a draggable graph and branch stopping reasons, and calculate a separate stoichiometric ceiling for the bicarbonate example. | Search frequencies are not product probabilities. A missed product is not ruled out; a sampled structure is not a verified reaction mechanism or bulk yield. |
| 3D inspector and stereochemistry | Rotate/zoom saved molecular geometries; preserve aromatic and explicit input R/S or E/Z identity for intact mapped components. | New fragments may lack reliable bond order, charge, or stereo identity. R/S separation alone does not predict enantiomeric excess. |
| Engines | Run compatible MACE-MP, ANI-2x, AIMNet2, and AIMNet2-RXN energy/force comparisons on CUDA; attempt bounded NEB electronic-path refinement. | Raw electronic energies are not solution free energies, barriers, rates, or calibrated confidence. Model agreement is not independent evidence. |
| State Kinetics and Efficiency | With complete user-supplied rates or sourced barriers, solve conditional state populations, reaction-network conversion, selectivity and target yield under declared conditions. | They cannot invent missing channels, activation data, solvent effects or true industrial yield. State populations are not bulk multi-reactant yield. |
| Protocol | Execute ordered reaction stages, mid-run additions, changed conditions, holds, sampling and imperfect transfers; check material balances; export SBML and cross-check with RoadRunner. | A staged glucose or Krebs-cycle **example** does not infer enzyme pathways, missing cofactors, or real rates. |
| Origins | Run a bounded closed boil/spark/condense loop; solve supported mixture EEDFs and collision-rate diagnostics; couple selected neutral/charged gas chemistry and gas/water capture; inspect compartment inventories, conservation, controls and numerical sensitivities. | Neither the historical experiment nor a modern wet-spark HCN yield has been quantitatively reproduced. Some electron-impact product states, spark exposure, low-temperature chemistry, transport and heat flows are unresolved. |

Origins also offers a **conditional reduced HCN route** that accepts
externally supplied effective fragment hazards, runs GRI gas chemistry and
one-compartment HCN contact at fixed pH, and rejects atom/charge-imbalanced
or direct C–N fragment injection. Its template source hazard is zero and
it supplies no nitrogen-fragment source; the total N2 collision table cannot
justify ground-state atoms. A separate element-count ceiling `min(C, N)` for HCN+CN− (and
`min(C, N, H)` for molecular HCN) is available before or after apparatus
simulation. Neither is a validated yield. See
[the reduced-route guide](REDUCED_HCN_AND_ATOM_CEILING.md).
An optional assay-free audit repeats the reduced calculation with a finer
solver step and reports source-off and atom-removal controls. It detects some
numerical or conservation problems; it does not establish real spark rates
or a measured HCN yield.

The active HCN benchmark has shifted from a wet Miller–Urey apparatus to a
[dry flowing CH4/N2 discharge](DRY_FLOW_HCN_GATE.md). A gas-packet probe can
report the complete simulated outlet composition, methane conversion and
carbon-normalized product amounts without any water chemistry. Under the
current state-safe N2 handling its HCN remains below the software reporting
floor; older positive dry-probe values are archived, not validated results.
An evaluator-only check against the 2024 dry-flow study now finds the model's
1% CH4 product ordering reversed: it predicts C2H6 above HCN, whereas the
experiment reports HCN above C2H6 in that feed range. The comparison is
qualitative because streamer exposure, N2 fragment states and the exact
flow/temperature basis are not yet matched; a 189–435 Td scenario sweep does
not remove the ordering failure.

Origins also displays saved full-week development scenarios. The UI now
labels those older numbers **archived** because they used a superseded N2
electron-impact fragment-state assignment; they are not current predictions.

Engine coverage is conditional: ANI-2x is used for neutral closed-shell
H/C/N/O/F/S/Cl systems; AIMNet2 covers its supported main-group elements
and charges but not sodium here; AIMNet2-RXN is restricted to neutral
closed-shell H/C/N/O reactions; MACE-MP has broader element support but
does not receive formal charge/spin as model inputs. Passing an input-domain
check still does not establish accuracy for that reaction or environment.

## Backend tools that require explicit evidence

- Dimensional and atom/charge checks on operation definitions; harmonic and
  conformer thermochemistry, reaction/activation Gibbs quantities, and
  reversible Eyring rate artifacts **only when** all required electronic,
  frequency, standard-state, solvation and transition-state evidence is
  supplied and internally consistent. The app does not produce a generally
  certified transition state or solution free energy on demand.
- Ideal-gas batch, steady CSTR and Lagrangian PFR simulations using a supplied
  Cantera mechanism; reaction rates are **inputs**, not products of the reactor
  solver. General semibatch and full heat-loss/surface-reactor operation are
  not implemented in this adapter.
- Single-neutral-solute Henry/kLa gas-liquid transfer with supplied parameters.
  The isolated Reaktoro equilibrium service is installed and started by
  `start.cmd`, but an API session without its configured service has no
  electrolyte-equilibrium fallback. Equilibrium is not aqueous kinetics.
- Explicit one-site-family mean-field catalyst microkinetics with sourced
  reversible barriers and fixed gas reservoirs. The CatMAP exporter makes
  setup text, not a running CatMAP calculation or catalyst discovery.
- One-reaction reactor/separator/purge/recycle flowsheet and bounded
  optimization, with an available IDAES/IPOPT algebraic cross-check. It does
  not size real equipment, compute rigorous phase behavior, or certify safety.
- Prediction freeze, separate calibration, and held-out evaluation with hash
  checks. This protects validation splits; it does not make a poor model
  accurate. Blind discovery rejects target products, recipes and yields as
  inputs; a user-declared protocol or a literature mechanism is a different,
  post-discovery task.

## Feasible projects with the current program

1. Teach or explore stoichiometry, reactant limits, named species and 3D
   structures in small supported systems.
2. Screen **hypotheses** for candidate molecules or reaction paths, including
   some organic/aromatic and stereochemically labelled inputs, then decide
   which need higher-level calculation or an experiment.
3. Compare multiple neural-energy engines on in-domain isolated molecules;
   record disagreement and geometry failures rather than report a consensus
   probability.
4. Simulate a **declared** multistep laboratory synthesis or biochemical
   network, including a reagent added halfway through, provided structures,
   balanced reactions and appropriate rates are supplied.
5. Explore conditional batch/flow-reactor, gas-liquid contact, electrolyte,
   catalyst-site and recycle scenarios with independently sourced parameters.
6. Study how assumed spark field, collision data, gas reactions and capture
   change a prebiotic **model**, test negative controls and conservation, and
   design the next measurement needed for HCN validation.

The [saved ten-case benchmark](../backend/data/benchmarks/reaction_suite_v01/latest.md)
is a cautionary example: AIMNet2 gives the correct reaction-energy *sign*
for all nine cases with a reference, yet raw electronic-energy deviations
against experimental 298 K enthalpies remain tens of kJ/mol and the
four-engine spread for benzene/pyridine hydrogenation is about
746–798 kJ/mol. Because electronic energy and reaction enthalpy are
different observables, those deviations are diagnostics rather than a
formal thermochemical accuracy score.

## Projects this build cannot defend yet

- Predicting an arbitrary mixture's complete products and calibrated product
  probabilities from its ingredient list alone.
- Reproducing a measured Miller-Urey or Barth wet-spark HCN amount, amino-acid
  yield, or seven-day energy efficiency from independent inputs.
- Inferring a realistic glucose or Krebs-cycle pathway and its enzyme-specific
  rates from molecular formulas, or predicting real chiral selectivity from
  R/S labels alone.
- Ranking complex aromatic, charged, metallic, catalytic or solvated systems
  as quantitatively accurate on the basis of the current four neural models.
- Designing or certifying an industrial plant, scale-up, runaway/safety
  envelope, emissions permit, operating procedure, cost or guaranteed yield.
- Running RMG-Py automatic mechanism generation, a CatMAP calculation, or a
  general-purpose reactive-transport/CFD/plasma-geometry solver as part of
  the current integrated workflow.

## Why the HCN gate remains open

The charged Pavan branch forms trace HCN from the declared inorganic Barth
feed in an assay-free model; the primary GRI result is below the reporting
floor. The charged branch shares all GRI HCN reaction records, so the two
are not independent kinetic validation. We do not know the real fraction of
gas exposed to the spark, product-state-resolved N2 dissociation, the wet
mixture CH+N2 association and HCNN/HNCN/NCN fate, several low-temperature
HCN/H2CN/HNC reactions, the experiment's gas-water contact coefficient,
or a complete thermal-energy budget. The time-dependent aqueous network is
also incomplete; even a current low-contact full-duration run reaches
pH below the published hydration-rate domain, so its old "source-bounded"
result has been retired. Recent assay-free rate stresses and 192/384-step checks
pass conservation and negative controls but quantify **model sensitivity**,
not a real HCN uncertainty or experimental yield. The held-out HCN assay has
not been used to tune these calculations. See `HCN_FORMATION_GATE.md` and
`BARTH_LONG_HORIZON_GATE.md` for numbers, sources and pending tests.

## Practical strength and weakness

The strongest part is **traceability**: explicit domains and inputs,
conservation ledgers, no silent invention of missing evidence, distinct
model outputs, reproducible numerical checks, and held-out validation
boundaries. Its weakness is not a missing GPU. GPU accelerates the neural
energy/force engines, while many kinetics/graph/transport tasks naturally
run on CPU. The decisive weakness for predictive yield is missing or
inapplicable **physical and kinetic evidence**, especially for wet plasma
and aqueous chemistry. A converged solver cannot correct that by itself.

Further scope details: `README.md`, `STAGED_LAB_PROTOCOLS.md`,
`PHASE_3_5_API.md`, `PHASE_6_8_API.md`,
`STEREOCHEMISTRY_AND_AROMATICS.md`, and the latest reaction benchmark.
