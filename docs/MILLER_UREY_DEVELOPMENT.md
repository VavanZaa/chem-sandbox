# Miller–Urey Computational Experiment

Status: Gates 1–2, product-agnostic closed-shell candidate enumeration, and
held-out post-discovery evaluation are implemented in Chem Sandbox v0.25.
This document does not claim that the current models can predict validated
spark-discharge yields.

## Implemented through v0.13

- Four explicitly tracked compartments: boiler liquid, gas loop, condensate
  trap, and wall inventory.
- Repeated boil, discharge, condensation, reflux, and wall-loss cycles.
- Mid-cycle feed additions and analytical sampling with a whole-run ledger.
- Separate gas/liquid mass-action and reduced electron-impact channel types.
- Pulsed voltage/current/duty/coupling electrical-energy accounting.
- Per-reaction atom and charge validation plus complete apparatus closure.
- Enantiomer-pair diagnostics without numerically forcing racemic output.
- A scientific firewall requiring a frozen manifest before a mechanism may
  claim blind-discovery provenance.
- A Cantera 3.2 capability probe that loads its packaged electron-collision
  plasma reference mechanism without transport.
- A strict cross-section endpoint that constructs a real Cantera
  `ElectronCollisionPlasmaRate` and fingerprints its provenance.
- Explicit Maxwellian-approximation or normalized tabulated electron-energy
  distributions (EEDFs), with no silent inference from voltage alone.
- Numerical evaluation of `<sigma*v>` in m3/s, an electron-density-dependent
  pseudo-first-order rate, and a target lifetime.
- Electron species, ion charge closure, and cross-section/EEDF channels coupled
  to the gas-loop discharge stage of the recirculating apparatus.
- Event-limited discharge integration: declared electron-channel energy costs
  stop further electron chemistry when the coupled electrical deposit is spent.
- Cantera's installed `lxcat2yaml` converter is available for user-obtained
  LXCATML files; source licensing and dataset selection remain external choices.
- A private local Morgan/LXCat complete-set download covers H2, H2O, CH4 and
  NH3 with 46 parsed processes. Its BOLSIG+ text is SHA-256 verified and is not
  redistributed with the source tree.
- BOLOS 0.2 solves the steady-state two-term electron Boltzmann equation from
  mixture mole fractions, gas temperature and reduced electric field. The API
  returns the normalized EEDF, mean energy, convergence diagnostics and 42
  inelastic rate coefficients for the full four-gas set.
- Cantera 3.2 independently solves a packaged CH4/H2O reference case when the
  submitted mixture lies inside that mechanism's electron-collision coverage.
- The effective-to-elastic reconstruction warning for H2, H2O and NH3 is
  elevated into the result rather than hidden in solver logs.
- UI and JSON export through the Origins workspace.

## Added in v0.14

- Energy-limited neutral electron-impact fragments feed an independent
  Cantera GRI heavy-gas network without treating raw collision opportunities
  as simultaneous reaction extents.
- Reviewed Henry-law records move supported HCN, aldehyde, alcohol, and
  ammonia species between gas and water while conserving each species.
- A product-agnostic graph enumerator combines only formulas actually present
  in the final aqueous inventory, generates neutral closed-shell acyclic CHNO
  structures, and expands unassigned stereocenters into explicit isomers.
- The generated discovery artifact is serialized and SHA-256 frozen before a
  separate Strecker evidence file or historical observation file is opened.
- A one-way cycle replay projects aminonitrile equilibration and hydrolysis
  sensitivity from each stored apparatus cycle. It does not feed back to the
  plasma, gas network, phase partition, or candidate generator.
- A versioned 1953 material/time profile represents 200 mL water, initial
  100/200/200 mmHg H2/CH4/NH3, and one week. Unknown discharge-current,
  coupling, pulse, reduced-field, and transport values remain explicitly
  hypothetical rather than fitted to products.
- The held-out benchmark reports condition matching, broad candidate identity
  coverage, precursor support, blockers, and an explicit `validated_reproduction:
  false` until complete pathways and evidence support that stronger claim.
- An optional background job ranks seven generic amino-carboxylic candidates
  formula-by-formula with MACE-MP, ANI-2x, AIMNet2, and AIMNet2-RXN. Absolute
  energies are never mixed across model families and rankings are not treated
  as probabilities.
- Exact R/S pairs begin from reflected coordinates and receive tied achiral
  ranks. Any raw mirror-energy difference is retained as a numerical
  diagnostic rather than interpreted as enantioselectivity.

## Added in v0.15

- Broad historical identity recovery is renamed `enumeration_coverage` and is
  explicitly forbidden from being interpreted as recall or prediction quality.
- A compact four-model prediction artifact is canonicalized and SHA-256 frozen
  before the 1953 observation file is opened. Only then is a class-conditional
  top-k score calculated. The live four-model artifact recovers glycine and
  alpha-alanine in the top three, while beta-alanine ranks fourth (2/3 top-3
  observation recall).
- Ranking artifacts, job state, and held-out evaluation artifacts persist across
  service restarts under `backend/data/runs/prebiotic-ranking/`.
- A versioned Stribling–Miller precursor benchmark records approximately 10
  nmol/J HCN and 2 nmol/J formaldehyde for the published CH4 case with
  H2/CH4 = 1. The code reports current simulated energy yields but refuses to
  score them when N2 is absent, NH3 is present, or the H2/CH4 ratio differs.
  Even with a matching feed it withholds an error score until the simulated
  observable uses the published zero-time extrapolation and the energy basis
  is tied to comparable calorimetry rather than a declared electrical deposit.
- A second collision-set option combines the private hash-verified Morgan
  H2/H2O/CH4/NH3 tables with 25 hash-verified Phelps N2 processes packaged by
  Cantera 3.2. BOLOS converges for N2-containing mixtures and the apparatus
  reports N2 excitation/ionization source terms. No N2 neutral-fragment source
  is invented: the current Phelps table and GRI coupling do not yet supply a
  validated N2-to-reactive-nitrogen plasma pathway.
- The historical 1953 development profile currently produces final-inventory
  diagnostics of 8.754e-5 nmol/J HCN and 0.5459 nmol/J formaldehyde at the
  declared hypothetical electrical settings. These are not compared as an
  error metric because the reference values were extrapolated to zero time in
  a different feed and apparatus protocol.

## Added in v0.16

- The public QST/Tabata recommended `e + N2 -> e + N + N` product cross
  section is stored with its source, units, retrieval date, and SHA-256. The
  16-point table spans 12–200 eV; the associated analytic fit reports a
  9.76 eV threshold.
- The neutral-dissociation cross section is folded through the converged EEDF
  after the 71-process Morgan/Phelps solve. It is intentionally excluded from
  the Phelps energy-loss assembly because counting it again alongside the
  electronic-excitation manifold would double count electron energy loss.
- The resulting two neutral N atoms pass atom/charge validation and enter the
  Cantera heavy-gas network only after target-inventory and deposited-energy
  limits are applied.
- The live 20-cycle N2-pathway preset with H2/CH4 = 1 closes its full material
  ledger, forms HCN, formaldehyde, methanol, HNCO, and other network species,
  and retains approximately 9.005 nanomol HCN and 2.614 micromol formaldehyde
  in the modeled water compartments. Its displayed final-inventory diagnostics
  are 0.01251 nmol/J HCN and 3.630 nmol/J formaldehyde. These remain unscored
  under hypothetical electrical controls and are not calibrated Miller–Urey
  yields.
- The deposited-energy cap is now allocated across every imported collision
  opportunity, including nonreactive excitation and unresolved channels,
  before resolved fragments enter chemistry. In the current last cycle only
  13.55 J of the 36.0 J collision allocation reaches resolved fragment
  chemistry; the UI exposes this allocation rather than hiding it.
- Optional resistance-heater calorimetry can override the electrical estimate
  without erasing either ledger. Cycle-resolved aqueous inventories spanning
  1 hour to 7 days are regressed to a zero-time nmol/J intercept with its
  standard error. A synthetic held-out fixture exactly recovers 10 nmol/J HCN
  and 2 nmol/J formaldehyde, while live runs remain unscored when the published
  weighting and uncertainty are unavailable.
- The Origins UI now includes a `Load N2 reactive pathway` preset and displays
  the reviewed reactive-nitrogen source separately from raw excitation and
  ionization opportunities.

## Added in v0.17

- Cantera's packaged Pavan 2023 plasma-assisted methane mechanism is loaded as
  a SHA-256-recorded downstream shadow (`3789dc895d9a...`). It contains 71
  species and 440 reactions and is based on a published extension of GRI-Mech.
- The same collision-wide energy scale used by the primary branch now supplies
  mapped N2(v1-v8), N2(A), N2(B), N2(C), and N2(a') populations to a 387-reaction
  neutral/heavy Pavan branch alongside the resolved neutral fragments.
- The 53 packaged electron-collision and two-temperature ion/electron reactions
  are excluded from this shadow because Cantera's zero-dimensional ideal-gas
  reactor cannot integrate the source `thermo: plasma` phase. They are not
  silently treated as zero-rate validated chemistry.
- GRI and Pavan per-cycle net changes are displayed side by side. Pavan
  extends GRI, so this is explicitly a correlated mechanism sensitivity rather
  than independent model confidence or a product probability.
- The shadow cannot mutate the primary apparatus or aqueous inventories. Its
  signed net changes expose both formation and consumption and remain
  mechanism-spread diagnostics only.
- A versioned Stribling–Miller timing profile now supplies 168 hourly apparatus
  cycles from 1 hour through 7 days with the published CH4/N2/H2O/H2 feed family
  and H2/CH4 = 1. Its absolute gas amounts and discharge remain explicitly
  hypothetical, and calorimetric deposited energy is deliberately absent, so
  the profile exercises zero-time extrapolation without fabricating a score.

## Added in v0.18

- The complete open supplementary data from Kamarinopoulou et al., *Science
  Advances* 10 (2024), adl4246, is stored with its CC BY-NC 4.0 provenance and
  source-PDF SHA-256. A reproducible importer extracts Data S1 and Data S2 and
  creates a stable, hash-verified Cantera mechanism.
- All 83 published thermodynamic species and 251 reaction records are
  preserved. The executable shadow retains 189 atom-balanced neutral/heavy
  reactions. The paper's 62 fixed 300 Td, electron-density-specific impact
  rates are archived but excluded, because the sandbox supplies electron
  sources from its own mixture EEDF.
- The CH4/N2 mechanism is integrated only for its published 18 ms microkinetic
  residence-time domain. It is a non-propagating local sensitivity branch, not
  an apparatus yield. It shares GRI/Glarborg ancestry with the other branches,
  so agreement does not become independent confidence. The 168-cycle replay
  completes all local solves; cycle 161 uses a bounded smaller-step retry that
  changes numerical controls only and is exposed in the result.
- Boiler and condenser gas inventories now enforce a pure-water saturation
  capacity using Cantera Water. Supersaturated water returns to the contacting
  liquid before plasma chemistry. Henry partitioning uses the resulting
  inventory-derived pressure; a declared pressure is retained as a reference
  and no hidden vent or pressure controller is applied.
- The corrected 168-hour Stribling-Miller scenario stays between 1.627 and
  1.701 bar instead of the former impossible roughly 74 bar gas state. Its
  zero-time diagnostics are 0.2449 +/- 0.0036 nmol/J HCN and 2.1166 +/- 0.0159
  nmol/J formaldehyde. The formaldehyde observable is near the published
  approximately 2 nmol/J value, while HCN remains roughly forty-fold low versus
  approximately 10 nmol/J. The run remains unscored because its energy is an
  electrical scenario rather than measured calorimetry and the public source
  does not supply the original regression weighting or uncertainty.

## Added in v0.19

- An evidence-labelled `discharge_zone_volume` is now part of the recirculating
  apparatus. It must be positive and cannot exceed the gas-loop volume.
- Electron-collision opportunities are evaluated from the gas amount and volume
  occupying the active zone. The primary and shadow heavy-gas networks use the
  whole-loop mean residence time `stage time * active volume / loop volume`.
- When the value is absent, the former full-loop assumption remains available
  only as an explicit upper-bound fallback with a warning. The shipped
  development profiles declare an illustrative 50 mL active zone inside the
  5 L loop and label it hypothetical rather than historical geometry.
- The Origins result view now auto-scrolls into view and reports the active-zone
  fraction and effective per-cycle chemistry exposure.
- The corrected 168-hour replay uses a 1% active-zone fraction and 12 s mean
  heavy-chemistry exposure per hourly cycle. Its zero-time diagnostics are
  0.24355 +/- 0.00360 nmol/J HCN and 2.11876 +/- 0.01606 nmol/J formaldehyde.
  The small change confirms that the current fragment source is dominated by
  the deposited-energy cap. HCN remains roughly forty-fold below the published
  diagnostic, so the residence correction does not justify a reproduction
  claim.

## Added in v0.20

- A separate Pavan plasma shadow now retains 406 reactions: 387 neutral/heavy
  reactions plus all 19 published two-temperature ion/electron reactions.
- The external, mixture-specific EEDF remains the only electron-collision
  driver. All 34 packaged Pavan electron-collision reactions are disabled so
  the same collision event cannot be counted twice.
- Energy- and inventory-limited collision extents enter as time-averaged,
  charge-balanced sources during the fixed-temperature/volume residence step.
  Positive ions receive matching electron sources; unsupported anions and
  fragment ions remain explicit unresolved-channel records.
- Direct SciPy BDF integration avoids Cantera's currently unavailable plasma
  internal-energy reactor equation while still evaluating Cantera's published
  reaction rates. Every solve checks both elemental and electrical charge
  closure.
- The full 168-cycle timing-profile replay completes 168/168 charged solves.
  Its largest final charge residual is 1.21e-16 mol charge equivalents and its
  largest atomic residual is 1.33e-10 mol against approximately 0.1-mol-scale
  inventories. The whole-apparatus primary ledger remains closed.
- The charged shadow cannot mutate primary gas/water inventories. It is shown
  beside the neutral branches as correlated mechanism sensitivity, not as an
  independent model vote, product probability, or validated yield.

## Added in v0.21

- The charge-balanced Pavan state is propagated as an isolated alternative
  apparatus branch instead of being restarted from the primary gas inventory
  on every cycle.
- Gas, boiler, condensate, wall, electron, ion, and excited-state inventories
  persist independently. Shared ground-state species use the same water-vapor,
  Henry-partition, reflux, sampling, and wall-transfer operations as the
  primary branch.
- Each propagated cycle solves a new two-term EEDF from the alternative
  branch's own neutral gas composition. The branch EEDF and Pavan chemistry do
  not alter the primary GRI inventory or its candidate-discovery artifact.
- The 168-hour replay completes 168/168 propagated charged cycles. Its
  end-to-end atom/charge ledger closes to 8.67e-16 mol and its largest local
  charge residual is 2.14e-16 mol charge equivalents.
- The propagated branch predicts zero-time diagnostics of 0.000291 +/- 0.000017
  nmol/J HCN and 0.7102 +/- 0.0085 nmol/J formaldehyde. These are lower than
  both the primary branch and the published diagnostics, so propagation
  falsifies rather than repairs the earlier optimistic local-shadow reading.
- Final-state pathway diagnostics identify `CH3 + N <=> H2 + HCN` as the
  leading HCN formation direction and `H + HCN (+M) <=> H2CN (+M)` as an almost
  exactly offsetting consumption direction. No literature yield was used to
  select or tune either reaction.

## Added in v0.22

- `continuous_overlap` apparatus operation divides each declared circulation
  cycle into configurable operator-split boil, discharge, condensation, and
  reflux substeps while preserving the declared total duration and electrical
  energy.
- Both historical profiles use two substeps per declared hour. This is an
  explicit provisional numerical setting, not inferred historical apparatus
  evidence.
- A 24-hour convergence study at 1, 2, 4, and 8 substeps closes both ledgers at
  every resolution. Between 4 and 8 substeps, primary formaldehyde changes by
  0.423% and propagated-branch HCN by 0.0122%, but primary HCN changes by 89.0%
  and propagated-branch formaldehyde by 10.1%; the model is therefore not
  timestep-converged as a whole.
- The full 168-hour, two-substep profile completes 336 plasma/kinetics/
  partition steps. Primary zero-time diagnostics are 0.13550 +/- 0.00554
  nmol/J HCN and 2.82967 +/- 0.00458 nmol/J formaldehyde. The propagated
  charge-balanced branch gives 0.0003071 +/- 0.0000218 nmol/J HCN and 0.99465
  +/- 0.01277 nmol/J formaldehyde.
- Both end-to-end ledgers close; the propagated branch's largest final
  conserved-quantity residual is 4.27e-14 mol. Quantitative yield remains
  disabled because HCN is far below the held-out diagnostic, energy is not
  calorimetrically matched, and the operator splitting is not converged.
- The UI reports operation mode, substeps per cycle, per-substep spark-zone
  exposure, and an explicit `convergence pending` label.

## Added in v0.23

- The primary gas mechanism can distribute its energy-limited, EEDF-derived
  neutral-fragment source through the spark-zone residence time instead of
  adding it as an instantaneous impulse. The shipped Stribling–Miller timing
  profile now selects this `time_averaged` source. The propagated charged
  branch already distributes its source through each residence step.
- The 24-hour source-forcing convergence study spans 4, 8, 16, and 32
  substeps/hour. At 16 to 32 substeps, primary and charged HCN each change by
  about 0.1%; primary formaldehyde changes by 1.54% and charged formaldehyde
  by 1.28%. This is only partial 24-hour convergence, not full-week yield
  convergence.
- An optional `condensate_nonwater_return_fraction` applies uniformly to
  **every** condensed non-water species. The 0%-return sensitivity case models
  a protective cold-trap scenario without naming or selecting expected
  products; it is not a measured collection efficiency. The default remains
  100% return to preserve the prior profile for comparison. The 0%-return
  24-hour result is not converged: from 16 to 32 substeps, isolated HCN rises
  by 84% (primary) and 74% (charged), while formaldehyde rises by 53% and 49%.
  The trap integration therefore needs a flow-resolved transport treatment
  before yield use.
- Three completed 168-hour scenarios and their zero-time fit quality are
  available in the Origins UI through a read-only validation endpoint. Their
  primary HCN intercepts are 0.13550 nmol/J (impulse, all condensate returned),
  0.00026393 nmol/J (time-averaged, all returned), and 0.00208505 nmol/J
  (time-averaged, hypothetical non-water holdback in an isolated collection
  reservoir). The large source-treatment spread is a model defect/sensitivity,
  not experimental uncertainty. The holdback formaldehyde fit has R² = 0.103,
  so its intercept is unsuitable as a benchmark value. The collection reservoir
  does not re-equilibrate with the gas or participate in active-boiler kinetics.
- Both branches conserve material in each recorded scenario. The backend has
  188 passing tests and the frontend build passes. Neither the published
  energy deposition nor transport coefficients are sufficiently constrained,
  and the apparent HCN deficit persists. No quantitative-yield claim is made.

## Added in v0.25

- A new, optional `coupled_flow_transport` mode advances boiler gas/water
  exchange and condenser extraction in the same stiff GRI and charged Pavan
  chemistry solve. Its flow contact is a mean-field Henry-equilibrium upper
  bound, not a resolved three-dimensional apparatus or a measured collection
  efficiency. The contact-water hold-up remains a declared hypothesis.
- The 24-hour 2/8/32-substep sensitivity run has both ledgers closed and
  changes HCN/formaldehyde amounts by no more than 0.0142% from 2 to 32
  substeps. This removes the previous artificial dependence on repeated
  boiler-equilibrium jumps for this scenario, but does not by itself establish
  seven-day convergence or experimental accuracy.
- The initial 168-hour two-substep result has primary HCN 0.0910 nmol/J and
  charged HCN 0.0881 nmol/J, still far below the published approximately
  10 nmol/J CH4 benchmark. The energy input is not the source experiment's
  calorimetric measurement, and the primary formaldehyde linear-fit R-squared
  is approximately 0.0004. A four-substep full-week replay agrees with the
  two-substep run to within 0.008% for the tracked final amounts and fitted
  intercepts. This is not a reproduced Miller-Urey result.
- Final-state instantaneous reaction diagnostics identify near-canceling
  GRI HCN formation from methyl radical plus atomic N and HCN consumption
  to H2CN. These are mechanism-level diagnostics, not independent product
  probabilities. The next scientific gate is measured apparatus power/contact
  constraints and a more appropriate, independently validated low-temperature
  plasma-afterglow/aqueous chemistry mechanism, tested without using expected
  products or benchmark yields to choose reactions.
- A three-point, unfitted 24-hour condenser contact-hold-up sensitivity
  (0.01/0.10/0.25 L, four substeps/hour) changes primary collected HCN by
  53.8%. This is a range across declared hypothetical inputs, not a
  probability distribution or calibrated confidence interval.

## Added in v0.24

- The optional `flow_limited_capture` condenser uses the measured-property
  Henry subset as a product-agnostic per-pass capture upper bound. It accounts
  for the declared gas circulation flow, contact-water hold-up, gas volume,
  elapsed condenser time, temperature and inventory-derived pressure. A 0.1 L
  hold-up is an explicit hypothetical sensitivity value, not a fitted or
  historically measured apparatus dimension.
- A standalone transfer test and no-spark collection test show stable capture
  when time is divided into smaller pieces. The coupled discharge/chemistry/
  collection computation still fails the resolution gate: at 24 hours, going
  from 16 to 32 substeps changes HCN by 52.2% primary / 46.3% charged and
  formaldehyde by 36.6% primary / 34.3% charged. A 168-hour two-substep
  replay is complete, conserves material, and is visible in Origins, but is
  **not a converged yield**. Its HCN intercept remains much lower than the
  held-out published diagnostic, and its formaldehyde regression fit is poor.
- The next numerical gate is simultaneous or demonstrably convergent
  integration of plasma fragment production, stiff neutral/charged reactions,
  and condenser extraction. The historical calorimetric energy and transport
  uncertainties remain external constraints; no literature yield is used to
  tune the internal solver.

The reproducible inputs-only dry run is stored at
`backend/data/examples/miller_urey_inputs_only.json`. It contains the starting
mixture and apparatus program but deliberately defines an empty reaction list.

## Private complete methane-set sensitivity

After the user accepted LXCat's data-use terms, a private, hash-pinned
Community/Song–Bouwman complete CH4 set was installed as an opt-in replacement
for Morgan CH4 in the mixture EEDF. The four source-labelled neutral methane
fragment channels are atom/charge checked and feed the existing energy-limited
apparatus chemistry without changing the blind product-discovery policy.
The [full sensitivity report](../backend/data/validation/prebiotic_latest/COMMUNITY_CH4_SENSITIVITY.md)
records 24-hour and 168-hour runs, 2/4/8-substep checks, local data provenance,
and the remaining benchmark mismatch. The 168-hour two-substep primary
zero-time diagnostics are 0.17169 nmol/J HCN and 11.9347 nmol/J formaldehyde;
both ledgers close, but the model energy and condenser parameters are still
hypothetical. This is not a defensible experimental-yield reproduction.

The EEDF is solved self-consistently and H2, H2O, CH4, NH3, and N2 have complete
contributor-labelled collision sets. One reviewed neutral N2 dissociation
outcome is coupled to the downstream network. The current result is still a
spatially uniform two-term plasma model, not a full discharge reactor.
Ion/electron chemistry is now represented in a propagated but correlated
alternative branch. Continuous overlapping boiler/discharge/condenser
operation has an operator-split implementation, but timestep convergence,
most excited-state heavy chemistry, pulse-resolved streamer/sheath physics, aqueous pH/activity,
transition paths, complete alpha-/beta-alanine kinetics, analytical response,
and held-out quantitative apparatus calibration remain incomplete, so
quantitative yields remain disabled by the API.

API endpoints:

- `GET /api/prebiotic/capabilities`
- `GET /api/prebiotic/validation/recorded`
- `GET /api/prebiotic/benchmarks/miller-urey-1953/profile`
- `GET /api/prebiotic/benchmarks/stribling-miller-1986/profile`
- `POST /api/prebiotic/apparatus/simulate`
- `POST /api/prebiotic/post-discovery/rank`
- `GET /api/prebiotic/post-discovery/rank/{job_id}`
- `POST /api/prebiotic/post-discovery/rank/{job_id}/held-out-evaluation`
- `POST /api/prebiotic/plasma/cross-sections/validate`
- `POST /api/prebiotic/plasma/rate-coefficient`
- `POST /api/prebiotic/plasma/eedf/solve`

The conditional rate fixture is stored at
`backend/data/examples/electron_collision_rate_fixture.json`. Dataset coverage
and missing-data status are explicit in
`backend/data/collision/prebiotic_coverage_manifest.json`.

## Scientific objective

Model a closed, recirculating gas–liquid apparatus in which water is boiled,
mixed with a reducing gas feed, exposed to a pulsed electrical discharge,
condensed, and returned to a collection flask. The experiment must predict a
distribution of products from computed chemistry rather than being supplied a
list of the historically reported products.

## Scientific firewall

Known Miller–Urey products and literature yields belong only in the held-out
evaluation data. They must not seed reaction discovery, candidate ranking, or
kinetic fitting. Apparatus conditions and starting materials are legitimate
inputs. Mechanisms from independent quantum calculations may be used when
their provenance is recorded.

## Development gates

1. Apparatus: conserve elements and charge across gas, liquid, condensate,
   walls, sampled material, and trapped nonvolatile products.
2. Discharge: represent electron-impact reactions separately from thermal
   reactions, calculate rates from cross sections plus an explicit EEDF, and
   account for deposited electrical energy. **Partial:** self-consistent EEDF
   and complete starting-gas collision coverage are implemented; sheath/
   streamer physics, excited-state populations and full heavy chemistry remain.
3. Discovery: enumerate radical, ionic, and closed-shell pathways without
   consulting the held-out product list.
4. Kinetics: solve the coupled multiphase network over repeated circulation
   cycles with uncertainty propagation.
5. Analysis: predict observable formulas, structures, stereoisomers, and
   virtual chromatography/mass-spectrometry signals.
6. Blind validation: reveal literature products only after a signed prediction
   artifact has been saved.

## Required state

- Gas and liquid temperature, pressure, and volume
- Gas flow and circulation rate
- Condenser and boiling temperatures
- Species amounts by phase
- Henry-law and vapour–liquid partition parameters
- Spark voltage, current, duration, frequency, and deposited energy
- Electron-energy distribution or a declared reduced-order surrogate
- Vessel material and optional wall-reaction model
- Sampling and transfer history

## Required provenance

Every rate or cross section records its source, method, temperature range,
pressure range, uncertainty, and applicability domain. Surrogate values are
visibly labelled and are excluded from quantitative yield claims.

## Stereochemistry

Stereoisomers remain distinct species. In an achiral apparatus, enantiomeric
pathways receive identical rate laws unless an explicit chiral influence is
present. A spontaneous non-racemic result is therefore a diagnostic warning,
not evidence of selectivity.

## Initial validation suite

- No-discharge control produces no discharge-only products.
- Zero-flow control prevents inter-compartment transfer.
- Condenser-off control does not accumulate condensate.
- Carbon, nitrogen, oxygen, hydrogen, and charge close at every time step.
- Splitting one circulation interval into smaller steps does not materially
  change the result.
- An achiral network remains racemic within numerical tolerance.
- Product-name and yield data cannot be read by discovery code before the
  prediction artifact is finalized.
