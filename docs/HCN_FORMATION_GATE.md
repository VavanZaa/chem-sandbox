# HCN formation gate (2026-09-23)

## Primary target changed — 2026-09-29

The user has retired wet Miller–Urey reproduction as the active development
goal. The primary gate is now the [dry flowing CH4/N2 discharge benchmark](DRY_FLOW_HCN_GATE.md).
This file preserves wet-apparatus research history and archived dry results;
its older positive dry amounts predate the corrected N2 product-state rule.
Use the new dry gate for current-code status and the next scientific step.

## Reduced-source route and independent atom ceiling — 2026-09-29

An opt-in [source-driven 0-D route](REDUCED_HCN_AND_ATOM_CEILING.md)
now skips the EEDF and spark geometry **only when effective fragment hazards
are supplied explicitly**. It couples GRI gas reactions to finite-rate HCN
Henry contact at declared fixed pH, rejects nonconserving/direct-C–N source
injection, records provenance, and exposes an independent C/N atom ceiling.
The UI template has zero source hazards. Synthetic source-on, zero-source,
carbon/nitrogen removal, energy-bound and solver-refinement tests pass, but
the synthetic ground-state N2 source is not physical evidence and no wet
HCN yield was validated. This is a shorter *conditional* route to an assay
comparison once independent source/contact measurements exist; it does not
resolve the state or thermal blockers below.

## Current scientific status — 2026-09-27

**The charged-branch comparison is not an independent HCN kinetics check.**
An assay-free, exact-record audit of the installed Cantera mechanisms finds
that the 440-reaction Pavan plasma mechanism contains all 325 GRI reaction
records, including all 30 HCN-containing reactions. Pavan adds plasma species
and reactions and evolves a different state, so its trajectory can differ;
nevertheless, agreement between the two branches cannot independently
validate the shared HCN rate constants. The application now exposes this
audit's counts and both installed-mechanism SHA-256 hashes in the
charged-branch result.
The earlier opt-in charged HCN endpoint trace reports *instantaneous net rates at each
solver substep endpoint*, not integrated pathway extents or product yields.
An endpoint ranking therefore cannot establish which reaction made most of
the accumulated aqueous HCN. These diagnostics read no HCN assay and do not
alter either propagating mechanism.
The focused overlap/long-horizon tests and the full backend test suite pass
after this change; the suite still emits third-party deprecation warnings.

**Integrated charged-branch pathway gate (assay-free).** An opt-in five-point
quadrature now integrates signed Pavan reaction progress over each adaptive
charged-BDF interval, with separate electron-source and coupled gas/water
terms in its HCN closure. In the hypothetical 925-minute, 192-substep Barth
feed run, it records `5.3080601e-9 mol` gross signed HCN formation and
`5.3038444e-9 mol` gross signed consumption, leaving
`4.2156767e-12 mol` observed net gas-chemical HCN; the quadrature closure
residual is `1.10e-20 mol`. The charged aqueous HCN+CN- endpoint remains
`3.9389551e-12 mol`, with closed material ledger and no failed substeps.
More than 99.9% of signed HCN formation is canceled by consumption in this
model. These are reaction-progress diagnostics, not independent chemical
evidence, molecular event counts, or measured yield. CVODES still emits
minimum-step warnings, and the complete thermal balance is not closed.
The exposed opposing-signed-reaction-contributions / absolute-net ratio is
about `2.52e3`. It is a warning that the modeled **gas-chemical net** is
the difference of much larger opposing progress terms, not a factor by
which to multiply any yield uncertainty: coupled gas-water capture and
nonlinear rate responses prevent that inference.

The dominant integrated charged HCN source is the *shared GRI* reaction
`HCNN + O <=> HCN + NO` (`4.3742e-9 mol`); the dominant sink is
`H + HCN (+M) <=> H2CN (+M)` (`5.3038e-9 mol`). In an assay-free 60-minute,
eight-substep precursor audit, `CH + N2 (+M) <=> HCNN (+M)` accounts for
the integrated HCNN formation (`1.57895e-8 mol`). The installed Pavan
mechanism contains HCNN but no NCN species. This is a serious mechanism-
applicability question: [Moskaleva and Lin's primary CH+N2 study](https://doi.org/10.1016/S0082-0784(00)80652-9)
finds a spin-conserved NCN+H channel rather than the older assumed
HCN+N(4S) channel, but its fitted rate covers **1500–4000 K**, not this
298.15 K wet flask. A [high-temperature shock-tube study](https://pubmed.ncbi.nlm.nih.gov/17958405/)
also supports NCN+H dominance over 2228–2905 K. Neither source justifies
swapping an NCN pathway into the present room-temperature solver or assigning
its yield. A [2026 primary multichannel master-equation study](https://doi.org/10.1021/acs.jpca.6c01935)
further distinguishes cyclic addition, HCNN chain addition, and HNCN
stabilization before H+NCN formation. Its reported appreciable stabilized
HNCN pertains to combustion-engine/turbine conditions, **not** a measured
298 K wet-spark product fraction. The current single-HCNN association
cannot claim to resolve those isomers. A separate assay-free comparison now
checks the **low-pressure,
Ar-bath CH+N2 association limit** against the direct 202-584 K experiment of
[Brownsword et al. (1996)](https://pubs.rsc.org/en/content/articlelanding/1996/ft/ft9969201087).
At 298.15 K, installed GRI gives `1.55885e-31` versus the experimental
fit's `1.69803e-31 cm6 molecule-2 s-1` (ratio `0.918`). This is a useful
partial rate check, **not** a validation of GRI's HCNN product assignment,
1-bar wet-gas falloff, subsequent `HCNN + O` HCN branch, or HCN yield.
At the other asymptote, [Fulle and Hippler's separate high-pressure
He-bath experiment](https://doi.org/10.1063/1.472383) yields a fitted
298.15 K limit of `4.10381e-11 cm3 molecule-1 s-1`, while installed GRI
gives `1.20997e-11` (ratio `0.29484`). This discrepancy further blocks
assuming the unmodified falloff is validated; neither asymptote alone
specifies the actual 1-bar wet-mixture rate. The source rates are not
substituted into the inventory-driving mechanism or fitted to HCN.
The charged amount still needs independently constrained mixture-dependent
and downstream low-temperature kinetics and explicit NCN/HCNN/HNCN state coverage
before it can be treated as a defensible estimate.
An [opt-in half-rate CH+N2 association stress](BARTH_LONG_HORIZON_GATE.md#ch--n2-association-partial-reference-audit-2026-09-27)
changes the charged 925-minute aqueous HCN+CN- amount from
`3.9389551e-12` to `2.7808761e-12 mol` at 192 substeps; at 384
substeps the stressed value is `2.7808041e-12 mol`. This reproducible
structural response does **not** calibrate wet-gas falloff, identify the
adduct, or turn an arbitrary rate multiplier into a physical error bar.

**Next falsifiable mechanism test:** obtain an independently usable 290–310 K,
mixture-/pressure-dependent CH+N2 association and adduct-isomer fate model,
then compare its *unfitted* effective rate against the current 1-bar wet-feed
value (`2.43836e-12 cm3 molecule-1 s-1`). The separate Ar low- and He
high-pressure measurements constrain opposite limits but do not provide
that wet-mixture curve. Keep the primary and charged propagating rates
unchanged until a candidate also passes thermochemistry, detailed-balance,
atom/charge, negative-control, and timestep checks. The HCNN+H and
HCNN+O downstream rates, spark exposure, and thermal budget are separate
requirements even if this association test succeeds.

An [independent CH+H2O competition audit](BARTH_LONG_HORIZON_GATE.md#competing-ch--water-vapour-loss-audit-2026-09-27)
adds a second source-bounded warning on the same carbon route: installed
GRI's 298.15 K water-reaction rate is `2.175` times a published
room-temperature CH-removal fit. At the declared initial feed, the model's
water-vapour CH-loss hazard is `0.243` of its CH+N2 hazard; its one-hour
signed-net CH trace shows comparable competition. The independent H-atom
product study does not resolve the CH2O-versus-CO split, so neither the
GRI product assignment nor an HCN correction is established. The
default propagating chemistry and held-out assay remain untouched. An
opt-in charged-shadow-only 0.45974 rate-factor stress (Blitz/GRI at 298 K)
raises its 925-minute aqueous HCN+CN- amount from `3.9389551e-12` to
`4.6957446e-12 mol` at 192 substeps; the stress gives
`4.6956284e-12 mol` at 384 substeps (`~2.5e-5` relative difference).
The unchanged primary branch, closed ledgers, zero failed steps, same-
process baseline restoration, and carbon/nitrogen negative controls make
this a useful *model-sensitivity* result, but the pressure extrapolation
and unmeasured product split prevent calling either amount a corrected
prediction or an experimental uncertainty bound. CVODES minimum-step
warnings persist. See the linked audit for its assumptions and traces.

**Current contact-domain correction:** The older Barth low-contact table
that called `kLa=1e-5 s^-1` source-bounded predates the N2 electronic-state
correction. In the current 925-minute model, that case reaches ideal
pH `4.71706`; the strict Winkelman pH-5-to-7 mode rejects it at pH
`4.9965` during cycle 4. The flagged extrapolation mode gives charged
aqueous HCN+CN- `1.5028622e-12` and `1.5028469e-12 mol` at 192/384
substeps, with closed ledgers and no assay access. These are stable
conditional outputs, **not** a source-bounded or measured yield. A
regression test now guards the strict rejection. The older table is marked
archived in the [Barth gate](BARTH_LONG_HORIZON_GATE.md#contact-rate-and-kinetic-domain-sensitivity).

The opt-in full-duration HCNN fate diagnostic sharpens that priority. Of
`3.6597240e-7 mol` integrated **signed net** HCNN consumption in the charged
mechanism, `85.040%` follows `H + HCNN <=> CH2 + N2`, `13.147%` follows
`HCNN + O <=> CO + H + N2`, and only `1.195%` follows
`HCNN + O <=> HCN + NO`; other routes account for the remainder. This
partitions the simulated mechanism's net progress, not microscopic forward
branching probabilities or measured selectivity. A clearly hypothetical
half-rate stress on the dominant `H + HCNN` sink changes the 925-minute
charged aqueous HCN+CN- amount from `3.9389551e-12` to
`6.2140679e-12 mol`, with closed material ledgers, no failed charged
substeps, and an unchanged primary branch. That response is a structural
sensitivity only, not a physical uncertainty bound or a better yield.
Doubling the full-duration resolution from 192 to 384 substeps leaves the
three HCNN net-fate fractions at `85.0403%`, `13.1475%`, and `1.19523%`,
respectively, versus `85.0403%`, `13.1475%`, and `1.19523%` at 192 steps
to the shown precision. The charged aqueous endpoint changes from
`3.9389551e-12` to `3.9388559e-12 mol` (`~2.52e-5` relative). Both
traces reconcile, both ledgers close, and neither assay is loaded. This
supports numerical stability of the *current mechanism*; CVODES still
emits minimum-step warnings and no physical parameter is thereby validated.
An opt-in [one-hour charged carbon-route audit](BARTH_LONG_HORIZON_GATE.md#one-hour-charged-carbon-route-audit-2026-09-27)
now traces the model's CO2+-to-CO, HCO, CH2O, CH2(S), CH, HCNN and HCN
links, with external electron-source extents kept separate from gas
reaction contributions. The resulting graph contains return cycles and
shared GRI/Pavan kinetics; it is a mechanism-internal dependency map,
not a measured molecular lineage or product probability. Its CH2O trace
also found and repaired a *diagnostic-only* missing methanediol term in
the gas-reaction closure; no species inventory or yield was changed.

Two explicitly **hypothetical, charged-branch-only** rate stresses quantify
this dependency without reading the assay or mutating the primary GRI
branch. At 925 minutes/192 steps, halving the Pavan/GRI
`H + HCN (+M) <=> H2CN (+M)` rate changes charged aqueous HCN+CN- from
`3.9389551e-12` to `7.3136207e-12 mol`; doubling
`HCNN + O <=> HCN + NO` changes it to `7.5322252e-12 mol`. Both cases
retain closed charged material ledgers, zero failed substeps and unchanged
primary HCN. The controls are **not measured rate intervals**, confidence
bounds, calibrated alternatives, or evidence that either stressed result is
closer to the real assay. A same-process one-hour baseline/stress/baseline
test confirms the cached Pavan phase restores default rates after a stress.
One-hour carbon-removed and nitrogen-removed controls under the half-sink
stress both return zero charged aqueous HCN+CN- and close both ledgers.

**Formation is not presently validated by the primary solver.** The QST/Tabata
N2 electron-impact table reports total `N + N` production, not a branching
fraction for two ground-state N(4S) atoms. The underlying crossed-beam study
finds N(2D)+N(4S) to be the primary mechanism over its measured energy range
([Cosby 1993](https://doi.org/10.1063/1.464385)); the
[Tabata 2006 compilation](https://doi.org/10.1016/j.adt.2006.02.002)
adopts Cosby's total production cross section. The model previously mapped
*all* these outcomes to `N + N` in GRI and Pavan. That mapping silently
erased electronic state and altered subsequent low-temperature chemistry.
The Origins UI now labels the post-EEDF N2 rate as a **total dissociation
outcome with unresolved electronic products**, rather than displaying it as
an injected `N2 -> 2 N` ground-state source. Its source table likewise
distinguishes atom/charge balance from electronic-state resolution.

The collision rate remains available as an EEDF-folded **total outcome** and
its threshold is still kept separate from the overlapping Phelps energy-loss
manifold. Until a justified state-selective branching and N(2D) fate network
are available, neither the primary nor charged branch injects that outcome as
two ground-state atoms. The one-cycle, two-substep CH4/N2/H2/H2O formation
check now returns `formation_gate_not_passed`: primary net HCN is
`8.56194e-17 mol`, below its numerical reporting floor. Its separate charged
shadow retains `9.27093e-11 mol` total HCN+CN- under other pathways, but is
not a quantitative yield or a passing primary gate. Spark-off, methane-free,
and nitrogen-free controls remain below the floor; both material ledgers
close. **All positive HCN quantities and "passes" below this section are
archived results of the superseded ground-state assignment, not predictions
of the current code.** No HCN assay was opened or used to set this correction.

[Peterson et al.'s N2+ dissociative-recombination experiment](https://ir.amolf.nl/pub/4894)
measured *zero-collision-energy, v=0 ion* product branches, not product
branches of direct electron-impact dissociation of neutral N2. A
[2023 N2 collision review](https://doi.org/10.1140/epjd/s10053-023-00687-5)
also lists recommended atomic-state fractions for that **ion
recombination** process. Importing either set into the neutral-dissociation
channel would repeat the state-assignment error. The review reports Cosby's
ground+ground neutral-dissociation channel as insignificant at 50 eV while
emphasizing that neutral product branching remains poorly measured. None of
these results supplies an EEDF-wide neutral branch law.

The separate [Barth thermal-capacity screen](BARTH_LONG_HORIZON_GATE.md#boundary-and-method)
compares the declared, unmeasured electrical estimate against the IAPWS-95
sensible enthalpy of the initial 50 mL of water. Its `5.23` ratio is a
conditional warning that thermal storage/transfer matters, not a measured
temperature or a closed heat balance.

Next: obtain independently measured or defensibly bounded N2 dissociation
product-state branching versus electron energy, add N(2D) excitation/
quenching/reactive kinetics with state and energy conservation, and rerun
the control/convergence gate before any experimental-yield comparison. In
parallel, vet the low-temperature HCN formation/consumption rates against
sources independent of the shared GRI/Pavan records; a second trajectory
using those same records cannot supply that check.

The source-energy applicability audit now partitions the EEDF-folded total
N2 dissociation event rate against Cosby's **18.5–148.5 eV** fragment-state
measurement window, without assigning its products. In the hypothetical
200 Td Barth initial feed, `56.014%` of the calculated rate comes from
electron energies **below 18.5 eV** and `43.986%` within the measured
window. Extending the EEDF grid from 60 eV/400 cells to 100 eV/700 cells
and 150 eV/1000 cells changes the total outcome rate by under `0.23%`
and the below-window fraction by under `0.005` percentage points. The
150 eV grid covers the full upper measurement endpoint; its calculated
rate above 148.5 eV is negligible. Across the 60-minute/48-step and
925-minute/192-step hypothetical wet-flask runs, opportunity-weighted
below-window fractions are `56.008%` and `55.975%`, respectively. These
fractions are **electron-energy support**, not N(2D), N(4S), or HCN
branching. Cosby's observed dominant N(2D)+N(4S) channel cannot be
extrapolated to the majority below-window part without independent evidence.
The [Walter–Cosby–Helm state-selective photofragment study](https://doi.org/10.1063/1.466152)
does measure atomic yields from *specific prepared N2 singlet rovibrational
levels*, so it is a useful candidate for a future
state-resolved intake. It does not itself give the electron-impact
excitation probabilities into each level. The current Phelps electron
collision set lumps this region into an `N2(SUM)` excitation at 13 eV;
that aggregate cannot be multiplied by the level-specific photofragment
fractions without an independently justified state-to-state mapping.
The assay-free `backend/scripts/audit_n2_source_energy_support.py` reproduces
the grid and 50/100/150/200 Td field checks from the feed card. At those
fields, the calculated total N2 outcome rate coefficients are respectively
`1.14003e-20`, `6.81885e-18`, `6.84094e-17`, and `2.37539e-16 m3 s^-1`;
the rate fractions below the product-state measurement window are 99.815%,
90.506%, 72.274%, and 56.014%. The approximately `2.08e4`-fold rate
range is sensitivity to an **assumed**, independently unmeasured reduced
field, not a measured uncertainty band. All four EEDF solves converged;
none determines fragment branching or HCN yield.
The same assay-free scan now includes a non-propagating source-form check
that evaluates the
[QST/Tabata official F3 fit](https://www.qst.go.jp/uploaded/attachment/8709.f)
against the currently imported
[recommended numeric table](https://www.qst.go.jp/uploaded/attachment/8710.dat)
on the same 200 Td initial-feed EEDF. The analytic fit folds to
`2.18626e-16 m3 s^-1`, versus `2.37539e-16 m3 s^-1` for linearly
interpolated tabular points (the latter is `8.65%` higher relative to the
fit). The analytic fit places `52.490%` of that rate below Cosby's
18.5 eV product-state measurement limit, versus `56.014%` from the
tabular interpolation. Thus the evidence gap survives this independent
source-form check, but its exact numerical share is representation-sensitive.
The fit and sparse table are distinct representations of the same
total outcome, not independent experiments or product-state evidence.
At the original tabulated points the fit is `31.5%` above the 14 eV
recommended value and `21.7%` below the 16 eV value, so replacing the
table solely because the fit is smooth would not be a source validation.
Across the 925-minute, 192-step feed-only run, tabular and fitted
*uncoupled outcome opportunities* are `0.00493504` and `0.00454248 mol`;
their opportunity-weighted below-window shares are `55.975%` and
`52.449%`. These are possible electron-collision events before energy
allocation and state-resolved chemistry, **not** formed N atoms or HCN.
The solver still uses the hash-verified table; no HCN rate was adjusted.
Before changing that representation, compare the energy-resolved fit to
the original measured cross sections. Keep this model-form sensitivity
separate from EEDF/field uncertainty.

A new **non-propagating energy-access check** uses the
[N2 dissociation limits](https://doi.org/10.1039/D3CP04854C): the first
neutral two-atom channel with *no* ground-state N(4S), N(2D)+N(2D), opens
at 14.526 eV above ground-state N2. In the nominal ground-target model,
events below that electron energy must have at least one N(4S) atom if the
imported total neutral-dissociation cross section is valid there. At the
hypothetical 200 Td initial feed, `11.098%` of the tabular EEDF-folded
total rate falls below 14.526 eV (`12.267%` for the read-only analytic
fit); the 60, 100 and 150 eV EEDF grids agree on the tabular fraction
within `0.007` percentage points. Over the 925-minute/192-step run, the
opportunity-weighted fractions are `11.088%` and `12.256%`, respectively.
The 50/100/150/200 Td initial-feed tabular fractions are `82.556%`,
`38.126%`, `19.122%`, and `11.098%`, illustrating sensitivity to an
unmeasured field. These are **conditional energy-access fractions**, not
measured N(4S) branching or atom yields: the imported near-threshold
Tabata total is below Cosby's 18.5 eV fragment-state measurement range,
the companion atom state and kinetic-energy distribution remain unknown,
and target internal excitation can shift the nominal threshold. No atom
source or HCN yield is changed by this audit.

### Conditional N(2D) fate in the wet feed (2026-09-27)

A new **read-only**, near-298 K first-collision audit uses the evaluated
bimolecular N(2D) removal rates and multiplicative error factors in
[Herron 1999, Tables 2-3](https://doi.org/10.1063/1.556043). It is attached
to EEDF snapshots but does not assign a branch to the total N2 dissociation
cross section or inject N(2D) products. For the Barth `220725001` initial
reported gas fractions at the explicitly hypothetical 1 bar, 298.15 K and
0.95 L headspace, the conditional lifetime of a **thermalized** N(2D) atom
is `5.89115e-8 s`. Its normalized first-collision hazards are H2O `91.593%`,
H2 `2.969%`, CO2 `2.942%`, N2 quenching `2.231%`, and O2 `0.265%`.
Independent extremes of the evaluator's multiplicative rate factors move
the H2O hazard fraction from `79.287%` to `96.861%`; these are **not** a
statistical confidence interval. The original HCN assay remains unopened.

The same evaluated source gives `k(298 K) = 1.4e-12 cm3 molecule^-1 s^-1`
with a factor-two error for N(2D)+O, which transfers the excitation to
`O(3P,1D)` without a resolved oxygen-state split. This channel is now
counted in the *read-only* removal audit; no oxygen state is injected into
the chemical solver. The tabulated expression was assessed over 300–400 K,
so use at 298.15 K is a two-kelvin extrapolation of an explicitly evaluated
298 K value, not an independently measured wet-spark rate.

For N(2D)+H2O, Herron also distinguishes a 1976 flash-photolysis estimate
of `2.5e-10` from a more sensitive 1998 laser-induced-fluorescence value
of `4.2e-11 cm3 molecule^-1 s^-1`. The recommended `4.0e-11` is based on
the latter, with a factor-two error reflecting water partial-pressure
measurement difficulty. Replacing it with the older larger number would
not be an evidence-based sensitivity correction. Neither measurement
provides quantitative product yields.

This changes the next scientific priority: even a resolved N2 -> N(2D)+N(4S)
source would mostly encounter water in this *initial wet feed*. The evaluated
H2O rate is a **total N(2D) removal** rate; its product branches are not
quantified by that number. Nascent N(2D) may react before thermalization,
and wall/radiative loss, evolving gas composition, and spark timing are not
captured by this first-collision audit. Do not turn the 91.6% into a water
reaction extent, HCN yield, or a ground-state N source. The next kinetic
intake must resolve N(2D)+H2O products and energy flow, alongside the
electron-energy-dependent N2 product-state branching.

A [2014 crossed-beam/QCT study](https://doi.org/10.1021/jz501757s)
identifies competing `HNO + H` and `HON + H` exits from N(2D)+H2O;
its [supporting methods](https://pstorage-acs-6854636.s3.amazonaws.com/3613440/jz501757s_si_001.pdf)
report collisions at 7.5 and 12.2 kcal mol^-1 under single-collision
beam conditions and QCT HON/(HNO+HON) values of 0.41 and 0.52,
respectively. Those ratios are **not** 298 K thermal-gas branching data,
do not cover all exits, and cannot be applied to the evaluated total
removal rate or this wet, approximately 1-bar reactor. They constrain
which product identities deserve a state-resolved kinetic intake, not
the current inventory, HCN amount, or an uncertainty interval.

In a 60-minute, 48-step wet-flask rerun, only the **initial** N(2D)
first-collision snapshot has complete evaluated bath coverage. The other
47 snapshots fail closed because the evolving gas includes H, OH, HO2,
or H2O2 without an evaluated rate in this audit. The largest combined
unreviewed mole fraction is `1.30877e-4`; H2O2 reaches
`2.63650e-7 mol` in the gas. At the most competitive modeled step, an
N(2D)+H2O2 rate of `5.34618e-9 cm3 molecule^-1 s^-1` would equal the
**sum** of the currently reviewed removal hazards. That is a calculated
competition threshold, **not** a measured rate or physical bound. A
targeted search found no sufficiently verified N(2D)+H2O2 rate/product
entry to install, so later-step fate remains unresolved rather than being
renormalized over the known gases. The same rerun retains closed material
ledgers and no injected total-N2 neutral fragment.

The reported-duration 925-minute, 192-step rerun extends the same audit:
only its initial snapshot has complete reviewed N(2D) bath coverage;
191 subsequent snapshots remain unavailable. The largest combined
unreviewed mole fraction is `1.92807e-4` after adding the evaluated
atomic-O removal rate. The remaining unreviewed partners are H, OH, HO2,
and H2O2. At the most competitive modeled
step, an N(2D)+H2O2 rate of `3.65007e-9 cm3 molecule^-1 s^-1` would equal
the sum of reviewed removal hazards. This is again an **equal-hazard
threshold**, not a measurement, an upper or lower bound, or an inferred
branching fraction. Primary aqueous HCN+CN- is `4.18575e-27 mol` and the
charged shadow is `3.93896e-12 mol`; both material ledgers close. Neither
result establishes a real HCN yield or authorizes using the held-out assay.

## Scope and result

The first gate asks whether an unseeded, reducing CH4/N2/H2/H2O discharge
scenario forms HCN in the simulated inventory, with atom balance and simple
negative controls. It does **not** require reconstructing a specific historical
apparatus and does **not** assert a measured or calibrated yield.
"Unseeded" here means no HCN was placed in the initial inventory. The
published GRI-Mech gas network already contains HCN species and elementary
reactions, so this formation gate is **not de novo product-identity discovery**
from first principles. Known HCN measurements and product yields remain held
out of the solver and are used only downstream for validation; the distinction
between a general pre-existing kinetic mechanism and a product-blind search
is explicit in the result metadata.
The gas-network outputs now report that known product identities are present
in their published GRI/Pavan mechanisms; only the separate EEDF collision
calculation is product-identity blind. Neither statement means an HCN recipe
or a measured yield was supplied as a starting inventory or fitted target.
The apparatus-wide `scientific_firewall` now also marks the Boltzmann/GRI
simulation as **not blind-discovery eligible**, flags pre-existing HCN
reactions as accessible to the gas solver, and emits a user-visible warning.
A frozen blind-discovery manifest alone cannot override a published mechanism
that already contains the target species and reactions.

The pre-carbonate one-hour and seven-day, two-substep control sweeps pass.
The one-hour five-case sweep also passes after carbonate coupling. With the
current NH3 transport extension, the 24-hour and seven-day five-case sweeps
now pass at two substeps. Seven-day spark-on checks pass at four and eight
substeps, and the three negative controls pass at four and eight. The
half-methane diagnostic also passes at four substeps. These separately run
subsets are not newly passed combined five-case sweeps at those resolutions.
With the opt-in, complete
IST-Lisbon H2 and Community/Song-Bouwman CH4 collision sets at the declared
200 Td, 400 K scenario point, the primary branch makes 2.668e-7 mol
HCN-plus-CN- equivalents;
2.606e-7 mol of HCN plus CN- is in water or collected condensate. This is
9.92e-6 mol cyanide-family product per
mol feed-carbon atoms. The no-spark, no-methane, and no-nitrogen controls are
below their explicitly reported numerical floors; all material ledgers close.
With half as much methane and the total starting gas amount held constant by
replacement with H2, HCN is 1.429e-7 mol (0.536 times the full-methane case).
That is an internal response check, not a fit to an experiment.

The control sweep exposed a BDF roundoff defect in the carbon-free case: a
carbon-containing intermediate could become slightly negative even though
the initial gas and injected fragments contained no carbon. The gas ODE now
keeps all species containing an absent element exactly zero. It does not add
an HCN-producing reaction or change a present-element feed.

The later fixed-volume transport audit found that the step-equilibrium and
flow-limited condenser paths could use a total-pressure snapshot taken before
condensation changed the headspace inventory. Both now calculate each gas
partial pressure from its **current** amount and the declared gas-loop volume
using the ideal-gas relation. The supplied pressure remains in diagnostics;
it no longer silently changes Henry capture at fixed volume. Direct tests
check material conservation, pressure-snapshot invariance, and stronger
partitioning at smaller gas volume. That earlier transport audit did not
include pH; the opt-in aqueous acid-base coupling below is a subsequent gate.

Run the check from `backend` with the installed environment:

`..\.venv\Scripts\python.exe -m scripts.run_hcn_formation_gate`

## Duration and numerical resolution checks

Run `..\.venv\Scripts\python.exe -m scripts.run_hcn_formation_gate
--resolution-sweep` from `backend` to repeat the same five input controls at
2, 4, 8 and 16 substeps per one-hour cycle. The 2026-09-23 sweep returned:

| Substeps | Net HCN + CN- formed (mol) | Aqueous HCN + CN- (mol) | Formation change from next-finer run |
| ---: | ---: | ---: | ---: |
| 2 | 2.6681906e-7 | 2.6059274e-7 | 0.01148% vs 4 |
| 4 | 2.6678844e-7 | 2.6056357e-7 | 0.00600% vs 8 |
| 8 | 2.6677244e-7 | 2.6054829e-7 | 0.00305% vs 16 |
| 16 | 2.6676431e-7 | 2.6054052e-7 | not tested further |

All no-spark, no-methane and no-nitrogen cases remained below their numerical
reporting floors, and all five material ledgers closed at every resolution.
This establishes stability at the **tested one-hour duration and scenario**
only. It does not establish a converged seven-day yield, measured uncertainty,
or experimental magnitude agreement.

The same five-case gate was separately run for six one-hour cycles at 2, 4, 8,
and 16 substeps per cycle. All cases retained closed material ledgers; the
no-spark, no-methane, and no-nitrogen controls stayed below their reporting
floors. The six-hour endpoint was:

| Substeps | Net HCN + CN- formed (mol) | Aqueous HCN + CN- (mol) | Formation change from next-finer run |
| ---: | ---: | ---: | ---: |
| 2 | 1.9385148e-6 | 1.9305220e-6 | 0.01749% vs 4 |
| 4 | 1.9381758e-6 | 1.9301852e-6 | 0.00879% vs 8 |
| 8 | 1.9380055e-6 | 1.9300160e-6 | 0.00440% vs 16 |
| 16 | 1.9379202e-6 | 1.9299313e-6 | not tested further |

This is a within-model six-hour resolution check, not a full-duration
convergence result or an experimental yield validation. The endpoint rises
more than sixfold from the one-hour amount, so extrapolating a one-hour rate
to the week-long apparatus is not justified.

A 24-hour five-case sweep at 2, 4, and 8 substeps also passed its input
controls. These reruns use a tightened per-conserved-quantity atom/charge
tolerance (the primary branch uses a maximum of 1e-12 mol absolute and
1e-12 relative to that quantity's expected amount), not a threshold scaled
to the large water pool.
Both the primary and alternative charged branch passed their material ledgers.
The largest primary/charged residuals in the 24-hour, 2-substep spark-on case
were 1.78e-14/5.33e-14 mol; at 8 substeps they were 1.07e-14/5.33e-14 mol.
The tightened check also exposed approximately 9.69e-10 mol N and
4.85e-10 mol charge drift in an older one-hour coupled-flow/Pavan case using
the earlier collision set. Tightening that branch's BDF relative/absolute
tolerances from 1e-6/1e-18 mol to 1e-9/1e-21 mol brought element and charge
residuals below 1e-11 mol in its 2- and 8-substep regression. Its global
per-quantity acceptance floor is now 1e-11 mol, about two orders tighter
than that detected drift; the regression also requires zero failed charged
substeps. A 24-hour/2-substep spark-on rerun with the tighter integrator
retained a closed charged ledger (6.52e-12 mol residual) and no failed
substeps. The later carbonate-coupled 24-hour *spark-on* case has now been
rerun at 2, 4 and 8 substeps with the final cation projection; full
five-case controls at those resolutions have not yet been repeated.

| Substeps | Net HCN + CN- formed (mol) | Aqueous HCN + CN- (mol) | Formation change from next-finer run |
| ---: | ---: | ---: | ---: |
| 2 | 8.6219142e-6 | 8.6139162e-6 | 0.02154% vs 4 |
| 4 | 8.6200577e-6 | 8.6120620e-6 | 0.01076% vs 8 |
| 8 | 8.6191302e-6 | 8.6111356e-6 | not tested further |

The 8-substep 24-hour scenario forms 0.1663 nmol HCN-plus-CN- per declared
electrical joule. The joules are a hypothetical electrical deposition
estimate during the discharge stage, not calorimetry or a closed thermal
energy balance. A 24-hour resolution result is still not a numerically
converged seven-day trajectory or a matched experimental comparison.
For longer endpoint refinement, `--case spark_on` runs only the spark-on
scenario and labels the result `partial_formation_diagnostic_only`; it never
passes the five-case control gate. The charged-branch cycle report now retains
every substep's solver status so an earlier ion-chemistry failure cannot be
hidden by a successful final substep. Complete negative-control runs must
still be repeated at the duration and resolution used for any gate claim.

The first 168-cycle (seven-day) **spark-on-only** diagnostic at 2 substeps
forms 4.51573e-5 mol HCN-plus-CN- (4.51532e-5 mol aqueous), equivalent to
0.12444 nmol per declared discharge joule. Primary and charged material
residuals are respectively 3.20e-13 and 2.63e-13 mol, and no charged substep
failed. At 4 substeps, the same spark-on scenario forms 4.51517e-5 mol
(4.51477e-5 mol aqueous), 0.01238% below the 2-substep net amount. Its
primary/charged ledger residuals are 6.04e-14/5.33e-13 mol, with no failed
charged substeps. At 8 substeps, it forms 4.51489e-5 mol (4.51449e-5 mol
aqueous), another 0.00618% change. The 8-substep primary/charged residuals
are 2.49e-14/6.68e-13 mol and no charged substep failed. This three-resolution
spark-on result establishes within-model endpoint stability for the tested
scenario, but by itself was not a full-week input-control gate, a measured
yield, or a matched experiment.

The complete seven-day, 2-substep **five-case** gate has since passed. The
full-feed case forms 4.51573037e-5 mol net HCN-plus-CN- and retains
4.51532481e-5 mol aqueous HCN-plus-CN- at the endpoint. Spark-off and
no-nitrogen form zero; no-methane produces 5.45e-33 mol, far below its
1e-12 mol numerical reporting floor. Halving methane while replacing it
with H2 gives 2.25580268e-5 mol, 0.49954 times the full-feed result. All
five primary material ledgers close. Every spark-on charged branch also
closes its material ledger and reports **zero failed substeps**. The largest
primary and charged residuals for the full-feed case are 3.20e-13 and
2.63e-13 mol. Its cumulative conservative ion clipping projection is
6.52e-11 mol, separately reported as a numerical correction rather than
chemical production. The full-feed endpoint is 0.12444 nmol net cyanide
family per **declared electrical joule**, not a measured/calorimetric
yield. Thus the full-week input-control and tested endpoint-resolution
checks for both propagated branches now pass, but the experiment-matched
magnitude, product-blindness,
aqueous kinetics and complete thermal-energy gates do **not**.

The gate now exposes the propagated ion/excited-state Pavan branch's final
HCN-plus-CN- inventory alongside the primary neutral GRI branch and requires
its own spark-on amount above the reporting floor and no-carbon/no-nitrogen
controls below that floor. At one hour and 2 substeps, it retains
2.12703e-7 mol total (2.07507e-7 mol aqueous), versus the primary
2.66819e-7 mol total: a charged-to-primary final-amount ratio of 0.7972.
Its no-methane/no-nitrogen artifacts are respectively 1.13e-22 and
1.35e-20 mol, both far below the 1e-12 mol-or-higher numerical floors;
the spark-on branch has no failed substeps and its material ledger closes.
This ~20% model-branch difference is **mechanism sensitivity**, not an
experimental uncertainty bar or two independent votes: Pavan extends GRI,
and both branches share the same collision data. The resolution-sweep
report now compares charged-branch HCN inventories too, so neutral-only
endpoint convergence cannot mask charged-branch drift.
The one-hour charged final totals at 2/4/8/16 substeps are respectively
2.1270258e-7, 2.1266758e-7, 2.1264910e-7 and 2.1263964e-7 mol; the
8-to-16 relative change is 4.45e-5. All five control cases pass at each
resolution. At the seven-day, 2-substep endpoint, the charged branch
retains 3.3342447e-5 mol, 0.7384 times the primary 4.5157304e-5 mol;
its no-methane/no-nitrogen controls are 4.24e-21/8.79e-15 mol, both below
their reporting floors, with zero failed charged substeps and a closed
material ledger. This strengthens the full-week input-control result, but
the 4- and 8-substep charged-branch runs were spark-on-only diagnostics,
not separate full-control sweeps. Their final HCN-plus-CN- totals are
3.3337428e-5 and 3.3334921e-5 mol, versus 3.3342447e-5 mol at 2 substeps;
the adjacent relative changes are 0.01506% and 0.00752%. All three charged
runs close their material ledgers and report zero failed substeps. This is
within-model seven-day endpoint refinement, not a measured uncertainty or
experimental validation. The spark-off negative control runs the primary no-discharge
apparatus; the alternative Pavan charged branch is not instantiated when the
spark is off, so its own spark-off output is not separately claimed.

The long-run charged audit found a numerical N2+ clipping defect. At one
72-hour failing point, the raw Pavan solution conserved atoms within
2.42e-15 mol, but setting a -1.28e-12 mol N2+ concentration to zero without
compensation created 2.56e-12 mol of apparent N and corresponding charge.
The branch now allows only tiny negative molecular cations (at most 2e-11
mol per event) to be projected to zero while adjusting the matching neutral
parent and electron together; this exactly preserves atoms and charge and is
reported separately from chemical reaction extents. The 72-hour/2-substep
and seven-day/2-substep spark-on diagnostics complete with no skipped charged
substeps, cumulative projection 3.71e-11 and 6.52e-11 mol respectively,
and strict whole-apparatus material closure. A larger negative excursion still
stops that charged substep and prevents a full-gate claim. These corrections
do not close the apparatus thermal energy budget.

## Separate dry-flow source-to-product probe

`..\.venv\Scripts\python.exe -m scripts.run_dry_flow_hcn_probe` runs an
unseeded **CH4/N2 gas packet** through the same EEDF and neutral GRI chemistry
used by the wet apparatus, without importing the 2024 paper's HCN mechanism or
measured yield into the solver. The paper supplies the nominal 1.57 cm3
reactor volume, 0.47 s residence point and 28.1 W mean discharge-power point.
The probe assumes 1% CH4, 1 atm, a uniform 555.15 K gas temperature inferred
from a reported 282 C outlet condition, and full-volume exposure. Its 200 Td
field and 1e20 m-3 electron-density ceiling are **hypothetical**; the latter is
reduced by the electron-drift power bound. The dry model uses 16 repeated
bounded discharge/chemistry steps, not an inferred plasma exposure time from
the paper's product fit.

The dry probe now consistently assigns the hot gas packet its ideal-gas
inventory, residence time and share of average power. With an explicitly
**assumed** 298.15 K sccm reference, the paper's 1.57 cm3/0.47 s nominal
point implies about 200.4 sccm and a 555.15 K, 1-atm packet residence of
0.25242 s. The packet receives 7.09298 J at 28.1 W, not the previous
13.207 J obtained by applying the nominal time to the smaller hot packet.
The former 16-step HCN amount of 9.4774e-8 mol and all dry-probe values
below it in earlier revisions were based on that inconsistent packet and
are superseded. The corrected 16-step impulse probe produces 5.35938e-8
mol gas HCN (0.1555 mol per mol inlet carbon; 7.556 nmol per declared J).
Its atom and allocated-collision upper-bound ledgers close. Spark-off,
no-methane and no-nitrogen controls stay below 1e-12 mol HCN. Impulse
substepping remains **not numerically converged** and should not be used for
a validation claim.

An alternative `--source-forcing inventory_limited` now integrates plasma
fragmentation simultaneously with the neutral reaction network. Its source
rate falls with the available parent molecule instead of removing a fixed
amount per second after that molecule has been consumed. This repairs the
negative-CH4 failure of the earlier fixed-rate continuous source. At the
same corrected flow scenario point, it predicts 4.64347e-8 mol gas HCN
per packet (0.13473 mol per mol inlet carbon). The maximum elemental residual
is 6.31e-20 mol, and 0.03635 J resolved fragment energy stays below
7.08285 J allocated collision energy and 7.09298 J declared discharge
energy. The
spark-off, no-methane and no-nitrogen controls remain below 1e-12 mol HCN.
This is **numerical** step refinement of a fixed-snapshot mean-field model,
not convergence with evolving EEDF, microdischarge exposure, reactor flow or
spatial temperature. The paper's measured HCN and published HCN mechanism
remain held out; no experimental magnitude agreement is claimed.

The `--source-forcing evolving_eedf` mode recomputes the EEDF and
collision-source rates from the **current** gas inventory at each bounded
time slice, then integrates inventory-limited source terms with the neutral
reaction network within that slice. In the corrected 200 Td flow scenario,
64 and 128 slices produce 4.6855171e-8 and 4.6857953e-8 mol gas HCN;
the relative endpoint change is 5.94e-5. At 128 slices, the elemental
residual is 1.89e-19 mol; 0.03629 J of resolved neutral fragmentation
remains below 7.06330 J allocated collision energy and 7.09298 J declared
discharge energy. The collision-set-covered fraction falls from 1.0000 to
0.9957
as products accumulate. All three feed/discharge negative controls stay
below 1e-12 mol HCN. These are **within-model** numerical checks only.

Field sensitivity is still large: with the corrected 64-slice method and
other scenario controls fixed, 189, 200 and 435 Td give respective HCN
fractions of feed carbon 0.1205, 0.1360 and 0.4146. The paper
reports 189 and 435 Td as average experimental field estimates, but used
300 Td in its fitted kinetic model; that fitted choice is not a blind
validation input here. The same study inferred about 1.3 ms of active
plasma exposure from its product-model comparison versus 0.47 s nominal
gas residence, indicating sparse streamers. Our current scenario exposes
the whole hot packet uniformly for the assumed 0.2524 s, so even a similar HCN
fraction would not establish independent quantitative reproduction. The
next gate is an independently constrained streamer exposure/flow model and
a held-out HCN measurement at a matching feed, sampled gas phase, energy,
and time basis; the wet gas--water apparatus remains a separate target.

The corrected flow/power basis is internally consistent, but does not show
that the published residence convention is wrong or resolve the real
temperature field. The probe accepts `--flow-reference-temperature-k` for a
declared reference-temperature sensitivity; it keeps the experimental
comparison flag false. Neither its HCN fraction nor its energy cost is
accepted as a matched measurement until the actual flow reference,
microdischarge exposure and gas-temperature profile are constrained.
The [primary 2024 methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10980260/)
state that the feed entered at ambient temperature and pressure, the
mass-flow controllers supplied 50-200 sccm, and the nominal residence was
computed from the approximately 1.57 cm3 reactor volume and flow. They do
not specify the standard-temperature convention in the accessible methods;
the probe's 298.15 K reference remains an assumption, not a reported value.
The 2026-09-23 targeted recheck of that primary article also confirms an
inline gas-chromatography product assay, 22 kHz operation with nanosecond
streamer-current pulses, and an outlet-adjacent temperature measurement used
to approximate the reactor average. The authors treated reduced field and
electron density as ranges because their experimental estimates were not
accurate enough to fix them; the active streamer exposure used in their
kinetic comparison was inferred by fitting products. These facts support the
gas-phase HCN observable but do **not** independently fix the exposure of
each gas parcel. The article links a supplementary PDF, but that file was
not retrievable through the available PMC route during this recheck; no
unread supplementary number was imported. This source is therefore still
unsuitable for a claim of independent magnitude agreement with the current
full-volume dry probe.

## Partial thermochemical energy inventory audit

The gas-only dry probe now evaluates the initial and final ideal-gas internal
energy from the **installed GRI species thermochemistry**, using Cantera's
partial molar internal energies in J/kmol and the actual species amounts in
kmol. This covers every nonzero dry-probe gas species but not excited-state
storage, walls, photons or a measured heat flux. At 200 Td with 128 evolving
EEDF slices, the corrected packet's modeled gas internal energy changes
from 0.076975 to 0.102786 J: +0.025811 J. The declared electrical input is
7.092979 J, of which 7.063298 J is allocated to imported collision
opportunities and 0.036292 J appears as resolved neutral-fragment source
cost. If **all** declared electrical input were deposited into a fixed-volume
modeled gas, the first-law remainder would be -7.067167 J of non-electrical
exchange into that gas. It is an **inferred, conditional remainder**, not a
measured or mechanistically allocated heat/radiation loss. The dry probe's
`complete_thermal_energy_balance_closed` flag remains false.

The wet apparatus now reports the corresponding GRI gas-stage inventory
change separately. Its one-hour, two-substep spark-on scenario gives
+25.1013 J across neutral gas stages, with all nonzero gas species covered
by the GRI thermochemical data. That sum includes coupled gas-water transfer
and therefore cannot be interpreted as chemical reaction heat or subtracted
from the electrical input as if the liquid, phase change and external
heating/cooling did not exist. Closing the full apparatus first law still
requires liquid/ionic enthalpies, latent heat, flow/work, wall/radiative
exchange and a constrained plasma-to-gas energy deposition model. The
[Cantera thermodynamic-property API](https://cantera.org/3.2/python/thermo.html)
specifies the partial molar internal-energy unit used in this audit.

The apparatus API now includes `hcn_formation`: initial/added/sampled HCN,
final gas/liquid/collected/wall inventory, net new HCN, feed-carbon fraction,
and scenario amount per declared joule. The Origins result panel displays
formation and fate separately from its historical-yield gate. The postprocessor
runs only after the blind discovery inventory is frozen. The Origins
**Run quick HCN scenario** button loads the declared one-hour opt-in inputs
and immediately solves them; it does not run the five-case control sweep.

## Independent numerical check

The direct EEDF result now folds each independently computed BOLOS and Cantera
distribution through the same imported collision curves. At the tested
five-gas 400 K, 200 Td point, all 77 material inelastic/outcome channels are
within a factor of ten. The exposed summary declares a 1e-20 m3/s rate floor
for this factor-ten count; it does not change the physical calculation or hide
the largest channel discrepancies. This is a **same-data numerical solver
comparison**, not an independent measurement of electron transport or HCN
production. Collision sets remain private and only derived rates/summary are
returned.

## Aqueous HCN/CN- and carbonate acid-base coupling

The apparatus has an opt-in `hcn_nh3_water_equilibrium` mode. It now carries
CN-, NH4+, H+, OH-, HCO3- and CO3-- as explicit **liquid** ledger species, and
solves HCN/CN-, NH3/NH4+, CO2/HCO3-/CO3-- and H+/OH- equilibrium at each
apparatus substep boundary. The mode name is retained for compatibility.
The new carbonate ions remain liquid; neutral CO2 and NH3 join the existing
Henry gas-water subset, which also includes HCN and other neutral solutes.
Each ion update checks atom and charge closure locally; the global
material ledger includes the new ions. Charged molecules are excluded from
the neutral blind candidate enumerator. The mode is off by default for
existing experiments but enabled in the unseeded HCN formation control gate.
This is a physical acid-base property model, not a known HCN product recipe
or yield fed into discovery.

The completed result also includes `aqueous_cyanide_speciation`, a post-run
equilibrium view of each liquid compartment that does not mutate the frozen
inventory. The HCN temperature trend is anchored to reviewed 298.15 K pKa
9.21; ammonium uses NIST pKa 9.245 and thermodynamic derivatives; water uses
IAPWS R11-24 with liquid-water density. The water calculation matches the
release's 300 K, 1 g/cm3 verification pKw of 13.906672. No measured pH is
available. The CO2 neutral pool uses NIST-reviewed 298.15 K Henry solubility
0.035 mol/(kg bar) and a 2400 K temperature coefficient. Neutral NH3 uses a
NIST literature-review value of 61 mol/(kg bar) and 4200 K; dissolved NH3
then enters the separate NH3/NH4+ equilibrium. Carbonate uses NIST
298.15 K pKa 6.351 and 10.329 plus selected ionization enthalpies and heat
capacities. The model still omits other buffers/ions, ionic activities, and
aqueous reaction kinetics. CO2 hydration and dissociation are assumed fast
enough to equilibrate at each boundary; the 373.15 K ammonium, carbonate and
CO2/NH3-solubility adjustments are extrapolations, not direct boiler measurements.
Acid-base equilibrium is updated at substep boundaries, not inside the stiff
gas/transfer solve.
The discharge ledger bounds electrical source terms, but acid-base reaction
enthalpy and the apparatus thermal balance are not included; material-ledger
closure must not be mistaken for complete energy conservation.

Before the carbonate extension, with this mode enabled, the one-hour
CH4/N2/H2/H2O formation gate produced
2.66819e-7 mol net **HCN plus CN-** (2.60593e-7 mol in aqueous compartments),
including 4.068e-11 mol CN- at the end. Spark-off, no-methane and no-nitrogen
controls remain below the numerical reporting floor, and all five control
ledgers close. The 2/4/8/16-substep resolution sweep passes all controls;
the final 8-to-16 changes are 3.05e-5 for net cyanide formation and
2.99e-5 for molecular aqueous HCN. These establish numerical stability
only for the declared one-hour hypothetical scenario. The seven-day endpoint
and full-duration controls reported above also predate carbonate coupling.
A measured pH/transport basis and an experimentally comparable HCN energy yield remain
unverified.

After adding CO2 transfer and carbonate buffering, but *before* adding NH3
transport, the one-hour five-case
gate still passes: the spark-on result is 2.6681906e-7 mol net HCN+CN-,
with 4.0650732e-11 mol final CN-. Its aqueous inventory includes
4.0848281e-8 mol neutral CO2, 2.2594609e-10 mol HCO3- and
2.3177947e-14 mol CO3--; the boiler's **ideal limited-species model pH** is 6.127.
These are simulated values, not measured concentrations. A separate
24-hour *spark-on only* refinement yielded 8.6219972e-6, 8.6201421e-6,
and 8.6192152e-6 mol net HCN+CN- at 2, 4 and 8 substeps, respectively;
the consecutive relative changes were 0.02152% and 0.01075%. Aqueous
HCO3- at 8 substeps was 1.6750393e-9 mol. Main and charged material
ledgers closed and no charged substep failed at all three resolutions.
This is **not** a repeat of the 24-hour five-case controls or seven-day
validation. In this baseline the carbonate amount is small and does not
materially change the predicted HCN total. A CO2-rich test does show the
expected pH and cyanide-ion response, but its magnitude is not calibrated.

With neutral NH3 Henry transport added, the **current** one-hour five-case
gate still passes. The spark-on case makes 2.6682164e-7 mol net HCN+CN-,
with 4.4687351e-11 mol final CN-. The aqueous NH3 and NH4+ inventories are
1.2027184e-7 and 2.9886280e-8 mol; the boiler's ideal limited-species pH
is 6.168. No-spark, no-methane and no-nitrogen controls still form no HCN
above the declared numerical floor; primary and charged material ledgers
close, with no failed charged substeps. The neutral HCN total changed little
from the carbonate-only run, but the ionized fraction did change; the pH is
still not calibrated. A current-model 24-hour *spark-on only* refinement
gave 8.6614592e-6, 8.6596260e-6 and 8.6587102e-6 mol net HCN+CN- at
2, 4 and 8 substeps (relative changes 0.02117% and 0.01058%). At 8
substeps, the final CN- amount was 3.0860351e-10 mol and the model boiler
pH was 6.892; all main and charged material ledgers closed with no failed
charged steps. A current-model 24-hour **five-case control** at two substeps
now passes: the spark-on branch forms 8.6614592e-6 mol net HCN+CN-, the
half-methane branch forms 3.9094423e-6 mol, and spark-off, no-methane and
no-nitrogen controls form none above the reporting floor. Both material
ledgers close and all charged substeps complete in the four spark-on cases.
The first 24-hour attempt did **not** pass: the propagated charged no-methane
control failed 44 of 48 substeps. The charged Pavan ODE lacked the
absent-element invariant already present in the primary GRI ODE; BDF roundoff
then produced a tiny negative CH4+ amount without enough neutral CH4 for a
conservative correction. Enforcing the absent-element invariant removed
those failures. The next rerun exposed 10 late charged no-nitrogen failures:
a negative N2+ amount of about 1.686e-14 mol exceeded its nearly equal
neutral N2 inventory. Trace nitrogen first appeared as 3.0e-40 mol of
N2(v1) after two cycles; the exact ODE derivative was zero, but the BDF
output retained roundoff that grew when carried between cycles. A bounded
post-solve projection now restores exactly absent species to zero, records
the correction, and rejects any correction above 1e-18 mol. A three-hour
nitrogen-free regression and the full 24-hour control sweep now pass. The
post-solve projection removes numerical residue, not a chemical source.
At this checkpoint, seven-day controls and current-model full-duration
refinement remained untested; earlier values cannot stand in for them. This
24-hour input-control gate is **not** an experimental yield validation.

The first current-model seven-day five-case run at two substeps did **not**
pass. The primary spark-on branch produced 4.6236284e-5 mol net HCN+CN-,
and its material ledger and input controls closed. The charged spark-on
branch, however, failed 102 of 336 substeps starting at cycle 118. The first
failure was a negative N2+ amount of 2.828e-11 mol: neutral N2 was abundant
(2.648e-2 mol), but the conservative-ion projection's 2e-11 mol bound was
exceeded. Failed charged steps retained the pre-discharge state, so their
closed global ledger did **not** validate that charged trajectory. The
projection bound was **not** raised. Tightening the charged Pavan BDF solve
from relative/absolute tolerances 1e-9/1e-21 mol and a residence/20 maximum
step to 1e-10/1e-22 mol and residence/40 removed the failure at 120 hours.
The rerun of all five cases through seven days then passed with zero failed
charged substeps, closed primary/charged material ledgers, and no HCN in the
spark-off, no-methane or no-nitrogen controls. At two substeps, the spark-on
primary result remains 4.6236284e-5 mol and the half-methane result is
2.2924350e-5 mol net HCN+CN-. This is numerical input-control evidence at
one resolution, not full-duration convergence, complete energy balance or
experimental yield validation.

A current-model seven-day **spark-on-only** run at four substeps also
completed with zero failed charged steps and closed material ledgers. It
formed 4.6230796e-5 mol net primary HCN+CN- (4.6226576e-5 mol aqueous),
0.01187% below the two-substep primary result. The charged HCN+CN- endpoints
were 3.4508517e-5 mol at two substeps and 3.4503479e-5 mol at four substeps,
a 0.01460% decrease; both had zero failed charged steps. At eight substeps,
the primary endpoint was 4.6228054e-5 mol (4.6223834e-5 mol aqueous),
0.00593% below the four-substep net amount. The charged endpoint was
3.4500962e-5 mol, 0.00729% below the four-substep amount, again with zero
failed charged steps and both material ledgers closed. The decreasing
2-to-4 and 4-to-8 changes support **within-model seven-day spark-on endpoint
stability** at this declared scenario; they are not an uncertainty interval.
The seven-day four-substep spark-off, no-methane and no-nitrogen controls were
then rerun as a **diagnostic subset**: all produced zero HCN, all primary
ledgers closed, and the two charged feed-removal controls had closed ledgers
and zero failed substeps. The runner reports
`partial_formation_diagnostic_only` for that subset; it does not combine
separate runs into a newly passed five-case gate. The four-substep
half-methane diagnostic has now also completed: it forms 2.2922865734e-5 mol
net primary HCN+CN- (2.2920463852e-5 mol aqueous), about 0.0065% below its
two-substep net amount. The charged endpoint is 1.7048427245e-5 mol total
(1.7046757550e-5 mol aqueous). Both material ledgers close, and no charged
substep failed. The separately run seven-day, eight-substep spark-off and
no-methane controls each report zero HCN+CN-; the no-nitrogen primary result
is 5.4346519e-33 mol, far below its 2.6890819e-12 mol reporting floor.
Both eight-substep charged feed-removal controls report zero HCN+CN-, closed
material ledgers, and zero failed substeps. This is a diagnostic subset, not
a newly passed combined five-case gate: the runner correctly labels it
`partial_formation_diagnostic_only`. These results complete the previously
missing seven-day resolution/control repeats for this declared scenario,
but they do not bound model-form error or establish an experimental yield.
An independent, comparable measured cyanide amount and complete thermal
energy balance are still required before a quantitative yield claim.

Provenance: [HCN acidity review](https://doi.org/10.1021/cr0780098),
[NIST ammonium thermodynamics](https://www.nist.gov/system/files/documents/srd/jpcrd615.pdf),
[NIST carbonate thermodynamics](https://www.nist.gov/system/files/documents/srd/jpcrd615.pdf),
[NIST CO2 Henry data](https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=10),
[NIST NH3 Henry data](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417&Mask=10),
[high-temperature ammonium measurements](https://doi.org/10.1016/0021-9614(75)90183-4),
and [IAPWS water ionization release](https://iapws.org/technical-guidance/release/Ionization.download).

## Evidence and use boundary

- [Kamarinopoulou, Wittreich & Vlachos (2024)](https://doi.org/10.1126/sciadv.adl4246)
  report a dry CH4/N2 dielectric-barrier-discharge flow reactor with gas
  products measured by inline chromatography. Their reactor has about 1.57
  cm3 active volume, 50-200 sccm feed flow, and a reported 0.47 s residence
  point at 28.1 W and 22 kHz; the paper reports a best HCN energy cost of
  44 kJ/mmol at some tested conditions. This is a promising **separate dry-flow
  validation branch** because it gives an explicit gas-phase observable and
  power basis. It cannot be compared numerically to our current wet,
  recirculating one-hour inventory: feed fraction, residence, sampled phase,
  and energy delivered to each processed parcel are different. The published
  HCN result is held out of the internal solver, and the previously imported
  mechanism from the same study cannot serve as independent validation. The
  paper's best 44 kJ/mmol HCN energy cost is not labelled as the **1% methane**
  point, so it is not used as an exact-condition pass threshold. Its
  supplementary PDF was not retrievable through the accessible PMC/publisher
  route during this audit; no exact 1% HCN value has been inferred from an
  unread plot.
  The paper's inline gas chromatography is a direct product measurement,
  but its inferred approximately 0.3% streamer-volume fraction and roughly
  1.3 ms kinetic-model residence were obtained by matching the model to the
  experiment. Neither fitted quantity can be imported into our *independent*
  validation scenario as though measured beforehand. The accessible main text
  reports 50-200 sccm feed at ambient inlet conditions, 282 C at the outlet,
  and an approximately 1.57 cm3 active volume; it does not establish a
  uniform reactor temperature or a unique standard-volume-flow temperature.

- [Pearce, He & Hoerst (2022)](https://doi.org/10.1021/acsearthspacechem.2c00138)
  measured mass-spectral changes in reducing H2-dominated cold-plasma
  experiments and reported an HCN interpretation sensitive to methane and
  water. Their gas mixture, flow reactor and
  measurement basis are not the current recirculating spark scenario, so its
  empirical relation is **not** a numerical target for this sweep. More
  importantly, their reported HCN mole fractions are *model-assisted
  estimates* from mass-spectrometer peaks, not independent direct HCN
  concentration measurements: the paper uses its CRAHCN-O kinetic model to
  estimate HCN's contribution to the m/z 26 signal. Its approximately 50%
  inferred reduction after adding 1% water is at most a qualitative trend;
  but the [HCN] = (0.13 +/- 0.01)[CH4] relation is not a clean independent
  order-of-magnitude validation datum for another gas-kinetic model.
- [Ravanbodshirazi et al. (2023)](https://doi.org/10.3390/life13112201)
  showed that spark-generator choice, temperature, pressure, and ammonia
  alter a Miller-Urey-type product mixture. Its reported total extracted
  product weights are **not** HCN yields and are not used as such.
- [Schlesinger & Miller (1983)](https://doi.org/10.1007/BF02101643)
  reported HCN yields up to 10% of methane carbon for some electric-discharge
  feed mixtures. “Up to” is not a universal acceptance threshold and is not
  compared with this one-hour scenario as if residence and measurement basis
  matched.

The previously recorded seven-day 400 K primary HCN diagnostic was about
0.165 nmol per declared electrical joule, versus the Stribling-Miller
conference abstract's approximately 10 nmol per *calorimetric* joule. These
are different energy and collection bases and differ by about a factor 61
numerically; an order-of-magnitude magnitude gate is **not** passed. The
historical comparator remains strict, while the new formation gate is a
separate developmental result.

## Bounded research and build strategy

1. Search a primary experimental paper only when it can supply the next
   measurable trend or comparable quantity. Record feed, discharge, sampled
   phase, units, timing and uncertainty once; stop if the observable cannot
   be mapped to the model without inventing conditions.
2. For each missing reaction or transport parameter, use a primary data
   source or a reviewed mechanism with explicit provenance. Do not tune a
   rate against the HCN validation observation and then call it a blind test.
3. Test one-cycle formation, zero-input controls and ledger closure before
   spending time on a long apparatus run. Run a long trajectory only after
   those checks and numerical convergence pass.
4. Next chemistry priority: state-resolved excited-H/O fate, supported ion
   reactions in the propagated branch, and aqueous reaction kinetics. The
   opt-in HCN/CN-/NH3/carbonate acid-base and gas-water coupling now exists, but
   its boundary-only split, omitted remaining buffers/activities, and thermal enthalpy
   remain unresolved. The seven-day negative controls now pass separate
   four- and eight-substep diagnostic subsets with the charged solver's
   correction bound unchanged; spark-on has 2/4/8-substep refinement.
   Prioritize independently constrained spark exposure and state/aqueous
   kinetics before testing HCN gross formation, destruction and collected
   amount against a truly comparable measurement. Historical-apparatus
   fidelity is a later, optional target.
5. Keep measurement independence explicit: the Pearce 2022 HCN relation is
   model-assisted, while the Kamarinopoulou 2024 HCN assay is direct but its
   streamer exposure is inferred by fitting the output. Do not use either
   derived quantity to calibrate an allegedly blind source model; first
   reconcile independently reported feed, residence, phase and energy bases.

## Theory-first checkpoint (2026-09-26; no kinetic constants imported)

The current weekly Codex usage window was reset before this audit. A one-hour,
two-substep spark-on diagnostic was rerun without changing the solver. It
produced 2.6682164e-7 mol net HCN+CN-, including 2.6059517e-7 mol in
aqueous compartments. Its ideal limited-species boiler pH was 6.16827; the
primary and charged material ledgers closed, and the charged branch had no
failed substeps. This is a partial diagnostic, **not** a repeated five-case
gate or experimental validation.

Before planning another long run, the earlier failure history was checked:
the charged methane-free and nitrogen-free numerical leaks were repaired by
the absent-element invariant and bounded projection; the later seven-day
charged N2+ failure was repaired by tighter BDF tolerances **without** raising
the 2e-11 mol conservative-ion bound. The seven-day two-substep five-case
gate, four/eight-substep spark-on refinements, and four-substep three-control
diagnostic subset already passed. They were **not rerun** in this checkpoint.
The then-missing repeats were the four-substep half-methane case and
eight-substep negative controls, not another replay of the known failures;
both have since completed as diagnostic subsets above. Further identical
numerical repeats are lower priority than independently constraining
discharge exposure and acquiring a comparable aqueous-cyanide benchmark.

### What the physical model must distinguish

1. **Discharge exposure is not reactor residence.** A gas packet's radical
   source is set by the space-time integral of electron density, target
   density and EEDF-folded collision rates over the *active* discharge,
   followed by chemistry in and outside that zone. The present wet scenario's
   50 mL discharge zone in a 5 L gas loop, 200 Td, electron-density input,
   voltage/current/pulse schedule and 0.1 L condenser contact hold-up are
   illustrative, not measured apparatus properties. The dry-flow study of
   [Kamarinopoulou et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10980260/)
   directly measured gas products, but its approximate 0.3% streamer-volume
   fraction and 1.3 ms effective plasma time came from fitting the product
   model. They cannot independently set our exposure parameter.
2. **HCN partition and destruction are separate.** For neutral HCN in dilute
   water, Henry partition constrains equilibrium concentration; HCN/CN-
   acid-base equilibrium sets the volatile neutral fraction. Neither process
   destroys cyanide-family atoms. Hydrolysis and other aqueous reactions do.
   Therefore a measured gas-phase HCN amount, total HCN+CN- formed, and
   collected aqueous HCN+CN- are three different observables and must not be
   interchanged.
3. **Aqueous kinetics may matter on the seven-day scale.** The peer-reviewed
   [Miyakawa, Cleaves & Miller (2002)](https://doi.org/10.1023/a:1016514305984)
   measured HCN/formamide hydrolysis over 30–150 C and pH 0–14. A newer
   [Murali, White & Rimmer (2026) preprint](https://arxiv.org/html/2608.00025v1)
   reports an acid/base-catalyzed rate law for formate production,
   `r = (k_acid[H+] + k_base[OH-])[HCN]`, with temperature-dependent rates.
   Its full kinetic analysis includes the neutral HCN fraction. Using **only
   its tabulated pure-water base term** as a screening calculation at 373.15 K,
   the current model water-ionization relation gives pKw=12.2534,
   HCN pKa=8.0121, and k_base=4.76 M^-1 s^-1. At *fixed assumed* pH 6.17 or
   6.89, including the neutral-HCN fraction gives pseudo-first-order
   base-only `k*t` of about 2.34 or 11.6 over seven days.
   This is **not** a predicted apparatus loss: HCN is produced continuously,
   moves among compartments, and sees changing pH and temperature. It shows
   that silently treating all aqueous HCN as inert could be consequential.
4. **Do not yet import that rate into the production solver.** The 2026
   source is a preprint, and its abstract lists a pure-water acid activation
   energy of 67.1 kJ/mol whereas Table 1 lists 60.0 kJ/mol. The discrepancy
   needs resolution against the PDF/SI or authors' correction. Its experiments
   are buffered and relatively concentrated compared with this hypothetical
   apparatus, whose pH and liquid residence distribution are unmeasured.
   An atom/charge-balanced hydrolysis network also needs explicit formamide,
   formate and ammonium fate, rather than deleting HCN from inventory.
5. **Henry-law temperature sensitivity is unresolved outside measured
   ranges.** The current property subset uses NIST's literature-review
   298 K HCN solubility of 12 mol/(kg bar) and inverse-temperature coefficient
   5000 K. [Ma et al. (2010)](https://doi.org/10.1021/es1001192) directly
   measured dilute HCN from 287–311 K and found 9.86 M/atm at 298 K and a
   fitted 8205.7 K inverse-temperature slope. These are useful independent
   partition checks near ambient conditions, but extrapolating that fit to
   the 373 K boiler (or even the 280 K condenser) would not be measured
   evidence. Do not silently replace the existing coefficient with it.

### Decision before implementation

The next small gate is an **aqueous-fate and residence sensitivity audit**,
not an asserted new HCN yield: record each liquid compartment's HCN/CN-
inventory, temperature, pH, and time exposure; compare accepted hydrolysis
rate sources over their measured domains; then decide whether a constrained,
atom/charge-balanced kinetic extension is justified. Independently, keep
discharge active-volume and plasma exposure as explicit uncertainties until
an experiment reports them without product fitting. Any resulting HCN amount
must be compared against an independently measured observable with matching
feed, phase, time, flow and energy basis. Do not calibrate a missing rate or
exposure using the HCN output and then call the same output validation.

### Read-only aqueous exposure diagnostic (2026-09-26)

The apparatus now reports endpoint-integrated HCN and total cyanide-family
mol-seconds for each aqueous compartment, with each step's temperature and
ideal-model pH, for the primary and charged branches. This is a diagnostic
derived after the predictive trajectories are frozen: it does not apply a
hydrolysis loss, change the reaction network, measure individual molecular
residence time, or establish a yield uncertainty. Synthetic two-step
quadrature and invalid-timeline tests pass, as does the existing HCN test
module (14 tests total). A one-hour spark-on integration also completes with
the new diagnostic attached; its formation amount is unchanged at
2.6682164e-7 mol HCN+CN- generated. The next decision is whether a
source-supported, atom/charge-balanced aqueous kinetic model can be justified
over this scenario's temperature and pH range.

### Modern wet-spark benchmark candidate: Barth et al. (2024)

[Barth et al.](https://doi.org/10.1051/0004-6361/202347286) report a
closed, stirred 1 L spark flask with 50 mL of water, a 49 kV-rated Tesla-coil
generator, and N2/CO2/H2 feed mixtures. They measured aqueous cyanide by a
colorimetric assay and inferred gas HCN from that measurement using Henry's
law; the latter is **not** an independent gas HCN assay. Their
[CC BY experimental dataset](https://research-portal.st-andrews.ac.uk/en/datasets/spark-experiments-and-photochemical-simulations-to-identify-light/)
includes individual feed compositions, spark durations, and the reported
aqueous `HCN` concentration in micromolar. Its README distinguishes 2 h wet,
2 h dry, and overnight wet experiments. For example, overnight wet record
`220725001` reports H2 fraction 0.00943, CO2 fraction 0.0571, duration
925 min, and aqueous assay result 777 uM. At the paper's 50 mL water volume,
that concentration corresponds arithmetically to 38.85 umol of assayed
cyanide-family material. This is a **held-out measurement candidate**, not a
source term, rate, fit target, or current model prediction.
The read-only `barth_2024_aqueous_cyanide_v1.json` benchmark card records all
13 nonblank-HCN rows (not just the example); every field was checked against
the source CSV, and the HCN module's 15 tests pass. Blank source cells were
not converted to zero or nondetection claims.

This is substantially closer in sampled phase than the dry-flow 2024 DBD
paper, but it is **not yet like-for-like** with our heated CH4/N2/H2/H2O
recirculating apparatus: their initial carbon source is CO2 rather than CH4;
their flask is closed and room-temperature rather than a boiler/condenser
loop; and the published energy estimate uses `E = 1/2 U I t`, not our declared
coupled electrical energy or measured calorimetry. No active-spark volume,
electron density, or time-resolved current waveform is supplied in the paper.
Their paper itself notes that repeated sparking of the same gas over water
can change HCN retention compared with an atmospheric single-pass parcel.
Therefore neither the 777 uM record nor its energy-normalized value may be
used to declare the existing model accurate. The next implementation decision
is to determine whether a closed, stirred, isothermal gas/water contact mode
can be represented without fitting its spark exposure or mass-transfer rate
to aqueous HCN. The older Stribling-Miller calorimetric aqueous benchmark
already exists separately in `backend/data/benchmarks`; it was not rediscovered
as a new validation result here.
At the time of that benchmark review, the EEDF collision manifest covered H2,
H2O, CH4, NH3, and N2 but **not CO2 or O2**. Therefore merely changing the feed in the current
apparatus profile would leave the Barth discharge physically unsupported.
Before a Barth comparison, the remaining requirements are: reviewed complete
CO2/O2 collision coverage and compatible fragment outcomes; an independently
specified spark exposure; and a closed-flask gas/water mode with the measured
water volume and room-temperature phase condition. Unknown HCN assays in
other CSV rows must remain missing, not become negative controls.

The first apparatus part is now implemented as an opt-in
`closed_stirred_flask` operation mode. It requires zero boil/condense stages,
zero gas circulation and a common gas/water temperature; it keeps all water
in one stirred liquid compartment and applies Henry partition before and
after each spark integration substep. That is an **instantaneous-contact
limit**, not a measured interfacial transfer coefficient or residence-time
distribution. A no-spark CO2 absorption test and a CH4/N2/H2/H2O spark test
pass material closure and leave both condenser compartments empty. The HCN
and prebiotic test modules pass together after this change. The spark test
uses only supported gas species and is deliberately *not* compared with the
Barth CO2/O2 measurements. The closed mode does not replenish consumed water
vapor to saturation, resolve spark geometry, or close thermal energy.

The closed-flask kinetics check exposed two zero-flow/occupancy errors:
declared gas-phase reactions were multiplied by the circulation pass, which
is zero by construction in a closed vessel; and the coupled heavy-gas solver
used only spark-zone mean occupancy as its residence time. The solver now
integrates ordinary mass-action gas chemistry over the whole stirred flask, while declared
electron-impact channels use the discharge-zone volume fraction times the
electrical pulse duty. Coupled heavy-gas chemistry runs over the full closed-
flask interval, while electron source opportunities remain localized to the
declared spark-zone fraction. An active closed-flask spark must state its
discharge-zone volume explicitly; a hypothetical value remains flagged as hypothetical
and cannot validate yield. Analytical first-order fixture tests confirm both
scalings, and a supported-feed spark test confirms full-interval residence,
material closure, and rejection of an omitted active-zone volume.
The HCN module passes 21 tests after the change. This repairs a silent
nonreaction for *declared* channels; it does not independently determine
actual streamer exposure or the true active volume.

The [IST-Lisbon LXCat database](https://nl.lxcat.net/data/) lists a reviewed
O2 complete set (Phelps-derived and updated) and a CO2 complete set updated
with channel-resolved dissociation work in 2025. Its descriptions say the
CO2 set needs the CO2(010) superelastic channel for agreement at low E/N,
and the O2 set needs rotational transitions below 10 Td. Its O2 page also
lists a separate dissociative-ionization channel explicitly marked **not**
part of the complete set. These distinctions must survive any import and
source-term mapping; selecting every visible process would not be a justified
complete-set import. The site's terms also prohibit third-party
redistribution, so raw curves must remain private and uncommitted. The CO2/O2
curves are now imported privately as EEDF inputs, but the Barth
benchmark remains held out rather than a claimed validation.
An active `closed_stirred_flask` Boltzmann run now fails early if any initial
neutral gas lacks coverage in its selected collision set. This prevents the
CO2/O2 Barth feed from being silently renormalized to the supported gases
while the EEDF still reports a result. The regression first demonstrates
the present CO2 refusal, while a supported CH4/N2/H2/H2O feed still runs.
This is an input-coverage safeguard, **not** CO2/O2 plasma chemistry; newly
formed gas species can still be outside the collision set and are separately
reported through the covered-gas fraction.

### Private CO2/O2 collision-set intake (2026-09-26)

After the user accepted LXCat's data terms, the reviewed IST-Lisbon
**complete** CO2 (17 channels) and O2 (14 channels) sets were downloaded to
ignored local storage. The separate O2 dissociative-ionization channel
explicitly marked *not part of the complete set* was not selected. The raw
text is hash-pinned in the collision manifest and is **not** included in
project source or test fixtures. Two O2 ionic product labels are normalized
in memory without changing numerical tables. The CO2 dissociative-attachment
source explicitly states a 3.9 eV onset, which the BOLOS BOLSIG attachment
parser otherwise reports as zero; the onset is restored in memory for
collision-energy bookkeeping. The imported set is opt-in and limited to
10–200 Td, the overlap of the selected H2 and O2 contributor domains without
additional rotational channels. BOLOS does not currently consume the
CO2(010) statistical-weight flag needed for low-field reverse-superelastic
behavior, so this is not an all-field validated plasma mechanism.

At a **synthetic numerical check point** (H2 0.01, CO2 0.05, O2 0.001,
H2O 0.02, N2 0.919; 298.15 K, 200 Td), the 400-cell BOLOS EEDF had mean
energy 5.3339 eV and the 800-cell result 5.3212 eV (0.24% difference).
The CO2 -> CO(X) + O(1D) collision rate differed by about 0.28% across
those two grids. A same-curve Cantera calculation at the nominal 400-cell
resolution gave mean energy 5.0417 eV, 5.8% below BOLOS; all 84 material
matched channel rates were within a factor of ten, but this is **not** an
independent swarm benchmark. BOLOS still reports a floating-point iteration
warning and clips one reconstructed H2O elastic curve, so the solver marks
this point sensitive. These checks establish that the curves can be parsed
and a finite, moderately grid-stable EEDF can be solved; they do **not**
validate spark-source amounts or aqueous HCN yield.

State-label parsing now distinguishes real products from `+` inside an
excited-state name. The O2 set declares O(3P)+O(3P) as a ground-state neutral
fragment outcome. Its O(1D) channel is excited, and both identified CO2
neutral fragment outcomes contain O(1D) or CO(a3P). The current chemistry
does **not** silently treat those as ground-state O/CO sources; attachment
and ionization products likewise require charge-balanced downstream kinetics.
The next gate is a source-supported state/ion reaction and relaxation network,
plus independently constrained discharge exposure. Only after those can a
closed-flask run against the held-out Barth aqueous cyanide assay be meaningful.
The synthetic one-minute closed-flask CO2/O2/H2/H2O/N2 spark regression now
selects the new collision set, solves the EEDF, and closes the material ledger;
its discharge-zone volume and power remain test-only assumptions. The six
relevant backend test modules pass together (119 tests), the frontend builds,
and both the collision manifest and held-out benchmark card parse as valid
JSON. No experimental HCN amount or energy-normalized yield is claimed.

The next theory intake is **state fate**, not a fitted CO2-to-HCN shortcut.
The closed-flask test explicitly verifies that the CO2 -> CO(X) + O(1D)
collision is reported as excluded from ground-state neutral injection because
O(1D) is excited. Published [CO2–N2](https://arxiv.org/abs/2301.08938)
and [CO2–O2](https://arxiv.org/abs/2210.07778) kinetic studies provide
state-resolved, externally tested mechanisms, including vibrational transfer
and CO2 back reactions. Their validations used 0.4–5 Torr DC glow discharges
(the CO2–N2 subset 0.6–4 Torr), not an approximately atmospheric, pulsed
spark over water. They are candidates for checking species/reaction topology
and source provenance, **not** interchangeable validation of this apparatus.
Before adding O(1D), CO(a3Π), N2 excited, and ion relaxation reactions, each
adopted rate must have a stated temperature/pressure domain and branching
outcome, conserve atoms and charge, and bookkeep chemical/excitation energy.

### O(1D) first-collision theory audit (2026-09-27)

The first state-fate calculation is deliberately **read-only**. It evaluates
five competing O(1D) removal hazards, `k_i [M_i]`, from the [NASA/JPL
Evaluation 17 rate table](https://jpldataeval.jpl.nasa.gov/pdf/JPL%2010-6%20Final%2015June2011_20110620.pdf)
at 298 K: physical deactivation by N2, O2, and CO2 to O(3P), versus OH/H
formation with H2O and H2. The rates are held constant only within a
295–301 K diagnostic window, not extrapolated to the boiler. The [direct
H2O/H2 measurements](https://kinetics.nist.gov/kinetics/Detail?id=2010VRA%2FPEE9213-9221%3A0)
and [CO2 measurement](https://kinetics.nist.gov/kinetics/Detail?id=2018NUN%2FHIC4002-4008%3A0)
give useful cross-checks near room temperature; they are not HCN-yield
calibrations. Source products and competing rates were specified independently
of the Barth aqueous cyanide measurements. These tabulated rates describe
near-thermal O(1D); electron dissociation may produce translationally hot
fragments, so the fractions below are conditional on thermalization before
subsequent reaction. The small O(3P) yields measured from O(1D)+H2/H2O are
not resolved as separate branches in this five-channel screening audit.

For the synthetic five-gas test mixture `N2 0.919 / CO2 0.05 / O2 0.001 /
H2O 0.02 / H2 0.01` at 298.15 K and 1 mol in 25 L, the resulting
conditional fractions are 0.86744 physical O(3P) and 0.13256 OH/H-forming;
the mean first-collision lifetime is 1.058e-9 s. This is a strong timescale
separation from a millisecond discharge pulse, **not** a product amount.
The one-minute closed-flask regression contains a covered first interval,
but subsequently develops CO and trace H, H2O2, HNCO and other bath species.
Their O(1D) rates/branches were not in this initial five-partner audit, so the later
interval is explicitly *unavailable* rather than renormalizing the known
hazards to 100%. Those unreviewed species are below 0.1% of gas molecules in
this regression, but a small mole fraction alone does not bound their O(1D)
reaction hazard without rate coefficients. The UI shows the first conditional result and the later
coverage failure. No O(1D) products have been injected into the primary
inventory, and CO(a3Π) is still unresolved. The next kinetic gate is to
validate fate/branching for the newly appearing bath species, then couple
state-resolved reactions and excitation/thermal energy to the gas solver.
The seven relevant backend modules pass together (126 tests), and the
frontend production build succeeds after adding this diagnostic.

### Benchmark eligibility recheck (2026-09-27)

The user's LXCat approval does not require another download: the complete
IST-Lisbon CO2/O2 export is already held privately and hash-checked above.
It supplies electron collision probabilities, **not** the number, volume,
current waveform, or energy deposition of actual sparks. It also does not
provide the subsequent excited-state or aqueous reaction network.

A targeted check of [Jiang et al. (2024)](https://doi.org/10.1073/pnas.2400819121)
found valuable independent interfacial-discharge observations: gas-to-water
filament geometry and measured current, CN radical optical emission, and
quantified aqueous nitrate, nitrite, ammonium, and formate. However, its
reported stable aqueous-product set does not include an HCN+CN- assay; it
used continuously flowing N2/CO2 over an interfacial discharge, not the
closed, stirred Barth flask. It can test discharge or ancillary chemistry
assumptions only after matching those specific conditions. The optical CN
signal is **not** a measured amount of collected cyanide and must not be
used as the HCN validation target.

[Barth's thesis](https://research-repository.st-andrews.ac.uk/bitstream/handle/10023/28875/Thesis-Patrick-Barth-complete-version.pdf?isAllowed=y&sequence=6)
describes the BD-50E as producing a continuous spark series at 0.5 MHz,
while the [HCN paper](https://doi.org/10.1051/0004-6361/202347286) uses
the nominal `E = 1/2 U I t` with 49 kV and 1 mA. Neither item supplies the
needed measured voltage/current waveform, filament volume, or local electron
density for the *same* wet reactor. Treating 0.5 MHz as a calibrated count
of identical energy-deposition events would add an unsupported assumption.
Thus the held-out Barth aqueous cyanide records remain the closest phase-
matched target, but they cannot yet verify an absolute predicted amount.
The next practical comparison needs independent spark-exposure measurements
for a cyanide-assayed reactor, not additional LXCat cross-section files alone.

### CO contribution to the O(1D) fate audit (2026-09-27)

The initial audit failed closed as soon as CO appeared. Independent
[Davidson et al. measurements](https://ntrs.nasa.gov/citations/19780065438)
now supply its *total O(1D) deactivation* coefficient,
`(4.7 +/- 0.9)e-11 exp[(126 +/- 33) cal mol-1 / RT] cm3 molecule-1 s-1`
over 113–333 K. The implementation evaluates the central expression in the
audit's existing 295–301 K window; it does not extend the other 298 K rates
to the full Davidson range. CO therefore contributes to a conditional mean
first-collision lifetime. The reported lifetime interval varies only the
published CO coefficient and exponent to their stated extremes; it is not
a statistical confidence interval or the full uncertainty of the five other
rates and nascent-fragment thermalization. Crucially, this measurement alone does not specify
the split between O(3P)+CO and stabilized CO2 in this wet, near-atmospheric
reactor. When CO is present and all other bath partners are covered, the
diagnostic reports a *partially resolved* fate with a separate unresolved-CO
fraction, not a fabricated ground-O or CO2 product. If H, H2O2, HNCO, or
other unreviewed partners are also present, it still returns unavailable
without renormalizing the known hazards. No CO/O(1D) products enter the
primary inventory, and this does not yet constrain the absolute HCN yield.

### Aqueous HCN hydrolysis rate screen (2026-09-27)

The apparatus now computes a **read-only, time-dependent screening estimate**
of HCN/CN- hydrolysis loss from each simulated aqueous compartment's
right-endpoint inventory, temperature, and *ideal limited-species* pH. It
does not subtract that loss from the simulated amount, add formamide/formate,
or close the thermal ledger. Any cyanide-bearing snapshot without a supported
pH or outside the deliberately narrow 293.15–323.15 K screen window is marked
unavailable; the remaining partial sum must not be interpreted as a complete
loss. The response counts every unavailable snapshot but stores at most 20
details to keep long experiments' diagnostics bounded. Zero-cyanide controls
give exactly zero loss.
The same screen is attached separately to the primary and optional propagated
charged branch; it does not blend their inventories or imply a model-ensemble
yield.

The independent [Murali, White, and Rimmer 2026 experimental
preprint](https://arxiv.org/pdf/2608.00025), with [archived primary
data](https://doi.org/10.7910/DVN/WT6NKZ), gives the pure-water expression
`d[formate]/dt = (k_acid[H+] + k_base[OH-])[HCN]`. The screen evaluates its
Table 1 reference log rates at 273.15 K and Arrhenius activation energies,
then multiplies by the molecular HCN fraction to express a hazard on total
HCN+CN-. It reports the Table 1 log-rate spreads as a *rate-parameter-only*
range, not a confidence interval for this reactor. The paper's abstract
quotes acid activation energy 67.1 kJ/mol, while Table 1 quotes 60.0
kJ/mol; the detailed table is the nominal choice and 67.1 is a separately
reported sensitivity. The discrepancy is unresolved. The study monitored
formate as a proxy for HCN hydrolysis; the screen does not assume a resolved
intermediate yield in our apparatus.

For a **conditional fixed-pH example** at 298.15 K, pH 4.53, and 925 min,
the pure-water central rate is about `8.26e-10 s-1`, corresponding to about
`0.00459%` loss of an initially present cyanide pool. At fixed pH 12 it is
about `0.306%`. These pH values are not measurements of the Barth reactor;
the 4.53 value is an idealized synthetic-mixture estimate. Actual pH, salts,
transport, within-step arrival, other aqueous sinks, and source-model errors
are outside this range. The screen uses endpoint quadrature, not a coupled
kinetic integration, so its amount is neither an actual yield correction nor
a strict bound. It suggests pure-water hydrolysis is unlikely to be the
leading *15-hour, near-room-temperature* error under those assumptions, but
does not settle heated or long-duration cases. The priority remains an
independently constrained spark exposure and comparable measured collected
HCN+CN- target.

The first broader HCN test attempt exhausted Windows committed virtual
memory (`bad allocation` and Python `MemoryError`), despite ample free
physical RAM. This was an environment failure, not a chemistry result.
After commit headroom recovered, the previously failing integration test
and combined suite passed, as recorded in the next section. No kinetic
parameter was changed to fit a benchmark.

### Measured O(1D)+HCN fate and retest (2026-09-27)

HCN itself becomes a bath partner after it forms, and the earlier O(1D)
first-collision audit rejected that later mixture entirely. The
[Strekowski, Nicovich, and Wine measurement](https://doi.org/10.1002/cphc.201000550)
now supplies both the **total** O(1D)+HCN removal coefficient
`1.08e-10 exp(105/T) cm3 molecule-1 s-1` and the O(3P) quenching yield
`0.150 exp(200/T)` over 193–430 K. At 298.15 K the latter is about 0.293
per O(1D)+HCN collision. These expressions enter only the existing
295–301 K read-only state-fate audit. Its measured ground-O share is
separated from the remaining HCN reactive share, whose molecular products
are **not** assigned. The reported H-atom yield at 298 K includes a fast
secondary contribution and is not treated as a primary H-producing branch.
No HCN loss, ground O, or other product is injected into the predictive gas
inventory; excitation energy and heat from this channel remain unclosed.

The HCN and CO unresolved shares remain separate and the first-collision
fractions sum to one. An unsupported bath partner such as H2O2 still makes
the audit unavailable, not falsely complete. The prior Windows commit-space
failure was transient: the previously failing spark integration test now
passes, and the combined HCN formation, aqueous hydrolysis, and O(1D) suite
passes (39 tests). The separate opt-in CO2/O2 LXCat collision-set module
passes (7 tests), and the frontend production build passes. This restores
regression confidence for the new diagnostic, but it does **not** resolve
the absolute spark exposure, H2O2/H/HNCO fate, or comparable collected
cyanide validation.

### O(1D) remaining-bath literature triage (2026-09-27)

The first-collision audit still fails closed when materially present H,
H2O2, or HNCO lacks a supported channel. A bounded check found a primary
[O(1D)+HNCO experiment](https://doi.org/10.1021/j100166a017) reporting a
reaction coefficient of `(3.3 +/- 0.3)e-10 cm3 molecule-1 s-1` and a
`0.44 +/- 0.15` primary branch to electronically excited NH(a1Delta).
The accessible [paper abstract](https://hero.epa.gov/reference/7545706/)
also distinguishes primary from secondary NH/NCO production; those products
cannot be collapsed into a single ground-state reaction. The accessible
record did not establish its temperature/pressure dependence or a complete
primary branching sum for this audit's narrow 295-301 K domain. No HNCO
rate, product, or heat term was therefore inserted into the predictive
inventory or the diagnostic. A search for an independently measured
O(1D)+H2O2 total rate and product branching did not resolve that channel
either; an assumed 50/50 split in a secondary mechanism is not a measured
branching ratio. The next intake needs the original experimental conditions
and separate defensible H/H2O2 fate data, after which a read-only hazard
can be tested before any state-resolved propagation.

### Spark-exposure identifiability and local sensitivity (2026-09-27)

The [BD-50E manufacturer's manual](https://www.electrotechnicproducts.com/content/BD-50E%20High%20Frequency%20Generator.pdf)
specifies a nominal 20-49 kV peak setting (about +/-4 kV) and approximately
1 mA maximum output, but says actual voltage changes with electrode loading,
unit-to-unit variation and contact condition. Its 0.5 MHz operating frequency
does not give the pulse waveform, plasma-channel volume, or deposited energy
in Barth's gas/water flask. Thus the paper's nominal `E = 1/2 U I t` remains
an electrical estimate, not an independently measured spark exposure or
calorimetric datum for our model.

The new `backend/scripts/audit_hcn_exposure_sensitivity.py` varies one
**hypothetical** parameter at a time around the existing one-hour,
two-substep CH4/N2/H2/H2O scenario. It loads no Barth HCN measurement and
never changes the initial feed. These half/double factors are mathematical
probes, **not** measured parameter intervals or error bars. "Aqueous" here
is the sum across modeled liquid compartments, not solely material withdrawn
from a collection trap. Run from `backend` with
`..\.venv\Scripts\python.exe -m scripts.audit_hcn_exposure_sensitivity`:

| Input varied | Half input: aqueous HCN+CN- (mol) | Baseline (mol) | Double input (mol) | Central log elasticity |
| --- | ---: | ---: | ---: | ---: |
| Active zone volume (0.025/0.05/0.1 L) | 3.2309283e-7 | 2.6059517e-7 | 2.1390412e-7 | -0.2975 |
| Declared electron density (0.5/1/2 x 1e15 m-3) | 2.6059517e-7 | 2.6059517e-7 | 2.6059517e-7 | 0.0000 |
| Electrical coupling fraction (0.15/0.3/0.6) | 1.8913951e-7 | 2.6059517e-7 | 3.8695311e-7 | +0.5164 |

The density result is **not** evidence that electron density is physically
irrelevant: the model's electron-drift energy cap lowers all three declared
values to about 1.66e12 m-3 at the baseline 2160 J electrical estimate.
The other two axes therefore change the effective source exposure while
the declared-density axis is inactive in this range. Every tested primary
and charged material ledger closes and every charged case has zero failed
substeps; the complete thermal balance remains open. The corresponding
charged-branch total HCN+CN- has similar local elasticities (-0.2873 for
zone volume and +0.5723 for coupling). These one-hour responses cannot be
extrapolated to the seven-day endpoint or treated as an experimental HCN
uncertainty interval.

The [Barth dataset README](https://pure.st-andrews.ac.uk/ws/portalfiles/portal/300064695/Readme.txt)
shows that NO, nitrate, nitrite, ammonium, CH4 and N2O were also measured in
some cyanide-assayed runs. A future staged calibration could use a **separate
non-HCN observable** while keeping HCN held out, but that would constrain a
combination of discharge source and chemical mechanism, not measure the
spark waveform or filament volume directly. The present state/ion and
aqueous networks are not complete enough to fit those co-products without
confounding source strength with missing pathways. No such calibration was
performed. The next independent physical input remains a measured waveform,
active-channel geometry and/or electron-density history for a cyanide-assayed
wet discharge, plus a closed energy ledger.

### CO2-fed wet-flask pathway probe and charged-source correction (2026-09-27)

The new `backend/data/benchmarks/barth_2024_feed_only_220725001_v1.json`
records only the published initial gas fractions for one cyanide-assayed
Barth run (`220725001`): H2 0.00943, CO2 0.0571, O2 0.000356, and H2O
0.016, with N2 as the balance. It records the reported 1 L flask and
50 mL water and source/provenance, but contains **no measured HCN value**.
The separate `backend/scripts/probe_barth_feed_pathway.py` initializes a
0.95 L headspace at an assumed 1 bar and 298.15 K. Its 50 mL active
discharge zone, 0.2 microsecond pulse width, 30% electrical coupling, and
200 Td field are deliberately hypothetical, not measurements in that run.
Only one hour of the reported 925-minute spark duration is probed. This
is a pathway/control test, **not** an independent yield comparison.

That probe exposed a real source-mapping error: LXCat's parenthesized
`CO2(+,X)` and `O2(+,X)` labels had been mapped to neutral CO2/O2 while
the electron source was still added. The resulting charged branch failed
its material/charge integration. The Pavan mapper now retains the cation
and checks each mapped product's atoms and charge against the collision
identity before accepting the source. The O2 ground-cation state is no
longer mislabeled as electronic excitation. Both primary and charged
material ledgers now close in this diagnostic, with zero failed charged
substeps. A 60-second spark-off control gives exactly zero HCN.

At two integration substeps, the **60-second spark-on diagnostic** gives
about `1.08e-7 mol` gas CO but less than the `1e-12 mol` numerical reporting
floor for aqueous HCN+CN-. At **one hour**, it gives about `1.40e-6 mol`
gas CO and `6.33e-12 mol` aqueous HCN+CN- in the primary branch; the
separate charged branch gives `6.26e-12 mol` aqueous HCN+CN-. These
branches are alternatives, not amounts to add. The hour result exceeds
the numerical floor but is a small, unvalidated model signal; it must not
be compared to the held-out 925-minute assay or extrapolated to it.
The 60-second O2 and N2 source extents are nonzero, but **CO2 neutral
fragment injection is zero**. Two `CO2 -> CO(X) + O(1D)` collision
channels are recognized yet excluded because the O(1D) product cannot
currently be propagated with a defensible fate and energy allocation.
Their raw one-hour collision opportunities (each about `7.55e-6 mol`)
are *not* accepted source amounts or predicted additional CO/HCN yield.
The solver also emitted SciPy finite-difference Jacobian overflow warnings
in the charged runs even though they completed and their material ledgers
closed. A one-hour substep check shows that **aqueous HCN is not converged**:

| Substeps in one hour | Primary aqueous HCN+CN- (mol) | Charged-branch aqueous HCN+CN- (mol) | Primary gas CO (mol) |
| ---: | ---: | ---: | ---: |
| 2 | 6.33e-12 | 6.26e-12 | 1.396e-6 |
| 4 | 9.43e-12 | 1.08e-11 | 1.397e-6 |
| 8 | 1.70e-11 | 1.77e-11 | 1.394e-6 |
| 16 | 2.82e-11 | 2.85e-11 | 1.387e-6 |

CO is comparatively stable, but the trace HCN keeps increasing as the
step shrinks; closed ledgers alone do not validate this small numerical
signal. At 32 substeps, the apparatus refuses the result because the
carbon ledger has a `+1.1025e-12 mol` residual against its `1e-12 mol`
tolerance. This is a small fraction of the initial carbon inventory but
large enough relative to the trace HCN amount to matter. The failure now
reports the element and residual rather than an opaque error; **no 32-step
yield is reported**. The warning's impact, the carbon drift and the
substep dependence require a solver/source-splitting and phase-transfer
audit before quantitative use. The complete thermal energy
balance remains open. The 120-test focused HCN/O(1D)/aqueous/plasma/
CO2-O2/prebiotic suite passes after the cation fix, as do the two
post-audit targeted regressions. All 27 HCN-formation tests pass after
the later diagnostic-error improvement; the SciPy warnings remain.

Next: diagnose and repair the trace-HCN substep dependence and charged
Jacobian warning, resolve state-specific CO2/O(1D) production and fate
with reviewed rates and explicit energy disposition, then independently
constrain the wet discharge exposure. Preserve the
HCN assay as a held-out target until source and transport assumptions
are independently supportable; do not use its known amount to tune the
internal mechanism.

### Conservative trace-state repair and one-hour numerical screen (2026-09-27)

The earlier 32-step carbon failure was localized to **GRI gas integration**,
not an acid-base atom-balance error. One substep's raw BDF state conserved
gas elements to `4.13e-16 mol`, but contained `-9.70e-13 mol` CH2(S) plus
smaller negative CH and HCNN traces. Simply replacing those values by zero
created about `1.02e-12 mol` carbon and `1.99e-12 mol` hydrogen. Disabling
aqueous acid-base updates had removed the failure in an isolation control
because it changed the upstream gas state; the acid-base update's own
maximum residual in the failing run was only `1.57e-15 mol`.

The GRI branch now projects small negative states onto nonnegative amounts
while preserving the **raw solver's** element totals within a checked
numerical tolerance. It adjusts only already-positive species, reports
every clipped species and the total absolute correction, and rejects a
negative species exceeding `1e-12 mol` or a projection that cannot close.
This is numerical error correction, not a new reaction or a chemistry
parameter fitted to HCN. In the separate Pavan branch, when a negative
cation is too small for a same-neutral correction, at most `1e-17 mol`
per such ion is treated as *explicitly reported* unprojected roundoff;
atom/charge residual checks still apply. The feed-only probe additionally
rejects more than `1e-15 mol` cumulative unprojected ion roundoff. No
unreported ion source is added.

`backend/scripts/audit_barth_feed_convergence.py` now reruns the one-hour
feed-only scenario at 48, 96, and 192 integration steps. A one-cycle/48-step
run and a two-cycle/24-step run agreed exactly, confirming that splitting
the declared hour across equivalent cycles does not change this case.
The end-to-end audit produced:

| One-hour steps | Primary aqueous HCN+CN- (mol) | Charged aqueous HCN+CN- (mol) |
| ---: | ---: | ---: |
| 48 | 4.631319e-11 | 4.664046e-11 |
| 96 | 5.064515e-11 | 5.180348e-11 |
| 192 | 4.889739e-11 | 5.087909e-11 |

The 96-to-192 relative changes are `3.57%` and `1.82%`, respectively,
meeting the audit's **5% one-hour substep screen**. Both material ledgers
close and the charged branch has zero failed substeps at all three tested
resolutions. Maximum summed primary projection was `2.79e-12 mol` across
the 96-step run; maximum cumulative unprojected ion roundoff was
`1.72e-18 mol`. The screen's pass status does **not** establish physical
accuracy, full-duration convergence, an HCN yield comparison, or a closed
thermal ledger. It uses only the feed-only card and does not load the
held-out assay. The Pavan finite-difference Jacobian overflow warning
still appears in some runs and needs a separate stiffness audit. The
123-test focused HCN/O(1D)/aqueous/plasma/CO2-O2/prebiotic regression
suite passes after these numerical changes.

Next scientific gate: independently constrain the wet discharge source
(waveform, active volume/energy coupling), resolve CO2-to-CO(X)+O(1D)
and subsequent state fate/heat before injecting that channel, and test a
time-dependent gas-water/aqueous calculation at the measured duration and
sampled phase. Compare HCN to the held-out assay only after those inputs
and the apparatus energy basis are defensible.

### Measured minor O(1D)+H2O/H2 quenching and bath audit (2026-09-27)

The earlier read-only O(1D) audit assigned **every** O(1D)+H2O and
O(1D)+H2 collision to its dominant radical products. The primary
[Vranckx, Peeters, and Carl experiment](https://doi.org/10.1039/b923959f)
instead measured total removal rates of `1.70e-10 exp(36 K/T)` and
`(1.35 +/- 0.05)e-10 cm3 molecule-1 s-1`, respectively, with minor
O(3P) quenching fractions `0.010 +/- 0.003` and `0.007 +/- 0.007`
(the latter at room temperature). The diagnostic now allocates these
measured nominal fractions to ground-O quenching and the remainders to
the atom-balanced `2 OH` or `OH + H` paths. It records the branch source
and reported spreads but does **not** inject either product into the
predictive inventory. Its domain stays at 295-301 K because the other
included bath rates are still fixed near 298 K, despite the wider domain
of the Vranckx study.

For the initial Barth feed-only diagnostic at 298.15 K, the conditional
*thermalized* O(1D) first-collision audit gives `0.8899` ground-O fraction,
`0.1101` OH/H radical fraction, and a roughly `1.05 ns` mean removal
lifetime. This is not a measured nascent-fragment thermalization time.
In a separate 60-second, two-substep run the first substep is covered,
but by the next substep the audit correctly becomes unavailable. The
unreviewed bath totals about `4.56e-6` of gas molecules; its largest
reported species are H2O2 (`5.92e-9 mol`) and atomic N
(`2.99e-9 mol`). Small abundance alone does not establish a small
O(1D) hazard without reviewed rate bounds. The next intake should
prioritize these channels or a defensible collision upper bound, plus
nascent-state relaxation and energy disposition. The CO2 excited
fragment channels remain excluded, and no HCN amount or energy closure
claim follows from this first-collision calculation. The expanded
124-test focused regression suite passes; the Pavan/SciPy Jacobian
overflow warnings noted above persist in the short wet-flask test.

### Explicit finite-rate gas/water contact sensitivity (2026-09-27)

The closed-flask path previously repartitioned every supported neutral
species instantaneously before and after each spark substep. This is a
high-contact limit, not a measured transfer rate. A separate, opt-in
`finite_rate_contact` mode now requires an evidence-labelled liquid-side
volumetric coefficient `kLa` in `1/s`; it never derives that coefficient
from HCN yield. The default equilibrium mode remains unchanged. In the
finite-rate mode the pre-spark instantaneous repartition is skipped, and
each post-chemistry substep integrates the fixed-volume two-compartment
balance analytically:

`d n_liq / dt = kLa (K n_gas - n_liq)`, where
`K = kH(T) (RT/V_g) m_water` after consistent bar conversion. For a
substep with fixed temperature and water mass,
`n_liq(t+dt) = N K/(1+K) + [n_liq(t)-N K/(1+K)] exp[-kLa(1+K)dt]`.
Gas and liquid are updated by equal and opposite amounts. The same
declared `kLa` is presently applied to all neutral species in the local
Henry subset; this is a **lumped sensitivity assumption**, not evidence
that their interfacial transport coefficients are identical. A zero
coefficient is an exact no-transfer control. The aqueous acid-base
equilibrium still runs at substep boundaries, not simultaneously with
gas chemistry and contact. [A primary measurement of dilute HCN Henry
partition near 298 K](https://doi.org/10.1021/es1001192) supports the
importance of the equilibrium relation, but was made in phosphate
buffer and supplies no transfer rate for the modeled spark flask. The
project's selected Henry values retain their separate per-species
provenance in `backend/data/thermo/prebiotic_henry_water.json`.

The one-hour Barth **feed-only** pathway diagnostic was rerun with the
explicitly *hypothetical* `kLa = 1e-4 s^-1`; feed, source assumptions,
and chemistry were unchanged. No measured HCN amount was read or fitted.

| One-hour steps | Finite-rate aqueous HCN+CN- (mol) | Instantaneous-contact aqueous HCN+CN- (mol) |
| ---: | ---: | ---: |
| 48 | 1.126488e-11 | 4.631319e-11 |
| 96 | 1.116448e-11 | 5.064515e-11 |
| 192 | 1.072468e-11 | 4.889739e-11 |

The finite-rate 96-to-192 difference is `4.10%` relative to the 192-step
value, meeting the existing **5% one-hour numerical screen**. Primary
and charged material ledgers closed in all three finite-rate runs, with
zero failed charged substeps. A 60-second spark-off feed control formed
zero HCN in both phases and closed its ledger. At 60 seconds with the
spark on, predicted HCN remained below the model's `1e-12 mol` reporting
floor, so those trace values are not interpreted. The recurring Pavan/
SciPy finite-difference Jacobian overflow warning remains unresolved.

These runs demonstrate that the water-collected amount depends strongly
on an **unmeasured** transfer parameter and on its feedback through gas
chemistry. They do not validate `1e-4 s^-1`, the full-duration yield, or
the thermal ledger. A matched hydrodynamic/contact measurement or a
defensible independent transfer correlation is required before using
this mode for quantitative comparison. The next combined numerical gate
is full-duration convergence with gas chemistry, transfer, and aqueous
speciation coupled at the sampled phase; source exposure, state-specific
O(1D) chemistry, and spark/thermal energy accounting remain independent
blockers.

The 32 HCN-formation tests, including the focused finite-contact checks,
pass after this change. The broader regression run did not complete:
Cantera raised `RuntimeError: bad allocation` while loading its packaged
plasma mechanism in an unrelated dry-plasma test. The Windows committed
memory counter was `99.27/101.33 GB` during isolation, with only about
`10.34 GB` attributed to visible process private bytes. That is a
host-level memory-pressure finding, not a passing chemistry test; the
broader suite must be rerun once commit headroom is restored. The
finite-contact one-hour probes themselves completed with closed ledgers
before this separate test-run failure.

### Feed-only contact-rate sweep and rejected coarse run (2026-09-27)

`backend/scripts/probe_barth_feed_pathway.py` now accepts an optional,
explicitly hypothetical `--gas-water-kla-per-s` input while holding its
feed, spark assumptions, and HCN-assay firewall fixed. The reproducible
`backend/scripts/audit_barth_transfer_sensitivity.py` tests a declared
contact-rate grid plus instantaneous equilibrium. It checks matching feed,
duration, electrical basis, resolution, primary and charged material
ledgers, and charged substep failures. It marks a rejected numerical case
as a failure with **no accepted amount**, rather than skipping it in a
passing sweep. It does not load the experimental HCN assay.

The one-hour, feed-only observations include:

| Hypothetical liquid-side kLa | 48-step aqueous HCN+CN- | 96-step | 192-step | Interpretation |
| ---: | ---: | ---: | ---: | --- |
| `1e-5 s^-1` | `1.354e-12 mol` | not tested | not tested | Near the `1e-12 mol` reporting floor; no convergence claim |
| `1e-4 s^-1` | `1.126e-11 mol` | `1.116e-11 mol` | `1.072e-11 mol` | 96-to-192 change `4.10%` |
| `1e-3 s^-1` | **rejected** | `4.340e-11 mol` | `4.270e-11 mol` | 96-to-192 change `1.63%` |
| Instantaneous contact | `4.631e-11 mol` | `5.065e-11 mol` | `4.890e-11 mol` | Unmeasured high-contact limit |

The rejected `1e-3 s^-1` 48-step calculation produced
`-1.68755e-12 mol` CH2(S) in the raw time-averaged GRI solver state,
exceeding the conservative projection's `1e-12 mol` negative-species
limit. The same scenario completed at 96 and 192 steps with both
material ledgers closed. This is a **numerical resolution failure**,
not a measured transfer effect and not permission to loosen the clip
threshold. The Pavan/SciPy Jacobian overflow warning persists. The
33-test HCN-formation module passes with one OpenBLAS thread; the broader
suite still needs a retest when host committed-memory headroom permits.

These values are sensitivity to unmeasured kLa and source assumptions,
not an uncertainty interval or comparison to Barth's 925-minute assay.
The immediate experimental blocker remains a compatible wet-spark data
set with measured discharge exposure and gas-water contact, or direct
measurements of those quantities for the candidate reactor. A bounded
search of modern wet-spark HCN studies did not establish both in the
same cyanide-assayed run; nominal generator ratings and stirring alone
cannot supply them.

A one-thread rerun of the full six-module focused regression set still
stopped with pytest `MemoryError` while reporting a later test; Windows
committed memory was `100.63/101.33 GB`. The individually isolated
Cantera plasma test that had failed earlier now passes, and all 33
HCN-formation tests pass. The six-module run is **not** recorded as
passed. No system memory setting or unrelated user process was changed.

### Opt-in simultaneous gas/contact integration awaiting runtime validation (2026-09-27)

The finite-rate closed-flask mode now also accepts
`closed_flask_finite_contact_solver: coupled_bdf` when the Boltzmann
source and time-averaged forcing are selected. This maps the same
declared `kLa` and Henry amount ratio `K` into the existing stiff
gas/liquid ODE with a relaxation coefficient `kLa(1+K)`:

`kLa(1+K)[n_g-(n_g+n_l)/(1+K)] = kLa(K n_g-n_l)`.

That identity is the governing balance already used by the analytic
post-step option; no HCN assay or new reaction rate enters. The
simultaneous option is intended to test and reduce gas-reaction/contact
operator-splitting error, not to infer the unmeasured transfer rate.
The primary and charged branches share the declared coefficients,
and the closed-flask local atom/charge check now includes both gas and
liquid when coupled transport is active. The old post-step analytic
option remains the default; its previously recorded numbers still
require a post-change regression check. Tests were added for the
coefficient identity, unsupported spark-mode rejection, and a short
feed-only end-to-end material ledger.

**Runtime validation of the new coupled option is pending.** At this
point Windows committed memory rose to `101.28/101.33 GB`; a Python
import failed with `MemoryError` before the experiment could run.
The modified Python files pass `py_compile`, but that verifies syntax
only. Do not use the coupled option or its unmeasured `kLa` for an HCN
yield claim until the new tests, a feed-only negative control, and a
48/96/192-step numerical comparison pass after memory headroom returns.
No user application, page-file setting, or other system state was
changed to work around the host memory condition.

On continuation, committed memory remained about `101.11-101.17` of
`101.33 GB`. Three lightweight coupled-contact tests pass with one
OpenBLAS thread: algebraic coefficient equivalence, direct ODE versus
analytic two-compartment transfer, and rejection of an unsupported
spark mode. The end-to-end 60-second plasma/contact test exited before
producing a result under this host condition. Thus the full coupled
apparatus branch, its primary/charged ledger, and numerical convergence
remain **unverified**; these lightweight tests must not be promoted to
an HCN amount claim.

Static review then found that the shared coupled-transfer application
could clip up to `1e-11 mol` of negative liquid inventory, comparable to
the trace HCN signal. The new closed-flask coupled path now rejects a
negative liquid amount larger than `1e-18 mol` and records any smaller
clipped roundoff by species; the pre-existing recirculating-flow
tolerance remains unchanged pending its own audit. A focused regression
test for this strict behavior was added. It could not run: the next
pytest invocation failed while importing `conftest.py` with
`MemoryError`, and even a subsequent syntax-check process could not
start reliably. Treat the strict guard as **implemented but unverified**
until host commit headroom returns and all new tests pass.

### Coupled contact numerical screen and charged-ion repair (2026-09-27)

After Windows committed-memory headroom recovered, all five focused
coupled-contact tests passed, including the strict liquid-inventory guard
and a 60-second end-to-end material-ledger test. The full HCN module and
the six-module HCN/prebiotic/O(1D)/aqueous/CO2-O2 regression set passed.
The SciPy finite-difference Jacobian overflow warnings persist, so passing
tests do not settle that solver diagnostic.

The first one-hour, 48-step `coupled_bdf` run closed the primary ledger but
failed 44 charged substeps. The first failure was a raw `-1.865e-17 mol`
CH4+ amount while neutral CH4 was only `6.157e-21 mol`. The old charged
projection could use only the exact neutral parent or leave at most
`1e-17 mol` of explicitly reported roundoff. It now permits a secondary,
**at most `1e-15 mol` per ion**, correction across already-positive Pavan
species, only when every Cantera element and formal charge can be preserved
to the checked numerical tolerance. Larger or infeasible corrections still
fail closed. This is a numerical repair, not a CH4+ reaction, product
source, or yield fit. Direct tests cover a feasible carbon/hydrogen bath,
an unrepresentable bath, the existing neutral-parent projection, and the
upper bound. The same 48-step run then had zero failed charged substeps
and both whole-apparatus material ledgers closed.

The first 48/96/192-step coupled screen still rejected 96 steps: raw GRI
CH was `-1.91375e-11 mol`, exceeding the unchanged `1e-12 mol` clipping
bound. The coupled gas/liquid BDF now retries **only** a solver failure or
that negative-inventory failure with relative/absolute tolerances
`1e-10/1e-21 mol` and a maximum step of at most 1% of the substep,
instead of accepting a larger clip. Retry use and its initial failure are
reported. One substep required this retry in the completed screen.

`backend/scripts/audit_barth_coupled_contact.py` repeats the *feed-only*,
one-hour run at a **hypothetical** common liquid-side `kLa=1e-4 s^-1`.
Any rejected resolution fails the audit without publishing an accepted
HCN amount. With the bounded numerical repairs, both propagated branches
close their ledgers, have zero failed charged substeps, and meet the
predeclared <=5% 96-to-192-step endpoint screen:

| One-hour steps | Coupled primary aqueous HCN+CN- (mol) | Coupled charged aqueous HCN+CN- (mol) |
| ---: | ---: | ---: |
| 48 | 4.3704958e-12 | 4.6139891e-12 |
| 96 | 4.3699916e-12 | 4.6141333e-12 |
| 192 | 4.3701733e-12 | 4.6142065e-12 |

The 96-to-192 relative changes are `0.00416%` primary and `0.00159%`
charged. Maximum summed primary positive-species projection over a run is
`2.12e-12 mol`; direct HCN projection adjustment is **zero**. The former
is not an HCN error bar because corrections to other radical species can
affect later chemistry. Maximum cumulative explicitly unprojected ion
roundoff is `3.40e-18 mol`. A one-hour, 48-step *carbon-removed* control
with the **same coupled solver and spark** gives `0` gas HCN,
`1.74e-40 mol` primary aqueous HCN+CN-, and `0` charged aqueous HCN+CN-:
all below the `1e-12 mol` reporting floor, with both ledgers closed and
zero failed charged substeps. The short carbon-removed regression also
passes. The separate spark-off control remains on the post-step branch;
this schema does not run a Boltzmann-coupled BDF with the spark off.

The prior post-step analytic contact path still reproduces its 48-step
primary `1.1264883e-11 mol` result with both ledgers closed. At that same
resolution and declared `kLa`, simultaneous contact gives only
`4.3704958e-12 mol`. The methods disagree by about 2.6-fold despite each
having its own small 96-to-192-step change. This shows that within-method
endpoint refinement is not enough to declare **operator-splitting**
convergence or a physical transfer rate; fast gas HCN formation/loss and
boundary-only aqueous speciation need a separate timescale audit.

None of these numbers is an experimental yield prediction. The feed-only
card has no HCN assay value, the wet spark exposure and `kLa` remain
independently unmeasured, the run is one hour rather than the reported
925 minutes, O(1D) fragment channels remain excluded, the aqueous network
and thermal energy ledger are incomplete, and the Jacobian warnings remain.
The next numerical work is to diagnose the split/coupled discrepancy and
charged Jacobian warnings without relaxing ledgers. The next scientific
work is independent wet-discharge/contact constraints and justified
state-resolved/aqueous chemistry, followed by a comparable full-duration
held-out measurement; no known HCN yield may be used to tune the solver.

### Formaldehyde state-identity finding (2026-09-27)

A controlled property-subset diagnostic localized most of the 2.6-fold
post-step/coupled HCN difference to the CH2O gas-water channel. With only
the CH2O Henry entry temporarily withheld, the 48-step one-hour amounts
are `1.1471963e-11 mol` post-step and `1.1668220e-11 mol` coupled,
within 1.7%; both branch ledgers pass. This is **not** an acceptable
predictive removal of CH2O. At the declared scenario point, its effective
Henry ratio implies a 2.555-second relaxation time, faster than the
18.75-second finest apparatus step. More fundamentally, the selected
Sander CH2O constant is defined for free formaldehyde *plus methanediol*
in water, whereas the then-default inventory had only an aqueous `CH2O`
label. The property record and apparatus firewall flag that legacy default;
the later opt-in hydration state is summarized below. The governing
equations, primary kinetic provenance, pH
domain, and implementation tests are in
[the formaldehyde transfer gate](FORMALDEHYDE_TRANSFER_GATE.md).
No current HCN amount has been promoted to a calibrated yield.

The source-bounded formaldehyde hydration/contact component now passes
forward/reverse equilibrium, measured-domain rejection, and explicit
C/H/O/water conservation checks. It is also wired as an **opt-in**
closed-flask mode to both GRI and Pavan stiff integrations with separate
aqueous free CH2O and methanediol inventories. A seeded component stress
test and short carbon-removed control close both material ledgers. The
strict source-bounded one-hour feed-only attempt stops when modeled pH leaves
the measured rate-law range (pH 4.9967 after substep 17), so **no
source-bounded one-hour HCN amount** is available. A separately flagged
pH-4-to-5 extrapolation diagnostic passes the one-hour 48/96/192-step
numerical screen. Its 48-step primary/charged aqueous HCN+CN− amounts are
`1.1647704354e-11`/`1.1440061070e-11 mol`, with closed ledgers and zero
failed charged substeps. Carbon- and nitrogen-removed one-hour controls
return zero in both branches. This is **not** a measured yield validation:
the assay is held out, pH extension and contact rate are unvalidated,
and exposure is hypothetical. This one-hour screen alone does not test
the full 925-minute duration; a later feed-only stress run below does so
numerically, without validating its physical assumptions.
The charged source now responds to target depletion after the nitrogen-free
control exposed a negative-O2/atom-closure failure. Details are in
[the formaldehyde transfer gate](FORMALDEHYDE_TRANSFER_GATE.md). Source
exposure, contact-rate evidence, lower-pH kinetics, and thermal closure
remain missing; the HCN yield is not experimentally validated.

The opt-in inventory-limited primary source can now run with simultaneous
closed-flask gas/water contact, matching the charged branch's target-
depletion treatment. Its 48/96/192-step one-hour pH-extrapolation screen
passes with primary aqueous HCN+CN− amounts `1.16468060194e-11`,
`1.16468168396e-11`, and `1.16468408290e-11 mol`; the charged endpoints
are `1.14400610701e-11`, `1.14401137960e-11`, and
`1.14401387516e-11 mol`. Both ledgers close at each resolution, zero
charged substeps fail, and the separate one-hour carbon/nitrogen-removal
controls stay below the numerical reporting floor. This tests consistency
of source competition, **not** the HCN assay or physical `kLa`/spark
exposure. See [the formaldehyde transfer gate](FORMALDEHYDE_TRANSFER_GATE.md)
for the mode comparison and remaining limitations.

### Published-duration feed-only numerical stress (2026-09-27)

The [long-horizon gate](BARTH_LONG_HORIZON_GATE.md) now tests the same
Barth feed-only, closed-flask diagnostic over the reported 925 minutes
without loading the HCN assay. At hypothetical `kLa=1e-4 s^-1` and with
explicit pH-4-to-5 hydration extrapolation, 192/384/768-step primary
aqueous HCN+CN− amounts are `7.5562351e-11`, `7.5560521e-11`, and
`7.5559608e-11 mol`; charged amounts are `1.2877222e-10`,
`1.2875241e-10`, and `1.2874252e-10 mol`. Both material/charge ledgers
close, zero charged substeps fail, and full-duration carbon/nitrogen
removal controls remain below the `1e-12 mol` numerical floor. Tenfold
contact-rate changes move the 384-step primary endpoint by up to about
2.29-fold; the low-contact case stays within the strict measured pH
domain, while baseline and high contact use explicit extrapolation. This
is numerical/sensitivity evidence, **not** experimental yield validation:
spark exposure and contact remain unmeasured, the aqueous/thermal models
are incomplete, and the assay remains held out.

### Long-run numerical warning repair (2026-09-27)

The repeated SciPy Jacobian overflow warning was reproduced with warnings
treated as errors and traced to finite-difference perturbations of passive
or inactive coordinates in both charged Pavan and primary GRI integrations.
Both now use a bounded active-coordinate Jacobian without changing the
governing ODEs or balance tolerances. The one-hour 48-step and full 925-minute
384-step feed-only runs pass warning-as-error checks; the latter changes
the aqueous HCN+CN− endpoint by less than `3e-8` relative in each branch,
and both ledgers still close. See [the long-horizon gate](BARTH_LONG_HORIZON_GATE.md)
for before/after values. This removes a numerical warning, not the missing
wet-spark exposure, measured contact coefficient, low-pH aqueous kinetics,
thermal closure, or held-out assay comparison.

A bounded follow-up search found no directly transferable measurement of
wet-spark exposure and gas-water contact for the selected cyanide-assayed
Barth run; see the [long-horizon evidence check](BARTH_LONG_HORIZON_GATE.md).
Other experiments separately characterize generators or wet-discharge
electrical input, but lack a matching aqueous HCN assay and reactor feed.

The [full-duration source audit](BARTH_LONG_HORIZON_GATE.md) now confirms
that the primary CO2 feed has zero *direct accepted neutral-fragment*
electron-source extent in all 384 steps, while the CO2 -> CO(X)+O(1D)
collision channels accumulate `2.34473e-4 mol` of excluded hypothetical
opportunities. These are not product or yield moles. The O(1D) fate audit
is fully resolved only for the first substep; later bath partners lack
reviewed fate data. The next chemistry work is defensible state-resolved
CO2 dissociation and oxygen/CO fate, with explicit excitation/heat accounting.

An isolated source projection now conserves C/H/O and charge for a small,
fully covered thermalized-O(1D) event batch, using reviewed first-collision
branches and an explicit [NIST O(1D) electronic level](https://physics.nist.gov/PhysRefData/Handbook/Tables/oxygentable5_a.htm).
It does not feed the apparatus or close heat; incomplete later baths fail
closed. The [long-horizon gate](BARTH_LONG_HORIZON_GATE.md) also reports
two hypothetical unknown-rate stress points, not measured bounds.

An accepted-fragment thermochemical audit now passes on the two-substep feed
check and all 384 full-duration substeps. Ordinary H2O/O2 fragment source
costs exceed their GRI standard enthalpy floors. The N2 -> 2N outcome is
charged to its parent excitation manifold rather than directly; the smallest
ratio of energy-scaled high-energy N2 excitation *opportunities* to accepted
outcomes is `4.058`. This is a necessary consistency screen, not a traced
branching mechanism or thermal-energy closure. The [long-horizon gate](BARTH_LONG_HORIZON_GATE.md)
records the values and limitations; no HCN assay was loaded or fitted.

Reviewed NASA/JPL O(1D) fates for N2O and NH3 now enter the isolated,
atom-balanced first-collision projector; CH4 enters only the removal-hazard
audit because an isomer branch is unresolved. The full-duration maximum
unreviewed-bath mole fraction falls to `2.32424e-5`, but 383/384 snapshots
still fail closed, notably on HNCO and H2O2. No CO2/O(1D) product enters
the apparatus, and the HCN endpoint is unchanged. See the
[long-horizon gate](BARTH_LONG_HORIZON_GATE.md) for provenance and the
non-bounding stress comparison.

A full-duration, feed-preserving hypothetical electrical-coupling probe now
finds `3.596`-fold primary versus `1.229`-fold charged aqueous HCN+CN−
spans between coupling fractions 0.15 and 0.60. Both off-baseline cases
pass 192/384-step numerical checks and material ledgers. This is not a
measured coupling interval or HCN yield uncertainty: the primary branch
even decreases as coupling rises, and its reaction-flux cause has not yet
been isolated. The [long-horizon gate](BARTH_LONG_HORIZON_GATE.md) gives
the amounts, energy estimates, and limits. Independent spark exposure is
still required before a held-out assay comparison.

An opt-in full-duration GRI reaction-flux trace now reconciles signed-net
HCN reaction contributions with the observed gas chemical inventory at
all three coupling settings. At baseline, positive and negative signed
net contributions are `5.41971020e-7` and `5.41890535e-7 mol`, leaving
only `8.04848e-11 mol` net; reaction 275 (`CH3 + N <=> H2 + HCN`) and
reaction 236 (`H + HCN (+M) <=> H2CN (+M)`) dominate the opposite sides.
The small surplus shrinks as hypothetical coupling rises, explaining the
primary branch's inverse response *within GRI*, not validating its rates
or real wet-spark behavior. The [long-horizon gate](BARTH_LONG_HORIZON_GATE.md)
contains the three-case quadrature closure and limitations. This strong
near-cancellation adds kinetic-mechanism uncertainty to the still-missing
independent discharge/contact and thermal constraints.

The subsequent [full-duration low-temperature kinetics audit](BARTH_LONG_HORIZON_GATE.md#low-temperature-gri-hcn-kinetics-challenge-2026-09-27)
found a larger structural blocker: GRI's dominant direct `N + CH3 -> HCN + H2`
branch is contested by [state-resolved 298 K work](https://arxiv.org/pdf/1902.05574),
and its H+HCN association rate at the assumed feed state exceeds a
[published 300 K theoretical high-pressure limit](https://quantchem.kuleuven.be/cmms/Articles/JPCA/jp981542u.pdf)
by about `4.47e4`. GRI also omits the source-supported
`H2CN + H -> HCN + H2` conversion. In opt-in **diagnostic** 925-minute
384-step stresses, removing direct N+CH3 formation gives
`4.94095141e-12 mol` primary aqueous HCN+CN−, versus unchanged
`7.55605215e-11 mol`; multiplying H+HCN association by `1e-5` gives
`4.26850386e-7 mol`. All material ledgers close, but these are neither
validated alternative mechanisms nor measured uncertainty bounds. This
greater-than-one-order-of-magnitude kinetic instability, along with the still-open
thermal balance and unmeasured spark exposure, blocks quantitative
validation. Do not open the held-out HCN assay or fit these rates to it;
first build/check a state-resolved, pressure-consistent low-temperature
HCN/H2CN subnetwork from independent kinetic evidence.

A second [source-form audit](BARTH_LONG_HORIZON_GATE.md#independent-h--hcn-association-applicability-check-2026-09-27)
now compares GRI's feed-state `H + HCN` association with Tsang and Herron's
1991 NIST-evaluated estimate. The GRI effective rate is `105.56` times the
source's **estimated high-pressure expression** at the hypothetical 298.15 K
feed. The NIST authors reported no direct measurement, the source begins at
300 K and lists a large low-temperature uncertainty, and its N2/CO2
low-pressure terms omit other bath gases. This strengthens the diagnosis of
an unqualified HCN sink without supplying a replacement kinetic law or a
measured rate bound. The audit is read-only and assay-free.

An opt-in [H2CN fate audit](BARTH_LONG_HORIZON_GATE.md#h2cn-fate-and-partial-low-temperature-branch-overlay-2026-09-27)
now compares independent published low-temperature H2CN+H and H2CN+N
branches without changing the default mechanism or loading the held-out
assay. The reported-duration numerical result shifts by roughly four
orders of magnitude across those diagnostic kinetic profiles, despite
192/384-step convergence, closed atom/charge ledgers, and zero-cyanide
carbon-/nitrogen-removed controls. This is a mechanism-identifiability
failure, not an experimental yield estimate. The user-approved LXCat
collision files are already private and available; more copies of the
same files will not resolve the missing HNC fate, excited N(2D) pathways,
pressure-dependent H+HCN rate, independently measured discharge exposure,
or open apparatus thermal budget. The next gate is to build and test a
state-resolved low-temperature HCN/HNC/H2CN subnetwork from independent
kinetics, then bound the discharge and compare against a comparable
assayed reactor without using that assay to set model rates.
The first HNC thermochemistry intake check found a large discrepancy:
the imported 2024 supplementary HNC record evaluates to `638.148 kJ mol^-1`
at 298.15 K, versus [NIST's evaluated `194.33 ± 9.00 kJ mol^-1`](https://cccbdb.nist.gov/exp2x.asp?casno=6914074).
It remains outside the inventory-driving primary GRI phase. The
non-propagating 2024 neutral shadow does contain it and now flags
HNC-sensitive rates as unvalidated. Do not transfer it to the primary
solver for reverse rates or energy until the ground-state record is
reconciled. See the long-horizon gate's HNC thermochemistry intake check.
A NIST-grounded 290–310 K HNC candidate and isolated H2CN+H exit-channel
enthalpy audit now pass their reference-state and isomer-gap tests. This
resolves a narrow room-temperature **thermochemistry** intake problem, not
the kinetic or apparatus gate: the measured H2CN+H total is only a lower
bound, the HCN/HNC split is estimated, H+HNC isomerization is uncertain,
and HNC gas-water fate has not been validated. The candidate is therefore
not part of the default inventory-driving solver.
The reference entropy is converted from NIST's 1 bar to GRI's 1 atm
standard pressure before constructing the diagnostic species.
An opt-in [HCN/HNC competing-exit diagnostic](BARTH_LONG_HORIZON_GATE.md#opt-in-hcnhnc-competing-exit-diagnostic-2026-09-27)
now couples the NIST-grounded HNC species to the gas inventory with both
estimated Hébrard `H2CN + H` exits, restricted to 290–310 K. Its 192/384-step
full-duration structural-stress results are numerically stable and conserve
material, but HNC gas-water transfer, validated H+HNC conversion, a full HNC kinetic
network, independent spark exposure and complete thermal accounting remain
unresolved. The default mechanism and held-out HCN assay are unchanged;
this is not quantitative yield validation. The user-approved LXCat
collision files remain private and already available.
An additional [opt-in estimated H+HNC isomerization sensitivity](BARTH_LONG_HORIZON_GATE.md#estimated-relaxed-hnc-h-atom-isomerization-sensitivity-2026-09-27)
uses Hébrard et al.'s `4.03e-11 exp(-1200/T)` expression, whose reported
300 K uncertainty is a factor of ten. In the same full-duration structural
stress, the 384-step aqueous HCN+CN− amount rises by `4.68%` relative to
the no-isomerization HNC diagnostic, while net gas HNC falls about
`5.05e4`-fold. This strengthens the need to distinguish isomer fates but
does not resolve HNC water transport, nascent-state relaxation, overall
mechanism uncertainty, or the independent spark/thermal constraints.
An opt-in, read-only [glycolonitrile endpoint screen](BARTH_LONG_HORIZON_GATE.md#conditional-aqueous-glycolonitrile-endpoint-screen-2026-09-27)
now checks one omitted aqueous cyanohydrin reservoir. Under an explicitly
conditional total-formaldehyde and frozen-pH interpretation, the unchanged
925-minute feed-only endpoint would bind 32.83% of its aqueous HCN+CN−.
This is not a coupled product prediction or HCN-yield correction. The
original equilibrium pool convention is supported by indexed primary-paper
text, but its PDF was not visually verified; source-bounded trajectory
kinetics remain unavailable before any inventory-driving aqueous implementation.
The subsequent [rate-aware correction](BARTH_LONG_HORIZON_GATE.md#glycolonitrile-equilibrium-versus-timescale-correction-2026-09-27)
uses the original paper's distinct free-formaldehyde kinetic pool. Its
read-only 384-step gross forward opportunity is only `2.55507e-14 mol`
under the feed-only model, while **all 384 cyanide-bearing steps are
outside the measured pH range**. Thus the earlier 32.83% *equilibrium*
potential is not a 925-minute yield prediction or a reason to inject
glycolonitrile into the current solver. It does not remove the independent
spark, gas-kinetic, contact, or thermal blockers.

### Private LXCat N2+ dissociative-recombination intake (2026-09-27)

The user approved LXCat access. The **eMol-LeHavre** electron + N2+ record
was obtained from the official LXCat browse/download interface on 2026-09-27
and retained only under `backend/data/collision/private/`, which is ignored by
version control. It is hash-pinned (`SHA-256
00ad61fca159d92a5adb00f7f80fb2cb9cf463e7431c9bb7b6211747f40e2967`).
The contributor identifies its source as [Little et al., Phys. Rev. A 90,
052705 (2014)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.90.052705);
LXCat citation: *eMol-LeHavre database, www.lxcat.net, retrieved September 27,
2026*. The source file is **not redistributed**. In BOLSIG format this is
classified `ATTACHMENT`, but the contributor comment describes physical
dissociative recombination. Its numeric energy support is only
`0.00001–0.39999 eV`, and it gives no N(4S)/N(2D)/N(2P) product branch or
ion-vibrational population. The experimentally reported product-state
fractions in [Peterson et al. (1998)](https://ir.amolf.nl/pub/4894) are
near-zero-energy, v=0-specific observations, not a branch law for the whole
nonthermal spark EEDF.

The raw 39,999-row export contains **8,100 repeated rounded energy entries**,
8,098 of them with unequal cross sections. A strict increasing-grid import
correctly failed. The new read-only audit retains the untouched export and
integrates separate minimum/maximum cross-section curves at each rounded
energy. This is an export-precision envelope, **not** a statistical error bar.
It neither enters the EEDF energy-loss calculation nor injects N atoms into
the gas network. `backend/scripts/audit_n2plus_recombination.py` reads only
the previously documented feed-only gas composition and assumed 200 Td field;
it never reads the held-out HCN assay or a known-product recipe.

For that conditional feed, the BOLOS EEDF mean energy is `5.3245 eV` on 400
cells and `5.3062 eV` on 4,000 cells. Yet the N2+ record reaches only
`0.39999 eV`: at 4,000 cells the modeled EEDF begins at `0.0075 eV`, so even
part of the record's low-energy interval remains unresolved. The fraction of
the tabulated EEDF in the common `0.0075–0.39999 eV` interval is `0.02635`.
The **partial** integrated coefficient in that interval is
`[2.22295, 2.23816]e-15 m3/s` at 4,000 cells, versus
`[1.78857, 1.79045]e-15 m3/s` at 400 cells. The changing lower grid edge
explains why these are not a converged *total* rate; the unmeasured
`>0.39999 eV` collision curve and post-discharge EEDF are still missing.
These numbers are not an N-atom source or an HCN-yield estimate. The tested
next step is to source a wider, state-/vibration-resolved DR curve and a
time-resolved electron distribution, then add explicitly balanced excited-N
reactions only when their kinetics are independently supported. Independent
discharge exposure and the kinetic, transport, and thermal blockers above
remain open.

Verification: the new private-record and no-extrapolation tests pass, the
focused 14-test collision/plasma subset passes, and the full backend `pytest`
suite exited successfully on 2026-09-27. The read-only diagnostic did not
alter any gas, water, charge, or energy inventory state.

### Reduced-route conditional numerical audit (2026-09-29)

The separate [reduced HCN route](REDUCED_HCN_AND_ATOM_CEILING.md) now has an
optional assay-free audit exposed in Origins. It repeats an input scenario at
maximum time steps of at most duration/16 and duration/32, checks gas and
aqueous HCN-equivalent plus all source extents against 1% relative and
`1e-12 mol` absolute tolerances, and reports source-off, carbon-removal and
nitrogen-removal controls. A control with no remaining gas phase is **incomplete**,
not a pass. The synthetic CH4/N2 source fixture passes; a seeded HCN-only
fixture correctly reports incomplete. The audit adds no source data, does not
alter the wet plasma model, and does not validate a real HCN amount.

The next scientific input is unchanged: an independently constrained,
state-resolved fragment source for the wet spark plus contact/thermal data.
This project's HCN assay values are not inputs to the audit or to the reduced
runner. During development, a broad file search did expose assay-file text
to the developer; no value was used for parameter selection or comparison.
Any future claim of **blind** validation against that record should use a
frozen model and an independent evaluator or a fresh holdout, not just the
absence of assay reads inside the simulation code.
