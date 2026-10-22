# Formaldehyde transfer and hydration gate (2026-09-27)

## Why this is a separate scientific gate

The wet Barth feed-only probe gave `1.1264883e-11 mol` aqueous HCN+CN-
after one hour with post-chemistry Henry contact, but `4.3704958e-12 mol`
with contact simultaneous with gas chemistry (both at 48 steps and the same
**hypothetical** `kLa=1e-4 s^-1`). The diagnostic-only
`scripts.audit_barth_contact_species_timing --exclude CH2O` leaves every
other Henry entry and the feed unchanged; the amounts become
`1.1471963e-11` and `1.1668220e-11 mol`, respectively, a 1.7% difference.
Both branches and their material ledgers pass. Removing CH2O is **not** a
valid predictive model; this is evidence that formaldehyde handling drives
the 2.6-fold full-model discrepancy.

For 0.95 L gas, 50 mL water, 298.15 K, and the declared `kLa`, the current
effective CH2O Henry entry gives liquid/gas amount ratio `K=3913.55` and
linear relaxation time `1/[kLa(1+K)]=2.555 s`. The finest one-hour apparatus
step tested is `18.75 s`, so delayed transfer is not time-resolved against
that relaxation. HCN's analogous relaxation time is about `600 s`; this
helps isolate the issue to a fast precursor pool, not the HCN transfer
equation itself. The bare two-compartment coupled ODE already matches the
analytic transfer update; the difference arises with chemistry active.

More importantly, [Sander's 2023 Henry compilation methanal entry](https://henry.mpch-mainz.gwdg.de/henry/casrn/50-00-0)
explicitly marks the selected ∼`3.0e3 mol/(kg water bar)` value as
**effective**, defined by `([HCHO] + [CH2(OH)2])/p(HCHO)` (note 457).
The current apparatus stores that entire aqueous pool under the molecular
`CH2O` label and does not distinguish free formaldehyde from methanediol.
This is a species-identity and kinetic-model defect, not just uncertainty in
the time-step size. The property record and apparatus firewall now flag it.

## Governing theory for the next implementation

The reaction to represent explicitly is
`CH2O(aq) + H2O(l) <-> CH2(OH)2(aq)`. The
[primary stirred-cell absorption study by Winkelman et al. (2002)](https://doi.org/10.1016/S0009-2509(02)00358-5)
reports, for 293-333 K and pH 5-7,

`k_h(T) = 2.04e5 exp(-2936/T) s^-1`,
`K_h(T) = exp(3769/T - 5.494)`, and `k_d(T) = k_h(T)/K_h(T)`.

At 298.15 K these are approximately `10.79 s^-1`, `1271`, and
`0.00849 s^-1`. These are *aqueous* hydration/dehydration parameters, not
gas reaction rates or gas-liquid contact coefficients. The study measured
chemically enhanced absorption and explicitly separated interfacial
transfer from liquid reaction. A [later primary dilute-solution NMR study](https://doi.org/10.1021/jp513020y)
confirms that free molecular formaldehyde is only a trace in water and
studied pH 1.8-7.8 and 278-333 K; no numeric rate from that paper has yet
been imported into the solver.

In the dilute equilibrium limit, if the effective Sander constant includes
only free HCHO plus methanediol, an intrinsic free-molecule Henry constant
would satisfy `H_eff = H_free (1+K_h)`. Combining the project's `3000`
effective value with the Winkelman `K_h(298.15 K)` gives a *derived*,
cross-source estimate `H_free ≈ 2.36 mol/(kg water bar)`. This is a
model-construction estimate, **not** a new measurement; temperature-law and
source consistency require sensitivity analysis. The existing effective
value must not simply be relabeled as `H_free`.

For gas amount `n_g`, free aqueous amount `n_f`, hydrate amount `n_h`,
liquid-water amount `n_w`, and an *independently supplied* free-species
contact coefficient, the candidate mass balances are:

```
J = kLa_free (K_free n_g - n_f)
r_h = k_h n_f - k_d n_h
d n_g/dt = gas_chemistry - J
d n_f/dt = J - r_h
d n_h/dt = r_h
d n_w/dt = -r_h
```

`K_free = H_free RT m_water/(1e5 V_g)` in the project's bar-based units.
These equations conserve C, H, O and charge when the water term is retained.
They do **not** determine `kLa_free`: interfacial area, stirring, gas-side
resistance, and chemical enhancement remain reactor-dependent. Assigning the
same `1e-4 s^-1` used in the earlier lumped scenario would be a labeled
sensitivity assumption, not an experimentally measured contact rate.
The source paper reports hydration enthalpy values spanning roughly
`-21.4` to `-39.4 kJ/mol`; no single thermal value or heat-flow closure is
adopted here without a compatible source and apparatus heat balance.

The current feed-only simulated ideal pH changes from about `5.63` to
`4.81` over the one-hour coupled run. Thus its endpoint already falls
outside the 2002 study's stated pH 5-7 rate-correlation domain. An
opt-in implementation must expose this domain violation and must not
quietly extrapolate that correlation through the full 925-minute assay.

## Required implementation and validation before an HCN claim

1. Add distinct aqueous free-CH2O and methanediol states, with the latter's
   `C1H4O2` composition and water consumption, while retaining gas CH2O.
   Keep Sander's effective total-pool value out of the free-molecule flux.
2. Couple free-CH2O gas transfer and measured hydration/dehydration kinetics
   inside the stiff gas/liquid integration. Independently label any assumed
   `kLa_free` and check the temperature/pH evidence domain each step.
3. Verify isolated transfer, isolated hydration, equilibrium, reverse
   desorption, exact material and charge closure, and separate source-off and
   carbon-removed controls. Compare with independent free/total formaldehyde
   measurements if available; HCN yield must not tune this subsystem.
4. Repeat HCN resolution and contact-rate sensitivity with the full aqueous
   state. Account for hydration enthalpy in a genuinely closed thermal
   budget before interpreting product per deposited joule.

Until then, the current one-hour HCN amounts are **diagnostic scenario
outputs**, even though their numerical ledgers and within-method resolution
checks pass. Spark exposure and contact hydrodynamics remain independent
unmeasured blockers to a comparable experimental yield test.

## Source-bounded isolated implementation check (2026-09-27)

`backend/data/kinetics/formaldehyde_hydration_winkelman_2002_v1.json`
records the 2002 primary rate and equilibrium correlations and their
293-333 K, pH 5-7 domain. The new
`backend/app/prebiotic/formaldehyde_hydration.py` evaluates them with a
fail-closed domain check (or an explicit extrapolation flag), derives the
*free* CH2O Henry estimate from the Sander effective total-pool entry, and
calculates gas/free-CH2O/methanediol/water changes. It also declares a
reserved, zero-gas-transfer `CH4O2_AQ` inventory identity with C1H4O2
composition for subsequent solver wiring. It accepts an independently
supplied `kLa` only; no HCN assay or product amount is read or fitted.

At 298.15 K and pH 5.5, tests reproduce `k_h ≈ 10.79 s^-1`,
`K_h ≈ 1271`, `k_d ≈ 0.00849 s^-1`, and the cross-source
`H_free ≈ 2.36 mol/(kg water bar)`. With the **assumed** 0.95 L gas,
2.7 mol liquid water and `kLa = 1e-4 s^-1`, the candidate free-species
liquid/gas amount ratio is `≈ 2.994`, not the old effective-pool ratio of
roughly 3913. Both all-gas and all-hydrate isolated BDF starts approach
the expected free/gas and hydrate/free ratios by 100,000 s while conserving
C, H, and O. Isolated hydration also conserves water explicitly and supports
reverse dehydration. This long transient is a numerical component test, **not**
a measured apparatus contact time. The new tests and the HCN-formation module
pass together; the pre-existing SciPy finite-difference Jacobian overflow
warnings still occur in two wet/flask tests.

At this isolated-component checkpoint, the calculation was not yet active
in the apparatus. The opt-in coupled-apparatus extension is reported below.
The common `kLa` remains hypothetical, and source and thermal-energy
limitations are unaffected.

## Opt-in coupled-apparatus hydration check (2026-09-27)

The closed-stirred-flask apparatus now accepts
`closed_flask_formaldehyde_hydration: winkelman_2002_source_bounded` **only**
with simultaneous finite-rate BDF contact and aqueous acid-base pH enabled.
It is off by default. The aqueous `CH2O` state is now the free molecule in
this mode; the separate `CH4O2_AQ` state is methanediol. Inside both GRI and
Pavan stiff gas/liquid solves, free CH2O transfers using the cross-source
derived *free* Henry constant while CH2O + H2O reversibly hydrates. The
hydration extent changes the liquid-water inventory, and the distinct
hydrate participates in both whole-apparatus atom ledgers. The chosen
`kLa=1e-4 s^-1` remains a declared, **hypothetical** common contact-rate
sensitivity, not a measurement or HCN fit.

The opt-in 60-second Barth **feed-only** diagnostic completes with both
material ledgers closed and zero failed charged substeps, but its methanediol
and HCN amounts are below the project's `1e-12 mol` reporting floor and are
not interpreted as product yields. An intentionally *seeded aqueous CH2O*
component stress test (not an unseeded formation test) converts about
`9.992e-7 mol` of an initial `1e-6 mol` free-CH2O pool to distinct aqueous
methanediol in 60 s; the primary and charged ledgers both close, and the
charged branch has zero failed substeps. Repeating this component stress
test at 2, 4, and 8 substeps gives primary hydrate amounts
`9.99199676401246e-7`, `9.99199676401506e-7`, and
`9.99199676401675e-7 mol`; both branches retain closed ledgers and zero
charged failures. This is an isolated-component numerical screen, **not**
convergence of an unseeded HCN yield. The carbon-removed, spark-on 60 s
control remains at zero HCN and zero methanediol, with both ledgers closed.
These observations test state wiring and conservation, not a reactor rate.

The strict 1500 s/20-step feed-only run stops after substep 17 at simulated
ideal pH `4.9967`, just outside the primary study's pH 5-7 rate-law range;
it returns **no accepted longer-duration HCN amount**. Source-bounded rates
are fixed over each substep using boundary pH; they are not continuously
re-evaluated inside the stiff ODE. A compatible lower-pH kinetic source or
an independently justified model is needed before the 925-minute wet-spark
assay can be compared. The formaldehyde component/apparatus tests and a combined
HCN/prebiotic/aqueous-speciation regression set pass. SciPy's existing
finite-difference Jacobian overflow warnings persist in some wet runs.

This upgrade resolves the *identity* mismatch in the opt-in path but does
not establish physical `kLa`, discharge exposure, excited/ion product fate,
aqueous-network completeness, or hydration and apparatus heat balance.
The effective-pool default remains available for legacy numerical
diagnostics and is still explicitly flagged as unvalidated for quantitative
HCN use. No HCN measurement has been used as a solver input.

The opt-in energy report now separately gives net methanediol formation,
`delta_H = -R * 3769 K ≈ -31.34 kJ/mol` from the primary study's Van't
Hoff equilibrium slope, and the corresponding isothermal heat release.
For the intentionally seeded `1e-6 mol` component stress test, this term is
about `0.0313 J`. It is **not** fed back into temperature and does not close
electrical input, gas reaction enthalpy, evaporation/condensation, wall heat,
radiation or other aqueous-reaction heat. The full thermal-energy flag stays
false.

### Lower-pH evidence intake, not yet a solver parameter

The [Winkelman et al. primary paper](https://pure.rug.nl/ws/files/3619450/2002ChemEngSciWinkelman.pdf)
reports that Schecker and Schulz (1969) measured
`k_h = 7800 exp(-1913/T) s^-1` over 298-333 K and pH 4-7, which would
give about `12.75 s^-1` at 298.15 K versus Winkelman's `10.79 s^-1`.
The [original 1969 study](https://doi.org/10.1524/zpch.1969.65.1_4.221)
has been identified but its method and complete data have not yet been
inspected directly. A [2015 NMR primary study](https://doi.org/10.1021/jp513020y)
measured pH-4 hydration and dehydration in **D2O**, not the modeled H2O;
it cites older H2O values with a substantial isotope effect. None of these
lower-pH numbers has been imported into the predictive solver. The immediate
research task is to verify the original H2O measurements, equilibrium
compatibility and uncertainty before replacing the fail-closed pH-5 limit.

A bounded follow-up checked the [author-uploaded 2015 primary NMR
article](https://www.researchgate.net/publication/273147900_NMR_Studies_of_the_Equilibria_and_Reaction_Rates_in_Aqueous_Solutions_of_Formaldehyde).
Its own numeric pH-4 hydration/dehydration pair was measured in D2O;
its table quotes the 1972 Sutton–Downes H2O values (`9.8 s^-1` hydration,
`5.1e-3 s^-1` dehydration) as **earlier work**, not new H2O measurements.
The [1972 original article's accessible publisher abstract](https://pubs.rsc.org/en/content/articlelanding/1972/c3/c39720000001)
confirms an irradiated aqueous-methanol flow method but does not show
the numeric conditions. Moreover, a [2013 literature table](https://acp.copernicus.org/articles/13/11625/2013/acp-13-11625-2013.pdf)
labels that `9.8 s^-1` point pH 7-9 and 295 K, not pH 4 and 298 K as
the 2015 table states. This pH/temperature provenance conflict prevents
using the quoted H2O number to lift the strict pH-5 guard. The full 1972
and 1969 primary methods remain the specific evidence needed. Bell and
Evans' [1966 primary abstract](https://doi.org/10.1098/rspa.1966.0097)
also reports general acid/base catalysis of dehydration, so buffer identity
cannot be assumed irrelevant when transferring low-pH rates to the flask.
Even taking the two quoted H2O rates together would imply an equilibrium
hydrate/free ratio near `1922`, versus `1271` from the current Winkelman
fit at 298.15 K; those mismatched conditions are not a justified way to
retune the free-species Henry estimate.

### One-hour pH-extrapolation diagnostic and source-depletion repair

An explicitly named diagnostic mode extends the unchanged Winkelman H2O
rate expression only through pH 4-5, outside this fit's measured domain.
The strict `winkelman_2002_source_bounded` mode still stops at pH 5. The
diagnostic fails closed below pH 4 or outside 293-333 K, records minimum
modeled ideal pH and marks extrapolation. This is a sensitivity probe, not
a newly validated kinetic parameter or an HCN yield claim.

At hypothetical `kLa=1e-4 s^-1`, the one-hour 48/96/192-step normal-feed
diagnostic passes its <=5% 96-to-192 numerical screen in both branches,
with closed ledgers and zero failed charged substeps. The 48-step primary
aqueous HCN+CN− endpoint is `1.1647704354e-11 mol`; the charged endpoint
after the source-depletion repair is `1.1440061070e-11 mol`. Minimum
modeled ideal pH is `4.81260`, so extrapolation is flagged. The HCN assay
was never loaded. Neither spark exposure nor contact rate is independently
constrained, and this is not experimental yield validation.

The one-hour 48-step carbon-removed control remains at zero HCN in both
branches. A nitrogen-removed control (N2 replaced by H2 at fixed total
gas moles) initially failed charged-branch atom closure: a fixed O2
collision sink consumed O2 already depleted by the evolving chemistry,
and clipping would have created about `1.13e-9 mol` O atoms. Tighter BDF
steps and tolerances did not reduce the `5.66e-10 mol` negative O2 state.
The charged collision source now scales with remaining target inventory,
and an integrated state records the realized extent separately from the
inventory-capped nominal opportunity. A bounded retry still rejects
material neutral negativity; no atom threshold was relaxed. The
nitrogen-removed control now has zero aqueous HCN+CN− in both branches,
both ledgers closed, and zero charged failures. A dedicated long-control
regression test covers this. SciPy finite-difference Jacobian overflow
warnings persist. Collision energy allocated to events prevented by
target depletion is not redistributed to other channels or converted into
a closed gas/wall heat account; the full heat balance remains open.
The charged shadow now also reports the independent threshold-energy lower
bound of its **realized** electron-impact events and rejects a value above
the shared collision allocation. Outcome-only fragment labels contribute
zero additional threshold cost because their parent excitation manifold
already carries it. This is a collision-budget check, not calorimetry or a
closed plasma/thermal-energy equation.

### Competing-source timing check in the wet flask

The earlier normal-feed run applied a fixed, cycle-averaged primary GRI
fragment sink even as heavy chemistry depleted its targets. An opt-in
`primary_source_forcing: inventory_limited` mode is now permitted with the
coupled wet-flask contact solver. It applies each published collision
opportunity as a first-order target-limited rate, integrates its realized
extent, and rejects a primary source that exceeds the allocated collision
energy. The hazard uses the **accepted inventory-capped** nominal extent,
not the larger pre-cap opportunity. The old `time_averaged` mode remains
available for controlled comparison; neither mode changes the EEDF,
spark geometry, contact coefficient, or held-out assay.

With the same diagnostic pH extension and hypothetical `kLa=1e-4 s^-1`,
the inventory-limited primary/charged one-hour aqueous HCN+CN− amounts are:

| Steps | Primary (mol) | Charged (mol) |
| ---: | ---: | ---: |
| 48 | `1.16468060194e-11` | `1.14400610701e-11` |
| 96 | `1.16468168396e-11` | `1.14401137960e-11` |
| 192 | `1.16468408290e-11` | `1.14401387516e-11` |

The 96-to-192 relative changes are `2.06e-6` primary and `2.18e-6`
charged. Both material ledgers close, no charged substep fails, and the
primary resolved-fragment energy stays within its allocation. At 48
steps the primary endpoint differs from the fixed-rate result by about
`0.0077%`, which is a *within-model source-timing sensitivity*, not an
experimental uncertainty estimate. The one-hour nitrogen-removed control
has zero aqueous cyanide in both branches. The carbon-removed control has
zero charged cyanide and `1.21e-40 mol` primary cyanide, far below the
project's `1e-12 mol` numerical reporting floor; both control ledgers
close. The model still extrapolates below pH 5 (minimum ideal pH `4.8126`),
retains the unresolved Jacobian warnings, and lacks a measured wet spark
exposure/contact coefficient and a closed thermal budget.
