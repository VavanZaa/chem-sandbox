# Conditional reduced HCN route and atom ceiling (2026-09-29)

This is a **separate, source-driven alternative** to the Origins EEDF/spark
calculation, not a corrected Miller–Urey yield. It deliberately skips the
electron Boltzmann solve and spatial spark geometry. It can only calculate a
conditional HCN trajectory after the user supplies effective gas-fragment
conversion hazards, their product-state identities, gas/water contact, Henry
solubility and pH with explicit provenance. It does not derive these from an
HCN assay. The UI's methane-source hazard is **zero by default** and it has
no nitrogen-fragment source. In particular, the reviewed N2 total collision
outcome cannot be entered as a ground-state `N2 -> 2 N` source. No new reaction rate or yield from
the literature has been imported into the predictive mechanism here.

## Governing approximation

The gas is an isothermal, fixed-volume, well-mixed Cantera GRI-3.0 phase.
SciPy BDF advances species amounts and integrated source extents. For each
externally supplied source channel `parent -> products`, the effective
whole-reactor source rate is `k_effective * n_parent` in mol/s; its declared
products must be GRI-compatible, charge/atom balanced, and cannot directly
inject HCN or any C–N species. This conditional hazard is **not** an electron
cross section, measured pulse waveform, or source rate automatically inferred
from deposited joules. Gas pressure follows `n_total RT/V_gas`.

Only neutral HCN moves between gas and one water compartment:

`J = kLa V_water (Hcp p_HCN - alpha_HCN n_cyanide,aq/V_water)`,

`alpha_HCN = 1/(1 + 10^(pH-pKa_HCN))`, and
`p_HCN = n_HCN,gas RT/V_gas`.

The `pKa_HCN(T)` comes from the project's reviewed aqueous subset. The
Henry coefficient, `kLa`, and **fixed** pH are explicit inputs, not fitted by
the tool. The aqueous amount is one HCN-equivalent pool containing molecular
HCN and CN−; its dissociation proton is lumped into that equivalent. This
allows an atom ledger but is **not** an explicit aqueous charge or reaction
network. Other gases do not dissolve in this reduced route. The model holds
temperature constant without resolving heat flow or evaporation.

If every source supplies a minimum event-energy cost and a deposited-energy
amount is supplied, the runner rejects a source inventory whose minimum
fragment energy exceeds that deposit. This is a necessary lower-bound
screen, not a closed first-law balance. It does not silently scale source
rates to match energy. The result reports mechanism hash, Cantera version,
source extents, gas/aqueous HCN, atom residuals, numerical controls and all
input provenance, including which fields are hypothetical. Every output is
explicitly conditional and unvalidated. A single run does not itself pass a
solver-step convergence check; the synthetic regression has a separate
5 s/2.5 s refinement comparison.

The separate atom ceiling requires no gas or plasma calculation. One HCN or
CN− equivalent requires one C and one N, so the total-family ceiling is
`min(n_C_atoms, n_N_atoms)`; molecular HCN additionally needs one H, giving
`min(n_C_atoms, n_N_atoms, n_H_atoms)`. The apparatus preflight includes all
scheduled additions; sampling and competing products only lower what can
actually be collected. This is an intentionally generous **upper bound**,
not predicted conversion, collection, probability, or an energetic bound.
The reduced-route result applies the ceiling to its **initial reactive gas
only** because its water is an inert solvent; the full apparatus preflight
counts atoms in every initial compartment and scheduled addition. These
different scopes are labelled in the returned records.
For the existing Barth `220725001` **feed-only card** at its declared
0.95 L, 1 bar and 298.15 K initialization, the full-apparatus preflight
counts `0.00218822 mol` carbon atoms and `0.0702923 mol` nitrogen atoms;
its HCN+CN− element ceiling is therefore `0.00218822 mol`, limited by
carbon. This does not use the held-out HCN assay, and it says nothing about
whether any measurable fraction of that carbon reacts or reaches water.

## Verified software gates, not chemistry validation

The synthetic regression fixture uses hypothetical CH4 and N2 cleavage
hazards of `1e-3 s^-1` each, with a **hypothetical** all-ground-state N2
product. At 298.15 K and 60 s it exercises the numerical pipeline and
produces approximately `1.899e-6 mol` HCN-equivalent in gas plus water,
including `6.706e-7 mol` in water. Those numbers are not experimental
predictions. Atom residuals are below `1e-15 mol`; carbon-removed,
nitrogen-removed and zero-source controls remain below the `1e-12 mol`
software reporting floor. The 5 s versus 2.5 s maximum-solver-step test
agrees to within 1%. Direct C–N source injection, atom-imbalanced sources,
HCN-assay-as-input declarations and minimum-source-energy violations are
rejected. The endpoints are `/api/prebiotic/hcn/reduced-mechanistic` and
`/api/prebiotic/hcn/atom-ceiling`; the Origins panel exposes both.

### Optional assay-blind numerical audit (2026-09-29)

The Origins panel now offers **Audit solver and controls**, backed by
`/api/prebiotic/hcn/reduced-mechanistic/audit`. It solves the exact supplied
scenario twice, with maximum steps no larger than duration/16 and duration/32
(or tighter if the user supplied a limit). It compares final gas HCN, aqueous
HCN+CN−, total HCN-equivalent, and each integrated source extent against a
declared 1% relative tolerance plus a `1e-12 mol` absolute reporting floor.
It then reruns without effective sources to expose any seeded or thermal
background, and separately removes all carbon-bearing or nitrogen-bearing
initial gas species. The latter controls must leave a positive gas inventory;
otherwise the audit reports **incomplete**, never a pass. Every run uses the
same declared conditions and still enforces the atom ledger and source-energy
lower-bound check. This is a numerical/conservation audit, not a source-rate
measurement or chemical-model validation. In particular, a near-zero
source-off background does not prove that a real spark produces the supplied
fragments, and no HCN assay is read by the audit.

The synthetic CH4/N2 fixture passes these controls and the step comparison;
an HCN-only seeded fixture correctly returns incomplete because neither
atom-removal control retains a gas phase. The UI labels any result
**conditional**, whether the audit passes or not.

## What remains before an experimental prediction

1. Independently establish effective, **state-resolved** fragment source
   hazards for the actual wet discharge. The total N2 dissociation cross
   section cannot supply a ground-state `N2 -> 2 N` hazard. A measured stable
   co-product may constrain a source only with an identifiable mechanism;
   fitting it to HCN and then validating on the same HCN value is circular.
2. Measure or independently bound the target run's `kLa`, water pH and
   temperature history; extend this reduced model to moving pH, major
   co-solutes and aqueous losses if their timescales are material.
3. Qualify the gas mechanism's low-temperature CH+N2, HCNN/NCN/HNCN,
   H+HCN/H2CN and other key rates. GRI/Pavan sharing is not independent
   validation. A few-key-channel reduced mechanism may ultimately be more
   appropriate than an unqualified full GRI network.
4. Close the delivered-energy/thermal balance and compare a frozen run
   against a held-out HCN+CN− assay with matching feed, sampled phase and
   duration. Until then, neither the conditional result nor its ratio to
   the atom ceiling is a real yield.

The current audit does not reduce these scientific blockers. Its purpose is
to prevent a numerically unstable or conservation-violating conditional run
from being mistaken for a useful scenario while independent physical inputs
are sought.

See [HCN_FORMATION_GATE.md](HCN_FORMATION_GATE.md) for the current wet-run
blockers and [BARTH_LONG_HORIZON_GATE.md](BARTH_LONG_HORIZON_GATE.md) for the
corrected Barth-feed numerical record. The reduced route's GRI mechanism is
packaged by [Cantera](https://cantera.org/stable/examples/input/gri30.html);
the 298 K HCN solubility basis in the UI template is the
[NIST-reviewed entry](https://webbook.nist.gov/cgi/cbook.cgi?ID=C74908&Mask=10).
