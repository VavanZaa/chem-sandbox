# Barth feed-only long-horizon numerical gate (2026-09-27)

**Reading guide:** The current source-audited status is in the sections
above `Boundary and method` and in `HCN_FORMATION_GATE.md`. Older tables
retained below are an audit trail and may predate the correction that stopped
assigning both N2 electron-impact fragments to ground-state N atoms. In
particular, the old low-contact pH-5-to-7 run is superseded by the
[current low-contact check](#contact-rate-and-kinetic-domain-sensitivity).

## CH + N2 association partial reference audit (2026-09-27)

The feed-only long-horizon report now includes the read-only
`unmodified_gri_ch_n2_initial_feed_reference_audit`. The assay is never
opened. [Brownsword et al. (1996)](https://pubs.rsc.org/en/content/articlelanding/1996/ft/ft9969201087)
measured predominantly associative CH + N2 loss in Ar at 202-584 K and
reported an extrapolated low-pressure coefficient of
`1.7e-31 (T/298)^-2.3 cm6 molecule-2 s-1`. At the declared 298.15 K
feed temperature, packaged GRI's Ar low-pressure limit is
`1.55885e-31`, versus `1.69803e-31` from that fit (ratio `0.918`).
The [NIST kinetics record](https://kinetics.nist.gov/kinetics/Detail?id=1996BRO%2FHER1087-1094%3A2)
lists the source Ar pressure range as 0.00533-0.53 bar, below the
declared 1-bar wet-flask pressure.
This comparison does **not** establish the 1-bar wet-mixture effective
rate, HCNN identity or fate, the competing NCN branch, or product yield.
The [primary paper's abstract](https://pubs.rsc.org/en/content/articlelanding/1996/ft/ft9969201087)
also explicitly cautions that its global-fit high-pressure constants are
merely fitting parameters, appreciably smaller than separate relaxation
measurements. We therefore cannot promote its falloff fit into an
independent wet-flask rate law.
The actual model's conditional wet-feed effective rate is shown separately,
not treated as a measured quantity. The source mechanism and assay remain
unchanged.
At the assumed feed, GRI's own Troe falloff calculation is only `0.5616`
times its wet-mixture low-pressure projection. Thus even within the model,
the low-pressure coefficient alone does not determine the effective rate;
the unvalidated high-pressure/falloff and collider assumptions matter.

An independent check of the other asymptote is now attached to the same
assay-free audit. [Fulle and Hippler (1996)](https://doi.org/10.1063/1.472383)
measured CH+N2 association across 1–150 bar **helium** and 200–715 K,
then extracted `k∞=(4.1±0.8)e-11(T/300)^-0.15 cm3 molecule-1 s-1` from
falloff curves at 200, 250, 300, 400 and 500 K. At the declared 298.15 K,
the installed GRI high-pressure limit is `1.20997e-11`, versus their fitted
`4.10381e-11 cm3 molecule-1 s-1`, a ratio of `0.29484`. The displayed
`±0.8e-11` is the paper's coefficient width, **not** an HCN-yield uncertainty.
The installed mechanism's **actual wet-feed effective** rate is only
`2.43836e-12 cm3 molecule-1 s-1`, or `0.20152` of its own high-pressure
limit, confirming that the 1-bar feed is in falloff in this model.
This is evidence of a high-pressure-limit discrepancy in the mechanism,
but neither asymptote measures the intermediate-pressure effective rate in
the 1-bar wet N2/CO2/H2/H2O/O2 feed. It does not identify the stabilized
adduct isomer, specify HCNN/HCN branching, or justify replacing the model
rate without a source-supported mixture-dependent falloff law.

An opt-in **diagnostic only** charged-shadow multiplier now tests the
installed reversible `CH + N2 (+M) <=> HCNN (+M)` record without changing
the primary GRI mechanism or claiming a new wet-mixture rate. Halving that
record is arbitrary, **not** the ratio of a measured wet rate to GRI. At
one hour/eight substeps, the signed-net CH+N2 contribution falls from
`1.57895e-8` to `1.05767e-8 mol`; the HCN endpoints in that short run
are both below the `1e-12 mol` reporting floor and are not compared as
yields. At 925 minutes the charged aqueous HCN+CN- amount is
`2.7808761e-12 mol`
at 192 substeps and `2.7808041e-12 mol` at 384, versus unmodified
`3.9389551e-12 mol` at 192. The two stressed resolutions differ by
`2.59e-5` relatively. The 192-step integrated HCN trace has a
`7.85e-21 mol` quadrature closure residual; both material ledgers close,
there are zero failed charged substeps, and the primary result is unchanged.
Both one-hour and 925-minute/192-substep carbon-/nitrogen-removal stress
controls return zero aqueous HCN+CN- in both branches; at full duration
both material ledgers close and there are zero failed charged substeps.
The HCN assay remains unread. This exposes structural sensitivity to an
unvalidated association/falloff and adduct identity; the roughly 29%
modeled HCN change is **not** a measured-rate uncertainty, corrected
yield, or independent validation. CVODES minimum-step warnings remain.

## Competing CH + water-vapour loss audit (2026-09-27)

The same assay-free feed-only report now includes a **read-only** comparison
for GRI reaction 126, `CH + H2O <=> CH2O + H`, with the installed
Cantera version and mechanism SHA-256 recorded. At the declared 298.15 K
initial feed, packaged GRI gives `3.39081e-11 cm3 molecule-1 s-1`, versus
`1.55889e-11` from [Blitz et al.'s primary CH-removal measurement](https://doi.org/10.1021/jp9903716)
(`1.56e-11(T/298)^-1.42`), a ratio of `2.17515`. The paper reports
`±0.25e-11` on its 298 K coefficient; its room-temperature pressure test
covered 20–200 Torr, **not** the assumed 1-bar wet feed. This comparison
is neither a wet-flask rate calibration nor an HCN-yield error bar.

With the *initial* measured gas fractions (`xH2O=0.016`, N2 balance
`xN2=0.917114`), the installed GRI CH+H2O loss hazard is `0.24261` of
its CH+N2 association hazard. This is an initial-state, same-density
hazard ratio, **not** a reaction branching probability. In the separate
assay-free one-hour/eight-substep charged trace, the signed-net CH loss
contributions are `4.01026e-9 mol` for CH+H2O and `1.57895e-8 mol` for
CH+N2; the evolving trace is broadly consistent with water being a
non-negligible competitor, but neither quantity is a microscopic forward
event count or an experimental yield.

[Bergeat et al.'s independent product study](https://doi.org/10.1016/j.cplett.2009.07.080)
measured 100% H-atom production at 330 K, but its proposed exits include
both H+CH2O and routes leading to H2+H+CO. It does **not** experimentally
assign all CH+H2O loss to the GRI `CH2O+H` product pair. No product split
or revised rate is inserted into the propagating solver. The next
mechanism step needs a validated pressure-/temperature-applicable CH+H2O
rate and product split alongside the CH+N2 wet falloff and HCNN fate.

An **opt-in, assay-free rate stress**, not a replacement mechanism, scaled
only the charged-shadow CH+H2O reaction by `1/2.17515` (the ratio of the
Blitz 298 K fit to packaged GRI). It deliberately keeps GRI's unverified
`CH2O+H` product assignment, extrapolates beyond the paper's 20–200 Torr
room-temperature pressure test to the assumed 1-bar feed, and leaves the
primary mechanism untouched. In the one-hour/eight-substep diagnostic,
signed-net CH+H2O consumption fell from `4.01026e-9` to `2.26692e-9 mol`
while signed-net CH+N2 consumption rose from `1.57895e-8` to
`1.94145e-8 mol`. Re-running the baseline in the same process restored its
HCN result exactly; carbon- and nitrogen-removal negative controls remain
zero. Separate 925-minute/192-substep rate-stressed carbon-removal and
nitrogen-removal controls each also return **zero** aqueous HCN+CN- in
both branches, with both material ledgers closed, zero failed charged
substeps, and the experimental HCN assay unopened. One-hour HCN amounts
lie below the `1e-12 mol` numerical reporting floor and are not compared
as yields.

At the declared 925-minute horizon, this stress gave charged aqueous
HCN+CN- `4.6957446e-12 mol` with 192 substeps and `4.6956284e-12 mol`
with 384 substeps, versus the unmodified 192-substep
`3.9389551e-12 mol`. The two stressed resolutions differ by about
`2.5e-5` relatively. The stressed 192-step HCN/CH traces reconcile to
their inventories, both material ledgers close, there are zero failed
steps, the primary result is unchanged, and no experimental HCN assay was
opened. CVODES minimum-step warnings persist. The roughly 19% change is
**model sensitivity to an extrapolated single-channel rate assumption**,
not a calibrated correction, experimental uncertainty, or evidence that
the real wet-flask HCN yield is closer to either value.

The full-duration charged HCNN signed-net-fate audit places `85.040%` of
modeled HCNN consumption in `H + HCNN <=> CH2 + N2`, `13.147%` in
`HCNN + O <=> CO + H + N2`, and `1.195%` in
`HCNN + O <=> HCN + NO` (total net consumption `3.6597240e-7 mol`).
These fractions are neither forward event probabilities nor measured
branching ratios. Halving only the dominant H+HCNN rate as an arbitrary
diagnostic stress raises charged aqueous HCN+CN- from `3.9389551e-12`
to `6.2140679e-12 mol` at 925 minutes; primary HCN is unchanged,
both material ledgers close and the HCN assay stays unread. This identifies
the H+HCNN sink as a further low-temperature kinetic validation target,
not a correction to the predicted amount or an error bar.
At 384 instead of 192 substeps, the three signed-net HCNN fractions are
`0.85040317`, `0.13147515`, and `0.01195229`, differing by at most
`2.9e-7` in absolute fraction from the 192-step values. Charged aqueous
HCN+CN- changes by about `2.52e-5` relatively; both integrated HCN/HCNN
traces reconcile and both material ledgers close. CVODES minimum-step
warnings and the kinetic/thermal validation gaps remain.

## One-hour charged carbon-route audit (2026-09-27)

The opt-in charged quadrature now traces precursor gas species as well as
HCN/HCNN, and separately aggregates accepted external electron-source
extents. For the assay-free 60-minute, eight-substep feed scenario, the
following **signed net reaction contributions** are present in the
installed Pavan/GRI mechanism (moles, not a single-molecule trajectory):

| Modeled link | Integrated signed net amount |
| --- | ---: |
| CO2+ + e -> CO + O | `5.67110e-6` |
| CO + H + M -> HCO + M | `8.20201e-6` |
| HCO + H + M -> CH2O + M | `9.03969e-7` |
| CH2O + H + M -> CH2OH/CH3O + M | `8.67393e-8` / `7.12412e-8` |
| CH2OH/CH3O + H -> CH2(S) + H2O | `1.45082e-8` / `1.03649e-8` |
| CH2(S) + N2 -> CH2 + N2 | `1.31931e-8` |
| CH2 + H -> CH + H2 | `2.14554e-8` |
| CH + N2 + M -> HCNN + M | `1.57895e-8` |
| HCNN + O -> HCN + NO | `1.30445e-10` |

The charged external source trace separately records CO2 ionization and
H2O -> H + OH collision extents; it does not establish that one source
event caused any particular downstream HCN molecule. The network contains
strong return cycles (especially HCNN + H -> CH2 + N2) and competing sinks,
so multiplying these numbers or calling them branching probabilities would
be wrong. CO2 is the supplied carbon-bearing gas, but the quantitative
route and yield still depend on unvalidated source, gas-rate, transfer,
spark and thermal assumptions. No HCN assay or reaction recipe was used to
construct this trace.

Tracing CH2O exposed a diagnostic closure defect: its gas+free-liquid
inventory excluded the coupled methanediol pool. The quadrature audit now
adds the methanediol delta **only when separating gas reaction progress
from aqueous hydration**. In the one-hour case, CH2O's gas-chemical net
is `6.07679118e-9 mol` from the corrected inventory decomposition and
`6.07679124e-9 mol` from reaction quadrature, residual `5.97e-17 mol`.
The propagating chemistry and HCN amount are unchanged; this repairs the
diagnostic ledger, not the missing aqueous/thermal mechanism.

## Current-state correction — 2026-09-27

The charged Pavan shadow is **not independent HCN-kinetic validation** of the
primary GRI branch. An exact-record audit of the Cantera-packaged files finds
325/325 GRI reactions embedded in Pavan's 440 reactions, including all 30
HCN-containing reactions. The application exposes the audit and both
mechanism hashes in the charged report. The earlier opt-in charged HCN
endpoint trace records
net reaction rates only at substep endpoints, not time-integrated
contributions; it cannot be used to apportion the final HCN amount. In the
925-minute/192-step
diagnostic, the final endpoint's largest listed HCN consumption rate is
H + HCN (+M) <=> H2CN (+M), but this single endpoint cannot establish an
integrated bottleneck or a corrected yield. The held-out assay was not read.

An opt-in integrated charged trace now supersedes that *endpoint-only*
ranking for the same 925-minute/192-step scenario: signed HCN formation
`5.3080601e-9 mol`, signed consumption `5.3038444e-9 mol`, observed net
gas-chemical change `4.2156767e-12 mol`, and quadrature closure residual
`1.10e-20 mol`. HCNN+O is the largest integrated HCN source; H+HCN is the
largest sink. The charged aqueous endpoint is still `3.9389551e-12 mol`.
This strong cancellation and GRI/Pavan ancestry make the amount sensitive to
unverified shared low-temperature rates. A shorter precursor trace finds
CH+N2 as the modeled HCNN source, while the mechanism omits NCN. See the
[formation gate](HCN_FORMATION_GATE.md#current-scientific-status--2026-09-27)
for the temperature-limited primary literature and why no NCN substitution
was made. This remains neither a calibrated uncertainty nor a passing yield
validation; the HCN assay was not loaded.
Assay-free full-duration charged rate stresses give `7.3136207e-12 mol`
aqueous HCN+CN- when the H+HCN sink is halved and `7.5322252e-12 mol`
when the HCNN+O-to-HCN source is doubled, versus the baseline
`3.9389551e-12 mol`. These are arbitrary mathematical factors, not
experimental uncertainty bounds; both ledgers still close and the primary
result is unchanged.

The total N2 electron-impact dissociation cross section no longer injects
two ground-state N atoms by default: [Cosby's product-state measurement](https://doi.org/10.1063/1.464385)
indicates N(2D)+N(4S) is the primary channel, while the imported
[Tabata total N+N table](https://doi.org/10.1016/j.adt.2006.02.002)
does not provide state-selective branching. A rerun of the feed-only 60-minute,
48-step scenario gives primary aqueous HCN+CN- `1.16285e-46 mol` and charged
aqueous HCN+CN- `9.17022e-14 mol`, both below the `1e-12 mol` numerical
reporting floor. Accepted primary N2 neutral-fragment extent is zero;
both material ledgers close, the held-out assay remains unopened, and the
thermal balance remains incomplete. This is **not** evidence that a real
spark cannot make HCN. It means the current state-resolved source/chemistry
chain is incomplete. All older positive quantities in this document were
computed with the superseded all-ground-state N2 assignment and must not be
used as current quantitative predictions, validation, or bounds. The original
tables remain as an audit trail of earlier numerical/kinetic work.

For the reported 925-minute duration, a 192-step rerun gives primary aqueous
HCN+CN- `4.18575e-27 mol` and charged aqueous HCN+CN- `3.93895508e-12 mol`;
the 384-step rerun gives `3.25350e-26 mol` and `3.93885589e-12 mol`,
respectively. The primary results are both many orders below the `1e-12 mol`
reporting floor, so their relative difference has no physical meaning. The
charged endpoint differs by about `2.52e-5` between these resolutions and
has zero reported failed substeps, but CVODE emitted repeated minimum-step
warnings; numerical conditioning is not fully cleared. Both material ledgers
close and the primary N2 neutral-fragment extent is zero. The charged branch
contains other excited/ion pathways from its published mechanism, but is
correlated with the same collision input and **is not a separately validated
yield**. The held-out HCN assay remains unopened and no rate was fit to it.

The initial-gas, conditional [N(2D) first-collision audit](HCN_FORMATION_GATE.md#conditional-n2d-fate-in-the-wet-feed-2026-09-27)
finds H2O responsible for 91.6% of evaluated thermalized-N(2D) removal
hazard at the declared 298.15 K, 1 bar geometry. Its source-evaluated
rate-factor extremes span 79.3-96.9%, not a confidence interval. This is
not a product yield: N2 electron-dissociation state branching and the
N(2D)+H2O product split are still unresolved.
Across a 60-minute 48-step rerun, the first N(2D) fate snapshot is covered
but 47 subsequent snapshots are incomplete because of evolving H/OH/HO2/
H2O2 bath species. The evaluated N(2D)+O removal rate is now included,
so O is no longer an unreviewed partner in this audit. See the formation
gate for the remaining unreviewed-inventory and equal-hazard-rate audit;
neither provides a product branch or HCN yield.
The 925-minute 192-step rerun likewise has 1 covered initial snapshot and
191 unavailable later N(2D) bath snapshots. Its minimum calculated
N(2D)+H2O2 equal-hazard rate is `3.65007e-9 cm3 molecule^-1 s^-1`;
this is not measured kinetics or a physical bound. No N(2D) products
were added to the propagating mechanisms.
The EEDF-folded total N2 outcome has an additional applicability gap:
`55.975%` of its 925-minute uncoupled event opportunities arise from
electron energies below the 18.5 eV lower edge of Cosby's fragment-state
measurement. The corresponding 60-minute fraction is `56.008%`; the
150 eV extended-grid initial-feed check changes the total rate by under
`0.23%` relative to the 60 eV grid. These are energy-support fractions,
not state branching or product yields. The branch for the majority of
computed N2 dissociation events remains experimentally unassigned.
The QST/Tabata official analytic fit, evaluated on the same EEDF as a
**read-only** source-form sensitivity, gives `52.449%` below that window
over the 925-minute run, versus `55.975%` from the imported tabular
interpolation. The corresponding `0.00454248` versus `0.00493504 mol`
figures are uncoupled collision opportunities, not N-atom amounts.
Neither representation resolves product electronic states or changes
the propagating chemistry.
The new conditional energy-access audit puts `11.088%` of the tabular
925-minute uncoupled N2 outcome opportunities below the nominal 14.526 eV
N(2D)+N(2D) opening, where a ground-state N2 target could only dissociate
to neutral atoms with at least one N(4S). The read-only analytic-fit share
is `12.256%`. This does not supply the missing N(2D)/N(4S) branching or
justify adding either atom to the solver: near-threshold total cross
sections are not fragment-state measurements, and the spark field is
unmeasured. See the formation gate for source and grid sensitivity.
An assay-free 50–200 Td scan of this feed changes the calculated *total*
N2 outcome rate by about `2.08e4`-fold and the below-measurement-window
share from 99.815% to 56.014%. These are model sensitivities to the
unmeasured reduced field, not a statistical or physical range for HCN.

## Boundary and method

The [published Barth feed table](https://pure.st-andrews.ac.uk/ws/portalfiles/portal/296987606/SparkExperiments.csv)
and [dataset README](https://pure.st-andrews.ac.uk/ws/portalfiles/portal/300064695/Readme.txt)
identify run `220725001` as a 925-minute wet, closed, stirred spark with
1 L flask and 50 mL water. The local feed-only card
`backend/data/benchmarks/barth_2024_feed_only_220725001_v1.json`
contains the initial H2/CO2/O2/H2O fractions with N2 balance but **no
HCN assay value**. `backend/scripts/probe_barth_long_horizon.py` reuses
that card and the existing coupled gas/water apparatus; it divides the
reported duration into 16 equal cycles and varies 12/24/48 integration
substeps per cycle. It never reads, fits, or scores the held-out assay.

The 0.95 L headspace, 1 bar, 298.15 K, active spark volume, pulse width,
electrical coupling, reduced field and common liquid-side `kLa` are the
same declared assumptions as the one-hour feed-only probe, **not**
measurements of run `220725001`. The primary GRI and correlated charged
Pavan branches use target-inventory-limited electron sources. The aqueous
quantity below is final HCN+CN− summed over the model's aqueous
compartments; it is not an assay concentration or inferred experimental
yield. The gas mechanisms already contain HCN pathways, so this is not
blind product discovery.
The scenario's `49 kV × 1 mA × 0.1` duty × `0.3` coupling arithmetic
allocates `88.2 J/min`, or `81,585 J` over 925 minutes, as a **declared
electrical estimate**. It is not measured delivered spark energy and is
not a closed thermal balance. The runner reports that basis and the
source-threshold budget separately so a numerical endpoint cannot be
silently labeled a measured yield per joule.

An assay-free [IAPWS-95 water-property screen](https://www.iapws.org/relguide/IAPWS-95.html)
is reproducible with `backend/scripts/audit_barth_thermal_capacity.py`.
At the declared 298.15 K, 1 bar and 50 mL initial liquid, Cantera's
water equation of state gives `0.049856 kg` water, an initial heat
capacity of `208.52 J K^-1`, and `15,586 J` of liquid sensible enthalpy
to reach saturation at `372.81 K` and the same pressure. The scenario's
**hypothetical**, not measured, `81,585 J` electrical-deposit estimate
is `5.23` times that isolated liquid sensible capacity. This does not
show that the experiment boiled: the real energy deposition, gas/glass
storage, chemical enthalpy, phase change, and environmental heat transfer
are unknown. It does show that treating the entire declared input as
heat in an isolated 50 mL liquid is incompatible with the fixed 298.15 K
state. The screen infers neither actual water temperature nor heat loss;
the full thermal ledger remains open and a measured power/temperature
history is needed before any yield-per-energy claim.
The long-horizon runner now embeds the same conditional screen at its
*selected* duration and checks that its declared input agrees with the
runner's electrical estimate; a shorter numerical probe is not silently
assigned the full 925-minute thermal exposure.

## Reported-duration numerical stress

At hypothetical `kLa=1e-4 s^-1`, using the explicitly flagged
Winkelman pH-4-to-5 extrapolation diagnostic:

| Integration steps | Primary aqueous HCN+CN− (mol) | Charged aqueous HCN+CN− (mol) |
| ---: | ---: | ---: |
| 192 | `7.55623512109e-11` | `1.28772216959e-10` |
| 384 | `7.55605214661e-11` | `1.28752412929e-10` |
| 768 | `7.55596080766e-11` | `1.28742519554e-10` |

The 384-to-768 relative changes are `1.209e-5` primary and `7.685e-5`
charged. Both material/charge ledgers close at all three resolutions and
the charged branch has zero failed substeps. Minimum modeled ideal pH is
about `4.714`, so this baseline **uses an extrapolated hydration rate**;
it is not source-bounded kinetic validation. The 768-step amount is not
promoted to an experimental yield. The primary and charged endpoints
also differ by roughly 1.70-fold, and the branches share source and
mechanism assumptions; their agreement or disagreement is not an
independent probability of products.

At 384 steps, full-duration carbon removal gives `8.79e-43 mol` primary
aqueous cyanide and zero charged cyanide. Full-duration nitrogen removal
gives zero in both branches. All are below the project's `1e-12 mol`
**numerical reporting floor** (not an instrumental detection limit), with
both ledgers closed and zero charged failures. The two-hour pilot also
completed with closed ledgers and minimum ideal pH `4.733`.

## Contact-rate and kinetic-domain sensitivity

The contact coefficient is **not** inferred from the HCN result. The next
table is an **archived, superseded** 384-step sweep from the former
all-ground-state N2-fragment assignment. Its pH/amounts and especially its
"strict source-bounded" low-contact row are **not current-code results**:

| Hypothetical common kLa (s⁻¹) | Hydration mode | Minimum ideal pH | Primary aqueous HCN+CN− (mol) | Charged aqueous HCN+CN− (mol) |
| ---: | --- | ---: | ---: | ---: |
| `1e-5` | strict source-bounded | `5.1857` | `3.29563945074e-11` | `5.52255050014e-11` |
| `1e-4` | flagged pH extrapolation | `4.7136` | `7.55605214661e-11` | `1.28752412929e-10` |
| `1e-3` | flagged pH extrapolation | `4.5531` | `7.42996905681e-11` | `1.25249742912e-10` |

In that superseded model, each contact-rate case also passed a 192-step run. Its 192-to-384
relative endpoint changes stay below `2.0e-4` in both branches. The
primary amount spans about a factor of 2.29 across these assumptions;
the response is not monotonic because contact affects both precursor
loss/formation and HCN collection. The old low-contact case stayed within
Winkelman's measured pH 5-7 domain for all 925 minutes, but its `kLa`
is just as hypothetical as the others. These values are controlled
model sensitivities, **not** a confidence interval, a fitted rate, or a
comparison to the assay.

The **current** corrected code has been rechecked for `kLa=1e-5 s^-1`.
At 925 minutes with the explicitly flagged pH-4-to-5 extrapolation mode,
the 192/384-step charged aqueous HCN+CN- endpoints are
`1.5028622e-12` / `1.5028469e-12 mol`; both primary endpoints are below
the `1e-12 mol` numerical floor. Minimum ideal pH is `4.71706` at both
resolutions, so the current path **does** leave Winkelman's pH-5-to-7
measured domain. Both ledgers close, there are zero failed charged steps,
and the assay remains unopened. The 192/384 charged difference is about
`1.02e-5` relatively, a numerical check only. Running the strict
source-bounded hydration mode at 192 steps correctly **rejects** the
current baseline after cycle 4/substep 8 when pH reaches `4.9965`.
Therefore the archived source-bounded row must not be cited as evidence
that the present model is source-bounded, nor may the flagged current
value be compared with the experimental HCN assay as a validated yield.

## Still required for quantitative validation

1. Independently measured or defensibly bounded wet-spark waveform,
   active filament volume/exposure, deposited energy, and gas/water
   contact for a cyanide-assayed run. Nominal BD-50E settings and
   stirring are insufficient.
2. A source-supported hydration model for any pH below 5 reached by the
   real contact regime, plus a more complete time-dependent aqueous
   network and phase-sampling model. Do not choose `kLa` from the assay.
3. A closed thermal-energy balance (reaction enthalpies, gas/water heat,
   evaporation/condensation, wall and radiation losses). The current
   electrical/collision budget and partial hydration heat do **not**
   establish actual deposited energy or yield per joule.
4. Independent sensitivity to spark assumptions at the full duration.
   The former SciPy finite-difference Jacobian overflow warning has been
   resolved as a numerical issue (see below), but three-point timestep
   stability alone does not validate the mechanism or transport parameters.

Only after these inputs are frozen independently should the held-out
aqueous HCN+CN− assay be opened for a like-for-like, order-of-magnitude
validation test. No known product amount has entered the solver in this
gate.

## Passive-coordinate Jacobian repair (2026-09-27)

A warning-as-error reproduction localized the repeated overflow to SciPy's
adaptive finite-difference multiplier in the charged Pavan BDF solve; the
longer run exposed the same issue in the primary GRI BDF solve. Captured
liquid and integrated source extents are write-only ledger coordinates,
and absent-element gas species are identically inactive. Their Jacobian
columns are exactly zero. Both branches now provide a bounded finite-
difference Jacobian for the remaining active columns, keeping the ODE,
solver tolerances, and conservation tests unchanged.

The one-hour 48-step run and the full 925-minute 384-step run complete with
`RuntimeWarning` treated as an error. The full-run primary amount is now
`7.55605211630e-11 mol` versus `7.55605214661e-11 mol` before the repair;
charged is `1.28752409501e-10 mol` versus `1.28752412929e-10 mol`.
Both material/charge ledgers close and no charged substep fails. These
differences are numerical, not a change in the chemistry model or evidence
for an experimental yield. The spark/contact, aqueous-kinetics and thermal
limitations above remain.

At the same 384-step resolution with warnings treated as errors, separate
full-duration carbon- and nitrogen-removed controls each report zero aqueous
HCN+CN− in both branches, closed ledgers, and zero failed charged substeps.
The new one-hour warning regression and the entire backend test suite pass
with `RuntimeWarning` treated as an error. A separate Starlette deprecation
warning from the test client's installed dependency remains unrelated to the
chemical solvers.

## Narrow independent-discharge evidence check (2026-09-27)

The [Barth cyanide-assayed wet-flask study](https://doi.org/10.1051/0004-6361/202347286)
reports nominal generator settings but no simultaneous measured wet-flask
waveform, active filament occupancy, or gas-water contact coefficient for
its cyanide-assayed run. A separate [Miller–Urey-style generator study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10672138/)
reports generator voltage-time characterization and GC/MS composition,
but not a quantitative HCN assay in the characterized runs. A separate
[gas–water lightning-electrochemistry experiment](https://pmc.ncbi.nlm.nih.gov/articles/PMC11317556/)
normalizes measured nitrogen fixation to recorded current profiles, but
measures nitrate, nitrite, and ammonium rather than HCN in a different
reactor/feed. None supplies the missing Barth exposure or `kLa` by
transfer; this was a bounded search, not proof no such dataset exists.

An [NMR formaldehyde study](https://doi.org/10.1021/jp513020y) spans
pH 1.8–7.8, but its accessible quantitative pH-4 rate concerns D2O,
whereas the model's water is H2O. The [Winkelman H2O correlation](https://doi.org/10.1016/S0009-2509(02)00358-5)
remains source-bounded to pH 5–7. The pH-4-to-5 option must retain its
explicit extrapolation label; mixing the D2O rate into the H2O model
would not be a measured-domain repair. The highest-value missing evidence
is a cyanide-assayed wet run with recorded discharge current/voltage,
filament geometry or optical occupancy, and independent gas-water contact
characterization under the same feed and duration.

There is a narrower **unimplemented** pH-4-to-5 lead: the primary
[Winkelman et al. (2002) paper](https://pure.rug.nl/ws/files/3619450/2002ChemEngSciWinkelman.pdf)
reports that [Schecker and Schulz (1969)](https://doi.org/10.1524/zpch.1969.65.1_4.221)
measured H2O formaldehyde hydration at pH 4–7 and 298–333 K by a
temperature-jump/UV method, fitting `k_h=7800 exp(-1913/T) s^-1`.
At 298.15 K this gives `12.7515 s^-1`, versus `10.7885 s^-1` from
Winkelman's own pH-5-to-7 H2O fit, an 18.2% difference relative to
Winkelman. This is a cross-study consistency lead, **not** yet an adopted
low-pH rate: the 1969 primary methods/uncertainties have not been checked
directly, and its reactor/buffer conditions have not been mapped to the
Barth liquid. Keep the existing extrapolation flag until that source and
the model's equilibrium/reverse-rate pairing are audited; do not change
the propagating HCN result on the strength of a cited formula alone.

A further bounded candidate check does not supply a substitute assay. The
[Scherer et al. wet Miller-type study](https://doi.org/10.1007/s11084-016-9528-8)
reports an approximately 20 W gas discharge and spectroscopic bands from
cyanide-containing material, but not a quantified aqueous HCN+CN- amount
matched to its electrical input. The
[Meng et al. microdroplet-discharge study](https://doi.org/10.1126/sciadv.adt8979)
reports HCN detection by mass spectrometry in a sprayed-droplet/gas system;
that geometry, collection basis, and observable are not the closed wet-flask
cyanide assay here. Neither observation may be inserted as a Barth spark
exposure, mass-transfer coefficient, or order-of-magnitude amount benchmark.

## Full-duration excluded-source audit (2026-09-27)

The same feed-only 925-minute, 384-step baseline now reports cumulative
accepted neutral electron-source extents and excluded O(1D) opportunities.
The primary branch accepts `0.00493621 mol` N2, `0.000305215 mol` H2O,
`5.21137e-8 mol` O2, and **zero direct CO2** fragment extent. Every one
of the 384 steps records a CO2 -> CO(X)+O(1D) channel as excluded; summed
uncoupled and energy-scaled opportunities are both `2.34473e-4 mol`
under this hypothetical discharge. They are possible electron events
before excited-state fate, target competition, and coupled chemistry, **not**
formed CO, oxygen, HCN, or measured product moles. Summing them cannot
repair the missing chemistry or validate an HCN amount.

The near-298 K O(1D) audit resolves the first substep but returns
`unavailable` for the other 383 as unreviewed bath partners appear. The
maximum simultaneous unreviewed-bath mole fraction reported is
`5.22e-5`; this small material fraction alone does not bound its O(1D)
reaction hazard without partner-specific rate evidence. Among individual
maximum amounts are HNCO `4.42e-8 mol`, N2O `4.39e-8 mol`, H2O2
`3.27e-8 mol`, and CH4 `1.56e-8 mol` (maxima need not occur at the same
time). The assay remains unloaded. This audit narrows the next chemistry
gate to an atom/charge/energy-conserving CO2 dissociation and O(1D)/CO(a3P)
fate model under this wet gas mixture, not another increase in search runs.
The source-reporting and O(1D) focused tests pass (16 tests) with runtime
warnings treated as errors.

To rank the unknown-bath issue without claiming a rate bound, the read-only
O(1D) audit now assigns *every* unreviewed species either a hypothetical
`1e-9` or `1e-8 cm3 molecule-1 s-1` total removal coefficient while
retaining its `unavailable` status. Against the reviewed hazards in each
step, the largest resulting hypothetical unknown shares across 925 minutes
are `0.001339` and `0.013234`, respectively. These are stress scenarios,
**not** measured rates, physical upper limits, confidence intervals, or
resolved product branches. They suggest where to prioritize rate research;
they do not authorize injecting the excluded CO2 fragments. The O(1D)
and long-horizon focused tests pass (17 tests), and the HCN endpoint is
unchanged because this audit is read-only.

An isolated `project_thermalized_co2_o1d_events` kernel now uses the
reviewed first-collision branches only when bath coverage is complete and
the event batch changes each consumed reactant by at most 1%. It maps
CO2 -> CO(X)+O(1D) to signed CO, ground-O, OH, H, H2O and H2 changes
with a closed C/H/O and charge ledger. [NIST's O I level table](https://physics.nist.gov/PhysRefData/Handbook/Tables/oxygentable5_a.htm)
places O(1D2) at `15867.862 cm-1` above ground, approximately `1.967 eV`
per atom. The kernel reports that initial electronic energy separately;
it **does not** call it thermal heat, claim a complete energy balance, or
inject any product into the apparatus. Nascent electron-fragment energy
and the 383 later unresolved baths still prevent full-run coupling. The
focused O(1D)/long-horizon suite passes (18 tests) with runtime warnings
treated as errors.

The imported private CO2 -> CO(X)+O(1D) cross section labels a `7.0 eV`
threshold and has its first nonzero tabulated point at `8.0 eV` incident
electron energy. From reviewed [CO2](https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=2781&Units=SI),
[CO](https://webbook.nist.gov/cgi/cbook.cgi?ID=C630080&Mask=FFFFFF), and
[ground O](https://webbook.nist.gov/cgi/cbook.cgi?ID=C17778802&Mask=29)
formation enthalpies plus the NIST O(1D) level, the standard 298.15 K
CO2 -> CO+O(1D) enthalpy requirement is about `7.48 eV` per event.
The `8.0 eV` point is an *incident-energy grid point*, not proof that
`8.0 eV` was deposited; the `7.0 eV` label alone would underbudget the
thermochemical change. The isolated projector reports the `7.48 eV`
requirement and explicitly says no electron-collision allocation has yet
been checked. Any production coupling must reconcile source-energy loss,
nascent fragment energy, chemical enthalpy, and heat before these events
can enter a claimed closed thermal budget.

After this diagnostic and isolated-kernel addition, the entire backend
test suite again passes with `RuntimeWarning` treated as an error. The
remaining Starlette test-client deprecation is unrelated to the chemical
calculations.

## Accepted-fragment thermochemical floor (2026-09-27)

The feed-only `audit_source_thermochemistry.py` now computes the standard
reaction enthalpy of every **accepted ground-state neutral electron fragment**
from GRI-Mech at the 298.15 K scenario temperature. It rejects an ordinary
accepted channel if its assigned per-event source-energy cost falls below
that endothermic floor by more than 0.02 eV. In the two-substep feed check,
H2O -> H+OH requires `5.174 eV` versus its `7.5 eV` assigned cost, and
O2 -> O+O requires `5.165 eV` versus `6.1 eV`. No ordinary accepted channel
has a deficit. This check uses no HCN assay and does not modify the solver.

The N2 -> N+N outcome is different: its direct energy charge is zero because
the parent electronic-excitation manifold carries the modeled electron loss.
Its GRI standard reaction enthalpy is `9.798 eV` per event. A separate,
weaker necessary check compares accepted N2 outcomes in each apparatus
substep to the energy-scaled opportunities for N2 electronic excitations
whose threshold is at least that enthalpy. All 384 substeps of the 925-minute
feed-only diagnostic pass; the minimum opportunity/outcome ratio is `4.058`.
These are **opportunities, not traced parent events**. The test does not
establish state-specific dissociation branching or a closed heat budget.
It protects the current source mapping from an obvious event-count or
enthalpy-floor contradiction, not from all plasma-chemistry error.
The full backend suite passed with runtime warnings treated as errors before
the final failure-path regression was added; the focused six-test audit and
long-horizon suite passes afterward. The new regression confirms that an
artificially underfunded H2O fragment or removed N2 parent opportunities
causes an explicit audit failure.

The full-duration aqueous HCN+CN− endpoint remains `7.5560521e-11 mol`
in the primary branch and `1.2875241e-10 mol` in the charged branch under
the explicitly hypothetical exposure/contact and pH extrapolation. The
carbon, nitrogen and thermal limitations above remain. The next mechanistic
gate is state-resolved source branching and nascent-energy accounting,
especially the excluded CO2/O(1D) channel; independently measured wet-spark
exposure and contact are still required before an absolute yield claim.

## Additional reviewed O(1D) bath partners (2026-09-27)

The [NASA/JPL Evaluation 17, Table 1-1, notes A8, A9 and A11](https://science.jpl.nasa.gov/documents/1489/JPL_10-6_Final_15June2011_20110620.pdf)
provides near-298 K O(1D) removal and product branches for N2O and NH3:
N2O has `4.95e-11 cm3 molecule-1 s-1` to N2+O2 and `7.75e-11` to 2NO;
NH3 has `2.50e-10` to OH+NH2. The isolated first-collision projector now
includes these explicit branches and checks C/H/N/O and charge closure.
The same review assigns CH4 total removal `1.75e-10`, but one product branch
leaves CH3O versus CH2OH unresolved; CH4 therefore contributes only a
reviewed removal hazard and makes product projection fail closed. These
are elementary collision data, not an HCN-target recipe or HCN fit.

On the same 925-minute, 384-step feed-only diagnostic, the maximum
simultaneous *unreviewed* bath mole fraction falls from `5.21769e-5` to
`2.32424e-5`. At the unchanged **hypothetical**, non-bounding common
unknown rates of `1e-9` and `1e-8 cm3 molecule-1 s-1`, maximum unknown
hazard shares fall from `0.001339`/`0.013234` to `0.000596`/`0.005932`.
This does not establish those rates as physical limits. Fate is still
unavailable for 383 later steps; among the largest remaining unknown
partners are HNCO `4.42e-8 mol` and H2O2 `3.27e-8 mol` at their separate
maxima. A bounded search did not identify source-sufficient near-room-
temperature O(1D) total rates **and product branches** for those partners.
Neither their material abundance nor the hypothetical stress points justify
renormalizing the known branches or injecting the excluded CO2 event.

The full-run primary and charged aqueous HCN+CN− endpoints remain
`7.5560521163e-11` and `1.2875240950e-10 mol`, respectively, with both
material ledgers closed and zero failed charged steps. This read-only
extension does not change the predictive chemistry, missing thermal ledger,
wet-spark exposure, gas-water contact evidence, or held-out assay status.
The complete backend suite passes with `RuntimeWarning` treated as an error;
the only emitted warning is the unrelated installed Starlette test-client
deprecation.

## Full-duration spark-coupling sensitivity (2026-09-27)

The one-hour spark sensitivity in `HCN_FORMATION_GATE.md` uses a different
CH4-rich feed and cannot be extrapolated to this 925-minute CO2-fed flask.
The long-horizon probe now accepts an explicit `--spark-coupling-factor`
while preserving its measured initial feed, duration, hypothetical common
`kLa=1e-4 s^-1`, 200 Td field, active-zone volume, and flagged pH-4-to-5
hydration extrapolation. The 0.5/1/2 factors correspond to hypothetical
electrical coupling fractions 0.15/0.30/0.60, **not** measured bounds.
The 384-step outcomes are:

| Coupling fraction | Declared electrical estimate (J) | Primary aqueous HCN+CN− (mol) | Charged aqueous HCN+CN− (mol) |
| ---: | ---: | ---: | ---: |
| 0.15 | 40,792.5 | `1.19253267585e-10` | `1.07986825518e-10` |
| 0.30 | 81,585.0 | `7.55605211630e-11` | `1.28752409501e-10` |
| 0.60 | 163,170.0 | `3.31624269246e-11` | `1.32693695771e-10` |

The primary endpoint spans `3.596`-fold; the charged endpoint spans
`1.229`-fold. The primary response is inverse over these points, despite
the larger electrical estimate, and the two gas-chemistry branches respond
differently. An opt-in gas-reaction trace below now identifies the within-
model near-cancellation behind the primary response; it does not validate
those reaction rates or the experimental spark. These three arbitrary
points are not a probability distribution or uncertainty interval.

Separate 192-step half/double runs reproduce the 384-step primary amounts
to `4.0e-5` and `4.9e-6` relative, and charged amounts to `7.7e-5` and
`2.5e-4` relative. Each run has closed primary/charged material ledgers,
zero failed charged substeps, and passes the accepted-source enthalpy and
N2 parent-opportunity screens. Both off-baseline runs still use the
explicitly extrapolated pH region; the thermal ledger remains open. The
assay was not loaded. This full-duration sensitivity strengthens the case
that independently measured deposited energy/exposure is necessary before
testing the HCN amount against a published assay.
The full backend suite passes after this probe was added, with
`RuntimeWarning` treated as an error; the only emitted warning is the
unrelated Starlette test-client deprecation.

## Full-duration primary HCN reaction-flux trace (2026-09-27)

The existing first-substep adaptive-interval pathway quadrature is now
available **opt-in** across all apparatus substeps with
`--trace-hcn-chemistry`. It integrates each GRI reaction's *signed net*
HCN contribution by two-point Gauss quadrature on the BDF intervals,
then compares the sum with the observed gas chemical change after the
modeled electron source and gas/water exchange are separated. This is not
microscopic forward and reverse event counting, and the trace never feeds
back into the ODE. A one-minute result below the `1e-15 mol` numerical
trace floor is explicitly unranked rather than interpreted.

At 925 minutes and 384 steps, all three hypothetical coupling cases pass
the quadrature-versus-inventory reconciliation. The largest total closure
residual is `9.07e-16 mol` (baseline), at most `1.13e-5` of that case's
observed net gas chemical HCN change. The primary aqueous endpoints are
unchanged by tracing. Signed net contributions are:

| Coupling fraction | Positive HCN contributions (mol) | Negative HCN contributions (mol) | Net gas chemical HCN (mol) | Net/positive |
| ---: | ---: | ---: | ---: | ---: |
| 0.15 | `2.56027753e-7` | `2.55900602e-7` | `1.27151243e-10` | `0.04966%` |
| 0.30 | `5.41971020e-7` | `5.41890535e-7` | `8.04847882e-11` | `0.01485%` |
| 0.60 | `7.32286842e-7` | `7.32251526e-7` | `3.53159761e-11` | `0.004823%` |

The largest positive modeled pathway at baseline is GRI reaction 275,
`CH3 + N <=> H2 + HCN` (`5.06e-7 mol` signed contribution); the largest
negative one is reaction 236, `H + HCN (+M) <=> H2CN (+M)`
(`5.42e-7 mol`). The smaller `HCNN + O <=> HCN + NO` pathway also forms
HCN. Increasing coupling raises both opposing sums, but the tiny positive
surplus shrinks. This explains the **internal** inverse primary response
without asserting that a real wet spark has the same branch rates or that
more electrical input physically lowers HCN yield. Gas-to-water capture,
aqueous CN− fate, and charged Pavan chemistry are outside this GRI gas-
reaction trace. Because the net is the difference of much larger terms,
key kinetic errors could easily dominate it; mechanism and discharge
constraints remain essential before an order-of-magnitude assay claim.
The focused long-horizon/core-chemistry tests and complete backend suite
pass with `RuntimeWarning` treated as an error; only the unrelated
Starlette test-client deprecation is emitted. The kinetic audit below
checks the provenance and validity domain of GRI reactions 275 and 236;
their flux dominance is a reason to audit them, not permission to fit
them to HCN.

## Low-temperature GRI HCN kinetics challenge (2026-09-27)

The full-run trace points to a **mechanism-domain problem**, not merely an
uncertain gas/water contact coefficient. The packaged [GRI-3.0 rate file](https://www.cantera.org/stable/examples/input/gri30.html)
contains Cantera zero-based reaction 275 (paper/CHEMKIN reaction 276),
`CH3 + N <=> H2 + HCN`, and reaction 236 (reaction 237),
`H + HCN (+M) <=> H2CN (+M)`. The primary collision projection has only
ground-state `N`, because unresolved excited products are excluded. Yet
[Pearce, Ayers & Pudritz (2019)](https://arxiv.org/pdf/1902.05574)
find direct `CH3 + N(2D) -> HCN + H2` below `1e-21 cm3 molecule^-1 s^-1`
at 298 K and omit it from their low-temperature network. They identify
`CH3 + N(4S) -> H2CN + H` as the supported triplet-surface route; the
historical ~10% direct-HCN branch inferred from H/H2 products could
instead proceed indirectly through H2CN. Their calculation does **not**
prove the GRI direct channel is exactly zero under every spark condition,
but it prevents treating that combustion-network branch as established
ground-state chemistry at this temperature.

At the feed card's assumed 298.15 K, 100 kPa and reported initial gas
fractions with N2 balance, local Cantera 3.2 reports forward effective
coefficients `1.06153e-10 cm3 molecule^-1 s^-1` for the H2CN branch,
`1.68107e-11` for direct HCN (13.671% of the two GRI N+CH3 branches),
and `1.63147e-12` for H+HCN association. These are **GRI evaluations**,
not measured wet-flask rates. The 1998 primary QRRK study by
[Sumathi & Nguyen](https://quantchem.kuleuven.be/cmms/Articles/JPCA/jp981542u.pdf)
calculates the 300 K high-pressure total H+HCN disappearance limit as
`2.2e7 cm3 mol^-1 s^-1` (`3.65e-17 cm3 molecule^-1 s^-1`), mainly adduct
stabilization. The GRI feed-state value is about `4.47e4` times larger.
These estimates have different methods/pressure limits and no direct
like-for-like rate measurement was found; this is a **conflict requiring
resolution**, not authority to replace GRI with that single theoretical
value. The main [2019 HCN study](https://arxiv.org/pdf/1902.05574)
also explicitly distinguishes N(4S) from N(2D), which the GRI primary
source presently does not model separately.
GRI also has **no** `H2CN + H -> HCN + H2` elementary step among its four
H2CN reactions, even though the [2019 low-temperature network](https://arxiv.org/pdf/1902.05574)
calculates `1.8e-11 cm3 molecule^-1 s^-1` for it at 298 K. Merely
suppressing one GRI rate is therefore not a complete alternative network:
the intermediate's subsequent fate and the N electronic state both need
source-consistent treatment.

An internal, opt-in, equation-guarded GRI rate-stress switch changes only
these two primary gas reactions; the ordinary cached mechanism and the
charged Pavan shadow are not changed. It is exposed solely in the feed-only
probe as `--gri-direct-hcn-rate-factor` and
`--gri-hcn-h-association-rate-factor`, **not** as calibrated parameters.
Removing the direct channel (`0`) and lowering H+HCN association by
`1e-5` are structural/rate stress cases motivated by the cited conflict,
not measured uncertainty bounds or adopted chemistry. At one hour and
48 steps, primary aqueous HCN+CN− is `1.16467775e-11 mol` at baseline,
`8.02495450e-13 mol` with the direct branch off, and
`1.48479402e-8 mol` with slower association. Both modified cases close
the primary/charged material ledgers; the Pavan aqueous amount stays
`1.14400613e-11 mol` because that shadow intentionally retains its own
rates. The direct-off gas reaction trace reconciles with the observed
chemical HCN change. These differences are **model sensitivity**, not
predicted experimental yield changes. Full-duration stress and numerical
checks follow below; no held-out HCN assay was loaded.

At the reported 925-minute duration and 384 integration steps, the same
unfitted structural/rate stresses give:

| Primary GRI scenario | Aqueous HCN+CN− (mol) | Ratio to unchanged GRI | Closed primary/charged ledgers? |
| --- | ---: | ---: | --- |
| Unchanged GRI | `7.55605215e-11` | `1` | yes |
| Direct N+CH3 to HCN channel off | `4.94095141e-12` | `0.06539` | yes |
| H+HCN association rate × `1e-5` | `4.26850386e-7` | `5,649` | yes |

Independent 192-step stress reruns give `4.94102923e-12 mol` and
`4.26850568e-7 mol`, respectively: 192-to-384 relative changes of
`1.575e-5` and `4.273e-7`. Both reruns close the primary and charged
material ledgers. The focused long-horizon/core prebiotic tests and the
complete backend suite pass with `RuntimeWarning` treated as an error;
unrelated dependency deprecations are emitted.

The charged Pavan result remains `1.28752410e-10 mol` in both stresses,
because its separate mechanism was not modified. The result is
**structurally unstable by much more than one order of magnitude** under
source-motivated low-temperature alternatives. Numerical conservation
holds but cannot resolve the kinetic conflict. The complete apparatus
thermal ledger is still open, and this test has not used the held-out
HCN assay. Resolving N(4S)/N(2D), H+HCN pressure-dependent kinetics,
and H2CN fate with an internally consistent low-temperature mechanism is
now a prerequisite to any claimed one-order-of-magnitude validation.

### Independent H + HCN association applicability check (2026-09-27)

The read-only `backend/scripts/audit_hcn_association_kinetics.py` now checks
the dominant GRI `H + HCN (+M) <=> H2CN (+M)` loss against a second,
independently evaluated kinetic source, using only the held-out-assay-free
Barth feed card. [Tsang and Herron (1991), pp. 654-655](https://doi.org/10.1063/1.555890)
explicitly state that no direct measurements were available for this
reaction. Their recommended **estimated** high-pressure expression is
`k_inf = 5.5e-11 exp(-2438/T) cm3 molecule^-1 s^-1`; they also give
N2- and CO2-specific low-pressure coefficients. The printed source table and
continuation were visually checked in the original NIST PDF. Its stated
temperature range begins at 300 K, so application to the hypothetical
298.15 K feed is a 1.85 K extrapolation, not a within-domain validation.

At the declared 100 kPa initial feed, GRI evaluates to
`1.63146997e-12 cm3 molecule^-1 s^-1`, while Tsang and Herron's estimated
`k_inf` is `1.54559864e-14` in the same units: a factor `105.56` lower.
Repeating the rate comparison at the source's 300 K lower boundary while
holding the same hypothetical pressure and composition gives a factor
`99.64`; the approximate hundredfold conflict is not created by the 1.85 K
extrapolation.
The N2 and CO2 terms of their low-pressure expression sum to
`4.73360153e-17 cm3 molecule^-1 s^-1` when multiplied by the calculated
initial number densities. That is **not** a complete mixture rate: H2, H2O,
and O2 collider terms are not supplied by this intake, and no falloff rate
is inferred. The source describes a factor-100 uncertainty at 500 K for
its provisional estimate; this is not a statistical 298 K uncertainty
interval, nor a measured upper bound. The earlier Sumathi-Nguyen 1998
theoretical limit is methodologically different again. These independent
discrepancies establish that the current GRI low-temperature HCN sink is
not qualified, but do not select a replacement rate or justify calibrating
it against the held-out HCN assay.

The audit is included in the feed-only long-horizon report as
`unmodified_gri_h_hcn_initial_feed_reference_audit`; it is omitted from
altered-feed controls and marked separately when an opt-in GRI rate stress
means it no longer matches the active mechanism. It does not modify the ordinary or
charged mechanisms. Focused source-form and one-minute integration tests
pass, including no-assay and missing-collider guards. A source-supported
mixture-dependent falloff law, or direct low-temperature measurements in a
relevant bath, is still required before the rate enters a predictive model.

## H2CN fate and partial low-temperature branch overlay (2026-09-27)

The full-duration baseline H2CN gas trace shows why suppressing direct
HCN formation is not a complete alternative: GRI `CH3 + N -> H + H2CN`
contributes about `3.194e-6 mol` H2CN, and `H + HCN -> H2CN` about
`5.419e-7 mol`, while GRI `H2CN + N -> CH2 + N2` removes about
`3.736e-6 mol`. Its small net is ill-conditioned: historical two-point
pathway quadrature differs from the observed net by `6.04e-15 mol`, only
`1.62e-9` of gross H2CN turnover but about `11%` of the tiny net.
The diagnostic now reports separate net-scale and gross-scale closure
errors, labeling a gross-reconciled but net-ill-conditioned trace
explicitly. Five-point Gauss quadrature is used when tracing H2CN; it
does not remove the one-hour net ill-conditioning. The adaptive BDF
history and rate evaluation remain possible error sources. These gross
rankings are useful for finding missing fates, not for claiming a
precise H2CN net yield.

Two primary studies provide independent 298 K kinetic constraints:

| H2CN reaction | Current GRI at assumed feed state | Independent low-temperature value |
| --- | ---: | ---: |
| `H2CN + H -> HCN + H2` | absent | `1.8e-11 cm3 molecule^-1 s^-1` from [Pearce et al. (2019)](https://arxiv.org/pdf/1902.05574), calculated singlet route; [Hébrard et al. (2012)](https://doi.org/10.1051/0004-6361/201218837) recommend `6e-11` |
| `H2CN + N(4S) -> CH2 + N2` | `5.072e-11 cm3 molecule^-1 s^-1` | `3e-11`, estimated branch in [Hébrard et al. (2012)](https://doi.org/10.1051/0004-6361/201218837) |
| `H2CN + N(4S) -> HCN + NH` | absent | `1e-11`, estimated branch in [Hébrard et al. (2012)](https://doi.org/10.1051/0004-6361/201218837); [Pearce et al. (2019)](https://arxiv.org/pdf/1902.05574) calculate `9.4e-13` |

The 2012 N+H2CN split sums to `4e-11 cm3 molecule^-1 s^-1`, close
to the reported overall room-temperature experiment (`4.4e-11`), but
its individual branches are estimates with substantial uncertainty.
The two sources' H+H2CN-to-HCN rates differ by a factor of 3.3, and
their N+H2CN-to-HCN rates by about 10.6. Hébrard et al. also allow
`H2CN + H -> HNC + H2` (`1.2e-11`); GRI has no HNC species, so that
competing fate is **not** represented in the inventory-driving branch.
The papers describe different environments and are not wet-spark
product-yield calibration data.

The 298 K [Nesbitt, Marston & Stief experimental study](https://doi.org/10.1021/j100375a035)
reports only a lower bound of `7e-11 cm3 molecule^-1 s^-1` for the
overall H2CN+H loss (inferred using D2CN as well, with no observed
isotope effect); it did **not** measure the HCN/HNC branching fraction.
Hébrard's estimated HCN and HNC branches sum to `7.2e-11`, just above
that bound. Pearce's `1.8e-11` HCN-only theoretical channel is not a
measurement of total loss. These facts support testing the competing
HNC exit but do not determine its amount or fate in the wet reactor.
Hébrard et al. identify H+HNC -> H+HCN isomerization as a major
kinetic uncertainty, so assigning instantaneous HNC conversion to HCN
would be unjustified.

An opt-in Cantera overlay now adds the Pearce 2019 H+H2CN route with
its reported constant 50–400 K forward rate and separately offers
the Hébrard 2012 **room-temperature** N(4S)+H2CN split. The latter
replaces GRI reaction 237's forward coefficient with the `3e-11`
branch and adds the `1e-11` NH+HCN branch; Cantera derives reverse
rates from the existing GRI thermochemistry. The branch profile is
guarded to 290–310 K. No modified reaction feeds the ordinary GRI
run, the charged Pavan shadow, or the held-out HCN assay. These
overlays are partial kinetic comparisons, not validated models.

At the reported 925-minute duration, 384 steps, fixed hypothetical
spark coupling and `kLa`, all primary/charged material ledgers close:

| Diagnostic primary profile | Aqueous HCN+CN− (mol) | Relative to unchanged GRI |
| --- | ---: | ---: |
| Unchanged GRI | `7.55605215e-11` | 1 |
| Add Pearce H+H2CN route only | `8.48481825e-11` | 1.123 |
| Direct N+CH3 HCN off; H+HCN association × `1e-5`; add Pearce route | `8.31563234e-8` | 1,101 |
| Same, plus Hébrard N+H2CN split | `7.52147835e-7` | 9,954 |

For the last profile, 192 steps give `7.52147602e-7 mol`, a
`3.095e-7` relative difference from 384. The 384-step HCN and H2CN
gas-reaction traces reconcile with their observed chemical changes;
the largest signed HCN-forming contributions are `7.870e-7 mol` from
the new N+H2CN branch, `7.259e-8 mol` from H+H2CN, and
`2.498e-8 mol` from HCNN+O. The charged branch remains
`1.28752410e-10 mol` because it was not modified. This extreme
source-sensitive span is decisive evidence **against** claiming that
the current simulated HCN amount is quantitatively predictive.
For the last profile, independent full-duration carbon-removed and
nitrogen-removed controls each return zero primary and charged aqueous
HCN+CN− at 192 steps with closed material ledgers and the source-energy
floor check passed. Those controls exclude spontaneous cyanide creation
without the requisite elements; they do not validate absolute yield.
The pathway aggregator now records per-substep quadrature methods: the
first apparatus substep traces H2CN at five-point order even when later
HCN-only substeps use two-point order, and this valid mixture is reported
rather than rejected.
The focused long-horizon regressions and complete backend test suite pass
with `RuntimeWarning` treated as an error; separate dependency
deprecation warnings remain.
The overlay still omits HNC, excited N(2D), other low-temperature
competing fates, independent spark exposure, and a closed apparatus
thermal budget. The next gate is a state-resolved, thermochemically
consistent HCN/HNC/H2CN network and independent discharge constraints;
do not fit any of these rates to the held-out HCN assay.

### HNC thermochemistry intake check

The already imported [Kamarinopoulou et al. (2024) supplementary
mechanism](https://pmc.ncbi.nlm.nih.gov/articles/PMC10980260/) contains
an HNC NASA7 record, but it is **not** suitable for direct transfer into
GRI without reconciliation. Cantera evaluation of that record gives
`H°(HNC, 298.15 K) = 638.148 kJ mol^-1`; its same-file HCN record gives
`128.538 kJ mol^-1`, so the implied isomer gap is `509.611 kJ mol^-1`.
GRI HCN gives `130.815 kJ mol^-1` at the same temperature. The
[NIST CCCBDB evaluated HNC formation enthalpy](https://cccbdb.nist.gov/exp2x.asp?casno=6914074)
is `194.33 ± 9.00 kJ mol^-1`, implying a roughly `63.5 kJ mol^-1`
HNC–GRI-HCN gap. The supplementary HNC record is about `444 kJ mol^-1`
above that evaluated HNC value. Its 298.15 K `Cp = 34.708 J mol^-1 K^-1`
and `S = 227.292 J mol^-1 K^-1` also differ from NIST's
`40.03 J mol^-1 K^-1` and `205.21 J mol^-1 K^-1`, respectively. Thus
changing only the NASA7 enthalpy offset would not establish a consistent
ground-state record. The source of this discrepancy has not
been established; it is not automatically corrected or blended with GRI.
Visual inspection of page 26 of the original supplementary PDF confirms
that its HNC NASA7 enthalpy coefficient is printed as `7.55555405E+04`;
the mismatch was not introduced by the local text extraction or CHEMKIN
conversion. Whether the source used a different state/reference or
contains an upstream data error remains unresolved. The Data S1 header
names several possible thermochemistry libraries but does not identify
which one supplied HNC, so it does not reconcile the ground-state value.
The record is absent from the inventory-driving GRI phase but is present
in the non-propagating Kamarinopoulou neutral shadow; that shadow now
explicitly flags HNC-sensitive rates as unvalidated. Any HNC branch must
first use a provenance-checked, temperature-valid
thermochemical record and pass a ground-state enthalpy/isomer-gap test.

A separate room-temperature HNC species constructor now
uses NIST's evaluated 298.15 K enthalpy (`194.33 ± 9.00 kJ mol^-1`),
entropy (`205.21 J mol^-1 K^-1`) and heat capacity
(`40.03 J mol^-1 K^-1`) in Cantera's constant-`Cp` form. The allowed
290–310 K interval is an engineering restriction on this local
approximation, **not** a source-validated thermal range. The NIST
entropy is defined at [1 bar standard pressure](https://cccbdb.nist.gov/thermox.asp);
the constructor converts it by the ideal-gas `-R ln(p2/p1)` relation to
`205.10056 J mol^-1 K^-1` at GRI/Cantera's 1 atm reference pressure.
The enthalpy and heat capacity are unaffected by this ideal-gas pressure
change. The implied HNC–GRI-HCN enthalpy gap is `63.515 kJ mol^-1`, within the NIST
uncertainty of the `62.903 kJ mol^-1` gap from the independent
[ATcT v1.122 HNC and HCN formation enthalpies](https://atct.anl.gov/Thermochemical%20Data/version%201.122/reaction/?c1=-1&c2=1&species_number%5B%5D=428&species_number%5B%5D=46).
These reference properties and cross-check do not specify an HNC kinetic
fate, gas-water transfer, or quantitative yield; HNC remains excluded from
the **default** primary apparatus. The separate opt-in gas comparison below
does not promote it to validated chemistry.
An isolated thermochemical audit combining this HNC candidate with GRI
H, H2, H2CN and HCN gives `ΔH(298.15 K) = -271.000 kJ mol^-1` for
`H2CN + H -> HNC + H2` and `-334.515 kJ mol^-1` for the HCN exit.
Their difference is exactly the `63.515 kJ mol^-1` isomer gap. These
values are close to the approximate `-276` and `-331 kJ mol^-1`
reaction energies in [Hébrard et al. (2012)](https://doi.org/10.1051/0004-6361/201218837),
but this is a consistency check, not a validation of their estimated
branching kinetics. HNC+H isomerization, gas-water fate, and a complete
independently justified HNC network remain absent from the default solver;
an estimated, opt-in HNC+H diagnostic is tested below.
The four HNC thermochemistry tests and complete backend suite pass with
`RuntimeWarning` treated as an error; existing dependency deprecations
remain. No held-out HCN assay was loaded for this intake.

### Opt-in HCN/HNC competing-exit diagnostic (2026-09-27)

An opt-in, **non-default** GRI extension now appends NIST-grounded HNC
thermochemistry and both [Hébrard et al. (2012)](https://doi.org/10.1051/0004-6361/201218837)
room-temperature estimates for `H2CN + H`: `6.0e-11 cm3 molecule^-1 s^-1`
to `HCN + H2` and `1.2e-11` to `HNC + H2`. They are reversible in Cantera;
reverse rates use GRI plus the narrow 290–310 K NIST HNC reference state.
The split is **estimated, not experimentally measured**. The switch is
mutually exclusive with the Pearce HCN-only `H2CN + H` overlay and does not
alter the default GRI mechanism or charged Pavan shadow. HNC is added to
the atom ledger and ideal-gas internal-energy inventory, but no HNC Henry
transfer, aqueous conversion, or H+HNC isomerization is assigned.

With the previously declared structural stress (direct GRI `N + CH3`
HCN exit off; GRI `H + HCN` association multiplied by `1e-5`; Hébrard
`N + H2CN` split on), the reported-duration feed-only run gives:

| Integration steps | Primary aqueous HCN+CN− (mol) | Net gas-chemical HNC (mol) | Charged aqueous HCN+CN− (mol) |
| ---: | ---: | ---: | ---: |
| 192 | `8.23245642967e-7` | `4.45140730150e-8` | `1.28772211652e-10` |
| 384 | `8.23246404297e-7` | `4.45144243045e-8` | `1.28752409501e-10` |

The primary aqueous amount changes by `9.25e-7` relatively between
resolutions; net gas-chemical HNC changes by `7.89e-6`. HCN and HNC
reaction-flux quadratures reconcile at both resolutions; H2CN's large
opposing turnover reconciles only on a gross-flux scale, not its tiny net.
Primary and charged atom/charge ledgers close, the charged branch has no
failed substeps, and the accepted source thermochemical-floor check passes.
The diagnostic ideal-gas internal-energy inventory includes HNC; a
one-minute regression confirms coverage of every nonzero primary gas species
in that test. This is not a closed gas–liquid–wall–radiation thermal balance.
At 192 steps, separate full-duration carbon-removed and nitrogen-removed
controls each yield zero primary aqueous HCN+CN−, zero charged aqueous
HCN+CN−, and zero net gas-chemical HNC, with both material ledgers closed.
At 384 steps this profile gives about `1.095` times the aqueous amount of
the previous Pearce-HCN-only plus Hébrard-N diagnostic; the difference
conflates the changed HCN forward rate with the new HNC exit and is **not**
an isolated HNC branching effect. The charged amount is unchanged because
the charged shadow does not use these rate overlays.

This is numerical stability for an **incomplete mechanism**, not an HCN
yield prediction or validation. HNC's missing water and H+HNC fates could
change both the collected HCN+CN− and gas inventory. Its estimated
branching is from an atmospheric-chemistry network, not a wet-spark
measurement. The full apparatus thermal balance and independent spark
exposure remain open, and the held-out HCN assay remains unloaded. A bounded
search for HNC-specific aqueous uptake/kinetics found HCN and HNCO data,
which are different species and were not transferred to HNC. The next
mechanism step needs source-supported HNC gas and aqueous fates; otherwise
the diagnostic must remain explicitly incomplete.

### Estimated relaxed-HNC H-atom isomerization sensitivity (2026-09-27)

The same opt-in HNC profile can additionally enable `H + HNC <=> H + HCN`
using [Hébrard et al. (2012), Appendix A.2](https://doi.org/10.1051/0004-6361/201218837):
`k = 4.03e-11 exp(-1200/T) cm3 molecule^-1 s^-1`, or `7.20026e-13`
at 298.15 K. The paper recommends this relaxed-HNC rate from a
[theoretical QRRK study](https://doi.org/10.1021/jp981542u) and barrier
estimates, and states a factor-of-ten uncertainty at 300 K. It is **not**
a wet-spark rate measurement. The reversible Cantera reaction uses the
same NIST-grounded HNC/GRI HCN thermochemistry, is restricted to 290–310 K,
and is absent from the default mechanism and charged Pavan shadow.

In the same declared full-duration structural-stress profile as the table
above, adding only this isomerization diagnostic gives:

| Integration steps | Primary aqueous HCN+CN− (mol) | Net gas-chemical HNC (mol) | Charged aqueous HCN+CN− (mol) |
| ---: | ---: | ---: | ---: |
| 192 | `8.61785179568e-7` | `8.81630407165e-13` | `1.28772211652e-10` |
| 384 | `8.61786171109e-7` | `8.81615746248e-13` | `1.28752409501e-10` |

The 192-to-384 primary aqueous difference is `1.15e-6` relatively.
At 384 steps, enabling isomerization raises collected aqueous HCN+CN−
by `4.68%` compared with the otherwise identical opt-in profile and
lowers net gas-chemical HNC by about `5.05e4`-fold. HCN and HNC gas
reaction traces reconcile; H2CN remains gross-flux reconciled but
net-ill-conditioned. All primary gas species are covered by the ideal-gas
energy inventory in these runs, material ledgers close, the source-energy
floor check passes, and charged substeps do not fail. Separate full-duration
carbon-removed and nitrogen-removed 192-step controls both give zero
aqueous cyanide and HNC with closed material ledgers. The matched one-hour
48-step comparison shows `2.85197656e-8` versus `2.98984153e-8 mol`
aqueous HCN+CN− without/with isomerization, with both gas-reaction traces
reconciled.

These percentages are **conditional model sensitivities**, not an error
bar or a real experimental yield. No factor-of-ten rate sweep has yet been
performed, and neither nascent hot-HNC relaxation nor its uptake into
water is resolved. The complete thermal balance and independently measured
spark exposure remain open; the held-out aqueous HCN assay was not loaded.
The focused mechanism tests, a one-hour run with `RuntimeWarning` treated
as an error, and the complete backend suite pass after this addition.
Remaining third-party deprecations are unrelated to the chemistry solver.

### Conditional aqueous glycolonitrile endpoint screen (2026-09-27)

A separate **read-only** screen now tests whether an omitted cyanohydrin
reservoir deserves a coupled aqueous mechanism. [Schlesinger and Miller's
measured aqueous glycolonitrile equilibrium](https://doi.org/10.1021/ja00792a043)
reports `log10(K_dissociation/M) = 1.66 - 2184/T` over 25–70 °C; its
association rate was separately measured only in dilute acetate-buffered
solution at pH 3.8–4.7 and 25–45 °C. The screen holds the final modeled pH
and liquid volume fixed, conditionally treats aqueous `CH2O + CH4O2_AQ`
as the formaldehyde pool, and solves `HCN + formaldehyde ⇌ HOCH2CN` at
the endpoint. It changes **no** simulated inventory, reaction network,
source term, or assay comparison. Full-text verification of the equilibrium
constant's formaldehyde-pool convention is still required; the paper's
abstract alone does not settle its precise free/hydrated basis. The bounded
source card is
[`glycolonitrile_schlesinger_miller_1973_v1.json`](../backend/data/kinetics/glycolonitrile_schlesinger_miller_1973_v1.json).

The unchanged 925-minute, 16-cycle, 24-substep feed-only baseline passes
its numerical screen and material ledger. Its endpoint has
`7.55605212e-11 mol` aqueous HCN+CN− and `5.30120022e-8 mol` free plus
hydrated aqueous formaldehyde. With the modeled final ideal pH frozen,
the conditional equilibrium screen gives `2.48082777e-11 mol`
glycolonitrile, or **32.83% of the endpoint cyanide pool**. The 60-second
two-substep check gives only `1.80276440e-31 mol` conditional
glycolonitrile from `4.76682286e-20 mol` aqueous cyanide, below the
apparatus reporting floor. Both figures are *scenario diagnostics*, not
predicted glycolonitrile product amounts or a corrected HCN yield.
Separate 60-second carbon-removed and nitrogen-removed feed controls
return zero aqueous cyanide and zero conditional glycolonitrile, with
their existing numerical screens passing.
They omit time-dependent production/consumption, gas-water feedback,
liquid pH evolution during association, other aqueous sinks, and any
possible glycolonitrile decomposition. The published `k2` is not applied
through a modeled trajectory whose minimum ideal pH is 4.714, slightly
outside its 3.8–4.7 measured range. No held-out HCN assay enters the
screen.

The next step at this point was to verify the equilibrium-pool definition
from the primary text, then establish source-supported reversible kinetics over
the model's temperature/pH history and add a material-conserving,
time-dependent liquid state **only if** those checks pass. This aqueous
branch cannot by itself address the larger independently unmeasured spark
exposure, GRI HCN kinetic instability, or open thermal balance.
The subsequent rate-aware check below supersedes this equilibrium-only
priority assessment.
The focused glycolonitrile and long-horizon tests pass with warnings
treated as errors; the full backend suite passes in normal warning mode.
The full strict-warning collection currently stops on an unrelated
Starlette `httpx`/`httpx2` deprecation, rather than a chemistry failure.

### Glycolonitrile equilibrium-versus-timescale correction (2026-09-27)

The [indexed original-paper text](https://www.electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20the%20American%20Chemical%20Society%20US/1973%20%20%28vol%20095%29/11%20%20%283445-3824%29/3729-3735.pdf)
on pp. 3731 and 3734 distinguishes `(H2CO)aq`, predominantly hydrated,
from **unhydrated** formaldehyde in the kinetic law. This supports using
free `CH2O` plus `CH4O2_AQ` for the **equilibrium** pool but free aqueous
`CH2O` alone for the measured forward rate. The original PDF could not be
rendered locally for visual verification, so the source card identifies
this as indexed primary-paper OCR evidence. The [1974 correction](https://doi.org/10.1021/ja00812a605)
changes benzaldehyde table entries and the preferred glycolonitrile name,
not the constants used here.

The opt-in read-only screen now additionally applies the paper's
`k2 = 3.5e5 M^-1 s^-1` at 298.15 K to each stored right-endpoint liquid
state as `k2 [CN-] [free CH2O]`; no reaction or product enters the solver.
This measured coefficient came from dilute acetate-buffered pH 3.8–4.7
water. The modeled 925-minute run's terminal pH is just above that range;
**all** of its 192/384 cyanide-bearing steps are outside the measured pH
domain. Values below are explicit extrapolations, not source-bounded
kinetics or bounds on a real experiment.

| Integration steps | Aqueous HCN+CN− at end (mol) | Conditional gross one-way association opportunity (mol) | Outside measured pH |
| ---: | ---: | ---: | ---: |
| 192 | `7.55623508e-11` | `2.56430402e-14` | 192/192 steps |
| 384 | `7.55605212e-11` | `2.55507032e-14` | 384/384 steps |

The step-opportunity changes by about 0.36% between these two resolutions,
and is below the apparatus's `1e-12 mol` numerical reporting floor. At the
terminal state, only `4.17778494e-11 mol` of the `5.30120022e-8 mol`
aqueous formaldehyde pool is free `CH2O`, and the modeled cyanide-ion
fraction is `6.00503e-5`. Holding **that endpoint state** for all 925
minutes would give an initial one-way rate of `1.75127e-8 s^-1`, a
characteristic time of `5.71e7 s`, and only `9.71484e-4` conditional
forward conversion fraction. These are not trajectory-integrated yields;
the step-opportunity also assumes reactants were present throughout each
step and neglects reverse reaction. It is not a strict upper or lower
bound. Short carbon-removed and nitrogen-removed controls retain zero
aqueous cyanide and report no association opportunity. Material ledgers
remain closed, and focused tests pass with warnings treated as errors.

**Priority correction:** the earlier **32.83% equilibrium potential**
must not be read as a 925-minute sink. The rate-aware screen suggests it
is much slower *under these modeled conditions*, but because the entire
trajectory leaves the paper's measured pH range, it cannot establish the
real glycolonitrile amount. Do not couple or fit this sink to an HCN assay
from this evidence alone. A wider-pH kinetic source or direct measurement
is needed before it enters the inventory-driving network. For the present
HCN gate, independently constrained spark exposure, low-temperature HCN
gas kinetics, gas-water contact, and thermal-energy accounting remain the
higher-priority quantitative blockers.
