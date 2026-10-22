# Dry-flow CH4/N2 HCN benchmark gate (2026-09-29)

This is the project's **primary HCN development target**. The wet
Miller–Urey-style apparatus is no longer the validation goal. The target is
the dry, flowing CH4/N2 dielectric-barrier-discharge experiment of
[Kamarinopoulou, Wittreich & Vlachos (2024)](https://doi.org/10.1126/sciadv.adl4246).
That paper reports controlled feed, flow, electrical power and inline gas
product measurements. It is a more directly comparable target than a wet
boil/spark/condense loop, but still has unresolved microdischarge exposure.
Its [published microkinetic methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10980260/)
represent an isothermal plug-flow
reactor by 100 stirred cells. Our evolving-packet integration is a
**Lagrangian plug-flow idealization**, not a single stirred tank, so adding
one well-mixed reactor would not improve methodological agreement. The packet
still assumes a fixed volume and uniform hot temperature and does not
resolve axial/radial gradients or streamer occupancy. In the 32-step
scenario, the final ideal-gas pressure inferred from its inventory is
101340.63 Pa, only 0.0154% above the 101325 Pa inlet; this particular
fixed-volume versus constant-pressure difference is small. That narrow
check does not validate the temperature, flow convention or plasma exposure.

## What would count as success

Build an input-frozen gas-flow model at one or more explicitly matched paper
conditions. Confirm the feed fraction, actual flow/residence convention,
pressure, temperature, power and sampled outlet basis; independently bound
the active discharge exposure and state-resolved fragment sources. Verify
numerical convergence, atom/charge accounting, source-energy bounds,
spark-off/no-carbon/no-nitrogen controls, and outlet CH4 conversion and
competing products. Then compare held-out HCN values and trends across feed,
flow and power conditions. An order-of-magnitude agreement may be reported
only with the precise matched condition and an uncertainty/sensitivity
envelope; it must not be inferred from a parameter fitted to that same HCN
measurement. The model already contains HCN reactions in GRI-3.0, so this
is **formation prediction**, not blind product discovery.

Published HCN yields, reaction-path conclusions and the paper's effective
~1.3 ms plasma exposure are external validation/context, **not** independent
source inputs. The paper's ~0.47 s nominal gas residence and much shorter
model-equivalent plasma time should not be conflated. The 0.3% streamer-volume
estimate is inferred through model/product comparison, not an independent
measurement of the streamer's space-time occupancy.
The [paper's methods](https://rcastoragev2.blob.core.windows.net/ea712d57e86c0bac7a14d5ec18c0b746/sciadv.adl4246.PMC10980260.pdf)
also explicitly simplify its neutral chemistry to ground electronic states
and assume direct N2 dissociation supplies N atoms; its E/N, electron density,
gas temperature, and effective residence are selected partly against product
trends. Reproducing those assumptions could be a **literature-model
replication**, not independent validation of the underlying fragment states
or a prediction with HCN held out.

## Current active-code result, not archived wet or dry results

`backend/scripts/run_dry_flow_hcn_probe.py` takes a 1% CH4/N2 gas packet
through the installed EEDF/collision and GRI neutral network. Its current
scenario uses the paper's 1.57 cm3 reactor volume and 28.1 W nominal point,
but treats a **downstream probe reading** of 282 C as a **uniform** 555.15 K
gas temperature,
uses an **assumed** 298.15 K flow reference, and sets an **assumed** 200 Td
reduced field and whole-packet discharge exposure. The hot-packet residence
is 0.25242 s under that flow convention; this is a scenario calculation,
not a recovered experimental exposure.
The [paper's thermometry method](https://pmc.ncbi.nlm.nih.gov/articles/PMC10980260/)
placed an optical-fiber probe **1 cm above the reactor**, with a reported
0.09–0.38 s gas transit from reactor outlet to probe. Supplementary Figure S3
shows a roughly 282 C steady reading after at least 10 minutes in a **separate
28.6 W, 1% CH4 run**; Figure 1's nominal comparison uses 28.1 W. Thus our
555.15 K uniform temperature is a proxy, not a directly measured in-reactor
field or an exactly matched Figure 1 temperature. The evaluator now reports
the thermometry-run mismatch explicitly even when feed, nominal residence,
and power match Figure 1.

After the correction that prevents unresolved total N2 dissociation from
becoming `2 N(4S)`, the active 16-step runs yield:

| Source forcing | Outlet HCN per packet |
| --- | ---: |
| Impulse substeps | 5.34189e-16 mol |
| Inventory-limited | 2.61610e-16 mol |
| Evolving EEDF | 2.61734e-16 mol |

All are **below** the project's `1e-12 mol` HCN reporting floor. Earlier
dry-flow positive amounts elsewhere in the historical formation gate are
archived, superseded by the N2 product-state correction, and must not be
used as a validation claim. This is not an experimental zero or a detection
limit; it is only the active model result for the declared scenario.

The new, read-only outlet diagnostics report every gas species amount and
mole fraction, methane conversion, feed-carbon-normalized product yield, and
converted-carbon-normalized selectivity. The latter two follow the algebra of
the paper's supplementary equations S4 and S7, but **are not identical to its
GC-based assay**: our all-species accounting includes radicals and has no
wall-deposited aromatics/coke. At 16 evolving-EEDF steps, the scenario converts 0.172315 of inlet
methane and has 2.19634e-8 mol C2H6 per packet, while HCN remains at trace
level. Thus the model activates methane but does not have a qualified primary
N2 fragment supply. The active neutral source IDs are CH4 channels only;
`N2_dissociation_outcome_000` is explicitly excluded with reason
`unresolved_electronic_product_state`. No published HCN or competing-product
measurement was used to produce these diagnostics.

At 32 and 64 evolving-EEDF steps, methane conversion changes from
0.172317785 to 0.172319074 and C2H6 from 2.19637464e-8 to
2.19639315e-8 mol. The material ledger and allocated collision-energy
upper bound close in both cases. The HCN trace remains below the reporting
floor, so a small relative difference between its tiny numerical endpoints
does **not** demonstrate a converged HCN yield. The full thermal balance,
microdischarge exposure and flow-temperature basis remain open.
The paper's 22 kHz AC drive has a period of about `45.45 us`. The assumed
`0.25242 s` hot packet spans about `5553` AC periods: 16 evolving-EEDF
chemistry substeps average about `347` periods each, and even 64 substeps
about `87` periods each. The probe now reports these ratios and explicitly
marks individual nanosecond streamer pulses **unresolved**. Substep
convergence of methane and ethane is therefore convergence of this
time-averaged source approximation, not proof that pulse/afterglow kinetics
or local streamer exposure are physically resolved.

The assay-free `backend/scripts/audit_dry_feed_convergence.py` extends this
temporal check to the 0.1%, 1%, and 7% feed cases, at 16, 32, and 64
evolving-EEDF steps. Predicted methane conversion at 64 steps is
`0.17138919`, `0.17231907`, and `0.17526937`, respectively. The 32-to-64
changes are at most `8.94e-6` absolute conversion fraction (under `0.001`
percentage point); ethane-selectivity changes are at most `1.48e-6`, and
hydrogen-selectivity changes at most `4.51e-7`. The high-minus-low-feed
conversion slope remains **positive** at all three resolutions
(`0.003854`, `0.003871`, `0.003880`), opposite the published feed trend.
Every run retains closed material and allocated-collision-energy upper-bound
ledgers and excludes the unresolved N2 product-state source. Hence ordinary
chemistry substep refinement does not repair the discrepancy. This is not
convergence of streamer pulse physics or of the trace HCN yield.

### Frozen qualitative external check

An evaluator-only card at
`backend/data/benchmarks/kamarinopoulou_2024_dry_flow_qualitative_v1.json`
records the [paper's explicit observation](https://rcastoragev2.blob.core.windows.net/ea712d57e86c0bac7a14d5ec18c0b746/sciadv.adl4246.PMC10980260.pdf)
that HCN is favored over C2H6 for feed methane fractions up to 7% at the
Figure 1 conditions. The 1% feed, nominal 0.47 s and average 28.1 W case is
within that range. The new held-out check runs the assay-free model first and
hashes the frozen prediction before reading the observation. Its 16-step evolving-EEDF result has
HCN carbon yield `7.59430e-10` versus C2H6 carbon yield `0.127455`.
The **product ordering is opposite to experiment**. This is a directional
failure under a scenario with assumed uniform field, full-packet exposure,
and unresolved N2 states; it is not a like-for-like quantitative error ratio.
The article's maximum 72% HCN yield belongs to an unspecified point in this
comparison and is **not** assigned to this 1% case. No figure values or
uncertainty bars were digitized, and no experimental yield was used to tune
the model.

An additional assay-free feed sweep, `backend/scripts/audit_dry_feed_sweep.py`,
freezes 0.1%, 1%, and 7% methane cases before an evaluator reads the paper's
Figure 1A/B trends. At constant assumed 200 Td and uniform 555.15 K, the
model's C2H6 carbon selectivity rises with methane fraction, while its H2
hydrogen selectivity falls. Both directions agree with the
[paper's Figure 1 and text](https://www.osti.gov/servlets/purl/2471801):

| Feed CH4 | CH4 conversion | Model C2H6 carbon selectivity | Model H2 hydrogen selectivity | Model HCN carbon yield, below reporting floor |
| ---: | ---: | ---: | ---: | ---: |
| 0.1% | 0.171389 | 0.533183 | 0.359333 | 3.56110e-9 |
| 1% | 0.172315 | 0.739663 | 0.301089 | 7.59430e-10 |
| 7% | 0.175243 | 0.835849 | 0.273648 | 8.34706e-11 |

These **match the two reported directions**, but not necessarily their
measured magnitudes: model selectivity includes radicals and omits deposited
coke, whereas the paper includes measured stable products and coke. Hydrogen
selectivity is separately normalized by hydrogen atoms in converted CH4; its
modeled species sum to one under the closed H ledger. These are limited
non-HCN trend checks, not validation of
the discharge exposure or HCN mechanism. At every feed point HCN remains
below the 1e-12 mol reporting floor and below C2H6, contradicting the
paper's HCN-dominant ordering. Material and collision-energy upper-bound
ledgers close, and the unresolved N2 product-state source stays excluded.

The same [Figure 1A](https://rcastoragev2.blob.core.windows.net/ea712d57e86c0bac7a14d5ec18c0b746/sciadv.adl4246.PMC10980260.pdf)
shows CH4 conversion **falling substantially** as feed CH4 rises. Our model
instead predicts `0.171389 -> 0.172315 -> 0.175243` across 0.1%, 1%, and
7% feed. This is a second, independent directional failure, not just an HCN
source issue: the model's methane activation responds incorrectly to feed
composition. The evaluator now reads that observation only after the
assay-free sweep is frozen. Figure 1A's plotted conversion symbols also
support two deliberately broad, held-out **visual bounds**: conversion at
0.1% feed is above 70%, while at 7% feed it is below 30%. The 200 Td model
fails the low-feed bound (`17.1%`) and passes the high-feed bound (`17.5%`);
the 189 Td scenario likewise fails low feed (`16.7%`) and passes high feed
(`17.1%`). These thresholds are not digitized point estimates, measured
error bars, or statistical uncertainty; the temperature, flow and local
discharge exposure are not like-for-like, so no calibrated quantitative
error ratio is claimed. Exact source values and experimental SD/error bars
have **not** been digitized.

An additional assay-free **field-by-feed interaction check** at the paper's
`435 Td` optical estimate predicts CH4 conversions `0.17683`, `0.17459`,
and `0.17366` for 0.1%, 1%, and 7% methane. Unlike the `189` and `200 Td`
cases, the slope is negative, so field choice can reverse the *sign* of
this model's feed trend. But its total drop is only `0.00317` in conversion
fraction across the entire feed range, while the plotted low-feed conversion
is above the deliberately coarse 70% bound; the predicted low-feed value is
only `17.68%`. This is not a recovered trend *magnitude*. The 435 Td number
comes from a separate 28.6 W optical run, and the optical audit below does
not reproduce its diagnostic inference; it is not a verified Figure 1 local
field. All three modeled ledgers close, and the unresolved N2 atom branch
remains excluded. No HCN assay selected the field.
Internally, the CH4-target electron-fragment loss falls only from `0.24164`
to `0.22790` of feed CH4 over that seventyfold feed increase; net neutral
regeneration weakens from `0.06481` to `0.05424` of feed and partly cancels
the source change. This source-versus-regeneration balance explains the
model's nearly flat high-field trend without treating either term as a
measured reactor mechanism.

A new separate, frozen
[`Figure 1A conversion plot-band card`](../backend/data/benchmarks/kamarinopoulou_2024_figure1a_conversion_plot_bands_v1.json)
now records all seven gray-curve x positions as **visual, pixel-assisted readings**
from the original paper PDF rendered at 300 dpi. Each approximate center is
surrounded by a deliberately broad +/-5 percentage-point reading band; the
band is **not** an experimental confidence interval or recovered replicate
SD. Figure 1A's caption says 0.1%, 0.5%, 2% and 7% feeds had triplicates,
while the displayed SD for the other points was set to the highest value.
The figure's HCN and other product values were not transcribed into this
card. `backend/scripts/evaluate_dry_conversion_plot_bands.py` computes the
seven assay-free model cases, freezes their hash, then reads the held-out
card through the evaluation-only module:

| Feed CH4 | Plot-read CH4 conversion center (band), % | 200 Td model CH4 conversion, % | Inside visual band? |
| ---: | ---: | ---: | :---: |
| 0.1% | 80 (75–85) | 17.14 | No |
| 0.25% | 68 (63–73) | 17.15 | No |
| 0.5% | 58 (53–63) | 17.17 | No |
| 1% | 44 (39–49) | 17.23 | No |
| 2% | 34 (29–39) | 17.34 | No |
| 3.5% | 28 (23–33) | 17.43 | No |
| 7% | 18 (13–23) | 17.52 | Yes |

Six of seven model results miss even these wide reading bands, and the
model's feed trend points the wrong way. This is stronger evidence against
the **current uniform-exposure scenario**, not a like-for-like numeric model
error or an independent experimental uncertainty estimate. The same
unmeasured in-reactor temperature and local discharge exposure remain open.

The separate, frozen
[`Figure 1A ethane-selectivity plot-band card`](../backend/data/benchmarks/kamarinopoulou_2024_figure1a_ethane_selectivity_plot_bands_v1.json)
reads the teal C2H6 *carbon-selectivity* segment at the same seven feed
points from the 300-dpi source PDF, using interior bar pixels away from the
overlaid black energy curve. The centers are rounded to whole percentages
and assigned broad plot-reading bands, not measured confidence intervals.
`backend/scripts/evaluate_dry_ethane_plot_bands.py` freezes the same
assay-free predictions before reading this separate reference:

| Feed CH4 | Plot-read C2H6 carbon selectivity center (band), % | 200 Td model, % |
| ---: | ---: | ---: |
| 0.1% | 1 (0–6) | 53.32 |
| 0.25% | 1 (0–6) | 62.02 |
| 0.5% | 3 (0–8) | 68.37 |
| 1% | 8 (3–13) | 73.97 |
| 2% | 14 (9–19) | 78.50 |
| 3.5% | 19 (14–24) | 81.29 |
| 7% | 29 (24–34) | 83.58 |

The **increasing direction** still agrees, but all seven model values
miss the wide visual bands. The paper's GC/coke-based carbon inventory and
our all-species/no-wall-coke inventory do not match exactly, and the
conversion and discharge conditions also fail the matching gate. Thus this
is a strong conditional magnitude discrepancy for the current scenario,
not a validated quantitative selectivity-error ratio or evidence that one
specific missing pathway fixes it. No experimental HCN amount enters the
source or mechanism calculation.

The assay-free `backend/scripts/evaluate_dry_joint_ethane_yield_bands.py`
now freezes the same seven predictions, then loads the *separate* conversion
and ethane-selectivity plot cards and multiplies their nonnegative visual
band endpoints. This derives a broad Figure 1A **feed-carbon ethane-yield
screen** without transcribing an HCN value or treating the product as a new
measured assay. At 1% feed, the two plot readings imply `0.0117–0.0637`
of inlet carbon as C2H6; the model predicts `0.12746`, more than twice the
upper visual reading. All seven model yields lie above their derived upper
bands. Thus changing only the *selectivity denominator* cannot reconcile
the current model's ethane amount with these two non-HCN readings. The
bands are not experimental confidence intervals, the model and GC/coke
bases are not fully matched, and the reactor's temperature and discharge
exposure remain unknown: this is a conditional failure of the current
scenario, not a calibrated quantitative prediction error.

The separately frozen
[`Figure 1B hydrogen-selectivity plot-band card`](../backend/data/benchmarks/kamarinopoulou_2024_figure1b_hydrogen_selectivity_plot_bands_v1.json)
reads the pale-pink H2 *hydrogen-selectivity* segment at the same seven
feed points. As above, these are broad visual reading bands, not experimental
confidence intervals. `backend/scripts/evaluate_dry_hydrogen_plot_bands.py`
evaluates predictions frozen before the reference is read:

| Feed CH4 | Plot-read H2 hydrogen selectivity center (band), % | 200 Td model, % |
| ---: | ---: | ---: |
| 0.1% | 64 (59–69) | 35.93 |
| 0.25% | 57 (52–62) | 33.48 |
| 0.5% | 55 (50–60) | 31.69 |
| 1% | 55 (50–60) | 30.11 |
| 2% | 50 (45–55) | 28.83 |
| 3.5% | 45 (40–50) | 28.04 |
| 7% | 45 (40–50) | 27.36 |

All seven model values miss the wide visual bands, although the decreasing
direction agrees. The model's hydrogen accounting and the paper's reported
selectivity basis are not yet established as exactly like-for-like, and
discharge exposure and temperature remain unresolved. This is a conditional
failure of the current scenario, not a calibrated experimental error bar.
The paper's supplementary Eq. S10 also has a printed-symbol ambiguity:
its hydrogen-selectivity numerator shows `N_CX`, while the accompanying
text defines `N_HX` as the number of hydrogen atoms and Figure 1B includes
H2 (which has zero carbon atoms). The model uses hydrogen-atom counts, the
chemically coherent reading, but this is one more reason not to claim an
exact paper-basis match without the authors' calculation details.

To avoid changing the mechanism on guesswork, the source-versus-neutral diagnostic is
now recorded directly in the packet and feed-sweep output.
For each gas packet, net CH4 consumption equals the integrated CH4-target
electron-fragment source consumption **plus** net neutral-chemistry CH4
consumption; a negative neutral term means net methane regeneration:

| Feed CH4 | Net conversion / feed CH4 | Electron-fragment CH4 loss / feed CH4 | Net neutral CH4 loss / feed CH4 | Resolved-fragment energy (J of 7.093 J declared) |
| ---: | ---: | ---: | ---: | ---: |
| 0.1% | 0.171389 | 0.213986 | -0.042597 | 0.0060 |
| 1% | 0.172315 | 0.213389 | -0.041074 | 0.0597 |
| 7% | 0.175243 | 0.208508 | -0.033265 | 0.4083 |

The power-capped effective electron density remains near `2.27–2.31e14
m^-3` across these inputs, and the model's **per-CH4 source loss remains
about 21%** even as feed concentration rises seventyfold. This arithmetic
explains the simulated flat conversion: the neutral network regenerates
roughly 3–4% of inlet CH4, while the electron-fragment source scarcely
changes as a feed fraction. It does **not** establish why the experiment's
conversion falls; possible changes in true streamer exposure, local EEDF,
or omitted chemistry still require independent evidence. Resolved-fragment
energy is only one allocation within the electrical budget, not a complete
plasma or thermal energy closure. No experimental conversion value or HCN
amount was used to obtain these internal source diagnostics.
An optional, read-only integrated GRI-3.0 reaction-rate trace now resolves
that neutral term across all 16 packet substeps. At 0.1%, 1%, and 7% CH4,
the largest single methane-regenerating contribution is the reversible
`CH3 + H (+M) <=> CH4 (+M)` reaction; leading methane-consuming contributions
include `CH + CH4 <=> C2H4 + H` and `CH2(S) + CH4 <=> 2 CH3`.
For the 1% case, the integrated recombination contribution is
`+2.72392e-8 mol CH4`, while net neutral chemistry is
`+1.41559e-8 mol CH4` after the other reactions. The reaction-integral
versus observed neutral-chemistry balance differs by less than `1e-15 mol`
in that case; regression checks require less than `1e-12 mol` at all three
feeds and confirm that enabling the trace does not change conversion.
These are rates **within the chosen mechanism and uniform-packet scenario**,
not evidence that the reaction is experimentally dominant or that its
coefficient should be changed. The optional trace is exposed by
`trace_methane_neutral_pathways=True` in the dry-flow probe and feed-sweep
audit; it is off for ordinary runs.

An [independent CH4/N2 DBD experiment/model study by Snoeckx et al.](https://medialibrary.uantwerpen.be/oldcontent/container2642/files/ijhe13influence.pdf)
finds N2(A) and N2(a') metastables important in its *different* H2-producing
reactor and uses time-integrated reaction-rate analysis to explain methane
loss and regeneration. This is a specific research lead, not a transferable
rate law or validation of our reactor. Our EEDF includes N2 electronic
excitation as collision energy demand, and a separate Pavan **shadow** can
map some bands to excited species, but the primary GRI dry-flow inventory
does not propagate those N2 excited states or their CH4 quenching chemistry.
The dry-flow probe now reports the **uncoupled, read-only** electron-impact
opportunities by N2 excitation-band label. At the default 1% CH4 case, its
three N2(A3) bands total about `7.65e-7 mol` of opportunities and the
N2(a'1) band about `2.52e-7 mol` across the modeled packet. These are
neither accepted state populations nor product yields: multiple collisions
per molecule, incomplete band-to-state assignments, hypothetical whole-
packet exposure, and shared energy allocation prevent that interpretation.
The assay-free seven-feed sweep now exposes these uncoupled counts **per
inlet CH4** without accepting them as reactions. N2(A3) excitation
opportunities decline from `22.42` at 0.1% methane to `2.22` at 1% and
`0.296` at 7%; N2(a') opportunities decline from `7.40` to `0.732` and
`0.0974`. This is an *availability screen*: a feed-dependent excited-state
route is arithmetically possible, but the counts do not show that any such
route reacts efficiently. It is not
an upper bound on methane conversion or a product yield: the same N2 can
be excited repeatedly, state populations are unpropagated, and reactive
quenching versus back-transfer is unresolved. In particular, these ratios
must not be multiplied by an assumed CH4 branching fraction to fit the
Figure 1 conversion curve.
Indeed, [Mignogna et al. (2022)](https://doi.org/10.1088/1361-6595/ac99fc)
measured time-resolved N2(A3,v=0,1) populations in N2-CH4 nanosecond-pulse
plasmas and found that the initial reactant CH4 quenches those low-v
electronic states slowly, while CH4 rapidly relaxes their vibration;
accumulated fragments change later decay. Their 150-Torr diffuse reactor
does not determine rates at this approximately 1-atm DBD, nor does a
low-v observation resolve our broad `v0-4`, `v5-9`, and `v10-` excitation
bands. It does, however, argue against treating every N2(A3) opportunity
as a methane-breaking collision. Singlet and higher-state fate, afterglow
accumulation, and CH4 versus fragment competition must be represented
before this lead becomes a source term.
The next chemistry step requires reviewed, state-resolved excitation,
quenching, and afterglow evidence appropriate to this reactor before such
states can be added to the propagating solver.
A bounded rate check found a
[NIST-evaluated N2(A)+CH4 total-removal entry](https://kinetics.nist.gov/kinetics/Detail?id=1999HER1453-1483%3A75),
but it labels the outcome only “Products” and its recommended temperature
range is 300–360 K. It cannot supply a product-resolved 555 K quenching
channel. Snoeckx et al. assign such channels in their model, but that is not
an independent hot-gas product measurement for this reactor. No new
metastable reaction was inserted into the primary solver on that basis.
[Sharipov et al. (2016)](https://doi.org/10.1021/acs.jpca.6b04244)
provide a theoretical reason not to assume every N2(A)+CH4 encounter is
irreversible CH4 destruction: their low-barrier entrance produces
N2H+CH3, but its subsequent chemistry can return N2(X)+CH4 as well as
quench dissociatively. That study does not give a measured, product-resolved
branch at this reactor's temperature and pressure. A single ad hoc
`N2(A)+CH4 -> N2+CH3+H` term could therefore misstate both methane
conversion and hydrogen selectivity; the channel remains a research lead,
not a promoted source term.
There is a more specific independent measurement for the N2(a', v=0)
channel: [Umemoto et al. (2002)](https://doi.org/10.1063/1.1502642)
report a CH4 quenching coefficient of `(2.9 +/- 0.2)e-10 cm3 molecule-1 s-1`
and an H-atom production yield of `0.7 +/- 0.2` under their stated
assumption about H2-producing exits. This constrains a prepared singlet
state's removal and H-atom fate, **not** the full carbon-product branching,
the mixed vibrational-state population, or its kinetics at the dry probe's
uniform 555 K scenario. It is a candidate for a temperature-bounded
conditional sensitivity, not justification to assert CH3+H as the unique
propagating outcome. The earlier private
[N2+ recombination intake](HCN_FORMATION_GATE.md#private-lxcat-n2-dissociative-recombination-intake-2026-09-27)
likewise remains a partial, read-only diagnostic: its cross-section export
ends below 0.4 eV and cannot provide a complete recombination rate or
state-resolved N atoms for this spark. No ion or metastable shortcut is
promoted into the primary HCN prediction.
The paper's [supplementary Figure S6](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10980260/supplementaryFiles)
plots average-field estimates roughly around 180–195 Td, average power near
28 W, and similar steady temperature trajectories across the same feed
range. These are chart-level observations, not precisely digitized values or
resolved streamer measurements. They do not support explaining the sharp
experimental conversion decline by a large fall in the **reported average**
field or input power. An assay-free sensitivity sweep at a fixed 189 Td,
near that reported average, still predicts CH4 conversion rising with feed:
`0.166643 -> 0.167682 -> 0.170992` at 0.1%, 1%, and 7% CH4. The
independently held-out conversion-direction check therefore still fails;
the C2H6 and H2 selectivity directions still pass, and material and
collision-energy ledgers still close. The paper's estimated average field
does not resolve the discrepancy, but this test does not establish the
unknown local streamer field, exposure, or chemistry.

### Held-out residence-time check

The same [article's Figure 1D](https://rcastoragev2.blob.core.windows.net/ea712d57e86c0bac7a14d5ec18c0b746/sciadv.adl4246.PMC10980260.pdf)
varies nominal residence time at 1% feed CH4 and 28.1 W. The new
assay-free `backend/scripts/audit_dry_residence_sweep.py` varies the inferred
volumetric flow while holding the reactor volume, assumed 298.15 K sccm
reference, uniform 555.15 K gas scenario and power fixed. Its default field
is the approximately 184 Td **average estimate** in supplementary Figure
S7B, not a resolved streamer field:

| Nominal residence (s) | Inferred flow (sccm) | Model CH4 conversion |
| ---: | ---: | ---: |
| 1.88 | 50.106 | 51.86% |
| 0.94 | 100.213 | 30.66% |
| 0.63 | 149.524 | 21.64% |
| 0.47 | 200.426 | 16.52% |

The frozen, external check finds the **correct declining direction** as
residence shortens, but fails two deliberately broad visual Figure 1D
thresholds: its mean conversion symbols at 1.88 and 0.47 s are visibly
above 70% and 35%, respectively. These are plot-reading bounds on the
displayed means, **not** exact data, replicate intervals, digitized error
bars, or quantitative validation. All four model cases close material and
collision-energy upper-bound ledgers, exclude unresolved N2 product-state
sources, and keep HCN below the model's `1e-12 mol` reporting floor. The
flow-reference and thermal profile remain unverified, and the model does
not resolve streamer pulses. Thus the residence sweep exposes a second
non-HCN magnitude failure rather than validating discharge exposure.

A separate frozen
[`Figure 1D conversion plot-band card`](../backend/data/benchmarks/kamarinopoulou_2024_figure1d_conversion_plot_bands_v1.json)
extends this check to all four gray-curve points. As for Figure 1A, these
are pixel-assisted 300-dpi visual centers with intentionally broad +/-5 percentage-
point reading bands, **not** experimental confidence intervals or extracted
replicate SD. The caption reports triplicates at all four residence
conditions. No HCN values are transcribed. The
`backend/scripts/evaluate_dry_residence_plot_bands.py` runner freezes the
assay-free residence predictions before the evaluation-only module reads
the card:

| Nominal residence (s) | Plot-read CH4 conversion center (band), % | 184 Td model CH4 conversion, % | Inside visual band? |
| ---: | ---: | ---: | :---: |
| 1.88 | 82 (77–87) | 51.86 | No |
| 0.94 | 70 (65–75) | 30.66 | No |
| 0.63 | 56 (51–61) | 21.64 | No |
| 0.47 | 46 (41–51) | 16.52 | No |

All four miss the broad bands despite the correctly declining direction.
This is an external failure of the **current uniform-field/exposure model**,
not an exact measurement-error ratio or proof of which omitted physical
process is responsible. Reactor temperature, actual flow reference,
streamer encounters, and nitrogen-product states remain unresolved.

The same run's read-only plasma diagnostic integrates `4.46748e-7 mol` of
**uncoupled total N2 outcome opportunities** at the hypothetical power-capped
electron density, with `0.420701 J` of overlapping threshold-energy estimate,
against `7.09298 J` declared electricity per gas packet. The inlet methane is
`3.44645e-7 mol`. This shows that the *calculated total collision opportunity*
is not obviously negligible in the assumed whole-packet field; it does **not**
establish actual N2 fragmentation, N(4S)/N(2D) branching, exposure, or HCN
production. It is not injected into chemistry or double-charged to the energy
ledger. Scaling it by the paper's product-fitted streamer fraction would be
invalid as an independent prediction.
Within that total, `5.00227e-8 mol` lies below 14.526 eV and is therefore a
**conditional ground-N(4S) opportunity**: for a genuine neutral dissociation
event of ground-state N2 below the first two-excited-atom threshold, at least
one fragment must be N(4S). This does not turn the modeled opportunity into
an observed event or validated atom amount. The near-threshold total cross
section lies below Cosby's state-measurement window, the partner fragment and
actual discharge exposure remain unresolved, and no atom source is injected.

The assay-free `backend/scripts/audit_dry_field_sweep.py` repeats the same
16-step packet at four declared field assumptions. The paper's 189 and 435 Td
are average-field estimates, 300 Td is its product-trend-selected model
setting, and 200 Td is our hypothetical base point; none resolves local
streamer field or exposure:

| E/N (Td) | CH4 conversion | C2H6 carbon yield | HCN carbon yield, all below reporting floor | Uncoupled N2 opportunity (mol) |
| ---: | ---: | ---: | ---: | ---: |
| 189 | 0.167682 | 0.123490 | 7.37228e-10 | 3.92249e-7 |
| 200 | 0.172315 | 0.127455 | 7.59430e-10 | 4.46748e-7 |
| 300 | 0.185954 | 0.142476 | 8.33350e-10 | 9.36387e-7 |
| 435 | 0.174595 | 0.138371 | 8.01220e-10 | 1.48617e-6 |

The carbon-product ordering failure persists across this field range. Each
run closes the material ledger and collision-energy upper bound; none closes
the full thermal balance or qualifies N2 fragment states. The figures are
scenario sensitivities, **not** measured yield estimates.

The independently controlled
`backend/scripts/audit_dry_temperature_sweep.py` varies a *uniform* gas
temperature while retaining the declared ambient flow reference and 200 Td
field. Only 555.15 K equals the paper's steady **downstream probe** reading;
450 and 700 K are scenarios, not measured alternatives or a reconstruction
of the reactor profile:

| Uniform scenario T (K) | Ideal-gas packet residence (s) | CH4 conversion | C2H6 carbon yield | HCN carbon yield, below reporting floor |
| ---: | ---: | ---: | ---: | ---: |
| 450 | 0.311401 | 0.170203 | 0.109303 | 3.91765e-12 |
| 555.15 | 0.252419 | 0.172315 | 0.127455 | 7.59430e-10 |
| 700 | 0.200186 | 0.198844 | 0.158874 | 2.75886e-8 |

All runs close the material and collision-energy upper-bound ledgers and keep
the unresolved N2 state source excluded. Temperature changes the nominal
time and hydrocarbon distribution, but does not reverse the HCN/C2H6
ordering. The paper's 0.47 s residence is the `V/Q` nominal value at
200 sccm. Its [supplementary Eq. S21](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10980260/supplementaryFiles)
prints `t = 60 V / Q`, with `V = 1.57 cm3` and the controller flow in
`cm3/min` for this calculation: `60 * 1.57 / 200 = 0.471 s`.
The same letter `Q` is described in `L/s` for the separate energy formula;
those units must not be mixed. Eq. S21 does not measure in-reactor gas
temperature or transit time. Converting a 555.15 K uniform inventory to the assumed 298.15 K
flow reference gives 0.2524 s. If, **additionally**, every point in the
reactor were between those two temperatures at fixed pressure, the ideal-gas
residence would lie between 0.2524 and 0.47 s. The downstream probe reading
does not establish that profile bound, so no actual residence time is
claimed from it.

### Mass-flow-controller reference convention

The article identifies **Brooks GF40** mass-flow controllers but does not
state their installed reference-temperature code. The manufacturer's
[GF40 installation manual](https://cdn.brooksinstrument.com/-/media/brooks/documentation/products/mass-flow-controllers/elastomer-sealed/gf40/installation-manual-gf40.pdf)
(2023, Table 5-1) lists 0, 15, 20 C, and 70 F reference options; its listed
standard-configuration bins use 0 C, but that does **not** establish which
code was used in this experiment. The assay-free
`backend/scripts/audit_dry_flow_reference.py` compares those listed options
with the project's earlier, separately labeled **assumed 25 C** reference,
holding the uniform 555.15 K reactor scenario, nominal 200 sccm and 28.1 W
fixed:

| Assumed flow reference | Hot-packet residence (s) | Declared J per packet | CH4 conversion | C2H6 carbon selectivity |
| ---: | ---: | ---: | ---: | ---: |
| 0 C, listed option | 0.231254 | 6.49823 | 0.158720 | 0.731578 |
| 15 C, listed option | 0.243953 | 6.85508 | 0.166902 | 0.736557 |
| 20 C, listed option | 0.248186 | 6.97403 | 0.169613 | 0.738130 |
| 70 F, listed option | 0.249127 | 7.00046 | 0.170214 | 0.738474 |
| 25 C, legacy assumption | 0.252419 | 7.09298 | 0.172315 | 0.739663 |

All runs retain the unresolved N2 state-source exclusion, close the material
and collision-energy upper-bound ledgers, and keep HCN below the model's
1e-12 mol reporting floor. They do not identify the installed controller
configuration. Choosing the option that best matches an HCN yield would be
product calibration, not an independent flow measurement. The remaining
temperature profile and discharge-exposure uncertainties are larger physical
questions than this reference-option sweep, but the convention must still be
fixed before claiming a like-for-like residence time or energy cost.

### Non-product electrical exposure audit

The paper's [supplementary electrical method](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10980260/supplementaryFiles)
(eqs. S23–S26) prints an emission-derived field of 435 Td, half-cycle mean
current `2.66 mA`, annular area `0.314 cm2`, current density `407 A/m2`,
BOLSIG+ drift speed `3.46e5 m/s`, and inferred electron density
`7.3e9 cm-3`. The new assay-free
`backend/scripts/audit_dry_electrical_exposure.py` checks these independently:

- The diameters give `0.314159 cm2`, agreeing with the printed area, but
  `2.66 mA / 0.314159 cm2 = 84.67 A/m2`, not `407 A/m2` (factor `4.81`).
  The latter implies `7.34e9 cm-3` with the printed drift speed; the
  current-and-area arithmetic instead implies `1.53e9 cm-3`. This may be a
  reporting/definition issue; the printed values cannot all be used
  simultaneously as one homogeneous-current measurement.
- Our separately solved 1% CH4/N2 EEDF at 435 Td and 550 K gives drift
  speed `3.347e5 m/s`, within about 3.3% of the paper's quoted value.
  The two calculations share Phelps N2 collision ancestry, so this is a
  transport cross-check, not fully independent confirmation; it does not
  validate local streamer density or product chemistry.
- At the paper's printed electron density, a **uniform, continuously active**
  435 Td field in the whole annulus would require about `3568 W` of electron
  drift power; using the current-and-area-recomputed density gives `746 W`.
  With `28.1 W` average electrical input, the respective *conditional*
  upper limits on a joint space-time filling fraction are `0.00788` and
  `0.0376` if those densities represent local streamer values and **all**
  measured power drives electron drift. If only some power does so, the
  allowed fractions are smaller.

These are **not** measured streamer occupancies. The paper estimated density
by assuming a homogeneous annular area; AC current may include non-electron
components, while true local field, current and duty cycle are unresolved.
The paper's product-fitted 0.3% streamer-volume estimate is not used as an
input or validated by these different conditional power bounds. The audit
does not change any fragment source or fit HCN. Before using electrical
estimates to drive spatially resolved chemistry, the current-density
definition and local-versus-averaged electron density must be reconciled.
A uniform power-capped electron density can match the *linear* total
electron-target source rate of a localized density times a filling fraction;
it cannot match nonlinear local radical chemistry, gradients, or repeated
streamer encounters. Entering either printed or recomputed electron density
into the existing whole-packet model simply activates the same power cap,
so it would not resolve the physical exposure problem.

The same supplement's **microdischarge-count estimate (S27)** has a further
arithmetic/definition ambiguity, checked by the assay-free electrical audit.
It states a product-fitted `0.3%` plasma volume (`4.2 mm3`), a *typical*
filament **diameter** of `0.1 mm`, length `2 mm`, and about `67`
microdischarges per cycle. As typeset, S27 reads
`N_MD = pi d_f^2 l_f / V_plasma`, which evaluates to `0.01496`, not 67.
Inverting that ratio to count cylindrical volumes gives `267.38` if the
stated `0.1 mm` is a **diameter**; it gives `66.85` only if that value is
treated as a **radius**. The discrepancy could reflect notation or
typesetting, but neither interpretation is an independent filament count:
the 4.2 mm3 plasma volume was derived from a product-fitted reactor fraction
and the diameter was a literature-typical assumption, not a direct
reactor-specific measurement. Neither `67`, `267`, nor the `0.3%` fraction
is used to select encounter count or calibrate chemistry.

### Independent optical field-ratio check fails with the installed collision set

The paper's [supplementary Eq. S24 and Fig. S12A](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10980260/supplementaryFiles)
report a measured `N2+(B-X,0-0) / N2(C-B,0-0)` emission-peak-area ratio
of `3.47e-3` from a **separate 28.6 W, 1% CH4** optical run, and infer
`435 Td` with their Phelps/BOLSIG+ workflow. This is a non-product
diagnostic, so it can constrain field without importing HCN yields.
`backend/scripts/audit_dry_optical_field.py` independently applies the
published radiative rates and N2 quenching constants to the installed
Cantera/Phelps `N2+(B)` and `N2(C)` electron-impact channels and our
two-term CH4/N2 EEDF at the paper's spectroscopic `550 K, 1 atm` setting:

| Field (Td) | Computed emission ratio | Computed / observed |
| ---: | ---: | ---: |
| 100 | 2.15e-6 | 0.000619 |
| 189 | 4.09e-5 | 0.0118 |
| 300 | 1.98e-4 | 0.0570 |
| 435 | 5.97e-4 | 0.172 |
| 500 | 8.84e-4 | 0.255 |

At `435 Td` the computed ratio is **5.81-fold below** the reported optical
ratio. The `60/100/150 eV` EEDF-grid results at that field are
`5.927e-4 / 5.969e-4 / 5.969e-4`, so grid truncation is not the primary
cause. The published quench/radiative constants give a N2+(B)-to-N2(C)
survival-factor ratio of `0.01081` in the stated gas. Our EEDF direct
excitation-rate ratio is `0.05522`; reproducing the measured optical ratio
under that same emission formula would require `0.3210`. These are
diagnostic decompositions, **not** permission to rescale a collision channel.
This is a failed reproduction of the *diagnostic inference*, not
evidence that the actual field exceeds 500 Td: the installed
[Cantera Phelps file](https://github.com/Cantera/cantera-example-data/blob/main/air-plasma-Phelps.yaml)
states an **August 2025** LXCat retrieval, whereas the paper cites
**March 2023**, so these are not verified identical curve versions; spectral response,
quenching, spatial averaging and two-term EEDF assumptions are not
independently validated here. The paper's optical run also differs from
Figure 1's 28.1 W nominal setting. The supplementary optical-method section
names the spectrometer and analysis software but does not state a relative
337/391 nm spectral-response correction; this is an unresolved measurement
definition, **not** evidence that a correction would remove the mismatch.
The [spectrometer manufacturer's calibration guidance](https://www.avantes.com/support/services-and-calibrations/)
specifies wavelength-dependent irradiance calibration for the assembled
spectrometer/fiber/collection-optic system; a generic detector curve cannot
substitute for the experimental instrument's calibration record.
An [independent atmospheric-pressure DBD optical-diagnostic study](https://d-nb.info/1261240057/34)
also identifies calibration, cross sections, quenching and model assumptions
as uncertainty sources for this nitrogen-band-ratio method. Its different
N2/O2 discharges provide **no numerical correction** for this CH4/N2 reactor.
No field was retuned, and no HCN or product assay was loaded. Reconcile
the cross-section and emission
definitions before treating `435 Td` as a verified model input.
The date difference alone does **not** establish a changed N2(C) or
N2+(B) curve. In [Phelps' revised ELECTRON.TXT note](https://nl.lxcat.net/notes/index.php?download=phelps5)
(last revised 2008), the listed N2(C) cross section is `0.2500` at 20 eV
and `0.1110` at 30 eV, and the N2+(B) cross section is `0.0720` at 30 eV
and `0.1470` at 50 eV (table units of 1e-20 m2). The corresponding installed
Cantera/Phelps channel values match at these sampled energies. The note says
the tabulated N2(C) curve already includes the historical 0.67 correction;
its `QSCALE` is for recovering older input data, **not** a second multiplier
for these tabulated values. This limited spot check weakens the simple
"export vintage changed the two channels" explanation, but does not compare
the complete 2023 and 2025 sets or validate the emission model. A matching
March 2023 Phelps export has not been verified; no restricted LXCat
cross-section data are copied into this repository.
The same note defines its N2 momentum-transfer curve as *effective*, including
inelastic contributions. BOLOS reconstructs the elastic curve by subtracting
the inelastic channels, and the installed N2/CH4 set reports no negative
elastic clipping across the audit's 100–500 Td field grid (covered by a
regression test); treating the effective curve as an
additional elastic process would be a double count, not a justified repair.
Using a **pure-N2** rather than a 1% CH4/99% N2 EEDF at `435 Td` changes
the computed ratio only from `5.969e-4` to `6.120e-4` (about `2.5%`);
the mixture-choice ambiguity does not close the observed factor-of-5.8 gap.
Changing the audit's uniform gas temperature from the paper's `550 K` to
`298.15 K` yields `5.725e-4` (about `4.1%` lower); this hypothetical
room-temperature limit likewise does not close the gap, and is **not** a
claim that the actual plasma gas was at room temperature.
The independent **Cantera two-term EEDF**, folded through the *same*
collision curves at `435 Td`, predicts `4.64e-4` (about `7.47` times below
the optical observation). Its N2+(B) and N2(C) rates differ from BOLOS by
factors `1.34` and `1.04`, with a `4.73%` mean-energy difference. Thus an
ordinary swap between these two solvers does not close the optical gap;
the cross-section vintage, optical model, and spatially averaged measurement
remain unresolved possibilities. This is not an
independent cross-section validation, because both solvers use the same file.
The next falsifiable comparison is the paper's *two direct excitation-rate
coefficients* at 435 Td (or its original BOLSIG+ input/output) and the
relative 337/391 nm instrument calibration. Those separate an EEDF/collision
discrepancy from an emission-measurement discrepancy without using HCN
product data to choose the field.

### Conditional two-zone localization test: encounter count not identified

`backend/scripts/audit_dry_localized_exposure.py` tests one explicit
limiting interpretation of those electrical bounds: during each chemistry
step, only a selected gas-volume fraction sees 435 Td at one of the two
current-derived local-density interpretations; the rest sees neutral gas
chemistry, and both zones are completely remixed afterward. Each local
source is constrained by the same packet electrical energy. **Equating the
power-limited joint space-time fraction to a continuously active volume is
an extra, unmeasured assumption**, not a recovery of the actual streamer
geometry. The 435 Td is an emission-derived estimate, not a resolved local
filament field. The two fraction/density pairs are linked by the respective
printed-current interpretations, not chosen against products:

| Complete remix encounters per packet | 0.00788 active fraction, printed-density interpretation: model CH4 conversion | 0.03765 active fraction, current/area-recomputed density: model CH4 conversion |
| ---: | ---: | ---: |
| 16 | 6.20% | 14.69% |
| 32 | 11.93% | 15.47% |
| 64 | 13.97% | 15.86% |
| 128 | 14.90% | 16.08% |

For reference, active fraction one reproduces the existing uniform 435 Td
packet calculation **exactly** at the same declared density and steps
(`17.46%` conversion at 16 steps). The localized values are not a
time-step convergence series: increasing “substeps” also increases the
assumed number of complete gas-remixing encounters, an unknown **physical**
parameter. The large change in the narrow-fraction series proves that
electrical power alone does not identify a unique outlet conversion under
this localized representation. Neither zone resolves nanosecond pulses,
actual filament paths, gas exchange rates, spatial temperature, or N2
product states. All checked cases close atom and collision-energy ledgers;
spark-off, carbon-free and nitrogen-free controls remain below the HCN
reporting floor. The result is therefore a negative **identifiability**
finding, not a validated localization correction, a numeric yield, or a
reason to choose a fraction that matches Figure 1. The primary predictive
solver remains unchanged.

An independent **molecular-diffusion scale** now sharpens the mixing
question. `backend/scripts/audit_dry_diffusive_mixing.py` uses the installed
Cantera 3.2 GRI-3.0 [mixture-averaged diffusion coefficient](https://cantera.org/stable/python/transport.html)
for 1% CH4 in N2 at the declared uniform `555.15 K`, `1 atm` scenario:
`D_CH4,mix = 6.61e-5 m2/s`. The paper's 7 mm quartz inner diameter and 3 mm
rod diameter imply a 2 mm radial annular gap. The one-dimensional RMS
diffusive distance `sqrt(2Dt)` is about `78 um` per 22 kHz AC period,
`1.44 mm` per 16-encounter interval, and `0.51 mm` per 128-encounter
interval. A 2 mm RMS displacement takes about `30 ms` on this molecular
diffusion estimate, roughly `666` AC periods. Across the whole assumed
`0.252 s` packet residence, the unbounded one-dimensional RMS length is
about `5.78 mm`, larger than the gap; this shows why slow bulk mixing and
very short-lived filament gradients must be distinguished. Molecular
diffusion **alone** does not establish
complete gap-scale remixing each period or each modeled encounter. These
RMS scales are not first-passage probabilities or complete-mixing times;
convection, thermal gradients, moving filaments and electrohydrodynamics
could change actual mixing. Their rates are not measured here. This supports
keeping encounter frequency as an unresolved physical parameter rather than
selecting one of the two-zone conversion values as a prediction.

### Independent pure-nitrogen atom-cost screen, not an HCN yield

A bounded [atmospheric-pressure N2 DBD analysis by Popov (2013)](https://doi.org/10.1134/S1063780X13050085)
provides a different check on the **total N2 collision source**, not a
state-branch law. It used a planar, homogeneous 2 mm Townsend discharge in
pure N2 near 300 K. From an earlier TALIF measurement of ground N(4S)
density, electrical power and modeled bulk/wall atom losses, it inferred
`337 +/- 35 eV` of discharge energy per N(4S) atom at an effective field
near `140 Td`. That is not a direct measurement of fragment branching, and
the neutral-dissociation cross-section comparison shares Cosby ancestry
with our imported Tabata total outcome. It is not the coaxial, filamentary
1% CH4/N2 reactor or its downstream 555.15 K temperature proxy.

The read-only `backend/scripts/audit_dry_nitrogen_atom_energy_cost.py`
computes an **intrinsic electron-drift energy per total N2 outcome** from
the current EEDF without specifying electron density, reactor occupancy or
residence. With reduced mobility `mu*N`, reduced field `E/N`, total-outcome
rate coefficient `k`, and N2 mole fraction `x_N2`, this is
`(mu*N)(E/N)^2/(k*x_N2)` eV per event. Pure N2 at 300 K/140 Td gives
`430.11 eV` per total outcome, or an optimistic `215.06 eV` per ground atom
**if both fragments were N(4S)**; if exactly one ground atom survives per
event, the conditional electron-drift cost is `430.11 eV` per ground atom.
The 60/100/150 eV EEDF grids give `430.11/430.66/430.11 eV` per outcome,
so this particular integral is numerically stable to about `0.13%`.

Popov's `337 +/- 35 eV/atom` lies between those two idealized numbers, but
that is **not** a fit or validation: the real electronic branching, atom
survival, wall losses and discharge-versus-electron power partition differ,
and the reference reactor is different. For the project's 1% CH4/N2,
uniform 555.15 K, **assumed average** 189 Td scenario, the corresponding
intrinsic cost is `187.74 eV` per total outcome, or `93.87 eV` per atom only
under the deliberately optimistic two-ground-atom assumption. These are
source-efficiency screens, **not** actual atom costs, HCN yield estimates or
permission to inject N2 fragments. The audit loads no HCN assay and leaves
the propagating solver unchanged. The comparison supports retaining the
total collision channel as a useful plausibility diagnostic while leaving
state-resolved products and actual streamer exposure unresolved.

### External N2 electron-swarm transport residual, not reactor calibration

`backend/scripts/audit_dry_n2_swarm_transport.py` compares the imported N2
collision set and two-term EEDF mobility against a small, separately
transcribed subset of [Dutton's NIST electron-swarm survey](https://srd.nist.gov/JPCRD/jpcrd69.pdf),
printed page 621, table 1.17. The table identifies the original Wagner,
Frommhold, Blevin and Wagner/Raether measurements. Its `E/N` unit is one Td
and its drift-speed unit `10^7 cm/s` is `10^5 m/s`. The audit solves **pure
N2 at a declared 295 K**, converts reduced mobility by
`v_d = (mu*N)(E/N)`, and never passes the measurements, HCN assay, or a
product-selected field into the solver.

| E/N (Td) | Historical measured drift (km/s) | Current model drift (km/s) | Model residual |
| ---: | ---: | ---: | ---: |
| 164.2 | 158 | 158.7 | +0.5% |
| 208.0 | 196 | 191.6 | -2.2% |
| 212.2 | 224 | 194.7 | -13.1% |
| 242.5 | 253 | 216.2 | -14.5% |
| 243.0 | 223 | 216.6 | -2.9% |
| 330.0 | 330 | 273.9 | -17.0% |
| 333.0 | 310 | 275.7 | -11.1% |
| 446.0 | 450 | 341.8 | -24.0% |
| 450.0 | 450 | 344.0 | -23.6% |

The 450 Td model drift changes only from `343.48` to `343.99 km/s` across
60 eV/400-cell, 100 eV/700-cell, and 150 eV/1000-cell EEDF grids, so the
high-field residual is not explained by this discretization. This is an
**external residual check**, not an independent validation: Phelps cross
sections may have been fitted to some of these historical swarm studies,
the studies differ somewhat among themselves, their drift-speed definition
may differ from BOLOS mobility-derived speed in ionizing conditions, and their pure
N2 ambient gas is not the paper's hot 1% CH4/N2 flow. Thus neither the
individual residuals nor a correction based on them is transferable to a
435 Td streamer without more evidence. The audit does **not** revise the
paper's local-density or conditional space-time filling bounds above. It
shows that the shared-ancestry ~3% BOLSIG+ agreement does not establish
transport accuracy to a few percent; transport-model uncertainty remains
in the electrical-exposure gate. No N2 fragments or HCN were introduced.

### Flow energy accounting, still conditional

The dry probe now distinguishes fixed-temperature **internal-energy change**
from a steady-flow **inlet-to-outlet ideal-gas enthalpy change**. At the 200 Td
base scenario, the latter is `0.264569 J` per hot gas packet or `1.04814 W`
at the assumed molar throughput, compared with `28.1 W` declared electrical
input. The conditional `27.0519 W` difference is **not** a measured wall-loss
rate or a closed first-law balance. It depends on the assumed 298.15 K feed
reference, modeled outlet composition, uniform 555.15 K outlet proxy, and
omits electrode coupling, radiation, plasma excitation/storage, and actual
temperature gradients. The companion Cantera-unit test checks that ideal-gas
enthalpy minus internal energy equals `RT` for one mole of methane. This
separation prevents the older fixed-temperature `0.002828 J` internal-energy
change from being mistaken for the gas-flow heating requirement.

## Next smallest defensible gate

Resolve the **electronic product states and fate** of N2 dissociation under
the experiment's electron-energy distribution before injecting N atoms into
the primary solver. A total N2 dissociation cross section is not a measured
`N2 -> 2 N(4S)` branch. Investigate independent state-resolved data and
reaction/relaxation rates with stated temperature and pressure domains; if
no defensible branch exists, retain the exclusion and report the blocker.
After a defensible nitrogen source exists, the GRI HCN subnetwork also needs
independent scrutiny rather than a rate fit to HCN. The earlier wet-gate audit
identified GRI's direct `CH3 + N -> HCN + H2` and pressure-dependent
`H + HCN (+M) -> H2CN (+M)` as potentially influential. A
[primary quantum-kinetics study](https://arxiv.org/pdf/1902.05574) finds the
direct CH3+N channel inefficient and favors an H2CN-mediated interpretation;
its reported temperature dependence is mainly 50–400 K, not a verified
replacement law for this 555 K dry-flow scenario. Thus even a qualified
N-atom source alone would not establish a credible quantitative HCN yield.
In parallel, keep the paper's feed/flow/power/outlet data as external
benchmark records, with measurements separated from model inputs. Do not
replace the missing source by fitting its strength to HCN or by importing
the published HCN reaction recipe.

### Bounded state-fate intake, 2026-09-29

[Cosby's neutral N2 electron-impact experiment](https://doi.org/10.1063/1.464385)
identifies N(2D)+N(4S) as the primary product interpretation over its
18.5–148.5 eV measured fragment window; it does not supply a branching law
over the full calculated EEDF. A separate
[2023 nitrogen-collision review](https://doi.org/10.1140/epjd/s10053-023-00687-5)
reports about 20–30% uncertainty even for the *total* neutral-dissociation
cross section and says the ground-state-pair channel is not significant at
50 eV. It identifies neutral dissociation as especially under-studied; the
review does not provide the energy-dependent product-state cross sections
needed to turn our total collision outcome into a qualified atom source. A
separate
[Zipf and McLaughlin near-threshold study](https://www.sciencedirect.com/science/article/pii/0032063378900661)
attributes much of 11–24 eV N2 dissociation to singlet-state predissociation
and models the **total** cross section, but does not supply the branch law
needed here. A separate
[Ratliff et al. state-specific excitation/predissociation experiment](https://spacewx.net/wp-content/uploads/2020/07/P_HEL_PSSD_AMP_JGR_1991.pdf)
measured the N2 `b(1Pi_u)` excitation cross section at 60 and 100 eV and
inferred roughly 96% **molecular-state predissociation** at 100 eV. That
percentage is the chance that the excited `b` state dissociates instead of
emitting, **not** the fraction yielding N(4S)+N(2D), and it covers neither
the full low-energy EEDF nor all neutral-dissociation channels. It cannot
be reused as an atom-product branch or an HCN source multiplier. A separate
[NIST-evaluated N(2D)+CH4 entry](https://kinetics.nist.gov/kinetics/Detail?id=1999HER1453-1483%3A13)
reports a total loss-rate expression but labels the products only
“Products.” A [later direct rate experiment](https://pubmed.ncbi.nlm.nih.gov/30849164/)
measured N(2D)+CH4 loss from 127–296 K. Those entries do not identify a
550 K product branch for the paper's nominal 555 K outlet condition. The
[Umemoto et al. direct rate experiment](https://doi.org/10.1021/jp9833146)
covers 223–298 K and likewise does not establish a 555 K rate or branching.
A [product-resolved crossed-beam experiment](https://doi.org/10.1021/jp904302g)
does identify methanimine and methylnitrene from the H-displacement channel,
but at 22.2–65.1 kJ/mol collision energies. It reports energy-dependent,
nonstatistical branching; it is not a thermal 555 K total-rate or branching
measurement. Importing those beam ratios as the hot-flow branching law would
not be defensible.
The existing first-collision N(2D) audit is intentionally limited to near
298 K; the dry probe now surfaces its exact status as **unavailable: outside
reviewed near-298 K diagnostic domain**. Extrapolating a total-loss rate to
555 K would still not identify whether it makes NH, CH2NH, CH3N or another
product.
Accordingly, no N(2D)+CH4 reaction or N2 fragment yield was added to the
propagating solver. A defensible extension needs energy-dependent N2
branching and state-specific reaction outcomes over the actual hot-gas
temperature range, not just a total disappearance rate.
As a deliberately **room-temperature-only** control, the installed
[Herron (1999) evaluated rates](https://doi.org/10.1063/1.556043) for
thermalized N(2D) in a hypothetical 1% CH4/99% N2, 1-atm, 298.15 K
homogeneous bath give first-removal fractions `0.704` by CH4 and `0.296`
by N2, with a conditional lifetime `7.15e-7 s`. The evaluated multiplicative
rate-factor extremes put the CH4 fraction in `0.514–0.842`; these are not
confidence limits. Thus simple N2 quenching is **not automatically** the
dominant first fate even in a nitrogen-rich *ambient* bath. The same audit
explicitly returns **unavailable** at 550 K; none of these ambient fractions
are used as hot-reactor chemistry or as HCN product branching.
There is a [primary theoretical N(2D)+N2 quenching study](https://doi.org/10.1021/jz401012t)
covering `240–1000 K`, but an accessible numerical 550 K rate and a
comparably supported **N(2D)+CH4** hot-gas removal/product law were not
established in this bounded intake. Neither the ambient N2 fraction nor a
secondary compiled high-temperature rate was extrapolated into the solver.

### Existing CH + N2 route is not a validated HCN bypass

An N2 atom source is **not a logical prerequisite for every HCN pathway**:
the installed GRI-3.0 mechanism already includes both thermal
`CH + N2 <=> HCN + N` and pressure-dependent
`CH + N2 (+M) <=> HCNN (+M)`. The read-only
`backend/scripts/audit_dry_ch_n2_routes.py` evaluates their **forward rate
constants**, without using the paper's HCN measurements. At the declared
uniform `555.15 K`, `1 atm`, `1% CH4/99% N2` gas state, the installed GRI
values are `9.6533` and `2.37461e8 m3 kmol-1 s-1`, respectively; direct
HCN formation divided by association is only `4.07e-8`. That ratio is
**not a measured product branching fraction** or a final-yield fraction:
association can reverse or lead to other products, and the model does not
resolve the competing NCN/HNCN adduct network. The simulated direct-route
outlet net rate is only `4.05e-16 mol/s` in the 16-step evolving-EEDF
scenario; the outlet HCN remains `2.62e-16 mol`, below the project's
`1e-12 mol` reporting floor. An outlet instantaneous rate cannot by itself
establish the pathway's integrated contribution.

The primary [Brownsword et al. CH + N2 experiment](https://doi.org/10.1039/FT9969201087)
covers `202–584 K` and reports association-dominated loss, with the direct
HCN+N route too slow to matter in **its** conditions. Its argon-supported
measurements [extend only to about `0.53 bar`](https://kinetics.nist.gov/kinetics/Detail?id=1996BRO%2FHER1087-1094%3A2),
so they do not independently
validate the model's 1-atm nitrogen-bath falloff or a local streamer.
This audit therefore identifies a real alternative route to scrutinize,
not an excuse to bypass state-resolved nitrogen chemistry or to increase
the direct GRI rate until HCN appears. Independently constraining CH
production and the HCNN/NCN/HNCN fates is a parallel scientific gate.

The read-only `backend/scripts/audit_dry_n2_source_energy_support.py` now
applies the installed CH4/N2 collision set to the **dry 1% CH4/99% N2,
555.15 K scenario**. It partitions the EEDF-folded *total* N2 outcome by
electron energy; it never creates N atoms or reads HCN observations:

| Declared reduced field | Total N2 rate coefficient (m3/s) | Fraction below Cosby's 18.5 eV measurement start |
| ---: | ---: | ---: |
| 189 Td, paper average-field estimate | 1.80231e-16 | 59.17% |
| 200 Td, project hypothetical scenario | 2.27148e-16 | 55.88% |
| 300 Td, paper model-selected setting | 9.78913e-16 | 34.06% |
| 435 Td, paper average-field estimate | 2.97458e-15 | 19.10% |

At 200 Td the 60/100/150 eV grid checks give total coefficients of
`2.27148e-16`, `2.26577e-16`, and `2.27148e-16 m3/s` (less than 1% spread).
This supports numerical stability of the **total collision-rate integral**,
not its unknown product-state branch. The field dependence is large, and
even the paper's reported average-field estimates do not specify a local
streamer EEDF. In particular, neither the measured-window fraction nor the
conditional energy threshold below 14.526 eV may be used as a ground-state
N(4S) branching probability. The next steps are independent state-resolved
N2 evidence and fate kinetics **and** a non-product-fitted explanation for
the failed CH4 conversion-versus-feed direction and the conversion
**magnitude** failures across feed and residence. Reported average E/N and
power are insufficient to determine local streamer exposure, so the primary
solver stays unchanged pending a defensible constraint. The two-zone
sensitivity further shows that a power-derived space-time bound cannot
specify a gas-packet encounter count or complete-remix rate; independent
filament/flow diagnostics or a spatial transport model are needed.
Before claiming a more precise electrical occupancy from drift velocity,
also reconcile the external pure-N2 swarm residual and the measured-versus-
calculated drift-speed definitions under the relevant ionizing field and gas
conditions.

The older [HCN formation gate](HCN_FORMATION_GATE.md) preserves the wet
research history and archived dry probes. The separate
[reduced-source route](REDUCED_HCN_AND_ATOM_CEILING.md) remains a conditional
sensitivity tool, not the primary quantitative benchmark.
