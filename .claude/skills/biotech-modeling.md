# Bioprocess modeling judgment

## 1. When to use this skill

Use this for mechanistic modeling work on batch or fed-batch bioprocesses:
choosing species/reaction-network resolution for a kinetic fit, deciding
whether an existing fitted-parameter artifact is still the right one to
build on, scaling a process up or down between vessel sizes, and validating
a mechanistic model by running it forward against a real batch it wasn't fit
to. It complements, and should be read alongside, `parameter-estimation.md`
(how to actually run the fits and judge fit quality) and
`ode-stiff-kinetics.md` (how to keep the resulting ODE/PDE numerically
well-behaved) — this file is about the modeling judgment calls that sit
around those two: how much structure to claim, when to trust an artifact,
how to move between scales, and how to know if a model is actually right.

## 2. Species resolution: pooled vs. full regiochemistry

**The core question is not "how many species could I write balances for" —
it's "how many species does the analytical method actually distinguish, and
of those, how many are identity-confirmed."** A reaction network can be
written at full site-resolved (regiochemistry) resolution whenever the
chemistry supports multiple isomers per degree of substitution, but that
doesn't mean the data can identify separate rate constants for each isomer.

**Decision criteria, in order:**
1. **What's confirmed vs. inferred.** If only a subset of peaks are
   identified against a reference standard (mass spec, authenticated
   standard) and the rest are assigned to a tier by an indirect method
   (retention-time clustering, stoichiometric peak-count matching), any
   rate constant that claims to distinguish *within* an unconfirmed tier is
   making a regiochemistry claim the data cannot support — the assignment
   of "which peak is isomer A vs. isomer B" is not trustworthy from kinetics
   fitting alone, only the tier (degree of substitution) is.
2. **Data-points-to-parameters ratio.** Count usable data points against
   the number of free parameters the network structure implies. A full
   site-resolved lattice with every isomer pathway explicit can easily reach
   parameter counts that are a large fraction of the usable data points —
   a red flag independent of whether the fit converges cleanly. A pooled
   representation that lumps same-tier isomers into one balance cuts the
   edge count (and so the parameter count) roughly in proportion, giving a
   much healthier ratio for the same data.
3. **Whether lumping is even valid, physically.** Pooling by molarity is
   only correct if every species being lumped together shares a molecular
   weight (or whatever normalizing property the balance uses) — check this
   explicitly, don't assume it. When it holds (e.g. every isomer at a given
   degree-of-substitution tier differs from the parent only by the same
   fixed mass addition, regardless of which site reacted), pooled
   concentration = sum of tier-member concentrations, and this is not an
   approximation, it's exact.
4. **AIC/BIC comparison, and read the boundary-pinning check.** Fit both a
   simpler baseline (e.g. integer/first-order kinetics) and the richer
   free-parameter version, and compare with an information criterion, not
   raw SSE (which always improves with more parameters). A decisive ΔAIC in
   favor of the richer model is necessary but not sufficient — also check
   whether any freed parameter landed pinned against its box bound. A
   parameter sitting exactly at its fence means the *direction* of the
   effect is supported by the data but the specific fitted value is not —
   report the qualitative finding, don't quote the boundary value as if it
   were a converged point estimate. Widening the bound and refitting (warm-
   started from the pinned result, not from scratch) tells you whether it
   keeps sliding (data wants an even more extreme value — reconsider whether
   the rate-law functional form itself is right, e.g. switch from a
   power-law order to a saturation/Michaelis-Menten form) or settles
   (the pin was a bound problem, not a form problem).

**The pooling convention** (concrete pattern that generalizes): define each
lumped species as the sum of every confirmed-or-inferred peak that shares
both (a) the same degree-of-substitution tier and (b) the same normalizing
property (MW), keep any individually-confirmed species un-pooled and
tracked on its own, and build the reaction network's edges at the pooled
level (pooled species -> pooled species), not by collapsing a fully
resolved network after the fact. **Pooling must happen before the fit, and
the same pooling must be applied again, explicitly, before any downstream
model consumes the fit** — never let a downstream reactor/scale-up model
implicitly assume it has finer species resolution than the upstream fit
actually resolved. If the downstream model's natural state vector is finer
than the fit (e.g. it wants individual isomers), it has to be pooled down to
match at the point of loading, not silently run at mismatched resolution.

**Illustrative case (generalize the lesson, not the numbers):** in this
project, an 8-peak site-resolved kinetic network was collapsed to 5 pooled
species — two peaks individually confirmed against a reference standard,
the rest pooled by degree-of-substitution tier — cutting the reaction
network from a full isomer-resolved lattice (a 12-edge, ~15-parameter
network) to a 5-edge, 8-parameter one, which fit substantially better
per-parameter against a data set of a few hundred usable points. The general
lesson: prefer the coarsest species resolution that (a) respects what's
actually identity-confirmed and (b) leaves a healthy data-to-parameter
ratio, and treat any finer resolution as a hypothesis to justify explicitly,
not a default to assume because the chemistry could in principle support it.

## 3. Fit versioning: when is a kinetic fit "current enough" to build on

**Before reusing an existing fit artifact, check three things:** has the
underlying experimental data changed (new batches, corrected data
processing), has the model structure changed (different species resolution,
different rate-law functional form, different reaction network), and does
the artifact's own metadata (R², AIC/BIC, boundary-pinning flags) still read
as a defensible fit rather than a placeholder or a superseded intermediate
step. Any one of these being true means refit; none being true means the
existing artifact is still current.

**What actually triggers a fit becoming stale, concretely:** a change in
species resolution (the pooled 5-species representation is a structurally
different fit target than an older full-regiochemistry representation — its
parameters are not interchangeable, one cannot substitute for the other even
if both nominally describe "the same reaction"); a widened parameter bound
after a boundary-pinning finding (the refined fit is not optional polish, it
supersedes the pinned one because the pinned one's headline parameter value
was flagged as untrustworthy); and a change in what's held fixed vs. free
(freeing a previously-fixed reaction order is a different optimization
problem, not a continuation of the old one, even though it's warm-started
from it).

**What does NOT automatically trigger a refit:** reusing a fit's parameters
in a different downstream model (e.g. a different compartment count, a
different reactor scale) is fine as long as the species resolution and rate
law are unchanged — that's a legitimate reuse, not staleness, provided you
still separately validate that the reused parameters make physical sense in
the new context (see §4 on carrying transport parameters across scale, which
generalizes the same idea to physical rather than kinetic parameters).

**Naming and lineage discipline:** date-stamp fit artifact filenames (so
`fit_results_YYYY-MM-DD` variants coexist rather than silently overwriting
each other across iterations), and never let a filename alone be the only
record of what changed between versions — save a self-describing bundle
(the raw optimizer result, fit-quality metrics like R²/RMSE, and enough
metadata to know which data and which model structure produced it) rather
than a bare parameter vector. When a refinement supersedes an earlier
result, say so explicitly in whatever loads it next (e.g. "refined,
supersedes the `_results.pkl` from the same date") rather than leaving two
similarly-named files with no stated precedence — a downstream notebook or
model that picks the wrong one silently is a much worse failure mode than
one that fails to find a file at all.

**Illustrative case:** this project moved from a 9-species full-
regiochemistry `.dill` fit to a 5-species DS-Pool `.pkl` fit specifically
because of the species-resolution argument in §2 — the older fit claimed
resolution the analytical method didn't support. Within the DS-Pool line
itself, a first fit (`n=m=1` baseline, then free `n,m`) found the freed
reaction order pinned at its lower bound; a second, explicitly-named
"refined" pass widened that bound and warm-started from the pinned result,
producing the fit now treated as current. Both transitions are documented
in the notebooks' own markdown, not just inferred from filenames — that's
the standard to match: the *reason* a fit was superseded should be
recoverable from the artifact trail, not just the fact that it was.

**When reuse across scale/condition is valid vs. when a new fit is
mandatory:** a rate-constant kinetic fit (chemistry-intrinsic) generally
transfers across reactor scale and geometry unchanged, since the reaction
itself doesn't know what vessel it's in — but check whether the process
conditions it was fit under (temperature, pH, solvent composition) are
still representative; empirical rate constants carry no explicit dependence
on a variable unless that variable was varied in the fitting data, so
operating outside the fitted range on such a variable is an extrapolation
to flag, not silently trust. A *transport* parameter (a flow number, a mass-
transfer coefficient) is different: it is generally NOT scale-invariant and
should be treated as vessel- and condition-specific — reuse across scale
only via an explicit, checked scaling law (see §4), never by direct
carry-over, and reuse across RPM/agitation condition within the *same*
vessel only after confirming the flow regime (turbulent plateau vs. not,
vortexing vs. not) is comparable at both conditions.

## 4. Scale-up / scale-down

**Dimensionless groups to check at the target scale, and what each one
governs:**
- **Impeller Reynolds number** (`Re_i = ρND²/μ`) — governs whether the flow
  number itself can be treated as a geometry-only constant (turbulent
  plateau, conventionally `Re_i >= ~10,000`) or is still Reynolds-dependent.
  Check this first at *both* scales — a scale-down to a much smaller vessel
  at the same RPM can fall out of the turbulent regime even if the parent
  vessel was safely in it.
- **Froude number** (`Fr_i = N²D/g`) — governs vortex onset in an unbaffled
  vessel; relevant only for unbaffled geometries, but check it explicitly
  rather than assuming baffling status carries over between scales.
- **Power per unit volume** (`P/V`, `P = N_p·ρ·N³·D⁵`) — a common
  scale-down invariant when the process is shear- or mixing-energy-
  sensitive.
- **Impeller tip speed** (`π·N·D`) — a common invariant when the process is
  sensitive to peak local shear at the blade rather than bulk energy input.
- **Mixing time** (`θ_95`, time to reach ~95% homogeneity — for a linear
  compartment system this is `ln(20)/λ_slow`, the slowest transport-matrix
  eigenvalue) — the group that actually matters when comparing to the
  reaction's own timescale (§2 of `parameter-estimation.md` covers deriving
  this from fitted compartment flow rates).

**The scale-down decision ladder — constant N, constant P/V, constant tip
speed are mutually exclusive; you cannot hold more than one fixed while
changing vessel size, so pick based on what the process is actually
sensitive to, not by default:**
- **Constant N (impeller speed)** preserves relative mixing pattern/timing
  best but lets P/V and tip speed both drift, usually upward as vessel size
  shrinks (smaller D needs more N to preserve tip speed, but constant-N
  scale-down doesn't even try) — reasonable only when the process is known
  to be insensitive to shear and to specific power input, or as a first
  baseline candidate to compare against the other two.
- **Constant P/V** preserves specific power input (good proxy for bulk
  mixing/mass-transfer intensity per unit volume) — the standard choice for
  processes where bulk mixing energy density is believed to be the
  controlling variable (e.g. mass-transfer-limited operations), at the cost
  of tip speed and mixing time generally not matching the parent scale.
- **Constant tip speed** preserves peak shear at the impeller blade — the
  standard choice for shear-sensitive processes (e.g. protein/cell
  integrity concerns), at the cost of P/V typically dropping sharply at
  smaller scale (power input falls off faster than volume as you shrink
  under this criterion).
- Report Re, Fr, tip speed, and mixing time together for *all three*
  candidate criteria side by side before picking one — the right choice
  depends on which downstream number moves into a concerning regime (e.g.
  Re dropping out of turbulent, Fr crossing into vortexing) under each
  criterion, which is exactly the kind of thing that's easy to miss if you
  only compute the one number the criterion is nominally about.

**Geometry changes that need special treatment, not a naive linear
scale-down:** a torispherical (dished) bottom head needs its own exact
volume formula (function of crown radius, knuckle radius, and head depth,
solved per standard — ASME/F&D vs. DIN/Klopper use different knuckle
ratios) rather than being approximated as a flat-bottom cylinder, especially
at small scale where the head is a large fraction of total volume. A
different impeller-to-tank diameter ratio (`D_impeller/D_tank`) between
scales can produce genuine *structural* changes, not just parameter
changes: if the impeller's nominal swept zone, computed the same way at the
smaller scale, extends past the vessel's physical bottom, that's a real
finding that a compartment believed to exist at the parent scale (e.g. a
quiescent dead zone below the impeller) may not exist at all at the new
scale — the compartment count itself needs to be re-derived at the new
scale's geometry, not assumed to carry over (see §2 of
`parameter-estimation.md` on the residence-time/data-feature-count decision
ladder — apply it fresh at each scale, not once). A different baffle
configuration between scales invalidates any Froude-number-derived vortex
ceiling computed at one scale being reused at the other — baffled and
unbaffled vessels have qualitatively different `Nq(Fr)` behavior, not just
different constants.

**Illustrative case:** a 5 L reference vessel's validated 3-compartment
transport model, scaled down to 250 mL with *exact* measured geometry and a
torispherical head, was found to structurally degenerate to 2 compartments
— the impeller's swept-zone boundary, computed by the same rule as at the
parent scale, fell below the physical vessel base at the smaller scale's
impeller-to-tank ratio. The general lesson: don't port a compartment count
across scales as a fixed fact about "how many zones this kind of vessel
has" — the number is a joint function of geometry ratios and impeller
position that must be recomputed, and can legitimately come out different,
at each new scale.

## 5. Batch validation: forward prediction against real data

**Workflow:** freeze every fitted/validated parameter from the template
model (kinetics, transport coefficients, numerical safeguards) — change
nothing about them — and run the forward model against a batch's actual,
fully-specified process conditions (its own charge, feed schedule, RPM,
timing) that were not used to fit anything. This is a "zero optimization"
pass: no peak detection, no sweep, no root-finding for a target endpoint.
The point is to see what the already-fitted model predicts when handed a
real, independent set of conditions, not to tune anything to match the
outcome.

**What counts as a passing validation vs. a red flag:**
- **Passing:** predicted trajectory shape and endpoint are consistent with
  the batch's observed outcome within run-to-run scatter, *and* mass balance
  (or whatever hard conserved-quantity check applies) holds to numerical
  precision, *and* any operating condition that falls outside the range the
  template was validated at (e.g. a much higher RPM than the validated SOP)
  is explicitly flagged and checked against its own relevant physical limit
  (e.g. re-running the vortex/turbulence checks at the new condition), not
  silently assumed to still be safe.
- **Red flag:** a prediction that requires silently reusing a numerical
  safeguard outside the context it was designed for. E.g. a terminal
  "reaction exhausted" stopping event built for a template that always ends
  the simulation once a reagent is depleted is the wrong tool for a batch
  protocol that has a fixed-duration hold *after* depletion — reusing it
  unmodified would silently truncate real hold time from the simulation.
  The fix is to adapt the safeguard's *mechanism* (e.g. gate the reaction
  rate terms to zero at the same physical threshold, but keep integrating
  in real time) while keeping its physical criterion unchanged — and to
  explicitly verify what happens if you *don't* adapt it (in this project,
  leaving the near-zero-order kinetics unguarded past depletion was
  confirmed to drift a reactant to non-physical negative values over a
  realistic hold duration) so the adaptation is justified by evidence, not
  just plausible-sounding caution.

**Handling run-to-run variability:** distinguish deviations that fall within
plausible experimental scatter (different batches, same nominal recipe) from
deviations that reveal a missing mechanism. A useful discipline: when a run
uses conditions genuinely outside the validated envelope (different RPM,
different impeller height, different stock concentrations), predict and
report the *consequences* of that specific difference explicitly (e.g. "at
this much higher circulation flow, the feed duration is now short enough
relative to circulation that the reagent is predicted to be numerically
exhausted long before the next scheduled addition") rather than treating an
unexpected-looking trajectory as automatically suspect — a large deviation
that's directly explainable by a stated, checkable difference in conditions
is not the same failure mode as an unexplained one.

**A forward prediction that "looks good" on yield/purity can still be
suspect** if it gets there via the wrong trajectory or for the wrong reason
— check the mechanism behind an endpoint number, not just the number.
Concretely: if a downstream reagent addition (e.g. a quench step) is
predicted to have negligible effect on the endpoint, check *why* — "the
dose was too small to matter" and "the reagent it was supposed to react with
was already fully consumed before it arrived" are different findings with
different implications, and only inspecting the trajectory (not just the
endpoint table) distinguishes them. Similarly, a mass-balance check passing
confirms the model is self-consistent, not that its trajectory is right —
use it as a necessary, not sufficient, validation criterion.

**Illustrative case:** a template model validated at one RPM/impeller
configuration was run forward, unmodified, against two real batches
specified at a much higher RPM and different impeller height. The forward
run correctly flagged that the new RPM exceeded the vortex-safety ceiling
the template had been deliberately capped below, and separately revealed
(by inspecting the full trajectory, not just the endpoint) that the faster
circulation and much shorter feed duration at this condition caused the key
reagent to be predicted numerically exhausted tens of minutes before a
scheduled downstream addition — making that addition's dose and timing
irrelevant to the predicted outcome regardless of how it was specified. That
finding came from tracing the mechanism, not from noticing that an endpoint
number looked odd.

## 6. Portability notes

**Changes between projects:** the specific reaction/species network and
which peaks or analytes are identity-confirmed vs. inferred; reactor
geometry, scale, and impeller configuration; the analytical method used to
track species over time (HPLC peaks here, but equally could be spectroscopic
signals, mass-spec transitions, or any other quantitative per-species
signal) and what it does and doesn't resolve; the specific dimensionless-
group targets (turbulent Reynolds threshold, vortex Froude threshold) which
are geometry- and baffle-configuration-dependent constants, not universal
numbers.

**Stays the same across projects:** the species-resolution decision
criteria (confirmed-vs-inferred identity, data-points-to-parameters ratio,
AIC/BIC plus boundary-pinning as the test for "does freeing this parameter
group earn its keep," pooling as an exact operation under a stated validity
condition rather than an approximation); the fit-versioning discipline
(dated, self-describing artifact bundles with explicit supersession
statements, treating a structural model change as fit-invalidating even
when a filename doesn't change); the scale-up dimensionless-group hierarchy
and the mutually-exclusive nature of constant-N / constant-P/V / constant-
tip-speed scale-down criteria, chosen by what the process is sensitive to
rather than by default; re-deriving compartment count fresh at each new
scale's actual geometry rather than porting a count; and the batch-
validation workflow structure (freeze all parameters, run forward against
independently-specified real conditions, flag out-of-envelope conditions
explicitly, trace mechanism behind an endpoint rather than trusting the
endpoint alone).
