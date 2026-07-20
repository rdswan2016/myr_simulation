# Parameter estimation for mechanistic bioprocess models

## 1. When to use this skill

Use this when you need to back out unknown physical parameters — rate
constants, reaction orders, mass-transfer/circulation flow rates, mixing
coefficients — from experimental time-series data by fitting a mechanistic
(ODE/PDE) forward model to it. Concretely: fitting a kinetic rate law to a
depletion/formation curve, fitting an inter-compartment flow rate to a
tracer or titration curve, fitting a mass-transfer coefficient to a probe
signal (pH, conductivity, absorbance, concentration). It does not cover
purely statistical/black-box regression — the point of this workflow is that
the model structure is physically motivated and the fitted parameters are
meant to be reused in a *different*, downstream mechanistic simulation.

## 2. Compartment number determination for STR problems

This is the part Claude must drive, not the user. Most users asking for a
stirred-tank-reactor (STR) model do not know a priori whether 1, 2, or more
compartments are appropriate — that is an engineering judgment call to
propose and argue for explicitly in the notebook/report, not a question to
punt back to the user.

### Engineering reasoning that feeds the decision

- **Reynolds number regime** (`Re_i = N·D²/ν`, impeller Reynolds number).
  Below the fully turbulent plateau (rule of thumb `Re_i >= ~10,000`, per
  standard transport-phenomena correlations, e.g. Welty et al.), the flow
  number `Nq` is not a constant — it depends on `Re_i` itself, and any
  `Q = Nq·N·D³` pumping-rate estimate built on a fixed `Nq` is invalid until
  you've checked you're actually in the plateau. Check this first, before
  anything else — it gates whether "flow rate scales linearly with RPM" is
  even a valid assumption in the rest of the analysis.
- **Baffling / Froude number** (`Fr_i = N²D/g`). An unbaffled vessel develops
  a free-surface vortex once `Fr_i` exceeds an empirical ceiling
  (`Fr_crit ≈ 0.15` is a reasonable literature default for unbaffled tanks).
  Past that ceiling, the fluid co-rotates with the impeller, effective
  `Nq` *collapses* (can drop an order of magnitude, not gently roll off), and
  `Q = Nq·N·D³` overestimates true pumping. A single angled probe or
  thermowell does not count as a baffle. If the vessel is unbaffled, check
  `Fr_i` across the RPM range of interest before trusting any Nq-derived flow
  number as an RPM-independent impeller property — it usually isn't, and the
  fitted value should be reported as apparent/effective at that specific RPM.
- **Mixing/residence timescale vs. process timescale.** This is the single
  most decisive check for compartment count. Compute each candidate
  compartment's residence time (`V_zone / Q_exchange`) and compare it to the
  characteristic timescale of the process you actually care about (reaction
  half-life, feed duration, dosing time). If residence time is far shorter
  (order-of-magnitude rule: ~50-200x) than the process timescale, the vessel
  is numerically indistinguishable from perfectly mixed *for that process* —
  regardless of how hydrodynamically complex the real vessel is. This can
  make a compartment model that is genuinely necessary for fitting the mixing
  data itself (e.g. a titration curve, which unfolds on the same
  seconds-to-minutes timescale the compartments resolve) collapse back to
  well-mixed once you plug its fitted flow rates into a *different*, slower
  downstream process (e.g. a reaction that takes tens of minutes). Compartment
  count is a property of the (model, process) pair, not of the vessel alone —
  don't reuse a compartment count across processes without redoing this
  comparison.
- **Data feature count (the direct experimental signature).** A chain of `N`
  well-mixed compartments has exactly `N-1` non-zero relaxation
  eigenvalues/timescales. If your experimental trace shows one clean
  exponential relaxation, one exchange coefficient (2 compartments) is the
  ceiling of what's identifiable from that data — don't fit more. If it shows
  two visually and temporally separated relaxation features (e.g. a fast
  rebound within seconds followed by a much slower drift over many minutes),
  a single lumped compartment *cannot* reproduce two well-separated time
  constants from one exchange coefficient — that is direct evidence you need
  at least 3 compartments. Let the data's own feature count set the floor on
  compartment number; don't add compartments beyond what distinguishable
  features justify.
- **Geometry / mechanism heuristics for where boundaries go.** Anchor
  compartment boundaries to physical features when you have them: impeller
  height as the natural zone boundary generator, impeller diameter as an
  order-of-magnitude proxy for the thickness of its actively-swept discharge
  zone. When there's no such anchor available (e.g. liquid fill level doesn't
  even reach the impeller-based boundary), fall back to the simplest
  defensible geometric split — equal-height bands — and say so explicitly
  rather than force-fitting an anchor that doesn't apply at that fill level.
  Specific impeller mechanism matters for *where* to suspect an extra
  compartment is needed: an up-pumping axial impeller (discharge jet toward
  the surface, diffuse return flow down the walls) is prone to leaving a
  poorly-swept dead zone directly under/around the impeller near the tank
  floor, especially when off-bottom clearance-to-diameter (`C/D`) is on the
  high side (order 1 or more). A down-pumping impeller has the opposite
  bias. Use this to generate a specific, testable hypothesis about where a
  4th (or Nth) compartment belongs — not to justify adding compartments
  generically.

### The decision ladder

1. **1 compartment (perfectly mixed)** — default when `Re_i` is turbulent and
   the mixing/residence timescale is much shorter than the process timescale
   you're actually simulating. This is a defensible default even for a
   hydrodynamically complex real vessel, *for that specific process* — state
   the residence-time-vs-process-time comparison that justifies it rather
   than asserting it.
2. **2 compartments** — the minimum whenever there's a physically distinct
   feed point and/or measurement point separate from the bulk (e.g. a feed/
   probe zone near the surface vs. everything below the impeller), and the
   experimental data shows a single relaxation timescale. This is also the
   right default when you don't yet have mixing data to justify anything
   finer — a 2-zone feed/bulk split is the standard first mechanistic model
   for an STR before you have a titration or tracer curve to fit.
3. **N > 2 compartments (chain topology)** — warranted when the data shows
   multiple relaxation timescales that one exchange coefficient cannot
   reproduce, or when a specific physical mechanism (axial gradient, feed-
   point asymmetry, a suspected dead zone from impeller pumping direction)
   generates a testable hypothesis for an additional zone. Use an
   adjacent-zone-only (chain) topology, not all-pairs exchange, when mass
   physically must pass through an intermediate zone to reach a farther one
   (e.g. a vertical axial stack under one impeller) — this is a topology
   decision to make from the physical layout, not something to fit.

### Validating the chosen number

Don't stop at "the data has two features, so 3 compartments." Validate:

- **Sensitivity test.** Refit with one more compartment than you settled on
  and check whether previously-fitted parameters move by more than ~5%, or
  whether the added compartment's own flow rate is even identifiable
  (profile-likelihood, see §4) rather than pinned to the edge of its bounds
  from every starting guess. A parameter that lands at essentially the same
  value regardless of how many downstream compartments exist, or where the
  optimizer starts, is good evidence it reflects a real physical exchange
  rate rather than an optimizer artifact — that stability *is* the
  validation, actively look for it.
- **Test the hypothesis, don't just assert it.** If you added a compartment
  to represent a specific physical mechanism (e.g. a dead zone), don't stop
  once the richer model exists — check whether forcing that compartment's
  flow rate to the value the mechanism predicts helps or hurts the fit
  (profile-likelihood comparison of cost). A mechanism that seems physically
  plausible can still be wrong-signed: e.g. a genuinely stagnant zone
  *withholds* buffering/reactive capacity from the fast-responding pool,
  which is the wrong direction if the data feature you're trying to explain
  needs *more* capacity reaching the fast zone over time, not less. Report
  the result plainly whichever way it comes out — a well-motivated
  compartment hypothesis that the data rejects is still a useful, reportable
  finding, not a failure to hide.
- **Structural checks independent of the fit.** Before trusting an overlay
  plot, ask what the best possible fit could *ever* look like independent of
  the specific fitted parameters — e.g. a closed system's fully-mixed
  equilibrium value (set purely by conserved total amounts and final volume)
  is a hard ceiling no choice of exchange rates in a transport-only model can
  exceed. If the data exceeds that ceiling, no amount of adding compartments
  within that model class will close the gap — that's evidence the model
  class itself is missing a mechanism, not that you haven't found the right
  compartment count yet.

### Default recommendation when geometry/mixing data is unavailable

Start at 2 compartments (feed/probe zone vs. bulk) for any STR problem where
you don't yet have experimental mixing data to justify more, and 1
(well-mixed) if you've separately confirmed the process of interest is much
slower than any plausible mixing timescale for that vessel scale. Only go to
3+ once you have a tracer/titration/probe curve in hand that actually shows
multiple relaxation features — don't pre-emptively build a finer compartment
model than the data can identify.

## 3. Parameter estimation workflow

1. **Set up the forward model** as an ODE (or PDE, method-of-lines) with the
   unknown parameters as explicit function arguments, not hardcoded — you'll
   call it repeatedly from the optimizer's residual function.
2. **Handle discontinuous forcing explicitly.** If the process has a feed or
   dosing schedule that turns on/off (a discontinuous source term), integrate
   in separate segments split exactly at the discontinuity rather than
   letting the adaptive solver step across it — an adaptive-step solver can
   silently smear over or miss a sharp transition otherwise.
3. **Build the residual function**: `predicted(params, t) - measured(t)`, fed
   unweighted (or explicitly weighted, if you have a reason) into
   `scipy.optimize.least_squares`. Prefer `least_squares` with `method='trf'`
   over unconstrained `curve_fit` whenever you have physical bounds — a
   bounded fit prevents the optimizer from wandering into unphysical
   parameter territory mid-search, which matters for stiff/singular
   forward-model evaluations too.
4. **Derive bounds from physics, not arbitrary numbers.** An arbitrarily
   large box bound (e.g. `1e4`) is a trap: if the optimizer lands *at* it,
   that's not a converged answer, it's a sign the parameter is poorly
   identified and the optimizer walked to the edge of the search space. Treat
   an optimizer landing at an arbitrary bound as a red flag requiring
   investigation, not a result to report. Replace arbitrary bounds with a
   physically-derived ceiling wherever one exists (e.g. an impeller's total
   plausible pumping capacity from a generous literature `Nq_max`, a
   stoichiometric ceiling from total charged mass) — then check whether the
   refit cost is statistically indistinguishable from the unconstrained fit.
   If so, the physically-capped value is strictly more interpretable at no
   cost in fit quality, and should become the number you report.
5. **Choose initial guesses deliberately.** When extending a validated
   simpler model (fewer compartments, fewer parameters) to a richer one,
   seed the new fit from the simpler model's already-fitted values for shared
   parameters, and use a physically-informed guess (e.g. from an independent
   analytical cross-check, see §4) for genuinely new parameters — not an
   arbitrary default.
6. **Run a multi-start robustness check.** Refit from several different
   initial guesses spanning orders of magnitude and confirm convergence to
   the same point (or, if not, that it's converging to the same bound every
   time — itself informative about identifiability). Report this, don't just
   run the fit once and trust it.

## 4. Fit quality evaluation

- **R²/RMSE tell you goodness-of-fit, not identifiability — check both.** A
  good R² can coexist with parameters that are not individually meaningful.
- **Check the Jacobian condition number at the solution** (SVD of `fit_result.jac`,
  or the ratio of largest to smallest singular value). A large condition
  number (order `1e6+`) signals near-singular parameter correlation — two
  parameters are trading off against each other along the data's information
  content, not independently determined. When this happens, the naive joint
  covariance (`inv(J^T J)`, what `curve_fit`'s default stderr computes) is
  numerically invalid and can produce a reported uncertainty *larger than the
  fitted value itself* — a strong tell that something's wrong with the error
  bars, not the fit.
- **Use profile-likelihood instead of naive covariance when the Jacobian is
  ill-conditioned.** Fix the poorly-identified parameter at a grid of values,
  refit everything else at each point, and compare cost against a chi-square
  threshold (`cost_best + 0.5 * chi2.ppf(0.95, df=1) * sigma²`, with
  `sigma² = 2·cost_best/(n_points - n_params)`) to get a defensible
  confidence interval — which may be one-sided (only a lower or upper bound
  identifiable, not both). Report a bound as a bound, not as a fabricated
  symmetric `value ± stderr` that implies more precision than the data
  supports.
- **Use structural/analytical cross-checks independent of the nonlinear
  optimizer wherever the model structure allows one.** Two useful patterns
  observed in practice: (a) a conserved-quantity ceiling — a value the system
  must approach regardless of the fitted parameters, computable directly from
  input amounts and final volume/state, giving a hard bound on achievable fit
  quality that isn't an optimization artifact; (b) an eigenvalue/timescale
  cross-check — for a linear (or locally-linearized, e.g. post-feed) system,
  the achievable relaxation timescales are fixed functions of the fitted
  transport parameters alone, computable with no fitting involved, and can be
  directly compared against what a data feature's own empirically-fit
  timescale requires. When the cross-check's required parameter value falls
  outside the profile-likelihood-supported range, that's independent,
  fit-optimizer-free confirmation (or refutation) of a hypothesis.
- **When is a fit good enough?** When (a) cost is statistically indistinguishable
  between a looser and a physically-constrained bound (prefer the constrained
  one), and (b) remaining residual structure (e.g. a gap between data and a
  hard analytical ceiling) is attributable to a stated, out-of-scope
  mechanism rather than unexplained. When residual structure instead tracks
  a data feature the current model class cannot produce even in principle
  (see the conserved-quantity-ceiling check), that's the signal to
  reparametrize — but reparametrizing means testing a specific new mechanism
  as a hypothesis (see §2's validation subsection), not blindly adding
  degrees of freedom until the residual shrinks.

## 5. Connecting fitted parameters to downstream models

- **Export a self-describing bundle, not bare numbers.** Pickle/dill a dict
  containing the raw optimizer result (or its parameter vector), plus fit
  quality metrics (R², RMSE) and enough metadata to know which run produced
  it. Use a dated, descriptive filename so successive fit iterations don't
  silently overwrite each other or get confused downstream. Prefer `dill`
  over stdlib `pickle` when the bundle holds a non-trivially-picklable object
  like a `scipy.optimize.OptimizeResult`.
- **Parameterize for positivity during fitting, un-transform at load time.**
  If rate constants must be positive, fit in log-space
  (`x = log(k)`, optimizer sees unconstrained/loosely-bounded `x`) and
  `exp()` only at the point the downstream model loads the bundle — document
  this convention next to both the export and the load code, since it's a
  silent source of a 100x-off bug if one side forgets the transform.
- **Write the downstream loader defensively.** Wrap the pickle/dill load in
  try/except with a hardcoded fallback of the last-known-good fitted values
  baked directly into the loader function. This keeps the downstream
  simulation runnable even if the artifact file moves, is regenerated with a
  different name, or the fitting notebook hasn't been rerun yet — and print
  explicitly which path was used (loaded-from-file vs. fallback) so it's
  never ambiguous which parameters actually produced a given simulation
  output.
- **Reconcile state resolution before consuming a fit.** If the downstream
  model tracks a different (usually coarser) set of species/states than what
  the upstream fit resolved, pool the fit's states to match *before* feeding
  it in — don't let the downstream model implicitly pretend it has finer
  resolution than the fit actually supports.
- **Carry transport parameters into a different downstream model only after
  checking they still matter there.** A flow rate fitted from a fast mixing
  experiment can be reused directly in a slower downstream reaction model's
  transport terms — but explicitly redo the residence-time-vs-process-time
  comparison from §2 for the new (model, process) pair. It's entirely
  possible, and worth showing rather than hiding, that a compartment
  structure real and necessary for the upstream fit turns out to be
  numerically indistinguishable from well-mixed once carried into the
  downstream process — that's a legitimate finding, not a wasted step.
- **When a carried-over parameter is only a bound, say so downstream too.**
  If the upstream fit could only establish a profile-likelihood bound (not a
  point value) for a transport parameter, state explicitly which bound you
  chose to use in the downstream deterministic simulation and why (e.g. "used
  the 95% lower bound as the most conservative/most-differentiated scenario
  the data doesn't rule out") — the caveat needs to propagate with the
  number, not get silently dropped at the handoff.

## 6. Portability notes

**Changes between projects:** reactor geometry and scale; impeller type and
pumping direction (up- vs down-pumping changes where a dead zone is
suspected); the specific reaction/species network; which physical quantity
is actually measured (pH via Henderson-Hasselbalch here, but equally could be
absorbance, conductivity, or direct concentration elsewhere — whatever the
probe reports, and the corresponding forward-model observable used in the
residual); literature `Nq`/`Np` ranges, which are impeller-type-specific;
feed/dosing schedule shape.

**Stays the same across projects:** the Re/Fr regime check as a gate before
trusting any flow-number-based estimate; using residence-time-vs-process-
time as the decisive compartment-count heuristic, redone per (model,
process) pair; using the experimental data's own relaxation-feature count as
a hard floor on identifiable compartment number; preferring physically-
derived optimizer bounds over arbitrary numeric ones, and treating an
optimizer parked at an arbitrary bound as an identifiability red flag;
checking Jacobian condition number before trusting any reported parameter
uncertainty, and falling back to profile-likelihood when it's ill-
conditioned; using structural cross-checks (conserved-quantity ceilings,
linear-system eigenvalue timescales) that don't depend on the nonlinear
optimizer having converged correctly; testing an added-compartment hypothesis
by seeing whether the data actually rewards it, not just building the richer
model and stopping; and the log-parameterize-plus-dilled-bundle-with-fallback
pattern for handing fitted parameters to a downstream mechanistic simulation.
