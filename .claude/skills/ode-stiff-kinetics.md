# Debugging stiff / ill-conditioned ODE kinetics

Use this when a mechanistic ODE or PDE simulation (reaction kinetics, transport,
population dynamics, anything with `solve_ivp`/`odeint`/similar) produces
implausible output: mass-balance violations, states crossing a physical bound
(negative concentrations, >100% yield, etc.), solver stalls, or blow-up. This is
especially likely when the rate law has a fractional or near-singular exponent
(e.g. rate ∝ c^n with n < 1 near c = 0), when rate constants span many orders of
magnitude (classic stiffness), or near a depletion/exhaustion boundary.

## Diagnostic-first, not patch-first

Do not apply a fix on the first implausible result. Root-cause it with concrete
printed evidence before changing anything:

1. **Look at raw solver nodes, not the dense/interpolated output.** Dense output
   (`sol.sol(t)` / interpolation) can mask or fabricate excursions past a
   physical bound that never occurred at the actual integration steps, or vice
   versa — hide a real one. Print/plot `sol.t`, `sol.y` directly.
2. **Check solver step sizes.** A step size collapsing toward zero is the
   signature of the solver fighting stiffness or a near-singular Jacobian, not
   a model bug — it tells you *why* a naive fix (e.g. a hard clamp) will trap
   the solver rather than resolve it.
3. **Check exact conserved quantities** (total mass, total moles, charge
   balance) at every step, not just at the end. A conserved quantity that holds
   globally but is wrong at the species level (e.g. one species pinned
   negative while the total still balances) tells you the bug is in how an
   individual state is regularized, not in the overall mass accounting.
4. **Locate the mechanism**, in order of likelihood: (a) a rate law with a
   fractional/negative exponent going singular as a reactant depletes, (b) rate
   constants spanning >3-4 orders of magnitude, (c) a depletion asymptote that
   is analytically unbounded in integration time (approaches zero but never
   exactly reaches it), (d) a stiff solver (BDF/Radau/LSODA) locked onto a
   spurious flat branch after a clamp discontinuity.

## Fix patterns, and why the naive ones fail

Tried in roughly this order of naivety; each failure mode below was observed
directly, not hypothesized:

- **Hard clamp** (`max(0, c)` / `np.maximum`): creates a discontinuity in the
  derivative at the bound. A stiff implicit solver can lock onto this as a
  spurious equilibrium — the state gets absorbed at exactly the clamp value
  and stops evolving, even though the "true" trajectory should continue past
  it. Aggregate/conserved quantities can look fine while the trapped species is
  silently wrong.
- **Sign-preserving continuation** (let the state go negative and continue
  evaluating the rate law with `sign(c)*|c|^n` or similar): avoids the
  discontinuity but does nothing about the underlying singular curvature. With
  large rate constants this commonly causes runaway blow-up (states reaching
  many orders of magnitude past the physical scale of the system) instead of
  a stall.
- **Smooth non-negative floor** (what actually worked): replace `c` with a
  smooth approximation to `max(c, 0)` everywhere it feeds a fractional-power or
  otherwise singular term, e.g.
  `c_floored = 0.5 * (c + sqrt(c^2 + eps^2))`
  with `eps` small relative to the system's characteristic concentration scale
  (pick it near or below solver tolerance, e.g. `eps ~ 1e-9` in the same units
  as `c`). This keeps the derivative continuous (no clamp discontinuity for
  the solver to lock onto) while keeping the rate law finite near zero (no
  blow-up). Apply it at the point of evaluation inside the rate law, not to
  the state vector itself.
- **Terminal event for the true asymptote**: if the underlying kinetics have a
  reactant that only approaches zero asymptotically (never exactly reaching
  it), the smooth floor alone won't stop the solver from grinding through an
  unbounded tail. Add a terminal `solve_ivp` event that stops integration once
  the relevant quantity drops below a small fraction (e.g. 0.01%) of its
  initial value. Without this, integration time/step count becomes effectively
  unbounded for no gain in accuracy.
- **A smooth floor on the rate law is not the same as a floor on the state.**
  If every reaction consuming a species has a rate ≥ 0, the ODE for that
  species is monotonically non-increasing everywhere — there is no restoring
  force once a stiff step overshoots past zero. The floored rate law keeps
  things numerically well-behaved but the *state itself* can still drift
  slightly negative (a few numerical-tolerance-scale units) after the
  species is nominally exhausted, particularly in a later integration stage
  that lacks its own terminal exhaustion event. Confirm this is cosmetic (check
  it doesn't move mass balance/yield/purity outputs, which should be computed
  from the raw unclamped state) before deciding whether a display-only clip is
  acceptable or whether the dynamical fix (extend the terminal event, or add
  an explicit restoring term) is warranted.

## When to stop patching and ask

If a second numerically-motivated fix doesn't obviously work on inspection, stop
and surface the tradeoff to the user before trying a third. The two concrete
options are usually:

- Keep debugging the numerically stiff/singular formulation (accurate to the
  fitted kinetics, more fragile integration), or
- Fall back to a simpler, well-behaved kinetic model (e.g. integer reaction
  orders instead of a fitted fractional exponent) that trades some fit
  accuracy for robustness.

Don't unilaterally pick one — this is a modeling judgment call, not a purely
numerical one. Always validate whichever fix is chosen against a concrete,
checkable criterion (exact mass balance to N decimals, no blow-up over the
full integration horizon, conserved quantity held at every raw solver step)
before treating it as done.

## Portability notes

The specific pattern (smooth floor + terminal depletion event) generalizes to
any ODE/PDE system with a rate law that goes singular or non-smooth as a state
approaches a physical bound — not just reaction kinetics. Population models
near extinction, transport models near a saturation limit, and control systems
near an actuator bound all show the same failure signature (clamp → spurious
frozen branch; naive continuation → blow-up) and the same fix.
