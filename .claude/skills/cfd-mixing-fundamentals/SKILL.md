---
name: cfd-mixing-fundamentals
description: CFD-based mixing analysis for stirred-tank-reactor (STR) compartment modeling — turbulence-model selection (k-epsilon variants, RSM, LES), MRF vs. sliding-mesh impeller modeling, CFD-derived Nq/Np/Q_ij extraction vs. correlation-based values, zone/compartment-boundary identification, RTD-from-CFD analysis, the macro/meso/micro mixing timescale hierarchy with Damkohler-number feed-point screening, practical snappyHexMesh/transient-scalar-transport pitfalls (STL unit scaling, non-manifold geometry, net-vs-gross flux, solver/timestep choice for transient species transport on a frozen flow field), and moving-mesh transient CFD for a growing liquid volume with a real inlet (prescribed mesh motion vs. VOF, inlet-patch under-resolution, the adjustPhi/pressure-reference trap for a closed domain with net inflow, and mesh-motion-diffusivity corner-shear failures). Use before constructing or revising a compartment model, before accepting a pumping-number correlation at face value for a non-standard (e.g. unbaffled) vessel, before assuming feed-point concentration gradients are negligible, before building/meshing a parametric CFD case or running transient scalar transport on a converged flow field, or before modeling liquid-volume growth/a real feed inlet via mesh motion.
---

# CFD-based mixing analysis for STR compartment modeling

<!-- Layout note: str_cfd/ is flat (no skills/ or knowledge/ subdirs), matching the existing
     str/ convention rather than cadet_chromatography_simulation's split skills/+knowledge/
     layout. Rationale: str_cfd/, like str/, is explicitly meant to be portable across more
     than one project (myr_simulation, str_mass_transfer) — the "Portability notes" section
     this file and mixing-time-correlations.md both carry is a cross-project convention that
     only makes sense in a flat, shared skill directory. cadet's split layout exists because
     its content is single-project and ingestion-heavy; that rationale doesn't apply here. -->

Source: Paul, Atiemo-Obeng, Kresta (eds.), *Handbook of Industrial Mixing* (Wiley, 2003),
Ch.5 "Computational Fluid Mixing," Ch.6 "Mechanically Stirred Vessels" §6-2.3/§6-3, Ch.9
"Blending of Miscible Liquids" §9-2/§9-5, Ch.13 "Mixing and Chemical Reactions" §13-1/§13-2/§13-4/§13-5.
Every formula and correlation below traces to a specific page in the extracted notes
(`/tmp/handbook_concepts.md` at authoring time) — items without a page-verified source are
marked `[TO VERIFY IN HANDBOOK]` rather than stated as fact.

## 1. When to use this skill

Consult before: constructing or revising a compartment model for any STR; accepting
Nq-correlation-derived pumping rates at face value for a non-standard vessel (unbaffled,
non-standard C/D, H/T far from 1, or an impeller/vessel combination outside a correlation's
tested range); deciding whether feed-point concentration gradients are negligible for
reaction selectivity. Read alongside `parameter-estimation.md` (compartment fitting
workflow — this skill supplies physics-based *priors* for that workflow's optimizer bounds
and compartment-count decision ladder, it does not replace either) and `biotech-modeling.md`
(scale-up judgment).

## 2. Turbulence modeling for STR CFD

- **Standard k–ε**: the default two-equation RANS model. Robust and cheap, but semi-empirical
  and tuned on high-Re flows; its dissipation-rate prediction is known to be unreliable enough
  that the Handbook recommends computing power draw from **blade-surface torque integration**
  (`P = 2πNτ`), never from integrating the turbulence model's own dissipation field, since
  different turbulence models can predict very different dissipation rates even when they agree
  on the mean flow pattern.
- **RNG k–ε**: modifies the dissipation equation for high-strain regions — flow around a bend
  or reattachment after a recirculation zone, i.e. exactly the region downstream of an impeller
  discharge jet. Extends acceptably into the transitional-Re range where standard k–ε is weaker.
  Recovers `Cµ≈0.0845` vs. the empirical 0.09 (within 7%) at high Re.
- **Realizable k–ε**: uses a variable `Cµ` (function of local strain/rotation) instead of the
  fixed 0.09, specifically to prevent unphysical negative normal stresses under high strain.
  Explicitly documented as **superior for predicting the spreading rate of round jets — including
  jets emitted from a rotating impeller blade**. Preferred over standard/RNG k–ε whenever the
  discharge-jet shape itself (not just bulk Nq/Np) is the thing being validated against
  experimental data.
- **RSM (Reynolds Stress Model)**: drops the isotropic-eddy-viscosity (Boussinesq) assumption
  entirely and solves the Reynolds stresses directly (6 extra transport equations in 3D).
  Explicitly recommended for **high swirl, rapid strain-rate change, and substantial streamline
  curvature — the Handbook names unbaffled stirred vessels by name as a case where RSM
  outperforms the k–ε family**. For an unbaffled vessel, RSM is the better-justified default,
  not standard k–ε — treat standard k–ε here as the fallback if RSM doesn't converge, not the
  first choice.
- **LES**: needed only to capture genuinely transient, sub-blade-passing-frequency behavior —
  low-frequency asymmetric flow instabilities (the slow "wobble" observed in real vessels) that
  no steady RANS model (any k–ε variant, RSM) can produce at all. Reserve for diagnosing an
  *observed* instability in fitted or experimental data, not as a routine tool for Nq/Q_ij
  extraction — it is transient and far more expensive than any RANS option above.
- **k–ω**: not covered by this Handbook edition; do not cite this skill file as a source for
  k–ω guidance.
- **Practical selection rule**: MRF steady-state Nq/Np extraction is not very sensitive to which
  RANS turbulence model is used — a validated benchmark (§3 below) shows all four RANS variants
  landing within ~10% of experimental Nq/Np on the same geometry. Turbulence-model choice matters
  far more for *local* accuracy (discharge-jet shape, dead-zone prediction, feed-point mixing)
  than for bulk pumping/power numbers. For an unbaffled vessel specifically, prefer realizable
  k–ε (cheap, good jet-shape fidelity) as a first pass and RSM as the accuracy check if the
  swirl/curvature at the vessel walls looks significant in the first-pass solution.

## 3. Impeller modeling: MRF vs. sliding mesh, and CFD-derived Nq vs. correlation Nq

**Impeller modeling.**
- **MRF (Multiple Reference Frames)**: steady-state; a rotating frame around the impeller, a
  stationary frame elsewhere, joined at a surface-of-revolution interface. The impeller's
  angular position relative to any asymmetric internal (baffle, probe, dip tube, feed line) is
  frozen for the whole solve — valid whenever that interaction is weak. **Good practice: solve
  MRF twice with the impeller at two different angular positions relative to the internal, and
  average the macroscopic results (Nq, Np)** rather than trusting one orientation. The
  "mixing-plane" MRF variant (circumferential averaging at the interface) is explicitly **not
  recommended for stirred tanks** — it produces unphysical results whenever inflow/outflow is
  not axisymmetric, which includes any single-point feed or dip-tube geometry.
- **Sliding mesh**: time-dependent; the grid around the impeller physically rotates. Required
  only for (a) resolving periodic impeller–internal interaction in the flow field itself, or
  (b) capturing/predicting low-frequency instabilities below blade-passing frequency — neither
  of which MRF can do. Initialize from a converged MRF solution rather than from rest, to skip
  the multi-revolution startup transient.
- **When MRF is sufficient**: for steady-state bulk-flow characterization — Nq, Np, zone
  boundaries, Q_ij — MRF is the appropriate, cheaper tool whenever impeller–internal interaction
  is weak. Escalate to sliding mesh only if a specific transient/periodic effect or instability
  is suspected in the data, not by default.

**CFD-derived Nq vs. correlation Nq.**
- **Extraction**: `Nq = Q/(N·D³)`, with `Q` obtained by integrating the outflow velocity across
  a discharge surface at the impeller — a **disk** normal to the shaft for an axial impeller, a
  **cylindrical shell section** for a radial impeller. `Np = P/(ρN³D⁵)`, with `P` from
  **blade-surface torque integration** (`P=2πNτ`), not from a volume integral of the turbulence
  model's own dissipation field (see §2).
- **When correlation Nq/Np break down** — each of these invalidates a standard correlation table
  and is grounds to extract Nq/Np from CFD instead: unbaffled or partially baffled vessel
  (published turbulent-plateau Nq tables are for baffled systems only — there is no unbaffled Nq
  table in this Handbook to fall back to); off-bottom clearance C/D far from standard
  (`Np ∝ (C/D)^-0.25` for a PBT); H/T far from 1; multiple impellers at nonstandard spacing
  (combined Np is *not* simply additive of single-impeller values); Re_i below ~1000 (the
  Handbook explicitly says not to trust turbulent-plateau Nq/Np charts below this).
- **Drop-in substitution**: `Q = Nq·N·D³` and `P = Np·ρN³D⁵` are the same functional forms
  whether Nq/Np come from a correlation table or a CFD flux/torque integral — substituting a
  CFD-derived value changes nothing about the downstream Re/Fr regime-check workflow in
  `parameter-estimation.md` §2, only the *provenance* of the number.

## 4. Zone identification and compartment boundary placement from CFD

- Visualization tools the Handbook uses for locating flow structure: velocity vectors,
  streamlines/stream function (2D), contours, isosurfaces (3D contour analog — used specifically
  to locate vortex cores via helicity isosurfaces), and the strain-rate-tensor modulus (used for
  identifying high-shear/dispersion regions).
- **Operational boundary rule**: the source text does not give a single named formula for "where
  to draw a compartment boundary" — the closest documented equivalent is reading the
  streamline/contour structure to locate the discharge jet, the recirculation loop(s) it drives,
  and any region the jet clearly does not reach (a candidate dead zone). Treating the locus of
  `v_z = 0` (time-averaged axial velocity sign change) as the boundary, as used in
  `parameter-estimation.md`, is a reasonable operational translation of this but is **not itself
  a verbatim Handbook rule** — `[TO VERIFY IN HANDBOOK Ch.6 §6-3, not extracted this pass]`.
- **Q_ij extraction**: same integration principle as Nq — integrate the net normal velocity
  across the zone-interface plane. The Handbook's own passage on this (§5-6.4.2) does not
  separately flag a net-vs-gross flux distinction; if the interface plane sees flow crossing
  both directions (recirculation straddling the boundary), decide explicitly whether the
  downstream compartment model wants net or gross exchange and document the choice —
  `[TO VERIFY]`.
- **Worked precedent for feed-point-vs-bulk zoning** (Bakker & Fasano 1993, using the Middleton
  et al. 1986 competitive-consecutive system A+B→R, R+B→S): moving the feed inlet from the bulk
  to directly above the impeller discharge — i.e. changing which zone the feed enters — changed
  by-product selectivity `X_s` by roughly a factor of 2 in a real CFD study, with no change to
  flow rate or RPM. This is direct precedent for treating a feed zone as a distinct compartment
  from the bulk on the strength of *where the feed enters relative to the impeller discharge*,
  independent of any Da-number calculation — feed-point geometry is itself evidence for or
  against a dedicated feed compartment, not just the reaction kinetics.

## 5. RTD analysis: linking CFD tracer simulation to compartment validation

- **Workflow**: converge the background-fluid flow field first (from experimental data or MRF).
  If the tracer's fluid properties match the background exactly, **disable momentum/continuity/
  turbulence equations** and solve only the transient scalar-transport equation for the tracer —
  cheaper and more numerically robust, since it isn't coupled to actively-changing variables.
  Track mean concentration and standard deviation across the vessel over time.
- **Two cases where the flow-field solve cannot be frozen during tracer transport**: (1) sliding
  mesh — the flow field is required at every time step regardless of the tracer; (2) the tracer
  enters through an inlet/dip tube carrying **significant momentum** (a real jet) — the flow
  solve must resume once that inlet becomes active, since the incoming stream perturbs the
  velocity field itself. Before assuming case (2) doesn't apply to a given feed line, check the
  feed-line's jet Reynolds number/momentum flux against the bulk circulation rather than assuming
  a low volumetric rate implies negligible momentum.
- **E(t)/F(t) formal definitions and the N-tanks-in-series formula belong to Handbook Chapter 1**
  (Residence Time Distributions), which was not extracted for this skill file — treat the
  existing E(t) formula already in `parameter-estimation.md` §2 as the authoritative source for
  that specific equation until Ch.1 is separately page-verified. `[TO VERIFY IN HANDBOOK Ch.1]`.
  Use RTD fitting (measured tracer/titration curve in hand) in preference to the eigenvalue
  cross-check in `parameter-estimation.md` §4 whenever experimental data is available; the
  eigenvalue approach remains the right tool when it isn't.

## 6. Macro, meso, and micro mixing: the three-scale hierarchy and Da screening

**Macro-mixing (θ_95).** Two independently-sourced correlations agree closely:
`Po^(1/3)·N·θ_95·D²/(T^1.5·H^0.5) = 5.20` (±10%, turbulent regime, Ch.9 eq. 9-1) and
`N·τ_B = 5.4·(T/D)²/Po^(1/3)` for Re>6400 (Ch.13 eq. 13-8, same underlying Grenville data,
restated). **Both are explicitly validated on baffled vessels only** — apply to an unbaffled
vessel only via the Fr-based correction path documented in `mixing-time-correlations.md` §4, not
directly. Converting between homogeneity thresholds is an *exact* log relation, not an
approximation: `θ_z = θ_95 · ln[(100−%homog)/100] / ln(0.05)` (e.g. `θ_99 = 1.537·θ_95`). From
CFD (no correlation needed): track max-minus-min tracer concentration spread vs. time in the
frozen-flow-field tracer simulation of §5.

**Meso-mixing.** *"Mesomixing effects most typically occur when the feed rate is greater than
the local mixing rate, allowing a plume of higher concentration to spread from the feed point"*
— this is the mechanism to check for any pumped/fed reagent, not just a Da-number check at a
single point. Two competing timescale estimates (Baldyga & Bourne):
`τ_D = Q_B/(U·D_t)` (`Q_B`=feed volumetric rate, `U`=local surrounding-fluid velocity at the feed
point, `D_t=0.1k²/ε`=local turbulent diffusivity) and
`τ_S = A·(L_s²/ε)^(1/3)` (`A`≈1–2, `L_s`=concentration macroscale set by feed geometry, e.g. feed
pipe diameter). Both are explicitly flagged in the source as limited — neither captures a
deliberately high-momentum feed jet correctly. **Empirical test for meso- vs. micro-mixing
control**: if yield/selectivity data from two different scales do *not* collapse onto one curve
when plotted against P/V, but *do* collapse (better, if imperfectly) when replotted against
impeller speed N, that is direct evidence of mesomixing control — this is the Handbook's own
worked diagnostic (Middleton et al. 1986 data), not a hypothetical. Minimum feed duration to
avoid a mesomixing-driven yield penalty is itself a mesomixing question — a critical-addition-time
relation `τ_crit·N^n = constant` was found empirically (Bourne & Hilber 1990); below the critical
addition time, selectivity depends on *absolute* impeller speed, not just local ε — meaning
constant-P/V scale-up alone is not sufficient once a feed is shorter than this critical duration.

**Micro-mixing.** Kolmogorov length `η=(ν³/ε)^(1/4)`; Batchelor length
`λ_B = η/√Sc = (ν·D_AB²/ε)^(1/4)` (for Sc≈1000, λ_B can be ~30× smaller than η — the relevant
scale for a typical aqueous solute is Batchelor, not Kolmogorov). Engulfment micromixing time
(Baldyga & Bourne): `τ_E = 17·(ν/ε)^(1/2)` — valid for Sc<~4000 (turbulent engulfment regime);
above that, viscous-stretching mixing dominates instead (different mechanism, not detailed here).
Corrsin's full mixing time `τ_M = (L_s²/ε)^(1/3) + 0.5·(ν/ε)^(1/2)·ln(Sc)` makes explicit that the
*mesoscale* term (`L_s`-dependent) usually dominates over the micromixing term for Sc≫1 liquids —
i.e. for most aqueous feed problems, the mesomixing contribution, not the micromixing one, is
where the real risk usually lives. Use `ε ≈ P/(ρV)` as a volume-averaged estimate only when a
spatially-resolved CFD field is unavailable — this specific framing is a standard identity
consistent with the Np/P relations in `mixing-time-correlations.md` §2, but was not found verbatim
in the extracted passages: `[TO VERIFY IN HANDBOOK]`.

**Da screening.** `Da_M = τ_M/τ_R` — the *mixing* Damköhler number (distinct from the
"traditional" Da = residence time/reaction time; don't conflate the two across sources).
Reaction time for the screening: `τ_R = 1/(k_R2·C_B0)`, with **`C_B0` = the feed-stream
concentration, not the bulk concentration** — the Handbook itself calls this "the worst
condition in the reactor," i.e. the deliberately conservative choice, confirming this is standard
practice rather than a project-specific convention. By-product selectivity rises continuously
(roughly log-linearly) with Da_M across ~3 orders of magnitude in the source data — Da is a
ranking/screening tool, not a sharp binary threshold. A single-point Da_M is itself a
simplification (real local mixing time varies as fed material moves through zones of different
ε) — the Handbook's own stated remedy for this, short of full spatially-resolved CFD, is exactly
the **zone/compartment model** approach (Patterson 1975) already in use for STR transport models
here — i.e. the compartment-modeling approach in `parameter-estimation.md` is not a workaround,
it is the Handbook's own recommended practical answer to Da_M's known point-estimate limitation.
Before accepting a well-mixed feed-zone assumption, compute Da_feed at the operating RPM range
using the feed-stream concentration and the local (or volume-averaged, if CFD isn't available) ε;
independently check the feed-pipe-exit-to-tip-speed velocity ratio (`v_f/v_t`, recommended
minimum ~0.5, some geometries need >2) as a second, purely-geometric red flag for backmixing —
a low ratio is grounds for a feed compartment even if the Da screening alone looks marginal.

## 7. Building parametric STL geometry for snappyHexMesh — pitfalls hit in practice

Lessons from actually building a parametric vessel+impeller mesh from scratch (tank
wall, torispherical head, shaft, blades) rather than starting from a hand-built
CAD/tutorial case — each of these produced a mesh that looked superficially fine
(`checkMesh` reported "Mesh OK", no solver error) but was silently wrong.

- **STL/case unit consistency.** A geometry script working in cm (a natural unit for
  vessel dimensions) writing STL vertices directly, dropped into an OpenFOAM case built
  in metres (`blockMeshDict`, `MRFProperties`, `locationInMesh`), makes the STL ~100x
  too large relative to the mesh domain. The background box then sits entirely *inside*
  the (apparently enormous) wall surface — nothing outside it is ever excluded, so
  `checkMesh` reports the *entire* background box as the fluid domain, with no error at
  any stage. Fix: `scale 0.01` (or whatever the actual ratio is) on each
  `triSurfaceMesh` entry's import in `snappyHexMeshDict`, not a unit conversion inside
  the STL-writing code itself (keeps the geometry script's native unit and the mesh
  case's unit independently sane, converting only at the one interface between them).
  **Diagnostic**: if the meshed domain's `checkMesh` bounding box/volume matches the
  *background box* dimensions rather than the intended vessel dimensions, suspect this
  first, before suspecting the castellation/snapping settings.
- **Non-manifold coincident vertices break the inside/outside flood fill.** A rotating
  shaft modeled as a surface-of-revolution cone down to a single point on the axis
  (instead of a constant-radius cylinder) put that point exactly coincident with the
  vessel bottom cap's own centre vertex — a different surface's vertex landing exactly
  on top of another. This kind of single-point non-manifold junction let
  `castellatedMesh`'s locationInMesh-based region split leak between what should be two
  disconnected regions. Model rotating shafts as genuine constant-radius cylinders (both
  end profile points at the same radius), never tapering to a point that could coincide
  with another surface's vertex.
- **Zero-thickness "sheet" obstacles are unreliable for a plain `wall`
  refinementSurface.** A flat impeller blade built as a single 2-triangle rectangle (no
  back face, no defined inside/outside) caused the mesh's inside/outside region count to
  fragment into tens of thousands of tiny regions mid-refinement, later merging back into
  one region that included the entire background box. Model any internal obstacle
  (blade, baffle meant as a real solid, not a genuine zero-thickness baffle patch) as a
  closed thin box (6 faces), not a single flat surface — a true zero-thickness baffle
  needs OpenFOAM's dedicated `createBaffles` mechanism, not a plain `wall`
  refinementSurface.
- **Mesh domain extent for a vessel with a curved bottom head must reach the true apex,
  not an intermediate reference plane.** Using the head/cylinder tangent line (a
  natural-seeming "where the straight wall starts" reference) as the mesh's lower z
  bound instead of the true geometric bottom silently clipped out nearly the entire head
  volume — the background box simply never extended down that far, so nothing there was
  ever castellated. **Diagnostic that catches this directly**: compare the meshed
  volume (`checkMesh`'s reported total volume) against the vessel's known liquid volume
  before trusting the mesh for anything — a ~10% shortfall for a vessel with a
  proportionally-sized head is exactly this bug's signature. This is the mass-balance-
  validation practice applied to geometry, before ever getting to a transport quantity.
- **Refinement levels don't port between scales at face value.** The same absolute
  surface-refinement levels, on the same background-block resolution (~40 cells across
  the vessel regardless of its physical size), gave a 1.1M-cell mesh for a 5 L vessel and
  a 3.7M-cell mesh for a geometrically similar 250 mL vessel — the smaller vessel's
  background cells are proportionally smaller in absolute terms, so the *same* level
  count over-resolves it. Scale refinement levels down (and/or the global cell-count cap)
  for smaller vessels rather than reusing a larger scale's settings verbatim; a rough
  rule that worked here was one fewer level across the board below roughly half the
  reference vessel's characteristic radius.

## 8. Transient scalar transport on a frozen MRF flow field — numerics and setup

Extends §5's frozen-flow-field shortcut (converge the flow once, disable
momentum/turbulence, solve only the scalar) with the concrete validation and numerics
that made it trustworthy and fast enough to actually run in practice.

- **Validate with the exact conserved quantity, every time, not just "it didn't
  crash."** Check the transported scalar's volume integral over the whole domain against
  (source strength × elapsed time) at multiple points during the run. This caught real
  setup mistakes early (see below) and, once passing to <1%, was strong enough evidence
  to trust the frozen-field shortcut's other outputs (probe time series, 3D snapshots)
  without further scrutiny.
- **Confirm the source region actually contains fluid cells before running, don't
  assume it.** A feed/source region (e.g. `fvOptions`' `scalarSemiImplicitSource`) that
  lands on a solid internal (a shaft occupying the vessel centreline) or straddles the
  domain boundary silently selects zero cells — the solver prints a `No cells selected!`
  warning and then crashes with a floating-point exception the moment it tries to divide
  the source strength by that region's (zero) volume. Query the mesh's actual cell
  centres directly (e.g. via a mesh-reading library) against the intended source
  geometry *before* running, not after hitting the crash. This is also the reason an
  on-axis feed point is often not a safe default in a stirred vessel: the impeller shaft
  usually occupies exactly that region.
- **Net vs. gross flux, worked concrete recipe** (extends §4's warning with the actual
  method, since the Handbook itself doesn't give one): for a closed, recirculating
  vessel, the *net* flux through an internal interface plane is ~zero by continuity
  (confirmed directly: computed net flux ~1e-5 relative to a resolved circulation of
  several L/min) — fluid crosses the plane in both directions as part of the
  recirculation loop. The compartment model's exchange flow `Q_ij` is the **gross**,
  one-directional flow: slice the domain at the interface plane, integrate *only* the
  positive (or only the negative) normal-velocity component over it, and use that as
  `Q_ij` — never the net integral a generic flux-integration function object reports by
  default.
- **Solver choice matters far more here than for a steady RANS solve.** A local smoother
  (GaussSeidel) needed on the order of 1000 sweeps per implicit step to converge a
  transient scalar equation on a ~500k-1M-cell mesh — a local smoother propagates
  information too slowly across a large mesh within one implicit step. Switching to a
  proper Krylov solver (PBiCGStab with a DILU preconditioner) cut that to ~20-100
  iterations at the same tolerance. Benchmark wall-time-per-simulated-second directly
  across a few candidate timesteps rather than assuming smaller is always more efficient
  — in one case a 2 s step was more wall-time-efficient than a 0.2 s step overall, since
  the Krylov iteration count did not grow proportionally with step size, despite the
  local Courant number reaching the thousands near the impeller. `adjustTimeStep` was
  tried first and behaved unreliably for this solver (jumped straight to a much larger
  step than the requested Courant target, repeatably) — prefer a benchmarked fixed step
  over trusting automatic step control blindly, at least until its behavior for the
  specific solver in use has been separately verified.
- **The transport diffusivity for a species riding on an already-resolved velocity field
  should come from the local turbulence field, not a lumped model's exchange-flow
  parameter.** `DT ≈ nut/Sc_t` (turbulent Schmidt number ~0.7, from the converged flow's
  own turbulent viscosity field) is the right scale for the sub-grid/molecular mixing
  still needed on top of a resolved velocity field. A lumped compartment model's
  `D_eff` (derived from an exchange flow rate, e.g. `D_eff = Q·Δz/A`) is typically an
  order of magnitude or more larger, because it was standing in for the *entire*
  unresolved advective exchange — reusing it once the real velocity field already
  provides that advection double-counts the same transport.

## 9. Moving-mesh transient CFD with a real inlet — modeling a growing liquid volume

Lessons from replacing a frozen-flow-field + interior-source-term setup (§8) with a
single transient run that resolves actual liquid-volume growth (a rising free surface)
and a real advective inlet, without going as far as a full VOF two-phase solve.

- **Prescribed mesh motion, not VOF, is the right first move when only the liquid phase
  matters.** If nothing about the air phase itself needs modeling (no interest in its
  velocity/pressure field, just where the liquid goes), a single-phase domain whose top
  boundary rises via a *prescribed* kinematic law (`dynamicMotionSolverFvMesh` +
  `displacementLaplacian`, driven by a known `V(t)`/`dV/dt`) is far cheaper and lower-risk
  than resolving two phases and a captured interface — at the cost of not modeling the
  actual free-surface shape (sloshing, meniscus). Reach for VOF only once the air phase's
  own behavior, or the interface shape itself, actually matters to the question being asked.
- **Don't pre-mesh headroom the moving-mesh case doesn't need.** A tempting-looking shortcut
  — extend the initial mesh domain up into the vessel's reserved freeboard "so the rising
  boundary has somewhere to go" — is wrong for this approach specifically: prescribed
  mesh motion *deforms existing cells*, it doesn't add new ones, so the mesh only needs to
  start at the true initial liquid height and stretch a few cm over the run. Extending it
  upfront silently solves the whole freeboard volume as liquid from t=0, inflating the
  starting liquid volume by whatever the freeboard fraction adds (this is exactly the kind
  of thing a VOF air-phase region would need, but a single-phase prescribed-motion domain
  does not). **Diagnostic**: compare mesh volume against the *known initial* liquid volume
  (not the final one) at t≈0 — a many-percent excess right at the start, not building up
  over time, is this bug specifically (contrast with real transport losses, which usually
  show up as a trend, not an offset present from the first timestep).
- **A real inlet patch's actual meshed area can differ substantially from its nominal
  geometric area at coarse refinement — verify it, don't assume it.** A small feed-inlet
  patch sized and given a `fixedValue` velocity based on its intended geometric
  area/flow-rate came out, once actually meshed at a coarsened refinement level, at only
  ~36% of that nominal area — silently injecting the tracked species at ~36% of the real
  rate, with no error or warning anywhere (the BC values themselves were exactly as
  specified; only the *count of faces actually carrying them* was short). **Diagnostic**:
  after meshing, read the actual patch back (e.g. via a mesh-reading library, summing face
  areas) and compare to the geometric area used to derive the prescribed velocity — this
  is a direct, cheap check worth doing for *any* small/localized boundary feature on a
  coarsened mesh, and it's exactly the same category of check as verifying total meshed
  volume against known liquid volume (§7), just applied to a boundary patch instead of the
  whole domain. **Fix**: refine just that one small feature's surface level independently
  of the rest of the (deliberately coarsened) mesh — refining a tiny, localized patch barely
  changes total cell count, unlike refining a large region such as the impeller zone.
- **A closed domain with a genuine net inflow (no real outlet) needs a real pressure
  reference — this is a distinct issue from moving vs. static mesh.** With every boundary
  patch's velocity fully prescribed (no-slip walls, rotating walls, a fixed-value inlet)
  and every pressure patch left `zeroGradient`, OpenFOAM's `adjustPhi` continuity check
  (triggered whenever `p.needReference()` is true, i.e. no patch fixes an absolute pressure
  value) has no way to reconcile a real net inflow against a domain with no outlet, and
  hard-fails on the very first timestep ("Continuity error cannot be removed by adjusting
  the outflow") — **regardless of whether the mesh is moving**, and this failure mode is
  easy to misdiagnose as a moving-mesh-specific problem when it isn't. The fix is to give
  **one** patch (here, the rising top boundary) a genuine `fixedValue` pressure (physically
  reasonable for a free surface, which is close to atmospheric/kinematic-zero pressure) —
  once `p` has a real reference, the *velocity* type at that same patch can still correctly
  be a true moving no-slip wall (`movingWallVelocity`, zero relative flow), which is both
  the physically correct choice (100% of the volume growth is then correctly attributed to
  the real inlet, matching reality: nothing actually flows out of a rising, closed vessel)
  and avoids a subtler trap: a first attempt used a bidirectional "vent" velocity BC
  (`pressureInletOutletVelocity`) specifically to give `adjustPhi` a patch to balance
  against — it ran without error, but let species mass genuinely leak out through whatever
  fraction of that vent had locally-outward flow at any given moment (confirmed directly:
  a mass-balance check that should read ~100% instead read ~35-45% and was *still falling*
  over time, not a one-time startup transient). A working fix must both run without error
  *and* pass the same mass-balance validation as everything else in this skill file —
  "no error" alone did not mean "correct" here.
- **`correctPhi yes` is required in `PIMPLE{}` for any `dynamicFvMesh` case** — without it,
  the flux isn't corrected for the mesh's own motion at the first sub-step of every mesh
  update. This is a one-line, easy-to-miss setting distinct from the pressure-reference fix
  above; both were needed together.
- **Mesh-motion diffusivity can concentrate deformation right at a rigid/moving patch
  junction, tearing a cell there — this one was not fully solved.** Using
  `inverseDistance(<moving patch>)` diffusivity (to protect an already-refined region
  elsewhere in the mesh from motion-related distortion) concentrates nearly all the
  prescribed deformation into the single cell layer nearest that moving patch. Exactly at
  the corner where that moving patch meets a *fixed* (zero-displacement) patch, this
  produced a severely sheared cell and a turbulence-field blowup (`epsilon` diverging past
  1e6, eventually 1e32) within about one second of simulated time for one vessel's mesh,
  and a much-delayed but likely related stall for a different vessel's mesh at the same
  scaling of that corner. Two follow-up attempts — lowering the Courant cap plus adding a
  non-orthogonality corrector, then switching to spatially-`uniform` diffusivity (spreading
  deformation over the whole domain instead of concentrating it) — each ran further but
  neither fully resolved it. Per this project's diagnostic-first practice: after a second
  fix attempt didn't obviously work, this was surfaced as an open tradeoff rather than
  guessed at a third time. **If revisiting this**: a spatially-varying custom diffusivity
  that specifically tapers deformation to zero right at a rigid-patch/moving-patch shared
  edge (rather than either a globally uniform or a single-patch-referenced
  inverse-distance field), or local mesh refinement concentrated at that specific edge, or
  accepting VOF (which has no such edge at all, since the interface itself is resolved
  rather than a kinematic constraint imposed on the mesh) are the more promising directions
  — not a third guess at a stock diffusivity keyword.
- **Wall-clock feasibility for transient MRF + mesh motion is a real, separate constraint
  from getting the physics right.** A mesh sized for steady-state accuracy (§7's guidance)
  can be 2-3 orders of magnitude too expensive for a transient run needing thousands of
  timesteps — a first attempt at steady-state-equivalent refinement measured multiple
  *days* of wall-clock to reach a feed duration of ~145 s of simulated time. Coarsening
  every refinement level (roughly by half) and the background mesh cut cell count ~6x and
  improved throughput ~40x. Benchmark actual physical-seconds-per-wall-clock-second
  directly on the coarsened mesh before committing to a multi-hour run, the same way §8
  already recommends benchmarking solver/timestep choice — don't assume a scaling factor.

## 10. Portability notes

**Changes between projects:** vessel geometry (baffled/unbaffled, H/T, C/D), impeller type and
pumping direction, RPM range, feed-stream concentration relative to bulk, feed-line velocity
relative to tip speed, and which specific turbulence model is warranted (unbaffled vessels shift
the recommendation toward RSM over standard k–ε; unbaffled also removes access to standard
baffled-vessel Nq/θ_95 correlation tables, forcing either CFD extraction or a Fr-based correction).

**Stays the same across projects:** the three-scale mixing hierarchy (macro/meso/micro) and the
diagnostic that mesomixing — not micromixing — is usually the dominant risk for a slowly fed
reagent in an aqueous (Sc≫1) system; the empirical meso-vs-micro diagnostic (replot yield against
N, not P/V, and see which collapses multi-scale data); Da_M evaluated at feed-stream (not bulk)
concentration as the conservative, standard convention; treating compartment/zone models as the
Handbook-endorsed practical substitute for full spatially-resolved CFD, not an approximation to
apologize for; MRF as the default steady-state tool for Nq/Np/Q_ij extraction, escalating to
sliding mesh only for a specific suspected transient/periodic effect; and using CFD-derived
Nq/Np/Q_ij as a drop-in replacement wherever a standard correlation's validity conditions
(baffling, C/D, H/T, impeller spacing, Re_i) are not met, without changing anything else about
the downstream Re/Fr-regime-check or compartment-fitting workflow; and defaulting to an
**executed Jupyter notebook** (numerical results and plots baked in as actual cell outputs,
not just described in prose) as the deliverable for CFD/simulation work — a written summary
alone is not the expected format for this class of task.
