---
name: cfd-mixing-fundamentals
description: CFD-based mixing analysis for stirred-tank-reactor (STR) compartment modeling — turbulence-model selection (k-epsilon variants, RSM, LES), MRF vs. sliding-mesh impeller modeling, CFD-derived Nq/Np/Q_ij extraction vs. correlation-based values, zone/compartment-boundary identification, RTD-from-CFD analysis, and the macro/meso/micro mixing timescale hierarchy with Damkohler-number feed-point screening. Use before constructing or revising a compartment model, before accepting a pumping-number correlation at face value for a non-standard (e.g. unbaffled) vessel, or before assuming feed-point concentration gradients are negligible.
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

## 7. Portability notes

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
the downstream Re/Fr-regime-check or compartment-fitting workflow.
