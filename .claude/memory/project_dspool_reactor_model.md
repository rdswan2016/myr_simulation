---
name: project-dspool-reactor-model
description: Status of the DS-Pool kinetics integration into the 5L mechanistic reactor model and the resulting SOP operating point
metadata: 
  node_type: memory
  type: project
  originSessionId: 99165efa-e454-4d46-833e-8127ba9c7a8c
---

`myr_mechanistic_model_5L_optimized_FIXED.ipynb` was rebuilt (2026-07-02) to use the DS-Pool refined kinetic fit
(`02jul2026_ds_pool_refined_results.pkl`: 6 lumped rate constants + free reaction orders n=0.119 on NHS-Myr,
m=1.89 on substrate, R^2=0.774) instead of the older 9-species full-regiochemistry `.dill` fit. Species had to be
pooled to match what the fit actually resolves: `other_DS1`=A1+B1, `DS2_total`=A1_B1+A1_B29+B1_B29, `DS3`=A1_B29_B1.

**Why:** user wanted the reactor model to simulate with the "current" kinetics and find pump speed / NHS-Myr stock
conditions giving a Detemir peak at ~30 min with minimal side products, subject to the shared NHS-Myr/glycine-quench
peristaltic pump range (10-200 mL/min) and 5 L reactor volume ceiling.

**Key finding:** the fitted NHS-Myr reaction order (n=0.119) makes the rate law nearly singular as NHS-Myr
concentration approaches zero (which happens constantly — natural depletion and by design during the glycine
quench). A hard `np.maximum(0, c)` clamp on concentrations traps stiff ODE solvers (BDF/Radau) on a spurious
flat/frozen branch (protein mass balance still summed correctly but individual species went unphysically negative,
e.g. desb30 stuck at -0.54 mmol against a 0.54 mmol charge). A sign-preserving continuation caused runaway blow-up
instead (species reaching >1 mol from a ~0.0005 mol system) because of the large rate constants (k up to ~13,700).
The fix that worked: a **smooth non-negative floor** `0.5*(c + sqrt(c^2 + eps^2))` with eps=1e-9 M applied before
every fractional-power term, combined with a terminal ODE event that stops integration once total NHS-Myr drops
below 0.01% of the charge (the true asymptotic decay toward zero is otherwise computationally unbounded). This
gives exact mass balance (protein conserved to 6 decimals) and fast (~1s/sim), stable integration.

**Chosen SOP operating point** (from the Section 8 sweep + root-find in the notebook): keep NHS-Myr stock
unchanged at 0.482 mg/mL (diluting further only bought ~1 selectivity percentage point, within the fit's noise),
set pump flow rate to ~12.96 mL/min (T_feed≈56.3 min nominal, pump stopped early at t=30.0 min, delivering ~53%
of the planned 2.0x charge ≈1.06x actual excess). Predicted outcome: 64.8% Detemir yield, 71.9% purity, batch
volume 3.475 L (well under the 5 L limit).

**How to apply:** if this kinetic fit or reactor model is revisited, don't reintroduce hard clamping on
concentrations near zero — use the smooth-floor + terminal-event pattern. See [[feedback_numerical_debugging]].

**2026-07-10 — spatially-resolved PDE version added:** built `transfer_mass_5L_pde.ipynb`, a 1D axial
(vertical) reaction-diffusion PDE version of Stage 2 (dosed NHS-Myr feed + reaction only; Stage 1/3 not
spatially resolved), reusing `reaction_derivs_mM_s` unchanged and replacing the 2-compartment `Q`-exchange
with a derived axial eddy-dispersion coefficient `D_eff = Q*Δz_centroid/A_tank ≈ 7.0e-5 m²/s`, matching the
box-exchange flux to an equivalent Fickian flux. Solved via method-of-lines (41 grid points, `solve_ivp`
`Radau`, ~75s wall clock on this machine). Confirmed the negative-NHS artifact above is present but *smaller*
on the finer 41-point grid (min ≈ -0.015 mM) than in the original coarse 2-box model (min ≈ -0.036 mM,
re-verified by rerunning the original ODE standalone) — i.e. finer spatial resolution doesn't make this
worse. PDE's domain-averaged Detemir peak lands at t≈31.5 min vs. the lumped model's 30.0 min, a good
macroscopic consistency check. Known limitations documented in the notebook itself: fixed domain length
(doesn't track the ~24% volume growth during dosing), `D_eff` is derived not measured, no grid-refinement
study done yet, Stage 3 (glycine quench) not included.

**Known residual artifact (found 2026-07-03, fixed for display):** the NHS-Myr trace in the Section 11 plot dips a
few uM below 0 mM. Confirmed real at raw (un-interpolated) solver nodes in both Stage 2 (~t=5 min, -0.013 mM) and
Stage 3 (drifts to ~-0.014 mM by end of the 30-min hold) — not a dense-output plotting artifact. Root cause:
`smooth_floor()` regularizes the rate-law *evaluation* but not the ODE state itself, and because every reaction
consuming NHS-Myr yields a rate >=0, `d(nhs_myr)/dt <= 0` everywhere — so once a stiff-solver step overshoots past
zero (driven by the near-singular curvature from n=0.12), there's no restoring force and it drifts further
negative, especially in Stage 3 which (unlike Stage 2) has no terminal "NHS exhausted" event. Fixed by clamping
only the plotted NHS-Myr line to `np.clip(..., 0.0, None)` in the Section 11 plotting cell — the underlying state,
mass balance, yield (64.8%), and purity (71.9%) numbers are computed from the raw (unclamped) state elsewhere and
are unaffected. If revisited, a real dynamical fix (terminal event in Stage 3, or a smooth restoring term) was
considered but not applied — user chose the display-only fix.
