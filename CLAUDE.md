# myr_simulation

Mechanistic modeling of the insulin desB30 -> Detemir myristoylation reaction
in a 5 L stirred-tank reactor: mass-transfer/mixing characterization,
kinetic parameter fitting, and reactor-scale predictive simulation.

## Environment

- Conda environment: `ph_adjustment_env` — activate with
  `conda activate ph_adjustment_env`.
- Working directory: `/mnt/databiotek/Surya/Code/Python/myr_simulation`
  (also reachable as `Z:\Surya\Code\Python\myr_simulation` from Windows).
- Python: environment reports 3.9.23, but notebook kernel metadata says
  3.10.13 — this mismatch has not been resolved [verify] before assuming
  which one governs a given notebook run.
- Key packages present in the env: numpy 2.0.2, scipy 1.13.1, pandas 2.3.3,
  matplotlib 3.9.4, dill 0.4.1. The `webapp/` subproject has its own
  `requirements.txt` (streamlit, plotly, cryptography) and is meant to run
  in its own virtualenv, not this conda env.

## Project structure

- `kinetika3b.ipynb` — 2-compartment (feed/probe zone vs. bulk) titration-
  curve fit; establishes the Re/Fr regime analysis and the theoretical
  Nq/Np/Q/P scale-up framework, fit at 450 RPM.
- `5L_tank_transfer_mass_3_compartments.ipynb` — 3- (and, as an extension,
  4-) compartment titration fit at 160 RPM; the current validated source of
  fitted inter-zone flow rates `Q_1->2`, `Q_2->3` for the 5 L vessel.
- `DS_Pool_Kinetic_Model_Refined.ipynb` — the kinetic parameter fit that
  produces the active `.pkl` (see below).
- `myr_mechanistic_model_5L_optimized_FIXED.ipynb` — the current 2-
  compartment mechanistic reactor model (smooth-floor fix applied); source
  of the SOP operating point (see Hard rules); the Da_feed screening from
  `cfd-mixing-fundamentals.md` §6 should be confirmed at 160 and 450 RPM
  before extending the well-mixed feed-zone assumption to operating points
  outside the validated SOP.
- `myr_mechanistic_model_5L_3_compartments.ipynb` — 3-compartment version
  of the same reactor model, reusing the fitted `Q` values above; finding:
  compartments are numerically indistinguishable from well-mixed for this
  reaction's timescale.
- `transfer_mass_5L_pde.ipynb` — 1D axial reaction-diffusion PDE version of
  Stage 2 (dosed NHS-Myr feed + reaction).
- `scale_down_250mL.ipynb` — 250 mL scale-down transfer-mass/mixing model
  (mixing physics only, no reaction kinetics), torispherical head geometry.
- `pra5_myr_model.ipynb` — forward-predictive (non-optimized) simulation
  against two real batch datasets (Run 1 / Run 2).
- `dataprep_rx_kinetic.ipynb` — experimental data preparation feeding the
  kinetic fitting notebooks.
- Older/superseded notebooks (`Kinetika.ipynb`, `kinetika2.ipynb`,
  `HydrolysisKineticModel*.ipynb`, `myr_mechanistic_model_5L.ipynb`,
  `myr_mechanistic_model_5L_optimized.ipynb`, `myr_mechanistic_model_5L_refined_kinetics.ipynb`,
  `DS_Pool_Kinetic_Model.ipynb`) are kept for history; do not build on them
  without checking whether a newer notebook already supersedes them.
- Portable skills directory:
  `/mnt/databiotek/Surya/Code/claude_portable_skills/str` — read the
  relevant skill file there before doing ODE-debugging or parameter-
  estimation work (see "Skill files in use" below).
- `webapp/` is a standalone Streamlit tool with its **own** ODE-simulation
  and least-squares-fitting engine (`webapp/engine/transport_model.py`,
  `model_selection.py`) that independently fits 2/3/4-compartment
  mass-transfer models to uploaded titration CSVs and recommends a
  compartment count. It does not read outputs from any notebook in this
  repo and no notebook reads its outputs. Do not change its CSV input
  format or its output/result schema without being explicitly asked.

## Kinetic model in use

- Active fit: `02jul2026_ds_pool_refined_results.pkl` — DS-Pool refined
  kinetic fit: 6 lumped rate constants + free reaction orders (n=0.119 on
  NHS-Myr, m=1.89 on substrate), R²=0.774. Loaded via `dill`, with rate
  constants stored as `log(k)` (un-transform with `exp()` at load time).
- Species pooling in use (to match what this fit resolves): `other_DS1` =
  A1 + B1, `DS2_total` = A1_B1 + A1_B29 + B1_B29, `DS3` = A1_B29_B1.
- `02jul2026_ds_pool_results.pkl` is an earlier (non-refined) version of the
  same fit — not the active one.
- Do **NOT** revert to the older 9-species full-regiochemistry `.dill` fits
  (`02jul2026_optimized_hydrolysis_kinetic_model.dill`,
  `30jul2026_hydrolysis_kinetic_model.dill`) as the kinetics source for any
  reactor model — the pooled DS-Pool `.pkl` fit above is what current
  reactor notebooks are built on.

## Hard rules (non-negotiable, never override)

- **Smooth-floor rule.** Never clamp concentrations near zero with a hard
  `np.maximum(0, c)` before a fractional-power rate-law term — it traps
  stiff solvers on a spurious frozen branch. Use the smooth non-negative
  floor (`0.5*(c + sqrt(c^2+eps^2))`) instead. See the `ode-stiff-kinetics`
  skill for the full pattern.
- **Mass balance validation requirement.** Before treating any numerical
  fix to a reactor ODE/PDE as done, verify it against an exact conserved
  quantity (e.g. total protein/mass balance) — not just "it no longer
  crashes."
- **The SOP operating point is fixed.** NHS-Myr stock 0.482 mg/mL, pump
  flow ~12.96 mL/min, pump stopped at t=30.0 min (predicted 64.8% Detemir
  yield, 71.9% purity). Do not re-optimize or shift this operating point
  unless explicitly asked to.
- **Diagnostic-first rule.** When a mechanistic ODE/PDE model produces
  implausible output (mass-balance violation, a state past a physical
  bound, solver stall/blow-up), root-cause it with concrete evidence (raw
  solver-node trajectories, step sizes, conserved quantities) before
  applying a fix — do not patch blindly. If a second fix attempt doesn't
  obviously work, stop and surface the tradeoff instead of trying a third.
  See the `ode-stiff-kinetics` skill.
- **Da-number feed-point screening rule.** Before treating the NHS-Myr feed
  zone as kinetically equivalent to the bulk, compute Da_feed =
  t_micro / t_reaction_local at the feed point using the Kolmogorov estimate
  from `mixing-time-correlations.md` §3 and the active kinetic fit
  (`02jul2026_ds_pool_refined_results.pkl`) evaluated at feed-stream
  concentration (not bulk concentration). If Da_feed > 0.1, a separate
  feed-point compartment with elevated local concentration is warranted.
  This check is mandatory whenever: (a) pump flow rate or NHS-Myr stock
  concentration changes from the current SOP, (b) RPM moves outside the
  validated 160–450 RPM range, or (c) a new compartment count is being
  explored. Document the computed Da_feed value alongside any simulation
  result produced outside the current SOP conditions.

## Skill files in use

The Paul et al. (2003) Handbook of Industrial Mixing (North American Mixing
Forum, Wiley-Interscience) is the primary reference for mixing hydrodynamics
in this project. Key content from Chapter 5 (Computational Fluid Mixing) and
the mixing-time/correlations chapters has been extracted into the two
`str_cfd` skill files below.

**At session start, read all skill files listed below before doing any work.**
Located at `/mnt/databiotek/Surya/Code/claude_portable_skills/`:
Skills in `str/` (bioprocess judgment, parameter estimation, ODE numerics)
and `str_cfd/` (CFD mixing, RTD analysis, Da screening):

- `ode-stiff-kinetics.md` — diagnosing and fixing stiff/singular ODE or PDE
  kinetics (mass-balance violations, solver blow-up/stall, fractional
  reaction orders near a depletion boundary). Read before touching any
  reactor ODE/PDE solve in this repo.
- `parameter-estimation.md` — methodology for fitting rate constants and
  mass-transfer/compartment flow rates to experimental data, including how
  to determine and validate STR compartment count, fit-quality/
  identifiability diagnostics, and how to hand a fitted parameter set off
  to a downstream mechanistic model. Read before any new kinetic or
  transport-parameter fitting work, or before choosing a compartment count
  for a new STR model.
- `biotech-modeling.md` — read before any kinetic modeling or scale work.
- `cfd-mixing-fundamentals.md` — CFD-based zone identification, RTD analysis,
  macro/meso/micro mixing timescale hierarchy, Damköhler number screening for
  the NHS-Myr feed point, and how to derive compartment exchange rates Q_ij
  from CFD rather than from a pumping-number correlation. Read before any new
  compartment model construction, before accepting Nq-correlation-based Q
  estimates for the unbaffled 5 L vessel, and before concluding that feed-point
  concentration gradients are negligible for the reaction.
- `mixing-time-correlations.md` — Reference equations for macro-mixing time
  θ_95, power number Np, pumping number Nq (with Froude correction for the
  unbaffled case), and Kolmogorov/Batchelor micro-mixing timescale estimates.
  Use for Da_feed screening and to cross-check CFD-derived or fitted Q values
  against order-of-magnitude correlation predictions.
