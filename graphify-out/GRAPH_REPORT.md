# Graph Report - .  (2026-07-27)

## Corpus Check
- 23 files · ~54,029 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 167 nodes · 255 edges · 24 communities (16 shown, 8 thin omitted)
- Extraction: 87% EXTRACTED · 13% INFERRED · 0% AMBIGUOUS · INFERRED: 32 edges (avg confidence: 0.79)
- Token cost: 165,000 input · 12,829 output

## Community Hubs (Navigation)
- Project Notebooks & Modeling References
- Webapp Engine Architecture
- Chemistry & Transport Model Internals
- DS-Pool Reactor Debugging & Fixes
- DS-Pool Kinetic Fit & SOP
- Geometry & Zone Layout Module
- Model Selection (AIC/BIC/F-test)
- CFD Zone ID & Mixing-Time Studies
- Nq vs RPM Chart (Vortex Correction)
- Chemistry Params & Mass-Balance Ceiling
- CSV I/O Utilities
- Compartment-Count Decision Ladder
- TLS Cert Generation
- Failed Fix Attempts (Clamp/Sign-Preserving)
- Macro-Mixing Time Correlations
- Moving-Mesh CFD Pitfalls
- Micro-Mixing Timescales
- Impeller Modeling (MRF vs Sliding Mesh)
- Power & Pumping Numbers
- Webapp Engine Init
- Identifiability Diagnostics
- LES Turbulence Model
- RNG k-epsilon Turbulence Model
- Standard k-epsilon Turbulence Model

## God Nodes (most connected - your core abstractions)
1. `myr_simulation project CLAUDE.md` - 23 edges
2. `run_model_selection()` - 14 edges
3. `FeedParams` - 13 edges
4. `ChemistryParams` - 11 edges
5. `webapp/README.md - Mass-Transfer & Kinetic Fitting Tool` - 11 edges
6. `simulate_ph()` - 10 edges
7. `FitResult` - 8 edges
8. `fit_n_zone_model()` - 8 edges
9. `DS-Pool Reactor Model Status` - 8 edges
10. `cfd-mixing-fundamentals skill (CFD-based STR mixing analysis)` - 8 edges

## Surprising Connections (you probably didn't know these)
- `webapp/README.md - Mass-Transfer & Kinetic Fitting Tool` --semantically_similar_to--> `parameter-estimation skill (fitting rate constants and transport parameters)`  [INFERRED] [semantically similar]
  webapp/README.md → .claude/skills/parameter-estimation/SKILL.md
- `9-species regiochemistry + hydrolysis kinetics network` --semantically_similar_to--> `Species pooling convention (other_DS1, DS2_total, DS3)`  [INFERRED] [semantically similar]
  myr_mechanistic_model_5L_optimized.html → CLAUDE.md
- `Species resolution decision criteria (pooled vs. full regiochemistry)` --semantically_similar_to--> `Species pooling convention (other_DS1, DS2_total, DS3)`  [INFERRED] [semantically similar]
  .claude/skills/biotech-modeling/SKILL.md → CLAUDE.md
- `Stage 3: dynamic glycine quench` --semantically_similar_to--> `SOP Operating Point (0.482 mg/mL, 12.96 mL/min, stop t=30 min, 64.8% yield, 71.9% purity)`  [INFERRED] [semantically similar]
  myr_mechanistic_model_5L_optimized.html → CLAUDE.md
- `Forward batch-validation workflow (freeze params, run forward)` --semantically_similar_to--> `pra5_myr_model.ipynb (forward-predictive validation vs. Run 1/Run 2)`  [INFERRED] [semantically similar]
  .claude/skills/biotech-modeling/SKILL.md → CLAUDE.md

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Portable str_cfd/str skill set for STR mixing and kinetics modeling** — claude_skills_biotech_modeling_skill_biotech_modeling, claude_skills_cfd_mixing_fundamentals_skill_cfd_mixing_fundamentals, claude_skills_mixing_time_correlations_skill_mixing_time_correlations, claude_skills_ode_stiff_kinetics_skill_ode_stiff_kinetics, claude_skills_parameter_estimation_skill_parameter_estimation [EXTRACTED 1.00]
- **DS-Pool reactor model numerical debugging fix sequence** — hard_clamp_fix_attempt, sign_preserving_continuation_fix_attempt, smooth_floor_fix, terminal_ode_event, claude_memory_feedback_numerical_debugging_principle [EXTRACTED 1.00]
- **Macro/meso/micro mixing timescale hierarchy for Da screening** — macro_mixing_theta95, meso_mixing, micro_mixing_kolmogorov, da_m_screening, da_feed_number [EXTRACTED 1.00]
- **Vortex-suppressed Nq falling below literature PBT range, calibrated to single experimental anchor** — grafik_nq_vs_rpm_nq_vortex_corrected, grafik_nq_vs_rpm_literature_pbt_range, grafik_nq_vs_rpm_nq_experimental_cell3, grafik_nq_vs_rpm_unbaffled_vessel [INFERRED 0.85]

## Communities (24 total, 8 thin omitted)

### Community 0 - "Project Notebooks & Modeling References"
Cohesion: 0.10
Nodes (31): 5L_tank_transfer_mass_3_compartments.ipynb (3/4-compartment fit, 160 RPM), Forward batch-validation workflow (freeze params, run forward), Hard rules (smooth-floor, mass-balance validation, fixed SOP, diagnostic-first, Da-feed screening), myr_simulation project CLAUDE.md, claude_portable_skills source repo (str/ + str_cfd/), biotech-modeling skill (bioprocess modeling judgment), cfd-mixing-fundamentals skill (CFD-based STR mixing analysis), mixing-time-correlations skill (reference formulas) (+23 more)

### Community 1 - "Webapp Engine Architecture"
Cohesion: 0.16
Nodes (15): Henderson-Hasselbalch <-> transport ODE coupling (moles-based, single closure on probe zone), Mass-balance ceiling check (closed-system equilibrium vs. data peak), N-compartment ODE engine (moles-based state, Fickian junction exchange), Nested extra-sum-of-squares F-test for compartment-count selection, Structural cross-check: conserved-quantity ceiling (independent of optimizer), webapp/app.py (Streamlit UI), webapp/engine/chemistry.py (Henderson-Hasselbalch + titrant feed model), webapp/engine/geometry.py (vessel/zone volume construction) (+7 more)

### Community 2 - "Chemistry & Transport Model Internals"
Cohesion: 0.19
Nodes (12): FeedParams, initial_split(), pH_from_moles(), ndarray, Inverse Henderson-Hasselbalch: fraction of total analyte in the CONJUGATE-ACID, Henderson-Hasselbalch: pH = pKa + log10([base]/[acid]). Concentrations cancel to, fit_n_zone_model(), ndarray (+4 more)

### Community 3 - "DS-Pool Reactor Debugging & Fixes"
Cohesion: 0.23
Nodes (12): Numerical Debugging Feedback: diagnostic-first before patching, Project Memory Index (MEMORY.md), DS-Pool Reactor Model Status, D_eff derived axial eddy-dispersion coefficient (~7.0e-5 m2/s), Diagnostic-first, not patch-first workflow, myr_mechanistic_model_5L_optimized_FIXED.ipynb (current 2-compartment reactor model), Residual negative NHS-Myr display artifact (display-only clip fix), Smooth non-negative floor regularization (accepted fix) (+4 more)

### Community 4 - "DS-Pool Kinetic Fit & SOP"
Cohesion: 0.18
Nodes (12): DS_Pool_Kinetic_Model_Refined.ipynb (produces active kinetic fit), DS-Pool Refined Kinetic Fit (02jul2026_ds_pool_refined_results.pkl), Fit versioning / staleness discipline, Log-parameterize + dilled-bundle-with-fallback handoff pattern, Stage 3: dynamic glycine quench, myr_mechanistic_model_5L_optimized.ipynb (superseded, 9-species regiochemistry), 9-species regiochemistry + hydrolysis kinetics network, SOP Operating Point (0.482 mg/mL, 12.96 mL/min, stop t=30 min, 64.8% yield, 71.9% purity) (+4 more)

### Community 5 - "Geometry & Zone Layout Module"
Cohesion: 0.29
Nodes (8): Standalone web application: automated multi-compartment mass-transfer parameter, build_zone_layout(), probe_zone_check(), Vessel/zone geometry for the N-compartment mass-transfer model.  Zone-splitting, True if the pH probe height falls inside Zone 0 (the feed zone), which is the, Build zone boundaries/volumes for n_zones in {2, 3, 4}, reusing the exact     st, VesselGeometry, ZoneLayout

### Community 6 - "Model Selection (AIC/BIC/F-test)"
Cohesion: 0.31
Nodes (10): _aic_bic(), _build_rationale(), ModelSelectionResult, _nested_f_test(), ndarray, Automatic compartment-count selection.  Fits the 2-, 3-, and 4-compartment model, run_model_selection(), FitResult (+2 more)

### Community 7 - "CFD Zone ID & Mixing-Time Studies"
Cohesion: 0.25
Nodes (8): Bakker & Fasano 1993 feed-point-vs-bulk zoning CFD study, Bourne & Hilber 1990 critical addition time relation (tau_crit*N^n=const), CFD-to-ODE coupling interface (fixed Q_ij / zone-averaged epsilon, no tight coupling), Corrsin full mixing time tau_M, Meso-mixing (feed-plume spreading, tau_D/tau_S timescales), Middleton et al. 1986 competitive-consecutive reaction system (A+B->R, R+B->S), RTD-based experimental validation (E(t), N-tanks-in-series), CFD zone/compartment boundary identification (v_z=0 locus, Q_ij flux integration)

### Community 8 - "Nq vs RPM Chart (Vortex Correction)"
Cohesion: 0.48
Nodes (7): Nq vs RPM: Vortex-Froude Theory Chart (Cell 3 anchored), Impeller Froude Number (Fr_i), Literature PBT Pumping Number Range (Rentang literatur PBT, ~0.5-0.87), Experimental Nq Anchor Point (Cell 3) = 0.027 at ~450 RPM, Vortex-Corrected Pumping/Flow Number Nq(N), Unbaffled Stirred-Tank Vessel (vortexing regime), Unbaffled Vortex Onset Limit (Fr=0.15, ~266 RPM)

### Community 9 - "Chemistry Params & Mass-Balance Ceiling"
Cohesion: 0.38
Nodes (6): ChemistryParams, feed_rate_mol_min(), mass_balance_ceiling_pH(), Henderson-Hasselbalch thermodynamics + continuous-titrant feed model.  This modu, Moles of titrant delivered per minute at time t (0 once the pump has stopped)., The one pH value a CLOSED system (titrant feed finished, no more mass entering)

### Community 10 - "CSV I/O Utilities"
Cohesion: 0.40
Nodes (5): DataFrame, _guess_column(), load_ph_csv(), CSV ingestion helpers: flexible column detection for non-engineer users who may, Read an uploaded CSV and return a DataFrame with columns guaranteed to be     na

### Community 11 - "Compartment-Count Decision Ladder"
Cohesion: 0.40
Nodes (5): Compartment number decision ladder (1 / 2 / N>2 compartments), 5L->250mL compartment degeneration (3 compartments -> 2), Reynolds/Froude regime check gate (Re_i turbulent plateau, Fr_i vortex ceiling), scale_down_250mL.ipynb (250 mL scale-down mixing model), Torispherical (dished) bottom head geometry handling

### Community 12 - "TLS Cert Generation"
Cohesion: 0.50
Nodes (4): Path, generate_self_signed_cert(), main(), Generate a self-signed TLS certificate/key pair for serving this app over HTTPS.

### Community 13 - "Failed Fix Attempts (Clamp/Sign-Preserving)"
Cohesion: 0.50
Nodes (4): Hard clamp fix attempt (rejected: absorbing negative state), Hard clamp fix pattern (generic): fails via spurious frozen branch, Sign-preserving continuation fix attempt (rejected: runaway blow-up), Sign-preserving continuation pattern (generic): fails via runaway blow-up

### Community 14 - "Macro-Mixing Time Correlations"
Cohesion: 0.67
Nodes (3): Cooke et al. 1988 multi-impeller/aspect-ratio mixing-time correction (H^2.43), Grenville (1992) turbulent-regime blend time correlation, Macro-mixing time theta_95

### Community 15 - "Moving-Mesh CFD Pitfalls"
Cohesion: 0.67
Nodes (3): Moving-mesh transient CFD for a growing liquid volume with a real inlet, snappyHexMesh parametric geometry pitfalls (STL scaling, non-manifold, refinement), Transient scalar transport on frozen MRF flow field: numerics/setup

## Ambiguous Edges - Review These
- `ph_adjustment_env conda environment` → `myr_simulation top-level requirements.txt dependencies`  [AMBIGUOUS]
  requirements.txt · relation: conceptually_related_to

## Knowledge Gaps
- **41 isolated node(s):** `D_eff derived axial eddy-dispersion coefficient (~7.0e-5 m2/s)`, `Forward batch-validation workflow (freeze params, run forward)`, `Scale-down decision ladder: constant N vs. constant P/V vs. constant tip speed`, `Torispherical (dished) bottom head geometry handling`, `Standard k-epsilon turbulence model` (+36 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **8 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **What is the exact relationship between `ph_adjustment_env conda environment` and `myr_simulation top-level requirements.txt dependencies`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **Why does `myr_simulation project CLAUDE.md` connect `Project Notebooks & Modeling References` to `Compartment-Count Decision Ladder`, `DS-Pool Reactor Debugging & Fixes`, `DS-Pool Kinetic Fit & SOP`?**
  _High betweenness centrality (0.123) - this node is a cross-community bridge._
- **Why does `webapp/README.md - Mass-Transfer & Kinetic Fitting Tool` connect `Webapp Engine Architecture` to `Project Notebooks & Modeling References`?**
  _High betweenness centrality (0.067) - this node is a cross-community bridge._
- **Why does `parameter-estimation skill (fitting rate constants and transport parameters)` connect `Project Notebooks & Modeling References` to `Webapp Engine Architecture`?**
  _High betweenness centrality (0.049) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `FeedParams` (e.g. with `ModelSelectionResult` and `FitResult`) actually correct?**
  _`FeedParams` has 2 INFERRED edges - model-reasoned connections that need verification._
- **Are the 2 inferred relationships involving `ChemistryParams` (e.g. with `ModelSelectionResult` and `FitResult`) actually correct?**
  _`ChemistryParams` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `D_eff derived axial eddy-dispersion coefficient (~7.0e-5 m2/s)`, `Forward batch-validation workflow (freeze params, run forward)`, `Scale-down decision ladder: constant N vs. constant P/V vs. constant tip speed` to the rest of the system?**
  _41 weakly-connected nodes found - possible documentation gaps or missing edges._