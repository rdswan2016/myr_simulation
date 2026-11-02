# Graph Report - Python/myr_simulation  (2026-11-02)

## Corpus Check
- 0 files · ~0 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 205 nodes · 314 edges · 24 communities (14 shown, 10 thin omitted)
- Extraction: 83% EXTRACTED · 17% INFERRED · 0% AMBIGUOUS · INFERRED: 52 edges (avg confidence: 0.82)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- Webapp App & Chemistry Engine
- Skills, Da Screening & Project Refs
- Kinetic Data Prep & 5L Reactor Notebooks
- Forward Validation & Compartment Ladder
- Webapp Transport & Validation Docs
- Mixing Timescales & CFD Zoning
- Numerical Debugging & Smooth Floor
- Nq vs RPM Vortex Chart
- Impeller Froude & Pumping Number
- CSV I/O Utilities
- TLS Cert Generation
- Failed Clamp/Sign-Preserving Fixes
- Macro-Mixing Time Correlations
- graphify MCP Registration
- Moving-Mesh & snappyHexMesh Pitfalls
- Micro-Mixing Timescales
- MRF vs Sliding Mesh
- Kinetika 3-Zone Fed-Batch
- Power & Pumping Numbers
- Webapp Engine Init
- Identifiability Diagnostics
- LES Turbulence Model
- RNG k-epsilon Model
- Standard k-epsilon Model

## God Nodes (most connected - your core abstractions)
1. `myr_simulation project CLAUDE.md` - 25 edges
2. `run_model_selection()` - 14 edges
3. `FeedParams` - 13 edges
4. `ChemistryParams` - 11 edges
5. `webapp/README.md - Mass-Transfer & Kinetic Fitting Tool` - 11 edges
6. `simulate_ph()` - 10 edges
7. `DS-Pool Refined Kinetic Fit (02jul2026_ds_pool_refined_results.pkl)` - 10 edges
8. `FitResult` - 8 edges
9. `fit_n_zone_model()` - 8 edges
10. `DS-Pool Reactor Model Status` - 8 edges

## Surprising Connections (you probably didn't know these)
- `Limitations and data needs (only 2 batches, 4/10 kinetic params loose, no pilot/3500 L validation, need 5-6 batches)` --conceptually_related_to--> `DS-Pool Refined Kinetic Fit (02jul2026_ds_pool_refined_results.pkl)`  [INFERRED]
  Python/myr_simulation/ringkasan_eksekutif_myristoylasi.pdf → Python/myr_simulation/CLAUDE.md
- `250 mL RANS fit vs experiment pH curve (Exp3; 250-mL-only C2=2.2337 vs combined C2=2.2469, RMSE 0.042)` --semantically_similar_to--> `scale_down_250mL.ipynb (250 mL scale-down mixing model)`  [INFERRED] [semantically similar]
  Python/myr_simulation/figs_ringkasan/rans_5L_combined_fit_comparison.png → Python/myr_simulation/CLAUDE.md
- `1 nM smooth floor on fractional-power terms (n=0.12 near-singular at NHS-Myr -> 0)` --semantically_similar_to--> `Hard rules (smooth-floor, mass-balance validation, fixed SOP, diagnostic-first, Da-feed screening)`  [INFERRED] [semantically similar]
  Python/myr_simulation/_notebook_md/myr_mechanistic_model_5L_optimized_FIXED.md → Python/myr_simulation/CLAUDE.md
- `Smooth non-negative floor pattern (generic working fix)` --semantically_similar_to--> `Smooth non-negative floor regularization (accepted fix)`  [INFERRED] [semantically similar]
  Python/myr_simulation/.claude/skills/ode-stiff-kinetics/SKILL.md → Python/myr_simulation/.claude/memory/project_dspool_reactor_model.md
- `9-species regiochemistry + hydrolysis kinetics network` --semantically_similar_to--> `Species pooling convention (other_DS1, DS2_total, DS3)`  [INFERRED] [semantically similar]
  Python/myr_simulation/myr_mechanistic_model_5L_optimized.html → Python/myr_simulation/CLAUDE.md

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Portable str_cfd/str skill set for STR mixing and kinetics modeling** — claude_skills_biotech_modeling_skill_biotech_modeling, claude_skills_cfd_mixing_fundamentals_skill_cfd_mixing_fundamentals, claude_skills_mixing_time_correlations_skill_mixing_time_correlations, claude_skills_ode_stiff_kinetics_skill_ode_stiff_kinetics, claude_skills_parameter_estimation_skill_parameter_estimation [EXTRACTED 1.00]
- **DS-Pool reactor model numerical debugging fix sequence** — hard_clamp_fix_attempt, sign_preserving_continuation_fix_attempt, smooth_floor_fix, terminal_ode_event, claude_memory_feedback_numerical_debugging_principle [EXTRACTED 1.00]
- **Macro/meso/micro mixing timescale hierarchy for Da screening** — macro_mixing_theta95, meso_mixing, micro_mixing_kolmogorov, da_m_screening, da_feed_number [EXTRACTED 1.00]
- **Vortex-suppressed Nq falling below literature PBT range, calibrated to single experimental anchor** — python_myr_simulation_grafik_nq_vs_rpm_nq_vortex_corrected, python_myr_simulation_grafik_nq_vs_rpm_literature_pbt_range, python_myr_simulation_grafik_nq_vs_rpm_nq_experimental_cell3, python_myr_simulation_grafik_nq_vs_rpm_unbaffled_vessel [INFERRED 0.85]

## Communities (24 total, 10 thin omitted)

### Community 0 - "Webapp App & Chemistry Engine"
Cohesion: 0.11
Nodes (36): Standalone web application: automated multi-compartment mass-transfer parameter, ChemistryParams, feed_rate_mol_min(), FeedParams, initial_split(), mass_balance_ceiling_pH(), pH_from_moles(), ndarray (+28 more)

### Community 1 - "Skills, Da Screening & Project Refs"
Cohesion: 0.13
Nodes (26): 5L_tank_transfer_mass_3_compartments.ipynb (3/4-compartment fit, 160 RPM), biotech-modeling skill (bioprocess modeling judgment), cfd-mixing-fundamentals skill (CFD-based STR mixing analysis), mixing-time-correlations skill (reference formulas), ode-stiff-kinetics skill (diagnosing stiff/singular ODE-PDE kinetics), parameter-estimation skill (fitting rate constants and transport parameters), Damkohler number at feed point (Da_feed = t_micro/t_reaction_local), Da_M mixing Damkohler screening (tau_M/tau_R) (+18 more)

### Community 2 - "Kinetic Data Prep & 5L Reactor Notebooks"
Cohesion: 0.10
Nodes (24): dataprep_rx_kinetic.ipynb (experimental data prep), 30jul2026_hydrolysis_kinetic_model.dill (superseded 9-species fit), Reactor model: 9-species regiochemistry + hydrolysis + glycine quench, 3-compartment (bulk/impeller/headspace), 3 stages (ACN mix, spike reaction, quench), 9-species regiochemistry + hydrolysis fit with n=m=1 hardcoded, log-param global DE + L-BFGS-B; saved as dill bundle, dataprep_rx_kinetic.md (notebook export: HPLC data prep), DS-pool refined kinetics coupled to 3-compartment 5 L reactor (n=0.12 NHS-Myr, m=1.89), moles-based state, Dosed NHS-Myr feed with time-varying volume; solve T_feed for 30-min Detemir peak (bound 1-30 min); dynamic glycine quench, FIXED fit: n,m freed (0.1-4), smooth quadratic mass-balance penalty replaces 1e12 cliff, boundary-pinning check (+16 more)

### Community 3 - "Forward Validation & Compartment Ladder"
Cohesion: 0.14
Nodes (17): Forward batch-validation workflow (freeze params, run forward), Compartment number decision ladder (1 / 2 / N>2 compartments), 5L->250mL compartment degeneration (3 compartments -> 2), Stage 3: dynamic glycine quench, gmp-3500l project (source of profitability figure), pra5_myr_model.ipynb (forward-predictive validation vs. Run 1/Run 2), Profitability scenarios chart (5 scenarios; annual revenue/OPEX/net profit; OPEX/batch Rp 3,941M vs BEP Rp 994M), 250 mL RANS fit vs experiment pH curve (Exp3; 250-mL-only C2=2.2337 vs combined C2=2.2469, RMSE 0.042) (+9 more)

### Community 4 - "Webapp Transport & Validation Docs"
Cohesion: 0.16
Nodes (15): Henderson-Hasselbalch <-> transport ODE coupling (moles-based, single closure on probe zone), Mass-balance ceiling check (closed-system equilibrium vs. data peak), N-compartment ODE engine (moles-based state, Fickian junction exchange), Nested extra-sum-of-squares F-test for compartment-count selection, webapp/README.md - Mass-Transfer & Kinetic Fitting Tool, webapp/requirements.txt dependencies (streamlit, plotly, cryptography, etc.), Structural cross-check: conserved-quantity ceiling (independent of optimizer), webapp/app.py (Streamlit UI) (+7 more)

### Community 5 - "Mixing Timescales & CFD Zoning"
Cohesion: 0.20
Nodes (14): DS_Pool_Kinetic_Model.ipynb (superseded, non-refined), DS_Pool_Kinetic_Model_Refined.ipynb (produces active kinetic fit), DS-Pool Refined Kinetic Fit (02jul2026_ds_pool_refined_results.pkl), Fit versioning / staleness discipline, Log-parameterize + dilled-bundle-with-fallback handoff pattern, AIC/BIC comparison: first-order baseline vs free n,m vs refinement, DS_Pool_Kinetic_Model.md (notebook export), DS_Pool_Kinetic_Model_Refined.md (notebook export) (+6 more)

### Community 6 - "Numerical Debugging & Smooth Floor"
Cohesion: 0.23
Nodes (12): Numerical Debugging Feedback: diagnostic-first before patching, Project Memory Index (MEMORY.md), DS-Pool Reactor Model Status, D_eff derived axial eddy-dispersion coefficient (~7.0e-5 m2/s), Diagnostic-first, not patch-first workflow, myr_mechanistic_model_5L_optimized_FIXED.ipynb (current 2-compartment reactor model), Residual negative NHS-Myr display artifact (display-only clip fix), Smooth non-negative floor regularization (accepted fix) (+4 more)

### Community 7 - "Nq vs RPM Vortex Chart"
Cohesion: 0.25
Nodes (8): Bakker & Fasano 1993 feed-point-vs-bulk zoning CFD study, Bourne & Hilber 1990 critical addition time relation (tau_crit*N^n=const), CFD-to-ODE coupling interface (fixed Q_ij / zone-averaged epsilon, no tight coupling), Corrsin full mixing time tau_M, Meso-mixing (feed-plume spreading, tau_D/tau_S timescales), Middleton et al. 1986 competitive-consecutive reaction system (A+B->R, R+B->S), RTD-based experimental validation (E(t), N-tanks-in-series), CFD zone/compartment boundary identification (v_z=0 locus, Q_ij flux integration)

### Community 8 - "Impeller Froude & Pumping Number"
Cohesion: 0.48
Nodes (7): Nq vs RPM: Vortex-Froude Theory Chart (Cell 3 anchored), Impeller Froude Number (Fr_i), Literature PBT Pumping Number Range (Rentang literatur PBT, ~0.5-0.87), Experimental Nq Anchor Point (Cell 3) = 0.027 at ~450 RPM, Vortex-Corrected Pumping/Flow Number Nq(N), Unbaffled Stirred-Tank Vessel (vortexing regime), Unbaffled Vortex Onset Limit (Fr=0.15, ~266 RPM)

### Community 9 - "CSV I/O Utilities"
Cohesion: 0.40
Nodes (5): DataFrame, _guess_column(), load_ph_csv(), CSV ingestion helpers: flexible column detection for non-engineer users who may, Read an uploaded CSV and return a DataFrame with columns guaranteed to be     na

### Community 10 - "TLS Cert Generation"
Cohesion: 0.50
Nodes (4): Path, generate_self_signed_cert(), main(), Generate a self-signed TLS certificate/key pair for serving this app over HTTPS.

### Community 11 - "Failed Clamp/Sign-Preserving Fixes"
Cohesion: 0.50
Nodes (4): Hard clamp fix attempt (rejected: absorbing negative state), Hard clamp fix pattern (generic): fails via spurious frozen branch, Sign-preserving continuation fix attempt (rejected: runaway blow-up), Sign-preserving continuation pattern (generic): fails via runaway blow-up

### Community 12 - "Macro-Mixing Time Correlations"
Cohesion: 0.67
Nodes (3): Cooke et al. 1988 multi-impeller/aspect-ratio mixing-time correction (H^2.43), Grenville (1992) turbulent-regime blend time correlation, Macro-mixing time theta_95

### Community 14 - "Moving-Mesh & snappyHexMesh Pitfalls"
Cohesion: 0.67
Nodes (3): Moving-mesh transient CFD for a growing liquid volume with a real inlet, snappyHexMesh parametric geometry pitfalls (STL scaling, non-manifold, refinement), Transient scalar transport on frozen MRF flow field: numerics/setup

## Ambiguous Edges - Review These
- `ph_adjustment_env conda environment` → `myr_simulation top-level requirements.txt dependencies`  [AMBIGUOUS]
  Python/myr_simulation/requirements.txt · relation: conceptually_related_to

## Knowledge Gaps
- **45 isolated node(s):** `graphify`, `D_eff derived axial eddy-dispersion coefficient (~7.0e-5 m2/s)`, `Forward batch-validation workflow (freeze params, run forward)`, `Scale-down decision ladder: constant N vs. constant P/V vs. constant tip speed`, `Torispherical (dished) bottom head geometry handling` (+40 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **10 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **What is the exact relationship between `ph_adjustment_env conda environment` and `myr_simulation top-level requirements.txt dependencies`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **Why does `myr_simulation project CLAUDE.md` connect `Skills, Da Screening & Project Refs` to `Kinetic Data Prep & 5L Reactor Notebooks`, `Forward Validation & Compartment Ladder`, `Mixing Timescales & CFD Zoning`, `Numerical Debugging & Smooth Floor`?**
  _High betweenness centrality (0.177) - this node is a cross-community bridge._
- **Why does `webapp/README.md - Mass-Transfer & Kinetic Fitting Tool` connect `Webapp Transport & Validation Docs` to `Skills, Da Screening & Project Refs`?**
  _High betweenness centrality (0.067) - this node is a cross-community bridge._
- **Why does `parameter-estimation skill (fitting rate constants and transport parameters)` connect `Skills, Da Screening & Project Refs` to `Webapp Transport & Validation Docs`?**
  _High betweenness centrality (0.048) - this node is a cross-community bridge._
- **Are the 2 inferred relationships involving `FeedParams` (e.g. with `ModelSelectionResult` and `FitResult`) actually correct?**
  _`FeedParams` has 2 INFERRED edges - model-reasoned connections that need verification._
- **Are the 2 inferred relationships involving `ChemistryParams` (e.g. with `ModelSelectionResult` and `FitResult`) actually correct?**
  _`ChemistryParams` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `graphify`, `D_eff derived axial eddy-dispersion coefficient (~7.0e-5 m2/s)`, `Forward batch-validation workflow (freeze params, run forward)` to the rest of the system?**
  _45 weakly-connected nodes found - possible documentation gaps or missing edges._