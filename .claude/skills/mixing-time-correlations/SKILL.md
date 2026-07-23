---
name: mixing-time-correlations
description: Reference correlations and formulas for STR mixing — macro-mixing time (theta_95, Grenville correlations), power number Np and pumping number Nq (turbulent-plateau values and geometry corrections), continuity/flux-integration basis for Nq/Q_ij, and Kolmogorov/Batchelor/Corrsin/engulfment micro- and meso-mixing timescale formulas, each with explicit validity ranges and when to prefer CFD instead. Use when computing or cross-checking a mixing-time, Nq/Np, or Damkohler-screening number.
---

# Mixing-time, power/pumping-number, and micro-mixing reference correlations

Reference data only — computable equations, coefficients, and validity ranges. For the
procedural judgment on when to trust these vs. CFD, see `cfd-mixing-fundamentals.md` §3/§6.
Source: Paul, Atiemo-Obeng, Kresta (eds.), *Handbook of Industrial Mixing* (Wiley, 2003) —
Ch.5 "Computational Fluid Mixing" §5-2.1.1, Ch.6 "Mechanically Stirred Vessels" §6-2.3,
Ch.9 "Blending of Miscible Liquids" §9-2/§9-5, Ch.13 "Mixing and Chemical Reactions" §13-1/§13-2.
Items not page-verified against the source are marked `[TO VERIFY IN HANDBOOK]`.

## 0. Governing equation: continuity (mass conservation)

Every Nq/Q/Q_ij value in §2 and in `cfd-mixing-fundamentals.md` §3–4 that is described as
"integrate the flow through a surface" rests on nothing more than the continuity equation —
this is the equation that makes flux integration valid, not an independent assumption layered
on top of it. As stated in the source (Ch.5 §5-2.1.1, PDF p. 314):
```
General (compressible) form:      ∂ρ/∂t + ∂(ρUᵢ)/∂xᵢ = 0                    (eq. 5-5)
Constant-density form used for    ∂Uᵢ/∂xᵢ = 0                                (eq. 5-4)
  liquid STR mixing:
```
Practical consequence: because density is constant for the liquid systems this file covers,
the net volumetric flow crossing *any* closed surface is exactly zero — that is precisely why
integrating the outflow velocity across one side of a control surface (the impeller discharge
disk for Nq, a zone-interface plane for Q_ij) gives a physically exact flow rate, not an
approximation. It is also why a **net** vs. **gross** flux distinction matters whenever a
surface has flow crossing both directions (see `cfd-mixing-fundamentals.md` §4): continuity
guarantees the *net* flux through a closed surface is zero, it does not say each open
sub-surface only carries flow one way.

The momentum equation (Navier–Stokes, Ch.5 §5-2.1.2, eq. 5-6) is the parent equation RANS
(§2 below, and `cfd-mixing-fundamentals.md` §2) is derived from by Reynolds-averaging — it is
not restated separately here because no formula in this file or in `cfd-mixing-fundamentals.md`
requires the instantaneous (laminar) form directly; every quantity used for STR compartment
work is already post-averaging (RANS) or an empirical correlation. Eq. 5-6 itself is at
Ch.5 §5-2.1.2, PDF p. 316, if ever needed verbatim.

## 1. Macro-mixing time correlations

**Correlation A — Grenville (1992) turbulent-regime blend time** (Ch.9 eq. 9-1, 9-2):
```
Po^(1/3) · N · θ_95 · D² / (T^1.5 · H^0.5) = 5.20     (± 10.0% std. dev.)
```
Reduces, for H = T, to `Po^(1/3)·N·θ_95·(D/T)² = 5.20`. Rearranged to Reynolds/Fourier form:
`Po^(1/3)·Re = 5.20/Fo`. **Validity**: turbulent regime, **baffled** vessels only — measured on
0.30–2.97 m diameter vessels with standard torispherical base, standard baffles, single impeller
at H/3 above the base, D/T = 1/3 to 1/2 (hydrofoils, pitched/flat-blade turbines, disk turbines
all tested). Transitional-regime companion correlation (eq. 9-6/9-7):
```
Po^(1/3) · Re = 183 / √Fo     (± 17.4%)
⟺  N·θ_95 = 1832·(T/D)² / (Po^(2/3)·Re)
```
Turbulent/transitional boundary: `Po^(1/3)·Re_TT = 6370` (equivalently `1/Fo_TT = 1225`).
Confidence bands: constant = `5.20 + 0.52s` (turbulent) or `183 + 31.1s` (transitional), where
s=1/2/3 gives 67%/95%/99% confidence.

**Correlation B — Grenville blend time, Re-binned restatement** (Ch.13 eq. 13-8, 13-9; same
underlying data set as Correlation A, restated with slightly different rounding — use as a
cross-check, cite Correlation A as primary):
```
N·τ_B = 5.4·(T/D)² / Po^(1/3)                          for Re > 6400
N·τ_B = (1/Re)·184.2·(T/D)² / Po^(2/3)                 for 500 < Re < 6400
```
Typical values quoted: vessel blend time ≈ 2 s in a 1 L vessel, ≈ 20 s in a 20 000 L vessel
(low-viscosity liquids).

**Multiple-impeller / aspect-ratio > 1 correction** (Cooke et al. 1988, Ch.9 §9-2.6):
```
θ_m · Po^(1/3) · N / D^3.3 ∝ H^2.43
```
Mixing time increases with liquid height at exponent **2.43**, materially steeper than the
single-impeller correlation's `H^0.5` — do not extrapolate the single-impeller correlation to a
tall/multi-impeller vessel by just plugging in a larger H. Staged/zoned mixing between multiple
radial-flow impellers roughly doubles blend time vs. a single impeller; axial (hydrofoil)
impeller pairs at aspect ratio 2 cut this zoning penalty by 30–60% vs. radial pairs.

**Conversion between homogeneity thresholds — exact, not approximate** (Ch.9 eq. 9-18):
first-order decay of concentration fluctuations gives
```
θ_z = θ_95 · ln[(100 − %homogeneity)/100] / ln(0.05)
θ_99 = 1.537 · θ_95
```

**Jet-mixing macro-mixing time** (for a pumped jet rather than an impeller, Ch.9 §9-5, Grenville
& Tilton 1996 model, validated against Fossett & Prosser 1949 form which fit best among
literature correlations):
```
θ = K_Z · Z² / (U·D)                    K_Z = 3.00 (± 11.0%)
```
`Z` = jet free-path length, `U` = nozzle velocity, `D` = nozzle diameter. **Validity**: Re>10 000
(turbulent), 0.2<H/T<2.0, 0.178<V<1200 m³, 1.32×10⁻²<(UD/Z)<0.137 m/s, 86<Z/D<753. Derived from
the Corrsin time-scale form `θ = K_Z·(Z/η)^x·(Z²/ε_Z)`; regression confirmed the Corrsin-predicted
exponent x=−1/3.

## 2. Power number Np and pumping number Nq

**Definitions and base formulas** (Ch.5 eq. 5-47/5-48/5-49; Ch.6 eq. 6-5):
```
Np = P / (ρ·N³·D⁵)                     Nq = Q / (N·D³)
P  = Np·ρ·N³·D⁵ / g_c                  Q  = Nq·N·D³
v_tip = π·N·D                          P/V = Np·ρ·N³·D⁵ / V
```
`P` from CFD: blade-surface torque integration, `P=2πNτ` — not a volume integral of the
turbulence model's dissipation field (model-dependent, unreliable for this purpose). `Q` from
CFD: integrate outflow velocity across a discharge surface (disk for axial, cylindrical-shell
section for radial impellers).

**Turbulent-plateau Nq by impeller type** (Ch.6 Table 6-3, baffled vessels):

| Impeller | Nq |
|---|---|
| Propeller | 0.4–0.6 |
| Pitched blade turbine (PBT) | 0.79 |
| Hydrofoil | 0.55–0.73 |
| Retreat curve blade | 0.3 |
| Flat-blade turbine | 0.7 |
| Rushton (disk flat-blade turbine) | 0.72 |
| Smith (hollow-blade turbine) | 0.76 |

**Np vs. Re regime behavior** (Ch.6 §6-2.3.2): `Np ∝ Re⁻¹` for Re<10 (laminar); Np ≈ constant
for Re>10⁴ (turbulent plateau); weak Re-dependence for 100<Re<10⁴ (transitional). Both Nq and Np
curves (e.g. for a PBT) only reach their turbulent plateau above **Re≈10 000** and should not be
trusted below **Re_i=1000** at all.

**Geometry-dependence corrections to Np** (turbulent regime, Ch.6 eqs. 6-10 to 6-16):
```
Blade width (6-blade Rushton):      Np ∝ (W/D)^1.45
Blade width (4-blade 45° PBT):      Np ∝ (W/D)^0.65
Blade count, 3–6 blades:            Np ∝ (n/D)^0.8
Blade count, 6–12 blades:           Np ∝ (n/D)^0.7
PBT blade angle θ:                  Np ∝ (sinθ)^2.6
PBT off-bottom clearance C:         Np ∝ (C/D)^-0.25
Propeller pitch p (1.0<p/D<2.0, Re>1000):  Np ∝ (p/D)^1.5
```
Baffle number×width (`N_b·B`) increases Np up to the standard 4-baffle, B=T/10 plateau, then Np
is constant thereafter (level depends on D/T). Multiple impellers on one shaft: combined Np is
**not simply additive**; for axial (PBT) pairs combined power is significantly less than 2×
single-impeller power; for flat-blade turbines it can exceed 2× depending on spacing S/D
(typical spacing ≈ 1 impeller diameter; closer spacing increases interaction).

**No Froude-number correction for unbaffled vessels was found with a numeric constant in the
extracted Ch.5/6/9 passages.** `[TO VERIFY IN HANDBOOK — the `Fr_crit≈0.15` value currently used
in `str/parameter-estimation.md` was not located verbatim in this Handbook edition's Ch.5, §6-2/
§6-3, or §9-2/§9-5; Chapter 2 "Turbulence in Mixing Applications" and Ch.6 §6-3 "Flow
Characteristics" were not fully extracted and may contain it.]` Do not cite this file as the
source for that constant until verified.

## 3. Kolmogorov scale and micro-mixing time estimates

All formulas as printed in the source (Ch.13 §13-2.1.3, eqs. 13-10 to 13-16); `ν`=kinematic
viscosity, `ε`=local turbulent energy dissipation rate per unit mass, `D_AB`=molecular
diffusivity, `Sc=ν/D_AB`:

```
Kolmogorov length:        η = (ν³/ε)^(1/4)
Kolmogorov time:          τ_K = (ν/ε)^(1/2)
Batchelor length:         λ_B = η/√Sc = (ν·D_AB²/ε)^(1/4)
Corrsin mixing time       τ_M = (L_s²/ε)^(1/3) + 0.5·(ν/ε)^(1/2)·ln(Sc)      [Sc ≫ 1, liquids]
  (Baldyga-Bourne restate the 2nd term as asinh(0.05·Sc)·(ν/ε)^(1/2) — more rigorous, same order)
Engulfment micromixing:   E = 0.06·(ε/ν)^(1/2)   ⟹   τ_E = 1/E = 17·(ν/ε)^(1/2)
                          (equivalently τ_E = 17·(µ/ρε)^(1/2))
```
`τ_E` scales as `√viscosity` at fixed ε/ρ — a 100× viscosity increase roughly doubles it. Valid
for **Sc < ~4000** (turbulent engulfment regime); above that, mixing proceeds by viscous
stretching instead (not detailed here). For Sc≈1000 (typical aqueous solute), the Batchelor
length can be **~30× smaller** than the Kolmogorov length — Batchelor, not Kolmogorov, is the
relevant smallest-striation scale for most aqueous mixing problems.

**Mesomixing time estimates** (Baldyga & Bourne; Ch.13 eqs. 13-18, 13-19) — the scale *between*
macro and micro, relevant whenever feed rate exceeds local mixing rate:
```
τ_D = Q_B / (U·D_t)                D_t = 0.1·k²/ε           (dispersion form)
τ_S = A·(L_s²/ε)^(1/3)             A ≈ 1–2                  (Corrsin/disintegration form)
```
`Q_B`=fed reagent's volumetric feed rate, `U`=local surrounding-fluid velocity at the feed
point, `L_s`=concentration macroscale set by feed geometry (e.g. feed pipe diameter). Both
forms are explicitly noted in the source as limited (neither captures a deliberately
high-momentum feed jet correctly).

**`ε ≈ P/(ρV)` as a named "volume-averaged dissipation rate" formula was not found verbatim** in
the extracted Ch.5/6/9/13 passages — it is dimensionally consistent with §2's P formula and
standard practice, but `[TO VERIFY IN HANDBOOK, likely Ch.2 "Turbulence in Mixing Applications,"
not extracted this pass]` before citing this file as its source.

## 4. Validity notes and when to prefer CFD over correlations

- **Turbulent regime required** for all §1/§2 correlations as stated: Re_i > ~10 000 for the
  Np/Nq turbulent plateau; Re > 6400 for the primary θ_95 blend-time correlation (transitional
  companion correlations exist for 500–6400, see §1).
- **Baffled vessel required** for every θ_95/Nq/Np correlation in §1/§2 — none of the tables or
  equations above were validated on an unbaffled geometry, and no unbaffled equivalent is given
  in this Handbook's extracted sections. **Do not apply §1/§2 numbers to an unbaffled vessel
  without an explicit, documented correction** — and note that the specific `Fr_crit` correction
  constant is currently unverified against this source (see §2's caveat).
- **H/T ≠ 1 invalidates** the simplified single-impeller θ_95 form (`Po^(1/3)Nθ_95(D/T)²=5.20`)
  — use the full form with the explicit H^0.5 dependence, and for H/T far from 1 with multiple
  impellers, prefer the Cooke et al. H^2.43 correction over the single-impeller correlation
  entirely.
- **Prefer CFD (see `cfd-mixing-fundamentals.md` §3–4) over these correlations** when: the vessel
  is unbaffled or partially baffled; C/D, D/T, or H/T fall outside the ranges validated above;
  multiple impellers are used at nonstandard spacing; or the process needs a spatially-resolved
  ε(x) field rather than a single volume-averaged number (e.g. a genuine feed-point Da_feed
  screening at the actual local ε near the feed line, rather than the vessel-average ε≈P/(ρV)).
