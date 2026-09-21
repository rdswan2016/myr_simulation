# myr_mechanistic_model_5L_optimized_FIXED

_Source notebook: `myr_mechanistic_model_5L_optimized_FIXED.ipynb`_

# Myristoylation (DesB30 -> Detemir) 5 L Tank Mechanistic Model — Process Optimization (DS-Pool Kinetics)
Reactor-coupled version of the **DS-Pool refined kinetic fit** (`02jul2026_ds_pool_refined_results.pkl`): 6 lumped rate
constants + free reaction orders (n=0.12 on NHS-Myr, m=1.89 on the acylation substrate), pooling regiochemistry into
`desb30 -> detemir (B29) / other_DS1 (A1+B1) -> DS2_total (all di-acyl) -> DS3 (tri-acyl)` + NHS-ester hydrolysis +
dynamic glycine quench, coupled to a 3-compartment hydrodynamic model of a 5 L reactor (3086.1768 mL working volume).

**Process control strategy implemented here:**
- NHS-Myr delivered as a **continuous timed feed** (peristaltic pump, volumetric-rate controlled) instead of an instant spike.
- **Time-varying tank volume** during dosing, tracked via a moles-based state so dilution is handled exactly.
- Reaction runs to its **Detemir concentration peak**, which triggers the glycine quench (also a timed feed, on the
  **same physical pump** as the NHS-Myr addition, so both are bounded to the pump's **10-200 mL/min** range).
- Section 8 **sweeps NHS-Myr pump flow rate and stock concentration** to target a Detemir peak at **~30 minutes**
  (fast enough to avoid wasting time, slow enough to stay controllable) while keeping side-product formation low and
  the total batch volume under the **5 L** reactor ceiling.
- Glycine dose is **dynamically sized** (1.1x molar excess over residual NHS-Myr at the peak).
- The full reaction network keeps running during the quench stage since NHS-Myr hasn't hit zero.

**Numerical note:** the fitted NHS-Myr reaction order (n=0.12) makes the rate law nearly singular as NHS-Myr
concentration approaches zero — which happens constantly in this process (natural depletion, and by design during
the quench). A small (1 nM) smooth floor is applied to every fractional-power term (Section 3) so the ODE stays
well-posed; without it, stiff solvers either stall chasing the singularity or produce unphysical negative
concentrations that silently break the mass balance.

```python
# --- TOP CELL: RUN THIS ONCE ---
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import dill
import os
# -------------------------------
```

## 1. Stoichiometry & process parameters
DesB30 at **2 mg/mL**. NHS-Myr stock concentration and pump flow rate are **process knobs** — a starting value is
set here and confirmed/refined by the Section 8 optimization (which explores stock concentration and pump speed
jointly, subject to the pump's 10-200 mL/min range and the 5 L reactor volume limit).
Impeller fixed at **450 rpm** (found to have no measurable effect on yield/selectivity across the 200-450 rpm range tested).

```python
MW_desb30 = 5706.0
MW_nhs_myr = 325.44
MW_glycine = 75.07

V_total_mix_mL = 3086.1768
V_buffer_L = (V_total_mix_mL / 2.0) / 1000.0
V_acn_L    = (V_total_mix_mL / 2.0) / 1000.0

DESB30_STOCK_MG_ML = 2.0
NHS_STOCK_MG_ML    = 0.482   # starting point; re-examined jointly with pump speed in Section 8
GLY_STOCK_G_ML     = 2.0

# Shared peristaltic pump (NHS-Myr addition AND glycine quench) and reactor hard limits
PUMP_MIN_ML_MIN, PUMP_MAX_ML_MIN = 10.0, 200.0
V_REACTOR_MAX_ML = 5000.0

mass_desb30_mg = DESB30_STOCK_MG_ML * (V_buffer_L * 1000.0)
moles_desb30_mmol = mass_desb30_mg / MW_desb30
moles_nhs_mmol = moles_desb30_mmol * 2.0   # 2.0x molar excess *charge target* (pump may be stopped early -- see Sec. 8)
V_nhs_add_L = ((moles_nhs_mmol * MW_nhs_myr) / NHS_STOCK_MG_ML) / 1000.0

N_RPM = 450.0

print(f'DesB30: {moles_desb30_mmol:.5f} mmol | NHS-Myr (2.0x charge target): {moles_nhs_mmol:.5f} mmol | '
      f'V_nhs_add: {V_nhs_add_L*1000:.1f} mL ({V_nhs_add_L/V_buffer_L*100:.1f}% of pre-dose buffer volume)')
```

## 2. Kinetics (DS-Pool refined fit: 6 lumped rate constants + free n, m)
Loaded from `02jul2026_ds_pool_refined_results.pkl` (`res_final.x`: first 6 entries are `log(k)`, last 2 are the
reaction orders `n` (NHS-Myr) and `m` (acylation substrate)). Species are **pooled**, not full regiochemistry:
`other_DS1` = A1 + B1, `DS2_total` = A1_B1 + A1_B29 + B1_B29, `DS3` = A1_B29_B1 — this is what the refined fit
actually resolves; the earlier full 9-species split isn't separable from this dataset.
Glycine quench reuses **k[1]** (desb30 -> other_DS1) as an SME-specified proxy: A1 (pooled inside other_DS1) is the
desb30 N-terminal glycine, structurally analogous to free glycine's amine.

```python
SPECIES = ['ACN', 'nhs_myr', 'desb30', 'detemir', 'other_DS1', 'DS2_total', 'DS3', 'glycine', 'gly_myr']
IDX = {name: i for i, name in enumerate(SPECIES)}
N_SPECIES = len(SPECIES)
N_COMPARTMENTS = 3
REACTIVE_SPECIES = ['nhs_myr', 'desb30', 'detemir', 'other_DS1', 'DS2_total', 'DS3', 'glycine', 'gly_myr']

def load_kinetics(path='02jul2026_ds_pool_refined_results.pkl'):
    try:
        with open(path, 'rb') as f:
            bundle = dill.load(f)
        x = bundle['res_final'].x
        k = list(np.exp(x[:6]))
        n, m = float(x[6]), float(x[7])
        print(f"SUCCESS: Loaded DS-Pool refined fit (R^2={bundle.get('r2'):.4f}, RMSE={bundle.get('rmse'):.4f})")
        print(f'k = {k}')
        print(f'n (NHS-Myr order) = {n:.4f}, m (substrate order) = {m:.4f}')
    except Exception as e:
        print(f'WARNING: Using fallback parameters ({e}).')
        k = [4824.6912685716925, 1134.9009076406996, 177.1701312057297,
             13710.007107598605, 2114.9020833218583, 0.12584014547014327]
        n, m = 0.11873756023404723, 1.8915245993718806
    return k, n, m

k, order_n, order_m = load_kinetics()
k_quench = k[1]  # desb30 -> other_DS1 proxy (A1-type N-terminal-glycine acylation)
```

## 3. Dimensional bridge & reaction network

```python
EPS_M = 1e-9  # 1 nM smoothing floor on every fractional-power term. The fitted NHS-Myr order (n=0.12) makes
              # d(rate)/d(NHS) diverge as NHS -> 0 (which happens constantly here). A hard clamp to zero
              # traps stiff solvers on a spurious flat branch (dydt=0 for any negative excursion, so the
              # implicit step equation is satisfied by *any* negative value); this smooth floor keeps the
              # ODE well-posed everywhere while leaving rates unchanged above ~uM concentrations.
def smooth_floor(c):
    return 0.5*(c + np.sqrt(c*c + EPS_M*EPS_M))

def calc_rate_mM_s(k_per_min, c_nhs_mM, c_sub_mM):
    c_nhs_M = smooth_floor(c_nhs_mM/1000.0)
    c_sub_M = smooth_floor(c_sub_mM/1000.0)
    return k_per_min * (c_nhs_M**order_n) * (c_sub_M**order_m) * (1000.0/60.0)

def calc_hydro_rate_mM_s(k_per_min, c_nhs_mM):
    c_nhs_M = smooth_floor(c_nhs_mM/1000.0)
    return k_per_min * c_nhs_M * (1000.0/60.0)

def reaction_derivs_mM_s(C, k, k_quench):
    """C: dict of one compartment's concentrations (mM). Pooled 5-reaction network + hydrolysis + quench."""
    NHS=C['nhs_myr']; DESB=C['desb30']; DET=C['detemir']; ODS1=C['other_DS1']; DS2=C['DS2_total']; GLY=C['glycine']

    r1 = calc_rate_mM_s(k[0], NHS, DESB)   # desb30 -> detemir (B29)
    r2 = calc_rate_mM_s(k[1], NHS, DESB)   # desb30 -> other_DS1 (A1+B1)
    r3 = calc_rate_mM_s(k[2], NHS, DET)    # detemir -> DS2_total
    r4 = calc_rate_mM_s(k[3], NHS, ODS1)   # other_DS1 -> DS2_total
    r5 = calc_rate_mM_s(k[4], NHS, DS2)    # DS2_total -> DS3
    r_hydro  = calc_hydro_rate_mM_s(k[5], NHS)     # NHS -> dead (hydrolysis)
    r_quench = calc_rate_mM_s(k_quench, NHS, GLY)  # NHS + glycine -> gly_myr

    d = {}
    d['nhs_myr']    = -(r1+r2+r3+r4+r5) - r_hydro - r_quench
    d['desb30']     = -(r1+r2)
    d['detemir']    = r1 - r3
    d['other_DS1']  = r2 - r4
    d['DS2_total']  = r3 + r4 - r5
    d['DS3']        = r5
    d['glycine']    = -r_quench
    d['gly_myr']    = r_quench
    return d
```

## 4. Hydrodynamics (3-compartment, impeller-driven mixing)

```python
def calc_compartments(V_current_L, N_rpm=N_RPM):
    Area_tank = np.pi*(0.205/2)**2
    H_fluid = (V_current_L/1000.0)/Area_tank
    baffle_factor = 0.1
    N_rps = N_rpm/60.0
    Q_L_s = (0.5*baffle_factor*N_rps*(0.075**3))*1000.0
    h_ratio = 0.02/H_fluid
    V2_L = (np.pi*(0.075/2)**2*(0.3*0.075))*1000.0
    V3_L = max(0.01, V_current_L*h_ratio - V2_L/2)
    V1_L = max(0.01, V_current_L - V2_L - V3_L)
    Q12 = Q_L_s*(h_ratio/0.33); Q23 = Q_L_s*((1.0-h_ratio)/0.67)
    return V1_L, V2_L, V3_L, Q12, Q23
```

## 5. Moles-based transport
Switched the state representation from **concentration** to **total moles per compartment**. This is what correctly handles a *time-varying* tank volume (the NHS-Myr feed adds ~24%+ extra volume during dosing) without needing extra dilution-correction terms: moles are conserved through transport and volume changes automatically; concentration is only reconstructed (moles / current volume) when evaluating rate laws.

```python
def transport_moles(row_moles, V1, V2, V3, Q12, Q23):
    N1, N2, N3 = row_moles
    C1, C2, C3 = N1/V1, N2/V2, N3/V3
    f12 = Q12*(C1-C2); f23 = Q23*(C2-C3)
    return np.array([-f12, f12-f23, f23])
```

## 6. Stage 1: ACN/buffer mixing

```python
def run_stage_1():
    V_tot = V_buffer_L + V_acn_L
    V1, V2, V3, Q12, Q23 = calc_compartments(V_tot)
    C_des_init = moles_desb30_mmol / V_tot
    y0 = np.zeros((N_SPECIES, N_COMPARTMENTS))
    y0[IDX['ACN']] = np.array([1.0, 0.0, 0.0]) * np.array([V1, V2, V3])
    y0[IDX['desb30']] = C_des_init * np.array([V1, V2, V3])

    def ode(t, yf):
        y = yf.reshape(N_SPECIES, N_COMPARTMENTS); dydt = np.zeros_like(y)
        dydt[IDX['ACN']] = transport_moles(y[IDX['ACN']], V1, V2, V3, Q12, Q23)
        return dydt.flatten()

    sol = solve_ivp(ode, [0, 60], y0.flatten(), max_step=0.5)
    return sol.y[:, -1].reshape(N_SPECIES, N_COMPARTMENTS), V_tot
```

## 7. Stage 2: dosed NHS-Myr feed + reaction (time-varying volume)
NHS-Myr is fed at a **constant molar rate** into zone 1 over `T_feed` minutes; tank volume grows continuously during
dosing and compartment geometry is recomputed at every ODE evaluation. Integrated with `Radau` (tight tolerances,
no artificial max-step) since the fractional-order kinetics are genuinely stiff; a **terminal event** stops the
integration once total NHS-Myr drops below 0.01% of the charge (remaining conversion is then negligible — and
without this cutoff the solver stalls chasing the n=0.12 singularity as NHS-Myr asymptotically approaches zero).

```python
def make_nhs_exhausted_event(nhs_threshold_mmol):
    def event(t, yf):
        y = yf.reshape(N_SPECIES, N_COMPARTMENTS)
        return y[IDX['nhs_myr']].sum() - nhs_threshold_mmol
    event.terminal = True
    event.direction = -1
    return event

def run_stage_2_dosed(y_initial_moles, V_start, T_feed_min, moles_nhs_mmol_local, V_nhs_add_L_local,
                       t_hold_s=3600.0, N_rpm=N_RPM, nhs_exhausted_frac=1e-4):
    T_feed_s = T_feed_min * 60.0
    feed_rate_mmol_s = moles_nhs_mmol_local / T_feed_s

    def get_geom(t):
        V_added = min(t, T_feed_s) / T_feed_s * V_nhs_add_L_local
        return calc_compartments(V_start + V_added, N_rpm)

    y0 = y_initial_moles.copy()  # moles carry straight through stage boundaries

    def ode(t, yf):
        y = yf.reshape(N_SPECIES, N_COMPARTMENTS)
        V1, V2, V3, Q12, Q23 = get_geom(t)
        Vvec = np.array([V1, V2, V3])
        dydt = np.zeros_like(y)
        for j in range(N_COMPARTMENTS):
            Cj = {name: y[IDX[name], j] / Vvec[j] for name in REACTIVE_SPECIES}
            d = reaction_derivs_mM_s(Cj, k, k_quench)
            for name, val_mM_s in d.items():
                dydt[IDX[name], j] += val_mM_s * Vvec[j]
        for name in REACTIVE_SPECIES:
            dydt[IDX[name]] += transport_moles(y[IDX[name]], V1, V2, V3, Q12, Q23)
        if t < T_feed_s:
            dydt[IDX['nhs_myr'], 0] += feed_rate_mmol_s
        return dydt.flatten()

    ev = make_nhs_exhausted_event(nhs_exhausted_frac * moles_nhs_mmol_local)
    sol = solve_ivp(ode, [0, t_hold_s], y0.flatten(), method='Radau', dense_output=True,
                     rtol=1e-7, atol=1e-10, events=ev)
    return sol, V_start

def detemir_peak_time_min(T_feed_min, y_mix, V1_, moles_nhs_mmol_local, V_nhs_add_L_local,
                           t_hold_s=None, grid_n=1200):
    if t_hold_s is None:
        t_hold_s = T_feed_min * 60.0 + 2400.0   # generous hold beyond feed end; terminal event usually cuts it short
    sol, V_start = run_stage_2_dosed(y_mix, V1_, T_feed_min, moles_nhs_mmol_local, V_nhs_add_L_local, t_hold_s)
    t_end = sol.t[-1]
    t_grid = np.linspace(0, t_end, grid_n)
    y_grid = sol.sol(t_grid).reshape(N_SPECIES, N_COMPARTMENTS, -1)
    T_feed_s = T_feed_min * 60.0
    det_conc = np.zeros(grid_n)
    for i, t in enumerate(t_grid):
        V_added = min(t, T_feed_s) / T_feed_s * V_nhs_add_L_local
        V1c, _, _, _, _ = calc_compartments(V_start + V_added)
        det_conc[i] = y_grid[IDX['detemir'], 0, i] / V1c
    pk = np.argmax(det_conc)
    return t_grid[pk] / 60.0, sol, V_start, t_grid, y_grid, pk
```

## 8. Process optimization: target a 30-minute Detemir peak, minimize side products
Sweeps **NHS-Myr pump flow rate** (10-200 mL/min — the same physical peristaltic pump used later for the glycine
quench) and **NHS-Myr stock concentration** jointly. For each combination the required nominal feed duration is
`T_feed = V_nhs_add / pump_rate`; Stage 2 is run to find where Detemir peaks, and the total volume after NHS-Myr
addition is checked against the 5 L reactor ceiling (with headroom left for the glycine dose).

**What the sweep shows:** selectivity (Detemir / total conjugated product) is only weakly sensitive to the
pump-rate/stock-concentration choice once dosing is slow enough to avoid a sharp initial NHS-Myr spike — it stays
in a narrow ~73-76% band across the whole explored range. The dominant lever is simply **how fast NHS-Myr is
delivered**, which sets *when* the Detemir peak lands (too fast: peak arrives early and hard to catch in time;
too slow: wastes process time). So the optimization reduces to: pick a stock concentration, then dial the pump
speed via a precise 1-D root-find so the peak lands at exactly 30 minutes.

The **existing NHS-Myr stock (0.482 mg/mL)** is kept unchanged (no reformulation needed) since diluting it further
only buys a marginal (~1 percentage point) selectivity gain in the sweep — within the noise of the kinetic fit
itself (R^2=0.77). Only the **pump flow rate** is adjusted.

```python
# --- Joint sweep: NHS-Myr pump flow rate x stock concentration ---
y_mix, V1_ = run_stage_1()

def V_nhs_add_L_for(stock_mg_ml):
    return ((moles_nhs_mmol * MW_nhs_myr) / stock_mg_ml) / 1000.0

stock_grid = np.linspace(0.2, 1.2, 6)                          # mg/mL, within the 0.1-1.5 mg/mL range considered
pump_grid  = np.linspace(PUMP_MIN_ML_MIN, PUMP_MAX_ML_MIN, 10)  # mL/min

sweep_rows = []
for stock in stock_grid:
    V_nhs = V_nhs_add_L_for(stock)
    if V_total_mix_mL + V_nhs*1000.0 > V_REACTOR_MAX_ML - 50:   # leave headroom for the (tiny) glycine dose
        continue
    for pump in pump_grid:
        T_feed = (V_nhs*1000.0) / pump
        t_pk, _, _, _, y_g, pk_i = detemir_peak_time_min(T_feed, y_mix, V1_, moles_nhs_mmol, V_nhs, grid_n=400)
        yp = y_g[:, :, pk_i]
        det = yp[IDX['detemir']].sum()
        side = sum(yp[IDX[s]].sum() for s in ('other_DS1', 'DS2_total', 'DS3'))
        sweep_rows.append(dict(stock=stock, pump=pump, T_feed=T_feed, t_peak=t_pk,
                                detemir=det, side=side, sel=100*det/(det+side)))

band = [r for r in sweep_rows if 25.0 <= r['t_peak'] <= 35.0]
print(f'Swept {len(sweep_rows)} (stock, pump) combinations -> {len(band)} land the Detemir peak within 25-35 min.')
print(f'{"stock(mg/mL)":>13} {"pump(mL/min)":>13} {"t_peak(min)":>12} {"selectivity(%)":>15}')
for r in sorted(band, key=lambda r: -r['sel'])[:6]:
    print(f"{r['stock']:13.2f} {r['pump']:13.2f} {r['t_peak']:12.2f} {r['sel']:15.2f}")

# --- Refine at the current NHS-Myr stock (Section 1): root-find the pump rate that lands the peak at exactly 30 min ---
def t_peak_minus_target(T_feed_min, target=30.0):
    t_pk, *_ = detemir_peak_time_min(T_feed_min, y_mix, V1_, moles_nhs_mmol, V_nhs_add_L)
    return t_pk - target

T_FEED_MIN = brentq(t_peak_minus_target, 20.0, 100.0, xtol=0.05)
pump_flow_rate_mL_min = (V_nhs_add_L * 1000.0) / T_FEED_MIN
assert PUMP_MIN_ML_MIN <= pump_flow_rate_mL_min <= PUMP_MAX_ML_MIN, 'pump rate out of the 10-200 mL/min range!'
assert V_total_mix_mL + V_nhs_add_L*1000.0 <= V_REACTOR_MAX_ML, 'exceeds 5 L reactor volume!'

t_peak_min, sol2, V_start2, t_grid2, y_grid2, pk = detemir_peak_time_min(T_FEED_MIN, y_mix, V1_, moles_nhs_mmol, V_nhs_add_L)
pct_delivered = t_peak_min / T_FEED_MIN * 100

print()
print(f'CHOSEN: NHS-Myr stock = {NHS_STOCK_MG_ML} mg/mL (unchanged), pump flow rate = {pump_flow_rate_mL_min:.2f} mL/min')
print(f'Nominal programmed feed duration = {T_FEED_MIN:.2f} min; Detemir peak at t = {t_peak_min:.2f} min')
print(f'>>> Pump must be stopped early, at t = {t_peak_min:.2f} min, delivering '
      f'{pct_delivered:.1f}% of the planned NHS-Myr charge ({pct_delivered/100*2:.2f}x molar excess instead of 2.0x) <<<')

y_peak_moles = y_grid2[:, :, pk]
T_feed_s = T_FEED_MIN * 60.0
V_added_at_peak = min(t_peak_min*60.0, T_feed_s) / T_feed_s * V_nhs_add_L
V_at_peak = V_start2 + V_added_at_peak   # ACTUAL tank volume at the quench trigger moment
```

## 9. Stage 3: dynamic glycine quench
Glycine is dosed as a **timed feed on the same peristaltic pump** used for the NHS-Myr addition (10-200 mL/min),
sized to **1.1x molar excess over residual NHS-Myr present at the Detemir peak**. Because the residual NHS-Myr at
that point is small, the required glycine volume is small too, so the feed duration is short regardless of which
pump speed is used here. The full acylation/hydrolysis network keeps running alongside the quench reaction.

```python
def run_stage_3_quench(y_state_moles, V_start, gly_pump_rate_mL_min, N_rpm=N_RPM, excess_factor=1.1, t_hold_after_s=1800.0):
    """Glycine dosed as a timed feed on the SAME pump used for NHS-Myr (must stay within 10-200 mL/min)."""
    N_nhs_total = y_state_moles[IDX['nhs_myr']].sum()
    dose_gly_mmol = excess_factor * N_nhs_total
    V_gly_add_L = ((dose_gly_mmol * MW_glycine) / GLY_STOCK_G_ML) / 1000.0
    T_gly_feed_min = max(1e-3, (V_gly_add_L * 1000.0) / gly_pump_rate_mL_min)
    T_gly_feed_s = T_gly_feed_min * 60.0
    gly_feed_rate_mmol_s = dose_gly_mmol / T_gly_feed_s
    t_hold_s = T_gly_feed_s + t_hold_after_s

    def get_geom(t):
        V_added = min(t, T_gly_feed_s) / T_gly_feed_s * V_gly_add_L
        return calc_compartments(V_start + V_added, N_rpm)

    y0 = y_state_moles.copy()

    def ode(t, yf):
        y = yf.reshape(N_SPECIES, N_COMPARTMENTS)
        V1, V2, V3, Q12, Q23 = get_geom(t)
        Vvec = np.array([V1, V2, V3])
        dydt = np.zeros_like(y)
        for j in range(N_COMPARTMENTS):
            Cj = {name: y[IDX[name], j] / Vvec[j] for name in REACTIVE_SPECIES}
            d = reaction_derivs_mM_s(Cj, k, k_quench)
            for name, val in d.items():
                dydt[IDX[name], j] += val * Vvec[j]
        for name in REACTIVE_SPECIES:
            dydt[IDX[name]] += transport_moles(y[IDX[name]], V1, V2, V3, Q12, Q23)
        if t < T_gly_feed_s:
            dydt[IDX['glycine'], 0] += gly_feed_rate_mmol_s
        return dydt.flatten()

    sol = solve_ivp(ode, [0, t_hold_s], y0.flatten(), method='Radau', dense_output=True, rtol=1e-7, atol=1e-10)
    V_tot = V_start + V_gly_add_L
    return sol, V_tot, dose_gly_mmol, N_nhs_total, T_gly_feed_min

sol3, V_tot3, dose_gly_mmol, nhs_residual, T_gly_feed_min = run_stage_3_quench(
    y_peak_moles, V_at_peak, gly_pump_rate_mL_min=pump_flow_rate_mL_min)
assert V_tot3*1000.0 <= V_REACTOR_MAX_ML, 'post-quench volume exceeds 5 L reactor!'
print(f'Residual NHS-Myr at peak (mmol): {nhs_residual:.5f}')
print(f'Glycine dose (1.1x excess, mmol): {dose_gly_mmol:.5f} '
      f'({dose_gly_mmol*MW_glycine/GLY_STOCK_G_ML*1000:.2f} uL of {GLY_STOCK_G_ML} g/mL stock, '
      f'fed over {T_gly_feed_min*60:.1f} s at {pump_flow_rate_mL_min:.2f} mL/min)')
print(f'Post-quench batch volume: {V_tot3*1000:.1f} mL (limit {V_REACTOR_MAX_ML:.0f} mL)')
```

## 10. Mass balance check

```python
V1f, V2f, V3f, _, _ = calc_compartments(V_tot3)
final_state = sol3.y[:, -1].reshape(N_SPECIES, N_COMPARTMENTS)
protein_species = ['desb30', 'detemir', 'other_DS1', 'DS2_total', 'DS3']
total_protein = sum(final_state[IDX[s]].sum() for s in protein_species)
print(f'Protein mass balance: {total_protein:.6f} mmol vs initial {moles_desb30_mmol:.6f} mmol '
      f'(ratio={total_protein/moles_desb30_mmol:.6f})')
print(f'Gly-Myr (quenched byproduct) formed: {final_state[IDX["gly_myr"]].sum():.5f} mmol')
```

## 11. Visualization

```python
t2_min = t_grid2[:pk+1] / 60.0
t3_min = (sol3.t / 60.0) + t_peak_min
t_all_min = np.concatenate([t2_min, t3_min])

def bulk_stage2(name): return y_grid2[IDX[name], 0, :pk+1]
def bulk_stage3(name): return sol3.y.reshape(N_SPECIES, N_COMPARTMENTS, -1)[IDX[name], 0, :]
def bulk_full(name): return np.concatenate([bulk_stage2(name), bulk_stage3(name)])

# Convert to concentration for stage 2 (volume was time-varying there)
def bulk_conc_stage2(name):
    T_feed_s = T_FEED_MIN * 60.0
    out = np.zeros(pk+1)
    for i, t in enumerate(t_grid2[:pk+1]):
        V_added = min(t, T_feed_s) / T_feed_s * V_nhs_add_L
        V1c, _, _, _, _ = calc_compartments(V_start2 + V_added)
        out[i] = y_grid2[IDX[name], 0, i] / V1c
    return out

def bulk_conc_stage3(name):
    return sol3.y.reshape(N_SPECIES, N_COMPARTMENTS, -1)[IDX[name], 0, :] / V1f

def bulk_conc_full(name):
    return np.concatenate([bulk_conc_stage2(name), bulk_conc_stage3(name)])

plt.figure(figsize=(13, 7))
plt.plot(t_all_min, bulk_conc_full('detemir'), label='Detemir (desired, B29)', color='green', linewidth=3)
plt.plot(t_all_min, bulk_conc_full('other_DS1'), label='other_DS1 (A1+B1, mono-acyl impurity)', color='purple')
plt.plot(t_all_min, bulk_conc_full('DS2_total'), label='DS2_total (A1_B1+A1_B29+B1_B29, di-acyl pooled)', color='magenta', linestyle='--')
plt.plot(t_all_min, bulk_conc_full('DS3'), label='DS3 (A1_B29_B1, tri-acyl)', color='black', linestyle=':')
# NHS-Myr dips a few uM below zero at the raw solver nodes (verified, not a dense-output artifact) --
# the near-zero-order fit (n=0.12) makes d(nhs_myr)/dt <= 0 everywhere, so a stiff-solver overshoot near
# c=0 has no restoring force and drifts slightly negative. Clamped for display only; all mass-balance,
# yield, and purity figures elsewhere are computed from the raw state and are unaffected.
plt.plot(t_all_min, np.clip(bulk_conc_full('nhs_myr'), 0.0, None), label='NHS-Myr (linker)', color='red', linestyle='-.')
plt.plot(t_all_min, bulk_conc_full('desb30'), label='DesB30 (precursor)', color='orange', linestyle='-.')
plt.plot(t_all_min, bulk_conc_full('gly_myr'), label='Gly-Myr (quenched byproduct)', color='steelblue', linewidth=2)
plt.axvline(x=t_peak_min, color='black', linestyle='--', linewidth=2,
            label=f'Detemir peak -- pump stopped early ({pct_delivered:.0f}% dosed) / quench triggered')
plt.title('5L Tank Bioconjugation: Dosed Feed + Dynamic Quench (DS-Pool Kinetics)', fontsize=13, fontweight='bold')
plt.xlabel('Process Time (Minutes)')
plt.ylabel('Bulk (Zone 1) Concentration (mM)')
plt.legend(fontsize=9, ncol=2)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

## 12. SOP Summary -- Optimum Operating Conditions
Prints the process parameters your scientists need to write the SOP, derived directly from the optimization above.

```python
nhs_delivered_mL = V_added_at_peak * 1000.0
nhs_delivered_mg = nhs_delivered_mL * NHS_STOCK_MG_ML
nhs_delivered_mmol = nhs_delivered_mg / MW_nhs_myr
actual_molar_excess = nhs_delivered_mmol / moles_desb30_mmol

gly_dose_mL = dose_gly_mmol * MW_glycine / GLY_STOCK_G_ML

detemir_final_mmol = final_state[IDX['detemir']].sum()
impurity_species = ['other_DS1', 'DS2_total', 'DS3']
impurity_final_mmol = sum(final_state[IDX[s]].sum() for s in impurity_species)
yield_pct = detemir_final_mmol / moles_desb30_mmol * 100.0
purity_pct = detemir_final_mmol / (detemir_final_mmol + impurity_final_mmol) * 100.0
desb30_unreacted_mmol = final_state[IDX['desb30']].sum()
conversion_pct = (1.0 - desb30_unreacted_mmol/moles_desb30_mmol) * 100.0

print('='*70)
print('  5 L TANK MYRISTOYLATION -- OPTIMUM OPERATING CONDITIONS (SOP INPUTS)')
print('='*70)
print()
print('CHARGE / STOCK PREPARATION')
print(f'  DesB30 stock concentration:        {DESB30_STOCK_MG_ML:.3f} mg/mL')
print(f'  DesB30 charge (in buffer half):     {mass_desb30_mg:.2f} mg  ({moles_desb30_mmol:.5f} mmol)')
print(f'  NHS-Myr stock concentration:        {NHS_STOCK_MG_ML:.3f} mg/mL  (unchanged; confirmed by Sec. 8 sweep)')
print(f'  NHS-Myr stock volume to prepare:     {V_nhs_add_L*1000:.1f} mL  (full 2.0x excess charge, {moles_nhs_mmol:.5f} mmol)')
print(f'  Glycine quench stock concentration: {GLY_STOCK_G_ML:.3f} g/mL')
print()
print('IMPELLER')
print(f'  Speed:                               {N_RPM:.0f} rpm  (45-deg PBT; confirmed no measurable')
print(f'                                        effect on yield/purity across 200-450 rpm)')
print()
print('NHS-MYR ADDITION (peristaltic pump, 10-200 mL/min range)')
print(f'  Pump flow rate:                       {pump_flow_rate_mL_min:.2f} mL/min  (constant rate)')
print(f'  Nominal programmed feed duration:     {T_FEED_MIN:.1f} min')
print(f'  >>> ACTUAL pump STOP time:             {t_peak_min:.2f} min  (stop early - do NOT run full program)')
print(f'  NHS-Myr volume actually delivered:    {nhs_delivered_mL:.1f} mL  ({pct_delivered:.1f}% of prepared stock)')
print(f'  Actual molar excess delivered:        {actual_molar_excess:.2f}x  (vs 2.0x if fully dosed)')
print()
print('GLYCINE QUENCH (same pump, switched to glycine reservoir)')
print(f'  Trigger time:                         {t_peak_min:.2f} min (immediately after pump stop, at Detemir peak)')
print(f'  Pump flow rate:                       {pump_flow_rate_mL_min:.2f} mL/min  (same setting as NHS-Myr addition)')
print(f'  Dose:                                 {gly_dose_mL*1000:.2f} uL  ({dose_gly_mmol:.5f} mmol, '
      f'1.1x molar excess over residual NHS-Myr), fed over {T_gly_feed_min*60:.1f} s')
print(f'  Residual NHS-Myr at trigger:          {nhs_residual:.5f} mmol')
print()
print('PREDICTED OUTCOME (at end of quench hold)')
print(f'  Detemir yield:                         {yield_pct:.1f}% of charged DesB30')
print(f'  Product purity (Detemir / total conj.): {purity_pct:.1f}%')
print(f'  DesB30 conversion:                      {conversion_pct:.1f}%')
print(f'  Total batch volume (post-quench):       {V_tot3:.3f} L  (limit 5.000 L)')
print('='*70)
```
