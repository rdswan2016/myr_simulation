# DS_Pool_Kinetic_Model

_Source notebook: `DS_Pool_Kinetic_Model.ipynb`_

# DS-Pool Kinetic Model

Collapses the 8-species site-resolved model down to 5 species by lumping unconfirmed
peaks into degree-of-substitution (DS) pools, while keeping the two **confirmed**
peaks (`desb30`, `detemir`) tracked individually.

**Why:** only `desb30` (starting material) and `detemir` (the commercial API) were
identified against a reference standard. The other 6 peaks were assigned to
mono/di/tri tiers by retention-time clustering + stoichiometric peak-count matching,
not by mass spec — so *within* a tier, which specific peak is A1 vs B1 vs A1_B29
etc. is not something we can trust from kinetics alone.

**Species tracked:**
- `nhs_myr` — myristoylation reagent
- `desb30` — starting material (DS0) — **confirmed**
- `detemir` — LysB29-myristoyl (DS1) — **confirmed**
- `other_DS1` — the other DS1 peak(s) lumped (A1 + B1) — pool, unconfirmed identity
- `DS2_total` — all three DS2 (di-substituted) peaks lumped — pool, unconfirmed identity
- `DS3` — the single DS3 (tri-substituted) peak — tier itself is unambiguous, only one peak possible

**Why lumping by molarity is valid here:** every species within a DS tier shares the
same MW (adding a myristoyl group changes mass by the same +210.36 Da regardless of
*which* site it's added to):
- DS1 members (A1, B1, detemir): all MW = 5916.885
- DS2 members (A1_B1, A1_B29, B1_B29): all MW = 6127.246
- DS3 member (A1_B29_B1): MW = 6337.607

So pool molarity = (sum of g/L across that tier's peaks) / (shared tier MW).

**Reaction network at the pool level** (5 myristoylation edges + hydrolysis, vs. 12
in the full site-resolved lattice — a much better match to 368 raw data points /
274 usable points in this pooled representation):

```
desb30    --k1-->  detemir      (trustworthy: detemir's identity is confirmed)
desb30    --k2-->  other_DS1    (lumped: desb30 -> either A1 or B1)
detemir   --k3-->  DS2_total    (detemir picks up a 2nd myristoyl, any site)
other_DS1 --k4-->  DS2_total    (A1 or B1 picks up a 2nd myristoyl, any site)
DS2_total --k5-->  DS3          (any DS2 species picks up the 3rd/last myristoyl)
nhs_myr   --k6-->  (hydrolysis, unchanged from the site-resolved model)
```

8 total free parameters (k1..k6, n, m) vs. 10 for the original site-resolved model
and 15 for the full 12-edge lattice — this is the version to trust as a headline
result, since it makes no claim about regiochemistry beyond the two confirmed peaks.

Run this **after** your other fits have finished (site-resolved notebook,
`pathway_search.py`) — don't run multiple DE optimizations at once on the same
4-core machine, they'll split cores and all get slower.

## 1. Imports

```python
import time, functools
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, differential_evolution
from sklearn.metrics import mean_squared_error, r2_score

print = functools.partial(print, flush=True)
```

## 2. Data loading + DS-pool aggregation

Loads `kinetic_data_processed.csv` (must be in the same folder as this notebook)
and aggregates the 8 raw peaks into the 5 DS-pool species described above.

```python
MW = {
    'desb30': 5706.524,
    'DS1': 5916.885,   # shared by A1, B1, detemir
    'DS2': 6127.246,   # shared by A1_B1, A1_B29, B1_B29
    'DS3': 6337.607,   # A1_B29_B1
    'nhs_myr': 325.449, 'nhs': 115.088
}

DS1_PEAKS = ['A1', 'B1', 'detemir']     # confirmed: detemir. unconfirmed: A1, B1
DS2_PEAKS = ['A1_B1', 'A1_B29', 'B1_B29']
DS3_PEAKS = ['A1_B29_B1']

df = pd.read_csv('kinetic_data_processed.csv')
df['Time'] = df['Time'].astype(float)
df['Perbandingan'] = df['Perbandingan'].astype(float)

ratios = sorted(df['Perbandingan'].unique())
time_points = sorted(df['Time'].unique())
time_exp = np.sort(df['Time'].unique())

def get_gL_series(df_ratio, peak_name):
    sub = df_ratio[df_ratio['Peak Identification'] == peak_name].sort_values('Time')
    return sub.set_index('Time')['Conc'].reindex(time_points, fill_value=0.0).values

experimental_dataset = {}
for r in ratios:
    df_ratio = df[df['Perbandingan'] == r]
    desb30_gL = get_gL_series(df_ratio, 'desb30')
    detemir_gL = get_gL_series(df_ratio, 'detemir')
    other_ds1_gL = get_gL_series(df_ratio, 'A1') + get_gL_series(df_ratio, 'B1')
    ds2_gL = sum(get_gL_series(df_ratio, p) for p in DS2_PEAKS)
    ds3_gL = get_gL_series(df_ratio, 'A1_B29_B1')

    experimental_dataset[r] = {
        'desb30': desb30_gL / MW['desb30'],
        'detemir': detemir_gL / MW['DS1'],
        'other_DS1': other_ds1_gL / MW['DS1'],
        'DS2_total': ds2_gL / MW['DS2'],
        'DS3': ds3_gL / MW['DS3'],
    }

pool_species_list = ['desb30', 'detemir', 'other_DS1', 'DS2_total', 'DS3']

N_DATAPOINTS = 0
for r in ratios:
    for sp in pool_species_list:
        vals_gL = experimental_dataset[r][sp] * MW['desb30' if sp == 'desb30' else
                                                    ('DS1' if sp in ('detemir', 'other_DS1') else
                                                     ('DS2' if sp == 'DS2_total' else 'DS3'))]
        N_DATAPOINTS += int(np.sum(~np.isnan(vals_gL) & (vals_gL > 1e-6)))
print(f"Total usable data points (DS-pool representation): {N_DATAPOINTS}")
print(f"Ratios: {ratios}")
```

## 3. The model

```python
class DSPoolKineticModel:
    species = ['nhs_myr', 'desb30', 'detemir', 'other_DS1', 'DS2_total', 'DS3']
    idx = {name: i for i, name in enumerate(species)}

    def __init__(self, params):
        self.k = params[:6]   # k1..k6 (k6 = hydrolysis)
        self.n = params[6]
        self.m = params[7]

    def _ode_system(self, t, y):
        y_safe = np.maximum(y, 1e-12)
        idx = self.idx
        C_nhsm = y_safe[idx['nhs_myr']]
        C_desb = y_safe[idx['desb30']]
        C_det = y_safe[idx['detemir']]
        C_oDS1 = y_safe[idx['other_DS1']]
        C_DS2 = y_safe[idx['DS2_total']]

        r1 = self.k[0] * (C_nhsm ** self.n) * (C_desb ** self.m)   # desb30 -> detemir
        r2 = self.k[1] * (C_nhsm ** self.n) * (C_desb ** self.m)   # desb30 -> other_DS1
        r3 = self.k[2] * (C_nhsm ** self.n) * (C_det ** self.m)    # detemir -> DS2
        r4 = self.k[3] * (C_nhsm ** self.n) * (C_oDS1 ** self.m)   # other_DS1 -> DS2
        r5 = self.k[4] * (C_nhsm ** self.n) * (C_DS2 ** self.m)    # DS2 -> DS3
        r_hydro = self.k[5] * C_nhsm

        dydt = np.zeros_like(y)
        dydt[idx['nhs_myr']] = -(r1 + r2 + r3 + r4 + r5) - r_hydro
        dydt[idx['desb30']] = -(r1 + r2)
        dydt[idx['detemir']] = r1 - r3
        dydt[idx['other_DS1']] = r2 - r4
        dydt[idx['DS2_total']] = r3 + r4 - r5
        dydt[idx['DS3']] = r5
        return dydt

    def simulate(self, y0, t_span, t_eval):
        return solve_ivp(self._ode_system, t_span, y0, t_eval=np.sort(t_eval),
                          method='Radau', atol=1e-9, rtol=1e-7)
```

## 4. Objective function

Same soft mass-balance penalty and weighted-SSE approach as the site-resolved
notebook, adapted to the 5 pool species. `detemir` keeps the highest weight since
it's both confirmed and the actual product of interest.

```python
weights = {
    'desb30': 2.0,
    'detemir': 3.0,       # this is the confirmed, business-relevant target
    'other_DS1': 0.75,
    'DS2_total': 0.75,
    'DS3': 0.5,
}
MASS_BALANCE_PENALTY_WEIGHT = 1e6

pool_mw = {'desb30': MW['desb30'], 'detemir': MW['DS1'], 'other_DS1': MW['DS1'],
           'DS2_total': MW['DS2'], 'DS3': MW['DS3']}

def objective_function_pool(params):
    k_vals = np.exp(params[:6])
    n, m = params[6], params[7]

    if np.any(k_vals > 1e8):
        return 1e12

    model = DSPoolKineticModel(list(k_vals) + [n, m])
    total_sse = 0.0

    for r in ratios:
        desb30_t0 = experimental_dataset[r]['desb30'][0]
        y0 = np.zeros(len(model.species))
        y0[model.idx['nhs_myr']] = desb30_t0 * r
        y0[model.idx['desb30']] = desb30_t0

        sol = model.simulate(y0, (0, 120), time_exp)
        if not sol.success or np.any(np.isnan(sol.y)) or np.any(np.isinf(sol.y)):
            return 1e12

        moles_in = desb30_t0
        moles_out = np.sum([sol.y[model.idx[sp]][-1] for sp in model.species if sp != 'nhs_myr'])
        rel_drift = abs(moles_in - moles_out) / moles_in
        total_sse += MASS_BALANCE_PENALTY_WEIGHT * rel_drift ** 2

        for sp in pool_species_list:
            sim_vals = sol.y[model.idx[sp]] * pool_mw[sp]
            exp_vals = experimental_dataset[r][sp] * pool_mw[sp]
            mask = ~np.isnan(exp_vals) & (exp_vals > 1e-6)
            if np.any(mask):
                err = sim_vals[mask] - exp_vals[mask]
                total_sse += np.sum((err * weights[sp]) ** 2)

    return total_sse

def objective_function_pool_first_order(log_k):
    return objective_function_pool(np.concatenate([log_k, [1.0, 1.0]]))
```

## 5. Fit driver

Same staged strategy as the site-resolved notebook: Stage 0 fits k's with n=m=1
fixed as a baseline, then Stage A frees n and m, warm-started from Stage 0's k's.
Reports AIC/BIC for both, plus a boundary-pinning check.

```python
def fit(popsize_stage0=30, maxiter_stage0=100, popsize_stageA=30, maxiter_stageA=120, workers=4, seed=42):
    bounds_k_only = [(np.log(1e-6), np.log(1e7))] * 6

    print("\nStage 0: DS-pool baseline fit (n = m = 1)...")
    t0 = time.time()
    res_de0 = differential_evolution(
        objective_function_pool_first_order, bounds_k_only,
        popsize=popsize_stage0, maxiter=maxiter_stage0, mutation=(0.5, 1.5),
        recombination=0.7, workers=workers, updating='deferred', seed=seed
    )
    res_first_order = minimize(
        objective_function_pool_first_order, res_de0.x, method='L-BFGS-B',
        bounds=bounds_k_only, options={'ftol': 1e-9, 'maxiter': 1000}
    )
    print(f"Baseline (n=m=1) SSE: {res_first_order.fun:.4f}  [{time.time()-t0:.1f}s]")

    n_m_bounds = [(0.1, 4.0), (0.1, 4.0)]
    bounds_full = bounds_k_only + n_m_bounds

    rng = np.random.default_rng(seed)
    n_dims = len(bounds_full)
    pop_size = popsize_stageA * n_dims
    init_center = np.concatenate([res_first_order.x, [1.0, 1.0]])
    init_pop = np.tile(init_center, (pop_size, 1))
    lo = np.array([b[0] for b in bounds_full]); hi = np.array([b[1] for b in bounds_full])
    noise = rng.normal(scale=0.5, size=init_pop.shape)
    noise[:, -2:] *= 0.6
    init_pop[1:] += noise[1:]
    init_pop = np.clip(init_pop, lo, hi)

    print("\nStage A: DS-pool fit with n, m free...")
    t1 = time.time()
    res_de = differential_evolution(
        objective_function_pool, bounds_full,
        popsize=popsize_stageA, maxiter=maxiter_stageA, mutation=(0.5, 1.5),
        recombination=0.7, workers=workers, updating='deferred', seed=seed, init=init_pop
    )
    res_final = minimize(
        objective_function_pool, res_de.x, method='L-BFGS-B',
        bounds=bounds_full, options={'ftol': 1e-9, 'maxiter': 1000}
    )
    print(f"Refined SSE (n,m free): {res_final.fun:.4f}   [{time.time()-t1:.1f}s]")
    print(f"Fitted n = {res_final.x[6]:.4f}, m = {res_final.x[7]:.4f}")

    for label, res, k_params in [('baseline (n=m=1)', res_first_order, 6),
                                  ('free n,m', res_final, 8)]:
        sse = res.fun
        aic = 2 * k_params + N_DATAPOINTS * np.log(sse / N_DATAPOINTS)
        bic = k_params * np.log(N_DATAPOINTS) + N_DATAPOINTS * np.log(sse / N_DATAPOINTS)
        print(f"  [{label}] SSE={sse:.4f} AIC={aic:.2f} BIC={bic:.2f}")

    print("\n--- Boundary check (free n,m fit) ---")
    names = ['k1(desb30->detemir)', 'k2(desb30->other_DS1)', 'k3(detemir->DS2)',
             'k4(other_DS1->DS2)', 'k5(DS2->DS3)', 'k6(hydro)', 'n', 'm']
    for name, val, (blo, bhi) in zip(names, res_final.x, bounds_full):
        span = bhi - blo
        if span > 0 and (val - blo) / span < 0.02:
            print(f"  WARNING: {name} pinned near LOWER bound")
        elif span > 0 and (bhi - val) / span < 0.02:
            print(f"  WARNING: {name} pinned near UPPER bound")

    k_vals = np.exp(res_final.x[:6])
    model = DSPoolKineticModel(list(k_vals) + [res_final.x[6], res_final.x[7]])
    y_true, y_pred = [], []
    for r in ratios:
        desb30_t0 = experimental_dataset[r]['desb30'][0]
        y0 = np.zeros(len(model.species))
        y0[model.idx['nhs_myr']] = desb30_t0 * r
        y0[model.idx['desb30']] = desb30_t0
        sol = model.simulate(y0, (0, 120), time_exp)
        for sp in pool_species_list:
            sim_vals = sol.y[model.idx[sp]] * pool_mw[sp]
            exp_vals = experimental_dataset[r][sp] * pool_mw[sp]
            mask = ~np.isnan(exp_vals)
            y_true.extend(exp_vals[mask]); y_pred.extend(sim_vals[mask])
    r2 = r2_score(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    print(f"\nGlobal R2 (DS-pool, free n,m)   : {r2:.4f}")
    print(f"Global RMSE (DS-pool, free n,m) : {rmse:.4f} g/L")

    print("\n================== DS-POOL RESULTS ==================")
    for name, val in zip(names[:6], k_vals):
        print(f"  {name:24}: {val:.4e}")
    print(f"  {'n':24}: {res_final.x[6]:.4f}")
    print(f"  {'m':24}: {res_final.x[7]:.4f}")
    print("=======================================================")

    return {'res_final': res_final, 'res_first_order': res_first_order,
            'r2': r2, 'rmse': rmse, 'k_vals': k_vals}
```

## 6. Run the fit

Sized for a 4-core machine (~8-12 min estimate given only 8 total parameters).
Adjust `popsize`/`maxiter` down for a quicker/rougher pass, or up for a more
thorough search if you have time/cores to spare.

```python
results = fit(popsize_stage0=30, maxiter_stage0=100,
              popsize_stageA=30, maxiter_stageA=120, workers=4)
```

## 7. Save results

```python
import pickle
with open('02jul2026_ds_pool_results.pkl', 'wb') as f:
    pickle.dump(results, f)
print("Saved 02jul2026_ds_pool_results.pkl")
```

## 8. Plot: simulated vs experimental (per ratio)

Visual check of fit quality across all 5 pool species and all 3 ratios.

```python
import matplotlib.pyplot as plt

plot_colors = {
    'desb30': '#1f77b4', 'detemir': '#2ca02c',
    'other_DS1': '#ff7f0e', 'DS2_total': '#d62728', 'DS3': '#9467bd'
}

k_vals = results['k_vals']
res_final = results['res_final']
model = DSPoolKineticModel(list(k_vals) + [res_final.x[6], res_final.x[7]])
time_smooth = np.linspace(0, 120, 1000)

fig, axes = plt.subplots(1, len(ratios), figsize=(7 * len(ratios), 6), sharey=True)
if len(ratios) == 1:
    axes = [axes]

for i, r in enumerate(ratios):
    desb30_t0 = experimental_dataset[r]['desb30'][0]
    y0 = np.zeros(len(model.species))
    y0[model.idx['nhs_myr']] = desb30_t0 * r
    y0[model.idx['desb30']] = desb30_t0

    sol_smooth = model.simulate(y0, (0, 120), time_smooth)
    sol_exp = model.simulate(y0, (0, 120), time_exp)

    for sp in pool_species_list:
        sim_gL_smooth = sol_smooth.y[model.idx[sp]] * pool_mw[sp]
        exp_gL = experimental_dataset[r][sp] * pool_mw[sp]
        axes[i].plot(time_smooth, sim_gL_smooth, color=plot_colors[sp], linewidth=1.8, alpha=0.9)
        axes[i].scatter(time_exp, exp_gL, color=plot_colors[sp], edgecolors='white',
                         linewidth=1, s=55, zorder=5)

    axes[i].set_title(f"Ratio 1:{r}", fontsize=12, fontweight='bold')
    axes[i].set_xlabel("Time (min)")
    axes[i].grid(True, linestyle='--', alpha=0.4)

axes[0].set_ylabel("Concentration (g/L)")
dummy_handles = [plt.Line2D([0], [0], color=plot_colors[sp], lw=2) for sp in pool_species_list]
fig.legend(dummy_handles, pool_species_list, loc='upper center', bbox_to_anchor=(0.5, -0.02),
           ncol=5, fontsize=11, frameon=True)
plt.suptitle(f"DS-Pool Model Fit  |  R2={results['r2']:.4f}  RMSE={results['rmse']:.4f} g/L",
             fontsize=13, fontweight='bold', y=1.03)
plt.tight_layout()
plt.show()
```
