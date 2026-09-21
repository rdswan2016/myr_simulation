# HydrolysisKineticModel

_Source notebook: `HydrolysisKineticModel.ipynb`_

```python
import sys
!{sys.executable} -m pip install scikit-learn
```

```python
# --- TOP CELL: RUN THIS ONCE ---
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy as np # Ensure numpy is imported too
from scipy.integrate import solve_ivp
import os
from sklearn.metrics import mean_squared_error, r2_score
# --- 2. SUPERCHARGED OPTIMIZER ---
from scipy.optimize import minimize, differential_evolution
import multiprocessing
# -------------------------------
```

```python
import pandas as pd
import numpy as np
import os

print("Cell 1: Memuat dataset hasil pemrosesan...")

MW = {
    'desb30': 5706.524, 'detemir': 5916.885, 'A1': 5916.885, 
    'A1_B1': 6127.246, 'A1_B29': 6127.246, 'B1_B29': 6127.246, 
    'A1_B29_B1': 6337.607, 'B1': 5916.885, 'nhs_myr': 325.449, 'nhs': 115.088
}

# Pastikan daftar ini sesuai dengan nilai unik di kolom 'Peak Identification' file kita
species_list = ['desb30', 'A1', 'detemir', 'A1_B29', 'B1_B29', 'A1_B29_B1', 'A1_B1', 'B1']

csv_filename = 'kinetic_data_processed.csv'
if not os.path.exists(csv_filename):
    raise FileNotFoundError(f"File '{csv_filename}' tidak ditemukan!")

df = pd.read_csv(csv_filename)

# Pastikan tipe data benar
df['Time'] = df['Time'].astype(float)
df['Perbandingan'] = df['Perbandingan'].astype(float)

experimental_dataset = {}
ratios = sorted(df['Perbandingan'].unique())
# Dapatkan daftar waktu unik yang tersedia dalam data (diurutkan)
time_points = sorted(df['Time'].unique())

for r in ratios:
    df_ratio = df[df['Perbandingan'] == r]
    experimental_dataset[r] = {}
    
    for sp in species_list:
        # Ambil data untuk spesies ini, urutkan berdasarkan Time
        sub_df = df_ratio[df_ratio['Peak Identification'] == sp].sort_values('Time')
        
        # Ekstrak konsentrasi (g/L) dan bagi dengan MW untuk mendapatkan Molaritas
        # Kita gunakan reindex untuk memastikan semua time_points ada (mengisi NaN dengan 0.0)
        conc_series = sub_df.set_index('Time').reindex(time_points, fill_value=0.0)['Conc']
        
        experimental_dataset[r][sp] = np.array(conc_series) / MW.get(sp, 1.0)

print("Cell 1 Selesai: Dataset telah dikonversi ke Molaritas.")
# Opsional: Cek hasil untuk satu perbandingan
print(f"Time points: {time_points}")
print(f"Contoh data (Ratio 1.0, desb30): {experimental_dataset[1.0]['desb30'][:5]}")
```

```python
import numpy as np
from scipy.integrate import solve_ivp

class FinalKineticModel:
    def __init__(self, params):
        # CHANGED: Now expects 8 k-values (index 0 to 7)
        self.k = params[:8] 
        # CHANGED: Indices for n and m shifted by 1
        self.n = params[8]
        self.m = params[9]
        
        self.species = ['nhs_myr', 'desb30', 'A1', 'detemir', 'B1', 'A1_B29', 'B1_B29', 'A1_B1', 'A1_B29_B1']
        self.idx = {name: i for i, name in enumerate(self.species)}

    def _ode_system(self, t, y):
        # 1. HARD CLIP: Force all concentrations to be at least 1e-12
        # This prevents the solver from ever trying to calculate negative numbers
        y_safe = np.maximum(y, 1e-12)
        idx = self.idx
        
        C_nhsm = y_safe[idx['nhs_myr']]
        C_desb = y_safe[idx['desb30']]
        C_a1   = y_safe[idx['A1']]
        C_det  = y_safe[idx['detemir']]
        C_a1b29= y_safe[idx['A1_B29']]
        
        r1 = self.k[0] * (C_nhsm**self.n) * (C_desb**self.m)
        r2 = self.k[1] * (C_nhsm**self.n) * (C_desb**self.m)
        r3 = self.k[2] * (C_nhsm**self.n) * (C_desb**self.m)
        
        r4 = self.k[3] * (C_nhsm**self.n) * (C_a1**self.m)
        r5 = self.k[4] * (C_nhsm**self.n) * (C_a1**self.m)
        
        r6 = self.k[5] * (C_nhsm**self.n) * (C_det**self.m)
        r7 = self.k[6] * (C_nhsm**self.n) * (C_a1b29**self.m)
        
        # NEW: Hydrolysis Rate 
        r_hydro = self.k[7] * C_nhsm 
        
        dydt = np.zeros_like(y)
        
        # CHANGED: nhs_myr is consumed by hydrolysis too
        dydt[idx['nhs_myr']]  = -(r1 + r2 + r3 + r4 + r5 + r6 + r7) - r_hydro
        dydt[idx['desb30']]   = -(r1 + r2 + r3)
        dydt[idx['A1']]       = r1 - (r4 + r5)
        dydt[idx['B1']]       = r2
        dydt[idx['detemir']]  = r3 - r6
        dydt[idx['A1_B1']]    = r4
        dydt[idx['A1_B29']]   = r5 - r7
        dydt[idx['B1_B29']]   = r6
        dydt[idx['A1_B29_B1']]= r7
        
        return dydt

    def simulate(self, y0, t_span, t_eval):
        # CHANGED: method='Radau' and tighter tolerances for the fast reaction
        return solve_ivp(self._ode_system, t_span, y0, t_eval=np.sort(t_eval), 
                         method='Radau', atol=1e-9, rtol=1e-7)

# Ensure time_exp is defined in the current cell
time_exp = np.sort(df['Time'].unique())
```

```python
def get_initial_conditions(ratio, experimental_dataset):
    # 1. Get desb30 molarity at t=0
    desb30_t0 = experimental_dataset[ratio]['desb30'][0]
    
    # 2. Calculate nhs_myr as requested
    nhs_myr_t0 = desb30_t0 * ratio
    
    # 3. Build the y0 array based on the class order
    # Note: Products are 0.0 at t=0
    y0 = np.zeros(9)
    y0[idx['nhs_myr']] = nhs_myr_t0
    y0[idx['desb30']]  = desb30_t0
    # Everything else (A1, detemir, etc.) stays 0.0
    
    return y0
```

```python
# --- 1. BALANCED ABSOLUTE WEIGHTS ---
weights = {
    'desb30': 2.0, 
    'A1': 1.0, 
    'detemir': 3.0,   # Strong priority, but not overwhelming
    'B1': 0.5, 'A1_B1': 0.5, 'A1_B29': 0.5, 'B1_B29': 0.5, 'A1_B29_B1': 0.5
}

# --- 2. THE OBJECTIVE FUNCTION ---
def objective_function_global(log_params):
    k_vals = np.exp(log_params[:8])
    n, m = 1.0, 1.0 # Locked to First-Order to prevent optimizer cheating
    
    # Reject ridiculously high k-values immediately
    if np.any(k_vals > 1e5): return 1e12
    
    full_params = list(k_vals) + [n, m]
    total_sse_weighted = 0
    model = FinalKineticModel(full_params)
    
    for r in ratios:
        desb30_t0_molL = experimental_dataset[r]['desb30'][0]
        y0 = np.zeros(len(model.species))
        y0[model.idx['nhs_myr']] = desb30_t0_molL * r
        y0[model.idx['desb30']] = desb30_t0_molL
            
        sol = model.simulate(y0, (0, 120), time_exp)
        
        # ODE Crash check
        if not sol.success or np.any(np.isnan(sol.y)) or np.any(np.isinf(sol.y)):
            return 1e12 
            
        # THE BACKBONE MASS BALANCE CHECK
        # Sum of all moles at the end (excluding nhs_myr) must equal starting desb30 moles
        moles_in = desb30_t0_molL
        moles_out = np.sum([sol.y[model.idx[sp]][-1] for sp in model.species if sp != 'nhs_myr'])
        if abs(moles_in - moles_out) / moles_in > 0.01:
            return 1e12 # Reject if backbone drifts by more than 1%

        # ERROR CALCULATION
        for sp in model.species:
            if sp == 'nhs_myr': continue
            
            sim_vals = sol.y[model.idx[sp]] * MW[sp]
            exp_vals = experimental_dataset[r][sp] * MW[sp]
            
            mask = ~np.isnan(exp_vals) & (exp_vals > 1e-6)
            if np.any(mask):
                abs_error = (sim_vals[mask] - exp_vals[mask])
                total_sse_weighted += np.sum((abs_error * weights[sp])**2)
            
    return total_sse_weighted

# Widen the bounds to prevent hitting invisible walls
bounds = [(np.log(1e-6), np.log(1e6))] * 8 

print("Stage 1: Global optimization (Differential Evolution)...")
print("This may take a few minutes. Searching 8-dimensional parameter space...")

# We are doubling the population, tripling iterations, and using all CPU cores
res_de = differential_evolution(
    objective_function_global, 
    bounds, 
    popsize=40,         # More search agents
    maxiter=150,        # Give them more time to search
    mutation=(0.5, 1.5),# Allow larger "jumps" to escape local minimums
    recombination=0.7, 
    workers=-1,         # Use all available CPU cores for parallel processing
    seed=42
)

print("Stage 2: Local refinement (L-BFGS-B Polishing)...")
res_final = minimize(
    objective_function_global, 
    res_de.x, 
    method='L-BFGS-B', 
    bounds=bounds,
    options={'ftol': 1e-9, 'maxiter': 1000} # Tighter tolerances for the final polish
)

final_params = list(res_final.x) + [1.0, 1.0] # Append the locked First-Order constants
print(f"Refined SSE: {res_final.fun:.4f}")
```

```python
opt_k = np.exp(final_params[:8])
opt_orders = final_params[8:]
fp = list(opt_k) + list(opt_orders)

print("\n================== HASIL KINETIKA MEKANISTIK GLOBAL ==================")
pathway_names = [
    'k1 (DesB30 -> A1)      ', 'k2 (DesB30 -> B1)      ', 'k3 (DesB30 -> Detemir) ',
    'k4 (A1 -> A1_B1)       ', 'k5 (A1 -> A1_B29)      ', 'k6 (Detemir -> B1_B29) ',
    'k7 (A1_B29 -> Tri)     ', 'k_hydro (NHS_myr hydrolisis)'
]
for name, val in zip(pathway_names, opt_k):
    print(f"  {name}: {val:.4e}")
print(f"  Order 'n' (NHS-myr)         : {opt_orders[0]:.3f}")
print(f"  Order 'm' (Substrat Insulin): {opt_orders[1]:.3f}")
print(f"  Final Minimum SSE (g/L)     : {res_final.fun:.5f}")
print("======================================================================\n")

# 1. Define the model with the final optimized parameters
model = FinalKineticModel(fp)

print("=== ANALISIS MASS CONVERSION & YIELD (MENIT KE-120) ===")
for r in ratios:
    # Use the same logic we used in the objective function
    desb30_t0_molL = experimental_dataset[r]['desb30'][0]
    
    y0 = np.zeros(len(model.species))
    y0[model.idx['nhs_myr']] = desb30_t0_molL * r
    y0[model.idx['desb30']] = desb30_t0_molL
    # Other products are 0.0 by default
        
    model = FinalKineticModel(final_params)
    sol = model.simulate(y0, (0, 120), time_exp)
    
    # Use model.idx to safely grab the final concentrations
    desb30_0_mass = desb30_t0_molL * MW['desb30']
    desb30_akhir_mass = sol.y[model.idx['desb30']][-1] * MW['desb30']
    detemir_akhir_mass = sol.y[model.idx['detemir']][-1] * MW['detemir']
    
    konversi = ((desb30_0_mass - desb30_akhir_mass) / desb30_0_mass) * 100
    yield_detemir = (detemir_akhir_mass / desb30_0_mass) * 100
    
    print(f"Rasio Skenario : {r}")
    print(f"  -> Konversi Massa DesB30 : {konversi:.2f}%")
    print(f"  -> Yield Hasil Detemir    : {yield_detemir:.2f}%")
```

```python
# 1. Transform the first 8 parameters from log-space to linear physical space
linear_k_vals = np.exp(res_final.x[:8])

# 2. Add the First-Order reaction orders (n=1.0, m=1.0)
linear_params = list(linear_k_vals) + [1.0, 1.0]

# 3. INSTANTIATE THE MODEL WITH THE CORRECT LINEAR PARAMETERS
model = FinalKineticModel(linear_params)

# --- Print Summary ---
print("\n================== HASIL KINETIKA MEKANISTIK GLOBAL ==================")
pathway_names = [
    'k1 (DesB30 -> A1)', 'k2 (DesB30 -> B1)', 'k3 (DesB30 -> Detemir)',
    'k4 (A1 -> A1_B1)', 'k5 (A1 -> A1_B29)', 'k6 (Detemir -> B1_B29)',
    'k7 (A1_B29 -> Tri)', 'k_hydro (NHS_myr hydrolisis)'
]

for name, val in zip(pathway_names, linear_k_vals):
    print(f"  {name:28}: {val:.4e}")
print(f"  Order 'n' (NHS-myr)         : 1.000")
print(f"  Order 'm' (Substrat Insulin): 1.000")
print(f"  Final Minimum SSE (g/L)     : {res_final.fun:.5f}")
print("======================================================================\n")

print("=== ANALISIS MASS CONVERSION & YIELD (MENIT KE-120) ===")

# --- Safe Simulation Loop ---
for r in ratios:
    desb30_t0_molL = experimental_dataset[r]['desb30'][0]
    y0 = np.zeros(len(model.species))
    y0[model.idx['nhs_myr']] = desb30_t0_molL * r
    y0[model.idx['desb30']] = desb30_t0_molL
    
    # Run the simulation safely
    sol_exp = model.simulate(y0, (0, 120), time_exp)
    
    # The "Kill Switch": If simulation fails, skip this loop so plotting doesn't crash
    if not sol_exp.success or np.any(np.isnan(sol_exp.y)):
        print(f"Simulation failed for Ratio {r}")
        continue 
        
    # Safely calculate and print mass conversion & yield
    konversi_desb = (desb30_t0_molL - sol_exp.y[model.idx['desb30']][-1]) / desb30_t0_molL * 100
    yield_detemir = (sol_exp.y[model.idx['detemir']][-1]) / desb30_t0_molL * 100
    print(f"Rasio Skenario : {r}")
    print(f"  -> Konversi Massa DesB30 : {konversi_desb:.2f}%")
    print(f"  -> Yield Hasil Detemir    : {yield_detemir:.2f}%")

    # -------------------------------------------------------------------------
    # DO NOT DELETE ANYTHING BELOW THIS LINE. 
    # Leave your existing code here that tracks RMSE, r_true, r_pred, and creates the matplotlib plots.
    # -------------------------------------------------------------------------
```

```python
import dill

# 1. Ensure your model is using the linear, optimized parameters
# (Important: Use np.exp for k-values, not the raw log-space res_final.x)
final_linear_params = list(np.exp(res_final.x[:8])) + [1.0, 1.0]
final_model = FinalKineticModel(final_linear_params)

# 2. Create a "Bundle" dictionary
# This saves the model, the specific params, and the SSE score for reference
model_bundle = {
    'model_instance': final_model,
    'parameters': final_linear_params,
    'final_sse': res_final.fun,
    'description': 'Hydrolysis-included mechanistic model',
    'date': '2026-06-30'
}

# 3. Save to disk
filename = '30jul2026_hydrolysis_kinetic_model.dill'
with open(filename, 'wb') as f:
    dill.dump(model_bundle, f)

print(f"Success! Model bundle saved to '{filename}'")
```

```python
import dill

# Load the file
with open('30jul2026_hydrolysis_kinetic_model.dill', 'rb') as f:
    loaded_bundle = dill.load(f)

# Extract your objects
my_model = loaded_bundle['model_instance']
my_params = loaded_bundle['parameters']

print(f"Loaded model version with SSE: {loaded_bundle['final_sse']:.4f}")

# Now you can use 'my_model.simulate(...)' immediately without re-optimizing!
```

```python
from sklearn.metrics import mean_squared_error, r2_score # Added r2_score

y_true_global, y_pred_global = [], []

plot_colors = {
    'desb30': '#1f77b4', 'A1': '#ff7f0e', 'detemir': '#2ca02c', 
    'A1_B29': '#d62728', 'B1_B29': '#9467bd', 'A1_B29_B1': '#8c564b', 
    'A1_B1': '#e377c2', 'B1': '#7f7f7f'
}

time_smooth = np.linspace(0, 120, 1000)
fig, axes = plt.subplots(1, 3, figsize=(22, 6.5), sharey=True)

print("=== ANALISIS ERROR PER SKENARIO RASIO ===")

# Ensure final_params exists (it was defined in your previous step)
# Then instantiate a fresh model object
model = FinalKineticModel(final_params)

for idx, r in enumerate(ratios):
    # DYNAMIC Y0: Matches your objective function logic
    desb30_t0_molL = experimental_dataset[r]['desb30'][0]
    y0 = np.zeros(len(model.species))
    y0[model.idx['nhs_myr']] = desb30_t0_molL * r
    y0[model.idx['desb30']] = desb30_t0_molL
        
    model = FinalKineticModel(final_params)
    sol_smooth = model.simulate(y0, (0, 120), time_smooth)
    sol_exp = model.simulate(y0, (0, 120), time_exp)
    
    # Use a list to store data for this ratio for local RMSE
    r_true, r_pred = [], []
    
    for sp in species_list:
        # DYNAMIC INDEXING: safer and cleaner
        ode_pos = model.idx[sp]
        
        sim_gL_smooth = sol_smooth.y[ode_pos] * MW[sp]
        sim_gL_exp = sol_exp.y[ode_pos] * MW[sp]
        exp_gL = experimental_dataset[r][sp] * MW[sp]
        
        # Masking NaNs to prevent calculation failure
        mask = ~np.isnan(exp_gL)
        r_true.extend(exp_gL[mask])
        r_pred.extend(sim_gL_exp[mask])
        
        # Add to global trackers
        y_true_global.extend(exp_gL[mask])
        y_pred_global.extend(sim_gL_exp[mask])
        
        axes[idx].plot(time_smooth, sim_gL_smooth, color=plot_colors[sp], linewidth=1.6, alpha=0.95)
        axes[idx].scatter(time_exp, exp_gL, color=plot_colors[sp], edgecolors='white', linewidth=1.1, s=65, zorder=5)

    rmse_ratio = np.sqrt(mean_squared_error(r_true, r_pred))
    print(f"Rasio 1 : {r:.1f} -> RMSE = {rmse_ratio:.4f} g/L")
    
    ax = axes[idx]
    ax.set_title(f"Rasio 1 : {r:.1f}\n(RMSE: {rmse_ratio:.3f} g/L)", fontsize=12, fontweight='bold', pad=10)
    ax.set_xlabel("Waktu Reaksi (Menit)", fontsize=12)
    if idx == 0:
        ax.set_ylabel("Konsentrasi Komponen (g/L)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)

# Hitung Statistik Global Akhir
global_r2 = r2_score(y_true_global, y_pred_global)
global_rmse = np.sqrt(mean_squared_error(y_true_global, y_pred_global))

# Penataan Legend & Title
handles, labels = axes[0].get_legend_handles_labels()
# Mengambil label unik tiruan untuk memunculkan legenda kromatogram tanpa tulisan "(Model)" berulang
dummy_handles = [plt.Line2D([0], [0], color=plot_colors[sp], lw=2) for sp in species_list]
fig.legend(dummy_handles, species_list, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=4, fontsize=11, frameon=True)

plt.suptitle(f"Kurva Kinetika Global DE (Super Smooth)\nGlobal R² = {global_r2:.4f} ({global_r2*100:.2f}%)  |  Global RMSE = {global_rmse:.4f} g/L", 
             fontsize=14, fontweight='bold', y=1.06)

plt.tight_layout()
plt.show()

print("\n==================== RINGKASAN PERFORMA MODEL ====================")
print(f" >> Akurasi Tren (Global R²) : {global_r2:.4f} ({global_r2*100:.2f}%)")
print(f" >> Deviasi Rata-rata (Global RMSE) : {global_rmse:.4f} g/L")
print("==================================================================")
```

```python
# Run this AFTER you define your model with final_params
# Let's see the rates at t=0 (the start of the reaction)
y0 = np.zeros(len(model.species))
y0[model.idx['desb30']] = experimental_dataset[1.0]['desb30'][0]
y0[model.idx['nhs_myr']] = y0[model.idx['desb30']] * 1.0

# Manually calculate the initial reaction rates (r1 to r7)
# This uses the same logic as your _ode_system
idx = model.idx
k = final_params[:8]
n, m = final_params[8], final_params[9]

# Use concentrations from y0
C_nhsm = y0[idx['nhs_myr']]
C_desb = y0[idx['desb30']]

r1 = k[0] * (C_nhsm**n) * (C_desb**m)
r2 = k[1] * (C_nhsm**n) * (C_desb**m)
r3 = k[2] * (C_nhsm**n) * (C_desb**m) # Detemir Pathway

print(f"Flux to A1 (r1): {r1:.4e}")
print(f"Flux to B1 (r2): {r2:.4e}")
print(f"Flux to Detemir (r3): {r3:.4e}")
```

```python
my_final_model = ExpandedKineticModel(final_params_full)
model_file_path = 'multi_ratio_kinetic_optimized_model.dill'

with open(model_file_path, 'wb') as f:
    dill.dump(my_final_model, f)
print(f"\nKinetic MODEL Multi-Rasio sukses disimpan ke '{model_file_path}'!")
```

```python
# Check total mass at t=0 and t=120
total_mass_t0 = np.sum([y0[model.idx[sp]] * MW[sp] for sp in model.species])
total_mass_t120 = np.sum([sol.y[model.idx[sp]][-1] * MW[sp] for sp in model.species])
print(f"Mass Balance Error: {abs(total_mass_t0 - total_mass_t120)/total_mass_t0 * 100:.2f}%")
```

```python
# 1. Provide a realistic set of starting k-values (not from the broken optimization)
# Let's start everything at a moderate rate (e.g., 0.1)
safe_k_values = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] 
n, m = 1.0, 1.0 # Locked to First-Order

safe_params = safe_k_values + [n, m]
test_model = FinalKineticModel(safe_params)

# 2. Setup initial conditions for Ratio 1.0
desb30_t0_molL = experimental_dataset[1.0]['desb30'][0]
y0 = np.zeros(len(test_model.species))
y0[test_model.idx['nhs_myr']] = desb30_t0_molL * 1.0
y0[test_model.idx['desb30']] = desb30_t0_molL

# 3. Simulate with LSODA (Best for stiff chemical networks)
sol = solve_ivp(test_model._ode_system, (0, 120), y0, method='LSODA', atol=1e-9, rtol=1e-7)

if sol.success:
    mass_in = np.sum([y0[test_model.idx[sp]] * MW[sp] for sp in test_model.species])
    mass_out = np.sum([sol.y[test_model.idx[sp]][-1] * MW[sp] for sp in test_model.species])
    print(f"Sanity Check: SUCCESS!")
    print(f"Mass Balance Error: {abs(mass_in - mass_out) / mass_in * 100:.5f}%")
else:
    print("Sanity Check: FAILED. ODEs are still exploding.")
```

```python
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np
import matplotlib.pyplot as plt
import dill

# ==============================================================================
# 1. LOAD & UNBOXING MODEL SECARA AMAN DARI FILE .DILL
# ==============================================================================
with open('30jul2026_hydrolysis_kinetic_model.dill', 'rb') as f:
    checkpoint = dill.load(f)

# Unboxing otomatis jika model dibungkus dalam dictionary
if isinstance(checkpoint, dict):
    print("INFO: File terbaca sebagai dictionary dengan keys:", list(checkpoint.keys()))
    if 'model' in checkpoint:
        model_loaded = checkpoint['model']
    elif 'model_loaded' in checkpoint:
        model_loaded = checkpoint['model_loaded']
    else:
        # Mencari objek di dalam dict yang memiliki method 'simulate'
        model_loaded = None
        for key, val in checkpoint.items():
            if hasattr(val, 'simulate'):
                print(f"INFO: Objek model ditemukan pada key: '{key}'")
                model_loaded = val
                break
        if model_loaded is None:
            raise AttributeError("Gagal mengekstrak objek model dari dictionary. Pastikan objek model memiliki method '.simulate'.")
else:
    model_loaded = checkpoint

# ==============================================================================
# 2. PROSES VISUALISASI GRAFIK & ANALISIS ERROR
# ==============================================================================
y_true_global, y_pred_global = [], []

plot_colors = {
    'desb30': '#1f77b4', 'A1': '#ff7f0e', 'detemir': '#2ca02c', 
    'A1_B29': '#d62728', 'B1_B29': '#9467bd', 'A1_B29_B1': '#8c564b', 
    'A1_B1': '#e377c2', 'B1': '#7f7f7f'
}

time_smooth = np.linspace(0, 120, 1000)
fig, axes = plt.subplots(1, 3, figsize=(22, 6.5), sharey=True)

print("\n=== ANALISIS ERROR PER SKENARIO RASIO ===")

for idx, r in enumerate(ratios):
    # DYNAMIC Y0: Menyelaraskan dengan logika fungsi objektif Anda
    desb30_t0_molL = experimental_dataset[r]['desb30'][0]
    
    # Ambil index dan ukuran spesies secara dinamis dari model yang di-load
    y0 = np.zeros(len(model_loaded.species))
    y0[model_loaded.idx['nhs_myr']] = desb30_t0_molL * r
    y0[model_loaded.idx['desb30']] = desb30_t0_molL
        
    # Jalankan simulasi menggunakan model dari file .dill
    sol_smooth = model_loaded.simulate(y0, (0, 120), time_smooth)
    sol_exp = model_loaded.simulate(y0, (0, 120), time_exp)
    
    # Penampung lokal untuk RMSE per rasio
    r_true, r_pred = [], []
    
    for sp in species_list:
        ode_pos = model_loaded.idx[sp]
        
        sim_gL_smooth = sol_smooth.y[ode_pos] * MW[sp]
        sim_gL_exp = sol_exp.y[ode_pos] * MW[sp]
        exp_gL = experimental_dataset[r][sp] * MW[sp]
        
        # --- FIX INDEX ERROR MASKING ---
        # Mengatasi inkonsistensi jumlah titik dari solver ODE dengan interpolasi linear ke time_exp
        if len(sim_gL_exp) != len(exp_gL):
            sim_gL_exp_fixed = np.interp(time_exp, time_smooth, sim_gL_smooth)
        else:
            sim_gL_exp_fixed = sim_gL_exp
            
        # Masking NaNs pada data eksperimen (hanya menguji data riil)
        mask = ~np.isnan(exp_gL)
        valid_exp = exp_gL[mask]
        valid_sim = sim_gL_exp_fixed[mask]
        
        r_true.extend(valid_exp)
        r_pred.extend(valid_sim)
        
        # Tambahkan ke global tracker
        y_true_global.extend(valid_exp)
        y_pred_global.extend(valid_sim)
        
        # Plotting hasil simulasi kurva mulus vs data eksperimen berupa titik scatter
        axes[idx].plot(time_smooth, sim_gL_smooth, color=plot_colors[sp], linewidth=1.6, alpha=0.95)
        axes[idx].scatter(time_exp, exp_gL, color=plot_colors[sp], edgecolors='white', linewidth=1.1, s=65, zorder=5)

    rmse_ratio = np.sqrt(mean_squared_error(r_true, r_pred))
    print(f"Rasio 1 : {r:.1f} -> RMSE = {rmse_ratio:.4f} g/L")
    
    ax = axes[idx]
    ax.set_title(f"Rasio 1 : {r:.1f}\n(RMSE: {rmse_ratio:.3f} g/L)", fontsize=12, fontweight='bold', pad=10)
    ax.set_xlabel("Waktu Reaksi (Menit)", fontsize=12)
    if idx == 0:
        ax.set_ylabel("Konsentrasi Komponen (g/L)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)

# Hitung Statistik Global Akhir
global_r2 = r2_score(y_true_global, y_pred_global)
global_rmse = np.sqrt(mean_squared_error(y_true_global, y_pred_global))

# Penataan Legend & Title
dummy_handles = [plt.Line2D([0], [0], color=plot_colors[sp], lw=2) for sp in species_list]
fig.legend(dummy_handles, species_list, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=4, fontsize=11, frameon=True)

plt.suptitle(f"Kurva Kinetika Global DE (Super Smooth)\nGlobal R² = {global_r2:.4f} ({global_r2*100:.2f}%)  |  Global RMSE = {global_rmse:.4f} g/L", 
             fontsize=14, fontweight='bold', y=1.06)

plt.tight_layout()
plt.show()

print("\n==================== RINGKASAN PERFORMA MODEL ====================")
print(f" >> Akurasi Tren (Global R²) : {global_r2:.4f} ({global_r2*100:.2f}%)")
print(f" >> Deviasi Rata-rata (Global RMSE) : {global_rmse:.4f} g/L")
print("==================================================================")
```
