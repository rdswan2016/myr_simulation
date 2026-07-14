"""
General N-compartment convective transport engine, coupled to Henderson-Hasselbalch
chemistry (engine/chemistry.py) to produce a simulated pH(t) trace at the probe zone.

============================================================================
HOW THE HENDERSON-HASSELBALCH LOOP IS LINKED TO THE TRANSPORT ODE ENGINE
============================================================================
State vector (2 * n_zones values): for each zone i, track moles of the
CONJUGATE-BASE form (n_base_i) and CONJUGATE-ACID form (n_acid_i) SEPARATELY.
Chemistry and transport are coupled by three rules, applied once per RHS evaluation:

  1. TRANSPORT moves n_base and n_acid between adjacent zones independently, each
     by the same chain-diffusion law:
         flux_base(i -> i+1)  = Q[i] * (C_base_i  - C_base_{i+1})
         flux_acid(i -> i+1)  = Q[i] * (C_acid_i  - C_acid_{i+1})
     (concentrations, not moles, drive the flux -- this is a Fickian/well-mixed-zone
     exchange law, not a directed net flow.)

  2. CHEMISTRY (Henderson-Hasselbalch, chemistry.py:pH_from_moles) never appears
     inside the ODE right-hand side at all -- equilibrium is IMPLICIT and instantaneous
     because titrant is fed directly as a conversion of n_base -> n_acid (or vice
     versa) in the feed zone; no separate "reaction rate" term exists. This mirrors
     the validated notebook models: neutralization is assumed fast relative to
     circulation, so only transport is rate-limiting.

  3. OUTPUT: after solve_ivp integrates the transport-only ODE for
     (n_base_i(t), n_acid_i(t)) at every zone, chemistry.py:pH_from_moles(...) is
     called ONCE, only on the probe zone's (n_base_0, n_acid_0) trajectory, to map
     moles -> pH for comparison against the uploaded data. This is the only place
     pH is computed -- everywhere else the state is pure mole bookkeeping.

So schematically, per RHS call:

    for each zone i:
        C_base_i, C_acid_i = n_base_i / V_i(t), n_acid_i / V_i(t)      # concentrations
    for each junction (i, i+1):
        flux_base, flux_acid = Q[i]*(C_base_i-C_base_{i+1}), Q[i]*(C_acid_i-C_acid_{i+1})
        dn_base_i -= flux_base;  dn_base_{i+1} += flux_base            # (and acid form)
    dn_base_0   -= feed_rate(t)   if titrant delivers acid  else  += feed_rate(t)
    dn_acid_0   += feed_rate(t)   if titrant delivers acid  else  -= feed_rate(t)

    # ... solve_ivp integrates this (transport only) ...
    pH_probe(t) = chemistry.pH_from_moles(n_base_0(t), n_acid_0(t), pKa)   # <- H-H loop closes HERE
============================================================================
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from .chemistry import ChemistryParams, FeedParams, feed_rate_mol_min, initial_split, pH_from_moles
from .geometry import ZoneLayout


def _zone0_volume(t_min: float, V0_initial_L: float, feed: FeedParams) -> float:
    added = feed.pump_rate_L_min * min(t_min, feed.t_feed_end_min)
    return V0_initial_L + added


def simulate_ph(Q_L_min: np.ndarray, layout: ZoneLayout, chem: ChemistryParams,
                 feed: FeedParams, t_eval_min: np.ndarray, pH0: float) -> np.ndarray:
    """Simulate probe-zone (Zone 0) pH(t) for an n-zone chain with n-1 flow rates.

    Q_L_min: array of length (n_zones - 1); Q_L_min[i] is the exchange flow between
             zone i and zone i+1.
    """
    n = layout.n_zones
    assert len(Q_L_min) == n - 1

    V_fixed = layout.fixed_volumes_L  # index 0 is Zone-0's INITIAL volume only
    frac_acid_0 = initial_split(pH0, chem.pKa)

    n_base0_state = np.array([(1 - frac_acid_0) * chem.C_analyte_stock_M * V for V in V_fixed])
    n_acid0_state = np.array([frac_acid_0 * chem.C_analyte_stock_M * V for V in V_fixed])
    y0 = np.empty(2 * n)
    y0[0::2] = n_base0_state
    y0[1::2] = n_acid0_state

    feed_sign = 1.0 if chem.titrant_delivers_acid else -1.0  # +1: feed converts base->acid in Zone 0

    def ode(t_min, y):
        n_base = np.maximum(y[0::2], 0.0)
        n_acid = np.maximum(y[1::2], 0.0)
        V = np.array(V_fixed, dtype=float)
        V[0] = _zone0_volume(t_min, V_fixed[0], feed)
        C_base = n_base / V
        C_acid = n_acid / V

        d_base = np.zeros(n)
        d_acid = np.zeros(n)
        for i in range(n - 1):
            flux_base = Q_L_min[i] * (C_base[i] - C_base[i + 1])
            flux_acid = Q_L_min[i] * (C_acid[i] - C_acid[i + 1])
            d_base[i] -= flux_base
            d_base[i + 1] += flux_base
            d_acid[i] -= flux_acid
            d_acid[i + 1] += flux_acid

        h_feed = feed_rate_mol_min(t_min, feed, chem)
        d_base[0] -= feed_sign * h_feed
        d_acid[0] += feed_sign * h_feed

        dydt = np.empty(2 * n)
        dydt[0::2] = d_base
        dydt[1::2] = d_acid
        return dydt

    t_feed_end = feed.t_feed_end_min
    t_end = t_eval_min[-1]
    if t_end <= t_feed_end:
        sol = solve_ivp(ode, (0, t_end), y0, t_eval=t_eval_min, method='BDF', rtol=1e-7, atol=1e-10)
        Y = sol.y
    else:
        mask1 = t_eval_min <= t_feed_end
        t1 = np.append(t_eval_min[mask1], t_feed_end) if t_feed_end not in t_eval_min[mask1] else t_eval_min[mask1]
        sol1 = solve_ivp(ode, (0, t_feed_end), y0, t_eval=t1, method='BDF', rtol=1e-7, atol=1e-10)
        y_mid = sol1.y[:, -1]
        t2 = t_eval_min[~mask1]
        sol2 = solve_ivp(ode, (t_feed_end, t_end), y_mid, t_eval=t2, method='BDF', rtol=1e-7, atol=1e-10)
        Y = np.hstack([sol1.y[:, :len(t_eval_min[mask1])], sol2.y])

    n_base0_t, n_acid0_t = Y[0], Y[1]
    return pH_from_moles(n_base0_t, n_acid0_t, chem.pKa)  # <-- Henderson-Hasselbalch loop closes here


def physical_flow_ceiling_L_min(N_rps: float, D_impeller_cm: float, Nq_max: float = 1.0) -> float:
    """Generous literature ceiling on ANY real circulation flow this impeller can drive,
    used as the upper optimizer bound for every Q (a conservative, engineering-defensible
    bound rather than an arbitrary large box constraint -- see notebook precedent)."""
    D_m = D_impeller_cm / 100.0
    return Nq_max * N_rps * (D_m ** 3) * 1000.0 * 60.0


@dataclass
class FitResult:
    n_zones: int
    Q_fit: np.ndarray
    cost: float
    rss: float
    rmse: float
    n_params: int
    n_points: int
    at_bound: np.ndarray  # bool per parameter: fitted value sits at its upper bound (non-identifiable flag)
    success: bool
    message: str


def fit_n_zone_model(layout: ZoneLayout, chem: ChemistryParams, feed: FeedParams,
                      t_sec: np.ndarray, pH_data: np.ndarray, Q_ceiling_L_min: float) -> FitResult:
    n = layout.n_zones
    k = n - 1
    V_total0 = sum(layout.fixed_volumes_L)
    x0 = np.array([V_total0 * 3.0 * (0.2 ** i) for i in range(k)])
    x0 = np.clip(x0, 1e-3, Q_ceiling_L_min * 0.9)
    lower = np.full(k, 1e-6)
    upper = np.full(k, Q_ceiling_L_min)

    def residuals(Q):
        pH0 = pH_data[0]
        pH_sim = simulate_ph(Q, layout, chem, feed, t_sec / 60.0, pH0)
        return pH_sim - pH_data

    result = least_squares(residuals, x0=x0, bounds=(lower, upper), method='trf', xtol=1e-13, ftol=1e-13)
    rss = float(np.sum(result.fun ** 2))
    rmse = float(np.sqrt(np.mean(result.fun ** 2)))
    at_bound = np.isclose(result.x, upper, rtol=1e-3)

    return FitResult(
        n_zones=n, Q_fit=result.x, cost=float(result.cost), rss=rss, rmse=rmse,
        n_params=k, n_points=len(pH_data), at_bound=at_bound,
        success=result.status > 0, message=result.message,
    )
