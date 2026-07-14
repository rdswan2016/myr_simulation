"""
Automatic compartment-count selection.

Fits the 2-, 3-, and 4-compartment models to the same dataset and picks the smallest
compartment count that is not significantly improved upon by the next tier, using:

  1. A nested extra-sum-of-squares F-test (Q_i+1 vs Q_i) -- the statistically correct
     way to ask "does adding a compartment explain significantly more variance," as
     opposed to just comparing raw RSS (which can only ever improve or tie as
     parameters are added, so it can never by itself justify added complexity).
  2. An "at-bound" identifiability flag -- if the added compartment's own flow rate(s)
     converge to the edge of their allowed range rather than an interior optimum, that
     compartment's rate is not actually pinned down by the data, regardless of what
     the F-test says about aggregate fit quality.
  3. A closed-system mass-balance ceiling check -- independent of compartment count
     entirely. If the data's peak/plateau pH exceeds the theoretical fully-mixed
     equilibrium implied by the total titrant/analyte charge, NO transport model of
     ANY compartment count can close that gap, and the rationale says so explicitly
     instead of recommending more compartments to chase an unreachable target.

AIC/BIC are also reported (for users who want the classic criteria alongside the more
rigorous nested test), but are not what the automatic recommendation is actually based
on -- with nested models fit by unweighted least squares, AIC/BIC can favor an
over-parameterized model whose extra flow rate is not statistically distinguishable
from a boundary value, exactly the failure mode the F-test + at-bound check are here
to catch.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.stats import f as f_dist

from .chemistry import ChemistryParams, FeedParams, mass_balance_ceiling_pH
from .geometry import VesselGeometry, build_zone_layout
from .transport_model import FitResult, fit_n_zone_model, physical_flow_ceiling_L_min

ALPHA = 0.05  # significance level for the nested F-test


@dataclass
class ModelSelectionResult:
    fits: dict  # {n_zones: FitResult}
    layouts: dict  # {n_zones: ZoneLayout}
    aic: dict
    bic: dict
    f_test_3v2: dict  # {'F':..., 'p':..., 'significant':...}
    f_test_4v3: dict
    recommended_n: int
    ceiling_pH: float
    data_peak_pH: float
    ceiling_exceeded: bool
    rationale: str = field(default="")


def _aic_bic(rss: float, n_points: int, n_params: int):
    rss = max(rss, 1e-15)
    aic = n_points * np.log(rss / n_points) + 2 * n_params
    bic = n_points * np.log(rss / n_points) + n_params * np.log(n_points)
    return aic, bic


def _nested_f_test(fit_reduced: FitResult, fit_full: FitResult):
    df_reduced = fit_reduced.n_points - fit_reduced.n_params
    df_full = fit_full.n_points - fit_full.n_params
    d_df = fit_full.n_params - fit_reduced.n_params
    rss_reduced, rss_full = fit_reduced.rss, fit_full.rss
    if rss_full <= 0 or df_full <= 0 or rss_reduced < rss_full:
        rss_full = max(rss_full, 1e-15)
    numerator = max(rss_reduced - rss_full, 0.0) / d_df
    denominator = rss_full / df_full
    F = numerator / denominator if denominator > 0 else 0.0
    p = float(f_dist.sf(F, d_df, df_full))
    return {"F": float(F), "p": p, "significant": p < ALPHA, "df1": d_df, "df2": df_full}


def run_model_selection(geom: VesselGeometry, N_rps: float, chem: ChemistryParams, feed: FeedParams,
                         t_sec: np.ndarray, pH_data: np.ndarray, Nq_max: float = 1.0) -> ModelSelectionResult:
    Q_ceiling = physical_flow_ceiling_L_min(N_rps, geom.D_impeller_cm, Nq_max)

    layouts, fits, aic, bic = {}, {}, {}, {}
    for n in (2, 3, 4):
        layout = build_zone_layout(geom, n)
        fit = fit_n_zone_model(layout, chem, feed, t_sec, pH_data, Q_ceiling)
        a, b = _aic_bic(fit.rss, fit.n_points, fit.n_params)
        layouts[n], fits[n], aic[n], bic[n] = layout, fit, a, b

    f_3v2 = _nested_f_test(fits[2], fits[3])
    f_4v3 = _nested_f_test(fits[3], fits[4])

    # --- mass-balance ceiling: independent of n_zones, computed once ---
    V_total_final_L = sum(build_zone_layout(geom, 3).fixed_volumes_L[1:]) + (
        build_zone_layout(geom, 3).fixed_volumes_L[0] + feed.V_titrant_total_mL / 1000.0
    )
    V_total_initial_L = sum(build_zone_layout(geom, 3).fixed_volumes_L)
    total_analyte_mol = chem.C_analyte_stock_M * V_total_initial_L
    total_titrant_mol = chem.C_titrant_stock_M * (feed.V_titrant_total_mL / 1000.0)
    ceiling_pH = mass_balance_ceiling_pH(total_analyte_mol, total_titrant_mol, V_total_final_L, chem)

    t_feed_end_s = feed.t_feed_end_min * 60.0
    mask_post_feed = t_sec > t_feed_end_s
    data_peak_pH = float(pH_data[mask_post_feed].max()) if mask_post_feed.any() else float(pH_data.max())
    if chem.titrant_delivers_acid:
        ceiling_exceeded = data_peak_pH > ceiling_pH + 0.01
    else:
        data_trough_pH = float(pH_data[mask_post_feed].min()) if mask_post_feed.any() else float(pH_data.min())
        ceiling_exceeded = data_trough_pH < ceiling_pH - 0.01

    # --- recommendation: smallest n not significantly (and identifiably) improved on ---
    recommended_n = 2
    if f_3v2["significant"] and not fits[3].at_bound.any():
        recommended_n = 3
        if f_4v3["significant"] and not fits[4].at_bound.any():
            recommended_n = 4

    rationale = _build_rationale(fits, f_3v2, f_4v3, recommended_n, ceiling_pH, data_peak_pH,
                                  ceiling_exceeded, chem)

    return ModelSelectionResult(
        fits=fits, layouts=layouts, aic=aic, bic=bic, f_test_3v2=f_3v2, f_test_4v3=f_4v3,
        recommended_n=recommended_n, ceiling_pH=ceiling_pH, data_peak_pH=data_peak_pH,
        ceiling_exceeded=ceiling_exceeded, rationale=rationale,
    )


def _build_rationale(fits, f_3v2, f_4v3, recommended_n, ceiling_pH, data_peak_pH, ceiling_exceeded, chem) -> str:
    lines = []
    lines.append(f"**Recommended layout: {recommended_n} compartments.**")
    lines.append("")

    if f_3v2["significant"]:
        if fits[3].at_bound.any():
            lines.append(
                f"- A 3-compartment layout explains significantly more of the pH trace than 2 "
                f"compartments (nested F-test p={f_3v2['p']:.4g}), **but** its added flow rate "
                f"converged to the edge of its allowed range rather than a stable interior value "
                f"-- meaning the data cannot actually pin down *how fast* that third zone exchanges, "
                f"only that treating it separately fits the overall shape better. This is reported as "
                f"a caveat, not treated as evidence for a genuine third, slow, identifiable zone."
            )
        else:
            lines.append(
                f"- A 3-compartment layout was selected over 2 because it captures a feature the "
                f"2-compartment model cannot: the nested F-test shows a statistically significant "
                f"improvement (p={f_3v2['p']:.4g}), and the added flow rate converged to a stable, "
                f"non-boundary value -- consistent with a real, identifiable secondary exchange "
                f"pathway (e.g. a stagnant zone below the impeller producing a secondary pH shock)."
            )
    else:
        lines.append(
            f"- A 2-compartment layout was selected: adding a 3rd compartment did **not** "
            f"significantly improve the fit (nested F-test p={f_3v2['p']:.4g}, not below 0.05). "
            f"The simplest model already explains the data as well as a more complex one can."
        )

    if recommended_n >= 3:
        if f_4v3["significant"] and not fits[4].at_bound.any():
            lines.append(
                f"- A 4th compartment (explicit dead zone) further improves the fit significantly "
                f"(p={f_4v3['p']:.4g}) with an identifiable (non-boundary) flow rate, and was kept."
            )
        elif f_4v3["significant"] and fits[4].at_bound.any():
            lines.append(
                f"- A 4th compartment shows a nominally significant improvement (p={f_4v3['p']:.4g}) "
                f"but its flow rate is not identifiable (converged to a bound) -- **not** adopted; "
                f"the 4-compartment fit is reported for reference only."
            )
        else:
            lines.append(
                f"- A 4th compartment (explicit dead zone) was tested and did **not** significantly "
                f"improve the fit beyond 3 compartments (p={f_4v3['p']:.4g}) -- not adopted."
            )

    lines.append("")
    if ceiling_exceeded:
        gap = data_peak_pH - ceiling_pH if chem.titrant_delivers_acid else ceiling_pH - data_peak_pH
        lines.append(
            f"- **Important:** the data's post-feed peak pH ({data_peak_pH:.3f}) exceeds this closed "
            f"system's own mass-balance ceiling ({ceiling_pH:.3f}, the fully-mixed equilibrium set "
            f"purely by total titrant/analyte charge and final volume) by {gap:.3f} pH units. "
            f"**No compartmental transport model, at any compartment count, can close this specific "
            f"gap** -- it is a hard analytical bound independent of flow rates or zone count. If your "
            f"data shows a late secondary drift that this tool's recommended model does not fully "
            f"track, treat it as a signal of a non-transport cause (CO2 outgassing, reagent-lot "
            f"concentration deviation, temperature drift, electrode drift) rather than adding more "
            f"compartments to chase it."
        )
    else:
        lines.append(
            f"- The mass-balance ceiling check passed: the data's post-feed extreme pH "
            f"({data_peak_pH:.3f}) does not exceed the closed-system equilibrium ({ceiling_pH:.3f}), "
            f"so a transport-only explanation remains plausible for the full trace."
        )

    return "\n".join(lines)
