"""
Henderson-Hasselbalch thermodynamics + continuous-titrant feed model.

This module is deliberately the ONLY place in the codebase that knows about
weak-acid/base equilibrium chemistry. transport_model.py calls it once per zone,
per ODE evaluation, to convert (moles of conjugate base, moles of conjugate acid)
into a pH -- the transport engine itself only ever moves "base form" and "acid form"
moles between zones by bulk flow. See engine/transport_model.py's module docstring
for the exact coupling diagram.
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class ChemistryParams:
    pKa: float
    C_analyte_stock_M: float      # stock concentration of the buffering species initially in the vessel
    C_titrant_stock_M: float      # stock concentration of the titrant being pumped in
    titrant_delivers_acid: bool   # True: titrant is an acid (adds H+, e.g. HCl into Tris base).
                                    # False: titrant is a base (adds OH-/consumes H+, e.g. NaOH into a weak acid).


@dataclass
class FeedParams:
    V_titrant_total_mL: float
    pump_rate_mL_min: float

    @property
    def pump_rate_L_min(self) -> float:
        return self.pump_rate_mL_min / 1000.0

    @property
    def t_feed_end_min(self) -> float:
        return self.V_titrant_total_mL / self.pump_rate_mL_min


def initial_split(pH0: float, pKa: float) -> float:
    """Inverse Henderson-Hasselbalch: fraction of total analyte in the CONJUGATE-ACID
    form at t=0, given a measured starting pH. Returns frac_acid_form in [0, 1]."""
    ratio_base_over_acid = 10 ** (pH0 - pKa)
    frac_acid_form = 1.0 / (1.0 + ratio_base_over_acid)
    return frac_acid_form


def pH_from_moles(n_base: np.ndarray, n_acid: np.ndarray, pKa: float, floor: float = 1e-12) -> np.ndarray:
    """Henderson-Hasselbalch: pH = pKa + log10([base]/[acid]). Concentrations cancel to a
    mole ratio directly (same zone volume in numerator and denominator), so this can be
    called on raw mole arrays without dividing by volume first."""
    b = np.maximum(n_base, floor)
    a = np.maximum(n_acid, floor)
    return pKa + np.log10(b / a)


def feed_rate_mol_min(t_min: float, feed: FeedParams, chem: ChemistryParams) -> float:
    """Moles of titrant delivered per minute at time t (0 once the pump has stopped)."""
    if t_min <= feed.t_feed_end_min:
        return feed.pump_rate_L_min * chem.C_titrant_stock_M
    return 0.0


def mass_balance_ceiling_pH(total_analyte_mol: float, total_titrant_mol: float,
                             V_total_final_L: float, chem: ChemistryParams) -> float:
    """The one pH value a CLOSED system (titrant feed finished, no more mass entering)
    must asymptote to, regardless of the number of compartments or their flow rates --
    set purely by total charge and final volume. Any compartmental transport model's
    achievable fit quality is upper-bounded by this value; a persistent gap between a
    data peak/plateau and this ceiling indicates a non-transport cause (e.g. CO2
    outgassing, reagent-lot deviation from the nominal stock concentration), not a
    transport-model shortfall -- no amount of extra compartments can close it."""
    if chem.titrant_delivers_acid:
        acid_form_mol = total_titrant_mol
        base_form_mol = total_analyte_mol - total_titrant_mol
    else:
        base_form_mol = total_titrant_mol
        acid_form_mol = total_analyte_mol - total_titrant_mol
    C_acid = max(acid_form_mol, 1e-12) / V_total_final_L
    C_base = max(base_form_mol, 1e-12) / V_total_final_L
    return chem.pKa + np.log10(C_base / C_acid)
