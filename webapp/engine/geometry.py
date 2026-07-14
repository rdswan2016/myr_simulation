"""
Vessel/zone geometry for the N-compartment mass-transfer model.

Zone-splitting rules are generalized, unchanged, from the validated 2/3/4-compartment
notebooks in this repository (5L_tank_transfer_mass_3_compartments.ipynb):

  - Zone 0 (top): always the feed-entry AND probe zone. Spans from the top structural
    boundary (defined below) up to the liquid surface. Its volume GROWS during the
    titrant feed (surface rises) and is fixed afterward.
  - For n_zones >= 3, Zone 1 is the impeller-swept zone: one impeller diameter tall,
    centered on the impeller centerline (h_impeller +/- D_impeller/2). This is a
    structural modeling choice (not fit from data) -- the free parameters are always
    the inter-zone FLOW RATES, never these boundaries.
  - For n_zones == 4, Zone 2 is the "swept lower bulk" and Zone 3 is an explicit dead
    zone at the tank floor, D_impeller/4 tall (again a structural proxy, not fit).
  - For n_zones == 2, there is no separate impeller zone: Zone 0 spans impeller-height
    to surface, Zone 1 is everything below the impeller (the classic 2-zone baseline).

All zones below Zone 0 have FIXED volume; only Zone 0 grows with the titrant feed.
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class VesselGeometry:
    D_tank_cm: float
    D_impeller_cm: float
    h_impeller_cm: float
    h_liquid_initial_cm: float

    @property
    def A_cross_cm2(self) -> float:
        return np.pi * (self.D_tank_cm / 2.0) ** 2


@dataclass
class ZoneLayout:
    n_zones: int
    boundaries_cm: list  # [(bottom_cm, top_cm), ...] index 0 = top/feed zone ... index n-1 = bottom zone
    fixed_volumes_L: list  # fixed volume of each zone (index 0's value is its INITIAL volume, before feed)
    labels: list


def build_zone_layout(geom: VesselGeometry, n_zones: int) -> ZoneLayout:
    """Build zone boundaries/volumes for n_zones in {2, 3, 4}, reusing the exact
    structural rules validated in the notebook analyses."""
    if n_zones not in (2, 3, 4):
        raise ValueError("n_zones must be 2, 3, or 4")

    A = geom.A_cross_cm2
    D_imp = geom.D_impeller_cm
    h_imp = geom.h_impeller_cm
    h_liq = geom.h_liquid_initial_cm

    if n_zones == 2:
        top0 = h_imp
        boundaries = [(top0, h_liq), (0.0, top0)]
        labels = ["Zone 1 (top, feed + probe)", "Zone 2 (lower bulk)"]
    elif n_zones == 3:
        z_lower = h_imp - D_imp / 2.0
        z_upper = h_imp + D_imp / 2.0
        boundaries = [(z_upper, h_liq), (z_lower, z_upper), (0.0, z_lower)]
        labels = ["Zone 1 (top, feed + probe)", "Zone 2 (impeller zone)", "Zone 3 (lower bulk)"]
    else:  # n_zones == 4
        z_lower = h_imp - D_imp / 2.0
        z_upper = h_imp + D_imp / 2.0
        z_dead_top = D_imp / 4.0
        if z_dead_top >= z_lower:
            raise ValueError(
                "Dead-zone proxy height (D_impeller/4) does not fit below the impeller "
                "zone floor for this geometry -- reduce impeller height or use n_zones=3."
            )
        boundaries = [
            (z_upper, h_liq),
            (z_lower, z_upper),
            (z_dead_top, z_lower),
            (0.0, z_dead_top),
        ]
        labels = [
            "Zone 1 (top, feed + probe)",
            "Zone 2 (impeller zone)",
            "Zone 3 (swept lower bulk)",
            "Zone 4 (dead zone, tank floor)",
        ]

    if boundaries[-1][0] != 0.0 or boundaries[0][1] != h_liq:
        raise ValueError("Zone layout must span the full liquid column from 0 to the surface.")
    for (b0, t0), (b1, t1) in zip(boundaries[:-1], boundaries[1:]):
        if abs(b0 - t1) > 1e-9:
            raise ValueError("Zone boundaries must be contiguous (chain topology).")
        if t0 <= b0:
            raise ValueError(f"Zone height must be positive: {b0}-{t0} cm")

    volumes_L = [(A * (top - bottom)) / 1000.0 for (bottom, top) in boundaries]
    return ZoneLayout(n_zones=n_zones, boundaries_cm=boundaries, fixed_volumes_L=volumes_L, labels=labels)


def probe_zone_check(layout: ZoneLayout, probe_height_cm: float) -> bool:
    """True if the pH probe height falls inside Zone 0 (the feed zone), which is the
    physical assumption this whole model relies on (probe reads Zone-0 concentration)."""
    bottom0, top0 = layout.boundaries_cm[0]
    return bottom0 < probe_height_cm < top0 + 50.0  # top0 grows with feed; generous upper check only warns, doesn't block
