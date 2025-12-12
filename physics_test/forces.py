from __future__ import annotations

import math

from physics_test import constants


def alpha_gravity(mass_kg: float) -> float:
    """
    Dimensionless gravitational coupling at a chosen mass scale:

        alpha_G(m) = G_N * m^2 / (hbar * c)
    """
    return constants.G_NEWTON * (mass_kg**2) / (constants.HBAR * constants.SPEED_OF_LIGHT)


def alpha_s_1loop(Q_GeV: float, *, Lambda_GeV: float = 0.2, n_f: int = 5) -> float:
    """
    Very rough 1-loop QCD running coupling (toy model):

        alpha_s(Q) = 12*pi / ((33 - 2 n_f) * ln(Q^2 / Lambda^2))

    This ignores thresholds and higher loops; it's only for broad exploratory scanning.
    """
    if Q_GeV <= 0 or Lambda_GeV <= 0:
        raise ValueError("Q_GeV and Lambda_GeV must be positive")
    b0 = 33.0 - 2.0 * float(n_f)
    x = (Q_GeV / Lambda_GeV) ** 2
    if x <= 1.0:
        # Below Lambda this blows up; return inf to signal non-perturbative region.
        return float("inf")
    return 12.0 * math.pi / (b0 * math.log(x))


