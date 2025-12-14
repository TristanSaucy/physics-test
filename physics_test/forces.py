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


def alpha_s_run_1loop_from_ref(
    Q_GeV: float,
    *,
    alpha_s_Q0: float,
    Q0_GeV: float,
    n_f: int = 5,
) -> float:
    """
    1-loop running for QCD using a *reference value* alpha_s(Q0) instead of a free Lambda:

      α_s^{-1}(Q) = α_s^{-1}(Q0) + (b0/(2π)) ln(Q/Q0),
      b0 = 11 - 2 n_f / 3.

    This is still approximate (thresholds + higher loops ignored), but it removes the
    "free Lambda" knob and is therefore a cleaner target-definition refinement.
    """

    if Q_GeV <= 0 or Q0_GeV <= 0:
        raise ValueError("Q_GeV and Q0_GeV must be positive")
    if alpha_s_Q0 <= 0:
        raise ValueError("alpha_s_Q0 must be positive")
    b0 = 11.0 - (2.0 / 3.0) * float(n_f)
    alpha_inv_Q = (1.0 / float(alpha_s_Q0)) + (b0 / (2.0 * math.pi)) * math.log(float(Q_GeV) / float(Q0_GeV))
    if alpha_inv_Q <= 0:
        # Would indicate hitting the Landau pole / non-perturbative region in this toy model.
        return float("inf")
    return 1.0 / alpha_inv_Q


