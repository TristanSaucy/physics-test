from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class BetaCoefficients:
    """
    One-loop beta function coefficients b_i for gauge couplings in:

      d(α_i^{-1}) / d ln μ = - b_i / (2π)

    Conventions vary; this module is only for exploratory convergence testing.
    """

    b1: float  # U(1)_Y (GUT-normalized)
    b2: float  # SU(2)
    b3: float  # SU(3)
    name: str = "custom"


SM_1LOOP = BetaCoefficients(b1=41.0 / 10.0, b2=-19.0 / 6.0, b3=-7.0, name="SM (1-loop)")
MSSM_1LOOP = BetaCoefficients(b1=33.0 / 5.0, b2=1.0, b3=-3.0, name="MSSM (1-loop)")


def run_alpha_inv(alpha_inv_mu0: float, mu0: float, mu: float, b: float) -> float:
    """1-loop running: α^{-1}(μ) = α^{-1}(μ0) - (b/(2π)) ln(μ/μ0)."""
    return float(alpha_inv_mu0) - (float(b) / (2.0 * math.pi)) * math.log(float(mu) / float(mu0))


def converge_score(alpha1_inv: float, alpha2_inv: float, alpha3_inv: float) -> float:
    """Simple convergence score: max pairwise difference in α^{-1}."""
    return max(abs(alpha1_inv - alpha2_inv), abs(alpha1_inv - alpha3_inv), abs(alpha2_inv - alpha3_inv))


def find_best_convergence(
    *,
    mu0: float,
    alpha1_inv_mu0: float,
    alpha2_inv_mu0: float,
    alpha3_inv_mu0: float,
    betas: BetaCoefficients,
    mu_min: float,
    mu_max: float,
    n: int = 2000,
) -> tuple[float, float, float, float, float]:
    """
    Scan log-spaced mu in [mu_min, mu_max] and return best point:

    (mu_best, score_best, a1_inv, a2_inv, a3_inv)
    """
    if mu_min <= 0 or mu_max <= 0 or mu_min >= mu_max:
        raise ValueError("mu_min/mu_max must be positive and mu_min < mu_max")
    if n < 2:
        raise ValueError("n must be >= 2")

    log_min = math.log(mu_min)
    log_max = math.log(mu_max)

    best_mu = None
    best_score = float("inf")
    best_vals = (float("nan"), float("nan"), float("nan"))

    for i in range(n):
        t = i / (n - 1)
        mu = math.exp(log_min + (log_max - log_min) * t)
        a1 = run_alpha_inv(alpha1_inv_mu0, mu0, mu, betas.b1)
        a2 = run_alpha_inv(alpha2_inv_mu0, mu0, mu, betas.b2)
        a3 = run_alpha_inv(alpha3_inv_mu0, mu0, mu, betas.b3)
        score = converge_score(a1, a2, a3)
        if score < best_score:
            best_score = score
            best_mu = mu
            best_vals = (a1, a2, a3)

    assert best_mu is not None
    a1, a2, a3 = best_vals
    return best_mu, best_score, a1, a2, a3


