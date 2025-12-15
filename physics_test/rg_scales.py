from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class QCDLambdaResult:
    n_f: int
    mu_GeV: float
    alpha_s_mu: float
    beta0: float
    beta1: float
    loops: int
    Lambda_GeV: float


def qcd_beta0(n_f: int) -> float:
    """QCD beta0 for SU(3): β0 = 11 - 2 n_f/3."""
    return 11.0 - (2.0 / 3.0) * float(n_f)


def qcd_beta1(n_f: int) -> float:
    """QCD beta1 for SU(3): β1 = 102 - 38 n_f/3 (MSbar)."""
    return 102.0 - (38.0 / 3.0) * float(n_f)


def lambda_qcd_from_alpha_s(
    *,
    alpha_s_mu: float,
    mu_GeV: float,
    n_f: int,
    loops: int = 2,
) -> QCDLambdaResult:
    """
    Dimensional transmutation (RG invariant) scale for QCD, Λ_QCD, from α_s(μ).

    Conventions:
      We use a common 1–2 loop approximation written in terms of α_s:

        dα/d ln μ = - (β0/(2π)) α^2 - (β1/(4π^2)) α^3 + ...

      Then the RG-invariant Λ (1-loop) is defined by:
        ln(μ/Λ) = 2π/(β0 α(μ))
        Λ = μ * exp(-2π/(β0 α(μ)))

      And the 2-loop improved expression:
        Λ = μ * exp(-2π/(β0 α)) * (β0 α/(2π))^{-β1/β0^2}

    This is still approximate (threshold matching + higher loops ignored).
    """

    if mu_GeV <= 0:
        raise ValueError("mu_GeV must be positive")
    if alpha_s_mu <= 0:
        raise ValueError("alpha_s_mu must be positive")
    if loops not in (1, 2):
        raise ValueError("loops must be 1 or 2")

    b0 = qcd_beta0(n_f)
    b1 = qcd_beta1(n_f)
    if b0 <= 0:
        # For QCD in the SM range this should be positive; guard anyway.
        raise ValueError("beta0 must be positive for asymptotic freedom")

    a = float(alpha_s_mu)
    mu = float(mu_GeV)

    # 1-loop
    Lambda = mu * math.exp(-(2.0 * math.pi) / (b0 * a))

    if loops == 2:
        # 2-loop correction factor
        x = (b0 * a) / (2.0 * math.pi)
        if x <= 0:
            raise ValueError("Invalid x for 2-loop Lambda computation")
        Lambda = Lambda * (x ** (-(b1 / (b0 * b0))))

    return QCDLambdaResult(
        n_f=int(n_f),
        mu_GeV=mu,
        alpha_s_mu=a,
        beta0=b0,
        beta1=b1,
        loops=int(loops),
        Lambda_GeV=float(Lambda),
    )


