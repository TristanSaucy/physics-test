from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class BetaCoefficients:
    """
    1-loop beta coefficients b_i for gauge couplings in convention:

        d(1/alpha_i)/d ln(Q) = - b_i / (2*pi)

    where alpha_i = g_i^2/(4*pi).
    """

    b1: float
    b2: float
    b3: float


SM_1LOOP = BetaCoefficients(b1=41.0 / 10.0, b2=-19.0 / 6.0, b3=-7.0)
MSSM_1LOOP = BetaCoefficients(b1=33.0 / 5.0, b2=1.0, b3=-3.0)


@dataclass(frozen=True)
class RunningPoint:
    Q_GeV: float
    inv_a1: float
    inv_a2: float
    inv_a3: float

    @property
    def a1(self) -> float:
        return 1.0 / self.inv_a1

    @property
    def a2(self) -> float:
        return 1.0 / self.inv_a2

    @property
    def a3(self) -> float:
        return 1.0 / self.inv_a3


def run_inv_alpha_1loop(inv_a0: float, *, b: float, Q0: float, Q: float) -> float:
    """1-loop running of 1/alpha."""
    if Q <= 0 or Q0 <= 0:
        raise ValueError("Q and Q0 must be positive")
    return inv_a0 - (b / (2.0 * math.pi)) * math.log(Q / Q0)


def run_gut_1loop(
    *,
    Q_GeV: float,
    Q0_GeV: float,
    inv_a1_0: float,
    inv_a2_0: float,
    inv_a3_0: float,
    betas: BetaCoefficients,
) -> RunningPoint:
    inv_a1 = run_inv_alpha_1loop(inv_a1_0, b=betas.b1, Q0=Q0_GeV, Q=Q_GeV)
    inv_a2 = run_inv_alpha_1loop(inv_a2_0, b=betas.b2, Q0=Q0_GeV, Q=Q_GeV)
    inv_a3 = run_inv_alpha_1loop(inv_a3_0, b=betas.b3, Q0=Q0_GeV, Q=Q_GeV)
    return RunningPoint(Q_GeV=Q_GeV, inv_a1=inv_a1, inv_a2=inv_a2, inv_a3=inv_a3)


def convergence_score(point: RunningPoint) -> float:
    """
    Simple convergence score: RMS of pairwise differences in inverse couplings.
    Lower is "more unified".
    """
    d12 = point.inv_a1 - point.inv_a2
    d23 = point.inv_a2 - point.inv_a3
    d13 = point.inv_a1 - point.inv_a3
    return math.sqrt((d12 * d12 + d23 * d23 + d13 * d13) / 3.0)


