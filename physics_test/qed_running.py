from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class ChargedParticle:
    """
    Minimal charged-matter entry for 1-loop QED running.

    We model the 1-loop beta function contribution from charged Dirac fermions:

      dα/d ln Q = (2 α^2 / (3π)) * Σ (N_c Q^2)

    which implies:

      d(α^{-1})/d ln Q = -(2 / (3π)) * Σ (N_c Q^2)

    Notes:
      - This is an intentionally simple toy prescription.
      - It does NOT model non-perturbative hadronic vacuum polarization precisely.
      - Threshold handling is done by a sharp turn-on at the particle mass.
    """

    name: str
    mass_GeV: float
    charge: float  # electric charge in units of e
    color: int = 1

    @property
    def weight(self) -> float:
        return float(self.color) * float(self.charge) * float(self.charge)


# Rough PDG-ish masses (GeV). These are only used for deterministic threshold placement.
# Quark masses here are current-quark / MSbar-ish order-of-magnitude values.
DEFAULT_CHARGED_PARTICLES: list[ChargedParticle] = [
    # leptons
    ChargedParticle("e", 0.00051099895, -1.0, 1),
    ChargedParticle("mu", 0.1056583745, -1.0, 1),
    ChargedParticle("tau", 1.77686, -1.0, 1),
    # quarks (color=3)
    ChargedParticle("u", 0.0022, +2.0 / 3.0, 3),
    ChargedParticle("d", 0.0047, -1.0 / 3.0, 3),
    ChargedParticle("s", 0.096, -1.0 / 3.0, 3),
    ChargedParticle("c", 1.27, +2.0 / 3.0, 3),
    ChargedParticle("b", 4.18, -1.0 / 3.0, 3),
    ChargedParticle("t", 172.76, +2.0 / 3.0, 3),
]


def qed_run_alpha_inv_1loop_from_ref(
    Q_GeV: float,
    *,
    alpha_inv_Q0: float,
    Q0_GeV: float,
    particles: list[ChargedParticle] | None = None,
) -> float:
    """
    Deterministic 1-loop QED running with sharp thresholds.

    Input:
      - alpha_inv_Q0 = 1/alpha(Q0)
      - Q0_GeV: reference scale in GeV (must be > 0)
      - Q_GeV: target scale in GeV (must be > 0)

    Output:
      - alpha_inv(Q) in the same 1-loop threshold approximation.
    """

    if Q0_GeV <= 0 or Q_GeV <= 0:
        raise ValueError("Q0_GeV and Q_GeV must be positive")

    if alpha_inv_Q0 <= 0:
        raise ValueError("alpha_inv_Q0 must be positive")

    parts = particles if particles is not None else DEFAULT_CHARGED_PARTICLES
    k = 2.0 / (3.0 * math.pi)  # coefficient for d(1/alpha)/d ln Q

    Q0 = float(Q0_GeV)
    Q = float(Q_GeV)
    out = float(alpha_inv_Q0)

    if Q == Q0:
        return out

    if Q > Q0:
        # Run upward: particles contribute only above their mass threshold.
        for p in parts:
            m = float(p.mass_GeV)
            if m <= 0:
                continue
            if Q <= m:
                continue
            start = max(Q0, m)
            out -= k * p.weight * math.log(Q / start)
        return out

    # Run downward
    for p in parts:
        m = float(p.mass_GeV)
        if m <= 0:
            continue
        if Q0 <= m:
            continue
        end = max(Q, m)
        out -= k * p.weight * math.log(end / Q0)  # log(end/Q0) < 0 -> increases out
    return out


