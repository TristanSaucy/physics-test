from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class EWChargedFermion:
    """
    Minimal charged-fermion entry for a toy γ–Z mixing / κ(Q) style running model.

    This is NOT a full electroweak calculation. The goal is a deterministic,
    threshold-aware log model that can be used as a principled next lever when
    comparing to low-Q effective weak mixing angle extractions (Qweak, E158, etc.).

    Fields:
      - charge: electric charge in units of e (Q_f)
      - T3: weak isospin of the left-handed component (±1/2). Right-handed has T3=0.
      - color: N_c
    """

    name: str
    mass_GeV: float
    charge: float
    T3: float
    color: int = 1

    @property
    def weight_QT3(self) -> float:
        return float(self.color) * float(self.charge) * float(self.T3)


# Rough masses (GeV) for deterministic threshold placement only.
# We intentionally skip the top quark here (m_t > mZ) by default.
DEFAULT_EW_FERMIONS: list[EWChargedFermion] = [
    # leptons (charged)
    EWChargedFermion("e", 0.00051099895, -1.0, -0.5, 1),
    EWChargedFermion("mu", 0.1056583745, -1.0, -0.5, 1),
    EWChargedFermion("tau", 1.77686, -1.0, -0.5, 1),
    # quarks (color=3) (current-quark / MSbar-ish order-of-magnitude masses)
    EWChargedFermion("u", 0.0022, +2.0 / 3.0, +0.5, 3),
    EWChargedFermion("d", 0.0047, -1.0 / 3.0, -0.5, 3),
    EWChargedFermion("s", 0.096, -1.0 / 3.0, -0.5, 3),
    EWChargedFermion("c", 1.27, +2.0 / 3.0, +0.5, 3),
    EWChargedFermion("b", 4.18, -1.0 / 3.0, -0.5, 3),
]


def kappa_gammaZ_1loop_from_ref(
    Q_GeV: float,
    *,
    sin2_ref: float,
    alpha_ref: float,
    Q0_GeV: float,
    fermions: list[EWChargedFermion] | None = None,
) -> float:
    """
    Toy low-Q κ(Q) model for the effective weak mixing angle, based on a 1-loop,
    threshold-aware logarithmic γ–Z mixing ansatz.

    This intentionally avoids trying to run unbroken SU(2) down to low Q (which is not
    the same object as the experimentally extracted effective weak angle).

    Model:
      κ(Q) = 1 + (α_ref/π) * Σ_f (N_c Q_f T3_f) * ln(Q0 / max(Q, m_f))

    Notes:
      - sin2_ref is included in the signature because more accurate weights depend on it;
        the current toy uses only Q*T3 weights (dominant structure in many leading-log
        treatments).
      - For Q > Q0, this returns κ < 1 via negative logs (not the intended regime).
    """

    _ = float(sin2_ref)  # reserved for future refinements
    if Q0_GeV <= 0 or Q_GeV <= 0:
        raise ValueError("Q0_GeV and Q_GeV must be positive")
    if alpha_ref <= 0:
        raise ValueError("alpha_ref must be positive")

    parts = fermions if fermions is not None else DEFAULT_EW_FERMIONS
    Q0 = float(Q0_GeV)
    Q = float(Q_GeV)

    s = 0.0
    for f in parts:
        m = float(f.mass_GeV)
        if m <= 0:
            continue
        # If a fermion threshold is above the reference scale, ignore it (e.g., top if Q0=mZ).
        if m >= Q0:
            continue
        end = max(Q, m)
        if end <= 0:
            continue
        s += float(f.weight_QT3) * math.log(Q0 / end)

    return 1.0 + (float(alpha_ref) / math.pi) * s


def sin2_eff_gammaZ_1loop_from_ref(
    Q_GeV: float,
    *,
    sin2_ref: float,
    alpha_ref: float,
    Q0_GeV: float,
    fermions: list[EWChargedFermion] | None = None,
) -> float:
    """
    Convenience wrapper: sin2_eff(Q) = κ(Q) * sin2_ref under the toy κ model.
    """

    kappa = kappa_gammaZ_1loop_from_ref(
        Q_GeV,
        sin2_ref=sin2_ref,
        alpha_ref=alpha_ref,
        Q0_GeV=Q0_GeV,
        fermions=fermions,
    )
    return float(kappa) * float(sin2_ref)


