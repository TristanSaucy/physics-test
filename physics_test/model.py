from __future__ import annotations

import math
from dataclasses import dataclass

from physics_test import constants


def phi() -> float:
    """Golden ratio φ."""
    return (1.0 + math.sqrt(5.0)) / 2.0


def gauge_G(C: float, m: float, *, phi_value: float | None = None) -> float:
    """
    Topological gauge (dimensionless): G = C / φ^m

    Notes:
      - If you interpret G as a *pure number* (dimensionless), then C must also be dimensionless.
    """
    p = phi_value if phi_value is not None else phi()
    # Use a log/exp form for better stability at large |m|.
    # G = C / phi^m = C * exp(-m * ln(phi))
    try:
        return C * math.exp(-float(m) * math.log(p))
    except OverflowError:
        return 0.0 if m > 0 else float("inf")


def frequency_F0(m: float, K: float, *, phi_value: float | None = None) -> float:
    """
    Frequency (Hz): F0 = φ^m * k_B * K / h

    where:
      - k_B*K is energy (J)
      - dividing by h (J*s) yields 1/s (Hz)
    """
    p = phi_value if phi_value is not None else phi()
    # F0 = phi^m * kB*K/h = exp(m*ln(phi)) * kB*K/h
    try:
        scale = math.exp(float(m) * math.log(p))
    except OverflowError:
        return float("inf") if m > 0 else 0.0
    return scale * constants.BOLTZMANN * K / constants.PLANCK


def temperature_K_from_frequency(m: float, F0: float, *, phi_value: float | None = None) -> float:
    """
    Solve the inverse of F0 = φ^m k_B K / h for temperature:

        K = F0 * h / (k_B * φ^m)
    """
    p = phi_value if phi_value is not None else phi()
    # K = F0*h/(kB*phi^m) = F0*h/kB * exp(-m*ln(phi))
    try:
        inv_scale = math.exp(-float(m) * math.log(p))
    except OverflowError:
        return 0.0 if m > 0 else float("inf")
    return F0 * constants.PLANCK / constants.BOLTZMANN * inv_scale


def coupling_invariant(C: float, K: float) -> float:
    """
    If G = C/φ^m and F0 = φ^m k_B K / h, then:

        G * F0 = C * k_B * K / h

    (independent of m).
    """
    return C * constants.BOLTZMANN * K / constants.PLANCK


@dataclass(frozen=True)
class FitResult:
    C: float
    m: float
    G: float
    target: float
    abs_err: float
    rel_err: float


def fit_C_for_target_G(target_G: float, m: float, *, phi_value: float | None = None) -> float:
    """Solve C = target_G * φ^m (so that G=C/φ^m matches target_G)."""
    p = phi_value if phi_value is not None else phi()
    return target_G * (p**m)


def evaluate_fit(target_G: float, C: float, m: float, *, phi_value: float | None = None) -> FitResult:
    """Compute G and errors vs a target."""
    G = gauge_G(C, m, phi_value=phi_value)
    abs_err = G - target_G
    rel_err = abs_err / target_G if target_G != 0 else float("nan")
    return FitResult(C=C, m=m, G=G, target=target_G, abs_err=abs_err, rel_err=rel_err)


