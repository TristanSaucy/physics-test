from __future__ import annotations

import math
from dataclasses import dataclass

from physics_test.model import phi


@dataclass(frozen=True)
class StepResult:
    """
    A C-independent "integer step" diagnostic for two targets assumed to lie on the same
    lattice (i.e., same C) so that their ratio should be a power of φ.
    """

    ratio: float
    dm_real: float
    dm_int: int
    dm_frac: float
    ratio_err_if_int: float


def step_from_targets(anchor_value: float, target_value: float, *, phi_value: float | None = None) -> StepResult:
    """
    Compute the integer-step diagnostic between two dimensionless target values.

    If both targets share the same C in:
        G = C / φ^m,
    then:
        (anchor / target) = φ^(m_target - m_anchor).

    We compute:
      - dm_real = log_phi(anchor/target)
      - dm_int  = nearest integer
      - ratio_err_if_int: relative error if you force the ratio to exactly φ^(dm_int)
    """

    if target_value == 0:
        raise ValueError("target_value must be non-zero")
    p = phi_value if phi_value is not None else phi()
    r = float(anchor_value) / float(target_value)
    dm_real = math.log(r) / math.log(p)
    dm_int = int(round(dm_real))
    dm_frac = dm_real - float(dm_int)
    # If the step were exactly integer: r_pred = φ^(dm_int)
    r_pred = p ** float(dm_int)
    ratio_err = abs((r_pred / r) - 1.0)
    return StepResult(ratio=r, dm_real=dm_real, dm_int=dm_int, dm_frac=dm_frac, ratio_err_if_int=ratio_err)


@dataclass(frozen=True)
class StepCRatioResult:
    """
    A "step signal" diagnostic that also allows a discrete C-ratio factor:

        (anchor/target) ≈ (C_a/C_b) * φ^(Δm),  with integer Δm and C_a,C_b in a fixed candidate set.
    """

    ratio: float
    C_a: float
    C_b: float
    C_ratio: float
    dm_real: float
    dm_int: int
    dm_frac: float
    ratio_err: float


def best_step_with_C_ratio(
    anchor_value: float,
    target_value: float,
    *,
    C_candidates: list[float],
    dm_min: int = -512,
    dm_max: int = 512,
    phi_value: float | None = None,
) -> StepCRatioResult:
    """
    Find the best (C_a, C_b, Δm_int) explanation of the ratio anchor/target using:

        ratio_pred = (C_a/C_b) * φ^(Δm_int)

    where Δm_int is constrained to [dm_min, dm_max].
    """

    if target_value == 0:
        raise ValueError("target_value must be non-zero")
    if dm_min > dm_max:
        raise ValueError("dm_min must be <= dm_max")
    if not C_candidates:
        raise ValueError("C_candidates must be non-empty")

    p = phi_value if phi_value is not None else phi()
    r = float(anchor_value) / float(target_value)

    best: StepCRatioResult | None = None

    for Ca in C_candidates:
        for Cb in C_candidates:
            if Cb == 0:
                continue
            c_ratio = float(Ca) / float(Cb)
            if c_ratio == 0:
                continue

            # Effective ratio that should be φ^Δm
            r_eff = r / c_ratio
            if r_eff <= 0:
                continue

            dm_real = math.log(r_eff) / math.log(p)
            dm_int = int(round(dm_real))
            if dm_int < dm_min or dm_int > dm_max:
                continue
            dm_frac = dm_real - float(dm_int)

            r_pred = c_ratio * (p ** float(dm_int))
            err = abs((r_pred / r) - 1.0)

            cand = StepCRatioResult(
                ratio=r,
                C_a=float(Ca),
                C_b=float(Cb),
                C_ratio=c_ratio,
                dm_real=dm_real,
                dm_int=dm_int,
                dm_frac=dm_frac,
                ratio_err=err,
            )
            if best is None or cand.ratio_err < best.ratio_err:
                best = cand

    if best is None:
        raise RuntimeError("No valid (C_a,C_b) candidate produced a finite step result")
    return best


