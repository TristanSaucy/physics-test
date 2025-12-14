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


