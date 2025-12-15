from __future__ import annotations

import math
from dataclasses import dataclass

from physics_test.gauge_groups import (
    GaugeGroup,
    casimir_C2_fundamental,
    dynkin_index_T_fundamental,
    standard_model_gauge_groups,
)
from physics_test.gut import SM_1LOOP


@dataclass(frozen=True)
class NormalizationFamily:
    """
    A principled normalization family to transform a "physical" target into the
    lattice coordinate used by G=C/phi^m.

    This is intentionally a *small, frozen menu* of well-motivated factors:
      - group theory factors (Casimirs, Dynkin index)
      - 1-loop beta-function coefficients (SM convention)

    NOTE: A multiplicative normalization is mathematically equivalent to rescaling C,
    but keeping it explicit makes the model assumptions auditable.
    """

    key: str
    note: str


def normalization_families() -> dict[str, NormalizationFamily]:
    return {
        "none": NormalizationFamily("none", "No normalization (status quo)."),
        "C2_fund": NormalizationFamily(
            "C2_fund",
            "Multiply by quadratic Casimir C2(F) (SU(N) fundamental).",
        ),
        "inv_C2_fund": NormalizationFamily(
            "inv_C2_fund",
            "Multiply by 1/C2(F) (SU(N) fundamental). Motivated by potentials ~ C2(F)*alpha.",
        ),
        "C2_adj": NormalizationFamily(
            "C2_adj",
            "Multiply by adjoint quadratic Casimir C2(adj)=h^vee (SU(N): N).",
        ),
        "inv_C2_adj": NormalizationFamily(
            "inv_C2_adj",
            "Multiply by 1/C2(adj)=1/h^vee.",
        ),
        "T_fund": NormalizationFamily(
            "T_fund",
            "Multiply by Dynkin index T(F) (SU(N) fundamental; common physics convention: 1/2).",
        ),
        "inv_T_fund": NormalizationFamily(
            "inv_T_fund",
            "Multiply by 1/T(F).",
        ),
        "sm_beta_abs": NormalizationFamily(
            "sm_beta_abs",
            "Multiply by |b_i| from SM 1-loop beta coefficients (gut.py convention).",
        ),
        "sm_beta_abs_over_2pi": NormalizationFamily(
            "sm_beta_abs_over_2pi",
            "Multiply by |b_i|/(2π) from SM 1-loop beta coefficients.",
        ),
        "two_pi_over_sm_beta_abs": NormalizationFamily(
            "two_pi_over_sm_beta_abs",
            "Multiply by (2π)/|b_i| from SM 1-loop beta coefficients.",
        ),
        "two_pi": NormalizationFamily("two_pi", "Multiply by 2π."),
        "inv_two_pi": NormalizationFamily("inv_two_pi", "Multiply by 1/(2π)."),
        "four_pi": NormalizationFamily("four_pi", "Multiply by 4π."),
        "inv_four_pi": NormalizationFamily("inv_four_pi", "Multiply by 1/(4π)."),
    }


def _group_for_force(force: str) -> GaugeGroup | None:
    force = force.lower().strip()
    if force == "em":
        want = "U(1)"
    elif force == "weak":
        want = "SU(2)"
    elif force == "strong":
        want = "SU(3)"
    else:
        return None

    for g in standard_model_gauge_groups():
        if g.name == want:
            return g
    return None


def normalization_factor_for_force(force: str, *, family: str) -> float:
    """
    Return a multiplicative normalization factor for the given force.

    This factor is applied as:
        target_norm = target * factor
    """

    fam = family.strip()
    if fam not in normalization_families():
        raise KeyError(f"Unknown normalization family: {family!r}")

    # Gravity: keep factor=1 for now (its "coupling" is already a derived dimensionless
    # combination; any alternative normalizations should be explicitly proposed and frozen).
    if force.lower().strip() == "gravity":
        return 1.0

    g = _group_for_force(force)
    if g is None:
        return 1.0

    if fam == "none":
        return 1.0

    if fam in ("two_pi", "inv_two_pi", "four_pi", "inv_four_pi"):
        if fam == "two_pi":
            return 2.0 * math.pi
        if fam == "inv_two_pi":
            return 1.0 / (2.0 * math.pi)
        if fam == "four_pi":
            return 4.0 * math.pi
        if fam == "inv_four_pi":
            return 1.0 / (4.0 * math.pi)

    # Group invariants (SU(N) only; for U(1) we fall back to 1.0)
    if fam in ("C2_fund", "inv_C2_fund"):
        cf = casimir_C2_fundamental(g)
        if cf is None or cf == 0.0:
            return 1.0
        return float(cf) if fam == "C2_fund" else (1.0 / float(cf))

    if fam in ("T_fund", "inv_T_fund"):
        tf = dynkin_index_T_fundamental(g)
        if tf is None or tf == 0.0:
            return 1.0
        return float(tf) if fam == "T_fund" else (1.0 / float(tf))

    if fam in ("C2_adj", "inv_C2_adj"):
        # For SU(N), C2(adj)=h^vee=N; for U(1) return 1.
        if g.dual_coxeter_h is None or g.dual_coxeter_h == 0:
            return 1.0
        ca = float(g.dual_coxeter_h)
        return ca if fam == "C2_adj" else (1.0 / ca)

    # SM 1-loop beta coefficients (gut.py convention uses b_i for α_i^{-1}):
    # force mapping: strong -> b3, weak -> b2, em -> (we use 1.0 here; EM is a mixture of b1/b2 below mZ).
    if fam in ("sm_beta_abs", "sm_beta_abs_over_2pi", "two_pi_over_sm_beta_abs"):
        if g.name == "SU(3)":
            b = abs(float(SM_1LOOP.b3))
        elif g.name == "SU(2)":
            b = abs(float(SM_1LOOP.b2))
        else:
            # U(1): EM is not identical to hypercharge α1_GUT; treat as no-op here.
            b = 1.0

        if b == 0.0:
            return 1.0
        if fam == "sm_beta_abs":
            return b
        if fam == "sm_beta_abs_over_2pi":
            return b / (2.0 * math.pi)
        if fam == "two_pi_over_sm_beta_abs":
            return (2.0 * math.pi) / b

    raise RuntimeError(f"Unhandled normalization family: {fam!r}")



