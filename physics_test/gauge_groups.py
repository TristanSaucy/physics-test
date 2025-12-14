from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class GaugeGroup:
    """
    Minimal gauge-group metadata for exploratory "topological constant" constructions.

    This is not a claim of correctness—just a structured way to generate non-arbitrary
    candidate numbers from group invariants.
    """

    name: str
    rank: int
    dim: int
    coxeter_h: int | None = None
    dual_coxeter_h: int | None = None


def standard_model_gauge_groups() -> list[GaugeGroup]:
    # U(1) has rank=1, dim=1 in the Lie algebra sense.
    # SU(N): rank=N-1, dim=N^2-1, Coxeter=dual_coxeter=N.
    return [
        GaugeGroup("U(1)", rank=1, dim=1, coxeter_h=None, dual_coxeter_h=None),
        GaugeGroup("SU(2)", rank=1, dim=3, coxeter_h=2, dual_coxeter_h=2),
        GaugeGroup("SU(3)", rank=2, dim=8, coxeter_h=3, dual_coxeter_h=3),
    ]


def _su_n(name: str) -> int | None:
    """
    Parse SU(N) from a name like "SU(3)".
    Returns None if not an SU(N) string.
    """

    if name.startswith("SU(") and name.endswith(")"):
        inner = name[3:-1]
        if inner.isdigit():
            return int(inner)
    return None


def roots_count(g: GaugeGroup) -> int:
    """
    Number of roots in the root system for a compact semisimple Lie algebra:
      dim(g) = rank(g) + |Phi|

    For an abelian factor like U(1), this yields 0.
    """

    n = int(g.dim) - int(g.rank)
    return max(0, n)


def positive_roots_count(g: GaugeGroup) -> int | None:
    """
    Number of positive roots |Phi^+| = |Phi|/2 (for reduced root systems).
    Returns None if |Phi| is odd.
    """

    r = roots_count(g)
    return (r // 2) if (r % 2 == 0) else None


def weyl_group_order(g: GaugeGroup) -> int | None:
    """
    Weyl group order |W|.

    Implemented minimally for:
      - U(1): 1
      - SU(N): N!
    """

    if g.name == "U(1)":
        return 1
    n = _su_n(g.name)
    if n is not None:
        return math.factorial(n)
    return None


def center_order(g: GaugeGroup) -> int | None:
    """
    Order of the center |Z(G)| for simply-connected SU(N): |Z| = N.

    For U(1), the center is continuous (not finite), so we return None.
    """

    n = _su_n(g.name)
    if n is not None:
        return n
    return None


def casimir_C2_fundamental(g: GaugeGroup) -> float | None:
    """
    Quadratic Casimir C2(F) for the fundamental representation of SU(N):
      C_F = (N^2 - 1) / (2N)

    This depends on a representation choice and normalization convention.
    We implement the common physics convention where long roots have length^2=2
    and for SU(N) the fundamental Dynkin index is T(F)=1/2.
    """

    n = _su_n(g.name)
    if n is None:
        return None
    return (float(n * n) - 1.0) / (2.0 * float(n))


def dynkin_index_T_fundamental(g: GaugeGroup) -> float | None:
    """
    Dynkin index T(F) for the fundamental representation.

    For SU(N) in the common physics normalization: T(F)=1/2.
    """

    n = _su_n(g.name)
    if n is None:
        return None
    return 0.5


def candidate_Cs_from_group(
    g: GaugeGroup,
    *,
    base: float = 360.0,
    include: tuple[str, ...] = ("base", "base/dim", "base/coxeter", "base/dual_coxeter", "base/(dim*coxeter)"),
) -> dict[str, float]:
    """
    Generate candidate C values from a base number (default 360) and group invariants.

    Returned keys describe the construction.
    """
    out: dict[str, float] = {}
    for key in include:
        if key == "base":
            out[key] = float(base)
        elif key == "base/rank":
            if g.rank != 0:
                out[key] = float(base) / float(g.rank)
        elif key == "base/dim":
            out[key] = float(base) / float(g.dim)
        elif key == "base/coxeter":
            if g.coxeter_h is not None:
                out[key] = float(base) / float(g.coxeter_h)
        elif key == "base/dual_coxeter":
            if g.dual_coxeter_h is not None:
                out[key] = float(base) / float(g.dual_coxeter_h)
        elif key == "base/C2_adj":
            # Under standard normalization for simple Lie algebras: C2(adj) = h^\vee
            if g.dual_coxeter_h is not None:
                out[key] = float(base) / float(g.dual_coxeter_h)
        elif key == "base/roots":
            rc = roots_count(g)
            if rc != 0:
                out[key] = float(base) / float(rc)
        elif key == "base/positive_roots":
            pr = positive_roots_count(g)
            if pr is not None and pr != 0:
                out[key] = float(base) / float(pr)
        elif key == "base/weyl_order":
            w = weyl_group_order(g)
            if w is not None and w != 0:
                out[key] = float(base) / float(w)
        elif key == "base/center_order":
            z = center_order(g)
            if z is not None and z != 0:
                out[key] = float(base) / float(z)
        elif key == "base/C2_fund":
            cf = casimir_C2_fundamental(g)
            if cf is not None and cf != 0.0:
                out[key] = float(base) / float(cf)
        elif key == "base/T_fund":
            tf = dynkin_index_T_fundamental(g)
            if tf is not None and tf != 0.0:
                out[key] = float(base) / float(tf)
        elif key == "base/(dim*coxeter)":
            if g.coxeter_h is not None:
                out[key] = float(base) / (float(g.dim) * float(g.coxeter_h))
        else:
            raise KeyError(f"Unknown construction key: {key}")
    return out


