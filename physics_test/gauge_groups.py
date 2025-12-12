from __future__ import annotations

from dataclasses import dataclass


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
        elif key == "base/dim":
            out[key] = float(base) / float(g.dim)
        elif key == "base/coxeter":
            if g.coxeter_h is not None:
                out[key] = float(base) / float(g.coxeter_h)
        elif key == "base/dual_coxeter":
            if g.dual_coxeter_h is not None:
                out[key] = float(base) / float(g.dual_coxeter_h)
        elif key == "base/(dim*coxeter)":
            if g.coxeter_h is not None:
                out[key] = float(base) / (float(g.dim) * float(g.coxeter_h))
        else:
            raise KeyError(f"Unknown construction key: {key}")
    return out


