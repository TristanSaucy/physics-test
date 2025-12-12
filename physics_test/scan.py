from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from physics_test.model import evaluate_fit, phi


@dataclass(frozen=True)
class ScanHit:
    C: float
    m: float
    G: float
    target: float
    abs_err: float
    rel_err: float


def frange(start: float, stop: float, step: float) -> list[float]:
    """Inclusive-ish float range for scanning."""
    if step == 0:
        raise ValueError("step must be non-zero")
    vals: list[float] = []
    x = start
    # include stop with a small tolerance
    if step > 0:
        while x <= stop + 1e-12:
            vals.append(x)
            x += step
    else:
        while x >= stop - 1e-12:
            vals.append(x)
            x += step
    return vals


def scan_candidates(
    *,
    Cs: Iterable[float],
    m_values: Iterable[float],
    target_G: float,
    phi_value: float | None = None,
) -> list[ScanHit]:
    p = phi_value if phi_value is not None else phi()
    hits: list[ScanHit] = []
    for C in Cs:
        for m in m_values:
            r = evaluate_fit(target_G, C, m, phi_value=p)
            hits.append(
                ScanHit(
                    C=r.C,
                    m=r.m,
                    G=r.G,
                    target=r.target,
                    abs_err=r.abs_err,
                    rel_err=r.rel_err,
                )
            )
    hits.sort(key=lambda h: abs(h.rel_err))
    return hits


def filter_hits_by_rel_err(hits: list[ScanHit], *, max_abs_rel_err: float) -> list[ScanHit]:
    """Return only hits with |relative error| <= threshold."""
    thr = abs(float(max_abs_rel_err))
    return [h for h in hits if abs(h.rel_err) <= thr]


