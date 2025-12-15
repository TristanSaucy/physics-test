from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class RGTarget:
    key: str
    Q_GeV: float
    note: str = ""


@dataclass(frozen=True)
class RGSuite:
    """
    A frozen RG+φ suite:
      - fit a lattice (C,m) to an anchor coupling target
      - hold C fixed (topological constant)
      - fit m for additional targets
      - compute an RG-generated scale (Lambda_QCD) implied by the lattice couplings

    This is meant to test whether adding RG / dimensional transmutation (exp(-const/alpha))
    provides a coherent bridge between discrete m-steps and physical scales.
    """

    key: str
    force: str
    anchor: RGTarget
    targets: list[RGTarget]
    n_f: int = 5
    loops: int = 2


def rg_suites() -> dict[str, RGSuite]:
    """
    Frozen suite menu.

    Note: We keep these small and interpretable; they are for *out-of-sample* style checks,
    not for tuning.
    """

    mW = 80.379
    mH = 125.0
    mt = 172.76
    return {
        "qcd-lambda-v1": RGSuite(
            key="qcd-lambda-v1",
            force="strong",
            anchor=RGTarget(
                "1/alpha_s_1loop_from_mZ(mH)",
                Q_GeV=mH,
                note="Strict strong anchor (1-loop-from-mZ at mH).",
            ),
            targets=[
                RGTarget("1/alpha_s_1loop_from_mZ(mW)", Q_GeV=mW, note="OOS: mW (nf=5)"),
                RGTarget("1/alpha_s_1loop_from_mZ(mt)", Q_GeV=mt, note="OOS-ish: mt (threshold-adjacent)"),
            ],
            n_f=5,
            loops=2,
        )
        ,
        "qcd-lambda-v2": RGSuite(
            key="qcd-lambda-v2",
            force="strong",
            anchor=RGTarget(
                "1/alpha_s(mZ)",
                Q_GeV=91.1876,
                note="Anchor at Z pole benchmark. Tests whether the lattice implies a consistent Lambda_QCD at another scale.",
            ),
            targets=[
                RGTarget("1/alpha_s_1loop_from_mZ(mH)", Q_GeV=mH, note="Cross-check at mH (nf=5)"),
            ],
            n_f=5,
            loops=2,
        ),
    }


