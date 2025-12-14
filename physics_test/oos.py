from __future__ import annotations

from dataclasses import dataclass

from physics_test.targets import TargetConstant, known_targets


@dataclass(frozen=True)
class OOSTarget:
    """
    A frozen out-of-sample target definition.

    The intent is: these are not used to choose the strict contract; they are used
    to test it afterward.
    """

    key: str
    rationale: str


def oos_targets_v1() -> list[OOSTarget]:
    """
    Frozen out-of-sample target menu (v1).

    Design principle:
      - standard, widely discussed dimensionless quantities
      - not used as part of the strict "fit" selection
      - a mix of: electroweak diagnostics, unification-style ratios, and one
        strong-running cross-check at a new scale (mt)
    """

    return [
        OOSTarget(
            "sin2thetaW(mZ)",
            "EW diagnostic (scheme-dependent but widely quoted). Not used in strict target freezing.",
        ),
        OOSTarget(
            "1/alpha2(alpha(mZ),sin2)",
            "EW inverse SU(2) coupling using the MSbar-ish sin^2thetaW(mZ) value (contrast vs on-shell strict choice).",
        ),
        OOSTarget(
            "alpha3/alpha2(mZ)",
            "Unification-style ratio probe at the Z scale.",
        ),
        OOSTarget(
            "alpha2/alpha1_GUT(mZ)",
            "Unification-style ratio probe at the Z scale (GUT normalization).",
        ),
        OOSTarget(
            "1/alpha_s_1loop_from_mZ(mt)",
            "Strong-running cross-check at a new scale (mt) using the same fixed 1-loop-from-mZ prescription (nf=5; no thresholds).",
        ),
        OOSTarget(
            "1/alpha(mZ)",
            "Effective EM inverse coupling at the Z pole (vacuum polarization). Included as an OOS check against the low-energy strict EM choice.",
        ),
    ]


def resolve_oos_targets(oos: list[OOSTarget]) -> list[TargetConstant]:
    """
    Resolve OOSTarget keys into TargetConstant entries using known_targets().
    """

    by_name = {t.name: t for t in known_targets()}
    missing = [x.key for x in oos if x.key not in by_name]
    if missing:
        raise KeyError(f"Missing OOS target keys: {missing}")
    return [by_name[x.key] for x in oos]


