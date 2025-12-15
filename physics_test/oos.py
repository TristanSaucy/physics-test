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


def oos_targets_v2() -> list[OOSTarget]:
    """
    Frozen out-of-sample target menu (v2).

    Design principle:
      - push harder on "predictive cross-checks" that are implied by the strict contract
        but not used to choose it
      - emphasize: (a) strong running cross-checks at new scales, (b) unification-style
        ratios computed using the on-shell EW definition
    """

    return [
        OOSTarget(
            "1/alpha_s_1loop_from_mZ(mW)",
            "Strong-running cross-check at mW using the same fixed 1-loop-from-mZ prescription (nf=5; no thresholds).",
        ),
        OOSTarget(
            "1/alpha_s_1loop_from_mZ(1TeV)",
            "Strong-running cross-check at 1 TeV using the same fixed 1-loop-from-mZ prescription (nf=5; no thresholds).",
        ),
        OOSTarget(
            "1/alpha_s_1loop_from_mZ(10TeV)",
            "Strong-running cross-check at 10 TeV using the same fixed 1-loop-from-mZ prescription (nf=5; no thresholds).",
        ),
        OOSTarget(
            "alpha3/alpha2(mZ,sin2_on_shell)",
            "Unification-style ratio probe at mZ, but using the on-shell EW definition for alpha2.",
        ),
        OOSTarget(
            "alpha2/alpha1_GUT(mZ,sin2_on_shell)",
            "Unification-style ratio probe at mZ (GUT normalization), using the on-shell EW definition.",
        ),
        OOSTarget(
            "1/alpha1_GUT(alpha(mZ),sin2_on_shell)",
            "Hypercharge inverse coupling derived from alpha(mZ) and on-shell sin2thetaW (GUT normalization).",
        ),
    ]


def oos_targets_v3() -> list[OOSTarget]:
    """
    Frozen out-of-sample target menu (v3).

    v3 focuses specifically on strong-running robustness: it uses a slightly more
    realistic 1-loop prescription above the top threshold by switching n_f from 5 to 6
    at mt (still no matching + no higher loops).
    """

    return [
        # High-scale cross-checks: nf changes above mt; this tests robustness of the
        # "discrete lattice" picture under a more realistic (still toy) running prescription.
        OOSTarget(
            "1/alpha_s_1loop_nf56_from_mZ(1TeV)",
            "Strong-running cross-check at 1 TeV using 1-loop running with nf=5 below mt and nf=6 above (no thresholds).",
        ),
        OOSTarget(
            "1/alpha_s_1loop_nf56_from_mZ(10TeV)",
            "Strong-running cross-check at 10 TeV using 1-loop running with nf=5 below mt and nf=6 above (no thresholds).",
        ),
        OOSTarget(
            "1/alpha_s_1loop_from_mZ(mW)",
            "Strong-running cross-check at mW (still nf=5). This tests whether the lattice prefers the Z/H anchoring vs W anchoring.",
        ),
    ]


def oos_targets_v4() -> list[OOSTarget]:
    """
    Frozen out-of-sample target menu (v4).

    v4 upgrades the strong-running cross-checks to a 2-loop beta function, while still
    avoiding any new free parameters by anchoring to alpha_s(mZ). We also keep the same
    minimal threshold story: switch n_f from 5 to 6 at mt (no matching).
    """

    return [
        OOSTarget(
            "1/alpha_s_2loop_from_mZ(mW)",
            "Strong-running cross-check at mW using 2-loop running from alpha_s(mZ) (nf=5; no thresholds).",
        ),
        OOSTarget(
            "1/alpha_s_2loop_nf56_from_mZ(1TeV)",
            "Strong-running cross-check at 1 TeV using 2-loop running with nf=5 below mt and nf=6 above (no thresholds).",
        ),
        OOSTarget(
            "1/alpha_s_2loop_nf56_from_mZ(10TeV)",
            "Strong-running cross-check at 10 TeV using 2-loop running with nf=5 below mt and nf=6 above (no thresholds).",
        ),
    ]


def oos_suites() -> dict[str, list[OOSTarget]]:
    return {
        "v1": oos_targets_v1(),
        "v2": oos_targets_v2(),
        "v3": oos_targets_v3(),
        "v4": oos_targets_v4(),
    }


def resolve_oos_targets(oos: list[OOSTarget]) -> list[TargetConstant]:
    """
    Resolve OOSTarget keys into TargetConstant entries using known_targets().
    """

    by_name = {t.name: t for t in known_targets()}
    missing = [x.key for x in oos if x.key not in by_name]
    if missing:
        raise KeyError(f"Missing OOS target keys: {missing}")
    return [by_name[x.key] for x in oos]


def predictive_force_suites() -> dict[str, tuple[dict[str, OOSTarget], dict[str, list[OOSTarget]]]]:
    """
    Predictive OOS suites.

    Unlike `oos_suites()` (which lets C vary per target), these suites are designed to
    freeze *one* best-fit C per force using a strict anchor target, then evaluate
    additional targets with C held fixed.
    """

    return {
        "v1": predictive_force_suite_v1(),
        "v2": predictive_force_suite_v2(),
        "v3": predictive_force_suite_v3(),
    }


def predictive_force_suite_v1() -> tuple[dict[str, OOSTarget], dict[str, list[OOSTarget]]]:
    """
    Predictive OOS suite v1.

    - Anchors are the frozen strict targets (one per force).
    - Predictive targets are a small, force-local menu of OOS checks.
    """

    anchors: dict[str, OOSTarget] = {
        "em": OOSTarget("1/alpha", "Strict EM anchor (vacuum, low-energy)."),
        "strong": OOSTarget("1/alpha_s_1loop_from_mZ(mH)", "Strict strong anchor (1-loop-from-mZ at mH)."),
        "weak": OOSTarget(
            "1/alpha2(alpha(mZ),sin2_on_shell)",
            "Strict weak anchor (alpha2 derived from alpha(mZ) and on-shell sin^2thetaW).",
        ),
        "gravity": OOSTarget("1/alpha_G(p)", "Strict gravity anchor (ordinary matter; proton mass)."),
    }

    targets: dict[str, list[OOSTarget]] = {
        "em": [
            OOSTarget(
                "1/alpha(mZ)",
                "OOS: effective EM inverse coupling at the Z pole (vacuum polarization) vs the low-energy anchor.",
            ),
        ],
        "strong": [
            OOSTarget(
                "1/alpha_s_1loop_from_mZ(mW)",
                "OOS: strong-running cross-check at mW (same fixed 1-loop-from-mZ prescription).",
            ),
            OOSTarget(
                "1/alpha_s_1loop_from_mZ(mt)",
                "OOS: strong-running cross-check at mt (same fixed 1-loop-from-mZ prescription).",
            ),
            OOSTarget(
                "1/alpha_s_1loop_from_mZ(1TeV)",
                "OOS: strong-running cross-check at 1 TeV (same fixed 1-loop-from-mZ prescription).",
            ),
            OOSTarget(
                "1/alpha_s_1loop_from_mZ(10TeV)",
                "OOS: strong-running cross-check at 10 TeV (same fixed 1-loop-from-mZ prescription).",
            ),
        ],
        "weak": [
            OOSTarget(
                "sin2thetaW(mZ)",
                "OOS: EW diagnostic (scheme-dependent). Tests whether the weak-lattice anchor predicts the quoted MSbar-ish value.",
            ),
            OOSTarget(
                "1/alpha2(alpha(mZ),sin2)",
                "OOS: inverse alpha2 using the MSbar-ish sin^2thetaW(mZ) value (contrast vs on-shell strict anchor).",
            ),
            OOSTarget(
                "1/alpha_w(mZ)",
                "OOS: legacy proxy 1/alpha_w(mZ) from g≈0.652 (comparison target).",
            ),
        ],
        "gravity": [
            OOSTarget(
                "1/alpha_G(e)",
                "OOS: mandatory cross-check using electron mass anchor (should imply a discrete Δm relative to proton anchor).",
            ),
            OOSTarget(
                "1/alpha_G(mP)",
                "OOS: extreme cross-check using Planck mass anchor (~1).",
            ),
        ],
    }

    return anchors, targets


def predictive_force_suite_v2() -> tuple[dict[str, OOSTarget], dict[str, list[OOSTarget]]]:
    """
    Predictive OOS suite v2.

    v2 extends v1 by adding **within-band running cross-checks for the weak sector**:
    evaluate SM 1-loop running predictions for alpha2^{-1}(Q) at additional scales,
    using the on-shell-defined alpha2(mZ) as the reference input.

    Anchors remain the strict v1 anchors.
    """

    anchors, targets = predictive_force_suite_v1()

    targets = dict(targets)
    targets["weak"] = [
        OOSTarget(
            "1/alpha2_1loop_from_mZ_on_shell(mW)",
            "OOS: SM 1-loop running cross-check for alpha2^{-1} at mW, using on-shell-defined alpha2(mZ) as reference.",
        ),
        OOSTarget(
            "1/alpha2_1loop_from_mZ_on_shell(mH)",
            "OOS: SM 1-loop running cross-check for alpha2^{-1} at mH, using on-shell-defined alpha2(mZ) as reference.",
        ),
        OOSTarget(
            "1/alpha2_1loop_from_mZ_on_shell(1TeV)",
            "OOS: SM 1-loop running cross-check for alpha2^{-1} at 1 TeV, using on-shell-defined alpha2(mZ) as reference.",
        ),
        OOSTarget(
            "1/alpha2_1loop_from_mZ_on_shell(10TeV)",
            "OOS: SM 1-loop running cross-check for alpha2^{-1} at 10 TeV, using on-shell-defined alpha2(mZ) as reference.",
        ),
    ]

    return anchors, targets


def predictive_force_suite_v3() -> tuple[dict[str, OOSTarget], dict[str, list[OOSTarget]]]:
    """
    Predictive OOS suite v3.

    v3 extends v2 with **hypercharge (GUT-normalized) within-band running** checks:
    evaluate SM 1-loop running predictions for alpha1_GUT^{-1}(Q) at additional scales,
    using the on-shell-derived alpha1_GUT(mZ) as the reference input.
    """

    anchors, targets = predictive_force_suite_v2()

    anchors = dict(anchors)
    anchors["hyper"] = OOSTarget(
        "1/alpha1_GUT(alpha(mZ),sin2_on_shell)",
        "Hypercharge (GUT-normalized) anchor derived from alpha(mZ) and on-shell sin^2thetaW.",
    )

    targets = dict(targets)
    targets["hyper"] = [
        OOSTarget(
            "1/alpha1_GUT_1loop_from_mZ_on_shell(mW)",
            "OOS: SM 1-loop running cross-check for alpha1_GUT^{-1} at mW, using on-shell-derived alpha1_GUT(mZ) as reference.",
        ),
        OOSTarget(
            "1/alpha1_GUT_1loop_from_mZ_on_shell(mH)",
            "OOS: SM 1-loop running cross-check for alpha1_GUT^{-1} at mH, using on-shell-derived alpha1_GUT(mZ) as reference.",
        ),
        OOSTarget(
            "1/alpha1_GUT_1loop_from_mZ_on_shell(1TeV)",
            "OOS: SM 1-loop running cross-check for alpha1_GUT^{-1} at 1 TeV, using on-shell-derived alpha1_GUT(mZ) as reference.",
        ),
        OOSTarget(
            "1/alpha1_GUT_1loop_from_mZ_on_shell(10TeV)",
            "OOS: SM 1-loop running cross-check for alpha1_GUT^{-1} at 10 TeV, using on-shell-derived alpha1_GUT(mZ) as reference.",
        ),
    ]

    return anchors, targets


