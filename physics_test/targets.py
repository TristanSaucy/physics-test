from __future__ import annotations

from dataclasses import dataclass

from physics_test import constants
from physics_test.forces import alpha_gravity, alpha_s_1loop, alpha_s_run_1loop_from_ref
from physics_test.units import mass_kg_from_GeV


@dataclass(frozen=True)
class TargetConstant:
    """A named dimensionless (pure-number) target constant."""

    name: str
    value: float
    note: str = ""


def known_targets() -> list[TargetConstant]:
    """
    A small curated list of *dimensionless* constants to test G against.

    Notes:
      - Fine-structure constant α is dimensionless; so are its inverse and 2π factors.
      - Many other famous constants (G_N, c, h, kB, etc.) have units and are not eligible here.
    """
    alpha0 = constants.FINE_STRUCTURE

    # Commonly cited effective EM coupling at the Z pole (vacuum polarization effects).
    # This is approximate and included for "all options" scanning.
    inv_alpha_mZ = 127.955  # ~ 1/alpha(mZ)
    alpha_mZ = 1.0 / inv_alpha_mZ

    # Strong coupling benchmark at mZ (commonly quoted); included as an exploratory target.
    alpha_s_mZ = 0.1179

    # Weak coupling (exploratory): alpha_W = g^2/(4*pi) at ~mZ scale.
    # Using a commonly cited SU(2)_L coupling g ≈ 0.652 at the electroweak scale.
    # This is approximate and scale-dependent; included for "all options" scanning.
    g_ew = 0.652
    alpha_w_mZ = (g_ew**2) / (4.0 * constants.PI)

    # Electroweak mixing angle (common on-shell-ish reference near mZ; approximate).
    # sin^2(theta_W) is dimensionless and often quoted; it is scheme/scale dependent.
    sin2_thetaW_mZ = 0.23122

    # On-shell weak mixing angle from pole masses:
    #   sin^2(theta_W)_OS = 1 - (mW^2 / mZ^2)
    # This is a different (standard) scheme than the MSbar-ish 0.23122 above.
    mW_GeV = 80.379
    mZ_GeV = 91.1876
    sin2_thetaW_on_shell = 1.0 - (mW_GeV * mW_GeV) / (mZ_GeV * mZ_GeV)
    alpha2_from_alpha_mZ_on_shell = alpha_mZ / sin2_thetaW_on_shell

    # Refined strong-coupling targets at other fixed scales via 1-loop running from mZ
    # (no free Lambda parameter).
    mH_GeV = 125.0
    alpha_s_mH_1loop_from_mZ = alpha_s_run_1loop_from_ref(mH_GeV, alpha_s_Q0=alpha_s_mZ, Q0_GeV=mZ_GeV, n_f=5)

    # Hypercharge-like coupling g' (approx near mZ; note normalization conventions vary).
    gprime_ew = 0.357
    alpha_y_mZ = (gprime_ew**2) / (4.0 * constants.PI)

    # Also derive alpha_2 and alpha_1 from alpha(mZ) and sin^2(thetaW) using:
    #   alpha_2 = alpha / sin^2(thetaW)
    #   alpha_1 = alpha / cos^2(thetaW)   (WITHOUT the GUT 5/3 normalization factor)
    cos2_thetaW_mZ = 1.0 - sin2_thetaW_mZ
    alpha2_from_alpha_mZ = alpha_mZ / sin2_thetaW_mZ
    alpha1_from_alpha_mZ = alpha_mZ / cos2_thetaW_mZ
    alpha1_gut_from_alpha_mZ = (5.0 / 3.0) * alpha1_from_alpha_mZ

    # Gravity: choose mass scale -> dimensionless coupling changes with scale.
    alphaG_e = alpha_gravity(constants.MASS_ELECTRON)
    alphaG_p = alpha_gravity(constants.MASS_PROTON)
    alphaG_P = alpha_gravity(constants.MASS_PLANCK)  # ~ 1 by construction

    # Additional gravity "types" via different mass scales (including TeV-scale for "quantum gravity" exploration).
    gravity_mass_scales_GeV: list[tuple[str, float]] = [
        ("W", 80.379),
        ("Z", 91.1876),
        ("H", 125.0),
        ("1TeV", 1_000.0),
        ("10TeV", 10_000.0),
        ("40TeV", 40_000.0),
        ("100TeV", 100_000.0),
        # Frozen GW-band / primordial gravity-type anchors (derived from strict inverse-gravity sweeps under CMB K).
        ("GW_CMB", 2.93012e4),
        ("GW_PTA", 1.58009e9),
        ("GW_LISA", 2.15524e12),
        ("GW_LIGO", 5.41086e13),
    ]
    extra_gravity_targets: list[TargetConstant] = []
    for label, gev in gravity_mass_scales_GeV:
        mkg = mass_kg_from_GeV(gev)
        a = alpha_gravity(mkg)
        if label.startswith("GW_"):
            band = label.split("_", 1)[1]
            note_a = f"Gravity (frozen GW type): alpha_G using GW-{band} mass anchor {gev} GeV"
            note_inv = f"Gravity (frozen GW type): inverse alpha_G using GW-{band} mass anchor {gev} GeV"
        else:
            note_a = f"Gravity: alpha_G using mass scale {gev} GeV"
            note_inv = f"Gravity: inverse alpha_G using mass scale {gev} GeV"
        extra_gravity_targets.append(
            TargetConstant(f"alpha_G({label})", a, note_a)
        )
        extra_gravity_targets.append(
            TargetConstant(f"1/alpha_G({label})", 1.0 / a, note_inv)
        )

    # Optional: toy running alpha_s(Q) values for a few Q points, just to scan patterns.
    alpha_s_10 = alpha_s_1loop(10.0)
    alpha_s_91 = alpha_s_1loop(91.1876)  # ~ mZ in GeV, toy model
    alpha_s_1 = alpha_s_1loop(1.0)
    alpha_s_2 = alpha_s_1loop(2.0)
    alpha_s_5 = alpha_s_1loop(5.0)
    alpha_s_1000 = alpha_s_1loop(1000.0)

    return [
        # EM (vacuum, low-energy)
        TargetConstant("alpha", alpha0, "EM: fine-structure constant alpha (vacuum, low-energy)"),
        TargetConstant("1/alpha", 1.0 / alpha0, "EM: inverse fine-structure (vacuum, low-energy)"),
        TargetConstant("2*pi/alpha", (2.0 * constants.PI) / alpha0, "EM: scaled variant"),
        TargetConstant("alpha/(2*pi)", alpha0 / (2.0 * constants.PI), "EM: scaled variant"),
        # EM (effective, higher scale)
        TargetConstant("alpha(mZ)", alpha_mZ, "EM: effective alpha at Z pole (approx)"),
        TargetConstant("1/alpha(mZ)", inv_alpha_mZ, "EM: effective inverse alpha at Z pole (approx)"),
        # Strong (benchmark)
        TargetConstant("alpha_s(mZ)", alpha_s_mZ, "Strong: alpha_s at Z pole (benchmark ~0.118)"),
        TargetConstant("1/alpha_s(mZ)", 1.0 / alpha_s_mZ, "Strong: inverse alpha_s at Z pole (benchmark)"),
        # Weak (multiple alternative dimensionless couplings)
        TargetConstant("alpha_w(mZ)", alpha_w_mZ, "Weak: alpha_2=g^2/(4*pi) near EW scale (approx)"),
        TargetConstant("1/alpha_w(mZ)", 1.0 / alpha_w_mZ, "Weak: inverse alpha_2 near EW scale (approx)"),
        TargetConstant("sin2thetaW(mZ)", sin2_thetaW_mZ, "Weak/EM: sin^2(theta_W) near mZ (approx)"),
        TargetConstant(
            "sin2thetaW(on-shell)",
            sin2_thetaW_on_shell,
            "Weak/EM: on-shell sin^2(theta_W)=1-mW^2/mZ^2 (pole-mass definition)",
        ),
        TargetConstant("alpha_Y(mZ)", alpha_y_mZ, "Weak/EM: alpha_Y=g'^2/(4*pi) near mZ (approx; normalization depends)"),
        TargetConstant("1/alpha_Y(mZ)", 1.0 / alpha_y_mZ, "Weak/EM: inverse alpha_Y near mZ (approx)"),
        TargetConstant(
            "alpha2(alpha(mZ),sin2)",
            alpha2_from_alpha_mZ,
            "Weak: alpha_2 derived from alpha(mZ)/sin^2(thetaW) (approx)",
        ),
        TargetConstant(
            "alpha2(alpha(mZ),sin2_on_shell)",
            alpha2_from_alpha_mZ_on_shell,
            "Weak: alpha_2 derived from alpha(mZ)/sin^2(thetaW)_on-shell (pole-mass definition)",
        ),
        TargetConstant(
            "alpha1(alpha(mZ),sin2)",
            alpha1_from_alpha_mZ,
            "Weak/EM: alpha_1 derived from alpha(mZ)/cos^2(thetaW) (no 5/3 GUT factor) (approx)",
        ),
        TargetConstant(
            "alpha1_GUT(alpha(mZ),sin2)",
            alpha1_gut_from_alpha_mZ,
            "Weak/EM: alpha_1 with 5/3 GUT normalization derived from alpha(mZ) and sin^2(thetaW) (approx)",
        ),
        TargetConstant("1/alpha2(alpha(mZ),sin2)", 1.0 / alpha2_from_alpha_mZ, "Inverse of alpha2(alpha(mZ),sin2)"),
        TargetConstant(
            "1/alpha2(alpha(mZ),sin2_on_shell)",
            1.0 / alpha2_from_alpha_mZ_on_shell,
            "Inverse of alpha2(alpha(mZ),sin2_on_shell)",
        ),
        TargetConstant("1/alpha1(alpha(mZ),sin2)", 1.0 / alpha1_from_alpha_mZ, "Inverse of alpha1(alpha(mZ),sin2)"),
        TargetConstant(
            "1/alpha1_GUT(alpha(mZ),sin2)",
            1.0 / alpha1_gut_from_alpha_mZ,
            "Inverse of alpha1_GUT(alpha(mZ),sin2)",
        ),
        # Strong (refined: 1-loop running from mZ without free Lambda; fixed scale)
        TargetConstant(
            "alpha_s_1loop_from_mZ(mH)",
            alpha_s_mH_1loop_from_mZ,
            "Strong: alpha_s at mH via 1-loop running from alpha_s(mZ) (nf=5; no thresholds)",
        ),
        TargetConstant(
            "1/alpha_s_1loop_from_mZ(mH)",
            (1.0 / alpha_s_mH_1loop_from_mZ) if alpha_s_mH_1loop_from_mZ != 0 else float("inf"),
            "Strong: inverse of alpha_s_1loop_from_mZ(mH)",
        ),
        # Unification probes (dimensionless ratios; still scale-dependent in reality)
        TargetConstant(
            "alpha2/alpha1_GUT(mZ)",
            alpha2_from_alpha_mZ / alpha1_gut_from_alpha_mZ,
            "Unification probe: alpha2 / alpha1_GUT at mZ (approx)",
        ),
        TargetConstant(
            "alpha3/alpha2(mZ)",
            alpha_s_mZ / alpha2_from_alpha_mZ,
            "Unification probe: alpha3(alpha_s) / alpha2 at mZ (approx)",
        ),
        # Strong (toy running examples)
        TargetConstant("alpha_s_1loop(1GeV)", alpha_s_1, "Strong: toy 1-loop alpha_s at 1 GeV"),
        TargetConstant("alpha_s_1loop(2GeV)", alpha_s_2, "Strong: toy 1-loop alpha_s at 2 GeV"),
        TargetConstant("alpha_s_1loop(5GeV)", alpha_s_5, "Strong: toy 1-loop alpha_s at 5 GeV"),
        TargetConstant("alpha_s_1loop(10GeV)", alpha_s_10, "Strong: toy 1-loop alpha_s at 10 GeV"),
        TargetConstant("alpha_s_1loop(mZ)", alpha_s_91, "Strong: toy 1-loop alpha_s at mZ"),
        TargetConstant("alpha_s_1loop(1TeV)", alpha_s_1000, "Strong: toy 1-loop alpha_s at 1 TeV"),
        TargetConstant(
            "1/alpha_s_1loop(1GeV)",
            (1.0 / alpha_s_1) if alpha_s_1 != 0 else float("inf"),
            "Strong: inverse toy 1-loop alpha_s at 1 GeV",
        ),
        TargetConstant(
            "1/alpha_s_1loop(2GeV)",
            (1.0 / alpha_s_2) if alpha_s_2 != 0 else float("inf"),
            "Strong: inverse toy 1-loop alpha_s at 2 GeV",
        ),
        TargetConstant(
            "1/alpha_s_1loop(5GeV)",
            (1.0 / alpha_s_5) if alpha_s_5 != 0 else float("inf"),
            "Strong: inverse toy 1-loop alpha_s at 5 GeV",
        ),
        TargetConstant(
            "1/alpha_s_1loop(10GeV)",
            (1.0 / alpha_s_10) if alpha_s_10 != 0 else float("inf"),
            "Strong: inverse toy 1-loop alpha_s at 10 GeV",
        ),
        TargetConstant(
            "1/alpha_s_1loop(mZ)",
            (1.0 / alpha_s_91) if alpha_s_91 != 0 else float("inf"),
            "Strong: inverse toy 1-loop alpha_s at mZ",
        ),
        TargetConstant(
            "1/alpha_s_1loop(1TeV)",
            (1.0 / alpha_s_1000) if alpha_s_1000 != 0 else float("inf"),
            "Strong: inverse toy 1-loop alpha_s at 1 TeV",
        ),
        # Gravity (scale depends on chosen mass)
        TargetConstant("alpha_G(e)", alphaG_e, "Gravity: alpha_G using electron mass scale"),
        TargetConstant("alpha_G(p)", alphaG_p, "Gravity: alpha_G using proton mass scale"),
        TargetConstant("alpha_G(mP)", alphaG_P, "Gravity: alpha_G using Planck mass scale (~1)"),
        TargetConstant("1/alpha_G(e)", 1.0 / alphaG_e, "Gravity: inverse alpha_G using electron mass scale"),
        TargetConstant("1/alpha_G(p)", 1.0 / alphaG_p, "Gravity: inverse alpha_G using proton mass scale"),
        TargetConstant(
            "1/alpha_G(mP)",
            (1.0 / alphaG_P) if alphaG_P != 0 else float("inf"),
            "Gravity: inverse alpha_G using Planck mass scale (~1)",
        ),
        *extra_gravity_targets,
    ]


