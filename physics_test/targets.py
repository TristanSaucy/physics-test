from __future__ import annotations

from dataclasses import dataclass

from physics_test import constants
from physics_test.target_registry import get_measurement
from physics_test.forces import (
    alpha_gravity,
    alpha_s_1loop,
    alpha_s_run_1loop_from_ref,
    alpha_s_run_1loop_from_ref_nf_switch,
    alpha_s_run_2loop_from_ref,
    alpha_s_run_2loop_from_ref_nf_switch,
)
from physics_test.units import mass_kg_from_GeV


@dataclass(frozen=True)
class TargetConstant:
    """A named dimensionless (pure-number) target constant."""

    name: str
    value: float
    note: str = ""
    sigma: float | None = None
    Q_GeV: float | None = None
    scheme: str = ""
    citation: str = ""


def known_targets() -> list[TargetConstant]:
    """
    A small curated list of *dimensionless* constants to test G against.

    Notes:
      - Fine-structure constant α is dimensionless; so are its inverse and 2π factors.
      - Many other famous constants (G_N, c, h, kB, etc.) have units and are not eligible here.
    """
    # Registry-backed measured inputs (value + 1σ + scheme/citation).
    m_alpha0 = get_measurement(
        "alpha0",
        default_value=constants.FINE_STRUCTURE,
        default_sigma=1.1e-12,  # CODATA 2018-ish; tiny
        default_scheme="vacuum (Thomson limit); CODATA 2018",
        default_citation="NIST CODATA 2018 (alpha)",
    )
    alpha0 = float(m_alpha0.value)
    sigma_alpha0 = m_alpha0.sigma

    m_mZ = get_measurement(
        "mZ_GeV",
        default_value=91.1876,
        default_sigma=0.0021,
        default_scheme="pole mass (benchmark)",
        default_citation="PDG RPP (mZ)",
    )
    mZ_GeV = float(m_mZ.value)
    sigma_mZ_GeV = m_mZ.sigma

    m_mW = get_measurement(
        "mW_GeV",
        default_value=80.379,
        default_sigma=0.012,
        default_scheme="pole mass (benchmark)",
        default_citation="PDG RPP (mW)",
    )
    mW_GeV = float(m_mW.value)
    sigma_mW_GeV = m_mW.sigma

    # Commonly cited effective EM coupling at the Z pole (vacuum polarization effects).
    # This is scheme-dependent and kept as an explicit target with metadata.
    m_inv_alpha_mZ = get_measurement(
        "alpha_inv_mZ",
        default_value=127.955,
        default_sigma=0.01,
        default_Q_GeV=mZ_GeV,
        default_scheme="effective EM coupling at mZ (scheme-dependent; commonly quoted)",
        default_citation="PDG RPP (alpha(mZ))",
    )
    inv_alpha_mZ = float(m_inv_alpha_mZ.value)
    sigma_inv_alpha_mZ = m_inv_alpha_mZ.sigma
    alpha_mZ = 1.0 / inv_alpha_mZ
    sigma_alpha_mZ = (sigma_inv_alpha_mZ / (inv_alpha_mZ**2)) if sigma_inv_alpha_mZ is not None else None

    # Strong coupling benchmark at mZ (commonly quoted; includes uncertainty).
    m_alpha_s_mZ = get_measurement(
        "alpha_s_mZ",
        default_value=0.1179,
        default_sigma=0.0009,
        default_Q_GeV=mZ_GeV,
        default_scheme="MSbar world average at mZ (approx)",
        default_citation="PDG RPP (alpha_s(mZ))",
    )
    alpha_s_mZ = float(m_alpha_s_mZ.value)
    sigma_alpha_s_mZ = m_alpha_s_mZ.sigma
    sigma_inv_alpha_s_mZ = (sigma_alpha_s_mZ / (alpha_s_mZ**2)) if sigma_alpha_s_mZ is not None else None

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
    sin2_thetaW_on_shell = 1.0 - (mW_GeV * mW_GeV) / (mZ_GeV * mZ_GeV)
    sigma_sin2_thetaW_on_shell: float | None = None
    if sigma_mW_GeV is not None and sigma_mZ_GeV is not None and mZ_GeV != 0:
        # First-order propagation from mW and mZ.
        # sin2 = 1 - mW^2/mZ^2
        d_dmW = -2.0 * mW_GeV / (mZ_GeV * mZ_GeV)
        d_dmZ = +2.0 * (mW_GeV * mW_GeV) / (mZ_GeV**3)
        sigma_sin2_thetaW_on_shell = (
            (d_dmW * float(sigma_mW_GeV)) ** 2 + (d_dmZ * float(sigma_mZ_GeV)) ** 2
        ) ** 0.5
    alpha2_from_alpha_mZ_on_shell = alpha_mZ / sin2_thetaW_on_shell
    sigma_alpha2_from_alpha_mZ_on_shell: float | None = None
    sigma_inv_alpha2_from_alpha_mZ_on_shell: float | None = None
    if sigma_alpha_mZ is not None or sigma_sin2_thetaW_on_shell is not None:
        rel2 = 0.0
        if sigma_alpha_mZ is not None and alpha_mZ != 0:
            rel2 += (float(sigma_alpha_mZ) / float(alpha_mZ)) ** 2
        if sigma_sin2_thetaW_on_shell is not None and sin2_thetaW_on_shell != 0:
            rel2 += (float(sigma_sin2_thetaW_on_shell) / float(sin2_thetaW_on_shell)) ** 2
        sigma_alpha2_from_alpha_mZ_on_shell = abs(float(alpha2_from_alpha_mZ_on_shell)) * (rel2**0.5)
        if alpha2_from_alpha_mZ_on_shell != 0:
            sigma_inv_alpha2_from_alpha_mZ_on_shell = float(sigma_alpha2_from_alpha_mZ_on_shell) / (
                float(alpha2_from_alpha_mZ_on_shell) ** 2
            )

    # Refined strong-coupling targets at other fixed scales via 1-loop running from mZ
    # (no free Lambda parameter).
    m_mH = get_measurement(
        "mH_GeV",
        default_value=125.0,
        default_scheme="benchmark scale",
        default_citation="project convention (mH benchmark)",
    )
    mH_GeV = float(m_mH.value)
    alpha_s_mH_1loop_from_mZ = alpha_s_run_1loop_from_ref(mH_GeV, alpha_s_Q0=alpha_s_mZ, Q0_GeV=mZ_GeV, n_f=5)

    # Out-of-sample strong running cross-checks (same fixed 1-loop prescription from mZ).
    # These are useful as "OOS targets" because they are *not* used to choose the strict targets.
    m_mt = get_measurement(
        "mt_GeV",
        default_value=172.76,
        default_sigma=0.3,
        default_scheme="pole mass (approx benchmark)",
        default_citation="PDG RPP (mt)",
    )
    mt_GeV = float(m_mt.value)
    alpha_s_mt_1loop_from_mZ = alpha_s_run_1loop_from_ref(mt_GeV, alpha_s_Q0=alpha_s_mZ, Q0_GeV=mZ_GeV, n_f=5)

    # Additional fixed-scale strong-running cross-checks (same 1-loop-from-mZ prescription, nf=5, no thresholds).
    # These are intended for the out-of-sample suite v2.
    alpha_s_mW_1loop_from_mZ = alpha_s_run_1loop_from_ref(mW_GeV, alpha_s_Q0=alpha_s_mZ, Q0_GeV=mZ_GeV, n_f=5)
    alpha_s_1TeV_1loop_from_mZ = alpha_s_run_1loop_from_ref(1_000.0, alpha_s_Q0=alpha_s_mZ, Q0_GeV=mZ_GeV, n_f=5)
    alpha_s_10TeV_1loop_from_mZ = alpha_s_run_1loop_from_ref(10_000.0, alpha_s_Q0=alpha_s_mZ, Q0_GeV=mZ_GeV, n_f=5)

    # Piecewise (nf=5 below mt, nf=6 above) variants for high-scale OOS checks.
    alpha_s_1TeV_1loop_nf56_from_mZ = alpha_s_run_1loop_from_ref_nf_switch(
        1_000.0, alpha_s_Q0=alpha_s_mZ, Q0_GeV=mZ_GeV, Q_switch_GeV=mt_GeV, n_f_below=5, n_f_above=6
    )
    alpha_s_10TeV_1loop_nf56_from_mZ = alpha_s_run_1loop_from_ref_nf_switch(
        10_000.0, alpha_s_Q0=alpha_s_mZ, Q0_GeV=mZ_GeV, Q_switch_GeV=mt_GeV, n_f_below=5, n_f_above=6
    )

    # 2-loop variants (still toy, but closer to standard QCD running). We keep the
    # same minimal threshold story: nf=5 below mt and nf=6 above mt.
    alpha_s_mW_2loop_from_mZ = alpha_s_run_2loop_from_ref(mW_GeV, alpha_s_Q0=alpha_s_mZ, Q0_GeV=mZ_GeV, n_f=5)
    alpha_s_1TeV_2loop_nf56_from_mZ = alpha_s_run_2loop_from_ref_nf_switch(
        1_000.0, alpha_s_Q0=alpha_s_mZ, Q0_GeV=mZ_GeV, Q_switch_GeV=mt_GeV, n_f_below=5, n_f_above=6
    )
    alpha_s_10TeV_2loop_nf56_from_mZ = alpha_s_run_2loop_from_ref_nf_switch(
        10_000.0, alpha_s_Q0=alpha_s_mZ, Q0_GeV=mZ_GeV, Q_switch_GeV=mt_GeV, n_f_below=5, n_f_above=6
    )

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

    # On-shell variants (use sin2_thetaW_on_shell)
    cos2_thetaW_on_shell = 1.0 - sin2_thetaW_on_shell
    alpha1_from_alpha_mZ_on_shell = alpha_mZ / cos2_thetaW_on_shell
    alpha1_gut_from_alpha_mZ_on_shell = (5.0 / 3.0) * alpha1_from_alpha_mZ_on_shell

    # Ratio-style unification probes using the on-shell EW definition
    alpha2_over_alpha1_gut_on_shell = alpha2_from_alpha_mZ_on_shell / alpha1_gut_from_alpha_mZ_on_shell
    alpha3_over_alpha2_on_shell = alpha_s_mZ / alpha2_from_alpha_mZ_on_shell

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

    inv_alpha0 = 1.0 / alpha0
    sigma_inv_alpha0 = (sigma_alpha0 / (alpha0 * alpha0)) if sigma_alpha0 is not None else None

    return [
        # EM (vacuum, low-energy)
        TargetConstant(
            "alpha",
            alpha0,
            "EM: fine-structure constant alpha (vacuum, low-energy)",
            sigma=sigma_alpha0,
            scheme=m_alpha0.scheme,
            citation=m_alpha0.citation,
        ),
        TargetConstant(
            "1/alpha",
            inv_alpha0,
            "EM: inverse fine-structure (vacuum, low-energy)",
            sigma=sigma_inv_alpha0,
            scheme=m_alpha0.scheme,
            citation=m_alpha0.citation,
        ),
        TargetConstant(
            "2*pi/alpha",
            (2.0 * constants.PI) * inv_alpha0,
            "EM: scaled variant",
            sigma=(2.0 * constants.PI) * sigma_inv_alpha0 if sigma_inv_alpha0 is not None else None,
            scheme=m_alpha0.scheme,
            citation=m_alpha0.citation,
        ),
        TargetConstant(
            "alpha/(2*pi)",
            alpha0 / (2.0 * constants.PI),
            "EM: scaled variant",
            sigma=(sigma_alpha0 / (2.0 * constants.PI)) if sigma_alpha0 is not None else None,
            scheme=m_alpha0.scheme,
            citation=m_alpha0.citation,
        ),
        # EM (effective, higher scale)
        TargetConstant(
            "alpha(mZ)",
            alpha_mZ,
            "EM: effective alpha at Z pole (scheme-dependent; commonly quoted)",
            sigma=sigma_alpha_mZ,
            Q_GeV=m_inv_alpha_mZ.Q_GeV,
            scheme=m_inv_alpha_mZ.scheme,
            citation=m_inv_alpha_mZ.citation,
        ),
        TargetConstant(
            "1/alpha(mZ)",
            inv_alpha_mZ,
            "EM: effective inverse alpha at Z pole (scheme-dependent; commonly quoted)",
            sigma=sigma_inv_alpha_mZ,
            Q_GeV=m_inv_alpha_mZ.Q_GeV,
            scheme=m_inv_alpha_mZ.scheme,
            citation=m_inv_alpha_mZ.citation,
        ),
        # Strong (benchmark)
        TargetConstant(
            "alpha_s(mZ)",
            alpha_s_mZ,
            "Strong: alpha_s at Z pole (benchmark; used as the reference input for alpha_s_1loop_from_mZ(mH))",
            sigma=sigma_alpha_s_mZ,
            Q_GeV=m_alpha_s_mZ.Q_GeV,
            scheme=m_alpha_s_mZ.scheme,
            citation=m_alpha_s_mZ.citation,
        ),
        TargetConstant(
            "1/alpha_s(mZ)",
            1.0 / alpha_s_mZ,
            "Strong: inverse alpha_s at Z pole (legacy strict benchmark; kept for comparison)",
            sigma=sigma_inv_alpha_s_mZ,
            Q_GeV=m_alpha_s_mZ.Q_GeV,
            scheme=m_alpha_s_mZ.scheme,
            citation=m_alpha_s_mZ.citation,
        ),
        # Weak (multiple alternative dimensionless couplings)
        TargetConstant("alpha_w(mZ)", alpha_w_mZ, "Weak: alpha_2=g^2/(4*pi) using g≈0.652 (legacy proxy; kept for comparison)"),
        TargetConstant("1/alpha_w(mZ)", 1.0 / alpha_w_mZ, "Weak: inverse of alpha_w(mZ) (legacy proxy; kept for comparison)"),
        TargetConstant("sin2thetaW(mZ)", sin2_thetaW_mZ, "Weak/EM: sin^2(theta_W) near mZ (legacy scheme-dependent value; comparison)"),
        TargetConstant(
            "sin2thetaW(on-shell)",
            sin2_thetaW_on_shell,
            "Weak/EM: on-shell sin^2(theta_W)=1-mW^2/mZ^2 (pole-mass definition; strict input)",
            sigma=sigma_sin2_thetaW_on_shell,
            Q_GeV=mZ_GeV,
            scheme="on-shell (pole masses): 1 - mW^2/mZ^2",
            citation=f"{m_mW.citation}; {m_mZ.citation}",
        ),
        TargetConstant("alpha_Y(mZ)", alpha_y_mZ, "Weak/EM: alpha_Y=g'^2/(4*pi) near mZ (approx; normalization depends)"),
        TargetConstant("1/alpha_Y(mZ)", 1.0 / alpha_y_mZ, "Weak/EM: inverse alpha_Y near mZ (approx)"),
        TargetConstant(
            "alpha2(alpha(mZ),sin2)",
            alpha2_from_alpha_mZ,
            "Weak: alpha_2 derived from alpha(mZ)/sin^2(thetaW) (legacy scheme-dependent; comparison)",
        ),
        TargetConstant(
            "alpha2(alpha(mZ),sin2_on_shell)",
            alpha2_from_alpha_mZ_on_shell,
            "Weak: alpha_2 derived from alpha(mZ)/sin^2(thetaW)_on-shell (pole-mass definition; strict weak definition)",
            sigma=sigma_alpha2_from_alpha_mZ_on_shell,
            Q_GeV=mZ_GeV,
            scheme="derived: alpha(mZ)/sin^2(thetaW)_on-shell",
            citation=f"{m_inv_alpha_mZ.citation}; {m_mW.citation}; {m_mZ.citation}",
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
        TargetConstant(
            "alpha1(alpha(mZ),sin2_on_shell)",
            alpha1_from_alpha_mZ_on_shell,
            "Weak/EM: alpha_1 derived from alpha(mZ)/cos^2(thetaW)_on-shell (pole-mass definition)",
        ),
        TargetConstant(
            "alpha1_GUT(alpha(mZ),sin2_on_shell)",
            alpha1_gut_from_alpha_mZ_on_shell,
            "Weak/EM: alpha_1 with 5/3 GUT normalization derived from alpha(mZ) and sin^2(thetaW)_on-shell (pole-mass definition)",
        ),
        TargetConstant("1/alpha2(alpha(mZ),sin2)", 1.0 / alpha2_from_alpha_mZ, "Inverse of alpha2(alpha(mZ),sin2)"),
        TargetConstant(
            "1/alpha2(alpha(mZ),sin2_on_shell)",
            1.0 / alpha2_from_alpha_mZ_on_shell,
            "Inverse of alpha2(alpha(mZ),sin2_on_shell) (strict weak inverse target)",
            sigma=sigma_inv_alpha2_from_alpha_mZ_on_shell,
            Q_GeV=mZ_GeV,
            scheme="derived: inverse(alpha(mZ)/sin^2(thetaW)_on-shell)",
            citation=f"{m_inv_alpha_mZ.citation}; {m_mW.citation}; {m_mZ.citation}",
        ),
        TargetConstant("1/alpha1(alpha(mZ),sin2)", 1.0 / alpha1_from_alpha_mZ, "Inverse of alpha1(alpha(mZ),sin2)"),
        TargetConstant(
            "1/alpha1_GUT(alpha(mZ),sin2)",
            1.0 / alpha1_gut_from_alpha_mZ,
            "Inverse of alpha1_GUT(alpha(mZ),sin2)",
        ),
        TargetConstant(
            "1/alpha1(alpha(mZ),sin2_on_shell)",
            1.0 / alpha1_from_alpha_mZ_on_shell,
            "Inverse of alpha1(alpha(mZ),sin2_on_shell)",
        ),
        TargetConstant(
            "1/alpha1_GUT(alpha(mZ),sin2_on_shell)",
            1.0 / alpha1_gut_from_alpha_mZ_on_shell,
            "Inverse of alpha1_GUT(alpha(mZ),sin2_on_shell)",
        ),
        # Strong (refined: 1-loop running from mZ without free Lambda; fixed scale)
        TargetConstant(
            "alpha_s_1loop_from_mZ(mH)",
            alpha_s_mH_1loop_from_mZ,
            "Strong: alpha_s at mH via 1-loop running from alpha_s(mZ) (nf=5; no thresholds) (strict strong definition)",
        ),
        TargetConstant(
            "1/alpha_s_1loop_from_mZ(mH)",
            (1.0 / alpha_s_mH_1loop_from_mZ) if alpha_s_mH_1loop_from_mZ != 0 else float("inf"),
            "Strong: inverse of alpha_s_1loop_from_mZ(mH) (strict strong inverse target)",
            sigma=sigma_inv_alpha_s_mZ,
            Q_GeV=mH_GeV,
            scheme="derived: 1-loop from alpha_s(mZ) (nf=5; no thresholds)",
            citation=m_alpha_s_mZ.citation,
        ),
        TargetConstant(
            "alpha_s_1loop_from_mZ(mW)",
            alpha_s_mW_1loop_from_mZ,
            "Strong (OOS): alpha_s at mW via 1-loop running from alpha_s(mZ) (nf=5; no thresholds)",
        ),
        TargetConstant(
            "1/alpha_s_1loop_from_mZ(mW)",
            (1.0 / alpha_s_mW_1loop_from_mZ) if alpha_s_mW_1loop_from_mZ != 0 else float("inf"),
            "Strong (OOS): inverse of alpha_s_1loop_from_mZ(mW)",
            sigma=sigma_inv_alpha_s_mZ,
            Q_GeV=mW_GeV,
            scheme="derived: 1-loop from alpha_s(mZ) (nf=5; no thresholds)",
            citation=m_alpha_s_mZ.citation,
        ),
        TargetConstant(
            "alpha_s_1loop_from_mZ(mt)",
            alpha_s_mt_1loop_from_mZ,
            "Strong (OOS): alpha_s at mt via 1-loop running from alpha_s(mZ) (nf=5; no thresholds)",
        ),
        TargetConstant(
            "1/alpha_s_1loop_from_mZ(mt)",
            (1.0 / alpha_s_mt_1loop_from_mZ) if alpha_s_mt_1loop_from_mZ != 0 else float("inf"),
            "Strong (OOS): inverse of alpha_s_1loop_from_mZ(mt)",
            sigma=sigma_inv_alpha_s_mZ,
            Q_GeV=mt_GeV,
            scheme="derived: 1-loop from alpha_s(mZ) (nf=5; no thresholds)",
            citation=m_alpha_s_mZ.citation,
        ),
        TargetConstant(
            "alpha_s_1loop_from_mZ(1TeV)",
            alpha_s_1TeV_1loop_from_mZ,
            "Strong (OOS): alpha_s at 1 TeV via 1-loop running from alpha_s(mZ) (nf=5; no thresholds)",
        ),
        TargetConstant(
            "1/alpha_s_1loop_from_mZ(1TeV)",
            (1.0 / alpha_s_1TeV_1loop_from_mZ) if alpha_s_1TeV_1loop_from_mZ != 0 else float("inf"),
            "Strong (OOS): inverse of alpha_s_1loop_from_mZ(1TeV)",
            sigma=sigma_inv_alpha_s_mZ,
            Q_GeV=1_000.0,
            scheme="derived: 1-loop from alpha_s(mZ) (nf=5; no thresholds)",
            citation=m_alpha_s_mZ.citation,
        ),
        TargetConstant(
            "alpha_s_1loop_from_mZ(10TeV)",
            alpha_s_10TeV_1loop_from_mZ,
            "Strong (OOS): alpha_s at 10 TeV via 1-loop running from alpha_s(mZ) (nf=5; no thresholds)",
        ),
        TargetConstant(
            "1/alpha_s_1loop_from_mZ(10TeV)",
            (1.0 / alpha_s_10TeV_1loop_from_mZ) if alpha_s_10TeV_1loop_from_mZ != 0 else float("inf"),
            "Strong (OOS): inverse of alpha_s_1loop_from_mZ(10TeV)",
            sigma=sigma_inv_alpha_s_mZ,
            Q_GeV=10_000.0,
            scheme="derived: 1-loop from alpha_s(mZ) (nf=5; no thresholds)",
            citation=m_alpha_s_mZ.citation,
        ),
        TargetConstant(
            "alpha_s_1loop_nf56_from_mZ(1TeV)",
            alpha_s_1TeV_1loop_nf56_from_mZ,
            "Strong (OOS): alpha_s at 1 TeV via 1-loop running with nf=5 below mt and nf=6 above (no thresholds)",
        ),
        TargetConstant(
            "1/alpha_s_1loop_nf56_from_mZ(1TeV)",
            (1.0 / alpha_s_1TeV_1loop_nf56_from_mZ) if alpha_s_1TeV_1loop_nf56_from_mZ != 0 else float("inf"),
            "Strong (OOS): inverse of alpha_s_1loop_nf56_from_mZ(1TeV)",
            sigma=sigma_inv_alpha_s_mZ,
            Q_GeV=1_000.0,
            scheme="derived: 1-loop from alpha_s(mZ) with nf=5 below mt and nf=6 above",
            citation=m_alpha_s_mZ.citation,
        ),
        TargetConstant(
            "alpha_s_1loop_nf56_from_mZ(10TeV)",
            alpha_s_10TeV_1loop_nf56_from_mZ,
            "Strong (OOS): alpha_s at 10 TeV via 1-loop running with nf=5 below mt and nf=6 above (no thresholds)",
        ),
        TargetConstant(
            "1/alpha_s_1loop_nf56_from_mZ(10TeV)",
            (1.0 / alpha_s_10TeV_1loop_nf56_from_mZ) if alpha_s_10TeV_1loop_nf56_from_mZ != 0 else float("inf"),
            "Strong (OOS): inverse of alpha_s_1loop_nf56_from_mZ(10TeV)",
            sigma=sigma_inv_alpha_s_mZ,
            Q_GeV=10_000.0,
            scheme="derived: 1-loop from alpha_s(mZ) with nf=5 below mt and nf=6 above",
            citation=m_alpha_s_mZ.citation,
        ),
        TargetConstant(
            "alpha_s_2loop_from_mZ(mW)",
            alpha_s_mW_2loop_from_mZ,
            "Strong (OOS): alpha_s at mW via 2-loop running from alpha_s(mZ) (nf=5; no thresholds)",
        ),
        TargetConstant(
            "1/alpha_s_2loop_from_mZ(mW)",
            (1.0 / alpha_s_mW_2loop_from_mZ) if alpha_s_mW_2loop_from_mZ != 0 else float("inf"),
            "Strong (OOS): inverse of alpha_s_2loop_from_mZ(mW)",
            sigma=sigma_inv_alpha_s_mZ,
            Q_GeV=mW_GeV,
            scheme="derived: 2-loop from alpha_s(mZ) (nf=5; no thresholds)",
            citation=m_alpha_s_mZ.citation,
        ),
        TargetConstant(
            "alpha_s_2loop_nf56_from_mZ(1TeV)",
            alpha_s_1TeV_2loop_nf56_from_mZ,
            "Strong (OOS): alpha_s at 1 TeV via 2-loop running with nf=5 below mt and nf=6 above (no thresholds)",
        ),
        TargetConstant(
            "1/alpha_s_2loop_nf56_from_mZ(1TeV)",
            (1.0 / alpha_s_1TeV_2loop_nf56_from_mZ) if alpha_s_1TeV_2loop_nf56_from_mZ != 0 else float("inf"),
            "Strong (OOS): inverse of alpha_s_2loop_nf56_from_mZ(1TeV)",
            sigma=sigma_inv_alpha_s_mZ,
            Q_GeV=1_000.0,
            scheme="derived: 2-loop from alpha_s(mZ) with nf=5 below mt and nf=6 above",
            citation=m_alpha_s_mZ.citation,
        ),
        TargetConstant(
            "alpha_s_2loop_nf56_from_mZ(10TeV)",
            alpha_s_10TeV_2loop_nf56_from_mZ,
            "Strong (OOS): alpha_s at 10 TeV via 2-loop running with nf=5 below mt and nf=6 above (no thresholds)",
        ),
        TargetConstant(
            "1/alpha_s_2loop_nf56_from_mZ(10TeV)",
            (1.0 / alpha_s_10TeV_2loop_nf56_from_mZ) if alpha_s_10TeV_2loop_nf56_from_mZ != 0 else float("inf"),
            "Strong (OOS): inverse of alpha_s_2loop_nf56_from_mZ(10TeV)",
            sigma=sigma_inv_alpha_s_mZ,
            Q_GeV=10_000.0,
            scheme="derived: 2-loop from alpha_s(mZ) with nf=5 below mt and nf=6 above",
            citation=m_alpha_s_mZ.citation,
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
        TargetConstant(
            "alpha2/alpha1_GUT(mZ,sin2_on_shell)",
            alpha2_over_alpha1_gut_on_shell,
            "Unification probe: alpha2 / alpha1_GUT using sin2thetaW(on-shell) (pole-mass definition)",
        ),
        TargetConstant(
            "alpha3/alpha2(mZ,sin2_on_shell)",
            alpha3_over_alpha2_on_shell,
            "Unification probe: alpha3(alpha_s) / alpha2 using sin2thetaW(on-shell) (pole-mass definition)",
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


