from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple


PHI = (1.0 + math.sqrt(5.0)) / 2.0


# ---------------------------------------------------------------------------
# 1-loop beta coefficients
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BetaCoefficients:
    """
    One-loop beta function coefficients b_i for gauge couplings in:

      d(α_i^{-1}) / d ln μ = - b_i / (2π)

    Conventions vary; this module is only for exploratory convergence testing.
    """

    b1: float  # U(1)_Y (GUT-normalized)
    b2: float  # SU(2)
    b3: float  # SU(3)
    name: str = "custom"


SM_1LOOP = BetaCoefficients(b1=41.0 / 10.0, b2=-19.0 / 6.0, b3=-7.0, name="SM (1-loop)")
MSSM_1LOOP = BetaCoefficients(b1=33.0 / 5.0, b2=1.0, b3=-3.0, name="MSSM (1-loop)")


# ---------------------------------------------------------------------------
# 2-loop beta coefficients  (Machacek & Vaughn 1983-84; Jones 1982)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Beta2Loop:
    """
    Full 2-loop beta coefficients for the three SM gauge couplings.

    Inverse-coupling convention:
      d(α_i⁻¹)/dt = -b_i/(2π) - Σ_j B_ij α_j/(8π²)

    Equivalently for the coupling:
      dα_i/dt = (b_i/(2π)) α_i²  +  Σ_j (B_ij/(8π²)) α_i² α_j

    where t = ln μ and i,j ∈ {1,2,3} for U(1)_Y, SU(2)_L, SU(3)_c.
    """

    b: Tuple[float, float, float]
    B: Tuple[Tuple[float, float, float],
             Tuple[float, float, float],
             Tuple[float, float, float]]
    name: str = "custom"


SM_2LOOP = Beta2Loop(
    b=(41.0 / 10.0, -19.0 / 6.0, -7.0),
    B=(
        (199.0 / 50.0, 27.0 / 10.0, 44.0 / 5.0),
        (9.0 / 10.0, 35.0 / 6.0, 12.0),
        (11.0 / 10.0, 9.0 / 2.0, -26.0),
    ),
    name="SM (2-loop)",
)

MSSM_2LOOP = Beta2Loop(
    b=(33.0 / 5.0, 1.0, -3.0),
    B=(
        (199.0 / 25.0, 27.0 / 5.0, 88.0 / 5.0),
        (9.0 / 5.0, 25.0, 24.0),
        (11.0 / 5.0, 9.0, 14.0),
    ),
    name="MSSM (2-loop)",
)


# ---------------------------------------------------------------------------
# GUT normalization factors for α₁ = k₁ · α_Y
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GUTNormalization:
    k1: float
    name: str
    note: str = ""


GUT_NORMALIZATIONS: List[GUTNormalization] = [
    GUTNormalization(5.0 / 3.0, "SU(5)",
                     "Standard GUT normalization; also applies to SO(10) standard embedding"),
    GUTNormalization(5.0 / 3.0, "SO(10)",
                     "Same k₁ as SU(5) for standard SU(5)⊂SO(10) decomposition"),
    GUTNormalization(5.0 / 3.0, "E6 (standard)",
                     "Standard E6 embedding inherits SU(5) normalization"),
    GUTNormalization(4.0 / 3.0, "SU(3)³ trinification",
                     "Trinification: SU(3)_c×SU(3)_L×SU(3)_R with democratic coupling"),
    GUTNormalization(1.0, "SU(4)_PS (Pati-Salam, direct)",
                     "Pati-Salam SU(4)×SU(2)×SU(2); k₁=1 treats α_Y directly as GUT coupling"),
    GUTNormalization(2.0, "Flipped SU(5)×U(1)_X (k=2)",
                     "Flipped SU(5) with alternative U(1) embedding"),
    GUTNormalization(8.0 / 3.0, "SU(6) minimal",
                     "SU(6) with 6-plet fundamental representation"),
    GUTNormalization(5.0 / 4.0, "E6 (U(1)_ψ mixed)",
                     "E6 with U(1)_ψ mixing rotated hypercharge generator"),
]


# ---------------------------------------------------------------------------
# 1-loop helpers
# ---------------------------------------------------------------------------

def run_alpha_inv(alpha_inv_mu0: float, mu0: float, mu: float, b: float) -> float:
    """1-loop running: α^{-1}(μ) = α^{-1}(μ0) - (b/(2π)) ln(μ/μ0)."""
    return float(alpha_inv_mu0) - (float(b) / (2.0 * math.pi)) * math.log(float(mu) / float(mu0))


def converge_score(alpha1_inv: float, alpha2_inv: float, alpha3_inv: float) -> float:
    """Simple convergence score: max pairwise difference in α^{-1}."""
    return max(abs(alpha1_inv - alpha2_inv), abs(alpha1_inv - alpha3_inv), abs(alpha2_inv - alpha3_inv))


def find_best_convergence(
    *,
    mu0: float,
    alpha1_inv_mu0: float,
    alpha2_inv_mu0: float,
    alpha3_inv_mu0: float,
    betas: BetaCoefficients,
    mu_min: float,
    mu_max: float,
    n: int = 2000,
) -> tuple[float, float, float, float, float]:
    """
    Scan log-spaced mu in [mu_min, mu_max] and return best point:

    (mu_best, score_best, a1_inv, a2_inv, a3_inv)
    """
    if mu_min <= 0 or mu_max <= 0 or mu_min >= mu_max:
        raise ValueError("mu_min/mu_max must be positive and mu_min < mu_max")
    if n < 2:
        raise ValueError("n must be >= 2")

    log_min = math.log(mu_min)
    log_max = math.log(mu_max)

    best_mu = None
    best_score = float("inf")
    best_vals = (float("nan"), float("nan"), float("nan"))

    for i in range(n):
        t = i / (n - 1)
        mu = math.exp(log_min + (log_max - log_min) * t)
        a1 = run_alpha_inv(alpha1_inv_mu0, mu0, mu, betas.b1)
        a2 = run_alpha_inv(alpha2_inv_mu0, mu0, mu, betas.b2)
        a3 = run_alpha_inv(alpha3_inv_mu0, mu0, mu, betas.b3)
        score = converge_score(a1, a2, a3)
        if score < best_score:
            best_score = score
            best_mu = mu
            best_vals = (a1, a2, a3)

    assert best_mu is not None
    a1, a2, a3 = best_vals
    return best_mu, best_score, a1, a2, a3


# ---------------------------------------------------------------------------
# 2-loop coupled RG running  (RK4 integrator)
# ---------------------------------------------------------------------------

def run_gut_2loop(
    *,
    Q0_GeV: float,
    Q_GeV: float,
    alpha1_0: float,
    alpha2_0: float,
    alpha3_0: float,
    betas2: Beta2Loop,
    steps_per_unit_log: int = 500,
) -> Tuple[float, float, float]:
    """
    2-loop coupled RG running for all three gauge couplings.

    Derived from d(α_i⁻¹)/dt = -b_i/(2π) - Σ_j B_ij α_j/(8π²),
    which gives (in terms of α_i directly):

      dα_i/dt = (b_i/(2π)) α_i²  +  Σ_j (B_ij/(8π²)) α_i² α_j

    where t = ln μ.  Returns (α₁(Q), α₂(Q), α₃(Q)).
    """
    if Q0_GeV <= 0 or Q_GeV <= 0:
        raise ValueError("Q0_GeV and Q_GeV must be positive")

    t0 = math.log(Q0_GeV)
    t1 = math.log(Q_GeV)
    if t0 == t1:
        return alpha1_0, alpha2_0, alpha3_0

    b = betas2.b
    B = betas2.B
    two_pi = 2.0 * math.pi
    eight_pi2 = 8.0 * math.pi * math.pi

    def deriv(a: List[float]) -> List[float]:
        da = [0.0, 0.0, 0.0]
        for i in range(3):
            ai2 = a[i] * a[i]
            one_loop = (b[i] / two_pi) * ai2
            two_loop = 0.0
            for j in range(3):
                two_loop += (B[i][j] / eight_pi2) * ai2 * a[j]
            da[i] = one_loop + two_loop
        return da

    total = abs(t1 - t0)
    n_steps = max(20, int(total * steps_per_unit_log))
    dt = (t1 - t0) / n_steps

    a = [alpha1_0, alpha2_0, alpha3_0]
    for _ in range(n_steps):
        if any(not math.isfinite(x) or x <= 0 for x in a):
            return (float("inf"), float("inf"), float("inf"))
        if any(x > 10.0 for x in a):
            return (float("inf"), float("inf"), float("inf"))

        k1 = deriv(a)
        a_tmp = [a[i] + 0.5 * dt * k1[i] for i in range(3)]
        k2 = deriv(a_tmp)
        a_tmp = [a[i] + 0.5 * dt * k2[i] for i in range(3)]
        k3 = deriv(a_tmp)
        a_tmp = [a[i] + dt * k3[i] for i in range(3)]
        k4 = deriv(a_tmp)
        a = [a[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) for i in range(3)]

    return (a[0], a[1], a[2])


def find_best_convergence_2loop(
    *,
    mu0: float,
    alpha1_mu0: float,
    alpha2_mu0: float,
    alpha3_mu0: float,
    betas2: Beta2Loop,
    mu_min: float = 1e3,
    mu_max: float = 1e19,
    n: int = 3000,
    steps_per_unit_log: int = 500,
) -> Tuple[float, float, float, float, float]:
    """
    Scan log-spaced mu and return (mu_best, score, inv_a1, inv_a2, inv_a3)
    using 2-loop coupled running.
    """
    if mu_min <= 0 or mu_max <= 0 or mu_min >= mu_max:
        raise ValueError("mu_min/mu_max must be positive and mu_min < mu_max")
    if n < 2:
        raise ValueError("n must be >= 2")

    log_min = math.log(mu_min)
    log_max = math.log(mu_max)

    best_mu: Optional[float] = None
    best_score = float("inf")
    best_inv = (float("nan"), float("nan"), float("nan"))

    for i in range(n):
        t = i / (n - 1)
        mu = math.exp(log_min + (log_max - log_min) * t)
        a1, a2, a3 = run_gut_2loop(
            Q0_GeV=mu0, Q_GeV=mu,
            alpha1_0=alpha1_mu0, alpha2_0=alpha2_mu0, alpha3_0=alpha3_mu0,
            betas2=betas2, steps_per_unit_log=steps_per_unit_log,
        )
        if any(not math.isfinite(x) or x <= 0 for x in (a1, a2, a3)):
            continue
        inv1, inv2, inv3 = 1.0 / a1, 1.0 / a2, 1.0 / a3
        score = converge_score(inv1, inv2, inv3)
        if score < best_score:
            best_score = score
            best_mu = mu
            best_inv = (inv1, inv2, inv3)

    assert best_mu is not None
    return best_mu, best_score, best_inv[0], best_inv[1], best_inv[2]


# ---------------------------------------------------------------------------
# Non-minimal GUT normalization scan  (Item 5)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NormScanResult:
    k1: float
    name: str
    mu_best: float
    score: float
    inv_a1: float
    inv_a2: float
    inv_a3: float


def scan_gut_normalizations(
    *,
    alpha_mZ: float,
    cos2_mZ: float,
    alpha_s_mZ: float,
    sin2_mZ: float,
    mu0: float,
    betas: BetaCoefficients,
    norms: Optional[List[GUTNormalization]] = None,
    mu_min: float = 1e3,
    mu_max: float = 1e19,
    n: int = 3000,
) -> List[NormScanResult]:
    """
    For each GUT normalization factor k₁, recompute α₁^GUT(m_Z) = k₁·α_Y(m_Z)
    and find the best convergence scale.  Returns results sorted by score (best first).
    """
    if norms is None:
        norms = GUT_NORMALIZATIONS

    alpha_Y = alpha_mZ / cos2_mZ
    inv_a2 = sin2_mZ / alpha_mZ  # 1/α₂ = sin²θ/α
    inv_a3 = 1.0 / alpha_s_mZ

    results: List[NormScanResult] = []
    for norm in norms:
        alpha1_gut = norm.k1 * alpha_Y
        inv_a1 = 1.0 / alpha1_gut

        mu_best, score, a1, a2, a3 = find_best_convergence(
            mu0=mu0,
            alpha1_inv_mu0=inv_a1,
            alpha2_inv_mu0=inv_a2,
            alpha3_inv_mu0=inv_a3,
            betas=betas,
            mu_min=mu_min,
            mu_max=mu_max,
            n=n,
        )
        results.append(NormScanResult(
            k1=norm.k1, name=norm.name,
            mu_best=mu_best, score=score,
            inv_a1=a1, inv_a2=a2, inv_a3=a3,
        ))

    results.sort(key=lambda r: r.score)
    return results


def find_optimal_k1(
    *,
    alpha_mZ: float,
    cos2_mZ: float,
    alpha_s_mZ: float,
    sin2_mZ: float,
    mu0: float,
    betas: BetaCoefficients,
    k1_min: float = 0.5,
    k1_max: float = 4.0,
    n_k1: int = 500,
    mu_min: float = 1e3,
    mu_max: float = 1e19,
    n_mu: int = 2000,
) -> Tuple[float, float, float, float, float, float]:
    """
    Brute-force scan over continuous k₁ values to find the one that minimizes
    the convergence score.

    Returns (k1_best, mu_best, score, inv_a1, inv_a2, inv_a3).
    """
    alpha_Y = alpha_mZ / cos2_mZ
    inv_a2 = sin2_mZ / alpha_mZ
    inv_a3 = 1.0 / alpha_s_mZ

    best_k1 = float("nan")
    best_mu = float("nan")
    best_score = float("inf")
    best_vals = (float("nan"), float("nan"), float("nan"))

    for ik in range(n_k1):
        k1 = k1_min + (k1_max - k1_min) * ik / (n_k1 - 1)
        inv_a1 = 1.0 / (k1 * alpha_Y)

        mu_best, score, a1, a2, a3 = find_best_convergence(
            mu0=mu0,
            alpha1_inv_mu0=inv_a1,
            alpha2_inv_mu0=inv_a2,
            alpha3_inv_mu0=inv_a3,
            betas=betas,
            mu_min=mu_min,
            mu_max=mu_max,
            n=n_mu,
        )
        if score < best_score:
            best_score = score
            best_k1 = k1
            best_mu = mu_best
            best_vals = (a1, a2, a3)

    return best_k1, best_mu, best_score, best_vals[0], best_vals[1], best_vals[2]


# ---------------------------------------------------------------------------
# Lattice-constrained GUT search  (Item 6)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LatticeGUTResult:
    Q_GeV: float
    C: float
    m: int
    lattice_value: float  # C / phi^m
    inv_a1: float
    inv_a2: float
    inv_a3: float
    max_dev: float  # max |α_i⁻¹(Q) - C/φ^m| across all 3 couplings
    rms_dev: float


def find_lattice_gut_point(
    *,
    mu0: float,
    alpha1_inv_mu0: float,
    alpha2_inv_mu0: float,
    alpha3_inv_mu0: float,
    betas: BetaCoefficients,
    Cs: List[float],
    m_range: List[int],
    mu_min: float = 1e3,
    mu_max: float = 1e19,
    n_mu: int = 5000,
    top: int = 10,
) -> List[LatticeGUTResult]:
    """
    Search for scales Q where all three running couplings simultaneously
    approach a single lattice point (C, m).

    For each Q in a log-spaced scan, run all 3 couplings to Q, then for each
    candidate (C, m), compute the max deviation of the three inverse couplings
    from C/φ^m.  Return the top results sorted by max_dev.
    """
    if mu_min <= 0 or mu_max <= 0 or mu_min >= mu_max:
        raise ValueError("mu_min/mu_max must be positive")

    lattice_vals: List[Tuple[float, int, float]] = []
    for C in Cs:
        for m_int in m_range:
            val = C / (PHI ** m_int)
            if 0.1 < val < 200.0:
                lattice_vals.append((C, m_int, val))

    log_min = math.log(mu_min)
    log_max = math.log(mu_max)

    results: List[LatticeGUTResult] = []

    for i in range(n_mu):
        t = i / (n_mu - 1)
        mu = math.exp(log_min + (log_max - log_min) * t)
        a1 = run_alpha_inv(alpha1_inv_mu0, mu0, mu, betas.b1)
        a2 = run_alpha_inv(alpha2_inv_mu0, mu0, mu, betas.b2)
        a3 = run_alpha_inv(alpha3_inv_mu0, mu0, mu, betas.b3)

        if any(x <= 0 for x in (a1, a2, a3)):
            continue

        for C, m_int, lv in lattice_vals:
            d1 = abs(a1 - lv)
            d2 = abs(a2 - lv)
            d3 = abs(a3 - lv)
            max_d = max(d1, d2, d3)
            rms_d = math.sqrt((d1 * d1 + d2 * d2 + d3 * d3) / 3.0)
            results.append(LatticeGUTResult(
                Q_GeV=mu, C=C, m=m_int, lattice_value=lv,
                inv_a1=a1, inv_a2=a2, inv_a3=a3,
                max_dev=max_d, rms_dev=rms_d,
            ))

    results.sort(key=lambda r: r.max_dev)
    return results[:top]


def find_lattice_gut_point_2loop(
    *,
    mu0: float,
    alpha1_mu0: float,
    alpha2_mu0: float,
    alpha3_mu0: float,
    betas2: Beta2Loop,
    Cs: List[float],
    m_range: List[int],
    mu_min: float = 1e3,
    mu_max: float = 1e19,
    n_mu: int = 3000,
    steps_per_unit_log: int = 300,
    top: int = 10,
) -> List[LatticeGUTResult]:
    """2-loop version of find_lattice_gut_point."""
    if mu_min <= 0 or mu_max <= 0 or mu_min >= mu_max:
        raise ValueError("mu_min/mu_max must be positive")

    lattice_vals: List[Tuple[float, int, float]] = []
    for C in Cs:
        for m_int in m_range:
            val = C / (PHI ** m_int)
            if 0.1 < val < 200.0:
                lattice_vals.append((C, m_int, val))

    log_min = math.log(mu_min)
    log_max = math.log(mu_max)
    results: List[LatticeGUTResult] = []

    for i in range(n_mu):
        t = i / (n_mu - 1)
        mu = math.exp(log_min + (log_max - log_min) * t)
        a1c, a2c, a3c = run_gut_2loop(
            Q0_GeV=mu0, Q_GeV=mu,
            alpha1_0=alpha1_mu0, alpha2_0=alpha2_mu0, alpha3_0=alpha3_mu0,
            betas2=betas2, steps_per_unit_log=steps_per_unit_log,
        )
        if any(not math.isfinite(x) or x <= 0 for x in (a1c, a2c, a3c)):
            continue
        inv1, inv2, inv3 = 1.0 / a1c, 1.0 / a2c, 1.0 / a3c

        for C, m_int, lv in lattice_vals:
            d1 = abs(inv1 - lv)
            d2 = abs(inv2 - lv)
            d3 = abs(inv3 - lv)
            max_d = max(d1, d2, d3)
            rms_d = math.sqrt((d1 * d1 + d2 * d2 + d3 * d3) / 3.0)
            results.append(LatticeGUTResult(
                Q_GeV=mu, C=C, m=m_int, lattice_value=lv,
                inv_a1=inv1, inv_a2=inv2, inv_a3=inv3,
                max_dev=max_d, rms_dev=rms_d,
            ))

    results.sort(key=lambda r: r.max_dev)
    return results[:top]


