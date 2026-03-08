"""
GUT validation suite: independent lattice predictions, tightened 2-loop
search, SU(5) threshold corrections, and Fibonacci structure analysis.
"""

from __future__ import annotations

import bisect
import math
from dataclasses import dataclass
from typing import List, Tuple

PHI = (1.0 + math.sqrt(5.0)) / 2.0
LN_PHI = math.log(PHI)


# ---------------------------------------------------------------------------
# Precomputed lattice for fast lookups
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LatticePoint:
    value: float
    C: float
    m: int


def build_sorted_lattice(Cs: List[float], m_range: List[int]) -> List[LatticePoint]:
    pts = []
    for C in Cs:
        for m_int in m_range:
            lv = C / (PHI ** m_int)
            if lv > 0 and math.isfinite(lv):
                pts.append(LatticePoint(lv, C, m_int))
    pts.sort(key=lambda p: p.value)
    return pts


def nearest_lattice(value: float, lattice: List[LatticePoint]) -> Tuple[float, int, float, float]:
    """Binary-search nearest lattice point. Returns (C, m, lattice_value, abs_err)."""
    if not lattice or value <= 0:
        return 0.0, 0, 0.0, float("inf")
    vals = [p.value for p in lattice]
    idx = bisect.bisect_left(vals, value)
    best_C, best_m, best_lv, best_err = 0.0, 0, 0.0, float("inf")
    for i in range(max(0, idx - 1), min(len(lattice), idx + 2)):
        p = lattice[i]
        err = abs(value - p.value)
        if err < best_err:
            best_err = err
            best_C, best_m, best_lv = p.C, p.m, p.value
    return best_C, best_m, best_lv, best_err


# ---------------------------------------------------------------------------
# 1. Independent scale predictions
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ScalePrediction:
    name: str
    value: float
    note: str
    C: float
    m: int
    lattice_value: float
    rel_err: float


def independent_scale_tests(lattice: List[LatticePoint]) -> List[ScalePrediction]:
    from physics_test.rg_scales import lambda_qcd_from_alpha_s

    mZ = 91.1876
    mW = 80.379
    alpha_s = 0.1179
    v_ew = 246.22
    m_e = 0.000510999
    m_mu = 0.105658
    m_tau = 1.77686
    m_p = 0.938272
    m_t = 172.76
    m_b = 4.18
    m_H = 125.25
    m_c = 1.27
    m_s = 0.093

    lqcd = lambda_qcd_from_alpha_s(alpha_s_mu=alpha_s, mu_GeV=mZ, n_f=5, loops=2)
    m_nu = 0.05e-9
    M_seesaw = v_ew ** 2 / (2.0 * m_nu)
    Lambda_CC_quarter = 2.4e-12

    ratios = [
        ("m_p/m_e", m_p / m_e, "proton-electron mass ratio"),
        ("m_μ/m_e", m_mu / m_e, "muon-electron mass ratio"),
        ("m_τ/m_μ", m_tau / m_mu, "tau-muon mass ratio"),
        ("m_τ/m_e", m_tau / m_e, "tau-electron mass ratio"),
        ("m_t/m_b", m_t / m_b, "top-bottom mass ratio"),
        ("m_t/m_τ", m_t / m_tau, "top-tau mass ratio"),
        ("m_b/m_τ", m_b / m_tau, "bottom-tau mass ratio"),
        ("m_c/m_s", m_c / m_s, "charm-strange mass ratio"),
        ("m_t/m_W", m_t / mW, "top-W mass ratio"),
        ("m_H/m_W", m_H / mW, "Higgs-W mass ratio"),
        ("mZ/mW", mZ / mW, "Z-W mass ratio"),
        ("v/mZ", v_ew / mZ, "EW VEV / Z mass"),
        ("mZ/Λ_QCD", mZ / lqcd.Lambda_GeV, "Z / QCD scale"),
        ("M_seesaw/mZ", M_seesaw / mZ, "seesaw / Z mass"),
        ("M_seesaw/M_GUT", M_seesaw / 1.6e16, "seesaw / GUT"),
        ("G_F·mZ²", 1.1663787e-5 * mZ ** 2, "Fermi × Z²"),
        ("G_F·v²", 1.1663787e-5 * v_ew ** 2, "Fermi × VEV²"),
        ("m_p/Λ_QCD", m_p / lqcd.Lambda_GeV, "proton / QCD scale"),
        ("mZ/m_t", mZ / m_t, "Z / top mass"),
        ("m_H/mZ", m_H / mZ, "Higgs / Z mass"),
    ]

    results = []
    for name, value, note in ratios:
        if value <= 0 or not math.isfinite(value):
            continue
        C, m_int, lv, _ = nearest_lattice(value, lattice)
        rel = (lv - value) / value if value != 0 else float("inf")
        results.append(ScalePrediction(name, value, note, C, m_int, lv, rel))
    return results


def lattice_coverage(lattice: List[LatticePoint], tol: float, n_samples: int = 5000) -> float:
    """Null-hypothesis: fraction of random log-uniform values within tol of a lattice point."""
    import random
    rng = random.Random(42)
    hits = 0
    for _ in range(n_samples):
        val = 10.0 ** rng.uniform(-3, 32)
        _, _, lv, _ = nearest_lattice(val, lattice)
        if lv > 0 and abs((lv - val) / val) < tol:
            hits += 1
    return hits / n_samples


# ---------------------------------------------------------------------------
# 2. Tightened 2-loop search around (C=15, m=-1)
# ---------------------------------------------------------------------------

def tightened_2loop_search(
    *,
    Q0_GeV: float,
    alpha1_0: float,
    alpha2_0: float,
    alpha3_0: float,
    C_target: float = 15.0,
    m_target: int = -1,
    Q_center: float = 2e16,
    Q_half_decades: float = 0.7,
    n_Q: int = 300,
    steps_per_unit_log: int = 300,
) -> dict:
    from physics_test.gut import run_gut_2loop, MSSM_2LOOP, converge_score

    lattice_value = C_target / (PHI ** m_target)
    Q_min = Q_center * 10 ** (-Q_half_decades)
    Q_max = Q_center * 10 ** (Q_half_decades)
    log_min = math.log(Q_min)
    log_max = math.log(Q_max)

    best = None
    best_conv = None
    trajectory: List[dict] = []

    for i in range(n_Q):
        t = i / (n_Q - 1) if n_Q > 1 else 0.5
        Q = math.exp(log_min + (log_max - log_min) * t)
        a1, a2, a3 = run_gut_2loop(
            Q0_GeV=Q0_GeV, Q_GeV=Q,
            alpha1_0=alpha1_0, alpha2_0=alpha2_0, alpha3_0=alpha3_0,
            betas2=MSSM_2LOOP, steps_per_unit_log=steps_per_unit_log,
        )
        if any(not math.isfinite(x) or x <= 0 for x in (a1, a2, a3)):
            continue

        inv1, inv2, inv3 = 1.0 / a1, 1.0 / a2, 1.0 / a3
        d1 = inv1 - lattice_value
        d2 = inv2 - lattice_value
        d3 = inv3 - lattice_value
        max_dev = max(abs(d1), abs(d2), abs(d3))
        rms_dev = math.sqrt((d1 ** 2 + d2 ** 2 + d3 ** 2) / 3.0)
        score = converge_score(inv1, inv2, inv3)

        entry = dict(Q=Q, inv1=inv1, inv2=inv2, inv3=inv3,
                     d1=d1, d2=d2, d3=d3,
                     max_dev=max_dev, rms_dev=rms_dev, score=score)
        trajectory.append(entry)

        if best is None or max_dev < best['max_dev']:
            best = entry
        if best_conv is None or score < best_conv['score']:
            best_conv = entry

    return dict(C=C_target, m=m_target, lattice_value=lattice_value,
                best=best, best_conv=best_conv, trajectory=trajectory)


# ---------------------------------------------------------------------------
# 3. SU(5) threshold corrections with lattice-quantized mass splittings
# ---------------------------------------------------------------------------

def scan_su5_thresholds(
    *,
    inv_a1_gut: float,
    inv_a2_gut: float,
    inv_a3_gut: float,
    delta_range: range = range(-5, 6),
    susy: bool = True,
) -> List[dict]:
    """
    Scan integer phi-rung mass splittings for X/Y bosons and colored Higgs.
    M_V = M_GUT * phi^delta_V,  M_HC = M_GUT * phi^delta_HC
    """
    from physics_test.gut import converge_score

    two_pi = 2.0 * math.pi

    if susy:
        b1_V, b2_V, b3_V = 10.0 / 3.0, 2.0, 2.0
        b1_HC, b2_HC, b3_HC = 2.0 / 5.0, 0.0, 1.0 / 2.0
    else:
        b1_V, b2_V, b3_V = 5.0 / 3.0, 1.0, 1.0
        b1_HC, b2_HC, b3_HC = 1.0 / 15.0, 0.0, 1.0 / 6.0

    no_corr_score = converge_score(inv_a1_gut, inv_a2_gut, inv_a3_gut)

    results = []
    for dV in delta_range:
        for dHC in delta_range:
            ln_MV = float(dV) * LN_PHI
            ln_MHC = float(dHC) * LN_PHI
            d1 = -(b1_V / two_pi) * ln_MV - (b1_HC / two_pi) * ln_MHC
            d2 = -(b2_V / two_pi) * ln_MV - (b2_HC / two_pi) * ln_MHC
            d3 = -(b3_V / two_pi) * ln_MV - (b3_HC / two_pi) * ln_MHC

            ca1 = inv_a1_gut + d1
            ca2 = inv_a2_gut + d2
            ca3 = inv_a3_gut + d3
            if any(x <= 0 for x in (ca1, ca2, ca3)):
                continue

            score = converge_score(ca1, ca2, ca3)
            results.append(dict(
                delta_V=dV, delta_HC=dHC,
                inv_a1=ca1, inv_a2=ca2, inv_a3=ca3,
                score=score, improvement=no_corr_score - score,
            ))

    results.sort(key=lambda r: r['score'])
    return results


# ---------------------------------------------------------------------------
# 4. Fibonacci / golden-ratio structure analysis
# ---------------------------------------------------------------------------

def fibonacci_up_to(n: int) -> List[int]:
    fibs = [1, 1]
    while fibs[-1] < n:
        fibs.append(fibs[-1] + fibs[-2])
    return fibs


def lucas_up_to(n: int) -> List[int]:
    seq = [2, 1]
    while seq[-1] < n:
        seq.append(seq[-1] + seq[-2])
    return seq


def zeckendorf(n: int) -> List[int]:
    """Non-consecutive Fibonacci summands (Zeckendorf's theorem)."""
    if n <= 0:
        return []
    fibs = fibonacci_up_to(n)
    result = []
    remaining = n
    for f in reversed(fibs):
        if f <= remaining:
            result.append(f)
            remaining -= f
        if remaining == 0:
            break
    return result


def fibonacci_index(n: int) -> int | None:
    """If n == F(k), return k (1-indexed). Else None."""
    fibs = fibonacci_up_to(n + 1)
    for i, f in enumerate(fibs):
        if f == n:
            return i + 1
    return None


def lucas_index(n: int) -> int | None:
    """If n == L(k), return k (0-indexed). Else None."""
    seq = lucas_up_to(n + 1)
    for i, v in enumerate(seq):
        if v == n:
            return i
    return None


def _test_ratios_for_base(
    base: float,
    include: tuple,
    m_range: List[int],
) -> Tuple[List[float], List[LatticePoint], List[ScalePrediction]]:
    """Build a lattice from a given base and run independent scale tests."""
    from physics_test.gauge_groups import candidate_Cs_from_group, standard_model_gauge_groups

    all_Cs: List[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=base, include=include)
        for _, v in cs.items():
            all_Cs.append(float(v))
    all_Cs = sorted(set(all_Cs))
    lattice = build_sorted_lattice(all_Cs, m_range)
    preds = independent_scale_tests(lattice)
    return all_Cs, lattice, preds


@dataclass(frozen=True)
class BaseEnrichmentResult:
    base: float
    n_Cs: int
    Cs: List[float]
    hits_1pct: int
    hits_3pct: int
    hits_5pct: int
    total: int
    null_1pct: float
    null_3pct: float
    null_5pct: float
    enrichment_1pct: float
    enrichment_3pct: float
    enrichment_5pct: float


def base_permutation_test(
    bases: List[float],
    include: tuple = ("base", "base/dim", "base/coxeter", "base/dual_coxeter", "base/(dim*coxeter)"),
    m_range: List[int] | None = None,
    null_samples: int = 3000,
) -> List[BaseEnrichmentResult]:
    """Strategy A: run independent predictions for multiple base values and compare enrichment."""
    if m_range is None:
        m_range = list(range(-80, 120))

    results = []
    for base in bases:
        all_Cs, lattice, preds = _test_ratios_for_base(base, include, m_range)
        total = len(preds)
        h1 = sum(1 for p in preds if abs(p.rel_err) < 0.01)
        h3 = sum(1 for p in preds if abs(p.rel_err) < 0.03)
        h5 = sum(1 for p in preds if abs(p.rel_err) < 0.05)

        n1 = lattice_coverage(lattice, 0.01, null_samples)
        n3 = lattice_coverage(lattice, 0.03, null_samples)
        n5 = lattice_coverage(lattice, 0.05, null_samples)

        obs_1 = h1 / total if total > 0 else 0
        obs_3 = h3 / total if total > 0 else 0
        obs_5 = h5 / total if total > 0 else 0

        results.append(BaseEnrichmentResult(
            base=base, n_Cs=len(all_Cs), Cs=all_Cs,
            hits_1pct=h1, hits_3pct=h3, hits_5pct=h5, total=total,
            null_1pct=n1, null_3pct=n3, null_5pct=n5,
            enrichment_1pct=obs_1 / n1 if n1 > 0 else float("inf"),
            enrichment_3pct=obs_3 / n3 if n3 > 0 else float("inf"),
            enrichment_5pct=obs_5 / n5 if n5 > 0 else float("inf"),
        ))
    return results


# ---------------------------------------------------------------------------
# Strategy B: C-value clustering test
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CClusterResult:
    C: float
    count: int
    ratio_names: List[str]
    m_values: List[int]
    p_value: float


def c_clustering_test(
    preds: List[ScalePrediction],
    all_Cs: List[float],
    n_trials: int = 100000,
) -> List[CClusterResult]:
    """
    Strategy B: test whether the observed clustering of ratios onto specific
    C values is more extreme than expected by chance.

    For each C value that has k hits, compute P(k or more on any single C)
    under the null hypothesis that each ratio independently picks a C value
    with probability proportional to the lattice density of that C.
    """
    import random
    from collections import Counter

    observed_C_counts: dict[float, list] = {}
    for p in preds:
        c_val = p.C
        if c_val not in observed_C_counts:
            observed_C_counts[c_val] = []
        observed_C_counts[c_val].append((p.name, p.m))

    n_ratios = len(preds)
    n_C = len(all_Cs)

    if n_C == 0:
        return []

    rng = random.Random(42)
    max_cluster_dist: list[int] = []
    for _ in range(n_trials):
        trial_counts = Counter(rng.randint(0, n_C - 1) for _ in range(n_ratios))
        max_cluster_dist.append(max(trial_counts.values()) if trial_counts else 0)

    results = []
    for C_val, entries in sorted(observed_C_counts.items(), key=lambda x: -len(x[1])):
        k = len(entries)
        p_val = sum(1 for mc in max_cluster_dist if mc >= k) / n_trials
        results.append(CClusterResult(
            C=C_val, count=k,
            ratio_names=[e[0] for e in entries],
            m_values=[e[1] for e in entries],
            p_value=p_val,
        ))

    results.sort(key=lambda r: -r.count)
    return results


# ---------------------------------------------------------------------------
# Strategy C: Pre-registered predictions (new dimensionless ratios)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PreRegisteredPrediction:
    name: str
    value: float
    note: str
    predicted_C: float
    predicted_m: int
    predicted_lattice: float
    predicted_rel_err: float


def pre_registered_predictions(lattice: List[LatticePoint]) -> List[PreRegisteredPrediction]:
    """
    Strategy C: dimensionless ratios NOT in the original test set.
    These serve as out-of-sample validation.
    """
    m_pi_charged = 0.13957  # GeV, charged pion
    m_pi_neutral = 0.13498  # GeV, neutral pion
    m_p = 0.938272
    m_e = 0.000510999
    m_mu = 0.105658
    m_n = 0.939565  # neutron
    m_t = 172.76
    m_W = 80.379
    m_Z = 91.1876
    m_H = 125.25
    m_c = 1.27
    m_b = 4.18
    m_tau = 1.77686
    m_s = 0.093
    m_d = 0.0047
    m_u = 0.0022
    alpha_em = 1.0 / 137.036

    sin2_cabibbo = 0.2253 ** 2
    sin2_12 = 0.307
    sin2_23 = 0.546
    sin2_13 = 0.0220
    jarlskog = 3.18e-5
    a_e = 1.15966e-3

    new_ratios = [
        ("sin(θ_C)", 0.2253, "Cabibbo angle sine"),
        ("sin²(θ_C)", sin2_cabibbo, "Cabibbo angle sin²"),
        ("sin²θ₁₂", sin2_12, "solar neutrino mixing"),
        ("sin²θ₂₃", sin2_23, "atmospheric neutrino mixing"),
        ("sin²θ₁₃", sin2_13, "reactor neutrino mixing"),
        ("J_CP", jarlskog, "Jarlskog invariant"),
        ("a_e", a_e, "electron (g-2)/2"),
        ("1/a_e", 1.0 / a_e, "inverse electron (g-2)/2"),
        ("m_π/m_p", m_pi_charged / m_p, "pion-proton ratio"),
        ("m_π⁺/m_π⁰", m_pi_charged / m_pi_neutral, "charged/neutral pion"),
        ("m_n/m_p", m_n / m_p, "neutron-proton ratio"),
        ("(m_n-m_p)/m_e", (m_n - m_p) / m_e, "n-p mass diff / electron"),
        ("m_c/m_τ", m_c / m_tau, "charm-tau mass ratio"),
        ("m_b/m_c", m_b / m_c, "bottom-charm mass ratio"),
        ("m_s/m_d", m_s / m_d, "strange-down mass ratio"),
        ("m_u/m_d", m_u / m_d, "up-down mass ratio"),
        ("m_t/m_H", m_t / m_H, "top-Higgs mass ratio"),
        ("m_W/m_H", m_W / m_H, "W-Higgs mass ratio"),
        ("m_t·m_b/m_W²", m_t * m_b / m_W ** 2, "Yukawa product"),
        ("α/(2π)", alpha_em / (2 * math.pi), "Schwinger coefficient"),
    ]

    results = []
    for name, value, note in new_ratios:
        if value <= 0 or not math.isfinite(value):
            continue
        C, m_int, lv, _ = nearest_lattice(value, lattice)
        rel = (lv - value) / value if value != 0 else float("inf")
        results.append(PreRegisteredPrediction(name, value, note, C, m_int, lv, rel))
    return results


# ---------------------------------------------------------------------------
# Extended analyses: CKM/PMNS, lattice operations, 360-family, n_gap check
# ---------------------------------------------------------------------------

def full_mixing_matrix_predictions(lattice: List[LatticePoint]) -> List[PreRegisteredPrediction]:
    """Full CKM and PMNS mixing matrix parameters."""
    import cmath

    # CKM parameters (PDG 2024)
    V_ud = 0.97373
    V_us = 0.2243   # = sin(theta_Cabibbo)
    V_ub = 0.00382
    V_cd = 0.221
    V_cs = 0.975
    V_cb = 0.0408
    V_td = 0.0080
    V_ts = 0.0388
    V_tb = 0.99917

    J_CKM = 3.18e-5  # Jarlskog CKM
    delta_CKM_rad = 1.144  # CP phase in radians

    # PMNS parameters (NuFIT 5.2, normal ordering)
    sin2_12_PMNS = 0.304
    sin2_23_PMNS = 0.573
    sin2_13_PMNS = 0.02219
    delta_CP_PMNS_rad = 3.42  # ~196 degrees

    ratios = [
        # CKM magnitudes
        ("|V_ud|", V_ud, "CKM ud element"),
        ("|V_us|", V_us, "CKM us = Cabibbo"),
        ("|V_ub|", V_ub, "CKM ub element"),
        ("|V_cd|", V_cd, "CKM cd element"),
        ("|V_cs|", V_cs, "CKM cs element"),
        ("|V_cb|", V_cb, "CKM cb element"),
        ("|V_td|", V_td, "CKM td element"),
        ("|V_ts|", V_ts, "CKM ts element"),
        ("|V_tb|", V_tb, "CKM tb element"),
        # CKM ratios (more lattice-friendly — further from 1)
        ("|V_us/V_cb|", V_us / V_cb, "Cabibbo/cb ratio"),
        ("|V_cb/V_ub|", V_cb / V_ub, "cb/ub ratio"),
        ("|V_td/V_ts|", V_td / V_ts, "td/ts Wolfenstein"),
        ("|V_us|²", V_us ** 2, "lambda² Wolfenstein"),
        ("|V_cb|²", V_cb ** 2, "A²lambda⁴ approx"),
        ("J_CKM", J_CKM, "Jarlskog CKM invariant"),
        ("δ_CKM/π", delta_CKM_rad / math.pi, "CKM CP phase / pi"),
        # PMNS
        ("sin²θ₁₂_PMNS", sin2_12_PMNS, "solar neutrino PMNS"),
        ("sin²θ₂₃_PMNS", sin2_23_PMNS, "atmospheric PMNS"),
        ("sin²θ₁₃_PMNS", sin2_13_PMNS, "reactor PMNS"),
        ("δ_PMNS/π", delta_CP_PMNS_rad / math.pi, "PMNS CP phase / pi"),
        # Cross-sector ratios
        ("θ_C/θ₁₃_CKM", math.asin(V_us) / math.asin(V_ub), "Cabibbo / theta13"),
        ("sin²θ₁₂_PMNS/sin²θ_C", sin2_12_PMNS / V_us ** 2, "PMNS12 / CKM-Cabibbo²"),
    ]

    results = []
    for name, value, note in ratios:
        if value <= 0 or not math.isfinite(value):
            continue
        C, m_int, lv, _ = nearest_lattice(value, lattice)
        rel = (lv - value) / value if value != 0 else float("inf")
        results.append(PreRegisteredPrediction(name, value, note, C, m_int, lv, rel))
    return results


def lattice_operation_analysis(Cs: List[float], m_range: List[int]) -> List[dict]:
    """
    Analyze how arithmetic operations (multiply/divide by 2pi, sqrt, etc.)
    map lattice addresses. Reveals internal lattice algebra.
    """
    lattice = build_sorted_lattice(Cs, m_range)

    known_points = [
        ("1/α", 137.036, "EM coupling"),
        ("α/(2π)", 0.001161, "Schwinger coeff"),
        ("1/α_s(mZ)", 8.475, "strong coupling"),
        ("1/α₂(mZ)", 28.54, "weak coupling"),
        ("15φ", 15.0 * PHI, "GUT coupling"),
    ]

    operations = [
        ("×2π", lambda x: x * 2 * math.pi, "multiply by 2pi"),
        ("÷2π", lambda x: x / (2 * math.pi), "divide by 2pi"),
        ("×π", lambda x: x * math.pi, "multiply by pi"),
        ("÷π", lambda x: x / math.pi, "divide by pi"),
        ("×φ", lambda x: x * PHI, "multiply by phi"),
        ("÷φ", lambda x: x / PHI, "divide by phi"),
        ("√", lambda x: math.sqrt(x), "square root"),
        ("²", lambda x: x * x, "square"),
    ]

    results = []
    for pt_name, pt_val, pt_note in known_points:
        C0, m0, lv0, _ = nearest_lattice(pt_val, lattice)
        for op_name, op_fn, op_note in operations:
            new_val = op_fn(pt_val)
            if new_val <= 0 or not math.isfinite(new_val):
                continue
            C1, m1, lv1, _ = nearest_lattice(new_val, lattice)
            rel = (lv1 - new_val) / new_val if new_val != 0 else float("inf")
            results.append(dict(
                source=pt_name, source_addr=f"({C0:g},{m0})",
                operation=op_name,
                result_val=new_val,
                result_addr=f"({C1:g},{m1})",
                delta_C=C1 - C0, delta_m=m1 - m0,
                rel_err=rel,
            ))
    return results


def base_family_analysis(perm_results: list) -> dict:
    """Analyze whether enrichment correlates with divisibility by 360's prime factors."""
    import math as _m

    def shared_prime_factor_count(a: int, b: int) -> int:
        factors_360 = {2: 3, 3: 2, 5: 1}
        count = 0
        for p, exp in factors_360.items():
            ea = 0
            temp = a
            while temp % p == 0:
                ea += 1
                temp //= p
            count += min(ea, exp)
        return count

    def gcd(a: int, b: int) -> int:
        while b:
            a, b = b, a % b
        return a

    analysis = []
    for r in perm_results:
        b = int(round(r.base))
        g = gcd(b, 360)
        shared = shared_prime_factor_count(b, 360)
        divides_360 = (360 % b == 0) if b > 0 else False
        multiple_360 = (b % 360 == 0) if b > 0 else False
        analysis.append(dict(
            base=r.base,
            gcd_with_360=g,
            shared_prime_factors=shared,
            divides_360=divides_360,
            multiple_of_360=multiple_360,
            enrichment_1pct=r.enrichment_1pct,
            hits_1pct=r.hits_1pct,
        ))

    # Correlation between shared factors and enrichment
    xs = [a['shared_prime_factors'] for a in analysis]
    ys = [a['enrichment_1pct'] for a in analysis]
    n = len(xs)
    if n > 1:
        mx = sum(xs) / n
        my = sum(ys) / n
        cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / n
        sx = _m.sqrt(sum((x - mx) ** 2 for x in xs) / n)
        sy = _m.sqrt(sum((y - my) ** 2 for y in ys) / n)
        corr = cov / (sx * sy) if sx > 0 and sy > 0 else 0.0
    else:
        corr = 0.0

    # Same for gcd
    xs_g = [a['gcd_with_360'] for a in analysis]
    if n > 1:
        mx_g = sum(xs_g) / n
        cov_g = sum((x - mx_g) * (y - my) for x, y in zip(xs_g, ys)) / n
        sx_g = _m.sqrt(sum((x - mx_g) ** 2 for x in xs_g) / n)
        corr_g = cov_g / (sx_g * sy) if sx_g > 0 and sy > 0 else 0.0
    else:
        corr_g = 0.0

    return dict(
        entries=analysis,
        corr_shared_factors=corr,
        corr_gcd=corr_g,
    )


def phi_power_decompose(exponents: dict) -> dict:
    """
    Analyze Fibonacci/golden-ratio structure in energy hierarchy exponents.
    exponents: name -> float  (e.g. {"GUT/mZ": 69.6, "Pl/GUT": 12.4, "Pl/mZ": 82.0})
    """
    fibs_set = set(fibonacci_up_to(500))
    lucas_set = set(lucas_up_to(500))

    results = {}
    for name, exp in exponents.items():
        exp_int = round(exp)
        z = zeckendorf(exp_int)
        fi = fibonacci_index(exp_int)
        li = lucas_index(exp_int)

        nearby_fibs = []
        for delta in range(-3, 4):
            if (exp_int + delta) in fibs_set:
                nearby_fibs.append((delta, exp_int + delta,
                                    fibonacci_index(exp_int + delta)))

        results[name] = dict(
            exact=exp, rounded=exp_int, deviation=exp - exp_int,
            zeckendorf=z, is_fibonacci=exp_int in fibs_set,
            is_lucas=exp_int in lucas_set,
            fib_index=fi, lucas_index=li,
            nearby_fibs=nearby_fibs,
        )
    return results


# ---------------------------------------------------------------------------
# Systematic lattice address table
# ---------------------------------------------------------------------------

def build_address_table(
    lattice: List[LatticePoint],
) -> List[dict]:
    """
    Combine ALL known dimensionless ratios and map them onto the lattice.
    Returns list of dicts sorted by (C, m) for pattern discovery.
    """
    m_e = 0.000510999
    m_mu = 0.105658
    m_tau = 1.77686
    m_p = 0.938272
    m_n = 0.939565
    m_t = 172.76
    m_b = 4.18
    m_c = 1.27
    m_s = 0.093
    m_d = 0.0047
    m_u = 0.0022
    m_W = 80.379
    m_Z = 91.1876
    m_H = 125.25
    v_ew = 246.22
    alpha_em = 1.0 / 137.036
    alpha_s = 0.1179
    G_F = 1.1663787e-5
    m_pi = 0.13957
    M_Planck = 1.2209e19
    M_GUT = 1.6e16

    from physics_test.rg_scales import lambda_qcd_from_alpha_s
    lqcd = lambda_qcd_from_alpha_s(alpha_s_mu=alpha_s, mu_GeV=m_Z, n_f=5, loops=2)
    L_QCD = lqcd.Lambda_GeV

    all_ratios = [
        # --- Gauge couplings ---
        ("1/α_em", 137.036, "gauge", "electromagnetic"),
        ("1/α_s(mZ)", 1.0 / alpha_s, "gauge", "strong"),
        ("1/α₂(mZ)", 29.587, "gauge", "weak SU(2)"),
        ("α/(2π)", alpha_em / (2 * math.pi), "gauge", "Schwinger"),
        ("α_s(mZ)", alpha_s, "gauge", "strong coupling"),
        # --- Lepton mass ratios ---
        ("m_μ/m_e", m_mu / m_e, "lepton", "mu-e mass"),
        ("m_τ/m_e", m_tau / m_e, "lepton", "tau-e mass"),
        ("m_τ/m_μ", m_tau / m_mu, "lepton", "tau-mu mass"),
        # --- Quark mass ratios ---
        ("m_t/m_b", m_t / m_b, "quark", "top-bottom"),
        ("m_t/m_τ", m_t / m_tau, "quark", "top-tau"),
        ("m_b/m_τ", m_b / m_tau, "quark", "bottom-tau"),
        ("m_c/m_s", m_c / m_s, "quark", "charm-strange"),
        ("m_t/m_u", m_t / m_u, "quark", "top-up"),
        ("m_b/m_d", m_b / m_d, "quark", "bottom-down"),
        ("m_s/m_d", m_s / m_d, "quark", "strange-down"),
        ("m_c/m_u", m_c / m_u, "quark", "charm-up"),
        # --- Baryon/meson ---
        ("m_p/m_e", m_p / m_e, "baryon", "proton-electron"),
        ("m_p/m_π", m_p / m_pi, "baryon", "proton-pion"),
        ("m_n/m_p", m_n / m_p, "baryon", "neutron-proton"),
        ("(m_n-m_p)/m_e", (m_n - m_p) / m_e, "baryon", "n-p split / e"),
        # --- EW sector ---
        ("m_t/m_W", m_t / m_W, "EW", "top-W"),
        ("m_H/m_W", m_H / m_W, "EW", "Higgs-W"),
        ("m_Z/m_W", m_Z / m_W, "EW", "Z-W"),
        ("v/m_Z", v_ew / m_Z, "EW", "VEV/Z"),
        ("m_H/m_Z", m_H / m_Z, "EW", "Higgs-Z"),
        ("m_Z/m_t", m_Z / m_t, "EW", "Z-top"),
        ("G_F·m_Z²", G_F * m_Z ** 2, "EW", "Fermi*Z²"),
        # --- Scale hierarchies ---
        ("m_Z/Λ_QCD", m_Z / L_QCD, "hierarchy", "EW/QCD"),
        ("m_p/Λ_QCD", m_p / L_QCD, "hierarchy", "proton/QCD"),
        ("M_GUT/m_Z", M_GUT / m_Z, "hierarchy", "GUT/EW"),
        ("M_Pl/m_Z", M_Planck / m_Z, "hierarchy", "Planck/EW"),
        ("M_Pl/M_GUT", M_Planck / M_GUT, "hierarchy", "Planck/GUT"),
        # --- Mixing angles ---
        ("|V_us|", 0.2243, "CKM", "Cabibbo"),
        ("|V_cb|", 0.0408, "CKM", "cb element"),
        ("|V_ub|", 0.00382, "CKM", "ub element"),
        ("|V_us/V_cb|", 0.2243 / 0.0408, "CKM", "Cabibbo/cb"),
        ("|V_cb/V_ub|", 0.0408 / 0.00382, "CKM", "cb/ub"),
        ("sin²θ₁₂_PMNS", 0.304, "PMNS", "solar mixing"),
        ("sin²θ₂₃_PMNS", 0.573, "PMNS", "atmos mixing"),
        ("sin²θ₁₃_PMNS", 0.02219, "PMNS", "reactor mixing"),
        ("J_CKM", 3.18e-5, "CKM", "Jarlskog"),
        # --- Anomalous magnetic moments ---
        ("a_e", 1.15966e-3, "anomalous", "electron g-2"),
        ("a_μ", 1.16592e-3, "anomalous", "muon g-2"),
    ]

    table = []
    for name, value, sector, description in all_ratios:
        if value <= 0 or not math.isfinite(value):
            continue
        C, m_int, lv, _ = nearest_lattice(value, lattice)
        rel = (lv - value) / value if value != 0 else float("inf")
        table.append(dict(
            name=name, value=value, sector=sector, description=description,
            C=C, m=m_int, lattice_val=lv, rel_err=rel,
        ))

    table.sort(key=lambda r: (r['C'], r['m']))
    return table


# ---------------------------------------------------------------------------
# n_gap self-consistency test
# ---------------------------------------------------------------------------

def ngap_self_consistency() -> dict:
    """
    Test: if n_gap = dim(G_SM) = 12, does that correctly predict the
    GUT-to-Planck ratio via M_Pl/M_GUT = phi^n_gap?
    And: starting from M_Z, does n_total = 82 predict M_Pl correctly?
    """
    n_gap = 12
    n_GUT = 70  # approximate phi-power from MZ to M_GUT

    M_Z = 91.1876  # GeV
    M_GUT_nominal = 1.6e16
    M_Pl = 1.2209e19

    # Forward predictions
    M_Pl_from_gap = M_GUT_nominal * PHI ** n_gap
    M_GUT_from_total = M_Z * PHI ** n_GUT
    M_Pl_from_total = M_Z * PHI ** (n_GUT + n_gap)

    # Actual ratios
    ratio_Pl_GUT_actual = M_Pl / M_GUT_nominal
    ratio_Pl_GUT_pred = PHI ** n_gap
    ratio_Pl_Z_actual = M_Pl / M_Z
    ratio_Pl_Z_pred = PHI ** 82

    # Backward: infer n_gap from M_Pl/M_GUT
    n_gap_inferred = math.log(M_Pl / M_GUT_nominal) / math.log(PHI)
    n_total_inferred = math.log(M_Pl / M_Z) / math.log(PHI)

    # Check if n_gap = dim(G_SM) = dim(SU(3)) + dim(SU(2)) + dim(U(1)) = 8+3+1
    dim_SM = 8 + 3 + 1  # = 12

    return dict(
        n_gap=n_gap,
        n_GUT=n_GUT,
        n_total=n_GUT + n_gap,
        dim_SM=dim_SM,
        n_gap_equals_dim_SM=(n_gap == dim_SM),
        M_GUT_predicted=M_GUT_from_total,
        M_GUT_actual=M_GUT_nominal,
        M_GUT_rel_err=(M_GUT_from_total - M_GUT_nominal) / M_GUT_nominal,
        M_Pl_from_gap_predicted=M_Pl_from_gap,
        M_Pl_actual=M_Pl,
        M_Pl_gap_rel_err=(M_Pl_from_gap - M_Pl) / M_Pl,
        M_Pl_from_total_predicted=M_Pl_from_total,
        M_Pl_total_rel_err=(M_Pl_from_total - M_Pl) / M_Pl,
        ratio_Pl_GUT_actual=ratio_Pl_GUT_actual,
        ratio_Pl_GUT_predicted=ratio_Pl_GUT_pred,
        ratio_Pl_GUT_rel_err=(ratio_Pl_GUT_pred - ratio_Pl_GUT_actual) / ratio_Pl_GUT_actual,
        n_gap_inferred=n_gap_inferred,
        n_total_inferred=n_total_inferred,
    )
