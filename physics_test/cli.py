from __future__ import annotations

import argparse
import math
import sys

from physics_test import constants
from physics_test.forces import alpha_gravity
from physics_test.model import (
    coupling_invariant,
    evaluate_fit,
    fit_C_for_target_G,
    frequency_F0,
    gauge_G,
    temperature_K_from_frequency,
)
from physics_test.scan import filter_hits_by_rel_err, frange, scan_candidates
from physics_test.toponumbers import candidate_sets, get_candidate_set
from physics_test.targets import known_targets
from physics_test.oos import ew_sin2_suites, oos_suites, predictive_force_suites, resolve_oos_targets
from physics_test.gravity_bands import bands as gravity_band_list
from physics_test.units import (
    energy_J_from_GeV,
    energy_J_from_MeV,
    energy_J_from_eV,
    frequency_Hz_from_energy_J,
    temperature_K_from_energy_J,
)
from physics_test.presets import em_frequency_presets, get_preset, particle_proxy_presets, thermal_presets
from physics_test.gauge_groups import candidate_Cs_from_group, standard_model_gauge_groups
from physics_test.units import mass_kg_from_GeV
from physics_test.gut import (
    MSSM_1LOOP, SM_1LOOP, MSSM_2LOOP, SM_2LOOP,
    converge_score, find_best_convergence, find_best_convergence_2loop,
    run_alpha_inv,
    scan_gut_normalizations,
    find_lattice_gut_point, find_lattice_gut_point_2loop,
)
from physics_test.normalization import normalization_factor_for_force, normalization_families
from physics_test.steps import step_from_targets
from physics_test.rg_scales import lambda_qcd_from_alpha_s
from physics_test.oos_rg import rg_suites
from physics_test.rg_within_band import QCDRunSpec, alpha_s_from_ref, qcd_run_spec_from_key, scale_GeV_from_target_key
from physics_test.qed_running import alpha_inv_mZ_from_delta_alpha, qed_run_alpha_inv_1loop_from_ref
from physics_test.ew_mixing import sin2_eff_gammaZ_1loop_from_ref
from physics_test.target_registry import get_measurement


def _try_configure_utf8_stdio() -> None:
    """
    Windows consoles often default to legacy encodings (e.g. cp1252), which can crash
    `print()` if notes contain symbols like α, ≈, etc. Reconfigure to UTF-8 when possible.
    """

    try:
        if hasattr(sys.stdout, "reconfigure"):
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")  # type: ignore[attr-defined]
    except Exception:
        pass
    try:
        if hasattr(sys.stderr, "reconfigure"):
            sys.stderr.reconfigure(encoding="utf-8", errors="replace")  # type: ignore[attr-defined]
    except Exception:
        pass


_try_configure_utf8_stdio()


def _add_common_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--m", type=int, default=2, help="Index m (integer harmonic step) (default: 2)")
    p.add_argument("--C", type=float, default=360.0, help="Topological constant C (default: 360)")
    p.add_argument("--K", type=float, default=300.0, help="Temperature K in Kelvin (default: 300)")


def cmd_calc(args: argparse.Namespace) -> int:
    G = gauge_G(args.C, args.m)
    F0 = frequency_F0(args.m, args.K)
    inv = coupling_invariant(args.C, args.K)

    print(f"phi = {((1+math.sqrt(5))/2):.15f}")
    print(f"G = C/phi^m = {G:.15g}  (C={args.C}, m={args.m})")
    print(f"F0 = phi^m*kB*K/h = {F0:.15g} Hz  (K={args.K} K)")
    print(f"G*F0 = {G*F0:.15g} Hz")
    print(f"C*kB*K/h = {inv:.15g} Hz  (should match G*F0)")
    return 0


def cmd_fits(args: argparse.Namespace) -> int:
    targets = known_targets()
    print(f"Using built-in CODATA values; alpha = {constants.FINE_STRUCTURE:.15g}")
    print(f"m = {args.m}\n")
    for t in targets:
        C_fit = fit_C_for_target_G(t.value, args.m)
        r = evaluate_fit(t.value, C_fit, args.m)
        print(f"{t.name:12s} target={t.value:.15g}  C_fit={C_fit:.15g}  rel_err={r.rel_err:.3g}")
    return 0


def cmd_check_example(args: argparse.Namespace) -> int:
    # Your example: C=360, m=2.
    G = gauge_G(360.0, 2.0)
    alpha = constants.FINE_STRUCTURE
    inv_alpha = 1.0 / alpha
    print(f"G(360,2) = {G:.15g}")
    print(f"alpha     = {alpha:.15g}")
    print(f"1/alpha   = {inv_alpha:.15g}")
    print(f"delta vs 1/alpha = {G - inv_alpha:.15g}")
    return 0


def _parse_Cs(csv: str) -> list[float]:
    if not csv.strip():
        return []
    out: list[float] = []
    for raw in csv.split(","):
        raw = raw.strip()
        if not raw:
            continue
        out.append(float(raw))
    return out


def cmd_list_sets(args: argparse.Namespace) -> int:
    sets = candidate_sets()
    for s in sets:
        print(f"{s.name:20s}  n={len(s.values):2d}  note={s.note}")
    return 0


def cmd_list_targets(args: argparse.Namespace) -> int:
    targets = known_targets()
    for t in targets:
        sigma_s = f"{t.sigma:.3g}" if getattr(t, "sigma", None) is not None else "NA"
        q_s = f"{t.Q_GeV:g}" if getattr(t, "Q_GeV", None) is not None else ""
        scheme = getattr(t, "scheme", "") or ""
        if q_s:
            q_s = f"Q={q_s}GeV"
        meta = " ".join(x for x in [f"sigma={sigma_s}", q_s, scheme] if x)
        print(f"{t.name:28s}  {t.value:.15g}  {meta}  note={t.note}")
    return 0


def cmd_list_norm_families(args: argparse.Namespace) -> int:
    fams = normalization_families()
    for k in sorted(fams.keys()):
        print(f"{k:24s}  {fams[k].note}")
    return 0


def cmd_rg_scales(args: argparse.Namespace) -> int:
    """
    Small RG/dimensional-transmutation helpers that make the role of e explicit.

    The key point: many physical scale hierarchies are not multiplicative constants,
    but exponentials like exp(-const/alpha). This is where 'e' naturally enters.
    """

    # QCD Lambda from alpha_s(mu)
    res1 = lambda_qcd_from_alpha_s(alpha_s_mu=args.alpha_s, mu_GeV=args.mu_GeV, n_f=args.n_f, loops=1)
    res2 = lambda_qcd_from_alpha_s(alpha_s_mu=args.alpha_s, mu_GeV=args.mu_GeV, n_f=args.n_f, loops=2)

    # Show equivalent temperatures/frequencies as a bridge back to the F0/K side.
    E1_J = energy_J_from_GeV(res1.Lambda_GeV)
    E2_J = energy_J_from_GeV(res2.Lambda_GeV)
    T1_K = temperature_K_from_energy_J(E1_J)
    T2_K = temperature_K_from_energy_J(E2_J)
    f1_Hz = frequency_Hz_from_energy_J(E1_J)
    f2_Hz = frequency_Hz_from_energy_J(E2_J)

    # How many φ-steps is μ/Λ ?
    p = (1.0 + math.sqrt(5.0)) / 2.0
    dm1 = math.log(res1.mu_GeV / res1.Lambda_GeV) / math.log(p) if res1.Lambda_GeV > 0 else float("nan")
    dm2 = math.log(res2.mu_GeV / res2.Lambda_GeV) / math.log(p) if res2.Lambda_GeV > 0 else float("nan")

    print("QCD dimensional transmutation (Lambda_QCD) from alpha_s(mu):")
    print(f"  inputs: alpha_s(mu)={args.alpha_s:.8g}, mu={args.mu_GeV:.8g} GeV, n_f={args.n_f}")
    print(f"  beta0={res2.beta0:.8g}, beta1={res2.beta1:.8g}")
    print("")
    print("  1-loop:")
    print(f"    Lambda = {res1.Lambda_GeV:.8g} GeV")
    print(f"    equiv  = {T1_K:.6g} K,  f=E/h ≈ {f1_Hz:.6g} Hz")
    print(f"    log_phi(mu/Lambda) ≈ {dm1:.6f}")
    print("")
    print("  2-loop (approx):")
    print(f"    Lambda = {res2.Lambda_GeV:.8g} GeV")
    print(f"    equiv  = {T2_K:.6g} K,  f=E/h ≈ {f2_Hz:.6g} Hz")
    print(f"    log_phi(mu/Lambda) ≈ {dm2:.6f}")
    print("")
    print("Interpretation:")
    print("  - 'e' enters through exp(-const/alpha). This is *not* a free multiplicative factor;")
    print("    it is the RG mechanism that generates large scale separations from dimensionless couplings.")
    return 0


def cmd_oos_rg(args: argparse.Namespace) -> int:
    """
    RG+φ out-of-sample style report.

    This fits an anchor coupling on the φ-lattice, then uses the resulting lattice
    couplings at additional scales to compute an RG-generated scale (Lambda_QCD).
    """

    suites = rg_suites()
    if args.suite not in suites:
        raise SystemExit(f"Unknown RG suite {args.suite!r}. Options: {', '.join(sorted(suites.keys()))}")
    suite = suites[args.suite]

    target_map = {t.name: t.value for t in known_targets()}
    if suite.anchor.key not in target_map:
        raise SystemExit(f"Unknown anchor key {suite.anchor.key!r}. Run `list-targets`.")
    for t in suite.targets:
        if t.key not in target_map:
            raise SystemExit(f"Unknown target key {t.key!r}. Run `list-targets`.")

    # Integer m grid
    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))

    # Gauge-derived C candidates
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    cand: list[tuple[str, float]] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for k, v in cs.items():
            cand.append((f"{g.name}:{k}", float(v)))

    seen: set[float] = set()
    Cs: list[float] = []
    label_by_C: dict[float, str] = {}
    for lab, c in cand:
        if c in seen:
            continue
        seen.add(c)
        Cs.append(c)
        label_by_C[c] = lab

    print(f"RG+phi suite: {suite.key}")
    print(f"force = {suite.force}")
    print(f"base = {args.base}   include = {','.join(include)}   |Cs| = {len(Cs)}")
    print(f"m range = [{min(m_values)}, {max(m_values)}]")
    print(f"Lambda_QCD: loops={suite.loops}, n_f={suite.n_f}\n")

    anchor_target = float(target_map[suite.anchor.key])

    # For each C candidate, fit anchor m and compute implied Lambda at anchor + targets.
    rows: list[tuple[float, float, float, float, float, list[tuple[str, float, float]]]] = []
    for C in Cs:
        # Fit anchor
        hits_anchor = scan_candidates(Cs=[C], m_values=m_values, target_G=anchor_target)
        best_a = hits_anchor[0]
        if args.max_rel_err is not None and abs(best_a.rel_err) > float(args.max_rel_err):
            continue

        # Interpret G as an inverse coupling for strong targets (our suite is strong-only for now)
        inv_alpha_anchor = float(best_a.G)
        alpha_anchor = 1.0 / inv_alpha_anchor if inv_alpha_anchor != 0 else float("inf")
        lam_a = lambda_qcd_from_alpha_s(alpha_s_mu=alpha_anchor, mu_GeV=float(suite.anchor.Q_GeV), n_f=suite.n_f, loops=suite.loops)

        per: list[tuple[str, float, float]] = []  # (key, rel_err, Lambda_GeV)
        lambdas: list[float] = [lam_a.Lambda_GeV]
        ok = True
        for t in suite.targets:
            tgt = float(target_map[t.key])
            hits = scan_candidates(Cs=[C], m_values=m_values, target_G=tgt)
            best = hits[0]
            if args.max_rel_err is not None and abs(best.rel_err) > float(args.max_rel_err):
                ok = False
                break
            inv_alpha = float(best.G)
            alpha = 1.0 / inv_alpha if inv_alpha != 0 else float("inf")
            lam = lambda_qcd_from_alpha_s(alpha_s_mu=alpha, mu_GeV=float(t.Q_GeV), n_f=suite.n_f, loops=suite.loops)
            per.append((t.key, float(best.rel_err), float(lam.Lambda_GeV)))
            lambdas.append(float(lam.Lambda_GeV))
        if not ok:
            continue

        lam_min = min(lambdas)
        lam_max = max(lambdas)
        lam_mean = sum(lambdas) / float(len(lambdas))
        spread = (lam_max - lam_min) / lam_mean if lam_mean != 0 else float("inf")

        rows.append((spread, abs(best_a.rel_err), C, int(best_a.m), lam_a.Lambda_GeV, per))

    rows.sort(key=lambda r: (r[0], r[1]))

    if not rows:
        print("No candidates survived filters.")
        return 0

    print(f"Anchor: {suite.anchor.key} at Q={suite.anchor.Q_GeV:g} GeV")
    print(f"  target={anchor_target:.12g}")
    if args.max_rel_err is not None:
        print(f"  filter: |rel_err| <= {float(args.max_rel_err):g}\n")
    else:
        print("")

    print(f"Top {min(args.top, len(rows))} candidates by Lambda spread:")
    for i, (spread, aerr, C, m, lam_anchor, per) in enumerate(rows[: int(args.top)]):
        lab = label_by_C.get(C, "")
        print(f"- #{i+1:02d} {lab:22s} C={C:g}, m_anchor={m:+d}, anchor_rel_err={aerr:.3e}")
        print(f"       Lambda(anchor)={lam_anchor:.6g} GeV   spread={(spread*100.0):.3g}%")
        for key, rel_err, lam in per:
            print(f"       {key:28s} rel_err={rel_err:.3e}  Lambda={lam:.6g} GeV")
    return 0


def cmd_oos_report(args: argparse.Namespace) -> int:
    """
    Out-of-sample report: evaluate a frozen list of targets against strict gauge-derived C
    under integer m and report pass/fail.
    """

    suites = oos_suites()
    if args.suite not in suites:
        raise SystemExit(f"Unknown OOS suite {args.suite!r}. Options: {', '.join(sorted(suites.keys()))}")
    oos = suites[args.suite]
    targets = resolve_oos_targets(oos)
    all_targets = {t.name: t for t in known_targets()}

    # Integer m grid
    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))

    # Gauge-derived C candidates
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    cand: list[tuple[str, float]] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for k, v in cs.items():
            cand.append((f"{g.name}:{k}", float(v)))

    # De-duplicate Cs (keep first label)
    seen: set[float] = set()
    Cs: list[float] = []
    label_by_C: dict[float, str] = {}
    for lab, c in cand:
        if c in seen:
            continue
        seen.add(c)
        Cs.append(c)
        label_by_C[c] = lab

    print(f"Out-of-sample target suite: {args.suite}")
    print(f"tol(|rel_err|) = {args.max_rel_err}")
    print(f"Gauge-derived Cs (unique) = {len(Cs)} from base={args.base}")
    print(f"include = {','.join(include)}")
    print(f"m range = [{min(m_values)}, {max(m_values)}]\n")

    n_pass = 0
    n_sigma = 0
    chi2 = 0.0
    for ot, t in zip(oos, targets):
        tgt = all_targets[t.name]
        hits = scan_candidates(Cs=Cs, m_values=m_values, target_G=tgt.value)
        best = hits[0]
        ok = abs(best.rel_err) <= args.max_rel_err
        status = "PASS" if ok else "FAIL"
        if ok:
            n_pass += 1
        lab = label_by_C.get(best.C, "")
        if tgt.sigma is not None and tgt.sigma > 0:
            z = best.abs_err / float(tgt.sigma)
            n_sigma += 1
            chi2 += float(z * z)
            z_s = f"  z={z:.3g}"
        else:
            z_s = ""
        print(
            f"[{status}] {t.name:28s} target={tgt.value:.12g}  "
            f"best: {lab:22s} C={best.C:g}, m={int(best.m):d}, G={best.G:.12g}, rel_err={best.rel_err:.3e}{z_s}"
        )
        print(f"       rationale: {ot.rationale}")
    if n_sigma:
        print(f"\nSummary: {n_pass}/{len(oos)} PASS at tol={args.max_rel_err}  (sigma-annotated: n={n_sigma}, chi2={chi2:.6g})")
    else:
        print(f"\nSummary: {n_pass}/{len(oos)} PASS at tol={args.max_rel_err}")
    return 0


def cmd_oos_predictive(args: argparse.Namespace) -> int:
    """
    Predictive OOS report:
      1) For each force, fit (C,m) for a strict anchor target using the strict gauge-derived C menu.
      2) Freeze C for that force.
      3) Evaluate additional targets for that force with C held fixed (only m can vary).
    """

    suites = predictive_force_suites()
    if args.suite not in suites:
        raise SystemExit(f"Unknown predictive OOS suite {args.suite!r}. Options: {', '.join(sorted(suites.keys()))}")

    anchors, by_force = suites[args.suite]
    all_targets = {t.name: t for t in known_targets()}

    # Integer m grid
    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))

    # Gauge-derived C candidates
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    cand: list[tuple[str, float]] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for k, v in cs.items():
            cand.append((f"{g.name}:{k}", float(v)))

    # De-duplicate Cs (keep first label)
    seen: set[float] = set()
    Cs: list[float] = []
    label_by_C: dict[float, str] = {}
    for lab, c in cand:
        if c in seen:
            continue
        seen.add(c)
        Cs.append(c)
        label_by_C[c] = lab

    # Which forces to run
    if args.force == "all":
        forces = ["em", "strong", "weak"]
        if "hyper" in anchors:
            forces.append("hyper")
        if "gravity" in anchors:
            forces.append("gravity")
    else:
        forces = [args.force]

    print(f"Predictive OOS suite: {args.suite}")
    print(f"tol(|rel_err|) = {args.max_rel_err}")
    print(f"Gauge-derived Cs (unique) = {len(Cs)} from base={args.base}")
    print(f"include = {','.join(include)}")
    print(f"m range = [{min(m_values)}, {max(m_values)}]")
    print(f"norm_family = {args.norm_family}")
    print("Rule: fit anchor with free C, then HOLD C fixed per force.\n")

    total_pass = 0
    total_n = 0
    total_n_sigma = 0
    total_chi2 = 0.0

    for force in forces:
        if force not in anchors:
            raise SystemExit(f"Unknown force {force!r}. Options: {', '.join(sorted(anchors.keys()))}")
        if force not in by_force:
            raise SystemExit(f"Predictive suite {args.suite!r} missing targets for force {force!r}")

        anchor = anchors[force]
        if anchor.key not in all_targets:
            raise SystemExit(f"Unknown anchor key {anchor.key!r} for force {force!r}. Run `list-targets`.")

        factor = normalization_factor_for_force(force, family=args.norm_family)

        # Fit the anchor with the strict gauge-derived menu (choose best C and m)
        anchor_target = all_targets[anchor.key].value * factor
        anchor_hits = scan_candidates(Cs=Cs, m_values=m_values, target_G=anchor_target)
        best_anchor = anchor_hits[0]
        C0 = best_anchor.C
        m0 = int(best_anchor.m)
        lab = label_by_C.get(C0, "")

        print(f"{force.upper():7s} anchor: {anchor.key:28s} target={all_targets[anchor.key].value:.12g}")
        if abs(factor - 1.0) > 1e-12:
            print(f"         norm: factor={factor:.12g}  target_norm={anchor_target:.12g}")
        if all_targets[anchor.key].sigma is not None:
            sigma_norm = abs(float(factor)) * float(all_targets[anchor.key].sigma)
        else:
            sigma_norm = None
        if sigma_norm is not None and sigma_norm > 0:
            z_a = (best_anchor.G - anchor_target) / sigma_norm
            z_s = f"  z={z_a:.3g}"
        else:
            z_s = ""
        print(f"         best: {lab:22s} C={C0:g}, m={m0:d}, G={best_anchor.G:.12g}, rel_err={best_anchor.rel_err:.3e}{z_s}")
        print(f"         note: {anchor.rationale}")

        n_pass = 0
        n = 0
        n_sigma = 0
        chi2 = 0.0
        if sigma_norm is not None and sigma_norm > 0:
            n_sigma += 1
            chi2 += float(z_a * z_a)
        for ot in by_force[force]:
            if ot.key not in all_targets:
                raise SystemExit(f"Unknown predictive target key {ot.key!r}. Run `list-targets`.")
            tgt0 = all_targets[ot.key].value
            tgt = tgt0 * factor
            hits = scan_candidates(Cs=[C0], m_values=m_values, target_G=tgt)
            best = hits[0]
            ok = abs(best.rel_err) <= args.max_rel_err
            status = "PASS" if ok else "FAIL"
            if ok:
                n_pass += 1
            n += 1
            dm = int(best.m) - m0
            if all_targets[ot.key].sigma is not None:
                sigma_norm = abs(float(factor)) * float(all_targets[ot.key].sigma)
            else:
                sigma_norm = None
            if sigma_norm is not None and sigma_norm > 0:
                z = (best.G - tgt) / sigma_norm
                n_sigma += 1
                chi2 += float(z * z)
                z_t = f"  z={z:.3g}"
            else:
                z_t = ""
            print(
                f"  [{status}] {ot.key:28s} target={tgt0:.12g}  "
                f"m={int(best.m):d} (dm={dm:+d})  G={best.G:.12g}  rel_err={best.rel_err:.3e}{z_t}"
            )
            if abs(factor - 1.0) > 1e-12:
                print(f"         target_norm={tgt:.12g}")
            print(f"         rationale: {ot.rationale}")
        if n_sigma:
            print(f"  Force summary: {n_pass}/{n} PASS at tol={args.max_rel_err}  (sigma-annotated: n={n_sigma}, chi2={chi2:.6g})\n")
        else:
            print(f"  Force summary: {n_pass}/{n} PASS at tol={args.max_rel_err}\n")

        total_pass += n_pass
        total_n += n
        total_n_sigma += n_sigma
        total_chi2 += chi2

    if total_n_sigma:
        print(
            f"Overall predictive summary: {total_pass}/{total_n} PASS at tol={args.max_rel_err}  "
            f"(sigma-annotated: n={total_n_sigma}, chi2={total_chi2:.6g})"
        )
    else:
        print(f"Overall predictive summary: {total_pass}/{total_n} PASS at tol={args.max_rel_err}")
    return 0


def cmd_oos_predictive_rg(args: argparse.Namespace) -> int:
    """
    Predictive OOS with deterministic within-band RG running (strong + EM + weak):

      1) Fit a strict anchor on the φ-lattice using the strict gauge-derived C menu.
      2) Freeze (C, m_anchor) -> implies an anchor inverse coupling inv(Q0)=C/φ^m.
      3) Predict inv(Q) at additional scales using a deterministic RG prescription,
         with no additional fit (no re-choosing m per target).

    This implements the "m = major transition, running happens within the band" lever:
    integer m labels the anchor band, and RG supplies the continuous within-band motion.

    Note: EW within-band running is currently implemented only for alpha2^{-1}(Q) in suite v2.
    """

    suites = predictive_force_suites()
    if args.suite not in suites:
        raise SystemExit(f"Unknown predictive OOS suite {args.suite!r}. Options: {', '.join(sorted(suites.keys()))}")

    anchors, by_force = suites[args.suite]
    all_targets = {t.name: t for t in known_targets()}

    # Integer m grid
    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))

    # Gauge-derived C candidates
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    cand: list[tuple[str, float]] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for k, v in cs.items():
            cand.append((f"{g.name}:{k}", float(v)))

    # De-duplicate Cs (keep first label)
    seen: set[float] = set()
    Cs: list[float] = []
    label_by_C: dict[float, str] = {}
    for lab, c in cand:
        if c in seen:
            continue
        seen.add(c)
        Cs.append(c)
        label_by_C[c] = lab

    # Build rationale map so custom target lists still print useful context when possible.
    rationale_by_key: dict[str, str] = {}
    try:
        for _suite_key, lst in oos_suites().items():
            for ot in lst:
                rationale_by_key.setdefault(ot.key, ot.rationale)
    except Exception:
        pass
    for f in sorted(anchors.keys()):
        for ot in by_force.get(f, []):
            rationale_by_key.setdefault(ot.key, ot.rationale)

    # Which forces to run
    if args.force == "all":
        # Preserve the historical meaning of v1 ("strong + EM only"), then extend:
        #  - v2 adds weak (alpha2^{-1}) running targets
        #  - v3 adds hypercharge (alpha1_GUT^{-1}) running targets
        if args.suite == "v1":
            forces = ["em", "strong"]
        elif args.suite == "v2":
            forces = ["em", "strong", "weak"]
        else:
            forces = ["em", "strong", "weak", "hyper"]
    else:
        forces = [args.force]

    total_pass = 0
    total_n = 0
    total_n_sigma = 0
    total_chi2 = 0.0

    for force in forces:
        if force == "weak" and args.suite == "v1":
            raise SystemExit("Suite v1 does not define weak RG-within-band targets. Use --suite v2 or v3 for --force weak.")
        if force == "hyper" and args.suite in ("v1", "v2"):
            raise SystemExit("Suite v3 defines hypercharge RG-within-band targets. Use --suite v3 for --force hyper.")
        if force not in anchors or force not in by_force:
            raise SystemExit(f"Unknown force {force!r} for suite {args.suite!r}. Options: {', '.join(sorted(anchors.keys()))}")

        anchor = anchors[force]
        if anchor.key not in all_targets:
            raise SystemExit(f"Unknown anchor key {anchor.key!r}. Run `list-targets`.")

        # Determine reference scale Q0
        if force == "strong":
            if getattr(args, "Q0_GeV", None) is not None:
                raise SystemExit("--Q0-GeV is not supported for strong (Q0 is defined by the anchor key scale).")
            Q0 = scale_GeV_from_target_key(anchor.key)
        elif force == "em":
            # α(0) is quoted at Q→0; we use m_e as a fixed reference scale for the 1-loop threshold model.
            Q0 = float(args.Q0_GeV) if getattr(args, "Q0_GeV", None) is not None else 0.00051099895
        elif force == "weak":
            # Weak anchor is defined at the Z scale; use mZ as the reference.
            Q0 = float(get_measurement("mZ_GeV", default_value=91.1876).value)
        elif force == "hyper":
            # Hypercharge (GUT-normalized) anchor is defined at the Z scale; use mZ as the reference.
            Q0 = float(get_measurement("mZ_GeV", default_value=91.1876).value)
        else:
            raise SystemExit("oos-predictive-rg currently supports only --force em|strong|weak|hyper|all")

        # Fit the anchor on the lattice (choose best C and integer m)
        anchor_target = float(all_targets[anchor.key].value)
        anchor_hits = scan_candidates(Cs=Cs, m_values=m_values, target_G=anchor_target)
        best_anchor = anchor_hits[0]
        C0 = float(best_anchor.C)
        m0 = int(best_anchor.m)
        inv0 = float(best_anchor.G)
        lab = label_by_C.get(C0, "")
        if all_targets[anchor.key].sigma is not None and all_targets[anchor.key].sigma > 0:
            z_anchor = (inv0 - anchor_target) / float(all_targets[anchor.key].sigma)
            z_anchor_s = f"  z={z_anchor:.3g}"
        else:
            z_anchor = None
            z_anchor_s = ""

        # Choose running spec (deterministic; no tuning knobs)
        runner = str(args.runner).strip().lower()
        spec_qcd: QCDRunSpec | None = None
        if force == "strong":
            if runner == "auto":
                spec_qcd = qcd_run_spec_from_key(anchor.key)
            elif runner == "1loop_nf5":
                spec_qcd = QCDRunSpec(loops=1, nf_mode="const", n_f=5)
            elif runner == "1loop_nf56":
                spec_qcd = QCDRunSpec(loops=1, nf_mode="nf56")
            elif runner == "2loop_nf5":
                spec_qcd = QCDRunSpec(loops=2, nf_mode="const", n_f=5, steps_per_unit_log=int(args.steps_per_unit_log))
            elif runner == "2loop_nf56":
                spec_qcd = QCDRunSpec(loops=2, nf_mode="nf56", steps_per_unit_log=int(args.steps_per_unit_log))
            elif runner == "qed_1loop":
                raise SystemExit("--runner qed_1loop is for EM only")
            elif runner == "ew_sm_1loop":
                raise SystemExit("--runner ew_sm_1loop is for weak/hyper only")
            else:
                raise SystemExit(f"Unknown --runner {args.runner!r}")

            # Allow overriding the (still principled) threshold switch scale used by nf56 variants
            spec_qcd = QCDRunSpec(
                loops=spec_qcd.loops,
                nf_mode=spec_qcd.nf_mode,
                n_f=spec_qcd.n_f,
                Q_switch_GeV=float(args.Q_switch_GeV),
                n_f_below=spec_qcd.n_f_below,
                n_f_above=spec_qcd.n_f_above,
                steps_per_unit_log=spec_qcd.steps_per_unit_log,
            )
        elif force == "em":
            # EM runners:
            #  - qed_1loop: simple 1-loop threshold model
            #  - qed_pdg_mZ: PDG-style Δα decomposition at mZ (uses external Δα inputs)
            #  - auto: defaults to qed_pdg_mZ
            em_runner = "qed_pdg_mZ" if runner == "auto" else runner
            if em_runner not in ("qed_1loop", "qed_pdg_mZ"):
                raise SystemExit("For EM, use --runner auto, --runner qed_pdg_mZ, or --runner qed_1loop")
        else:
            # EW: SM 1-loop running of inverse couplings (fixed beta coefficient).
            ew_runner = "ew_sm_1loop" if runner == "auto" else runner
            if ew_runner != "ew_sm_1loop":
                raise SystemExit("For weak/hyper, use --runner auto or --runner ew_sm_1loop")

        # Pretty header
        print(f"Predictive RG-within-band OOS suite: {args.suite}")
        print(f"force = {force}")
        print(f"tol(|rel_err|) = {args.max_rel_err}")
        print(f"Gauge-derived Cs (unique) = {len(Cs)} from base={args.base}")
        print(f"include = {','.join(include)}")
        print(f"m range = [{min(m_values)}, {max(m_values)}]")
        if force == "strong":
            assert spec_qcd is not None
            print(f"runner = {args.runner}  -> loops={spec_qcd.loops}, nf_mode={spec_qcd.nf_mode}, n_f={spec_qcd.n_f}, Q0={Q0:g} GeV\n")
            print(f"STRONG  anchor: {anchor.key:28s} target={anchor_target:.12g}  Q0={Q0:g} GeV")
        elif force == "em":
            if em_runner == "qed_1loop":
                runner_desc = "QED 1-loop (fermion thresholds)"
            else:
                runner_desc = "PDG-style Δα(mZ^2)=Δα_lept+Δα_had^(5)+Δα_top"
            print(f"runner = {args.runner}  -> {runner_desc}, Q0={Q0:g} GeV\n")
            print(f"EM      anchor: {anchor.key:28s} target={anchor_target:.12g}  Q0={Q0:g} GeV (fixed)")
        elif force == "weak":
            print(f"runner = {args.runner}  -> SM 1-loop EW running (alpha2^{-1}), Q0={Q0:g} GeV\n")
            print(f"WEAK    anchor: {anchor.key:28s} target={anchor_target:.12g}  Q0={Q0:g} GeV")
        else:
            print(f"runner = {args.runner}  -> SM 1-loop EW running (alpha1_GUT^{-1}), Q0={Q0:g} GeV\n")
            print(f"HYPER   anchor: {anchor.key:28s} target={anchor_target:.12g}  Q0={Q0:g} GeV")
        print(f"         best: {lab:22s} C={C0:g}, m={m0:d}, inv0={inv0:.12g}, rel_err={best_anchor.rel_err:.3e}{z_anchor_s}")
        print(f"         note: {anchor.rationale}")

        # Precompute phi/log(phi)
        p = (1.0 + math.sqrt(5.0)) / 2.0
        ln_phi = math.log(p)

        # Target list
        if force == "strong" and str(getattr(args, "targets", "")).strip():
            keys = [s.strip() for s in str(args.targets).split(",") if s.strip()]
            targets_to_eval: list[tuple[str, str]] = [(k, rationale_by_key.get(k, "(custom target list)")) for k in keys]
        else:
            targets_to_eval = [(ot.key, ot.rationale) for ot in by_force[force]]

        n_pass = 0
        n = 0
        n_sigma = 0
        chi2 = 0.0
        if z_anchor is not None:
            n_sigma += 1
            chi2 += float(z_anchor * z_anchor)
        for key, rationale in targets_to_eval:
            if key not in all_targets:
                raise SystemExit(f"Unknown target key {key!r}. Run `list-targets`.")
            tgt = float(all_targets[key].value)
            Q = scale_GeV_from_target_key(key)

            if force == "strong":
                assert spec_qcd is not None
                alpha0 = (1.0 / inv0) if inv0 != 0 else float("inf")
                aQ = alpha_s_from_ref(Q, alpha_s_Q0=alpha0, Q0_GeV=Q0, spec=spec_qcd)
                inv_pred = (1.0 / aQ) if aQ not in (0.0, float("inf")) else float("inf")
            elif force == "em":
                if em_runner == "qed_1loop":
                    inv_pred = qed_run_alpha_inv_1loop_from_ref(Q, alpha_inv_Q0=inv0, Q0_GeV=Q0)
                else:
                    # PDG-style Δα inputs are defined at mZ only (by construction).
                    mZ = get_measurement("mZ_GeV", default_value=91.1876).value
                    if abs(float(Q) - float(mZ)) > 1e-3:
                        raise SystemExit(f"runner qed_pdg_mZ only supports Q≈mZ; got Q={Q:g} GeV from target {key!r}")
                    da_lept = get_measurement("delta_alpha_lept_mZ2", default_value=0.0314977).value
                    da_had5 = get_measurement("delta_alpha_had5_mZ2", default_value=0.02764).value
                    da_top = get_measurement("delta_alpha_top_mZ2", default_value=-0.00007).value
                    inv_pred = alpha_inv_mZ_from_delta_alpha(
                        alpha_inv_0=inv0, delta_alpha_lept=da_lept, delta_alpha_had5=da_had5, delta_alpha_top=da_top
                    )
            else:
                # EW: SM 1-loop running of inverse couplings.
                inv_pred = run_alpha_inv(inv0, Q0, Q, SM_1LOOP.b2) if force == "weak" else run_alpha_inv(inv0, Q0, Q, SM_1LOOP.b1)

            rel_err = (inv_pred - tgt) / tgt if tgt != 0 else float("nan")
            if all_targets[key].sigma is not None and all_targets[key].sigma > 0:
                z = (inv_pred - tgt) / float(all_targets[key].sigma)
                z_s = f"  z={z:.3g}"
                n_sigma += 1
                chi2 += float(z * z)
            else:
                z_s = ""

            ok = abs(rel_err) <= float(args.max_rel_err)
            status = "PASS" if ok else "FAIL"
            if ok:
                n_pass += 1
            n += 1

            # Interpret the RG-shift as a real-valued Δm inside the same C-band:
            #   inv(Q) = C / φ^(m_eff(Q))  => m_eff = log_phi(C/inv(Q))
            if inv_pred > 0 and C0 > 0 and ln_phi != 0:
                m_eff = math.log(C0 / inv_pred) / ln_phi
                dm_real = m_eff - float(m0)
                dm_int = int(round(dm_real))
                dm_frac = dm_real - float(dm_int)
            else:
                dm_real = float("nan")
                dm_int = 0
                dm_frac = float("nan")

            print(
                f"  [{status}] {key:28s} target={tgt:.12g}  Q={Q:g} GeV  "
                f"pred={inv_pred:.12g}  rel_err={rel_err:.3e}  "
                f"dm_real={dm_real:.6f}  dm_int={dm_int:+d}  frac={dm_frac:+.6f}{z_s}"
            )
            print(f"         rationale: {rationale}")

        if n_sigma:
            print(f"\nForce summary: {n_pass}/{n} PASS at tol={args.max_rel_err}  (sigma-annotated: n={n_sigma}, chi2={chi2:.6g})\n")
        else:
            print(f"\nForce summary: {n_pass}/{n} PASS at tol={args.max_rel_err}\n")

        total_pass += n_pass
        total_n += n
        total_n_sigma += n_sigma
        total_chi2 += chi2

    if len(forces) > 1:
        if total_n_sigma:
            print(
                f"Overall RG-within-band predictive summary: {total_pass}/{total_n} PASS at tol={args.max_rel_err}  "
                f"(sigma-annotated: n={total_n_sigma}, chi2={total_chi2:.6g})"
            )
        else:
            print(f"Overall RG-within-band predictive summary: {total_pass}/{total_n} PASS at tol={args.max_rel_err}")
    return 0


def cmd_oos_ew_mix(args: argparse.Namespace) -> int:
    """
    Electroweak mixing cross-check (derived sin^2θW(Q) from α2 and α1_GUT):

      1) Fit lattice anchors for α2^{-1}(mZ) and α1_GUT^{-1}(mZ) independently.
      2) Run each inverse coupling with SM 1-loop running to other scales Q.
      3) Form a derived mixing angle:

            sin^2θW(Q) := αY(Q) / (α2(Q) + αY(Q)),

         where αY = (3/5) α1_GUT and inv(αY) = (5/3) inv(α1_GUT).

    This is not a substitute for a full scheme-aware EW analysis; it is a deterministic
    internal-consistency probe in the same “band + within-band RG” spirit.
    """

    all_targets = {t.name: t for t in known_targets()}

    key_a2 = "1/alpha2(alpha(mZ),sin2_on_shell)"
    key_a1 = "1/alpha1_GUT(alpha(mZ),sin2_on_shell)"
    key_sin2_os = "sin2thetaW(on-shell)"

    for k in (key_a2, key_a1, key_sin2_os):
        if k not in all_targets:
            raise SystemExit(f"Missing required target {k!r}. Run `list-targets`.")

    inv_a2_target0 = float(all_targets[key_a2].value)
    inv_a1_target0 = float(all_targets[key_a1].value)
    sin2_os_target = float(all_targets[key_sin2_os].value)

    # Reference scale (mZ)
    Q0 = float(get_measurement("mZ_GeV", default_value=91.1876).value)

    # Integer m grid
    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))

    # Gauge-derived C candidates
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    cand: list[tuple[str, float]] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for k, v in cs.items():
            cand.append((f"{g.name}:{k}", float(v)))

    # De-duplicate Cs (keep first label)
    seen: set[float] = set()
    Cs: list[float] = []
    label_by_C: dict[float, str] = {}
    for lab, c in cand:
        if c in seen:
            continue
        seen.add(c)
        Cs.append(c)
        label_by_C[c] = lab

    # Fit the two anchors independently (no shared constraints yet)
    best_a2 = scan_candidates(Cs=Cs, m_values=m_values, target_G=inv_a2_target0)[0]
    best_a1 = scan_candidates(Cs=Cs, m_values=m_values, target_G=inv_a1_target0)[0]

    inv_a2_0 = float(best_a2.G)
    inv_a1_0 = float(best_a1.G)

    inv_aY_0 = (5.0 / 3.0) * inv_a1_0
    sin2_pred_mZ = inv_a2_0 / (inv_a2_0 + inv_aY_0) if (inv_a2_0 + inv_aY_0) != 0 else float("nan")
    rel_err_mZ = (sin2_pred_mZ - sin2_os_target) / sin2_os_target if sin2_os_target != 0 else float("nan")

    print("EW mixing (derived) OOS check")
    print(f"tol(|rel_err|) = {args.max_rel_err}")
    print(f"Gauge-derived Cs (unique) = {len(Cs)} from base={args.base}")
    print(f"include = {','.join(include)}")
    print(f"m range = [{min(m_values)}, {max(m_values)}]")
    print(f"Q0 = {Q0:g} GeV (mZ)\n")

    print("Anchor fits (independent):")
    print(
        f"  alpha2^-1 @ mZ: target={inv_a2_target0:.12g}  "
        f"best: {label_by_C.get(best_a2.C, ''):22s} C={best_a2.C:g}, m={int(best_a2.m):d}, inv0={inv_a2_0:.12g}, rel_err={best_a2.rel_err:.3e}"
    )
    print(
        f"  alpha1_GUT^-1 @ mZ: target={inv_a1_target0:.12g}  "
        f"best: {label_by_C.get(best_a1.C, ''):22s} C={best_a1.C:g}, m={int(best_a1.m):d}, inv0={inv_a1_0:.12g}, rel_err={best_a1.rel_err:.3e}"
    )
    print(f"  derived sin2thetaW(on-shell) @ mZ: target={sin2_os_target:.12g}  pred={sin2_pred_mZ:.12g}  rel_err={rel_err_mZ:.3e}\n")

    # Scales to evaluate
    from physics_test.scales import scale_GeV  # local import to keep CLI imports stable

    scale_labels = [s.strip() for s in str(args.scales).split(",") if s.strip()]
    if not scale_labels:
        raise SystemExit("--scales must contain at least one label (e.g. mW,mH,1TeV,10TeV)")

    n_pass = 0
    n = 0
    for lab in scale_labels:
        Q = float(scale_GeV(lab))
        if Q <= 0:
            raise SystemExit(f"Invalid scale label {lab!r} -> Q={Q:g} GeV (must be > 0)")

        inv2_pred = run_alpha_inv(inv_a2_0, Q0, Q, SM_1LOOP.b2)
        inv1_pred = run_alpha_inv(inv_a1_0, Q0, Q, SM_1LOOP.b1)
        invY_pred = (5.0 / 3.0) * float(inv1_pred)
        sin2_pred = float(inv2_pred) / (float(inv2_pred) + invY_pred) if (float(inv2_pred) + invY_pred) != 0 else float("nan")

        inv2_tgt = run_alpha_inv(inv_a2_target0, Q0, Q, SM_1LOOP.b2)
        inv1_tgt = run_alpha_inv(inv_a1_target0, Q0, Q, SM_1LOOP.b1)
        invY_tgt = (5.0 / 3.0) * float(inv1_tgt)
        sin2_tgt = float(inv2_tgt) / (float(inv2_tgt) + invY_tgt) if (float(inv2_tgt) + invY_tgt) != 0 else float("nan")

        rel_err = (sin2_pred - sin2_tgt) / sin2_tgt if sin2_tgt != 0 else float("nan")
        ok = abs(rel_err) <= float(args.max_rel_err)
        status = "PASS" if ok else "FAIL"
        if ok:
            n_pass += 1
        n += 1

        print(f"[{status}] Q={Q:g} GeV ({lab})  sin2_target={sin2_tgt:.12g}  sin2_pred={sin2_pred:.12g}  rel_err={rel_err:.3e}")

    print(f"\nSummary: {n_pass}/{n} PASS at tol={args.max_rel_err}")
    return 0


def cmd_ew_sin2(args: argparse.Namespace) -> int:
    """
    Predict sin^2θW(Q) from lattice-quantized EW anchors + SM/MSSM 1-loop running.

    This is a convenience wrapper around the same logic used in `oos-ew-mix`, but:
      - it prints predictions at user-chosen scales, and
      - it can optionally compare against user-supplied measurements.
    """

    all_targets = {t.name: t for t in known_targets()}

    key_a2 = "1/alpha2(alpha(mZ),sin2_on_shell)"
    key_a1 = "1/alpha1_GUT(alpha(mZ),sin2_on_shell)"

    for k in (key_a2, key_a1):
        if k not in all_targets:
            raise SystemExit(f"Missing required target {k!r}. Run `list-targets`.")

    inv_a2_target0 = float(all_targets[key_a2].value)
    inv_a1_target0 = float(all_targets[key_a1].value)

    # Reference scale (mZ)
    Q0 = float(get_measurement("mZ_GeV", default_value=91.1876).value)

    # Integer m grid
    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))

    # Gauge-derived C candidates
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    cand: list[tuple[str, float]] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for k, v in cs.items():
            cand.append((f"{g.name}:{k}", float(v)))

    # De-duplicate Cs (keep first label)
    seen: set[float] = set()
    Cs: list[float] = []
    label_by_C: dict[float, str] = {}
    for lab, c in cand:
        if c in seen:
            continue
        seen.add(c)
        Cs.append(c)
        label_by_C[c] = lab

    # Fit lattice anchors (independent)
    best_a2 = scan_candidates(Cs=Cs, m_values=m_values, target_G=inv_a2_target0)[0]
    best_a1 = scan_candidates(Cs=Cs, m_values=m_values, target_G=inv_a1_target0)[0]
    inv_a2_0 = float(best_a2.G)
    inv_a1_0 = float(best_a1.G)

    # Choose beta set for running
    betas = SM_1LOOP if str(args.model).strip().lower() == "sm" else MSSM_1LOOP
    method = str(getattr(args, "method", "gammaZ_1loop")).strip()

    # Parse measurement overrides: label -> (value, sigma|None)
    meas: dict[str, tuple[float, float | None]] = {}
    for trip in args.measurement or []:
        if len(trip) != 3:
            raise SystemExit("--measurement expects: <label> <sin2> <sigma_or_0>")
        lab, v_s, s_s = trip
        v = float(v_s)
        s = float(s_s)
        meas[str(lab)] = (v, None if s <= 0 else s)

    # Scales to evaluate
    from physics_test.scales import scale_GeV  # local import to keep CLI imports stable

    scale_labels = [s.strip() for s in str(args.scales).split(",") if s.strip()]
    if not scale_labels:
        raise SystemExit("--scales must contain at least one label (e.g. mW,mH,1TeV,10TeV)")

    print("EW sin^2thetaW(Q) prediction (lattice anchors + 1-loop running)")
    print(f"model = {betas.name}")
    print(f"method = {method}")
    print(f"Q0 = {Q0:g} GeV (mZ)")
    print(f"Gauge-derived Cs (unique) = {len(Cs)} from base={args.base}")
    print(f"include = {','.join(include)}")
    print(f"m range = [{min(m_values)}, {max(m_values)}]\n")

    print("Anchor fits (independent):")
    print(
        f"  inv_alpha2(mZ): target={inv_a2_target0:.12g}  "
        f"best: {label_by_C.get(float(best_a2.C), ''):22s} C={best_a2.C:g}, m={int(best_a2.m):d}, inv0={inv_a2_0:.12g}, rel_err={best_a2.rel_err:.3e}"
    )
    print(
        f"  inv_alpha1_GUT(mZ): target={inv_a1_target0:.12g}  "
        f"best: {label_by_C.get(float(best_a1.C), ''):22s} C={best_a1.C:g}, m={int(best_a1.m):d}, inv0={inv_a1_0:.12g}, rel_err={best_a1.rel_err:.3e}\n"
    )

    # Predict at scales
    for lab in scale_labels:
        Q = float(scale_GeV(lab))
        if method == "gauge_1loop":
            inv2 = run_alpha_inv(inv_a2_0, Q0, Q, betas.b2)
            inv1 = run_alpha_inv(inv_a1_0, Q0, Q, betas.b1)
            invY = (5.0 / 3.0) * float(inv1)
            sin2_pred = float(inv2) / (float(inv2) + invY) if (float(inv2) + invY) != 0 else float("nan")
        elif method == "gammaZ_1loop":
            invY0 = (5.0 / 3.0) * float(inv_a1_0)
            sin2_ref = float(inv_a2_0) / (float(inv_a2_0) + float(invY0)) if (float(inv_a2_0) + float(invY0)) != 0 else float("nan")
            inv_alpha_mZ = float(all_targets["1/alpha(mZ)"].value) if "1/alpha(mZ)" in all_targets else 127.955
            alpha_ref = (1.0 / inv_alpha_mZ) if inv_alpha_mZ != 0 else float("nan")
            sin2_pred = sin2_eff_gammaZ_1loop_from_ref(Q, sin2_ref=sin2_ref, alpha_ref=alpha_ref, Q0_GeV=Q0)
        else:
            raise SystemExit(
                "Unknown --method. Use --method gauge_1loop, --method gammaZ_1loop, or --method zpole_kappa_approx"
            )

        if lab in meas:
            v, s = meas[lab]
            rel_err = (sin2_pred - v) / v if v != 0 else float("nan")
            if s is not None and s > 0:
                z = (sin2_pred - v) / s
                z_s = f"  z={z:.3g}"
            else:
                z_s = ""
            ok = abs(rel_err) <= float(args.max_rel_err)
            status = "PASS" if ok else "FAIL"
            print(f"[{status}] Q={Q:g} GeV ({lab})  meas={v:.12g}  pred={sin2_pred:.12g}  rel_err={rel_err:.3e}{z_s}")
        else:
            print(f"[PRED] Q={Q:g} GeV ({lab})  sin2_pred={sin2_pred:.12g}")

    return 0


def cmd_oos_ew_sin2(args: argparse.Namespace) -> int:
    """
    OOS check for externally provided sin^2θW(Q) targets.

    This command looks for registry-backed targets whose names start with "sin2thetaW("
    and whose notes start with "Registry target:" (i.e., added via keys like
    "tgt_sin2thetaW(Qweak)" in the measurement registry).

    For each such target, it predicts sin^2θW at the target's Q using lattice-quantized
    α2^{-1}(mZ) and α1_GUT^{-1}(mZ) anchors + SM/MSSM 1-loop running, then compares.
    """

    all_targets = {t.name: t for t in known_targets()}

    suite = str(getattr(args, "suite", "registry-all")).strip()
    scheme_prefix = str(getattr(args, "scheme_prefix", "") or "").strip()
    method = str(getattr(args, "method", "gammaZ_1loop")).strip()
    if (suite.startswith("ew-independent-v") or suite == "ew-exploratory-v1") and not scheme_prefix:
        scheme_prefix = "sin2thetaW_eff:"
    if suite == "ew-dis-exploratory-v1" and not scheme_prefix:
        scheme_prefix = "sin2thetaW_on_shell:"
    if suite == "ew-zpole-exploratory-v1" and not scheme_prefix:
        scheme_prefix = "sin2thetaW_eff_lept:"
    if suite == "ew-zpole-exploratory-v1" and method == "gammaZ_1loop":
        # Default the Z-pole suite to its dedicated mapping method.
        method = "zpole_kappa_approx"

    # Collect external sin2 targets:
    #  - suite="registry-all": auto-discover all registry-backed sin2thetaW(Q) targets
    #  - suite=<frozen>: evaluate an explicit, frozen list (does not change as registry grows)
    ext: list[tuple[str, float, float | None, float, str, str, str]] = []  # (name, v, sigma, Q, scheme, citation, rationale)
    if suite == "registry-all":
        for t in all_targets.values():
            if not t.name.startswith("sin2thetaW("):
                continue
            if not str(t.note).startswith("Registry target:"):
                continue
            if t.Q_GeV is None:
                continue
            scheme = str(getattr(t, "scheme", "") or "")
            # Scheme isolation: when evaluating the auto-discovered registry menu, treat
            # --scheme-prefix as a *filter* (skip non-matching targets) rather than an error.
            if scheme_prefix and not scheme.startswith(scheme_prefix):
                continue
            ext.append(
                (
                    t.name,
                    float(t.value),
                    t.sigma,
                    float(t.Q_GeV),
                    scheme,
                    str(getattr(t, "citation", "") or ""),
                    "Registry target (auto-discovered)",
                )
            )
        if not ext:
            if scheme_prefix:
                raise SystemExit(
                    f"No external sin2thetaW(Q) targets matched --scheme-prefix {scheme_prefix!r}. "
                    "Either widen the prefix, or omit --scheme-prefix to list all registry targets."
                )
            raise SystemExit(
                "No external sin2thetaW(Q) targets found. Add registry keys like "
                "'tgt_sin2thetaW(Qweak)' with fields value/sigma/Q_GeV/scheme/citation "
                "(you can use PHYSICS_TEST_TARGET_REGISTRY to point to a local registry file)."
            )
    else:
        suites = ew_sin2_suites()
        if suite not in suites:
            raise SystemExit(f"Unknown --suite {suite!r}. Options: registry-all, {', '.join(sorted(suites.keys()))}")
        if (
            (suite.startswith("ew-independent-v") or suite == "ew-exploratory-v1" or suite == "ew-dis-exploratory-v1")
            and method != "gammaZ_1loop"
        ):
            raise SystemExit(
                "Suites ew-independent-v*, ew-exploratory-v1, and ew-dis-exploratory-v1 are frozen to "
                "--method gammaZ_1loop (low-Q weak-angle toy model)."
            )
        if suite == "ew-zpole-exploratory-v1" and method != "zpole_kappa_approx":
            raise SystemExit(
                "Suite ew-zpole-exploratory-v1 is frozen to --method zpole_kappa_approx "
                "(Z-pole effective leptonic weak mixing angle mapping)."
            )
        for ot in suites[suite]:
            if ot.key not in all_targets:
                raise SystemExit(f"Missing target {ot.key!r} for suite {suite!r}. Run `list-targets`.")
            t = all_targets[ot.key]
            if t.Q_GeV is None:
                raise SystemExit(f"Target {ot.key!r} is missing Q_GeV metadata (required for oos-ew-sin2).")
            scheme = str(getattr(t, "scheme", "") or "")
            if scheme_prefix and not scheme.startswith(scheme_prefix):
                raise SystemExit(
                    f"Suite {suite!r} requires scheme prefix {scheme_prefix!r} for target {ot.key!r}; got scheme={scheme!r}"
                )
            ext.append(
                (
                    t.name,
                    float(t.value),
                    t.sigma,
                    float(t.Q_GeV),
                    scheme,
                    str(getattr(t, "citation", "") or ""),
                    ot.rationale,
                )
            )

    key_a2 = "1/alpha2(alpha(mZ),sin2_on_shell)"
    key_a1 = "1/alpha1_GUT(alpha(mZ),sin2_on_shell)"
    for k in (key_a2, key_a1):
        if k not in all_targets:
            raise SystemExit(f"Missing required target {k!r}. Run `list-targets`.")

    inv_a2_target0 = float(all_targets[key_a2].value)
    inv_a1_target0 = float(all_targets[key_a1].value)

    # Reference scale (mZ)
    Q0 = float(get_measurement("mZ_GeV", default_value=91.1876).value)

    # Integer m grid
    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))

    # Gauge-derived C candidates
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    cand: list[tuple[str, float]] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for k, v in cs.items():
            cand.append((f"{g.name}:{k}", float(v)))

    # De-duplicate Cs (keep first label)
    seen: set[float] = set()
    Cs: list[float] = []
    label_by_C: dict[float, str] = {}
    for lab, c in cand:
        if c in seen:
            continue
        seen.add(c)
        Cs.append(c)
        label_by_C[c] = lab

    # Fit lattice anchors (independent)
    best_a2 = scan_candidates(Cs=Cs, m_values=m_values, target_G=inv_a2_target0)[0]
    best_a1 = scan_candidates(Cs=Cs, m_values=m_values, target_G=inv_a1_target0)[0]
    inv_a2_0 = float(best_a2.G)
    inv_a1_0 = float(best_a1.G)

    # Choose beta set for running
    betas = SM_1LOOP if str(args.model).strip().lower() == "sm" else MSSM_1LOOP

    z_max = getattr(args, "z_max", None)
    if suite == "ew-independent-v1" and (z_max is None or float(z_max) <= 0):
        z_max = 3.0
    if suite == "ew-independent-v2" and (z_max is None or float(z_max) <= 0):
        z_max = 2.0
    if suite == "ew-independent-v3" and (z_max is None or float(z_max) <= 0):
        z_max = 2.0
    if suite == "ew-exploratory-v1" and (z_max is None or float(z_max) <= 0):
        z_max = 2.0
    if suite == "ew-zpole-exploratory-v1" and (z_max is None or float(z_max) <= 0):
        z_max = 2.0
    sigma_theory = getattr(args, "sigma_theory", None)
    if (
        (suite.startswith("ew-independent-v") or suite == "ew-exploratory-v1")
        and (sigma_theory is None or float(sigma_theory) < 0)
    ):
        # Upgraded κ(Q) model: default to zero additional theory sigma for effective-angle suites.
        sigma_theory = 0.0
    sigma_theory_f = float(sigma_theory) if sigma_theory is not None else 0.0
    report_thresh = bool(getattr(args, "report_theory_threshold", False))

    print("OOS: external sin2thetaW(Q) targets (registry-driven)")
    print(f"model = {betas.name}")
    print(f"method = {method}")
    if z_max is not None:
        print(f"PASS rule = |z| <= {float(z_max):g} (fallback: |rel_err| <= {float(args.max_rel_err):g} if sigma missing)")
    else:
        print(f"tol(|rel_err|) = {args.max_rel_err}")
    print(f"suite = {suite}")
    if scheme_prefix:
        print(f"scheme requirement: scheme startswith {scheme_prefix!r}")
    if sigma_theory_f > 0:
        print(f"theory sigma: {sigma_theory_f:g} (combined in quadrature with experimental sigma)")
    print(f"Q0 = {Q0:g} GeV (mZ)")
    print(f"Gauge-derived Cs (unique) = {len(Cs)} from base={args.base}")
    print(f"include = {','.join(include)}")
    print(f"m range = [{min(m_values)}, {max(m_values)}]\n")

    print("Anchor fits (independent):")
    print(
        f"  inv_alpha2(mZ): target={inv_a2_target0:.12g}  "
        f"best: {label_by_C.get(float(best_a2.C), ''):22s} C={best_a2.C:g}, m={int(best_a2.m):d}, inv0={inv_a2_0:.12g}, rel_err={best_a2.rel_err:.3e}"
    )
    print(
        f"  inv_alpha1_GUT(mZ): target={inv_a1_target0:.12g}  "
        f"best: {label_by_C.get(float(best_a1.C), ''):22s} C={best_a1.C:g}, m={int(best_a1.m):d}, inv0={inv_a1_0:.12g}, rel_err={best_a1.rel_err:.3e}\n"
    )

    n_pass = 0
    n = 0
    n_sigma = 0
    chi2 = 0.0
    # For threshold reporting: per-target minimum sigma_theory needed to satisfy |z|<=z_max,
    # assuming sigma_eff = sqrt(sigma_exp^2 + sigma_theory^2).
    thresh_rows: list[tuple[str, float, float, float, float, float | None]] = []  # (name, delta, sigma_exp, z0, req_sigma_theory, Q)
    for name, v, s, Q, scheme, citation, rationale in sorted(ext, key=lambda x: x[3]):
        if method == "gauge_1loop":
            inv2 = run_alpha_inv(inv_a2_0, Q0, Q, betas.b2)
            inv1 = run_alpha_inv(inv_a1_0, Q0, Q, betas.b1)
            invY = (5.0 / 3.0) * float(inv1)
            pred = float(inv2) / (float(inv2) + invY) if (float(inv2) + invY) != 0 else float("nan")
        elif method == "gammaZ_1loop":
            # Build a lattice-defined reference angle at mZ from the fitted inverse couplings,
            # then apply a toy low-Q κ(Q) model for γ–Z mixing.
            invY0 = (5.0 / 3.0) * float(inv_a1_0)
            sin2_ref = float(inv_a2_0) / (float(inv_a2_0) + float(invY0)) if (float(inv_a2_0) + float(invY0)) != 0 else float("nan")
            # Use alpha(mZ) as a fixed coefficient (deterministic; registry-backed).
            inv_alpha_mZ = float(all_targets["1/alpha(mZ)"].value) if "1/alpha(mZ)" in all_targets else 127.955
            alpha_ref = (1.0 / inv_alpha_mZ) if inv_alpha_mZ != 0 else float("nan")
            pred = sin2_eff_gammaZ_1loop_from_ref(Q, sin2_ref=sin2_ref, alpha_ref=alpha_ref, Q0_GeV=Q0)
        elif method == "zpole_kappa_approx":
            # Z-pole effective leptonic weak mixing angle mapping (toy but deterministic):
            #
            #   sin^2θ_eff^lept(mZ) ≈ κ_Z * sin^2θ_ref(mZ),
            #
            # where:
            #   - sin^2θ_ref is the lattice-derived angle from the independently fitted α2^{-1}(mZ) and α1_GUT^{-1}(mZ),
            #   - κ_Z is an approximate net electroweak correction factor dominated by vacuum polarization (Δα)
            #     with a compensating leading top-loop term (Δρ_top).
            #
            # We model:
            #   κ_Z ≈ 1 + Δα_total(mZ^2) − Δρ_top
            #
            # This is not a substitute for a full EW precision calculation; it is a principled, auditable
            # mapping layer for Z-pole pseudo-observables.
            if suite == "registry-all" and not scheme_prefix:
                raise SystemExit("--method zpole_kappa_approx requires --scheme-prefix (e.g. sin2thetaW_eff_lept:).")
            if abs(float(Q) - float(Q0)) > 1e-3:
                raise SystemExit(
                    f"--method zpole_kappa_approx is only defined at the Z pole (Q≈mZ). Got Q={Q:g} GeV for {name}."
                )
            invY0 = (5.0 / 3.0) * float(inv_a1_0)
            sin2_ref = (
                float(inv_a2_0) / (float(inv_a2_0) + float(invY0)) if (float(inv_a2_0) + float(invY0)) != 0 else float("nan")
            )
            key_da = "delta_alpha_total(mZ2)"
            key_drho = "delta_rho_top(GF,mt)"
            if key_da not in all_targets or key_drho not in all_targets:
                raise SystemExit(
                    f"--method zpole_kappa_approx requires targets {key_da!r} and {key_drho!r}. Run `list-targets`."
                )
            da = float(all_targets[key_da].value)
            drho = float(all_targets[key_drho].value)
            kappa_Z = 1.0 + da - drho
            pred = float(kappa_Z) * float(sin2_ref)
        else:
            raise SystemExit("Unknown --method. Use --method gauge_1loop or --method gammaZ_1loop")

        rel_err = (pred - v) / v if v != 0 else float("nan")
        delta = float(pred - v)
        z = None
        sigma_eff = None
        if s is not None and s > 0:
            sigma_exp = float(s)
            sigma_eff = sigma_exp
            if sigma_theory_f > 0:
                sigma_eff = math.sqrt((sigma_eff * sigma_eff) + (sigma_theory_f * sigma_theory_f))
            z = (pred - v) / float(sigma_eff)
            if report_thresh:
                if z_max is None or float(z_max) <= 0:
                    raise SystemExit("--report-theory-threshold requires --z-max (or a suite with a default z-max)")
                z0 = delta / sigma_exp if sigma_exp != 0 else float("inf")
                # Need sigma_eff >= |delta|/z_max
                need = (abs(delta) / float(z_max)) if float(z_max) != 0 else float("inf")
                req2 = (need * need) - (sigma_exp * sigma_exp)
                req = math.sqrt(req2) if req2 > 0 else 0.0
                thresh_rows.append((name, delta, sigma_exp, z0, req, Q))

        if z_max is not None and z is not None:
            ok = abs(float(z)) <= float(z_max)
        else:
            ok = abs(rel_err) <= float(args.max_rel_err)
        status = "PASS" if ok else "FAIL"
        if ok:
            n_pass += 1
        n += 1

        if z is not None:
            n_sigma += 1
            chi2 += float(float(z) * float(z))
            if sigma_eff is not None and sigma_eff != float(s):
                z_s = f"  z={float(z):.3g}  sigma_eff={sigma_eff:.3g}"
            else:
                z_s = f"  z={float(z):.3g}"
        else:
            z_s = ""

        print(f"[{status}] {name:28s} Q={Q:g} GeV  meas={v:.12g}  pred={pred:.12g}  rel_err={rel_err:.3e}{z_s}")
        if rationale:
            print(f"         rationale: {rationale}")
        if scheme:
            print(f"         scheme: {scheme}")
        if citation:
            print(f"         citation: {citation}")

    if n_sigma:
        if z_max is not None:
            print(
                f"\nSummary: {n_pass}/{n} PASS  (sigma-annotated: n={n_sigma}, chi2={chi2:.6g}, chi2/ndf={chi2/float(n_sigma):.6g})"
            )
        else:
            print(
                f"\nSummary: {n_pass}/{n} PASS at tol={args.max_rel_err}  (sigma-annotated: n={n_sigma}, chi2={chi2:.6g})"
            )
    else:
        print(f"\nSummary: {n_pass}/{n} PASS at tol={args.max_rel_err}")

    if report_thresh and thresh_rows:
        # The minimal sigma_theory that makes all targets satisfy |z|<=z_max is the max required across targets.
        limiting = max(thresh_rows, key=lambda row: row[4])
        req_global = float(limiting[4])
        lim_name = limiting[0]
        print("\nTheory-sigma threshold report:")
        print(f"  z_max = {float(z_max):g}")
        print("  Per-target required sigma_theory (absolute on sin^2thetaW):")
        for name, delta, sigma_exp, z0, req, Q in sorted(thresh_rows, key=lambda row: row[4], reverse=True):
            print(
                f"   - {name:24s} Q={Q:g}  delta={delta:+.6g}  sigma_exp={sigma_exp:.6g}  z0={z0:+.3g}  req_sigma_theory={req:.6g}"
            )
        print(f"  => Minimal sigma_theory to pass all (at |z|<={float(z_max):g}) is {req_global:.6g}  (limiting target: {lim_name})")

    # CI gate behavior:
    #  - ew-independent-v1 is treated as a "one-command gate" by default: non-zero exit if any target FAILs.
    #  - other suites can opt into this behavior via --ci.
    if (suite.startswith("ew-independent-v") or bool(getattr(args, "ci", False))) and n_pass != n:
        return 2
    return 0


def cmd_oos_steps(args: argparse.Namespace) -> int:
    """
    Step-signal OOS report (C-independent):

    For each force, take a strict anchor and compute whether each additional target
    is consistent with an *integer* Δm step under the assumption of the same C, i.e.:

        (anchor/target) ≈ φ^(Δm),   with integer Δm.
    """

    suites = predictive_force_suites()
    if args.suite not in suites:
        raise SystemExit(f"Unknown step-suite {args.suite!r}. Options: {', '.join(sorted(suites.keys()))}")

    anchors, by_force = suites[args.suite]
    target_map = {t.name: t.value for t in known_targets()}

    # Which forces to run
    if args.force == "all":
        forces = ["em", "strong", "weak", "gravity"]
    else:
        forces = [args.force]

    tol = float(args.max_ratio_err)
    if tol < 0:
        raise SystemExit("--max-ratio-err must be non-negative")

    print(f"Step-signal suite: {args.suite}")
    print(f"tol(|ratio_err_if_int|) = {tol}")
    print("Rule: for each force, test whether (anchor/target) is close to φ^integer.\n")

    total_pass = 0
    total_n = 0

    for force in forces:
        if force not in anchors:
            raise SystemExit(f"Unknown force {force!r}. Options: {', '.join(sorted(anchors.keys()))}")
        if force not in by_force:
            raise SystemExit(f"Step-suite {args.suite!r} missing targets for force {force!r}")

        anchor = anchors[force]
        if anchor.key not in target_map:
            raise SystemExit(f"Unknown anchor key {anchor.key!r}. Run `list-targets`.")
        anchor_val = target_map[anchor.key]

        print(f"{force.upper():7s} anchor: {anchor.key:28s} value={anchor_val:.12g}")
        print(f"         note: {anchor.rationale}")

        n_pass = 0
        n = 0
        for ot in by_force[force]:
            if ot.key not in target_map:
                raise SystemExit(f"Unknown target key {ot.key!r}. Run `list-targets`.")
            tgt_val = target_map[ot.key]
            sr = step_from_targets(anchor_val, tgt_val)
            ok = sr.ratio_err_if_int <= tol
            status = "PASS" if ok else "FAIL"
            if ok:
                n_pass += 1
            n += 1
            print(
                f"  [{status}] {ot.key:28s} ratio={sr.ratio:.12g}  "
                f"dm_real={sr.dm_real:.6f}  dm_int={sr.dm_int:+d}  frac={sr.dm_frac:+.6f}  "
                f"ratio_err={sr.ratio_err_if_int:.3%}"
            )
            print(f"         rationale: {ot.rationale}")
        print(f"  Force summary: {n_pass}/{n} PASS at tol={tol}\n")

        total_pass += n_pass
        total_n += n

    print(f"Overall step-signal summary: {total_pass}/{total_n} PASS at tol={tol}")

    # Rough null baseline: assume dm_frac is uniformly distributed in [-0.5,0.5).
    try:
        ln_phi = math.log((1.0 + math.sqrt(5.0)) / 2.0)
        if tol >= 1.0:
            delta_thr = 0.5
        else:
            # asymmetric bounds for +δ/-δ; take the larger magnitude as the threshold
            d_pos = -math.log(max(1e-12, 1.0 - tol)) / ln_phi
            d_neg = math.log(1.0 + tol) / ln_phi
            delta_thr = min(0.5, max(d_pos, d_neg))
        p_null = min(1.0, 2.0 * delta_thr)

        # Binomial tail probability for >= total_pass successes (independence assumed; rough).
        # We keep this simple to avoid external deps.
        def binom_tail(n: int, k: int, p: float) -> float:
            if k <= 0:
                return 1.0
            if k > n:
                return 0.0
            out = 0.0
            for i in range(k, n + 1):
                out += math.comb(n, i) * (p**i) * ((1.0 - p) ** (n - i))
            return out

        p_tail = binom_tail(total_n, total_pass, p_null) if total_n > 0 else float("nan")
        print("\nNull baseline (rough): assume dm fractional part is uniform in [-0.5,0.5).")
        print(f"  dm_frac threshold ≈ {delta_thr:.6f}  => expected pass prob ≈ {p_null:.3%} per pair")
        print(f"  binomial P(X >= {total_pass} | n={total_n}, p={p_null:.3%}) ≈ {p_tail:.3g}")
    except Exception:
        # Baseline is optional; never fail the command because of it.
        pass
    return 0


def cmd_list_gravity_bands(args: argparse.Namespace) -> int:
    for b in gravity_band_list():
        print(f"{b.name:8s}  {b.f_min_hz:.3e} .. {b.f_max_hz:.3e} Hz  note={b.note}")
    return 0


def cmd_list_frequency_presets(args: argparse.Namespace) -> int:
    print("EM presets:")
    for p in em_frequency_presets():
        print(f"- {p.key:22s}  F0={p.F0_hz:.6e} Hz  note={p.note}")
    print("\nParticle/energy proxies (f=E/h):")
    for p in particle_proxy_presets():
        print(f"- {p.key:22s}  F0={p.F0_hz:.6e} Hz  note={p.note}")
    print("\nThermal anchors (kB*T/h):")
    for p in thermal_presets():
        print(f"- {p.key:28s}  F0={p.F0_hz:.6e} Hz  note={p.note}")
    return 0


def cmd_list_gauge_candidates(args: argparse.Namespace) -> int:
    groups = standard_model_gauge_groups()
    for g in groups:
        include = tuple(s.strip() for s in args.include.split(",") if s.strip())
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        print(f"{g.name}: rank={g.rank}, dim={g.dim}, h={g.coxeter_h}, h*={g.dual_coxeter_h}")
        for k, v in cs.items():
            print(f"  - {k:18s} C={v:.15g}")
    return 0


def cmd_scan_gauge_Cs(args: argparse.Namespace) -> int:
    targets = {t.name: t.value for t in known_targets()}
    if args.target not in targets:
        raise SystemExit(f"Unknown --target {args.target!r}. Run `python -m physics_test.cli list-targets`.")

    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))

    # Build candidate Cs from gauge groups and constructions
    Cs: list[float] = []
    labels: list[str] = []
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for k, v in cs.items():
            Cs.append(float(v))
            labels.append(f"{g.name}:{k}")

    # De-duplicate while keeping a label for the first occurrence
    seen: set[float] = set()
    uniq_Cs: list[float] = []
    uniq_labels: list[str] = []
    for c, lab in zip(Cs, labels):
        if c in seen:
            continue
        seen.add(c)
        uniq_Cs.append(c)
        uniq_labels.append(lab)

    hits = scan_candidates(Cs=uniq_Cs, m_values=m_values, target_G=targets[args.target])
    hits = filter_hits_by_rel_err(hits, max_abs_rel_err=args.max_rel_err)

    print(f"Target {args.target} = {targets[args.target]:.15g}")
    print(f"Gauge-derived Cs (unique) = {len(uniq_Cs)} from base={args.base}")
    print(f"include = {','.join(include)}")
    print(f"m range = [{min(m_values)}, {max(m_values)}] step={args.m_step} (integerized)")
    print(f"Kept: {len(hits)} hits with |rel_err| <= {args.max_rel_err}\n")

    # Build reverse map for labeling
    label_by_C = {c: lab for c, lab in zip(uniq_Cs, uniq_labels)}

    print("rank  label                 C           m     G=C/phi^m           rel_err")
    for i, h in enumerate(hits[: max(1, args.top)], start=1):
        lab = label_by_C.get(h.C, "")
        print(f"{i:4d}  {lab:20s}  {h.C:10.6g}  {int(h.m):4d}  {h.G:18.12g}  {h.rel_err: .6e}")
    return 0


def cmd_sweep_quantum_gravity(args: argparse.Namespace) -> int:
    """
    Sweep a gravity mass scale (in GeV) and check for solutions where:
      - C is restricted to gauge-derived candidates (from base + invariants)
      - m is integer
      - |rel_err| <= tolerance for EM/strong/weak/gravity coupling targets
      - gravity's predicted F0 (from fixed K) lies in a chosen GW band

    This is designed to be fast by restricting gravity m to the band-implied window.
    """
    targets = {t.name: t.value for t in known_targets()}
    for k in (args.em_target, args.strong_target, args.weak_target):
        if k not in targets:
            raise SystemExit(f"Unknown target {k!r}. Run `python -m physics_test.cli list-targets`.")

    # Gauge-derived C candidates
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    cand: list[tuple[str, float]] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for k, v in cs.items():
            cand.append((f"{g.name}:{k}", float(v)))
    seen: set[float] = set()
    Cs: list[float] = []
    label_by_C: dict[float, str] = {}
    for lab, c in cand:
        if c in seen:
            continue
        seen.add(c)
        Cs.append(c)
        label_by_C[c] = lab

    # Integer m grid (for EM/strong/weak fits)
    m_all = list(range(int(args.m_min), int(args.m_max) + 1))

    # Fit EM/strong/weak once (best hit within tolerance)
    def _best_hit(target_key: str) -> object | None:
        hits = scan_candidates(Cs=Cs, m_values=m_all, target_G=targets[target_key])
        hits = filter_hits_by_rel_err(hits, max_abs_rel_err=args.max_rel_err)
        return hits[0] if hits else None

    best_em = _best_hit(args.em_target)
    best_s = _best_hit(args.strong_target)
    best_w = _best_hit(args.weak_target)
    if best_em is None or best_s is None or best_w is None:
        print("No baseline EM/strong/weak fits under gauge-derived Cs at this tolerance.")
        print(f"EM best: {best_em is not None}, strong best: {best_s is not None}, weak best: {best_w is not None}")
        return 0

    # Gravity band => implied m window for fixed K
    bands: dict[str, tuple[float, float]] = {
        "cmb": (1e-18, 1e-16),
        "pta": (1e-9, 1e-7),
        "lisa": (1e-4, 1e-1),
        "ligo": (10.0, 1000.0),
    }
    fmin, fmax = bands[args.gravity_band]
    phi_val = (1.0 + math.sqrt(5.0)) / 2.0
    base_hz = constants.BOLTZMANN * args.gravity_K / constants.PLANCK
    # m = log_phi(F0/base)
    m_lo = math.log(fmin / base_hz, phi_val)
    m_hi = math.log(fmax / base_hz, phi_val)
    m_g_min = max(int(math.ceil(min(m_lo, m_hi))), int(args.m_min))
    m_g_max = min(int(math.floor(max(m_lo, m_hi))), int(args.m_max))
    m_band = list(range(m_g_min, m_g_max + 1)) if m_g_min <= m_g_max else []

    if not m_band:
        print("No integer m values fall inside the GW band for the provided K and m range.")
        return 0

    # Build log-spaced GeV scales
    if args.scale_min_GeV <= 0 or args.scale_max_GeV <= 0:
        raise SystemExit("--scale-min-GeV/--scale-max-GeV must be > 0")
    if args.n_scales < 2:
        scales = [float(args.scale_min_GeV)]
    else:
        a = math.log10(float(args.scale_min_GeV))
        b = math.log10(float(args.scale_max_GeV))
        scales = [10 ** (a + (b - a) * i / (args.n_scales - 1)) for i in range(args.n_scales)]

    print("\nQuantum-gravity sweep (gauge-derived C only):")
    print(f"- C candidates: {len(Cs)} (base={args.base}, include={include})")
    print(
        f"- baseline fits: EM({args.em_target}) C={best_em.C:g}, m={int(best_em.m)}, err={best_em.rel_err:.3e}; "
        f"S({args.strong_target}) C={best_s.C:g}, m={int(best_s.m)}, err={best_s.rel_err:.3e}; "
        f"W({args.weak_target}) C={best_w.C:g}, m={int(best_w.m)}, err={best_w.rel_err:.3e}"
    )
    print(f"- gravity: K={args.gravity_K} K, band={args.gravity_band} => m window {m_g_min}..{m_g_max} (n={len(m_band)})")
    print(f"- scanning scales: {args.scale_min_GeV:g}..{args.scale_max_GeV:g} GeV (n={len(scales)})")
    print("")

    results: list[tuple[float, float, int, float, float, str]] = []
    # tuple: (scaleGeV, targetVal, m_g, rel_err_g, F0_g, C_label)

    for gev in scales:
        mass_kg = mass_kg_from_GeV(gev)
        aG = alpha_gravity(mass_kg)
        target_g = (1.0 / aG) if args.gravity_mode == "inverse" else aG
        best: tuple[float, int, float] | None = None  # (C, m, rel_err)
        for C in Cs:
            for m in m_band:
                Gp = gauge_G(C, m)
                rel_err = (Gp - target_g) / target_g if target_g != 0 else float("nan")
                if abs(rel_err) <= args.max_rel_err:
                    if best is None or abs(rel_err) < abs(best[2]):
                        best = (C, m, rel_err)
        if best is not None:
            C, m_g, rel_err_g = best
            F0_g = frequency_F0(m_g, args.gravity_K)
            results.append((gev, target_g, m_g, rel_err_g, F0_g, label_by_C.get(C, "")))

    results.sort(key=lambda r: abs(r[3]))
    print(f"Found {len(results)} passing scales (gravity coupling + GW band + |rel_err|<= {args.max_rel_err}).\n")
    print("rank  scale(GeV)     m_g   C-label              rel_err_g     F0_g(Hz)")
    for i, (gev, _tg, m_g, rel_err_g, F0_g, lab) in enumerate(results[: max(1, args.top)], start=1):
        print(f"{i:4d}  {gev:11.4g}  {m_g:4d}  {lab:18s}  {rel_err_g: .3e}  {F0_g: .3e}")
    return 0


def cmd_gut_run(args: argparse.Namespace) -> int:
    """
    1-loop GUT-style running test: run alpha1_GUT, alpha2, alpha3 from mZ upward and
    find the scale with best convergence (by inverse-coupling RMS).

    This is a standard "do couplings meet?" diagnostic, not a proof of a specific GUT.
    """
    # Reference scale (mZ)
    mZ = 91.1876
    targets = {t.name: t.value for t in known_targets()}
    # Use our derived targets for alpha1_GUT and alpha2 at mZ and alpha_s(mZ) as alpha3.
    a1 = targets["alpha1_GUT(alpha(mZ),sin2)"]
    a2 = targets["alpha2(alpha(mZ),sin2)"]
    a3 = targets["alpha_s(mZ)"]

    inv_a1_0 = 1.0 / a1
    inv_a2_0 = 1.0 / a2
    inv_a3_0 = 1.0 / a3

    betas = SM_1LOOP if args.model == "sm" else MSSM_1LOOP

    mu_best, score, inv1, inv2, inv3 = find_best_convergence(
        mu0=mZ,
        alpha1_inv_mu0=inv_a1_0,
        alpha2_inv_mu0=inv_a2_0,
        alpha3_inv_mu0=inv_a3_0,
        betas=betas,
        mu_min=args.Q_min_GeV,
        mu_max=args.Q_max_GeV,
        n=args.n,
    )

    print(f"Model={betas.name}")
    print("Inputs at mZ (approx):")
    print(f"  alpha1_GUT(mZ) = {a1:.12g}  (inv={inv_a1_0:.6g})")
    print(f"  alpha2(mZ)     = {a2:.12g}  (inv={inv_a2_0:.6g})")
    print(f"  alpha3(mZ)     = {a3:.12g}  (inv={inv_a3_0:.6g})\n")
    print(f"Best convergence in scan: Q ~ {mu_best:.6g} GeV; score(max delta alpha^-1)={score:.6g}")
    print(f"  inv_a1={inv1:.6g}, inv_a2={inv2:.6g}, inv_a3={inv3:.6g}")
    return 0


def cmd_gut_run_lattice(args: argparse.Namespace) -> int:
    """
    1-loop GUT-style running test, but with **lattice-quantized inputs**:

      1) Fit φ-lattice anchors for:
           - α1_GUT^{-1}(mZ)   (hypercharge, GUT-normalized)
           - α2^{-1}(mZ)       (weak SU(2))
           - α3^{-1}(mZ)       (QCD), obtained by fitting α3^{-1}(mH) then running to mZ
      2) Run the resulting inverse couplings from mZ to high scales and find the best
         convergence point, using either SM or MSSM 1-loop beta coefficients.

    This is an exploratory diagnostic: it asks whether the lattice-quantized couplings
    imply “more convergent” running than the raw measured inputs.
    """

    # Reference scale (mZ)
    mZ = float(get_measurement("mZ_GeV", default_value=91.1876).value)

    # Integer m grid
    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))

    # Gauge-derived C candidates
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    cand: list[tuple[str, float]] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for k, v in cs.items():
            cand.append((f"{g.name}:{k}", float(v)))
    # De-duplicate Cs (keep first label)
    seen: set[float] = set()
    Cs: list[float] = []
    label_by_C: dict[float, str] = {}
    for lab, c in cand:
        if c in seen:
            continue
        seen.add(c)
        Cs.append(c)
        label_by_C[c] = lab

    # Targets
    all_targets = {t.name: t for t in known_targets()}
    key_a1 = "1/alpha1_GUT(alpha(mZ),sin2_on_shell)"
    key_a2 = "1/alpha2(alpha(mZ),sin2_on_shell)"
    key_a3_anchor = "1/alpha_s_1loop_from_mZ(mH)"
    for k in (key_a1, key_a2, key_a3_anchor):
        if k not in all_targets:
            raise SystemExit(f"Missing required target {k!r}. Run `list-targets`.")

    # Fit lattice anchors
    best_a1 = scan_candidates(Cs=Cs, m_values=m_values, target_G=float(all_targets[key_a1].value))[0]
    best_a2 = scan_candidates(Cs=Cs, m_values=m_values, target_G=float(all_targets[key_a2].value))[0]
    best_a3 = scan_candidates(Cs=Cs, m_values=m_values, target_G=float(all_targets[key_a3_anchor].value))[0]

    inv_a1_mZ = float(best_a1.G)
    inv_a2_mZ = float(best_a2.G)

    # Run strong anchor to mZ (deterministic; no tuning knobs)
    Q0_strong = float(scale_GeV_from_target_key(key_a3_anchor))
    spec_qcd = qcd_run_spec_from_key(key_a3_anchor)
    alpha3_Q0 = (1.0 / float(best_a3.G)) if float(best_a3.G) != 0 else float("inf")
    alpha3_mZ = alpha_s_from_ref(mZ, alpha_s_Q0=alpha3_Q0, Q0_GeV=Q0_strong, spec=spec_qcd)
    inv_a3_mZ = (1.0 / float(alpha3_mZ)) if alpha3_mZ not in (0.0, float("inf")) else float("inf")

    betas = SM_1LOOP if args.model == "sm" else MSSM_1LOOP
    mu_best, score, inv1, inv2, inv3 = find_best_convergence(
        mu0=mZ,
        alpha1_inv_mu0=inv_a1_mZ,
        alpha2_inv_mu0=inv_a2_mZ,
        alpha3_inv_mu0=inv_a3_mZ,
        betas=betas,
        mu_min=args.Q_min_GeV,
        mu_max=args.Q_max_GeV,
        n=args.n,
    )

    print(f"Model={betas.name}")
    print("Lattice-quantized inputs (at mZ):")
    print(
        f"  inv_alpha1_GUT(mZ) = {inv_a1_mZ:.12g}  "
        f"(fit: {label_by_C.get(float(best_a1.C), ''):22s} C={best_a1.C:g}, m={int(best_a1.m):d}, rel_err={best_a1.rel_err:.3e})"
    )
    print(
        f"  inv_alpha2(mZ)     = {inv_a2_mZ:.12g}  "
        f"(fit: {label_by_C.get(float(best_a2.C), ''):22s} C={best_a2.C:g}, m={int(best_a2.m):d}, rel_err={best_a2.rel_err:.3e})"
    )
    print(
        f"  inv_alpha3(mZ)     = {inv_a3_mZ:.12g}  "
        f"(fit anchor @Q0={Q0_strong:g} GeV: {label_by_C.get(float(best_a3.C), ''):22s} C={best_a3.C:g}, m={int(best_a3.m):d}, rel_err={best_a3.rel_err:.3e}; "
        f"run to mZ with QCD {spec_qcd.loops}L, nf_mode={spec_qcd.nf_mode})\n"
    )
    print(f"Best convergence in scan: Q ~ {mu_best:.6g} GeV; score(max delta alpha^-1)={score:.6g}")
    print(f"  inv_a1={inv1:.6g}, inv_a2={inv2:.6g}, inv_a3={inv3:.6g}")
    return 0


def cmd_pair_forces_gaugeCs(args: argparse.Namespace) -> int:
    """
    Full pairing test using ONLY gauge-derived C candidates (non-arbitrary):
      - Coupling fit: G = C/phi^m matches each force target within tolerance
      - Quantum forces (EM/strong/weak): Option-2, so user supplies F0 and we solve K
      - Gravity: fixed K (default CMB), optional GW band filter on predicted F0
    """
    targets = {t.name: t.value for t in known_targets()}

    gravity_targets = [s.strip() for s in args.gravity_targets.split(",") if s.strip()]
    for k in (args.em_target, args.strong_target, args.weak_target, *gravity_targets):
        if k not in targets:
            raise SystemExit(f"Unknown target {k!r}. Run `python -m physics_test.cli list-targets`.")

    # Frequencies (Option-2)
    def _F0_from_arg(value_hz: float | None, preset_key: str | None, label: str) -> float:
        if preset_key:
            return float(get_preset(preset_key).F0_hz)
        if value_hz is None:
            raise SystemExit(f"Must provide either --{label}-F0 or --{label}-preset")
        return float(value_hz)

    F0_em = _F0_from_arg(args.em_F0, args.em_preset, "em")
    F0_s = _F0_from_arg(args.strong_F0, args.strong_preset, "strong")
    F0_w = _F0_from_arg(args.weak_F0, args.weak_preset, "weak")

    # Integer m grid
    m_values = sorted(set(int(round(x)) for x in frange(args.m_min, args.m_max, args.m_step)))

    # Gauge-derived C candidates
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    cand: list[tuple[str, float]] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for k, v in cs.items():
            cand.append((f"{g.name}:{k}", float(v)))
    # de-dupe C values, keep first label
    seen: set[float] = set()
    Cs: list[float] = []
    label_by_C: dict[float, str] = {}
    for lab, c in cand:
        if c in seen:
            continue
        seen.add(c)
        Cs.append(c)
        label_by_C[c] = lab

    # Optionally force include 360 only (strict mode)
    if args.C360_only:
        Cs = [c for c in Cs if abs(c - 360.0) < 1e-12]
        if not Cs:
            raise SystemExit("No C candidates left after --C360-only")

    # Build hit lists for each target
    def _hits_for(target_key: str) -> list[object]:
        hits = scan_candidates(Cs=Cs, m_values=m_values, target_G=targets[target_key])
        hits = filter_hits_by_rel_err(hits, max_abs_rel_err=args.max_rel_err)
        return hits[: max(1, args.max_hits)]

    em_hits = _hits_for(args.em_target)
    s_hits = _hits_for(args.strong_target)
    w_hits = _hits_for(args.weak_target)
    g_hits_by_key = {gk: _hits_for(gk) for gk in gravity_targets}

    # Gravity band preset (applies only when using fixed K)
    gravity_bands: dict[str, tuple[float | None, float | None]] = {
        "any": (None, None),
        "ligo": (10.0, 1000.0),
        "lisa": (1e-4, 1e-1),
        "pta": (1e-9, 1e-7),
        "cmb": (1e-18, 1e-16),
    }
    if args.gravity_band not in gravity_bands:
        raise SystemExit(f"Unknown gravity band {args.gravity_band!r}")
    g_f0_min, g_f0_max = gravity_bands[args.gravity_band]
    if args.gravity_f0_min is not None:
        g_f0_min = args.gravity_f0_min
    if args.gravity_f0_max is not None:
        g_f0_max = args.gravity_f0_max

    def _k_ok(K: float, kmin: float | None, kmax: float | None) -> bool:
        if kmin is not None and K < kmin:
            return False
        if kmax is not None and K > kmax:
            return False
        return True

    def _f0_ok(F0: float) -> bool:
        if g_f0_min is not None and F0 < g_f0_min:
            return False
        if g_f0_max is not None and F0 > g_f0_max:
            return False
        return True

    print("\nGauge-C pairing (Option-2 for EM/strong/weak):")
    print(f"- C candidates (unique) = {len(Cs)} from base={args.base} include={include}")
    print(
        f"- targets: EM={args.em_target}, S={args.strong_target}, W={args.weak_target}, "
        f"G in {gravity_targets}"
    )
    print(f"- tol(|rel_err|) = {args.max_rel_err}")
    print(f"- F0: EM={F0_em:.6e} Hz, S={F0_s:.6e} Hz, W={F0_w:.6e} Hz")
    print(f"- gravity: K={args.gravity_K} K, band={args.gravity_band}, F0 bounds=[{g_f0_min},{g_f0_max}]")
    print("")

    shown = 0
    for eh in em_hits:
        m_em = int(eh.m)
        if args.em_m_sign == "positive" and m_em <= 0:
            continue
        if args.em_m_sign == "negative" and m_em >= 0:
            continue
        K_em = temperature_K_from_frequency(m_em, F0_em)
        if not _k_ok(K_em, args.em_K_min, args.em_K_max):
            continue

        for sh in s_hits:
            m_s = int(sh.m)
            if args.strong_m_sign == "positive" and m_s <= 0:
                continue
            if args.strong_m_sign == "negative" and m_s >= 0:
                continue
            K_s = temperature_K_from_frequency(m_s, F0_s)
            if not _k_ok(K_s, args.strong_K_min, args.strong_K_max):
                continue

            for wh in w_hits:
                m_w = int(wh.m)
                if args.weak_m_sign == "positive" and m_w <= 0:
                    continue
                if args.weak_m_sign == "negative" and m_w >= 0:
                    continue
                K_w = temperature_K_from_frequency(m_w, F0_w)
                if not _k_ok(K_w, args.weak_K_min, args.weak_K_max):
                    continue

                for gk, ghits in g_hits_by_key.items():
                    for gh in ghits:
                        m_g = int(gh.m)
                        if args.gravity_m_sign == "positive" and m_g <= 0:
                            continue
                        if args.gravity_m_sign == "negative" and m_g >= 0:
                            continue

                        K_g = args.gravity_K
                        if not _k_ok(K_g, args.gravity_K_min, args.gravity_K_max):
                            continue
                        F0_g = frequency_F0(m_g, K_g)
                        if not _f0_ok(F0_g):
                            continue

                        shown += 1
                        print(
                            f"#{shown}: m=[EM {m_em}, S {m_s}, W {m_w}, G {m_g}] "
                            f"rel_err=[{eh.rel_err:.2e}, {sh.rel_err:.2e}, {wh.rel_err:.2e}, {gh.rel_err:.2e}]"
                        )
                        print(
                            f"  EM: {label_by_C.get(eh.C,''):18s} C={eh.C:g}  K={K_em:.6g} K"
                        )
                        print(
                            f"  S : {label_by_C.get(sh.C,''):18s} C={sh.C:g}  K={K_s:.6g} K"
                        )
                        print(
                            f"  W : {label_by_C.get(wh.C,''):18s} C={wh.C:g}  K={K_w:.6g} K"
                        )
                        print(
                            f"  G : {gk:18s} {label_by_C.get(gh.C,''):18s} C={gh.C:g}  "
                            f"K={K_g:.6g} K -> F0={F0_g:.3e} Hz"
                        )
                        print("")

                        if shown >= args.max_results:
                            print(f"Stopped after max_results={args.max_results}.")
                            return 0

    print(f"Done. Found {shown} configurations.")
    return 0


def cmd_scan(args: argparse.Namespace) -> int:
    targets = {t.name: t.value for t in known_targets()}
    if args.target not in targets:
        raise SystemExit(f"Unknown --target {args.target!r}. Options: {', '.join(targets.keys())}")

    Cs = _parse_Cs(args.Cs)
    if not Cs:
        s = get_candidate_set(args.set_name)
        Cs = list(s.values)
        print(f"Using candidate set: {s.name} ({len(Cs)} values) - {s.note}")
    else:
        print(f"Using custom Cs list: n={len(Cs)}")

    m_values = frange(args.m_min, args.m_max, args.m_step)
    if args.m_integer:
        m_values = sorted(set(int(round(x)) for x in m_values))
    hits = scan_candidates(Cs=Cs, m_values=m_values, target_G=targets[args.target])
    if args.max_rel_err is not None:
        hits = filter_hits_by_rel_err(hits, max_abs_rel_err=args.max_rel_err)

    print(f"Target {args.target} = {targets[args.target]:.15g}")
    print(f"Scanned: nC={len(Cs)} x nm={len(m_values)} = {len(Cs)*len(m_values)} combos")
    if args.max_rel_err is not None:
        print(f"Kept: {len(hits)} hits with |rel_err| <= {args.max_rel_err}\n")
    else:
        print(f"Ranked: {len(hits)} combos\n")
    print("rank  C           m        G=C/phi^m           rel_err")
    for i, h in enumerate(hits[: max(1, args.top)], start=1):
        print(f"{i:4d}  {h.C:10.6g}  {h.m:7.3f}  {h.G:18.12g}  {h.rel_err: .6e}")
    return 0


def cmd_scan_all(args: argparse.Namespace) -> int:
    all_targets = known_targets()

    Cs = _parse_Cs(args.Cs)
    if not Cs:
        s = get_candidate_set(args.set_name)
        Cs = list(s.values)
        print(f"Using candidate set: {s.name} ({len(Cs)} values) - {s.note}")
    else:
        print(f"Using custom Cs list: n={len(Cs)}")

    m_values = frange(args.m_min, args.m_max, args.m_step)
    if args.m_integer:
        m_values = sorted(set(int(round(x)) for x in m_values))

    print(f"Scanning all targets: nTargets={len(all_targets)}; nC={len(Cs)}; nm={len(m_values)}\n")
    for t in all_targets:
        hits = scan_candidates(Cs=Cs, m_values=m_values, target_G=t.value)
        best = hits[0]
        ok = abs(best.rel_err) <= args.max_rel_err
        status = "PASS" if ok else "FAIL"
        print(
            f"[{status}] {t.name:18s} target={t.value:.12g}  "
            f"best: C={best.C:g}, m={best.m:g}, G={best.G:.12g}, rel_err={best.rel_err:.3e}"
        )
    return 0


def _classify_sector(name: str) -> str:
    n = name.lower()
    if "alpha_g" in n or "gravity" in n:
        return "GRAVITY"
    if "alpha_s" in n or "strong" in n:
        return "STRONG"
    if "alpha2" in n or "alpha_w" in n or "weak" in n or "sin2theta" in n or "delta_r" in n or "delta_rho" in n:
        return "WEAK/EW"
    if "alpha1" in n or "alpha_y" in n or "hypercharge" in n:
        return "HYPERCHARGE"
    if "alpha" in n or "delta_alpha" in n:
        return "EM"
    if "over" in n or "mp_over" in n or "mmu_over" in n or "mtau_over" in n or "mt_over" in n or "mb_over" in n or "mw_over" in n:
        return "MASS RATIO"
    if "unification" in n:
        return "UNIFICATION"
    return "OTHER"


def cmd_spectrum(args: argparse.Namespace) -> int:
    """
    Phi-lattice spectrum: find the best (C, m) for each target and display
    organized by m-value, revealing the lattice topology.
    """
    all_targets = known_targets()

    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    Cs: list[float] = []
    labels: list[str] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for k, v in cs.items():
            Cs.append(float(v))
            labels.append(f"{g.name}:{k}")
    seen: set[float] = set()
    uniq_Cs: list[float] = []
    label_by_C: dict[float, str] = {}
    for c, lab in zip(Cs, labels):
        if c in seen:
            continue
        seen.add(c)
        uniq_Cs.append(c)
        label_by_C[c] = lab

    m_values = list(range(args.m_min, args.m_max + 1))
    max_err = args.max_rel_err

    filter_str = args.filter.lower() if args.filter else ""
    results: list[tuple[int, float, str, float, str, str]] = []

    for t in all_targets:
        if filter_str and filter_str not in t.name.lower():
            continue
        hits = scan_candidates(Cs=uniq_Cs, m_values=m_values, target_G=t.value)
        best = hits[0]
        if abs(best.rel_err) > max_err:
            continue
        m_int = int(best.m)
        lab = label_by_C.get(best.C, f"C={best.C:g}")
        sector = _classify_sector(t.name)
        results.append((m_int, best.C, t.name, best.rel_err, lab, sector))

    results.sort(key=lambda r: (r[0], abs(r[3])))

    print(f"PHI-LATTICE SPECTRUM: G = C / phi^m")
    print(f"Gauge-derived C menu: {sorted(set(uniq_Cs))}")
    print(f"m range: [{args.m_min}, {args.m_max}], max |rel_err| = {max_err}")
    if filter_str:
        print(f"Filter: '{args.filter}'")
    print(f"Hits: {len(results)}\n")

    header = f"{'m':>5s}  {'C':>6s}  {'target':40s}  {'rel_err':>10s}  {'sector':12s}  {'gauge label'}"
    print(header)
    print("-" * len(header) + "----------")

    prev_m = None
    for m_int, C, name, rel_err, lab, sector in results:
        if prev_m is not None and m_int != prev_m:
            print()
        prev_m = m_int
        err_pct = f"{rel_err*100:+.2f}%"
        print(f"{m_int:5d}  {C:6g}  {name:40s}  {err_pct:>10s}  {sector:12s}  {lab}")

    if not args.no_summary:
        print(f"\n{'='*80}")
        print("SECTOR SUMMARY BY m:")
        sector_ms: dict[str, list[int]] = {}
        for m_int, _, _, _, _, sector in results:
            sector_ms.setdefault(sector, []).append(m_int)
        for sector in ["EM", "STRONG", "WEAK/EW", "HYPERCHARGE", "GRAVITY", "MASS RATIO", "UNIFICATION", "OTHER"]:
            if sector not in sector_ms:
                continue
            ms = sorted(set(sector_ms[sector]))
            rng = f"[{min(ms)}, {max(ms)}]"
            print(f"  {sector:14s}  m in {rng:16s}  n_hits={len(sector_ms[sector]):3d}  m-values: {ms}")

    return 0


def cmd_gut_run_mssm(args: argparse.Namespace) -> int:
    """
    Compare SM vs MSSM GUT convergence at 1-loop and 2-loop, with:
      - Lattice-quantized inputs
      - Non-minimal GUT normalization scan
      - Lattice-constrained unification search
    """
    from physics_test.target_registry import get_measurement

    m_mZ = get_measurement("mZ_GeV", default_value=91.1880)
    mZ = float(m_mZ.value)

    m_inv_alpha_mZ = get_measurement("alpha_inv_mZ", default_value=127.930)
    inv_alpha_mZ = float(m_inv_alpha_mZ.value)
    alpha_mZ = 1.0 / inv_alpha_mZ

    m_sin2 = get_measurement("sin2thetaW_mZ_MSbar", default_value=0.23129)
    sin2_mZ = float(m_sin2.value)

    m_alpha_s = get_measurement("alpha_s_mZ", default_value=0.1180)
    alpha_s_mZ = float(m_alpha_s.value)

    sin2_on_shell = 1.0 - (float(get_measurement("mW_GeV", default_value=80.3692).value) ** 2) / (mZ ** 2)
    cos2_on_shell = 1.0 - sin2_on_shell
    alpha2_on_shell = alpha_mZ / sin2_on_shell
    alpha1_gut_on_shell = (5.0 / 3.0) * alpha_mZ / cos2_on_shell

    inv_a1_mZ = 1.0 / alpha1_gut_on_shell
    inv_a2_mZ = 1.0 / alpha2_on_shell
    inv_a3_mZ = 1.0 / alpha_s_mZ

    # -- Build lattice C candidates for quantized inputs + lattice-constrained search --
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for _, v in cs.items():
            all_Cs.append(float(v))
    all_Cs = sorted(set(all_Cs))
    m_range = list(range(-10, 20))

    # =====================================================================
    # SECTION 1:  SM vs MSSM  --  1-loop and 2-loop comparison
    # =====================================================================
    models_1l = [("SM", SM_1LOOP), ("MSSM", MSSM_1LOOP)]
    models_2l = [("SM", SM_2LOOP), ("MSSM", MSSM_2LOOP)]

    print("=" * 72)
    print("GUT CONVERGENCE: SM vs MSSM (1-loop AND 2-loop)")
    print("=" * 72)
    print(f"Inputs at mZ = {mZ} GeV (on-shell EW definition, SU(5) k1=5/3):")
    print(f"  alpha1_GUT^-1(mZ) = {inv_a1_mZ:.4f}")
    print(f"  alpha2^-1(mZ)     = {inv_a2_mZ:.4f}")
    print(f"  alpha3^-1(mZ)     = {inv_a3_mZ:.4f}")
    print()

    for (name, betas_1l), (_, betas_2l) in zip(models_1l, models_2l):
        # --- 1-loop ---
        mu_best, score, a1, a2, a3 = find_best_convergence(
            mu0=mZ,
            alpha1_inv_mu0=inv_a1_mZ,
            alpha2_inv_mu0=inv_a2_mZ,
            alpha3_inv_mu0=inv_a3_mZ,
            betas=betas_1l,
            mu_min=1e3, mu_max=1e19, n=2000,
        )
        print(f"  {name:6s} 1-LOOP  Q_GUT = {mu_best:.3e} GeV  score = {score:.3f}")
        print(f"         a1^-1 = {a1:.3f}   a2^-1 = {a2:.3f}   a3^-1 = {a3:.3f}")

        # --- 2-loop coupled (coarse scan -- pure-Python RK4 is slow) ---
        mu_best_2, score_2, a1_2, a2_2, a3_2 = find_best_convergence_2loop(
            mu0=mZ,
            alpha1_mu0=alpha1_gut_on_shell,
            alpha2_mu0=alpha2_on_shell,
            alpha3_mu0=alpha_s_mZ,
            betas2=betas_2l,
            mu_min=1e3, mu_max=1e19, n=200, steps_per_unit_log=80,
        )
        delta = score_2 - score
        arrow = "better" if delta < 0 else "worse"
        print(f"  {name:6s} 2-LOOP  Q_GUT = {mu_best_2:.3e} GeV  score = {score_2:.3f}  ({arrow}, delta={delta:+.3f})")
        print(f"         a1^-1 = {a1_2:.3f}   a2^-1 = {a2_2:.3f}   a3^-1 = {a3_2:.3f}")

        # --- Lattice-quantized 1-loop ---
        hits_a1 = scan_candidates(Cs=all_Cs, m_values=m_range, target_G=inv_a1_mZ)
        hits_a2 = scan_candidates(Cs=all_Cs, m_values=m_range, target_G=inv_a2_mZ)
        hits_a3 = scan_candidates(Cs=all_Cs, m_values=m_range, target_G=inv_a3_mZ)
        best_a1, best_a2, best_a3 = hits_a1[0], hits_a2[0], hits_a3[0]

        mu_lat, score_lat, a1_lat, a2_lat, a3_lat = find_best_convergence(
            mu0=mZ,
            alpha1_inv_mu0=best_a1.G,
            alpha2_inv_mu0=best_a2.G,
            alpha3_inv_mu0=best_a3.G,
            betas=betas_1l,
            mu_min=1e3, mu_max=1e19, n=2000,
        )
        print(f"  {name:6s} LATTICE-QUANTIZED (1-loop):")
        print(f"         anchors: a1^-1={best_a1.G:.3f}(C={best_a1.C:g},m={int(best_a1.m)})  "
              f"a2^-1={best_a2.G:.3f}(C={best_a2.C:g},m={int(best_a2.m)})  "
              f"a3^-1={best_a3.G:.3f}(C={best_a3.C:g},m={int(best_a3.m)})")
        print(f"         Q_GUT = {mu_lat:.3e} GeV  score = {score_lat:.3f}")
        print(f"         a1^-1 = {a1_lat:.3f}   a2^-1 = {a2_lat:.3f}   a3^-1 = {a3_lat:.3f}")
        print()

    # =====================================================================
    # SECTION 2:  Non-minimal GUT normalization scan  (discrete k1 values)
    # =====================================================================
    print("=" * 72)
    print("NON-MINIMAL GUT NORMALIZATION SCAN")
    print("=" * 72)
    print("Scanning k₁ in α₁^GUT = k₁·α_Y for known GUT groups.")
    print()

    for model_name, betas_1l in models_1l:
        print(f"  --- {model_name} (1-loop) ---")
        norm_results = scan_gut_normalizations(
            alpha_mZ=alpha_mZ,
            cos2_mZ=cos2_on_shell,
            alpha_s_mZ=alpha_s_mZ,
            sin2_mZ=sin2_on_shell,
            mu0=mZ,
            betas=betas_1l,
            mu_min=1e3, mu_max=1e19, n=500,
        )
        print(f"  {'Rank':>4s}  {'k1':>6s}  {'GUT group':<30s}  {'Q_GUT':>11s}  {'Score':>7s}  {'a_GUT^-1':>9s}")
        print(f"  {'----':>4s}  {'------':>6s}  {'-'*30:<30s}  {'-'*11:>11s}  {'-'*7:>7s}  {'-'*9:>9s}")
        for rank, nr in enumerate(norm_results, 1):
            a_avg = (nr.inv_a1 + nr.inv_a2 + nr.inv_a3) / 3.0
            print(f"  {rank:4d}  {nr.k1:6.3f}  {nr.name:<30s}  {nr.mu_best:11.3e}  {nr.score:7.3f}  {a_avg:9.3f}")
        print()

    # =====================================================================
    # SECTION 3:  Lattice-constrained unification
    # =====================================================================
    print("=" * 72)
    print("LATTICE-CONSTRAINED UNIFICATION")
    print("=" * 72)
    print("Searching for Q where all three α_i^-1(Q) ≈ C/φ^m (same C,m).")
    print()

    for model_name, betas_1l in models_1l:
        print(f"  {model_name} (1-loop):")
        lc_results = find_lattice_gut_point(
            mu0=mZ,
            alpha1_inv_mu0=inv_a1_mZ,
            alpha2_inv_mu0=inv_a2_mZ,
            alpha3_inv_mu0=inv_a3_mZ,
            betas=betas_1l,
            Cs=all_Cs, m_range=m_range,
            mu_min=1e3, mu_max=1e19, n_mu=1000, top=5,
        )
        _print_lattice_gut_results(lc_results)

    return 0


def _print_lattice_gut_results(results: list) -> None:
    if not results:
        print("    (no results found)")
        print()
        return
    print(f"    {'#':>3s}  {'C':>6s}  {'m':>3s}  {'C/φ^m':>8s}  {'Q (GeV)':>11s}  {'max_dev':>8s}  {'rms_dev':>8s}  {'a1^-1':>8s}  {'a2^-1':>8s}  {'a3^-1':>8s}")
    for i, r in enumerate(results, 1):
        print(f"    {i:3d}  {r.C:6g}  {r.m:3d}  {r.lattice_value:8.3f}  {r.Q_GeV:11.3e}  {r.max_dev:8.3f}  {r.rms_dev:8.3f}  {r.inv_a1:8.3f}  {r.inv_a2:8.3f}  {r.inv_a3:8.3f}")
    best = results[0]
    print(f"    BEST: (C={best.C:g}, m={best.m}) => {best.lattice_value:.3f} at Q = {best.Q_GeV:.3e} GeV  max_dev = {best.max_dev:.3f}")
    print()


def cmd_gut_trajectory(args: argparse.Namespace) -> int:
    """
    GUT trajectory diagnostic: show how the three couplings' lattice addresses
    evolve from mZ to the Planck scale, and test consistency of the GUT scale
    with the phi-lattice.
    """
    import math as _math
    from physics_test.target_registry import get_measurement
    from physics_test import constants

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0

    m_mZ = get_measurement("mZ_GeV", default_value=91.1880)
    mZ = float(m_mZ.value)
    m_inv_alpha_mZ = get_measurement("alpha_inv_mZ", default_value=127.930)
    inv_alpha_mZ = float(m_inv_alpha_mZ.value)
    alpha_mZ = 1.0 / inv_alpha_mZ
    m_alpha_s = get_measurement("alpha_s_mZ", default_value=0.1180)
    alpha_s_mZ = float(m_alpha_s.value)

    sin2_on_shell = 1.0 - (float(get_measurement("mW_GeV", default_value=80.3692).value) ** 2) / (mZ ** 2)
    cos2_on_shell = 1.0 - sin2_on_shell
    alpha2_on_shell = alpha_mZ / sin2_on_shell
    alpha1_gut_on_shell = (5.0 / 3.0) * alpha_mZ / cos2_on_shell

    inv_a1_mZ = 1.0 / alpha1_gut_on_shell
    inv_a2_mZ = 1.0 / alpha2_on_shell
    inv_a3_mZ = 1.0 / alpha_s_mZ

    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for _, v in cs.items():
            all_Cs.append(float(v))
    all_Cs = sorted(set(all_Cs))
    m_range = list(range(-50, 50))

    def nearest_lattice(value: float) -> tuple:
        best_C, best_m, best_lv, best_err = 0.0, 0, 0.0, float("inf")
        for C in all_Cs:
            for m_int in m_range:
                lv = C / (PHI ** m_int)
                err = abs(value - lv)
                if err < best_err:
                    best_err = err
                    best_C, best_m, best_lv = C, m_int, lv
        return best_C, best_m, best_lv, best_err

    # =====================================================================
    # SECTION 1: MSSM 1-loop trajectory from mZ to M_Planck
    # =====================================================================
    print("=" * 80)
    print("MSSM RG TRAJECTORY: lattice addresses from mZ to M_Planck")
    print("=" * 80)
    print()
    print(f"  {'Q (GeV)':>11s}  |  {'a1^-1':>8s} ->{'(C,m)':>10s}  |  {'a2^-1':>8s} ->{'(C,m)':>10s}  |  {'a3^-1':>8s} ->{'(C,m)':>10s}  |  {'score':>6s}")
    print(f"  {'-'*11}  |  {'-'*8}  {'-'*10}  |  {'-'*8}  {'-'*10}  |  {'-'*8}  {'-'*10}  |  {'-'*6}")

    log_min = _math.log(mZ)
    log_max = _math.log(1e19)
    n_steps = 30
    gut_Q = None
    gut_score = float("inf")

    for i in range(n_steps + 1):
        t = i / n_steps
        Q = _math.exp(log_min + (log_max - log_min) * t)
        a1 = run_alpha_inv(inv_a1_mZ, mZ, Q, MSSM_1LOOP.b1)
        a2 = run_alpha_inv(inv_a2_mZ, mZ, Q, MSSM_1LOOP.b2)
        a3 = run_alpha_inv(inv_a3_mZ, mZ, Q, MSSM_1LOOP.b3)

        if a1 <= 0 or a2 <= 0 or a3 <= 0:
            continue

        C1, m1, _, e1 = nearest_lattice(a1)
        C2, m2, _, e2 = nearest_lattice(a2)
        C3, m3, _, e3 = nearest_lattice(a3)
        score = converge_score(a1, a2, a3)

        marker = ""
        if score < gut_score:
            gut_score = score
            gut_Q = Q
        if abs(Q - 1.6e16) / 1.6e16 < 0.5:
            marker = "  <-- GUT region"

        print(f"  {Q:11.3e}  |  {a1:8.3f}  ({C1:g},{m1:d})  |  {a2:8.3f}  ({C2:g},{m2:d})  |  {a3:8.3f}  ({C3:g},{m3:d})  |  {score:6.2f}{marker}")

    print()

    # =====================================================================
    # SECTION 2: GUT scale consistency checks
    # =====================================================================
    from physics_test.gut import find_best_convergence
    mu_best, best_score, a1_gut, a2_gut, a3_gut = find_best_convergence(
        mu0=mZ,
        alpha1_inv_mu0=inv_a1_mZ,
        alpha2_inv_mu0=inv_a2_mZ,
        alpha3_inv_mu0=inv_a3_mZ,
        betas=MSSM_1LOOP,
        mu_min=1e3, mu_max=1e19, n=2000,
    )
    a_gut_avg = (a1_gut + a2_gut + a3_gut) / 3.0
    C_gut, m_gut, lv_gut, err_gut = nearest_lattice(a_gut_avg)

    print("=" * 80)
    print("GUT SCALE CONSISTENCY CHECKS")
    print("=" * 80)
    print()
    print(f"  MSSM best convergence: Q_GUT = {mu_best:.4e} GeV  (score = {best_score:.3f})")
    print(f"  Couplings at Q_GUT:  a1^-1 = {a1_gut:.4f}   a2^-1 = {a2_gut:.4f}   a3^-1 = {a3_gut:.4f}")
    print(f"  Average a_GUT^-1 = {a_gut_avg:.4f}")
    print(f"  Nearest lattice point: (C={C_gut:g}, m={m_gut}) => {lv_gut:.4f}  (dev = {err_gut:.4f})")
    print()

    # Energy hierarchy in phi-units
    M_Planck_GeV = constants.MASS_PLANCK * constants.SPEED_OF_LIGHT ** 2 / 1.6022e-10
    ratio_gut_mZ = mu_best / mZ
    ratio_Pl_gut = M_Planck_GeV / mu_best
    ratio_Pl_mZ = M_Planck_GeV / mZ

    m_gut_mZ = _math.log(ratio_gut_mZ) / _math.log(PHI)
    m_Pl_gut = _math.log(ratio_Pl_gut) / _math.log(PHI)
    m_Pl_mZ = _math.log(ratio_Pl_mZ) / _math.log(PHI)

    print("  Energy hierarchy in φ-units:")
    print(f"    Q_GUT / mZ       = {ratio_gut_mZ:.4e} = φ^{m_gut_mZ:.2f}  (nearest integer: {round(m_gut_mZ)})")
    print(f"    M_Planck / Q_GUT = {ratio_Pl_gut:.4e} = φ^{m_Pl_gut:.2f}  (nearest integer: {round(m_Pl_gut)})")
    print(f"    M_Planck / mZ    = {ratio_Pl_mZ:.4e} = φ^{m_Pl_mZ:.2f}  (nearest integer: {round(m_Pl_mZ)})")
    print(f"    Check: {round(m_gut_mZ)} + {round(m_Pl_gut)} = {round(m_gut_mZ) + round(m_Pl_gut)}  (should ≈ {round(m_Pl_mZ)})")
    print()

    # Gravitational coupling at GUT scale
    M_gut_kg = mu_best * 1.6022e-10 / (constants.SPEED_OF_LIGHT ** 2)
    alpha_G_gut = constants.G_NEWTON * M_gut_kg ** 2 / (constants.HBAR * constants.SPEED_OF_LIGHT)
    inv_alpha_G_gut = 1.0 / alpha_G_gut

    C_grav, m_grav, lv_grav, err_grav = nearest_lattice(inv_alpha_G_gut)
    rel_err_grav = (lv_grav - inv_alpha_G_gut) / inv_alpha_G_gut

    print("  Gravitational coupling at Q_GUT:")
    print(f"    M_GUT = {M_gut_kg:.4e} kg")
    print(f"    α_G(M_GUT) = {alpha_G_gut:.4e}")
    print(f"    1/α_G(M_GUT) = {inv_alpha_G_gut:.4e}")
    print(f"    Nearest lattice: (C={C_grav:g}, m={m_grav}) => {lv_grav:.4e}  (rel_err = {rel_err_grav:+.2%})")
    print()

    # The unified coupling vs gravity gap
    print("  Unification-gravity bridge:")
    gap_coupling = inv_alpha_G_gut / a_gut_avg
    gap_m = _math.log(gap_coupling) / _math.log(PHI)
    print(f"    1/α_G(M_GUT) / α_GUT^-1 = {gap_coupling:.4e}")
    print(f"    This gap = φ^{gap_m:.2f}  (nearest integer: {round(gap_m)})")
    print(f"    Coupling m-address of α_GUT:   m = {m_gut}")
    print(f"    Coupling m-address of α_G(GUT): m = {m_grav}")
    print(f"    Lattice distance: Δm = {m_grav - m_gut}")
    print()

    # Proton lifetime estimate from lattice-quantized values
    if lv_gut > 0:
        m_proton_GeV = 0.93827
        alpha_gut_lat = 1.0 / lv_gut
        tau_p_factor = mu_best ** 4 / (alpha_gut_lat ** 2 * m_proton_GeV ** 5)
        tau_p_years = tau_p_factor * 1.6022e-10 / (constants.SPEED_OF_LIGHT * constants.HBAR) * (1.0 / (3.156e7))
        log10_tau = _math.log10(abs(tau_p_factor))
        print("  Proton lifetime (dimensional estimate):")
        print(f"    Using α_GUT = 1/{lv_gut:.3f} and M_GUT = {mu_best:.3e} GeV:")
        print(f"    τ_p ∝ M_GUT^4 / (α_GUT^2 · m_p^5)")
        print(f"    log10(M_GUT^4 / (α_GUT^2 · m_p^5)) = {log10_tau:.1f}")
        print(f"    (Super-K bound: log10(τ_p/yr) > 34;  Hyper-K target: ~35)")
    print()

    return 0


def cmd_gut_validate(args: argparse.Namespace) -> int:
    """
    GUT validation suite — four independent tests to strengthen the (C=15, m=-1) finding.

    1. Independent lattice predictions from physical mass ratios / energy scales
    2. Tightened 2-loop search around (C=15, m=-1)
    3. SU(5) threshold corrections with lattice-quantized mass splittings
    4. Fibonacci / golden-ratio structure in energy hierarchy exponents
    """
    import math as _math
    from physics_test.target_registry import get_measurement
    from physics_test import constants
    from physics_test.gut import (
        MSSM_1LOOP, MSSM_2LOOP,
        find_best_convergence, converge_score, run_alpha_inv,
    )
    from physics_test.gut_validate import (
        build_sorted_lattice, nearest_lattice as nl,
        independent_scale_tests, lattice_coverage,
        tightened_2loop_search,
        scan_su5_thresholds,
        zeckendorf, fibonacci_index, lucas_index,
        phi_power_decompose, fibonacci_up_to, lucas_up_to,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0

    # ----- Common setup -----
    mZ = float(get_measurement("mZ_GeV", default_value=91.1880).value)
    inv_alpha_mZ = float(get_measurement("alpha_inv_mZ", default_value=127.930).value)
    alpha_mZ = 1.0 / inv_alpha_mZ
    alpha_s_mZ = float(get_measurement("alpha_s_mZ", default_value=0.1180).value)
    mW = float(get_measurement("mW_GeV", default_value=80.3692).value)

    sin2_os = 1.0 - (mW ** 2) / (mZ ** 2)
    cos2_os = 1.0 - sin2_os
    alpha2_os = alpha_mZ / sin2_os
    alpha1_gut_os = (5.0 / 3.0) * alpha_mZ / cos2_os

    inv_a1_mZ = 1.0 / alpha1_gut_os
    inv_a2_mZ = 1.0 / alpha2_os
    inv_a3_mZ = 1.0 / alpha_s_mZ

    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=args.base, include=include)
        for _, v in cs.items():
            all_Cs.append(float(v))
    all_Cs = sorted(set(all_Cs))
    m_range = list(range(-80, 120))

    lattice = build_sorted_lattice(all_Cs, m_range)

    # ===================================================================
    # SECTION 1: INDEPENDENT LATTICE PREDICTIONS
    # ===================================================================
    print("=" * 80)
    print("SECTION 1: INDEPENDENT LATTICE PREDICTIONS")
    print("=" * 80)
    print()
    print("Dimensionless ratios NOT used in lattice construction or fitting:")
    print()

    preds = independent_scale_tests(lattice)

    print(f"  {'Name':<18s}  {'Value':>14s}  {'C':>6s}  {'m':>4s}  {'C/φ^m':>14s}  {'rel%':>7s}  Note")
    print(f"  {'-'*18}  {'-'*14}  {'-'*6}  {'-'*4}  {'-'*14}  {'-'*7}  {'-'*25}")

    hits_1pct = 0
    hits_3pct = 0
    hits_5pct = 0
    for p in preds:
        pct = p.rel_err * 100.0
        marker = "***" if abs(pct) < 1.0 else "** " if abs(pct) < 3.0 else "*  " if abs(pct) < 5.0 else "   "
        if abs(pct) < 1.0:
            hits_1pct += 1
        if abs(pct) < 3.0:
            hits_3pct += 1
        if abs(pct) < 5.0:
            hits_5pct += 1

        if p.value > 1e6:
            val_str = f"{p.value:.4e}"
            lat_str = f"{p.lattice_value:.4e}"
        else:
            val_str = f"{p.value:.6f}"
            lat_str = f"{p.lattice_value:.6f}"
        print(f"  {p.name:<18s}  {val_str:>14s}  {p.C:6g}  {p.m:4d}  {lat_str:>14s}  {pct:+6.2f}%  {marker} {p.note}")

    total = len(preds)
    print()
    print(f"  Hits within 1%: {hits_1pct}/{total}")
    print(f"  Hits within 3%: {hits_3pct}/{total}")
    print(f"  Hits within 5%: {hits_5pct}/{total}")
    print()

    print("  Null-hypothesis (random log-uniform values):")
    for tol_pct, tol_frac in [(1, 0.01), (3, 0.03), (5, 0.05)]:
        null_rate = lattice_coverage(lattice, tol_frac, n_samples=5000)
        obs_rate = (hits_1pct if tol_pct == 1 else hits_3pct if tol_pct == 3 else hits_5pct) / total
        enrichment = obs_rate / null_rate if null_rate > 0 else float("inf")
        print(f"    {tol_pct}% tolerance:  null = {null_rate:.1%},  observed = {obs_rate:.1%},  enrichment = {enrichment:.2f}x")
    print()

    # ===================================================================
    # SECTION 2: TIGHTENED 2-LOOP SEARCH AROUND (C=15, m=-1)
    # ===================================================================
    print("=" * 80)
    print("SECTION 2: TIGHTENED 2-LOOP SEARCH AROUND (C=15, m=-1)")
    print("=" * 80)
    print()

    C_tgt = 15.0
    m_tgt = -1
    lv_tgt = C_tgt / (PHI ** m_tgt)
    print(f"  Target lattice point: C={C_tgt:g}, m={m_tgt}")
    print(f"  Lattice value: C/φ^m = {lv_tgt:.6f}")
    print()

    mu1l_best, _, a1_1l, a2_1l, a3_1l = find_best_convergence(
        mu0=mZ, alpha1_inv_mu0=inv_a1_mZ, alpha2_inv_mu0=inv_a2_mZ,
        alpha3_inv_mu0=inv_a3_mZ, betas=MSSM_1LOOP,
        mu_min=1e3, mu_max=1e19, n=3000,
    )
    score_1l = converge_score(a1_1l, a2_1l, a3_1l)
    dev_1l = [abs(x - lv_tgt) for x in (a1_1l, a2_1l, a3_1l)]

    print(f"  1-loop baseline (MSSM):")
    print(f"    Q_GUT = {mu1l_best:.4e} GeV")
    print(f"    α₁⁻¹ = {a1_1l:.4f}   α₂⁻¹ = {a2_1l:.4f}   α₃⁻¹ = {a3_1l:.4f}")
    print(f"    Convergence score = {score_1l:.4f}")
    print(f"    Deviations from {lv_tgt:.3f}:  {dev_1l[0]:.4f}  {dev_1l[1]:.4f}  {dev_1l[2]:.4f}  max = {max(dev_1l):.4f}")
    print()

    print("  Running 2-loop coupled RK4 (MSSM) — narrow scan ±0.7 decades around Q_GUT...")
    res2 = tightened_2loop_search(
        Q0_GeV=mZ,
        alpha1_0=alpha1_gut_os,
        alpha2_0=alpha2_os,
        alpha3_0=alpha_s_mZ,
        C_target=C_tgt, m_target=m_tgt,
        Q_center=mu1l_best,
        Q_half_decades=0.7,
        n_Q=200, steps_per_unit_log=200,
    )

    b = res2['best']
    bc = res2['best_conv']
    if b:
        print(f"  2-loop result (closest to lattice point):")
        print(f"    Q = {b['Q']:.4e} GeV")
        print(f"    α₁⁻¹ = {b['inv1']:.4f}   α₂⁻¹ = {b['inv2']:.4f}   α₃⁻¹ = {b['inv3']:.4f}")
        print(f"    Deviations from {lv_tgt:.3f}:  {b['d1']:+.4f}  {b['d2']:+.4f}  {b['d3']:+.4f}")
        print(f"    max_dev = {b['max_dev']:.4f}   rms_dev = {b['rms_dev']:.4f}")
        print(f"    Convergence score = {b['score']:.4f}")
    if bc and bc != b:
        print(f"  2-loop result (best convergence):")
        print(f"    Q = {bc['Q']:.4e} GeV")
        print(f"    α₁⁻¹ = {bc['inv1']:.4f}   α₂⁻¹ = {bc['inv2']:.4f}   α₃⁻¹ = {bc['inv3']:.4f}")
        print(f"    Convergence score = {bc['score']:.4f}")

    print()
    if b:
        improvement = "IMPROVES" if b['max_dev'] < max(dev_1l) else "does not improve"
        print(f"  Comparison: 2-loop {improvement} distance to (C=15, m=-1)")
        print(f"    1-loop max_dev = {max(dev_1l):.4f}")
        print(f"    2-loop max_dev = {b['max_dev']:.4f}")
        if max(dev_1l) > 0:
            ratio = b['max_dev'] / max(dev_1l)
            print(f"    Ratio (2L/1L) = {ratio:.3f}")

    # Sample trajectory points
    traj = res2['trajectory']
    if traj:
        print()
        print(f"  Trajectory sample (10 points):")
        print(f"    {'Q (GeV)':>12s}  {'α₁⁻¹':>8s}  {'α₂⁻¹':>8s}  {'α₃⁻¹':>8s}  {'max_dev':>8s}  {'score':>8s}")
        step = max(1, len(traj) // 10)
        for entry in traj[::step]:
            print(f"    {entry['Q']:12.3e}  {entry['inv1']:8.3f}  {entry['inv2']:8.3f}  {entry['inv3']:8.3f}  {entry['max_dev']:8.3f}  {entry['score']:8.3f}")
    print()

    # ===================================================================
    # SECTION 3: SU(5) THRESHOLD CORRECTIONS
    # ===================================================================
    print("=" * 80)
    print("SECTION 3: SU(5) THRESHOLD CORRECTIONS (LATTICE-QUANTIZED)")
    print("=" * 80)
    print()
    print("  Mass splittings: M_V = M_GUT × φ^δ_V,  M_HC = M_GUT × φ^δ_HC")
    print("  Correction: Δα_i⁻¹ = -(b_i^a / 2π) ln(M_a/M_GUT)")
    print()

    a_avg_1l = (a1_1l + a2_1l + a3_1l) / 3.0
    no_corr_score = converge_score(a1_1l, a2_1l, a3_1l)
    print(f"  Baseline (no corrections): score = {no_corr_score:.4f}")
    print(f"    α₁⁻¹ = {a1_1l:.4f}   α₂⁻¹ = {a2_1l:.4f}   α₃⁻¹ = {a3_1l:.4f}")
    print()

    for model, susy in [("SUSY SU(5)", True), ("Non-SUSY SU(5)", False)]:
        print(f"  --- {model} ---")
        thresh = scan_su5_thresholds(
            inv_a1_gut=a1_1l, inv_a2_gut=a2_1l, inv_a3_gut=a3_1l,
            delta_range=range(-5, 6), susy=susy,
        )
        print(f"  {'#':>3s}  {'δ_V':>4s}  {'δ_HC':>5s}  {'score':>7s}  {'impr':>7s}  {'α₁⁻¹':>8s}  {'α₂⁻¹':>8s}  {'α₃⁻¹':>8s}  {'near latt':>12s}")
        for i, r in enumerate(thresh[:10]):
            avg_corrected = (r['inv_a1'] + r['inv_a2'] + r['inv_a3']) / 3.0
            C_near, m_near, lv_near, err_near = nl(avg_corrected, lattice)
            lat_str = f"({C_near:g},{m_near})"
            impr_str = f"+{r['improvement']:.4f}" if r['improvement'] > 0 else f"{r['improvement']:.4f}"
            print(f"  {i+1:3d}  {r['delta_V']:4d}  {r['delta_HC']:5d}  {r['score']:7.4f}  {impr_str:>7s}  {r['inv_a1']:8.4f}  {r['inv_a2']:8.4f}  {r['inv_a3']:8.4f}  {lat_str:>12s}")

        best_t = thresh[0]
        M_V_ratio = PHI ** best_t['delta_V']
        M_HC_ratio = PHI ** best_t['delta_HC']
        print()
        print(f"  Best: δ_V={best_t['delta_V']}, δ_HC={best_t['delta_HC']}")
        print(f"    M_V/M_GUT = φ^{best_t['delta_V']} = {M_V_ratio:.4f}")
        print(f"    M_HC/M_GUT = φ^{best_t['delta_HC']} = {M_HC_ratio:.4f}")
        print(f"    Score improvement: {no_corr_score:.4f} → {best_t['score']:.4f} ({best_t['improvement']:+.4f})")

        # Check if corrected couplings hit a lattice point
        avg_best = (best_t['inv_a1'] + best_t['inv_a2'] + best_t['inv_a3']) / 3.0
        C_b, m_b, lv_b, err_b = nl(avg_best, lattice)
        rel_b = (lv_b - avg_best) / avg_best if avg_best != 0 else float("inf")
        print(f"    Corrected avg α_GUT⁻¹ = {avg_best:.4f}")
        print(f"    Nearest lattice: (C={C_b:g}, m={m_b}) = {lv_b:.4f}  (rel = {rel_b:+.4f})")
        print()

    # ===================================================================
    # SECTION 4: FIBONACCI STRUCTURE IN ENERGY HIERARCHY
    # ===================================================================
    print("=" * 80)
    print("SECTION 4: FIBONACCI / GOLDEN-RATIO STRUCTURE")
    print("=" * 80)
    print()

    M_Planck_GeV = constants.MASS_PLANCK * constants.SPEED_OF_LIGHT ** 2 / 1.6022e-10
    ratio_gut_mZ = mu1l_best / mZ
    ratio_Pl_gut = M_Planck_GeV / mu1l_best
    ratio_Pl_mZ = M_Planck_GeV / mZ

    exp_gut = _math.log(ratio_gut_mZ) / _math.log(PHI)
    exp_Pl_gut = _math.log(ratio_Pl_gut) / _math.log(PHI)
    exp_Pl_mZ = _math.log(ratio_Pl_mZ) / _math.log(PHI)

    exponents = {
        "Q_GUT/mZ": exp_gut,
        "M_Pl/Q_GUT": exp_Pl_gut,
        "M_Pl/mZ": exp_Pl_mZ,
    }

    print("  4a. Phi-unit exponents of the energy hierarchy:")
    for name, exp in exponents.items():
        print(f"    {name:15s} = φ^{exp:.4f}  (nearest int: {round(exp)})")
    print()

    print("  4b. Zeckendorf representations (non-consecutive Fibonacci sums):")
    for name, exp in exponents.items():
        n = round(exp)
        z = zeckendorf(n)
        fi = fibonacci_index(n)
        li = lucas_index(n)
        fib_str = f"F({fi})" if fi else "not Fibonacci"
        luc_str = f"L({li})" if li else "not Lucas"
        z_str = " + ".join(str(x) for x in z)
        print(f"    {n:4d} = {z_str}")
        print(f"         {fib_str};  {luc_str}")
    print()

    n_gut = round(exp_gut)
    n_gap = round(exp_Pl_gut)
    n_total = round(exp_Pl_mZ)

    print("  4c. Self-similar hierarchy structure:")
    ratio_total_gap = n_total / n_gap if n_gap != 0 else float("inf")
    phi_power = _math.log(ratio_total_gap) / _math.log(PHI) if ratio_total_gap > 0 else float("nan")
    print(f"    n_total / n_gap = {n_total}/{n_gap} = {ratio_total_gap:.4f}")
    print(f"    φ^4 = {PHI**4:.4f}")
    print(f"    Ratio ≈ φ^{phi_power:.4f}")
    print(f"    => M_Pl/mZ exponent ≈ (M_Pl/Q_GUT exponent) × φ^4")
    print()
    predicted_total = n_gap * PHI ** 4
    predicted_desert = n_gap * (PHI ** 4 - 1)
    print(f"    If n_gap = {n_gap}:")
    print(f"      Predicted n_total = {n_gap} × φ^4 = {predicted_total:.2f}  (actual: {n_total}, dev: {n_total - predicted_total:+.2f})")
    print(f"      Predicted n_desert = {n_gap} × (φ^4-1) = {predicted_desert:.2f}  (actual: {n_gut}, dev: {n_gut - predicted_desert:+.2f})")
    print()

    print("  4d. φ^4 decomposition (since φ^4 = 3φ+2 = F(4)φ+F(3)):")
    print(f"    φ^4 = {PHI**4:.6f}")
    print(f"    3φ + 2 = {3*PHI + 2:.6f}")
    print(f"    So: n_total ≈ {n_gap}(3φ+2) = {n_gap*3}φ + {n_gap*2} = {n_gap*3*PHI + n_gap*2:.2f}")
    print()

    print("  4e. F(12) = 144 = 12² — unique Fibonacci square:")
    fibs = fibonacci_up_to(200)
    fib_squares = [(i + 1, f) for i, f in enumerate(fibs) if _math.isqrt(f) ** 2 == f and f > 1]
    for idx, val in fib_squares:
        root = _math.isqrt(val)
        print(f"    F({idx}) = {val} = {root}²")
    print(f"    The Planck-GUT gap ({n_gap}) is the ONLY index k>1 where F(k) = k²")
    print()

    print("  4f. α_GUT⁻¹ ≈ dim(SU(5)):")
    lv_15_m1 = 15.0 * PHI
    print(f"    C/φ^m = 15/φ^(-1) = 15φ = {lv_15_m1:.4f}")
    print(f"    dim(SU(5)) = 5²-1 = 24")
    print(f"    Deviation: {lv_15_m1:.4f} - 24 = {lv_15_m1 - 24:.4f}  ({(lv_15_m1 - 24)/24*100:.2f}%)")
    print()

    print("  4g. C=15 group-theoretic derivations:")
    print(f"    15 = 360/dim(SU(5)) = 360/24")
    print(f"    15 = dim(SU(4)) = 4²-1")
    print(f"    15 = dim(fundamental_antisymmetric_2-tensor of SU(6)) = C(6,2)")
    print(f"    15 = dim(adjoint of SU(4) = Pati-Salam color)")
    print()

    print("  4h. Cross-scale Fibonacci multiples of n_gap:")
    for k in range(0, 6):
        phi_k = PHI ** k
        predicted = n_gap * phi_k
        predicted_int = round(predicted)
        scale_GeV = mZ * PHI ** predicted_int
        log10_GeV = _math.log10(scale_GeV) if scale_GeV > 0 else 0
        print(f"    {n_gap} × φ^{k} = {predicted:.2f} → {predicted_int}  (φ^{predicted_int} × mZ = 10^{log10_GeV:.1f} GeV)")
    print()

    print("  4i. Exponent relationships to SM group dimensions:")
    sm_dims = [("U(1)", 1), ("SU(2)", 3), ("SU(3)", 8), ("SM total", 12),
               ("SU(5)", 24), ("SO(10)", 45), ("E6", 78)]
    for gname, d in sm_dims:
        if d > 0:
            ratio_to_gap = n_gap / d
            ratio_to_desert = n_gut / d
            ratio_to_total = n_total / d
            print(f"    dim({gname:6s}) = {d:3d}:  gap/{d} = {ratio_to_gap:.2f},  desert/{d} = {ratio_to_desert:.2f},  total/{d} = {ratio_to_total:.2f}")
    print()

    return 0


def cmd_gut_significance(args: argparse.Namespace) -> int:
    """
    Statistical significance tests for the phi-lattice:
      A. Base-uniqueness permutation test
      B. C=15 clustering test
      C. Pre-registered predictions (out-of-sample validation)
    """
    import math as _math
    from physics_test.gut_validate import (
        build_sorted_lattice, nearest_lattice as nl,
        independent_scale_tests, lattice_coverage,
        base_permutation_test, c_clustering_test,
        pre_registered_predictions,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    m_range = list(range(-80, 120))

    # ===================================================================
    # STRATEGY A: BASE-UNIQUENESS PERMUTATION TEST
    # ===================================================================
    print("=" * 80)
    print("STRATEGY A: BASE-UNIQUENESS PERMUTATION TEST")
    print("=" * 80)
    print()
    print("For each base, generate C menu from SM gauge groups, build lattice,")
    print("run 20 independent predictions, and compare enrichment vs null.")
    print()

    bases = [60, 100, 120, 180, 200, 240, 270, 300, 315, 330, 350, 355, 359,
             360, 361, 365, 370, 380, 400, 420, 450, 500, 600, 720, 1000, 1080]

    perm_results = base_permutation_test(bases, include=include, m_range=m_range, null_samples=3000)

    print(f"  {'Base':>6s}  {'#Cs':>4s}  {'<1%':>4s}  {'<3%':>4s}  {'<5%':>4s}  {'null1%':>6s}  {'null3%':>6s}  {'null5%':>6s}  {'enr1%':>6s}  {'enr3%':>6s}  {'enr5%':>6s}  C values")
    print(f"  {'------':>6s}  {'----':>4s}  {'----':>4s}  {'----':>4s}  {'----':>4s}  {'------':>6s}  {'------':>6s}  {'------':>6s}  {'------':>6s}  {'------':>6s}  {'------':>6s}  --------")

    best_enr_1 = max(perm_results, key=lambda r: r.enrichment_1pct)
    for r in perm_results:
        marker = " <-- BEST" if r.base == best_enr_1.base else ""
        if r.base == 360:
            marker = " <-- 360" + (" (BEST)" if r.base == best_enr_1.base else "")
        cs_str = ",".join(f"{c:g}" for c in r.Cs[:8])
        if len(r.Cs) > 8:
            cs_str += "..."
        print(f"  {r.base:6g}  {r.n_Cs:4d}  {r.hits_1pct:4d}  {r.hits_3pct:4d}  {r.hits_5pct:4d}  "
              f"{r.null_1pct:5.1%}  {r.null_3pct:5.1%}  {r.null_5pct:5.1%}  "
              f"{r.enrichment_1pct:6.2f}  {r.enrichment_3pct:6.2f}  {r.enrichment_5pct:6.2f}  "
              f"{cs_str}{marker}")

    print()

    r360 = next(r for r in perm_results if r.base == 360)
    all_enr1 = [r.enrichment_1pct for r in perm_results]
    rank_360 = sorted(all_enr1, reverse=True).index(r360.enrichment_1pct) + 1

    print(f"  Summary:")
    print(f"    Base 360 enrichment (1%): {r360.enrichment_1pct:.2f}x  (rank {rank_360}/{len(perm_results)})")
    print(f"    Best base: {best_enr_1.base:g} with enrichment {best_enr_1.enrichment_1pct:.2f}x")
    print(f"    Mean enrichment (1%): {sum(all_enr1)/len(all_enr1):.2f}x")
    print(f"    Std dev: {(_math.sqrt(sum((x - sum(all_enr1)/len(all_enr1))**2 for x in all_enr1) / len(all_enr1))):.2f}")
    if r360.enrichment_1pct > 0:
        sigma_above = (r360.enrichment_1pct - sum(all_enr1) / len(all_enr1))
        std_dev = _math.sqrt(sum((x - sum(all_enr1) / len(all_enr1)) ** 2 for x in all_enr1) / len(all_enr1))
        if std_dev > 0:
            sigma_above /= std_dev
            print(f"    Base 360 is {sigma_above:.1f}σ above mean")
    print()

    # ===================================================================
    # STRATEGY B: C=15 CLUSTERING TEST
    # ===================================================================
    print("=" * 80)
    print("STRATEGY B: C-VALUE CLUSTERING TEST")
    print("=" * 80)
    print()

    all_Cs_360: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=360.0, include=include)
        for _, v in cs.items():
            all_Cs_360.append(float(v))
    all_Cs_360 = sorted(set(all_Cs_360))

    lattice_360 = build_sorted_lattice(all_Cs_360, m_range)
    preds_360 = independent_scale_tests(lattice_360)

    cluster_results = c_clustering_test(preds_360, all_Cs_360, n_trials=100000)

    print(f"  C menu ({len(all_Cs_360)} values): {', '.join(f'{c:g}' for c in all_Cs_360)}")
    print(f"  Each of {len(preds_360)} ratios maps to its nearest (C, m).")
    print()
    print(f"  Observed C-value clustering:")
    print(f"  {'C':>8s}  {'count':>5s}  {'p-value':>8s}  Ratios")
    print(f"  {'--------':>8s}  {'-----':>5s}  {'--------':>8s}  ------")

    for cr in cluster_results:
        pv_str = f"{cr.p_value:.4f}" if cr.p_value >= 0.0001 else f"{cr.p_value:.1e}"
        sig = ""
        if cr.p_value < 0.01:
            sig = "  ***"
        elif cr.p_value < 0.05:
            sig = "  **"
        elif cr.p_value < 0.10:
            sig = "  *"
        names_m = [f"{n}(m={m})" for n, m in zip(cr.ratio_names, cr.m_values)]
        print(f"  {cr.C:8g}  {cr.count:5d}  {pv_str:>8s}  {', '.join(names_m)}{sig}")

    print()
    c15 = next((cr for cr in cluster_results if cr.C == 15), None)
    if c15:
        print(f"  C=15 cluster analysis:")
        print(f"    {c15.count} independent ratios map to C=15")
        print(f"    p-value (seeing this many or more on any single C by chance): {c15.p_value:.4f}")
        if c15.p_value < 0.05:
            print(f"    This is statistically significant at p < 0.05")
        print(f"    Combined with C=15 appearing at GUT unification (an independent finding),")
        print(f"    the recurrence of C=15 across different physics sectors is notable.")
    print()

    # ===================================================================
    # STRATEGY C: PRE-REGISTERED PREDICTIONS (OUT-OF-SAMPLE)
    # ===================================================================
    print("=" * 80)
    print("STRATEGY C: PRE-REGISTERED PREDICTIONS (OUT-OF-SAMPLE)")
    print("=" * 80)
    print()
    print("20 new dimensionless ratios NOT in the original test set.")
    print("These are mixing angles, mass ratios, and coupling combinations")
    print("that were never used in any lattice construction or prior test.")
    print()

    pre_reg = pre_registered_predictions(lattice_360)

    print(f"  {'Name':<18s}  {'Value':>14s}  {'C':>6s}  {'m':>4s}  {'C/φ^m':>14s}  {'rel%':>7s}  Note")
    print(f"  {'-'*18}  {'-'*14}  {'-'*6}  {'-'*4}  {'-'*14}  {'-'*7}  {'-'*25}")

    pr_hits_1 = 0
    pr_hits_3 = 0
    pr_hits_5 = 0
    pr_total = len(pre_reg)

    for p in pre_reg:
        pct = p.predicted_rel_err * 100.0
        marker = "***" if abs(pct) < 1.0 else "** " if abs(pct) < 3.0 else "*  " if abs(pct) < 5.0 else "   "
        if abs(pct) < 1.0:
            pr_hits_1 += 1
        if abs(pct) < 3.0:
            pr_hits_3 += 1
        if abs(pct) < 5.0:
            pr_hits_5 += 1

        if p.value > 1e4:
            val_str = f"{p.value:.4e}"
            lat_str = f"{p.predicted_lattice:.4e}"
        elif p.value < 0.001:
            val_str = f"{p.value:.6e}"
            lat_str = f"{p.predicted_lattice:.6e}"
        else:
            val_str = f"{p.value:.6f}"
            lat_str = f"{p.predicted_lattice:.6f}"
        print(f"  {p.name:<18s}  {val_str:>14s}  {p.predicted_C:6g}  {p.predicted_m:4d}  {lat_str:>14s}  {pct:+6.2f}%  {marker} {p.note}")

    print()
    print(f"  Pre-registered hit rates:")
    print(f"    < 1%: {pr_hits_1}/{pr_total} ({pr_hits_1/pr_total:.0%})")
    print(f"    < 3%: {pr_hits_3}/{pr_total} ({pr_hits_3/pr_total:.0%})")
    print(f"    < 5%: {pr_hits_5}/{pr_total} ({pr_hits_5/pr_total:.0%})")
    print()

    null_1 = lattice_coverage(lattice_360, 0.01, 3000)
    null_3 = lattice_coverage(lattice_360, 0.03, 3000)
    null_5 = lattice_coverage(lattice_360, 0.05, 3000)

    print(f"  Comparison with null (base=360 lattice):")
    for tol_pct, pr_hits, null in [(1, pr_hits_1, null_1), (3, pr_hits_3, null_3), (5, pr_hits_5, null_5)]:
        obs = pr_hits / pr_total if pr_total > 0 else 0
        enr = obs / null if null > 0 else float("inf")
        print(f"    {tol_pct}%:  observed = {obs:.0%},  null = {null:.1%},  enrichment = {enr:.2f}x")
    print()

    # C-value clustering in pre-registered set
    pr_c_counts: dict[float, list] = {}
    for p in pre_reg:
        c_val = p.predicted_C
        if c_val not in pr_c_counts:
            pr_c_counts[c_val] = []
        pr_c_counts[c_val].append(p.name)

    print(f"  C-value distribution in pre-registered set:")
    for c_val in sorted(pr_c_counts.keys()):
        names = pr_c_counts[c_val]
        print(f"    C={c_val:6g}: {len(names)} hits — {', '.join(names)}")
    print()

    # ===================================================================
    # COMBINED SUMMARY
    # ===================================================================
    print("=" * 80)
    print("COMBINED SIGNIFICANCE SUMMARY")
    print("=" * 80)
    print()

    orig_1pct = r360.hits_1pct
    orig_total = r360.total
    combined_1pct = orig_1pct + pr_hits_1
    combined_total = orig_total + pr_total
    combined_rate = combined_1pct / combined_total if combined_total > 0 else 0
    combined_enr = combined_rate / null_1 if null_1 > 0 else float("inf")

    print(f"  Original test set (20 ratios):   {orig_1pct}/{orig_total} < 1% hits ({orig_1pct/orig_total:.0%})")
    print(f"  Pre-registered set (20 ratios):  {pr_hits_1}/{pr_total} < 1% hits ({pr_hits_1/pr_total:.0%})")
    print(f"  Combined (40 ratios):            {combined_1pct}/{combined_total} < 1% hits ({combined_rate:.0%})")
    print(f"  Combined enrichment (1%):        {combined_enr:.2f}x over null ({null_1:.1%})")
    print()
    print(f"  Base uniqueness:  360 ranks #{rank_360} of {len(perm_results)} bases tested")
    if c15:
        print(f"  C=15 clustering:  {c15.count} ratios on C=15, p-value = {c15.p_value:.4f}")
    print()

    return 0


def cmd_gut_deep(args: argparse.Namespace) -> int:
    """
    Deep analysis: address table, CKM/PMNS predictions, lattice operations,
    360-family analysis, n_gap self-consistency.
    """
    import math as _math
    from physics_test.gut_validate import (
        build_sorted_lattice, nearest_lattice as nl,
        lattice_coverage,
        full_mixing_matrix_predictions,
        lattice_operation_analysis,
        base_family_analysis,
        base_permutation_test,
        build_address_table,
        ngap_self_consistency,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    m_range = list(range(-80, 120))

    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=360.0, include=include)
        for _, v in cs.items():
            all_Cs.append(float(v))
    all_Cs = sorted(set(all_Cs))

    lattice = build_sorted_lattice(all_Cs, m_range)

    # ===================================================================
    # 1. SYSTEMATIC LATTICE ADDRESS TABLE
    # ===================================================================
    print("=" * 80)
    print("1. SYSTEMATIC LATTICE ADDRESS TABLE")
    print("=" * 80)
    print()

    table = build_address_table(lattice)
    print(f"  {len(table)} dimensionless ratios mapped to nearest (C, m)")
    print()
    print(f"  {'Name':<20s}  {'Value':>14s}  {'Sector':<10s}  {'C':>6s}  {'m':>4s}  {'C/φ^m':>14s}  {'rel%':>7s}")
    print(f"  {'-'*20}  {'-'*14}  {'-'*10}  {'-'*6}  {'-'*4}  {'-'*14}  {'-'*7}")

    hits_1 = hits_3 = hits_5 = 0
    c_distribution: dict[float, list] = {}
    sector_stats: dict[str, dict] = {}

    for row in table:
        pct = row['rel_err'] * 100.0
        marker = "***" if abs(pct) < 1.0 else "** " if abs(pct) < 3.0 else "*  " if abs(pct) < 5.0 else "   "
        if abs(pct) < 1.0:
            hits_1 += 1
        if abs(pct) < 3.0:
            hits_3 += 1
        if abs(pct) < 5.0:
            hits_5 += 1

        c_val = row['C']
        if c_val not in c_distribution:
            c_distribution[c_val] = []
        c_distribution[c_val].append(row['name'])

        sec = row['sector']
        if sec not in sector_stats:
            sector_stats[sec] = {'total': 0, 'h1': 0, 'h3': 0, 'h5': 0}
        sector_stats[sec]['total'] += 1
        if abs(pct) < 1.0:
            sector_stats[sec]['h1'] += 1
        if abs(pct) < 3.0:
            sector_stats[sec]['h3'] += 1
        if abs(pct) < 5.0:
            sector_stats[sec]['h5'] += 1

        if row['value'] > 1e4:
            val_str = f"{row['value']:.4e}"
            lat_str = f"{row['lattice_val']:.4e}"
        elif row['value'] < 0.001:
            val_str = f"{row['value']:.6e}"
            lat_str = f"{row['lattice_val']:.6e}"
        else:
            val_str = f"{row['value']:.6f}"
            lat_str = f"{row['lattice_val']:.6f}"

        print(f"  {row['name']:<20s}  {val_str:>14s}  {row['sector']:<10s}  {row['C']:6g}  {row['m']:4d}  {lat_str:>14s}  {pct:+6.2f}% {marker}")

    total = len(table)
    null_1 = lattice_coverage(lattice, 0.01, 3000)
    null_3 = lattice_coverage(lattice, 0.03, 3000)
    null_5 = lattice_coverage(lattice, 0.05, 3000)

    print()
    print(f"  Hit rate summary ({total} ratios):")
    for tol, h, null in [(1, hits_1, null_1), (3, hits_3, null_3), (5, hits_5, null_5)]:
        obs = h / total
        enr = obs / null if null > 0 else float("inf")
        print(f"    <{tol}%: {h}/{total} ({obs:.0%}),  null = {null:.1%},  enrichment = {enr:.2f}x")

    print()
    print(f"  C-value distribution (address population):")
    for c_val in sorted(c_distribution.keys()):
        names = c_distribution[c_val]
        print(f"    C={c_val:6g}: {len(names):2d} ratios — {', '.join(names[:6])}" + ("..." if len(names) > 6 else ""))

    print()
    print(f"  Sector-by-sector accuracy:")
    print(f"  {'Sector':<12s}  {'Total':>5s}  {'<1%':>4s}  {'<3%':>4s}  {'<5%':>4s}  {'Rate<1%':>7s}")
    for sec in sorted(sector_stats.keys()):
        s = sector_stats[sec]
        r1 = s['h1'] / s['total'] if s['total'] > 0 else 0
        print(f"  {sec:<12s}  {s['total']:5d}  {s['h1']:4d}  {s['h3']:4d}  {s['h5']:4d}  {r1:6.0%}")

    # Look for m-value patterns within each C
    print()
    print(f"  Lattice address patterns:")
    for c_val in sorted(c_distribution.keys()):
        names = c_distribution[c_val]
        if len(names) >= 2:
            m_vals = [row['m'] for row in table if row['C'] == c_val]
            m_vals_sorted = sorted(m_vals)
            if len(m_vals_sorted) >= 2:
                diffs = [m_vals_sorted[i+1] - m_vals_sorted[i] for i in range(len(m_vals_sorted)-1)]
                print(f"    C={c_val:6g}: m values = {m_vals_sorted}, gaps = {diffs}")

    # ===================================================================
    # 2. FULL CKM + PMNS MIXING MATRIX
    # ===================================================================
    print()
    print("=" * 80)
    print("2. FULL CKM + PMNS MIXING MATRIX PREDICTIONS")
    print("=" * 80)
    print()

    mixing = full_mixing_matrix_predictions(lattice)

    print(f"  {'Name':<22s}  {'Value':>14s}  {'C':>6s}  {'m':>4s}  {'C/φ^m':>14s}  {'rel%':>7s}  Note")
    print(f"  {'-'*22}  {'-'*14}  {'-'*6}  {'-'*4}  {'-'*14}  {'-'*7}  {'-'*25}")

    mx_h1 = mx_h3 = mx_h5 = 0
    ckm_c_dist: dict[float, list] = {}
    for p in mixing:
        pct = p.predicted_rel_err * 100.0
        marker = "***" if abs(pct) < 1.0 else "** " if abs(pct) < 3.0 else "*  " if abs(pct) < 5.0 else "   "
        if abs(pct) < 1.0:
            mx_h1 += 1
        if abs(pct) < 3.0:
            mx_h3 += 1
        if abs(pct) < 5.0:
            mx_h5 += 1

        cv = p.predicted_C
        if cv not in ckm_c_dist:
            ckm_c_dist[cv] = []
        ckm_c_dist[cv].append(p.name)

        if p.value > 1e4:
            val_str = f"{p.value:.4e}"
            lat_str = f"{p.predicted_lattice:.4e}"
        elif p.value < 0.001:
            val_str = f"{p.value:.6e}"
            lat_str = f"{p.predicted_lattice:.6e}"
        else:
            val_str = f"{p.value:.6f}"
            lat_str = f"{p.predicted_lattice:.6f}"
        print(f"  {p.name:<22s}  {val_str:>14s}  {p.predicted_C:6g}  {p.predicted_m:4d}  {lat_str:>14s}  {pct:+6.2f}% {marker} {p.note}")

    mx_total = len(mixing)
    print()
    print(f"  CKM+PMNS hit rates ({mx_total} parameters):")
    print(f"    < 1%: {mx_h1}/{mx_total} ({mx_h1/mx_total:.0%})")
    print(f"    < 3%: {mx_h3}/{mx_total} ({mx_h3/mx_total:.0%})")
    print(f"    < 5%: {mx_h5}/{mx_total} ({mx_h5/mx_total:.0%})")

    print()
    print(f"  C-value clustering in mixing sector:")
    for cv in sorted(ckm_c_dist.keys()):
        names = ckm_c_dist[cv]
        print(f"    C={cv:6g}: {len(names)} — {', '.join(names)}")

    # ===================================================================
    # 3. LATTICE OPERATION ANALYSIS (g-2, 2pi, phi maps)
    # ===================================================================
    print()
    print("=" * 80)
    print("3. LATTICE OPERATION ANALYSIS")
    print("=" * 80)
    print()
    print("  How do arithmetic operations (×2π, ÷2π, ×φ, √, etc.) map between")
    print("  lattice addresses? This reveals the internal algebra of the lattice.")
    print()

    ops = lattice_operation_analysis(all_Cs, m_range)

    print(f"  {'Source':<12s}  {'Addr':>10s}  {'Op':>4s}  →  {'Result':>14s}  {'Addr':>10s}  {'ΔC':>6s}  {'Δm':>4s}  {'rel%':>7s}")
    print(f"  {'-'*12}  {'-'*10}  {'-'*4}     {'-'*14}  {'-'*10}  {'-'*6}  {'-'*4}  {'-'*7}")

    # Find operations with small rel errors (clean lattice maps)
    clean_ops = [o for o in ops if abs(o['rel_err']) < 0.05]
    clean_ops.sort(key=lambda o: abs(o['rel_err']))

    for o in clean_ops:
        pct = o['rel_err'] * 100.0
        marker = "***" if abs(pct) < 1.0 else "** " if abs(pct) < 3.0 else "*  " if abs(pct) < 5.0 else "   "
        if o['result_val'] > 1e4:
            res_str = f"{o['result_val']:.4e}"
        elif o['result_val'] < 0.001:
            res_str = f"{o['result_val']:.6e}"
        else:
            res_str = f"{o['result_val']:.6f}"
        print(f"  {o['source']:<12s}  {o['source_addr']:>10s}  {o['operation']:>4s}  →  {res_str:>14s}  {o['result_addr']:>10s}  {o['delta_C']:+6g}  {o['delta_m']:+4d}  {pct:+6.2f}% {marker}")

    # Analyze delta-m patterns for each operation
    print()
    print(f"  δm patterns for operations on clean (<5%) mappings:")
    op_dm: dict[str, list] = {}
    for o in clean_ops:
        op = o['operation']
        if op not in op_dm:
            op_dm[op] = []
        op_dm[op].append(o['delta_m'])
    for op in sorted(op_dm.keys()):
        dms = op_dm[op]
        print(f"    {op:>4s}: δm values = {dms}, mean = {sum(dms)/len(dms):.1f}")

    # ===================================================================
    # 4. 360-FAMILY FACTORIZATION ANALYSIS
    # ===================================================================
    print()
    print("=" * 80)
    print("4. 360-FAMILY FACTORIZATION ANALYSIS")
    print("=" * 80)
    print()

    bases_for_family = [60, 100, 120, 180, 200, 240, 270, 300, 315, 330, 350, 355, 359,
                        360, 361, 365, 370, 380, 400, 420, 450, 500, 600, 720, 1000, 1080]
    perm_results = base_permutation_test(bases_for_family, include=include, m_range=m_range, null_samples=3000)
    family = base_family_analysis(perm_results)

    print(f"  360 = 2³ × 3² × 5")
    print(f"  Correlating enrichment with shared prime structure:")
    print()
    print(f"  {'Base':>6s}  {'GCD':>4s}  {'Shared':>6s}  {'Div360':>6s}  {'Mult':>4s}  {'Enr(1%)':>7s}  {'Hits':>4s}")
    print(f"  {'-'*6}  {'-'*4}  {'-'*6}  {'-'*6}  {'-'*4}  {'-'*7}  {'-'*4}")

    for e in sorted(family['entries'], key=lambda x: -x['enrichment_1pct']):
        print(f"  {e['base']:6g}  {e['gcd_with_360']:4d}  {e['shared_prime_factors']:6d}  "
              f"{'yes' if e['divides_360'] else 'no':>6s}  "
              f"{'yes' if e['multiple_of_360'] else 'no':>4s}  "
              f"{e['enrichment_1pct']:7.2f}  {e['hits_1pct']:4d}")

    print()
    print(f"  Correlation(shared_prime_factors, enrichment_1%): r = {family['corr_shared_factors']:+.3f}")
    print(f"  Correlation(gcd_with_360, enrichment_1%):         r = {family['corr_gcd']:+.3f}")

    # Group by divisibility relationship
    divides = [e for e in family['entries'] if e['divides_360']]
    multiples = [e for e in family['entries'] if e['multiple_of_360']]
    neither = [e for e in family['entries'] if not e['divides_360'] and not e['multiple_of_360']]

    if divides:
        avg_div = sum(e['enrichment_1pct'] for e in divides) / len(divides)
        print(f"  Mean enrichment for divisors of 360:    {avg_div:.2f}x  ({len(divides)} bases)")
    if multiples:
        avg_mult = sum(e['enrichment_1pct'] for e in multiples) / len(multiples)
        print(f"  Mean enrichment for multiples of 360:   {avg_mult:.2f}x  ({len(multiples)} bases)")
    if neither:
        avg_neither = sum(e['enrichment_1pct'] for e in neither) / len(neither)
        print(f"  Mean enrichment for unrelated bases:    {avg_neither:.2f}x  ({len(neither)} bases)")

    # ===================================================================
    # 5. n_gap SELF-CONSISTENCY TEST
    # ===================================================================
    print()
    print("=" * 80)
    print("5. n_gap SELF-CONSISTENCY: dim(G_SM) → ENERGY HIERARCHY")
    print("=" * 80)
    print()

    sc = ngap_self_consistency()

    print(f"  Key hypothesis: n_gap = dim(G_SM) = dim(SU(3)) + dim(SU(2)) + dim(U(1)) = 8 + 3 + 1 = {sc['dim_SM']}")
    print(f"  This predicts:  M_Planck / M_GUT = φ^{sc['n_gap']} = {PHI**sc['n_gap']:.4f}")
    print()
    print(f"  n_gap (Planck-to-GUT exponent):")
    print(f"    Assumed:    {sc['n_gap']}")
    print(f"    Inferred:   {sc['n_gap_inferred']:.4f}")
    print(f"    Match:      {'yes' if abs(sc['n_gap_inferred'] - sc['n_gap']) < 0.5 else 'no (>0.5 off)'}")
    print()
    print(f"  n_total (Planck-to-mZ exponent):")
    print(f"    Predicted:  {sc['n_total']} (= {sc['n_GUT']} + {sc['n_gap']})")
    print(f"    Inferred:   {sc['n_total_inferred']:.4f}")
    print()
    print(f"  Forward predictions:")
    print(f"    M_GUT = m_Z × φ^{sc['n_GUT']}:")
    print(f"      predicted:  {sc['M_GUT_predicted']:.4e} GeV")
    print(f"      actual:     {sc['M_GUT_actual']:.4e} GeV")
    print(f"      rel error:  {sc['M_GUT_rel_err']:+.1%}")
    print()
    print(f"    M_Planck = M_GUT × φ^{sc['n_gap']}:")
    print(f"      predicted:  {sc['M_Pl_from_gap_predicted']:.4e} GeV")
    print(f"      actual:     {sc['M_Pl_actual']:.4e} GeV")
    print(f"      rel error:  {sc['M_Pl_gap_rel_err']:+.1%}")
    print()
    print(f"    M_Planck = m_Z × φ^{sc['n_total']}:")
    print(f"      predicted:  {sc['M_Pl_from_total_predicted']:.4e} GeV")
    print(f"      actual:     {sc['M_Pl_actual']:.4e} GeV")
    print(f"      rel error:  {sc['M_Pl_total_rel_err']:+.1%}")
    print()
    print(f"  Dimensional analysis:")
    print(f"    n_gap = 12 = dim(G_SM)  ← gauge group dimensionality")
    print(f"    n_GUT ≈ 70              ← EW-to-GUT desert")
    print(f"    n_total = 82 ≈ 12 × φ⁴ ← self-similar hierarchy")
    print(f"    F(12) = 144 = 12²       ← unique Fibonacci square relation")
    print(f"    n_gap = n_total / φ⁴ = {82/PHI**4:.2f}")
    print()

    # ===================================================================
    # COMBINED SUMMARY
    # ===================================================================
    print("=" * 80)
    print("DEEP ANALYSIS SUMMARY")
    print("=" * 80)
    print()

    print(f"  Address Table:  {hits_1}/{total} ratios hit <1% ({hits_1/total:.0%}), enrichment {hits_1/total/null_1:.2f}x")
    print(f"  CKM+PMNS:       {mx_h1}/{mx_total} params hit <1% ({mx_h1/mx_total:.0%})")
    print(f"  Clean lattice operations: {len(clean_ops)} mappings with <5% error")
    print(f"  360-family:     correlation with shared prime factors = {family['corr_shared_factors']:+.3f}")
    if abs(sc['n_gap_inferred'] - sc['n_gap']) < 0.5:
        print(f"  n_gap:          dim(G_SM) = 12 CONFIRMED as Planck-GUT exponent (inferred: {sc['n_gap_inferred']:.2f})")
    else:
        print(f"  n_gap:          dim(G_SM) = 12 vs inferred {sc['n_gap_inferred']:.2f} — discrepancy")
    print()

    return 0


def cmd_gut_predict(args: argparse.Namespace) -> int:
    """
    Extended predictions: running-coupling lattice energies, m_W cross-consistency,
    neutrino sector, G↔1/G duality, vacancy catalog.
    """
    import math as _math
    from physics_test.gut_validate import (
        build_sorted_lattice, nearest_lattice as nl,
        running_coupling_lattice_energies,
        mw_cross_consistency,
        duality_analysis,
        vacancy_catalog,
        neutrino_predictions,
        build_address_table,
        binomial_p_value,
        address_coincidences,
        m_value_zeckendorf_analysis,
        extended_hadron_predictions,
        m_value_distribution,
        lattice_coverage,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    m_range = list(range(-80, 120))

    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=360.0, include=include)
        for _, v in cs.items():
            all_Cs.append(float(v))
    all_Cs = sorted(set(all_Cs))

    lattice = build_sorted_lattice(all_Cs, m_range)

    # ===================================================================
    # 1. RUNNING COUPLING LATTICE-HIT ENERGIES
    # ===================================================================
    print("=" * 80)
    print("1. RUNNING COUPLING LATTICE-HIT ENERGIES")
    print("=" * 80)
    print()
    print("  At what energy Q does each gauge coupling exactly hit a lattice point?")
    print("  Using 1-loop SM beta functions from m_Z.")
    print()

    rg_hits = running_coupling_lattice_energies(all_Cs, m_range)

    # Group by coupling
    by_coupling: dict[str, list] = {}
    for h in rg_hits:
        if h.coupling not in by_coupling:
            by_coupling[h.coupling] = []
        by_coupling[h.coupling].append(h)

    for coupling in sorted(by_coupling.keys()):
        hits = by_coupling[coupling]
        print(f"  {coupling}:")
        print(f"  {'C':>6s}  {'m':>4s}  {'C/φ^m':>10s}  {'Q (GeV)':>14s}  {'log₁₀Q':>8s}  Note")
        print(f"  {'-'*6}  {'-'*4}  {'-'*10}  {'-'*14}  {'-'*8}  ----")

        notable_scales = [
            (1.0, "Λ_QCD"),
            (2.0, "charm threshold"),
            (5.0, "bottom threshold"),
            (80.4, "m_W"),
            (91.2, "m_Z"),
            (125.3, "m_H"),
            (172.8, "m_t"),
            (1e3, "1 TeV"),
            (1e4, "10 TeV"),
            (1e16, "M_GUT"),
            (1.2e19, "M_Planck"),
        ]

        for h in hits:
            note = ""
            for scale, label in notable_scales:
                if abs(_math.log10(h.Q_GeV) - _math.log10(scale)) < 0.15:
                    note = f"≈ {label}"
                    break
            if h.Q_GeV > 1e6:
                q_str = f"{h.Q_GeV:.3e}"
            else:
                q_str = f"{h.Q_GeV:.2f}"
            print(f"  {h.C:6g}  {h.m:4d}  {h.lattice_val:10.4f}  {q_str:>14s}  {_math.log10(h.Q_GeV):8.2f}  {note}")
        print()

    # Summary: crossings between couplings at lattice points
    print("  Lattice-point crossings (two couplings at same C):")
    for C in all_Cs:
        hits_at_C: dict[str, list] = {}
        for h in rg_hits:
            if h.C == C:
                if h.coupling not in hits_at_C:
                    hits_at_C[h.coupling] = []
                hits_at_C[h.coupling].append(h)
        couplings_present = [k for k in hits_at_C if hits_at_C[k]]
        if len(couplings_present) >= 2:
            print(f"    C={C:g}: {', '.join(couplings_present)}")
    print()

    # ===================================================================
    # 2. m_W CROSS-CONSISTENCY
    # ===================================================================
    print("=" * 80)
    print("2. m_W CROSS-CONSISTENCY FROM MULTIPLE LATTICE PATHS")
    print("=" * 80)
    print()

    mw_paths = mw_cross_consistency(lattice)

    print(f"  {'Path':<20s}  {'Address':>12s}  {'Lattice ratio':>14s}  {'m_W pred':>10s}  {'m_W meas':>10s}  {'rel%':>7s}")
    print(f"  {'-'*20}  {'-'*12}  {'-'*14}  {'-'*10}  {'-'*10}  {'-'*7}")

    mw_values = []
    for p in mw_paths:
        pct = p['rel_err'] * 100.0
        marker = "***" if abs(pct) < 1.0 else "** " if abs(pct) < 3.0 else "*  " if abs(pct) < 5.0 else "   "
        print(f"  {p['path']:<20s}  {p['addr']:>12s}  {p['lattice_ratio']:14.6f}  {p['mW_predicted']:10.3f}  {p['mW_measured']:10.3f}  {pct:+6.2f}% {marker}")
        if _math.isfinite(p['mW_predicted']):
            mw_values.append(p['mW_predicted'])

    if mw_values:
        mean_mw = sum(mw_values) / len(mw_values)
        spread = max(mw_values) - min(mw_values)
        print()
        print(f"  Lattice m_W predictions: {', '.join(f'{v:.2f}' for v in mw_values)} GeV")
        print(f"  Mean: {mean_mw:.3f} GeV,  Spread: {spread:.3f} GeV")
        print(f"  Measured (PDG): 80.377 ± 0.012 GeV")
        print(f"  Measured (CDF): 80.434 ± 0.009 GeV")
    print()

    # ===================================================================
    # 3. NEUTRINO SECTOR
    # ===================================================================
    print("=" * 80)
    print("3. NEUTRINO SECTOR PREDICTIONS")
    print("=" * 80)
    print()

    nu = neutrino_predictions(lattice)

    print(f"  {'Name':<22s}  {'Value':>14s}  {'C':>6s}  {'m':>4s}  {'C/φ^m':>14s}  {'rel%':>7s}  Note")
    print(f"  {'-'*22}  {'-'*14}  {'-'*6}  {'-'*4}  {'-'*14}  {'-'*7}  {'-'*20}")

    for r in nu:
        pct = r['rel_err'] * 100.0
        marker = "***" if abs(pct) < 1.0 else "** " if abs(pct) < 3.0 else "*  " if abs(pct) < 5.0 else "   "
        if r['value'] > 1e4 or r['value'] < 1e-4:
            val_str = f"{r['value']:.4e}"
            lat_str = f"{r['lattice_val']:.4e}"
        else:
            val_str = f"{r['value']:.6f}"
            lat_str = f"{r['lattice_val']:.6f}"
        print(f"  {r['name']:<22s}  {val_str:>14s}  {r['C']:6g}  {r['m']:4d}  {lat_str:>14s}  {pct:+6.2f}% {marker} {r['note']}")
    print()

    # ===================================================================
    # 4. G ↔ 1/G DUALITY ANALYSIS
    # ===================================================================
    print("=" * 80)
    print("4. G ↔ 1/G DUALITY ANALYSIS")
    print("=" * 80)
    print()
    print("  For each ratio G → (C₁, m₁), what does 1/G → (C₂, m₂)?")
    print("  If the lattice has an inversion symmetry, ΔC and Δm should follow patterns.")
    print()

    dual = duality_analysis(lattice)

    print(f"  {'Name':<12s}  {'G addr':>10s}  {'err%':>6s}  {'1/G addr':>10s}  {'err%':>6s}  {'ΔC':>6s}  {'Δm':>4s}  {'C₁×C₂':>8s}  Clean?")
    print(f"  {'-'*12}  {'-'*10}  {'-'*6}  {'-'*10}  {'-'*6}  {'-'*6}  {'-'*4}  {'-'*8}  ------")

    product_counts: dict[float, int] = {}
    for d in dual:
        pct1 = d['rel_err'] * 100.0
        pct2 = d['inv_rel_err'] * 100.0
        clean = "YES" if d['both_clean'] else "no"
        prod = d['product_C']
        product_counts[prod] = product_counts.get(prod, 0) + 1
        print(f"  {d['name']:<12s}  {d['addr']:>10s}  {pct1:+5.1f}%  {d['inv_addr']:>10s}  {pct2:+5.1f}%  {d['delta_C']:+6g}  {d['delta_m']:+4d}  {prod:8g}  {clean}")

    print()
    print(f"  C₁ × C₂ distribution (inversion product):")
    for prod in sorted(product_counts.keys()):
        print(f"    C₁×C₂ = {prod:g}: {product_counts[prod]} ratios")

    # Check for m1 + m2 pattern
    m_sums = [d['m1'] + d['m2'] for d in dual if d['both_clean']]
    if m_sums:
        print(f"  m₁ + m₂ for clean pairs: {m_sums}")
    print()

    # ===================================================================
    # 5. VACANCY CATALOG
    # ===================================================================
    print("=" * 80)
    print("5. VACANCY CATALOG: EMPTY LATTICE SITES")
    print("=" * 80)
    print()

    addr_table = build_address_table(lattice)
    vacancies = vacancy_catalog(lattice, addr_table, m_range_display=range(-15, 25))

    filled = [v for v in vacancies if v['matches']]
    empty_count = sum(1 for v in vacancies if not v['matches'])

    print(f"  Scanning {len(vacancies)} lattice sites in m ∈ [-15, 25)")
    print(f"  Vacancies with candidate matches (<5% of a known constant):")
    print()
    print(f"  {'C':>6s}  {'m':>4s}  {'C/φ^m':>14s}  {'Candidate':>22s}  {'Value':>14s}  {'rel%':>7s}  Note")
    print(f"  {'-'*6}  {'-'*4}  {'-'*14}  {'-'*22}  {'-'*14}  {'-'*7}  {'-'*20}")

    for v in filled:
        for match in v['matches']:
            pct = match['rel_err'] * 100.0
            if v['lattice_val'] > 1e4:
                lv_str = f"{v['lattice_val']:.3e}"
                mv_str = f"{match['value']:.3e}"
            elif v['lattice_val'] < 0.001:
                lv_str = f"{v['lattice_val']:.5e}"
                mv_str = f"{match['value']:.5e}"
            else:
                lv_str = f"{v['lattice_val']:.6f}"
                mv_str = f"{match['value']:.6f}"
            print(f"  {v['C']:6g}  {v['m']:4d}  {lv_str:>14s}  {match['name']:>22s}  {mv_str:>14s}  {pct:+6.2f}%  {match['note']}")

    print()
    print(f"  Total vacancies scanned: {len(vacancies)}")
    print(f"  With candidate matches: {len(filled)}")
    print(f"  Truly empty: {empty_count}")
    print()

    # ===================================================================
    # 6. BINOMIAL SIGNIFICANCE
    # ===================================================================
    print("=" * 80)
    print("6. BINOMIAL SIGNIFICANCE TEST")
    print("=" * 80)
    print()

    null_1 = lattice_coverage(lattice, 0.01, 5000)
    null_3 = lattice_coverage(lattice, 0.03, 5000)
    null_5 = lattice_coverage(lattice, 0.05, 5000)

    for label, threshold, null_rate in [("1%", 0.01, null_1), ("3%", 0.03, null_3), ("5%", 0.05, null_5)]:
        hits = sum(1 for row in addr_table if abs(row['rel_err']) < threshold)
        n = len(addr_table)
        pval = binomial_p_value(hits, n, null_rate)
        sigma = -_math.log10(pval) if pval > 0 else float("inf")
        print(f"  Tolerance < {label}:")
        print(f"    Observed: {hits}/{n} ({hits/n:.0%})")
        print(f"    Null rate: {null_rate:.1%}")
        print(f"    P(X ≥ {hits} | Binomial(n={n}, p={null_rate:.4f})): {pval:.2e}")
        print(f"    Significance: {sigma:.1f} orders of magnitude")
        print()

    # ===================================================================
    # 7. ADDRESS COINCIDENCES
    # ===================================================================
    print("=" * 80)
    print("7. ADDRESS COINCIDENCE ANALYSIS")
    print("=" * 80)
    print()

    coinc = address_coincidences(addr_table)
    print(f"  {coinc['n_ratios']} ratios mapped to ~{coinc['L_effective']} effective lattice sites")
    print(f"  Expected collisions (birthday): {coinc['expected_collisions']:.1f}")
    print(f"  Observed collisions: {coinc['observed_collisions']}")
    print()

    if coinc['shared_addresses']:
        print(f"  Shared addresses:")
        for (C, m), rows in sorted(coinc['shared_addresses'].items()):
            names = [r['name'] for r in rows]
            sectors = list(set(r.get('sector', '?') for r in rows))
            cross = len(sectors) > 1
            errs = [f"{r['rel_err']*100:+.1f}%" for r in rows]
            flag = " ← CROSS-SECTOR" if cross else ""
            print(f"    ({C:g}, {m:2d}): {', '.join(names)} [{', '.join(errs)}]{flag}")
            if cross:
                print(f"             sectors: {', '.join(sectors)}")
        print()

    # ===================================================================
    # 8. EXTENDED HADRON MASS RATIOS
    # ===================================================================
    print("=" * 80)
    print("8. EXTENDED HADRON/MESON/BARYON MASS RATIOS")
    print("=" * 80)
    print()

    hadrons = extended_hadron_predictions(lattice)
    had_h1 = had_h3 = had_h5 = 0

    print(f"  {'Name':<16s}  {'Value':>10s}  {'Sector':<14s}  {'C':>6s}  {'m':>4s}  {'C/φ^m':>10s}  {'rel%':>7s}")
    print(f"  {'-'*16}  {'-'*10}  {'-'*14}  {'-'*6}  {'-'*4}  {'-'*10}  {'-'*7}")

    for h in hadrons:
        pct = h['rel_err'] * 100.0
        marker = "***" if abs(pct) < 1.0 else "** " if abs(pct) < 3.0 else "*  " if abs(pct) < 5.0 else "   "
        if abs(pct) < 1.0: had_h1 += 1
        if abs(pct) < 3.0: had_h3 += 1
        if abs(pct) < 5.0: had_h5 += 1
        print(f"  {h['name']:<16s}  {h['value']:10.4f}  {h['sector']:<14s}  {h['C']:6g}  {h['m']:4d}  {h['lattice_val']:10.4f}  {pct:+6.2f}% {marker}")

    had_total = len(hadrons)
    had_null_1 = lattice_coverage(lattice, 0.01, 3000)
    print()
    print(f"  Hadron hit rates ({had_total} ratios):")
    print(f"    < 1%: {had_h1}/{had_total} ({had_h1/had_total:.0%}), null = {had_null_1:.1%}, enrichment = {had_h1/had_total/had_null_1:.2f}x")
    print(f"    < 3%: {had_h3}/{had_total} ({had_h3/had_total:.0%})")
    print(f"    < 5%: {had_h5}/{had_total} ({had_h5/had_total:.0%})")
    print()

    # ===================================================================
    # 9. ZECKENDORF DECOMPOSITION OF m-VALUES
    # ===================================================================
    print("=" * 80)
    print("9. ZECKENDORF DECOMPOSITION OF m-VALUES")
    print("=" * 80)
    print()

    zeck = m_value_zeckendorf_analysis(addr_table)

    fib_count = sum(1 for z in zeck if z['is_fibonacci'])
    lucas_count = sum(1 for z in zeck if z['is_lucas'])

    print(f"  {'Name':<20s}  {'C':>6s}  {'m':>4s}  {'|m|':>4s}  {'Fib?':>5s}  {'Luc?':>5s}  Zeckendorf")
    print(f"  {'-'*20}  {'-'*6}  {'-'*4}  {'-'*4}  {'-'*5}  {'-'*5}  {'-'*20}")

    for z in sorted(zeck, key=lambda x: (x['C'], x['m'])):
        fib_str = "YES" if z['is_fibonacci'] else ""
        luc_str = "YES" if z['is_lucas'] else ""
        zeck_str = " + ".join(str(f) for f in z['zeckendorf'])
        print(f"  {z['name']:<20s}  {z['C']:6g}  {z['m']:4d}  {z['abs_m']:4d}  {fib_str:>5s}  {luc_str:>5s}  {zeck_str}")

    print()
    print(f"  Fibonacci m-values: {fib_count}/{len(zeck)} ({fib_count/len(zeck):.0%})")
    print(f"  Lucas m-values:     {lucas_count}/{len(zeck)} ({lucas_count/len(zeck):.0%})")

    # Zeckendorf length distribution
    len_counts: dict[int, int] = {}
    for z in zeck:
        l = z['zeckendorf_len']
        len_counts[l] = len_counts.get(l, 0) + 1
    print(f"  Zeckendorf length distribution: {dict(sorted(len_counts.items()))}")
    print()

    # ===================================================================
    # 10. m-VALUE DISTRIBUTION ANALYSIS
    # ===================================================================
    print("=" * 80)
    print("10. m-VALUE DISTRIBUTION ANALYSIS")
    print("=" * 80)
    print()

    mdist = m_value_distribution(addr_table)

    print(f"  {mdist['total_ratios']} ratios occupy {mdist['distinct_m']} distinct m-values")
    print()
    print(f"  Most populated m-values:")
    for m_val, count in mdist['hot_m']:
        names_at_m = [row['name'] for row in addr_table if row['m'] == m_val]
        print(f"    m = {m_val:4d}: {count} ratios — {', '.join(names_at_m[:5])}" + ("..." if len(names_at_m) > 5 else ""))

    print()
    print(f"  Band-by-band gap analysis:")
    fibs_list = [1, 2, 3, 5, 8, 13, 21, 34]
    for C in sorted(mdist['band_stats'].keys()):
        bs = mdist['band_stats'][C]
        if bs['n'] < 2:
            continue
        frac = bs['fib_gap_fraction']
        print(f"    C={C:6g}: {bs['n']:2d} values, m ∈ [{bs['m_range'][0]}, {bs['m_range'][1]}], "
              f"gaps = {bs['gaps']}, Fibonacci fraction = {frac:.0%}")

    print()

    # ===================================================================
    # 11. CROSS-BAND HORIZONTAL STRUCTURE
    # ===================================================================
    print("=" * 80)
    print("11. CROSS-BAND HORIZONTAL STRUCTURE")
    print("=" * 80)
    print()

    from physics_test.gut_validate import (
        cross_band_analysis,
        lattice_implied_relations,
        weinberg_angle_lattice_running,
        cosmological_constant_analysis,
        c_band_catalog,
    )

    horiz = cross_band_analysis(addr_table)
    print(f"  m-values occupied by 2+ ratios from different C-bands:")
    print()
    for m in sorted(horiz.keys()):
        h = horiz[m]
        if not h['cross_band']:
            continue
        cs_str = ", ".join(f"C={c:g}" for c in h['Cs'])
        cross = " ← CROSS-SECTOR" if h['cross_sector'] else ""
        print(f"  m = {m:4d}: {h['count']} ratios across [{cs_str}]{cross}")
        for e in h['entries']:
            pct = e['rel_err'] * 100
            print(f"    ({e['C']:g},{m:2d}) {e['name']:<20s}  sector={e['sector']:<10s}  err={pct:+.1f}%")

        for cr in h['c_ratios']:
            print(f"    → {cr[0]}/{cr[1]}: C ratio = {cr[2]:g}/{cr[3]:g} = {cr[4]:.4f}")
        print()

    # ===================================================================
    # 12. LATTICE-IMPLIED MASS RELATIONS
    # ===================================================================
    print("=" * 80)
    print("12. LATTICE-IMPLIED MASS RELATIONS")
    print("=" * 80)
    print()

    rels = lattice_implied_relations(addr_table)

    print(f"  Same-address relations (implied equality):")
    same_addr = [r for r in rels if r['type'] == 'same_address']
    for r in same_addr[:10]:
        pct = r['rel_err'] * 100
        print(f"    {r['implied_relation']:<40s}  actual ratio = {r['actual_ratio']:.6f}  ({pct:+.2f}%)")

    print()
    print(f"  Same-C, small Δm relations (implied φ-power ratio):")
    dm_rels = [r for r in rels if r['type'].startswith('same_C_dm') and abs(r['rel_err']) < 0.10]
    dm_rels.sort(key=lambda r: abs(r['rel_err']))
    for r in dm_rels[:15]:
        pct = r['rel_err'] * 100
        print(f"    {r['implied_relation']:<55s}  actual = {r['actual_ratio']:.4f}  ({pct:+.2f}%)")

    print()
    print(f"  Same-m, cross-band relations (implied C-ratio):")
    same_m = [r for r in rels if r['type'] == 'same_m' and abs(r['rel_err']) < 0.10]
    same_m.sort(key=lambda r: abs(r['rel_err']))
    for r in same_m[:15]:
        pct = r['rel_err'] * 100
        print(f"    {r['implied_relation']:<55s}  actual = {r['actual_ratio']:.4f}  ({pct:+.2f}%)")
    print()

    # ===================================================================
    # 13. WEINBERG ANGLE RUNNING ON THE LATTICE
    # ===================================================================
    print("=" * 80)
    print("13. WEINBERG ANGLE RG RUNNING ON THE LATTICE")
    print("=" * 80)
    print()

    weinberg = weinberg_angle_lattice_running(lattice)

    print(f"  {'Scale':<28s}  {'sin²θ_W':>10s}  {'C':>6s}  {'m':>4s}  {'C/φ^m':>10s}  {'rel%':>7s}")
    print(f"  {'-'*28}  {'-'*10}  {'-'*6}  {'-'*4}  {'-'*10}  {'-'*7}")

    for w in weinberg:
        pct = w['rel_err'] * 100
        marker = "***" if abs(pct) < 1 else "** " if abs(pct) < 3 else "*  " if abs(pct) < 5 else "   "
        print(f"  {w['label']:<28s}  {w['sin2w']:10.5f}  {w['C']:6g}  {w['m']:4d}  {w['lattice_val']:10.5f}  {pct:+6.2f}% {marker}")

    # Check if address changes with scale
    addrs = [(w['C'], w['m']) for w in weinberg]
    unique_addrs = list(set(addrs))
    print()
    if len(unique_addrs) == 1:
        print(f"  All scales map to the same address: ({unique_addrs[0][0]:g}, {unique_addrs[0][1]})")
    else:
        print(f"  Address changes with scale: {len(unique_addrs)} distinct addresses")
        for addr in unique_addrs:
            scales = [w['label'] for w in weinberg if (w['C'], w['m']) == addr]
            print(f"    ({addr[0]:g}, {addr[1]}): {', '.join(scales)}")
    print()

    # ===================================================================
    # 14. COSMOLOGICAL CONSTANT AND DARK SECTOR
    # ===================================================================
    print("=" * 80)
    print("14. COSMOLOGICAL CONSTANT AND DARK SECTOR")
    print("=" * 80)
    print()

    cosmo = cosmological_constant_analysis(lattice)

    print(f"  {'Name':<22s}  {'Value':>14s}  {'C':>6s}  {'m':>4s}  {'C/φ^m':>14s}  {'rel%':>7s}  Note")
    print(f"  {'-'*22}  {'-'*14}  {'-'*6}  {'-'*4}  {'-'*14}  {'-'*7}  {'-'*20}")

    for c in cosmo:
        pct = c['rel_err'] * 100
        marker = "***" if abs(pct) < 1 else "** " if abs(pct) < 3 else "*  " if abs(pct) < 5 else "   "
        if c['value'] > 1e4:
            val_str = f"{c['value']:.4e}"
            lat_str = f"{c['lattice_val']:.4e}"
        elif c['value'] < 1e-4:
            val_str = f"{c['value']:.4e}"
            lat_str = f"{c['lattice_val']:.4e}"
        else:
            val_str = f"{c['value']:.6f}"
            lat_str = f"{c['lattice_val']:.6f}"
        print(f"  {c['name']:<22s}  {val_str:>14s}  {c['C']:6g}  {c['m']:4d}  {lat_str:>14s}  {pct:+6.2f}% {marker} {c['note']}")
    print()

    # ===================================================================
    # 15. C-BAND PHYSICS CATALOG
    # ===================================================================
    print("=" * 80)
    print("15. C-BAND PHYSICS CATALOG")
    print("=" * 80)
    print()

    catalog = c_band_catalog(addr_table)

    for C in sorted(catalog.keys()):
        band = catalog[C]
        print(f"  ┌─ C = {C:g}  ({band['origin']})")
        print(f"  │  Role: {band['physics']}")
        print(f"  │  Occupants: {band['n_occupants']}, Sectors: {', '.join(band['sectors'])}")
        print(f"  │  Accuracy: {band['hits_1pct']}/{band['n_occupants']} < 1%, {band['hits_3pct']}/{band['n_occupants']} < 3%")
        print(f"  │")
        for occ in band['occupants']:
            pct = occ['rel_err'] * 100
            marker = "***" if abs(pct) < 1 else "** " if abs(pct) < 3 else "*  " if abs(pct) < 5 else "   "
            print(f"  │  m={occ['m']:4d}  {occ['name']:<20s}  [{occ['sector']:<10s}]  {pct:+6.2f}% {marker}")
        print(f"  └{'─'*70}")
        print()

    return 0


def cmd_gut_explore(args: argparse.Namespace) -> int:
    """
    Round 4 explorations: lattice arithmetic predictions, duality axis analysis,
    α_s comparison to data, lattice thermodynamics, 360 number theory.
    """
    import math as _math
    from physics_test.gut_validate import (
        build_sorted_lattice, nearest_lattice as nl,
        lattice_arithmetic_predictions,
        duality_axis_analysis,
        alpha_s_comparison,
        lattice_thermodynamics,
        base_360_number_theory,
        lattice_selection_rules,
        lattice_closure_analysis,
        error_budget_analysis,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    m_range = list(range(-80, 120))

    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=360.0, include=include)
        for _, v in cs.items():
            all_Cs.append(v)
    Cs = sorted(set(all_Cs))
    lattice = build_sorted_lattice(Cs, m_range)

    # ===================================================================
    # 1. LATTICE ARITHMETIC PREDICTIONS
    # ===================================================================
    print("=" * 80)
    print("1. LATTICE ARITHMETIC PREDICTIONS")
    print("=" * 80)
    print()

    result = lattice_arithmetic_predictions(lattice)

    print("  A) φ-power and C-ratio predictions from confirmed addresses")
    print("  " + "-" * 70)
    print(f"  {'Relation':<55s} {'Predicted':>10s} {'Actual':>10s} {'Error':>8s}")
    print("  " + "-" * 70)
    for p in result['lattice_relations'][:30]:
        pct = p['rel_err'] * 100
        marker = " ***" if abs(pct) < 0.1 else " ** " if abs(pct) < 1 else " *  " if abs(pct) < 3 else "    "
        print(f"  {p['relation']:<55s} {p['predicted']:10.4f} {p['actual']:10.4f} {pct:+7.2f}%{marker}")

    n_sub01 = sum(1 for p in result['lattice_relations'] if abs(p['rel_err']) < 0.001)
    n_sub1 = sum(1 for p in result['lattice_relations'] if abs(p['rel_err']) < 0.01)
    n_total = len(result['lattice_relations'])
    print()
    print(f"  Total: {n_total} relations, {n_sub01} < 0.1%, {n_sub1} < 1%")
    print()

    print("  B) Novel predictions (quantities predicted from lattice addresses)")
    print("  " + "-" * 70)
    for np_item in result['novel_predictions']:
        pred = np_item['predicted_value']
        meas = np_item['measured_value']
        rel = (pred - meas) / meas * 100
        unit = np_item.get('unit', '')
        print(f"  {np_item['name']}")
        print(f"    Formula:   {np_item['formula']}")
        print(f"    Predicted: {pred:.6g} {unit}")
        print(f"    Measured:  {meas:.6g} {unit}")
        print(f"    Error:     {rel:+.2f}%")
        print()

    # ===================================================================
    # 2. DUALITY AXIS ANALYSIS: m₁+m₂=23
    # ===================================================================
    print("=" * 80)
    print("2. DUALITY AXIS ANALYSIS: m₁ + m₂ = 23")
    print("=" * 80)
    print()

    da = duality_axis_analysis()
    print(f"  Duality sum:  {da['duality_sum']}")
    print(f"  dim(SU(5)):   {da['dim_SU5']}")
    print(f"  Is prime:     {da['is_prime']}")
    print(f"  Is Fibonacci: {da['is_fibonacci']}")
    print()
    print("  Connections:")
    for c in da['connections']:
        print(f"    • {c}")
    print()
    print("  Predictions for other GUT groups:")
    for group, pred in da['gut_predictions'].items():
        print(f"    {group}: duality sum → {pred}")
    print()

    # ===================================================================
    # 3. α_s COMPARISON TO EXPERIMENTAL DATA
    # ===================================================================
    print("=" * 80)
    print("3. α_s COMPARISON TO EXPERIMENTAL DATA AT MULTIPLE SCALES")
    print("=" * 80)
    print()

    alpha_results = alpha_s_comparison(lattice)

    print(f"  {'Measurement':<25s} {'Q [GeV]':>8s} {'α_s(exp)':>10s} {'±':>6s} "
          f"{'(C,m)':>8s} {'α_s(lat)':>10s} {'Err':>7s} {'σ':>5s}")
    print("  " + "-" * 85)
    for r in alpha_results:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        sig = r['sigma_tension']
        marker = " ***" if sig < 1 else " ** " if sig < 2 else " *  " if sig < 3 else "    "
        print(f"  {r['label']:<25s} {r['Q_GeV']:8.1f} {r['alpha_s_measured']:10.4f} "
              f"{r['uncertainty']:6.4f} {addr:>8s} {r['lattice_alpha_s']:10.4f} "
              f"{pct:+6.2f}% {sig:4.1f}σ{marker}")
    print()

    n_1sig = sum(1 for r in alpha_results if r['sigma_tension'] < 1)
    n_2sig = sum(1 for r in alpha_results if r['sigma_tension'] < 2)
    print(f"  Within 1σ: {n_1sig}/{len(alpha_results)}")
    print(f"  Within 2σ: {n_2sig}/{len(alpha_results)}")
    print()

    # ===================================================================
    # 4. LATTICE THERMODYNAMICS
    # ===================================================================
    print("=" * 80)
    print("4. LATTICE THERMODYNAMICS")
    print("=" * 80)
    print()

    thermo = lattice_thermodynamics(Cs)
    print(f"  Number of levels in core range: {thermo['n_levels']}")
    print(f"  Average level spacing (in log g): {thermo['avg_spacing']:.4f}")
    print(f"  Minimum spacing: {thermo['min_spacing']:.6f}")
    print(f"  log(φ) = {thermo['log_phi']:.6f}")
    ratio = thermo['avg_spacing'] / thermo['log_phi'] if thermo['log_phi'] > 0 else 0
    print(f"  Avg spacing / log(φ) = {ratio:.4f}")
    print()

    print(f"  {'β':>8s} {'Z':>12s} {'⟨log g⟩':>10s} {'⟨g⟩':>12s} {'Entropy':>10s}")
    print("  " + "-" * 60)
    for bs in thermo['beta_scan']:
        print(f"  {bs['beta']:8.2f} {bs['Z']:12.2f} {bs['avg_log_g']:10.4f} "
              f"{bs['avg_g']:12.4g} {bs['entropy']:10.4f}")
    print()

    # Identify "natural temperature"
    # <log g> ≈ -4.92 (log 1/137) for EM
    # <log g> ≈ -2.14 (log 0.118) for QCD
    log_1_over_alpha = -_math.log(137.036)
    log_alpha_s = _math.log(0.1179)
    print(f"  Target ⟨log g⟩ for EM:  {log_1_over_alpha:.4f} (log(1/137))")
    print(f"  Target ⟨log g⟩ for QCD: {log_alpha_s:.4f} (log(0.118))")
    print()

    # ===================================================================
    # 5. 360 NUMBER THEORY
    # ===================================================================
    print("=" * 80)
    print("5. 360 NUMBER THEORY")
    print("=" * 80)
    print()

    nt = base_360_number_theory()
    print(f"  Factorization: {nt['factorization']}")
    print(f"  Number of divisors τ(360) = {nt['n_divisors']}")
    print(f"  τ(360) = dim(SU(5)) = 24: {nt['tau_equals_dim_SU5']}")
    print(f"  {nt['totient_relation']}")
    print(f"  Sum of divisors σ(360) = {nt['sigma']}")
    print()

    print("  All 24 divisors of 360:")
    divs = nt['all_divisors']
    for row_start in range(0, len(divs), 8):
        chunk = divs[row_start:row_start+8]
        print("    " + "  ".join(f"{d:4d}" for d in chunk))
    print()

    print("  C values as 360/X:")
    for ca in nt['c_analysis']:
        print(f"    C = {ca['C']:4g}  →  360/{ca['complement']}  (complement = {ca['complement']})")
    print()

    print(f"  C complements used: {nt['divisors_used']}")
    print(f"  All C-values divide 360: {nt['C_all_divisors']}")
    print()

    print("  Contexts:")
    for name, val, category in nt['contexts']:
        print(f"    [{category:>14s}]  {name}: {val}")
    print()

    print("  360 mod primes:")
    for p, r in nt['mod_analysis'].items():
        print(f"    360 mod {p:2d} = {r}")
    print()

    # ===================================================================
    # 6. LATTICE SELECTION RULES FOR PARTICLE DECAYS
    # ===================================================================
    print("=" * 80)
    print("6. LATTICE SELECTION RULES FOR PARTICLE DECAYS")
    print("=" * 80)
    print()

    sr = lattice_selection_rules()
    for i, d in enumerate(sr['decays'], 1):
        print(f"  {i}. {d['decay']}")
        print(f"     Parent:  {d['parent']}")
        if 'products' in d:
            for prod in d['products']:
                print(f"     Product: {prod}")
        if 'delta_m' in d:
            print(f"     Δm = {d['delta_m']}")
        if 'mass_ratio' in d:
            print(f"     Mass ratio = {d['mass_ratio']:.4f}")
        print(f"     Notes: {d['notes']}")
        print()

    print("  EMERGENT RULES:")
    for rule in sr['rules_summary']:
        print(f"    → {rule}")
    print()

    # ===================================================================
    # 7. LATTICE PRODUCT/QUOTIENT CLOSURE
    # ===================================================================
    print("=" * 80)
    print("7. LATTICE PRODUCT/QUOTIENT CLOSURE")
    print("=" * 80)
    print()

    closure = lattice_closure_analysis(lattice)

    print("  A) Products of lattice quantities (nearest lattice hit)")
    print("  " + "-" * 70)
    print(f"  {'Pair':<35s} {'Value':>12s} {'(C,m)':>10s} {'Lattice':>12s} {'Error':>8s}")
    print("  " + "-" * 70)
    for p in closure['products']:
        addr = f"({p['C_actual']:g},{p['m_actual']})"
        pct = p['rel_err'] * 100
        alg = f"[{p['C_algebraic']:g},{p['m_algebraic']}]"
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {p['pair']:<35s} {p['value']:12.4g} {addr:>10s} {p['lattice_val']:12.4g} {pct:+7.2f}%{marker}")
    print()

    print("  B) Quotients of lattice quantities")
    print("  " + "-" * 70)
    print(f"  {'Pair':<35s} {'Value':>12s} {'(C,m)':>10s} {'Lattice':>12s} {'Error':>8s}")
    print("  " + "-" * 70)
    for q in closure['quotients']:
        addr = f"({q['C_actual']:g},{q['m_actual']})"
        pct = q['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {q['pair']:<35s} {q['value']:12.4g} {addr:>10s} {q['lattice_val']:12.4g} {pct:+7.2f}%{marker}")
    print()

    print("  C) C-algebra: which C-products map to lattice via φ-shift")
    n_alg = 0
    for key in sorted(closure['c_algebra'].keys()):
        entry = closure['c_algebra'][key]
        print(f"    {key:>10s} = {entry['C_prod']:>6d}  →  {entry['maps_to_C']:g} × φ^{entry['phi_power']}")
        n_alg += 1
    if n_alg == 0:
        print("    (none found with exact φ-power mapping)")
    print()

    # ===================================================================
    # 8. ERROR BUDGET: INDEPENDENT VS DERIVED
    # ===================================================================
    print("=" * 80)
    print("8. ERROR BUDGET: INDEPENDENT VS DERIVED PREDICTIONS")
    print("=" * 80)
    print()

    eb = error_budget_analysis()
    print(f"  Independent input quantities: {eb['n_independent_inputs']}")
    print(f"  Sub-1% error:  {eb['n_sub_1pct']}")
    print(f"  Sub-0.5% error: {eb['n_sub_05pct']}")
    print(f"  Unique (C,m) addresses: {eb['effective_dof']}")
    print(f"  m-value span: {eb['m_span']} (from {min(eb['m_values'])} to {max(eb['m_values'])})")
    print()

    print("  C-band distribution:")
    for c in sorted(eb['c_distribution'].keys()):
        n = eb['c_distribution'][c]
        bar = "█" * n
        print(f"    C={c:4g}: {n:2d}  {bar}")
    print()

    print("  All independent entries:")
    print(f"  {'Name':<25s} {'(C,m)':>10s} {'Value':>12s} {'Error%':>8s}")
    print("  " + "-" * 60)
    for e in sorted(eb['entries'], key=lambda x: x[4]):
        addr = f"({e[1]:g},{e[2]})"
        print(f"  {e[0]:<25s} {addr:>10s} {e[3]:12.4g} {e[4]:+7.2f}%")
    print()

    print("  Derived relations (not independent predictions):")
    for name, rtype, note in eb['derived_relations']:
        print(f"    {name}")
        print(f"      Type: {rtype}")
        print(f"      {note}")
        print()

    return 0


def cmd_gut_spectrum(args: argparse.Namespace) -> int:
    """
    Round 5: fermion mass spectrum, symmetry breaking, new ratios, RG flow, outliers.
    """
    import math as _math
    from physics_test.gut_validate import (
        build_sorted_lattice, nearest_lattice as nl,
        fermion_mass_spectrum,
        symmetry_breaking_analysis,
        new_ratio_predictions,
        rg_flow_on_lattice,
        outlier_analysis,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    m_range = list(range(-80, 120))

    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=360.0, include=include)
        for _, v in cs.items():
            all_Cs.append(v)
    Cs = sorted(set(all_Cs))
    lattice = build_sorted_lattice(Cs, m_range)

    # ===================================================================
    # 1. COMPLETE FERMION MASS SPECTRUM
    # ===================================================================
    print("=" * 80)
    print("1. COMPLETE FERMION MASS SPECTRUM")
    print("=" * 80)
    print()

    spec = fermion_mass_spectrum(lattice)

    print("  A) All fermion masses as m_X/m_e → lattice address")
    print("  " + "-" * 75)
    print(f"  {'Fermion':<8s} {'Mass [GeV]':>12s} {'m_X/m_e':>12s} "
          f"{'(C,m)':>10s} {'Lattice':>12s} {'Pred [GeV]':>12s} {'Error':>8s}")
    print("  " + "-" * 75)
    for r in sorted(spec['spectrum'], key=lambda x: -x['mass_GeV']):
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else " *  " if abs(pct) < 5 else "    "
        print(f"  {r['fermion']:<8s} {r['mass_GeV']:12.6g} {r['ratio_to_me']:12.4g} "
              f"{addr:>10s} {r['lattice_val']:12.4g} {r['predicted_mass']:12.6g} {pct:+7.2f}%{marker}")
    print()

    print("  B) Inter-generation mass ratios (same charge)")
    print("  " + "-" * 65)
    print(f"  {'Ratio':<8s} {'Value':>10s} {'(C,m)':>10s} {'Lattice':>10s} {'Error':>8s}")
    print("  " + "-" * 65)
    for r in spec['generation_ratios']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {r['ratio']:<8s} {r['value']:10.4f} {addr:>10s} {r['lattice_val']:10.4f} {pct:+7.2f}%{marker}")
    print()

    print("  C) Cross-sector ratios (quark/lepton)")
    print("  " + "-" * 65)
    print(f"  {'Ratio':<8s} {'Value':>10s} {'(C,m)':>10s} {'Lattice':>10s} {'Error':>8s}")
    print("  " + "-" * 65)
    for r in spec['cross_sector_ratios']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {r['ratio']:<8s} {r['value']:10.4f} {addr:>10s} {r['lattice_val']:10.4f} {pct:+7.2f}%{marker}")
    print()

    print("  D) Neutrino masses on the lattice")
    print("  " + "-" * 70)
    for r in spec['neutrinos']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        print(f"  {r['fermion']:<4s}  mass = {r['mass_eV']:.4g} eV → "
              f"{addr}  lattice = {r['predicted_mass_eV']:.4g} eV  ({pct:+.2f}%)")
    print()

    h = spec['total_hierarchy']
    print(f"  Total hierarchy m_t/m_ν₂ = {h['value']:.3g} → ({h['C']:g},{h['m']})  "
          f"lattice = {h['lattice_val']:.3g}  ({h['rel_err']*100:+.2f}%)")
    print()

    # ===================================================================
    # 2. SU(5) SYMMETRY BREAKING ON THE LATTICE
    # ===================================================================
    print("=" * 80)
    print("2. SU(5) SYMMETRY BREAKING ON THE LATTICE")
    print("=" * 80)
    print()

    sb = symmetry_breaking_analysis(lattice)

    print("  Breaking chain: SU(5) → SU(3) × SU(2) × U(1)")
    print(f"  Broken generators: {sb['dim_broken']} (X and Y gauge bosons)")
    print()

    print("  C-ratios encode group dimension ratios:")
    for name, val in sb['c_ratios'].items():
        print(f"    {name:<20s} = {val:.4g}")
    print()

    print("  Key identity: C_A/C_B = dim(B)/dim(A)")
    print("  → The INVERSE of group dimension ratio determines the C-band ratio")
    print()

    print("  Branching rules and lattice sector mapping:")
    for sm in sb['sector_mapping']:
        print(f"    C = {sm['C_band']:4d} ({sm['group']:>15s}) → {sm['sector']}")
        print(f"         Evidence: {sm['evidence']}")
        print()

    print("  SU(5) representation decompositions:")
    for br in sb['branching_rules']:
        print(f"    {br['su5_rep']:>3s}: {br['decomposition']}")
        print(f"         Particles: {br['particles']}")
    print()

    # ===================================================================
    # 3. NEW DIMENSIONLESS RATIO PREDICTIONS
    # ===================================================================
    print("=" * 80)
    print("3. NEW DIMENSIONLESS RATIO PREDICTIONS")
    print("=" * 80)
    print()

    new_preds = new_ratio_predictions(lattice)

    print(f"  {'Name':<30s} {'Value':>12s} {'(C,m)':>10s} {'Lattice':>12s} {'Error':>8s}")
    print("  " + "-" * 80)
    for p in new_preds:
        addr = f"({p['C']:g},{p['m']})"
        pct = p['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else " *  " if abs(pct) < 5 else "    "
        print(f"  {p['name']:<30s} {p['value']:12.4g} {addr:>10s} {p['lattice_val']:12.4g} {pct:+7.2f}%{marker}")
    print()

    n_sub1 = sum(1 for p in new_preds if abs(p['rel_err']) < 0.01)
    n_sub3 = sum(1 for p in new_preds if abs(p['rel_err']) < 0.03)
    print(f"  Sub-1%: {n_sub1}/{len(new_preds)}, Sub-3%: {n_sub3}/{len(new_preds)}")
    print()

    # ===================================================================
    # 4. RG FLOW PATHS ON THE LATTICE
    # ===================================================================
    print("=" * 80)
    print("4. RG FLOW PATHS ON THE LATTICE")
    print("=" * 80)
    print()

    rg = rg_flow_on_lattice(lattice)

    for name in ["α₁⁻¹", "α₂⁻¹", "α₃⁻¹"]:
        print(f"  {name} trajectory:")
        path = rg['paths'][name]
        for p in path:
            addr = f"({p['C']:g},{p['m']})"
            pct = p['rel_err'] * 100
            marker = " *" if abs(pct) < 1 else "  "
            print(f"    Q = {p['Q_GeV']:10.2e} GeV  α⁻¹ = {p['value']:7.2f}  → {addr:>10s}  "
                  f"lat = {p['lattice_val']:7.2f}  ({pct:+5.2f}%){marker}")
        print()

    print("  Transition energies (where nearest lattice site changes):")
    for name in ["α₁⁻¹", "α₂⁻¹", "α₃⁻¹"]:
        trans = rg['transitions'][name]
        if trans:
            print(f"  {name}:")
            for t in trans:
                print(f"    ({t['addr_from'][0]:g},{t['addr_from'][1]}) → "
                      f"({t['addr_to'][0]:g},{t['addr_to'][1]})  "
                      f"at Q ∈ [{t['Q_from']:.2e}, {t['Q_to']:.2e}] GeV")
            print()

    # ===================================================================
    # 5. OUTLIER ANALYSIS: WHAT DOESN'T FIT
    # ===================================================================
    print("=" * 80)
    print("5. OUTLIER ANALYSIS: WHAT DOESN'T FIT")
    print("=" * 80)
    print()

    out = outlier_analysis(lattice)

    print(f"  Tested {len(out['all_results'])} dimensionless ratios:")
    print(f"    < 1%:  {out['n_good']:3d}  (good)")
    print(f"    1-3%:  {out['n_marginal']:3d}  (marginal)")
    print(f"    3-5%:  {out['n_poor']:3d}  (poor)")
    print(f"    > 5%:  {out['n_bad']:3d}  (bad)")
    print()

    print("  Full results (sorted by error):")
    print("  " + "-" * 75)
    print(f"  {'Name':<25s} {'Value':>12s} {'(C,m)':>10s} {'Lattice':>12s} {'Error':>8s} {'Cat':>8s}")
    print("  " + "-" * 75)
    for r in out['all_results']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else " *  " if abs(pct) < 5 else "    "
        print(f"  {r['name']:<25s} {r['value']:12.4g} {addr:>10s} {r['lattice_val']:12.4g} "
              f"{pct:+7.2f}%{marker} {r['category']:>8s}")
    print()

    print("  WORST 5 OUTLIERS:")
    for r in out['worst_5']:
        pct = r['rel_err'] * 100
        print(f"    {r['name']:<25s}  {pct:+6.2f}%  [{r['category']}]")
    print()

    return 0


def cmd_gut_pheno(args: argparse.Namespace) -> int:
    """
    Round 6: phenomenological predictions — proton decay, dark matter,
    muon g-2, information theory, neutrino properties.
    """
    import math as _math
    from physics_test.gut_validate import (
        build_sorted_lattice, nearest_lattice as nl,
        proton_decay_prediction,
        dark_matter_predictions,
        muon_g2_analysis,
        lattice_information_theory,
        neutrino_lattice_predictions,
        build_address_table,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    m_range = list(range(-80, 120))

    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=360.0, include=include)
        for _, v in cs.items():
            all_Cs.append(v)
    Cs = sorted(set(all_Cs))
    lattice = build_sorted_lattice(Cs, m_range)

    # ===================================================================
    # 1. PROTON DECAY PREDICTION
    # ===================================================================
    print("=" * 80)
    print("1. PROTON DECAY LIFETIME PREDICTION")
    print("=" * 80)
    print()

    pd = proton_decay_prediction(lattice)

    print(f"  Lattice-predicted M_GUT = {pd['M_GUT_lattice']:.4e} GeV")
    print(f"  Conventional M_GUT     = {pd['M_GUT_conventional']:.4e} GeV")
    print(f"  Lattice 1/α_GUT        = {pd['alpha_GUT_inv_lattice']:.4f}")
    print(f"  Lattice α_GUT          = {pd['alpha_GUT_lattice']:.6f}")
    print()
    print(f"  Proton lifetime predictions (p → e⁺ + π⁰):")
    print(f"    Simple estimate (lattice):      τ_p ~ {pd['tau_simple_yr']:.2e} years")
    print(f"    Refined estimate (lattice):     τ_p ~ {pd['tau_refined_yr']:.2e} years")
    print(f"    Conventional (M_GUT=1.6e16):    τ_p ~ {pd['tau_conventional_yr']:.2e} years")
    print()
    print(f"  Experimental bound (Super-K):     τ_p > {pd['tau_exp_bound_yr']:.1e} years")
    print(f"  Hyper-K sensitivity:              τ_p ~ 1.0e+35 years")
    print()
    if pd['tau_refined_yr'] > pd['tau_exp_bound_yr']:
        print(f"  ✓ Lattice prediction CONSISTENT with current bound")
    else:
        print(f"  ✗ Lattice prediction EXCLUDED by current bound")
    print(f"  Observable at Hyper-K: {'YES' if pd['observable_at_hyperK'] else 'NO'}")
    print()

    # ===================================================================
    # 2. DARK MATTER MASS PREDICTIONS
    # ===================================================================
    print("=" * 80)
    print("2. DARK MATTER MASS PREDICTIONS")
    print("=" * 80)
    print()

    dm = dark_matter_predictions(lattice)

    print(f"  DM candidates from lattice × m_e: {dm['n_candidates']} in 1 GeV – 100 TeV range")
    print()
    print("  Targeted mass searches:")
    print(f"  {'Label':<25s} {'Target':>10s} {'(C,m)':>10s} {'Lattice':>10s} {'Error':>8s}")
    print("  " + "-" * 68)
    for t in dm['targeted']:
        addr = f"({t['C']:g},{t['m']})"
        pct = t['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {t['label']:<25s} {t['target_GeV']:10.1f} {addr:>10s} "
              f"{t['lattice_GeV']:10.2f} {pct:+7.2f}%{marker}")
    print()

    dfo = dm['dm_from_omega']
    print(f"  From Ω_DM address (360, 15):")
    print(f"    Lattice value:  {dfo['omega_lattice_val']:.4f}")
    print(f"    As coupling:    g_DM = {dfo['as_coupling']:.4f}")
    print(f"    As √Ω × v_H:   {dfo['as_mass_ratio_to_vH']:.1f} GeV")
    print(f"    As √Ω × m_Z:   {dfo['as_mass_ratio_to_mZ']:.1f} GeV")
    print()

    # ===================================================================
    # 3. MUON g-2 AND THE MUON MASS PROBLEM
    # ===================================================================
    print("=" * 80)
    print("3. MUON g-2 AND THE MUON MASS PROBLEM")
    print("=" * 80)
    print()

    g2 = muon_g2_analysis(lattice)

    print(f"  a_μ (experiment):         {g2['a_mu_exp']:.6e}")
    print(f"  a_μ (SM, data-driven):    {g2['a_mu_SM_dd']:.6e}")
    print(f"  Δa_μ (data-driven):       {g2['Delta_a_mu_dd']:.0e}")
    print(f"  Δa_μ (lattice QCD):       {g2['Delta_a_mu_lat']:.0e}")
    print()
    print(f"  Muon mass displacement:   δm_μ/m_μ(lat) = {g2['delta_mass_fraction']:.4f} ({g2['delta_mass_fraction']*100:.2f}%)")
    print(f"  g-2 anomaly fraction:     Δa_μ/a_μ      = {g2['frac_anomaly']:.4e} ({g2['frac_anomaly']*100:.4f}%)")
    print(f"  Ratio (Δa/a) / (δm/m):   {g2['ratio_frac']:.4f}")
    print()
    print(f"  If muon were at lattice mass (194.16 × m_e):")
    print(f"    Hadronic VP shift:      Δa_had ~ {g2['delta_a_had']:.0e}")
    print(f"    This is {g2['delta_a_had']/g2['Delta_a_mu_dd']*100:.0f}% of the data-driven anomaly")
    print()
    print(f"  a_μ/a_e = {g2['a_mu_over_a_e']:.6f} → nearest lattice: {g2['a_mu_ratio_addr']}")
    print()

    # ===================================================================
    # 4. INFORMATION THEORY: LATTICE ENTROPY
    # ===================================================================
    print("=" * 80)
    print("4. INFORMATION THEORY: LATTICE ENTROPY AND SPACING")
    print("=" * 80)
    print()

    addr_table = build_address_table(lattice)
    info = lattice_information_theory(lattice, addr_table)

    print(f"  Occupied sites:    {info['n_occupied']}")
    print(f"  C-bands used:      {info['n_bands']}")
    print(f"  Shannon entropy of C-distribution: H = {info['H_C']:.4f} bits")
    print(f"  Maximum entropy (uniform):         H_max = {info['H_max']:.4f} bits")
    print(f"  Efficiency H/H_max:                {info['H_ratio']:.4f}")
    print()
    print(f"  Global m-spacing statistics:")
    print(f"    Mean spacing:     {info['global_mean_spacing']:.2f}")
    print(f"    Variance:         {info['global_var_spacing']:.2f}")
    print(f"    Dispersion index: {info['dispersion_index']:.2f} (1=Poisson, <1=regular, >1=clustered)")
    print()
    print(f"  Fibonacci spacings: {info['n_fib_spacings']}/{info['n_total_spacings']} "
          f"({info['fib_fraction']*100:.0f}%)")
    print(f"  Combinatorial information: {info['bits_combinatorial']:.1f} bits")
    print()

    print("  Per-band spacing analysis:")
    for c in sorted(info['band_spacings'].keys()):
        bs = info['band_spacings'][c]
        print(f"    C={c:4g}: m-values = {bs['m_values']}")
        print(f"           spacings = {bs['spacings']}")
        print(f"           mean={bs['mean_spacing']:.1f}, min={bs['min_spacing']}, max={bs['max_spacing']}")
        print()

    # ===================================================================
    # 5. NEUTRINO MASS ORDERING AND PREDICTIONS
    # ===================================================================
    print("=" * 80)
    print("5. NEUTRINO MASS ORDERING PREDICTIONS")
    print("=" * 80)
    print()

    nu = neutrino_lattice_predictions(lattice)

    for label in ["NH_minimal", "IH_minimal"]:
        ordering = nu['orderings'][label]
        kind = "Normal Hierarchy" if "NH" in label else "Inverted Hierarchy"
        print(f"  {kind} (minimal, lightest ≈ 0):")
        for i, h in enumerate(ordering['hits']):
            if h['mass_eV'] > 0:
                addr = f"({h['C']:g},{h['m']})"
                pct = h['rel_err'] * 100
                marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
                print(f"    m_{i+1} = {h['mass_eV']:.4e} eV → {addr:>10s}  "
                      f"lattice = {h['lattice_mass_eV']:.4e} eV  ({pct:+.2f}%){marker}")
            else:
                print(f"    m_{i+1} = 0 (lightest)")
        print(f"    Σm_ν = {ordering['sum_masses_eV']*1000:.2f} meV  "
              f"(lattice: {ordering['sum_lattice_eV']*1000:.2f} meV)")
        print()

    print("  Mass ratio m₃/m₂ (NH):")
    mr = nu['mass_ratio_m3_m2']
    print(f"    m₃/m₂ = {mr['ratio']:.4f} → ({mr['C']:g},{mr['m']})  "
          f"lattice = {mr['lattice_val']:.4f}  ({mr['rel_err']*100:+.2f}%)")
    print()

    print("  Neutrinoless double beta decay effective mass <m_ββ>:")
    for mc in nu['mbb_candidates']:
        addr = f"({mc['C']:g},{mc['m']})"
        pct = mc['rel_err'] * 100
        print(f"    {mc['label']:<10s}: <m_ββ> = {mc['mbb_eV']*1000:.2f} meV → {addr:>10s}  "
              f"lattice = {mc['lattice_mbb_eV']*1000:.2f} meV  ({pct:+.2f}%)")
    print()

    return 0


def cmd_gut_struct(args: argparse.Namespace) -> int:
    """
    Round 7: lattice renormalization, cosmological phase transitions,
    anomaly matching, hadron spectrum, lattice topology.
    """
    import math as _math
    from physics_test.gut_validate import (
        build_sorted_lattice, nearest_lattice as nl,
        lattice_renormalization,
        cosmological_phase_transitions,
        anomaly_matching_analysis,
        full_hadron_spectrum,
        lattice_topology,
        build_address_table,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    m_range = list(range(-80, 120))

    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=360.0, include=include)
        for _, v in cs.items():
            all_Cs.append(v)
    Cs = sorted(set(all_Cs))
    lattice = build_sorted_lattice(Cs, m_range)

    # ===================================================================
    # 1. LATTICE RENORMALIZATION
    # ===================================================================
    print("=" * 80)
    print("1. LATTICE RENORMALIZATION: COUPLING RESIDENCE TIMES")
    print("=" * 80)
    print()

    lr = lattice_renormalization(lattice)

    for name in ["α₁⁻¹", "α₂⁻¹", "α₃⁻¹"]:
        sites = lr['residence_times'].get(name, [])
        print(f"  {name}: {len(sites)} lattice crossings")
        if sites:
            for s in sites[:20]:
                bar = "█" * max(1, int(s['residence_dlogQ'] * 5))
                print(f"    ({s['C']:4g},{s['m']:3d}) val={s['value']:7.2f}  "
                      f"Q={s['Q_GeV']:.2e}  Δlog₁₀Q={s['residence_dlogQ']:.2f} {bar}")
            if len(sites) > 20:
                print(f"    ... ({len(sites)-20} more)")
        print()

    print("  SLOW ZONES (longest residence at a single site):")
    for sz in lr['slow_zones']:
        print(f"    {sz['coupling']}: ({sz['C']:g},{sz['m']}) at Q={sz['Q_GeV']:.2e} GeV  "
              f"Δlog₁₀Q = {sz['residence']:.2f}")
    print()

    # ===================================================================
    # 2. COSMOLOGICAL PHASE TRANSITIONS
    # ===================================================================
    print("=" * 80)
    print("2. COSMOLOGICAL PHASE TRANSITIONS ON THE LATTICE")
    print("=" * 80)
    print()

    cpt = cosmological_phase_transitions(lattice)

    print(f"  {'Transition':<30s} {'E [GeV]':>10s} {'(C,m)':>10s} {'Lattice E':>10s} {'Error':>8s}")
    print("  " + "-" * 75)
    for t in cpt['transitions']:
        addr = f"({t['C']:g},{t['m']})"
        pct = t['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else " *  " if abs(pct) < 5 else "    "
        print(f"  {t['name']:<30s} {t['E_GeV']:10.3g} {addr:>10s} {t['lattice_E_GeV']:10.3g} {pct:+7.2f}%{marker}")
    print()

    if cpt['phi_relations']:
        print("  φ-power relations between transition energies:")
        for pr in cpt['phi_relations']:
            print(f"    {pr['high']:<30s} / {pr['low']:<20s} ≈ φ^{pr['phi_power']} "
                  f"({pr['ratio']:.3g} vs {pr['phi_predicted']:.3g}, {pr['frac_err']*100:+.1f}%)")
        print()

    # ===================================================================
    # 3. ANOMALY MATCHING
    # ===================================================================
    print("=" * 80)
    print("3. GAUGE ANOMALY CANCELLATION ON THE LATTICE")
    print("=" * 80)
    print()

    am = anomaly_matching_analysis()

    print("  SM fermion hypercharges:")
    for name, f in am['fermions'].items():
        print(f"    {name:<5s}: Y = {f['Y']:+6.3f}  SU(3)={f['SU3']}  SU(2)={f['SU2']}  [{f['desc']}]")
    print()

    print("  Anomaly cancellation checks:")
    print(f"    Tr[Y]         = {am['tr_Y']:.6f}  {'✓' if abs(am['tr_Y']) < 1e-10 else '✗'}")
    print(f"    Tr[Y³]        = {am['tr_Y3']:.6f}  {'✓' if abs(am['tr_Y3']) < 1e-10 else '✗'}")
    print(f"    Tr[SU(2)²×Y]  = {am['tr_SU2Y']:.6f}  {'✓' if abs(am['tr_SU2Y']) < 1e-10 else '✗'}")
    print(f"    Tr[SU(3)²×Y]  = {am['tr_SU3Y']:.6f}  {'✓' if abs(am['tr_SU3Y']) < 1e-10 else '✗'}")
    print()

    print(f"  Y in units of 1/6: {am['Y_units']}")
    print(f"  Weighted sum (should be 0): {am['sum_weighted']}")
    print(f"  Sum |Y|: {am['sum_abs_Y']:.4f} = 8/3")
    print(f"  {am['denominator_connection']}")
    print()

    # ===================================================================
    # 4. FULL HADRON SPECTRUM
    # ===================================================================
    print("=" * 80)
    print("4. HADRON SPECTRUM ON THE LATTICE")
    print("=" * 80)
    print()

    hs = full_hadron_spectrum(lattice)

    print(f"  {hs['n_hadrons']} hadrons tested (mass/m_π):")
    print(f"    Sub-1%: {hs['n_sub1_pi']}, Sub-3%: {hs['n_sub3_pi']}")
    print()

    print(f"  {'Hadron':<10s} {'Mass':>8s} {'m/m_π':>10s} {'(C,m)':>10s} {'Lattice':>10s} {'Error':>8s} {'Type':>8s}")
    print("  " + "-" * 75)
    for r in hs['pi_ratios']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {r['name']:<10s} {r['mass']:8.4f} {r['ratio']:10.4f} {addr:>10s} "
              f"{r['lattice_val']:10.4f} {pct:+7.2f}%{marker} {r['kind']:>8s}")
    print()

    if hs['meson_cross_ratios']:
        print("  Best inter-meson mass ratios (< 2%):")
        for r in hs['meson_cross_ratios'][:10]:
            addr = f"({r['C']:g},{r['m']})"
            pct = r['rel_err'] * 100
            print(f"    {r['pair']:<12s} = {r['ratio']:8.4f} → {addr:>10s}  ({pct:+.2f}%)")
        print()

    # Baryon ratios to proton
    if hs['p_ratios']:
        print("  Baryon/proton mass ratios:")
        for r in hs['p_ratios']:
            addr = f"({r['C']:g},{r['m']})"
            pct = r['rel_err'] * 100
            marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
            print(f"    {r['name']:<5s}: m/m_p = {r['ratio']:8.5f} → {addr:>10s}  ({pct:+.2f}%){marker}")
        print()

    # ===================================================================
    # 5. LATTICE TOPOLOGY
    # ===================================================================
    print("=" * 80)
    print("5. LATTICE TOPOLOGY")
    print("=" * 80)
    print()

    addr_table = build_address_table(lattice)
    topo = lattice_topology(addr_table)

    print(f"  Sites:          {topo['n_sites']}")
    print(f"  Same-C edges:   {topo['n_same_c_edges']} (Δm = 1 within a band)")
    print(f"  Same-m edges:   {topo['n_same_m_edges']} (cross-band at same m)")
    print(f"  Total edges:    {topo['n_total_edges']}")
    print(f"  Components:     {topo['n_components']}")
    print(f"  Largest:        {topo['largest_component']} sites")
    print()

    if topo['horizontal_lines']:
        print("  Horizontal lines (same m, multiple C-bands):")
        for m_val in sorted(topo['horizontal_lines'].keys()):
            cs = topo['horizontal_lines'][m_val]
            print(f"    m = {m_val:4d}: C = {cs}  ({len(cs)} bands)")
        print()

    if topo['vertical_chains']:
        print("  Vertical chains (consecutive m in same C):")
        for vc in topo['vertical_chains']:
            print(f"    C = {vc['C']:4g}: m = {vc['chain']}  (length {vc['length']})")
        print()

    return 0


def cmd_gut_precision(args: argparse.Namespace) -> int:
    """
    Round 8: electroweak precision, lattice action, CP violation, inflation.
    """
    import math as _math
    from physics_test.gut_validate import (
        build_sorted_lattice, nearest_lattice as nl,
        electroweak_precision,
        lattice_action_principle,
        cp_violation_analysis,
        inflation_analysis,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    m_range = list(range(-80, 120))

    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=360.0, include=include)
        for _, v in cs.items():
            all_Cs.append(v)
    Cs = sorted(set(all_Cs))
    lattice = build_sorted_lattice(Cs, m_range)

    # ===================================================================
    # 1. ELECTROWEAK PRECISION OBSERVABLES
    # ===================================================================
    print("=" * 80)
    print("1. ELECTROWEAK PRECISION OBSERVABLES")
    print("=" * 80)
    print()

    ew = electroweak_precision(lattice)

    print("  Mass ratios on the lattice:")
    print(f"  {'Quantity':<35s} {'Value':>10s} {'(C,m)':>10s} {'Lattice':>10s} {'Error':>8s}")
    print("  " + "-" * 78)
    for r in ew['mass_ratios']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {r['name']:<35s} {r['value']:10.6f} {addr:>10s} "
              f"{r['lattice']:10.6f} {pct:+7.2f}%{marker}")
    print()

    print(f"  ρ parameter: {ew['rho']:.6f} (should be 1.0000)")
    print(f"  1/(G_F·v²) = {ew['fermi_check']:.4f} (should be √2 = {_math.sqrt(2):.4f})")
    print()

    print("  Weinberg angle tests:")
    for a in ew['angle_tests']:
        addr = f"({a['C']:g},{a['m']})"
        pct = a['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"    {a['name']:<25s} = {a['value']:.6f} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    # ===================================================================
    # 2. LATTICE ACTION PRINCIPLE
    # ===================================================================
    print("=" * 80)
    print("2. LATTICE ACTION PRINCIPLE")
    print("=" * 80)
    print()

    la = lattice_action_principle(lattice)

    print(f"  Occupied sites: {la['n_occupied']}")
    print(f"  Total |m|:      {la['total_abs_m']}  (mean = {la['mean_abs_m']:.1f})")
    print(f"  C-weighted |m|: {la['c_weighted_m']}")
    print(f"  Pairwise dist:  {la['total_pairwise_dist']}")
    print()
    print(f"  C-band entropy:  {la['c_entropy']:.3f} bits (max = {la['max_entropy']:.3f})")
    print(f"  Efficiency:      {la['entropy_efficiency']*100:.1f}%")
    print(f"  Golden zone fraction (m ∈ [-5,15]): {la['frac_golden_zone']*100:.1f}%")
    print(f"  Fibonacci spacing fraction: {la['n_fib_gaps']}/{la['n_gaps']} = {la['fib_spacing_frac']*100:.1f}%")
    print()

    # ===================================================================
    # 3. CP VIOLATION
    # ===================================================================
    print("=" * 80)
    print("3. CP VIOLATION ON THE LATTICE")
    print("=" * 80)
    print()

    cp = cp_violation_analysis(lattice)

    print("  Jarlskog invariants:")
    for r in cp['jarlskog']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"    {r['name']:<20s} = {r['value']:12.4e} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    pf = cp['phase_fraction']
    pct = pf['rel_err'] * 100
    print(f"  δ_CKM = {cp['delta_CKM_rad']:.3f} rad (68.4°)")
    print(f"  δ_CKM/(2π) = {pf['value']:.4f} → ({pf['C']:g},{pf['m']})  ({pct:+.2f}%)")
    print()

    print("  Wolfenstein parameters:")
    for w in cp['wolfenstein']:
        addr = f"({w['C']:g},{w['m']})"
        pct = w['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"    {w['name']:<20s} = {w['value']:.6f} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    # ===================================================================
    # 4. INFLATION PARAMETERS
    # ===================================================================
    print("=" * 80)
    print("4. INFLATION PARAMETERS ON THE LATTICE")
    print("=" * 80)
    print()

    inf = inflation_analysis(lattice)

    print("  CMB observables:")
    for r in inf['observables']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"    {r['name']:<25s} = {r['value']:12.4e} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    print("  Slow-roll estimates:")
    for r in inf['slow_roll']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"    {r['name']:<25s} = {r['value']:12.4e} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    print(f"  N from n_s: {inf['N_from_ns']:.1f} e-folds")
    if inf['inflation_scale_GeV'] > 0:
        print(f"  Inflation energy scale: V^{{1/4}} ≈ {inf['inflation_scale_GeV']:.2e} GeV")
    print()

    return 0


def cmd_gut_nuclear(args: argparse.Namespace) -> int:
    """
    Round 9: nuclear physics, astrophysical scales, lattice bootstrap,
    coupling ratios, QCD observables.
    """
    import math as _math
    from physics_test.gut_validate import (
        build_sorted_lattice,
        nuclear_physics_analysis,
        astrophysical_scales,
        lattice_bootstrap,
        coupling_ratio_evolution,
        qcd_observables,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    m_range = list(range(-80, 120))

    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=360.0, include=include)
        for _, v in cs.items():
            all_Cs.append(v)
    Cs = sorted(set(all_Cs))
    lattice = build_sorted_lattice(Cs, m_range)

    # ===================================================================
    # 1. NUCLEAR PHYSICS
    # ===================================================================
    print("=" * 80)
    print("1. NUCLEAR PHYSICS ON THE LATTICE")
    print("=" * 80)
    print()

    nuc = nuclear_physics_analysis(lattice)

    print(f"  {'Ratio':<35s} {'Value':>12s} {'(C,m)':>10s} {'Lattice':>12s} {'Error':>8s}")
    print("  " + "-" * 82)
    for r in nuc['nuclear_ratios']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {r['name']:<35s} {r['value']:12.4f} {addr:>10s} "
              f"{r['lattice']:12.4f} {pct:+7.2f}%{marker}")
    print()

    print("  Nuclear magic numbers on the lattice:")
    for r in nuc['magic_numbers']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"    {r['magic_number']:>4d} → {addr:>10s}  lattice = {r['lattice']:8.2f}  ({pct:+.2f}%){marker}")
    print()

    if nuc['magic_ratios']:
        print("  Magic number ratios (< 5%):")
        for r in nuc['magic_ratios'][:8]:
            addr = f"({r['C']:g},{r['m']})"
            pct = r['rel_err'] * 100
            print(f"    {r['pair']:<8s} = {r['ratio']:8.3f} → {addr:>10s}  ({pct:+.2f}%)")
        print()

    print("  Bethe-Weizsäcker coefficients:")
    for r in nuc['bethe_weiszacker']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"    {r['name']:<20s} = {r['value']:10.5f} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    # ===================================================================
    # 2. ASTROPHYSICAL SCALES
    # ===================================================================
    print("=" * 80)
    print("2. ASTROPHYSICAL MASS SCALES")
    print("=" * 80)
    print()

    astro = astrophysical_scales(lattice)

    for r in astro['astro_ratios']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {r['name']:<30s} = {r['value']:12.4e} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    print(f"  M_Pl/m_p as φ-power: φ^{astro['MPl_over_mp_phi_power']:.2f}")
    print(f"  Nearest integer: φ^{astro['nearest_phi_power']}  (error: {astro['phi_power_err']*100:+.2f}%)")
    print()

    print("  Dirac large numbers:")
    for r in astro['large_numbers']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"    {r['name']:<30s} = {r['value']:12.4e} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    # ===================================================================
    # 3. LATTICE BOOTSTRAP
    # ===================================================================
    print("=" * 80)
    print("3. LATTICE BOOTSTRAP: DERIVING PHYSICS FROM SEEDS")
    print("=" * 80)
    print()

    boot = lattice_bootstrap(lattice)

    print(f"  Seeds: {boot['seeds']}")
    print()
    print("  Minimum hops from seed set to known sites:")
    for sd in boot['site_distances']:
        print(f"    {sd['name']:<25s} ({sd['C']:g},{sd['m']}) → {sd['hops']} hops")
    print()

    # ===================================================================
    # 4. COUPLING RATIO EVOLUTION
    # ===================================================================
    print("=" * 80)
    print("4. COUPLING RATIO EVOLUTION")
    print("=" * 80)
    print()

    cre = coupling_ratio_evolution(lattice)

    print(f"  {'E [GeV]':>10s}  {'α₁⁻¹/α₂⁻¹':>12s} {'(C,m)':>8s} {'err':>6s}  "
          f"{'α₂⁻¹/α₃⁻¹':>12s} {'(C,m)':>8s} {'err':>6s}")
    print("  " + "-" * 78)
    for row in cre['evolution']:
        r12 = row.get("α₁⁻¹/α₂⁻¹", {})
        r23 = row.get("α₂⁻¹/α₃⁻¹", {})
        if r12 and r23:
            a12 = f"({r12['C']:g},{r12['m']})"
            a23 = f"({r23['C']:g},{r23['m']})"
            print(f"  {row['E_GeV']:10.1e}  {r12['value']:12.4f} {a12:>8s} {r12['rel_err']*100:+5.1f}%  "
                  f"{r23['value']:12.4f} {a23:>8s} {r23['rel_err']*100:+5.1f}%")
    print()

    if cre['best_hits']:
        print("  Best coupling-ratio lattice hits (< 0.5%):")
        for h in cre['best_hits'][:10]:
            addr = f"({h['C']:g},{h['m']})"
            print(f"    {h['ratio']:<15s} at E={h['E_GeV']:.1e}: {h['value']:.4f} → {addr:>10s}  ({h['rel_err']*100:+.2f}%)")
        print()

    # ===================================================================
    # 5. QCD OBSERVABLES
    # ===================================================================
    print("=" * 80)
    print("5. QCD NON-PERTURBATIVE OBSERVABLES")
    print("=" * 80)
    print()

    qcd = qcd_observables(lattice)

    for r in qcd['qcd_ratios']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {r['name']:<30s} = {r['value']:10.4f} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    reg = qcd['regge']
    pct = reg['rel_err'] * 100
    print(f"  Regge slope: α'·m_p² = {reg['alpha_prime_mp2']:.4f} → ({reg['C']:g},{reg['m']})  ({pct:+.2f}%)")
    print()

    return 0


def cmd_gut_final(args: argparse.Namespace) -> int:
    """
    Round 10: Yukawa hierarchy, vacuum energy, grid symmetries, scorecard.
    """
    import math as _math
    from physics_test.gut_validate import (
        build_sorted_lattice,
        yukawa_hierarchy,
        vacuum_energy_analysis,
        grid_symmetries,
        comprehensive_scorecard,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    m_range = list(range(-80, 120))

    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=360.0, include=include)
        for _, v in cs.items():
            all_Cs.append(v)
    Cs = sorted(set(all_Cs))
    lattice = build_sorted_lattice(Cs, m_range)

    # ===================================================================
    # 1. YUKAWA COUPLING HIERARCHY
    # ===================================================================
    print("=" * 80)
    print("1. YUKAWA COUPLING HIERARCHY")
    print("=" * 80)
    print()

    yh = yukawa_hierarchy(lattice)

    print("  Inter-generation mass ratios:")
    print(f"  {'Sector':<20s} {'Ratio':<15s} {'Value':>10s} {'(C,m)':>10s} {'Error':>8s}")
    print("  " + "-" * 68)
    for r in yh['inter_generation']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        phi_p = f"  φ^{r.get('phi_power', 0):.1f}" if 'phi_power' in r else ""
        print(f"  {r['sector']:<20s} {r['ratio_name']:<15s} {r['ratio']:10.2f} {addr:>10s} "
              f"{pct:+7.2f}%{marker}{phi_p}")
    print()

    print("  Cross-sector same-generation ratios:")
    for r in yh['cross_sector']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"    Gen {r['generation']}: {r['ratio_name']:<12s} = {r['ratio']:10.2f} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    print("  Yukawa couplings as lattice quantities:")
    for r in yh['yukawa_lattice']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"    y_{r['fermion']:<3s} = {r['yukawa']:12.5e} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    # ===================================================================
    # 2. VACUUM ENERGY / COSMOLOGICAL CONSTANT
    # ===================================================================
    print("=" * 80)
    print("2. VACUUM ENERGY AND COSMOLOGICAL CONSTANT")
    print("=" * 80)
    print()

    ve = vacuum_energy_analysis(lattice)

    print(f"  ρ_vac/ρ_Pl = {ve['rho_vac_over_Pl']:.2e}")
    print(f"  log_φ(ρ_Pl/ρ_vac) = {ve['log_phi_ratio']:.2f}")
    print(f"  Nearest integer: φ^{ve['nearest_phi_power']}  (err: {ve['phi_power_err']*100:+.2f}%)")
    print()

    for r in ve['quantities']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {r['name']:<35s} = {r['value']:12.4e} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    # ===================================================================
    # 3. HIDDEN SYMMETRIES
    # ===================================================================
    print("=" * 80)
    print("3. HIDDEN SYMMETRIES OF THE (C,m) GRID")
    print("=" * 80)
    print()

    gs = grid_symmetries(lattice)

    print(f"  Sites: {gs['n_sites']}")
    print(f"  m-center: {gs['m_center']:.1f}")
    print(f"  m-reflection fraction: {gs['m_reflection_frac']*100:.1f}%")
    print(f"  C-duality (360/C at same m): {gs['c_dual_frac']*100:.1f}%")
    print(f"  Diagonal (C₂/C₁ = φ^(m₁-m₂)) hits: {gs['diagonal_hits']}")
    print()

    print("  Best m-translations (Δm → number of shared sites):")
    for dm, count in gs['best_translations']:
        print(f"    Δm = {dm:+3d}: {count} sites overlap")
    print()

    print("  Most frequent C-band pairings (same m):")
    for (c1, c2), count in gs['best_c_pairs']:
        print(f"    ({c1:g}, {c2:g}): {count} shared m-values")
    print()

    # ===================================================================
    # 4. COMPREHENSIVE SCORECARD
    # ===================================================================
    print("=" * 80)
    print("4. COMPREHENSIVE SCORECARD: ALL TESTED QUANTITIES")
    print("=" * 80)
    print()

    sc = comprehensive_scorecard(lattice)

    print(f"  Total quantities tested: {sc['n_total']}")
    print(f"  Sub-1%:  {sc['n_sub1']:3d}  ({sc['n_sub1']/sc['n_total']*100:.0f}%)")
    print(f"  Sub-3%:  {sc['n_sub3']:3d}  ({sc['n_sub3']/sc['n_total']*100:.0f}%)")
    print(f"  Sub-5%:  {sc['n_sub5']:3d}  ({sc['n_sub5']/sc['n_total']*100:.0f}%)")
    print(f"  Over 5%: {sc['n_over5']:3d}  ({sc['n_over5']/sc['n_total']*100:.0f}%)")
    print()

    print(f"  C-band distribution: {sc['c_distribution']}")
    print()

    print(f"  {'Rank':>4s}  {'Quantity':<30s} {'Value':>12s} {'(C,m)':>10s} {'Lattice':>12s} {'Error':>8s}")
    print("  " + "-" * 82)
    for i, r in enumerate(sc['all_results']):
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        tier = "***" if abs(pct) < 1 else "** " if abs(pct) < 3 else "*  " if abs(pct) < 5 else "   "
        print(f"  {i+1:4d}  {r['name']:<30s} {r['value']:12.4e} {addr:>10s} "
              f"{r['lattice']:12.4e} {pct:+7.2f}% {tier}")
    print()

    return 0


def cmd_gut_deep2(args: argparse.Namespace) -> int:
    """
    Round 11: black hole thermodynamics, lepton universality,
    base-360 derivation, fine-tuning integer decomposition.
    """
    import math as _math
    from physics_test.gut_validate import (
        build_sorted_lattice,
        black_hole_thermodynamics,
        lepton_universality,
        base_360_derivation,
        fine_tuning_integers,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    m_range = list(range(-80, 120))

    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=360.0, include=include)
        for _, v in cs.items():
            all_Cs.append(v)
    Cs = sorted(set(all_Cs))
    lattice = build_sorted_lattice(Cs, m_range)

    # ===================================================================
    # 1. BLACK HOLE THERMODYNAMICS
    # ===================================================================
    print("=" * 80)
    print("1. BLACK HOLE THERMODYNAMICS ON THE LATTICE")
    print("=" * 80)
    print()

    bh = black_hole_thermodynamics(lattice)

    for r in bh['bh_ratios']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {r['name']:<35s} = {r['value']:10.4f} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    se = bh['planck_bh_entropy']
    print(f"  Planck-BH entropy 4π = {se['value']:.4f} → ({se['C']:g},{se['m']})  ({se['rel_err']*100:+.2f}%)")
    ht = bh['hawking_T_ratio']
    print(f"  Hawking T/M_Pl = 1/(8π) = {ht['value']:.5f} → ({ht['C']:g},{ht['m']})  ({ht['rel_err']*100:+.2f}%)")
    print(f"  log_φ(M_Chandrasekhar) ≈ {bh['log_phi_Chandrasekhar']:.1f}")
    print()

    # ===================================================================
    # 2. LEPTON UNIVERSALITY
    # ===================================================================
    print("=" * 80)
    print("2. LEPTON UNIVERSALITY AND B-PHYSICS")
    print("=" * 80)
    print()

    lu = lepton_universality(lattice)

    print("  B/D-meson mass ratios:")
    for r in lu['meson_ratios']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"    {r['name']:<15s} = {r['value']:8.4f} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    print("  R-ratios (lepton universality):")
    for r in lu['r_ratios']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"    {r['name']:<25s} = {r['value']:8.5f} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    # ===================================================================
    # 3. WHY BASE 360?
    # ===================================================================
    print("=" * 80)
    print("3. WHY BASE 360? DERIVATION ATTEMPTS")
    print("=" * 80)
    print()

    b360 = base_360_derivation()

    print("  Candidate constructions yielding 360:")
    for name, val in b360['candidates'].items():
        hit = " ✓" if val == 360 else ""
        print(f"    {name:<40s} = {val:>6d}{hit}")
    print()

    print(f"  τ(360) = {b360['n_divisors']} = dim(SU(5))  ← number of divisors")
    print(f"  φ(360) = {b360['euler_totient']} = 4 × dim(SU(5))  ← Euler totient")
    print()

    print("  360 = Σdim(SM) × (dim(SU(5)) + h(SU(5)) + 1) = 12 × 30")
    print("      = (dim(SU(3)) + dim(SU(2)) + dim(U(1))) × (24 + 5 + 1)")
    print()

    # ===================================================================
    # 4. THE FINE-TUNING INTEGERS: 91 AND 588
    # ===================================================================
    print("=" * 80)
    print("4. NUMBER THEORY OF 91 (HIERARCHY) AND 588 (CC)")
    print("=" * 80)
    print()

    ft = fine_tuning_integers()

    for label in ["91 (hierarchy)", "588 (CC)"]:
        info = ft[label]
        print(f"  {label}:")
        print(f"    Factorization: {info['factorization']}")
        print(f"    Zeckendorf (Fibonacci sum): {info['n']} = {' + '.join(str(f) for f in info['zeckendorf'])}")
        print(f"    Is Fibonacci: {info['is_fibonacci']}")
        for c in info['connections']:
            print(f"    • {c}")
        print()

    rel = ft['relationship']
    print(f"  Relationship between 91 and 588:")
    print(f"    588 / 91 = {rel['ratio']:.4f}")
    print(f"    588 mod 91 = {rel['remainder']}")
    print(f"    gcd(91, 588) = {rel['gcd']}")
    print(f"    {rel['note_588_eq_91x6_plus_42']}")
    print(f"    {rel['note_gcd']}")
    print()

    return 0


def cmd_gut_predict2(args: argparse.Namespace) -> int:
    """
    Round 12: mathematical constants, condensed matter, sharp predictions.
    """
    import math as _math
    from physics_test.gut_validate import (
        build_sorted_lattice,
        mathematical_constants,
        condensed_matter_constants,
        sharp_predictions,
    )

    PHI = (1.0 + _math.sqrt(5.0)) / 2.0
    include = tuple(s.strip() for s in args.include.split(",") if s.strip())
    m_range = list(range(-80, 120))

    all_Cs: list[float] = []
    for g in standard_model_gauge_groups():
        cs = candidate_Cs_from_group(g, base=360.0, include=include)
        for _, v in cs.items():
            all_Cs.append(v)
    Cs = sorted(set(all_Cs))
    lattice = build_sorted_lattice(Cs, m_range)

    # ===================================================================
    # 1. MATHEMATICAL CONSTANTS
    # ===================================================================
    print("=" * 80)
    print("1. MATHEMATICAL CONSTANTS ON THE LATTICE")
    print("=" * 80)
    print()

    mc = mathematical_constants(lattice)

    for r in mc['constants']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {r['name']:<35s} = {r['value']:12.6f} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    # ===================================================================
    # 2. CONDENSED MATTER CONSTANTS
    # ===================================================================
    print("=" * 80)
    print("2. CONDENSED MATTER AND METROLOGY")
    print("=" * 80)
    print()

    cm = condensed_matter_constants(lattice)

    for r in cm['condensed_matter']:
        addr = f"({r['C']:g},{r['m']})"
        pct = r['rel_err'] * 100
        marker = " ***" if abs(pct) < 1 else " ** " if abs(pct) < 3 else "    "
        print(f"  {r['name']:<35s} = {r['value']:12.6e} → {addr:>10s}  ({pct:+.2f}%){marker}")
    print()

    # ===================================================================
    # 3. SHARP PREDICTIONS
    # ===================================================================
    print("=" * 80)
    print("3. SHARP TESTABLE PREDICTIONS FROM THE LATTICE")
    print("=" * 80)
    print()

    sp = sharp_predictions(lattice)

    for p in sp['predictions']:
        print(f"  {p['name']}:")
        if isinstance(p['predicted'], float) and p['predicted'] > 1e10:
            print(f"    Lattice predicts: {p['predicted']:.2e}")
        else:
            print(f"    Lattice predicts: {p['predicted']:.4f}")
        print(f"    Current exp:     {p['current_exp']}")
        if p['discrepancy_pct'] != 0:
            print(f"    Discrepancy:     {p['discrepancy_pct']:+.3f}%")
        print(f"    How to test:     {p['testable']}")
        print()

    return 0


def cmd_solve_K(args: argparse.Namespace) -> int:
    K = temperature_K_from_frequency(args.m, args.F0)
    print(f"m = {args.m}")
    print(f"F0 = {args.F0:.15g} Hz")
    print(f"K = {K:.15g} K")
    return 0


def _best_hit_for_target(
    *,
    target_name: str,
    target_value: float,
    Cs: list[float],
    m_values: list[float],
    max_rel_err: float,
):
    hits = scan_candidates(Cs=Cs, m_values=m_values, target_G=target_value)
    best = hits[0]
    ok = abs(best.rel_err) <= max_rel_err
    return ok, best


def cmd_pair_forces(args: argparse.Namespace) -> int:
    """
    "Pair" test:
      1) Find best (C,m) for each force's dimensionless coupling target via G=C/phi^m
      2) Use the same m in F0=phi^m*kB*K/h with K chosen per force:
         - gravity: K = CMB temperature
         - EM/strong/weak: K = E/kB for chosen characteristic energy scale E
    """
    targets = {t.name: t.value for t in known_targets()}

    # Candidate (C,m) grid
    Cs = _parse_Cs(args.Cs)
    if not Cs:
        s = get_candidate_set(args.set_name)
        Cs = list(s.values)
        print(f"Using candidate set: {s.name} ({len(Cs)} values) - {s.note}")
    else:
        print(f"Using custom Cs list: n={len(Cs)}")

    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))  # force integer stepping

    # Force definitions: (label, target_key, K_spec)
    # Note: EM uses inverse by default since your worked example is near 1/alpha.
    em_target = args.em_target
    strong_target = args.strong_target
    weak_target = args.weak_target
    gravity_target = args.gravity_target

    missing = [k for k in (em_target, strong_target, weak_target, gravity_target) if k not in targets]
    if missing:
        raise SystemExit(f"Unknown target(s): {missing}. Run `python -m physics_test.cli list-targets`.")

    # Temperatures from energy scales (E -> K)
    E_em_J = energy_J_from_eV(args.em_eV)
    E_strong_J = energy_J_from_MeV(args.strong_MeV)
    E_weak_J = energy_J_from_GeV(args.weak_GeV)

    K_em = temperature_K_from_energy_J(E_em_J)
    K_strong = temperature_K_from_energy_J(E_strong_J)
    K_weak = temperature_K_from_energy_J(E_weak_J)
    K_grav = args.gravity_K

    print(f"\nPairing settings (you can override these):")
    print(f"- tol(|rel_err|) = {args.max_rel_err}")
    print(f"- EM:    target={em_target},  E={args.em_eV} eV  -> K={K_em:.6g} K")
    print(f"- Strong: target={strong_target}, E={args.strong_MeV} MeV -> K={K_strong:.6g} K")
    print(f"- Weak:  target={weak_target}, E={args.weak_GeV} GeV -> K={K_weak:.6g} K")
    print(f"- Grav:  target={gravity_target}, K={K_grav} K (e.g., CMB)\n")

    rows = [
        ("EM", em_target, K_em, E_em_J),
        ("Strong", strong_target, K_strong, E_strong_J),
        ("Weak", weak_target, K_weak, E_weak_J),
        ("Gravity", gravity_target, K_grav, None),
    ]

    print("force   status  target_key           best(C,m)         G_pred              rel_err      K_used(K)        F0_pred(Hz)")
    for label, tkey, K_used, E_J in rows:
        ok, best = _best_hit_for_target(
            target_name=tkey,
            target_value=targets[tkey],
            Cs=Cs,
            m_values=m_values,
            max_rel_err=args.max_rel_err,
        )
        status = "PASS" if ok else "FAIL"
        F0 = frequency_F0(best.m, K_used)

        print(
            f"{label:7s} {status:5s}  {tkey:18s}  "
            f"C={best.C:8.6g}, m={int(best.m):4d}  "
            f"{best.G:18.12g}  {best.rel_err: .3e}  "
            f"{K_used:12.6g}  {F0: .6e}"
        )

        # For quantum forces, also show the baseline frequency f=E/h to make the scaling visible.
        if E_J is not None:
            f0_base = frequency_Hz_from_energy_J(E_J)
            # since K=E/kB, model gives F0 = phi^m * (E/h)
            print(f"         (baseline f=E/h={f0_base:.6e} Hz; model scales by phi^m)")
    return 0


def cmd_pair_forces_option2(args: argparse.Namespace) -> int:
    """
    Option 2 for quantum forces:
      - Choose observed/proxy frequencies F0 for EM/strong/weak
      - Fit (C,m) from coupling targets via G=C/phi^m (within tolerance)
      - Solve implied temperatures: K = F0*h/(kB*phi^m)
      - Gravity can still use a fixed K (e.g., CMB) or its own F0 if desired.
    """
    targets = {t.name: t.value for t in known_targets()}

    # Candidate (C,m) grid
    Cs = _parse_Cs(args.Cs)
    if not Cs:
        s = get_candidate_set(args.set_name)
        Cs = list(s.values)
        print(f"Using candidate set: {s.name} ({len(Cs)} values) - {s.note}")
    else:
        print(f"Using custom Cs list: n={len(Cs)}")

    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))  # integer stepping

    # Targets
    for k in (args.em_target, args.strong_target, args.weak_target, args.gravity_target):
        if k not in targets:
            raise SystemExit(f"Unknown target {k!r}. Run `python -m physics_test.cli list-targets`.")

    # Frequencies
    def _F0_from_arg(value_hz: float | None, preset_key: str | None, label: str) -> float:
        if preset_key:
            return float(get_preset(preset_key).F0_hz)
        if value_hz is None:
            raise SystemExit(f"Must provide either --{label}-F0 or --{label}-preset")
        return float(value_hz)

    F0_em = _F0_from_arg(args.em_F0, args.em_preset, "em")
    F0_s = _F0_from_arg(args.strong_F0, args.strong_preset, "strong")
    F0_w = _F0_from_arg(args.weak_F0, args.weak_preset, "weak")

    # Gravity: either fixed K (default CMB) OR a gravity F0 preset/value.
    F0_g: float | None
    if args.gravity_F0 is not None or args.gravity_preset is not None:
        F0_g = _F0_from_arg(args.gravity_F0, args.gravity_preset, "gravity")
    else:
        F0_g = None

    # Precompute hits within tolerance for each force target
    def _hits_for(target_key: str) -> list[object]:
        hits = scan_candidates(Cs=Cs, m_values=m_values, target_G=targets[target_key])
        hits = filter_hits_by_rel_err(hits, max_abs_rel_err=args.max_rel_err)
        return hits[: max(1, args.max_hits)]

    em_hits = _hits_for(args.em_target)
    s_hits = _hits_for(args.strong_target)
    w_hits = _hits_for(args.weak_target)
    g_hits = _hits_for(args.gravity_target)

    print("\nOption-2 pairing (solve K from F0 using fitted m):")
    print(f"- tol(|rel_err|) = {args.max_rel_err}")
    print(f"- EM:    target={args.em_target}, F0={F0_em:.6e} Hz")
    print(f"- Strong: target={args.strong_target}, F0={F0_s:.6e} Hz")
    print(f"- Weak:  target={args.weak_target}, F0={F0_w:.6e} Hz")
    if F0_g is None:
        print(f"- Gravity: target={args.gravity_target}, K={args.gravity_K} K (fixed)")
    else:
        print(f"- Gravity: target={args.gravity_target}, F0={F0_g:.6e} Hz (solve K)")
    print("")

    shown = 0
    for eh in em_hits:
        m_em = int(eh.m)
        if args.em_m_sign == "positive" and m_em <= 0:
            continue
        if args.em_m_sign == "negative" and m_em >= 0:
            continue

        K_em = temperature_K_from_frequency(m_em, F0_em)
        for sh in s_hits:
            m_s = int(sh.m)
            if args.strong_m_sign == "positive" and m_s <= 0:
                continue
            if args.strong_m_sign == "negative" and m_s >= 0:
                continue
            K_s = temperature_K_from_frequency(m_s, F0_s)
            for wh in w_hits:
                m_w = int(wh.m)
                if args.weak_m_sign == "positive" and m_w <= 0:
                    continue
                if args.weak_m_sign == "negative" and m_w >= 0:
                    continue
                K_w = temperature_K_from_frequency(m_w, F0_w)
                for gh in g_hits:
                    m_g = int(gh.m)
                    if args.gravity_m_sign == "positive" and m_g <= 0:
                        continue
                    if args.gravity_m_sign == "negative" and m_g >= 0:
                        continue

                    if F0_g is None:
                        K_g = args.gravity_K
                        F0g_pred = frequency_F0(m_g, K_g)
                        gravity_line = f"K={K_g:.6g}K -> F0_pred={F0g_pred:.3e}Hz"
                    else:
                        K_g = temperature_K_from_frequency(m_g, F0_g)
                        gravity_line = f"F0={F0_g:.3e}Hz -> K={K_g:.6g}K"

                    # Optional sanity bounds on K (global + per-force)
                    if args.K_min is not None and min(K_em, K_s, K_w, K_g) < args.K_min:
                        continue
                    if args.K_max is not None and max(K_em, K_s, K_w, K_g) > args.K_max:
                        continue

                    def _k_ok(K: float, kmin: float | None, kmax: float | None) -> bool:
                        if kmin is not None and K < kmin:
                            return False
                        if kmax is not None and K > kmax:
                            return False
                        return True

                    if not _k_ok(K_em, args.em_K_min, args.em_K_max):
                        continue
                    if not _k_ok(K_s, args.strong_K_min, args.strong_K_max):
                        continue
                    if not _k_ok(K_w, args.weak_K_min, args.weak_K_max):
                        continue
                    if not _k_ok(K_g, args.gravity_K_min, args.gravity_K_max):
                        continue

                    shown += 1
                    print(
                        f"#{shown}: m=[EM {m_em}, S {m_s}, W {m_w}, G {m_g}]  "
                        f"rel_err=[{eh.rel_err:.2e}, {sh.rel_err:.2e}, {wh.rel_err:.2e}, {gh.rel_err:.2e}]"
                    )
                    print(f"  EM:    C={eh.C:g}  K={K_em:.6g} K")
                    print(f"  Strong: C={sh.C:g}  K={K_s:.6g} K")
                    print(f"  Weak:  C={wh.C:g}  K={K_w:.6g} K")
                    print(f"  Grav:  C={gh.C:g}  {gravity_line}")
                    print("")
                    if shown >= args.max_results:
                        return 0

    print(f"Done. Found {shown} configurations.")
    return 0


def cmd_pair_forces_all(args: argparse.Namespace) -> int:
    """
    Explore combinations of target choices + energy scales and report configurations
    where all 4 forces achieve |rel_err| <= tolerance on the G fit.
    """
    targets = {t.name: t.value for t in known_targets()}

    # Apply gravity band presets if requested and explicit bounds not provided.
    # These are rough observational bands.
    gravity_bands: dict[str, tuple[float | None, float | None]] = {
        "any": (None, None),
        "ligo": (10.0, 1000.0),
        "lisa": (1e-4, 1e-1),
        "pta": (1e-9, 1e-7),
        "cmb": (1e-18, 1e-16),
    }
    if args.gravity_band not in gravity_bands:
        raise SystemExit(f"Unknown gravity band {args.gravity_band!r}")
    if args.gravity_f0_min is None and args.gravity_f0_max is None and args.gravity_band != "any":
        args.gravity_f0_min, args.gravity_f0_max = gravity_bands[args.gravity_band]

    # Candidate (C,m) grid
    Cs = _parse_Cs(args.Cs)
    if not Cs:
        s = get_candidate_set(args.set_name)
        Cs = list(s.values)
        print(f"Using candidate set: {s.name} ({len(Cs)} values) - {s.note}")
    else:
        print(f"Using custom Cs list: n={len(Cs)}")

    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))  # integer stepping

    # Candidate target keys for each force
    em_keys = [k.strip() for k in args.em_targets.split(",") if k.strip()]
    strong_keys = [k.strip() for k in args.strong_targets.split(",") if k.strip()]
    weak_keys = [k.strip() for k in args.weak_targets.split(",") if k.strip()]
    grav_keys = [k.strip() for k in args.gravity_targets.split(",") if k.strip()]

    for group, keys in [("EM", em_keys), ("Strong", strong_keys), ("Weak", weak_keys), ("Gravity", grav_keys)]:
        missing = [k for k in keys if k not in targets]
        if missing:
            raise SystemExit(
                f"Unknown {group} target(s): {missing}. Run `python -m physics_test.cli list-targets`."
            )

    # Candidate energy scales (for Option 1-style pairing): E -> K -> F0
    em_eVs = [float(x) for x in args.em_eVs.split(",") if x.strip()]
    strong_MeVs = [float(x) for x in args.strong_MeVs.split(",") if x.strip()]
    weak_GeVs = [float(x) for x in args.weak_GeVs.split(",") if x.strip()]

    print(f"\nSearching for all-4 PASS with tol={args.max_rel_err}")
    print(f"Targets:")
    print(f"- EM:     {em_keys}")
    print(f"- Strong:  {strong_keys}")
    print(f"- Weak:    {weak_keys}")
    print(f"- Gravity: {grav_keys}")
    print(f"E choices:")
    print(f"- EM eV:      {em_eVs}")
    print(f"- Strong MeV: {strong_MeVs}")
    print(f"- Weak GeV:   {weak_GeVs}")
    print(f"- Gravity K:  {args.gravity_K}\n")
    if args.gravity_f0_min is not None or args.gravity_f0_max is not None:
        print(f"Sanity: gravity F0 bounds = [{args.gravity_f0_min}, {args.gravity_f0_max}] Hz")
    if args.gravity_band != "any":
        print(f"Sanity: gravity band preset = {args.gravity_band}")
    if args.em_m_sign != "any" or args.strong_m_sign != "any" or args.weak_m_sign != "any" or args.gravity_m_sign != "any":
        print(
            "Sanity: m sign constraints = "
            f"EM {args.em_m_sign}, Strong {args.strong_m_sign}, Weak {args.weak_m_sign}, Gravity {args.gravity_m_sign}"
        )
    if args.require_order != "any":
        print(f"Sanity: require-order = {args.require_order}")
    print("")

    shown = 0
    tried = 0

    def _sign_ok(m: int, want: str) -> bool:
        if want == "any":
            return True
        if want == "positive":
            return m > 0
        if want == "negative":
            return m < 0
        if want == "nonnegative":
            return m >= 0
        if want == "nonpositive":
            return m <= 0
        raise ValueError(f"Unknown sign constraint: {want}")

    def _f0_ok(F0: float, *, min_v: float | None, max_v: float | None) -> bool:
        if min_v is not None and not (F0 >= min_v):
            return False
        if max_v is not None and not (F0 <= max_v):
            return False
        return True

    # Precompute all hits within tolerance per target key (capped).
    needed_keys = sorted(set(em_keys + strong_keys + weak_keys + grav_keys))
    hits_by_key: dict[str, list[object]] = {}
    for k in needed_keys:
        all_hits = scan_candidates(Cs=Cs, m_values=m_values, target_G=targets[k])
        tol_hits = filter_hits_by_rel_err(all_hits, max_abs_rel_err=args.max_rel_err)
        hits_by_key[k] = tol_hits[: max(1, args.max_hits_per_target)]

    for em_key in em_keys:
        for best_em in hits_by_key.get(em_key, []):
            for s_key in strong_keys:
                for best_s in hits_by_key.get(s_key, []):
                    for w_key in weak_keys:
                        for best_w in hits_by_key.get(w_key, []):
                            for g_key in grav_keys:
                                for best_g in hits_by_key.get(g_key, []):

                                    # Sanity: m sign constraints (macro/micro intuition)
                                    m_em = int(best_em.m)
                                    m_s = int(best_s.m)
                                    m_w = int(best_w.m)
                                    m_g = int(best_g.m)

                                    if not _sign_ok(m_em, args.em_m_sign):
                                        continue
                                    if not _sign_ok(m_s, args.strong_m_sign):
                                        continue
                                    if not _sign_ok(m_w, args.weak_m_sign):
                                        continue
                                    if not _sign_ok(m_g, args.gravity_m_sign):
                                        continue

                                    # Sanity: optional ordering constraints
                                    if args.require_order == "strong>weak>em":
                                        if not (m_s > m_w > m_em):
                                            continue
                                    elif args.require_order == "strong>em>weak":
                                        if not (m_s > m_em > m_w):
                                            continue
                                    elif args.require_order == "weak>strong>em":
                                        if not (m_w > m_s > m_em):
                                            continue
                                    elif args.require_order == "any":
                                        pass
                                    else:
                                        raise SystemExit(f"Unknown --require-order {args.require_order!r}")

                                    # Now sweep energy choices for the F0 side
                                    for em_eV in em_eVs:
                                        E_em_J = energy_J_from_eV(em_eV)
                                        K_em = temperature_K_from_energy_J(E_em_J)
                                        for s_MeV in strong_MeVs:
                                            E_s_J = energy_J_from_MeV(s_MeV)
                                            K_s = temperature_K_from_energy_J(E_s_J)
                                            for w_GeV in weak_GeVs:
                                                E_w_J = energy_J_from_GeV(w_GeV)
                                                K_w = temperature_K_from_energy_J(E_w_J)

                                                tried += 1
                                                # Compute F0 predictions using the fitted m values.
                                                F0_em = frequency_F0(best_em.m, K_em)
                                                F0_s = frequency_F0(best_s.m, K_s)
                                                F0_w = frequency_F0(best_w.m, K_w)
                                                F0_g = frequency_F0(best_g.m, args.gravity_K)

                                                # Sanity: F0 bounds (especially useful for gravity with CMB)
                                                if not _f0_ok(F0_em, min_v=args.em_f0_min, max_v=args.em_f0_max):
                                                    continue
                                                if not _f0_ok(F0_s, min_v=args.strong_f0_min, max_v=args.strong_f0_max):
                                                    continue
                                                if not _f0_ok(F0_w, min_v=args.weak_f0_min, max_v=args.weak_f0_max):
                                                    continue
                                                if not _f0_ok(F0_g, min_v=args.gravity_f0_min, max_v=args.gravity_f0_max):
                                                    continue

                                                shown += 1
                                                print(
                                                    f"PASS-ALL #{shown}: "
                                                    f"targets=[{em_key}, {s_key}, {w_key}, {g_key}] "
                                                    f"m=[EM {m_em}, S {m_s}, W {m_w}, G {m_g}]"
                                                )
                                                print(
                                                    f"  EM: C={best_em.C:g}, m={m_em}, rel_err={best_em.rel_err:.3e}, "
                                                    f"E={em_eV}eV -> F0={F0_em:.3e}Hz"
                                                )
                                                print(
                                                    f"  S : C={best_s.C:g}, m={m_s}, rel_err={best_s.rel_err:.3e}, "
                                                    f"E={s_MeV}MeV -> F0={F0_s:.3e}Hz"
                                                )
                                                print(
                                                    f"  W : C={best_w.C:g}, m={m_w}, rel_err={best_w.rel_err:.3e}, "
                                                    f"E={w_GeV}GeV -> F0={F0_w:.3e}Hz"
                                                )
                                                print(
                                                    f"  G : C={best_g.C:g}, m={m_g}, rel_err={best_g.rel_err:.3e}, "
                                                    f"K={args.gravity_K}K -> F0={F0_g:.3e}Hz"
                                                )
                                                print("")

                                                if shown >= args.max_results:
                                                    print(f"Stopped after max_results={args.max_results}.")
                                                    print(f"(Tried {tried} energy combos for passing target combos.)")
                                                    return 0

    print(f"Done. Found {shown} PASS-ALL configurations.")
    return 0


def _parse_int_list_csv(csv: str) -> list[int]:
    if not str(csv).strip():
        return []
    out: list[int] = []
    for raw in str(csv).split(","):
        raw = raw.strip()
        if not raw:
            continue
        out.append(int(raw))
    return out


def _num_divisors(n: int) -> int:
    """
    Simple divisor-count function for small/medium integers (used for base-scan diagnostics).
    """
    if n <= 0:
        return 0
    x = n
    p = 2
    total = 1
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            total *= e + 1
        p += 1 if p == 2 else 2  # 2, then odd primes
    if x > 1:
        total *= 2
    return total


def cmd_base_vs_alt_bases(args: argparse.Namespace) -> int:
    """
    Pre-registered style scan over candidate "base" values.

    For each base, we:
      - generate the gauge-derived C menu under the frozen invariant constructions,
      - measure how short it stays after de-duplication,
      - and check whether the strict anchor targets can be hit within a tolerance with integer m.

    This is meant to support the *selection rule* framing in `freezing_the_rules.md`,
    without turning base into an ad hoc tuning knob.
    """

    all_targets = {t.name: t for t in known_targets()}

    # Anchor targets (frozen choices; hyper is optional)
    anchors: list[tuple[str, str]] = [
        ("em", "1/alpha"),
        ("strong", "1/alpha_s_1loop_from_mZ(mH)"),
        ("weak", "1/alpha2(alpha(mZ),sin2_on_shell)"),
    ]
    if bool(getattr(args, "include_hyper", True)):
        anchors.append(("hyper", "1/alpha1_GUT(alpha(mZ),sin2_on_shell)"))

    missing = [k for _lab, k in anchors if k not in all_targets]
    if missing:
        raise SystemExit(f"Missing required anchor targets: {missing}. Run `list-targets`.")

    # Candidate bases: default frozen menu of highly-composite integers (overrideable explicitly).
    bases = _parse_int_list_csv(str(getattr(args, "bases", "")).strip())
    if not bases:
        bases = [
            12,
            24,
            36,
            48,
            60,
            120,
            180,
            240,
            360,
            420,
            480,
            600,
            720,
            840,
            1260,
            1680,
            2520,
        ]

    bases = sorted(set(int(b) for b in bases if int(b) > 0))
    if not bases:
        raise SystemExit("--bases must contain at least one positive integer base")

    tol = abs(float(getattr(args, "tol", 0.05)))
    max_nCs = int(getattr(args, "max_nCs", 10))

    # Integer m grid
    m_values = frange(args.m_min, args.m_max, args.m_step)
    m_values = sorted(set(int(round(x)) for x in m_values))

    include = tuple(s.strip() for s in str(args.include).split(",") if s.strip())

    def _sign_ok(m: int, want: str) -> bool:
        if want == "any":
            return True
        if want == "positive":
            return m > 0
        if want == "negative":
            return m < 0
        if want == "nonnegative":
            return m >= 0
        if want == "nonpositive":
            return m <= 0
        raise ValueError(f"Unknown sign constraint: {want}")

    m_sign_by_force: dict[str, str] = {
        "em": str(getattr(args, "em_m_sign", "nonnegative")),
        "strong": str(getattr(args, "strong_m_sign", "nonnegative")),
        "weak": str(getattr(args, "weak_m_sign", "nonnegative")),
        "hyper": str(getattr(args, "hyper_m_sign", "nonnegative")),
    }

    print("Base-vs-alt-bases scan (pre-registered experiment scaffolding)")
    print(f"bases = {bases}" if str(getattr(args, 'bases', '')).strip() else "bases = (default frozen highly-composite list)")
    print(f"tol(|rel_err|) = {tol:g}")
    print(f"max_nCs = {max_nCs}")
    print(f"include = {','.join(include)}")
    print(f"m range = [{min(m_values)}, {max(m_values)}] (integer grid)\n")
    print(
        "m sign constraints (micro anchors): "
        f"EM {m_sign_by_force['em']}, Strong {m_sign_by_force['strong']}, Weak {m_sign_by_force['weak']}"
        + (f", Hyper {m_sign_by_force['hyper']}" if any(f == 'hyper' for f, _k in anchors) else "")
        + "\n"
    )

    # Evaluate each base
    rows: list[tuple[int, int, int, bool, dict[str, tuple[float, float, int, float]]]] = []
    # (base, n_divisors, nCs, meets_shortness, best_by_force -> (rel_err, C, m, G))

    for base in bases:
        cand: list[tuple[str, float]] = []
        for g in standard_model_gauge_groups():
            cs = candidate_Cs_from_group(g, base=float(base), include=include)
            for k, v in cs.items():
                cand.append((f"{g.name}:{k}", float(v)))

        seen: set[float] = set()
        Cs: list[float] = []
        for _lab, c in cand:
            if c in seen:
                continue
            seen.add(c)
            Cs.append(c)

        Cs = sorted(Cs)
        nCs = len(Cs)
        meets_shortness = nCs <= max_nCs

        best_by_force: dict[str, tuple[float, float, int, float]] = {}
        for force, key in anchors:
            tgt = float(all_targets[key].value)
            want = m_sign_by_force.get(force, "any")
            m_allowed = [m for m in m_values if _sign_ok(int(m), want)]
            if not m_allowed:
                raise SystemExit(f"No m values allowed under sign constraint {want!r} for force {force!r}")
            best = scan_candidates(Cs=Cs, m_values=m_allowed, target_G=tgt)[0]
            best_by_force[force] = (float(best.rel_err), float(best.C), int(best.m), float(best.G))

        rows.append((base, _num_divisors(int(base)), nCs, meets_shortness, best_by_force))

    # Print summary table
    print("Summary (best anchor fits per base):")
    header_forces = [f for f, _k in anchors]
    print(f"{'base':>6s}  {'d(n)':>4s}  {'nCs':>3s}  {'short':>5s}  " + "  ".join(f"{f:>6s}" for f in header_forces))
    print(f"{'':>6s}  {'':>4s}  {'':>3s}  {'':>5s}  " + "  ".join(f"{'rel_err':>6s}" for _f in header_forces))

    best_candidate_base: int | None = None
    for base, dcnt, nCs, short_ok, best_by_force in rows:
        rels = [abs(best_by_force[f][0]) for f in header_forces]
        pass_all = all(r <= tol for r in rels) and short_ok
        if pass_all and best_candidate_base is None:
            best_candidate_base = base
        print(
            f"{base:6d}  {dcnt:4d}  {nCs:3d}  {('yes' if short_ok else 'no'):>5s}  "
            + "  ".join(f"{best_by_force[f][0]:+6.3g}" for f in header_forces)
        )

    print("\nDetails (per base):")
    for base, dcnt, nCs, short_ok, best_by_force in rows:
        print(f"- base={base}  d(n)={dcnt}  nCs={nCs}  short_ok={short_ok}")
        for force, key in anchors:
            rel, C, m, G = best_by_force[force]
            ok = abs(rel) <= tol
            status = "PASS" if ok else "FAIL"
            print(f"    {force:6s}  {status}  {key:36s}  rel_err={rel:+.6g}  best(C={C:g}, m={m}, G={G:.12g})")
        print("")

    if best_candidate_base is not None:
        print(f"First base meeting criteria (short + all anchors within tol): base={best_candidate_base}")
    else:
        print("No base met criteria (short + all anchors within tol) under current settings.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="topogauge")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_calc = sub.add_parser("calc", help="Compute G, F0, and the invariant relation")
    _add_common_args(p_calc)
    p_calc.set_defaults(func=cmd_calc)

    p_fits = sub.add_parser("fits", help="Fit C for several known dimensionless targets (given m)")
    p_fits.add_argument("--m", type=int, default=2, help="Index m (integer) (default: 2)")
    p_fits.set_defaults(func=cmd_fits)

    p_ex = sub.add_parser("check-example", help="Check the C=360,m=2 example vs alpha and 1/alpha")
    p_ex.set_defaults(func=cmd_check_example)

    p_scan = sub.add_parser("scan", help="Scan candidate topological C values across m and rank by target fit")
    p_scan.add_argument(
        "--target",
        default="1/alpha",
        help="Target constant name (alpha, 1/alpha, 2*pi/alpha, alpha/(2*pi)) (default: 1/alpha)",
    )
    p_scan.add_argument(
        "--set",
        dest="set_name",
        default="rotation-degrees",
        help="Candidate C set name (default: rotation-degrees). Use `list-sets` to see options.",
    )
    p_scan.add_argument(
        "--Cs",
        default="",
        help="Comma-separated override list of C values (e.g. 360,432,144). If provided, ignores --set.",
    )
    p_scan.add_argument("--m-min", type=float, default=-6.0, help="Min m (default: -6)")
    p_scan.add_argument("--m-max", type=float, default=6.0, help="Max m (default: 6)")
    p_scan.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_scan.add_argument(
        "--m-integer",
        action="store_true",
        help="Force m values to be integers (recommended for harmonic stepping).",
    )
    p_scan.add_argument(
        "--max-rel-err",
        type=float,
        default=None,
        help="Optional filter: only show hits with |rel_err| <= this value (e.g. 0.05 for 5%).",
    )
    p_scan.add_argument("--top", type=int, default=20, help="Show top N hits (default: 20)")
    p_scan.set_defaults(func=cmd_scan)

    p_solve = sub.add_parser("solve-K", help="Solve for temperature K from frequency F0 (given m)")
    p_solve.add_argument("--m", type=int, default=2, help="Index m (integer) (default: 2)")
    p_solve.add_argument("--F0", type=float, required=True, help="Frequency F0 (Hz)")
    p_solve.set_defaults(func=cmd_solve_K)

    p_targets = sub.add_parser("list-targets", help="List built-in dimensionless targets (EM/strong/gravity options)")
    p_targets.set_defaults(func=cmd_list_targets)

    p_nf = sub.add_parser("list-norm-families", help="List principled normalization-factor families (for oos-predictive)")
    p_nf.set_defaults(func=cmd_list_norm_families)

    p_rg = sub.add_parser("rg-scales", help="RG/dimensional-transmutation helper: compute Lambda_QCD from alpha_s(mu)")
    p_rg.add_argument("--alpha-s", dest="alpha_s", type=float, default=0.1179, help="Input alpha_s(mu) (default: 0.1179)")
    p_rg.add_argument("--mu-GeV", type=float, default=91.1876, help="Scale mu in GeV (default: mZ)")
    p_rg.add_argument("--n-f", type=int, default=5, help="Active flavors n_f (default: 5)")
    p_rg.set_defaults(func=cmd_rg_scales)

    p_oos_rg = sub.add_parser("oos-rg", help="RG+phi OOS: fit anchor on lattice, then test Lambda_QCD consistency across scales")
    p_oos_rg.add_argument("--suite", default="qcd-lambda-v1", help="RG suite key (default: qcd-lambda-v1)")
    p_oos_rg.add_argument("--base", type=float, default=360.0, help="Base constant to derive from (default: 360)")
    p_oos_rg.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include.",
    )
    p_oos_rg.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_oos_rg.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_oos_rg.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_oos_rg.add_argument(
        "--max-rel-err",
        type=float,
        default=0.05,
        help="Optional filter on each coupling fit before computing Lambda (default: 0.05)",
    )
    p_oos_rg.add_argument("--top", type=int, default=10, help="Show top N candidates (default: 10)")
    p_oos_rg.set_defaults(func=cmd_oos_rg)

    p_oos = sub.add_parser("oos-report", help="Out-of-sample report: run a frozen target suite against strict gauge-derived C")
    p_oos.add_argument("--base", type=float, default=360.0, help="Base constant to derive from (default: 360)")
    p_oos.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include.",
    )
    p_oos.add_argument(
        "--suite",
        choices=["v1", "v2", "v3", "v4", "v5", "v6", "v7"],
        default="v1",
        help="Frozen OOS suite to run (default: v1)",
    )
    p_oos.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_oos.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_oos.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_oos.add_argument("--max-rel-err", type=float, default=0.02, help="Tolerance on |rel_err| (default: 0.02)")
    p_oos.set_defaults(func=cmd_oos_report)

    p_oos_pred = sub.add_parser(
        "oos-predictive",
        help="Predictive OOS: freeze one best-fit C per force (from strict anchors), then evaluate OOS targets with C fixed",
    )
    p_oos_pred.add_argument("--base", type=float, default=360.0, help="Base constant to derive from (default: 360)")
    p_oos_pred.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include.",
    )
    p_oos_pred.add_argument("--suite", choices=["v1", "v2", "v3"], default="v1", help="Predictive OOS suite to run (default: v1)")
    p_oos_pred.add_argument(
        "--force",
        choices=["all", "em", "strong", "weak", "hyper", "gravity"],
        default="all",
        help="Which force group(s) to run (default: all)",
    )
    p_oos_pred.add_argument(
        "--norm-family",
        choices=sorted(normalization_families().keys()),
        default="none",
        help="Principled normalization family to apply to targets before lattice fitting (default: none).",
    )
    p_oos_pred.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_oos_pred.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_oos_pred.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_oos_pred.add_argument("--max-rel-err", type=float, default=0.02, help="Tolerance on |rel_err| (default: 0.02)")
    p_oos_pred.set_defaults(func=cmd_oos_predictive)

    p_oos_pred_rg = sub.add_parser(
        "oos-predictive-rg",
        help="Predictive OOS with deterministic within-band RG running (strong + EM + weak)",
    )
    p_oos_pred_rg.add_argument("--suite", choices=["v1", "v2", "v3"], default="v1", help="Predictive OOS suite to run (default: v1)")
    p_oos_pred_rg.add_argument(
        "--force",
        choices=["all", "em", "strong", "weak", "hyper"],
        default="all",
        help="Force group(s) to run (default: all)",
    )
    p_oos_pred_rg.add_argument("--base", type=float, default=360.0, help="Base constant to derive from (default: 360)")
    p_oos_pred_rg.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include.",
    )
    p_oos_pred_rg.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_oos_pred_rg.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_oos_pred_rg.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_oos_pred_rg.add_argument("--max-rel-err", type=float, default=0.02, help="Tolerance on |rel_err| (default: 0.02)")
    p_oos_pred_rg.add_argument(
        "--runner",
        choices=["auto", "qed_1loop", "qed_pdg_mZ", "ew_sm_1loop", "1loop_nf5", "1loop_nf56", "2loop_nf5", "2loop_nf56"],
        default="auto",
        help="Deterministic running prescription to use (strong: QCD; EM: QED; weak: EW) (default: auto).",
    )
    p_oos_pred_rg.add_argument(
        "--Q0-GeV",
        dest="Q0_GeV",
        type=float,
        default=None,
        help="EM only: reference scale for the QED threshold model (default: m_e). Strong ignores this (Q0 from anchor key).",
    )
    p_oos_pred_rg.add_argument(
        "--Q-switch-GeV",
        dest="Q_switch_GeV",
        type=float,
        default=172.76,
        help="Flavor-switch scale for nf56 variants (default: mt≈172.76 GeV).",
    )
    p_oos_pred_rg.add_argument(
        "--steps-per-unit-log",
        type=int,
        default=500,
        help="Integrator resolution for 2-loop running (default: 500).",
    )
    p_oos_pred_rg.add_argument(
        "--targets",
        default="",
        help="Optional comma-separated override list of STRONG target keys to evaluate (default: suite v1 strong targets).",
    )
    p_oos_pred_rg.set_defaults(func=cmd_oos_predictive_rg)

    p_oos_ew = sub.add_parser(
        "oos-ew-mix",
        help="EW mixing cross-check: fit alpha2^-1 and alpha1_GUT^-1 on the lattice, SM 1-loop run both, and compare derived sin^2thetaW(Q)",
    )
    p_oos_ew.add_argument("--base", type=float, default=360.0, help="Base constant to derive from (default: 360)")
    p_oos_ew.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include.",
    )
    p_oos_ew.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_oos_ew.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_oos_ew.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_oos_ew.add_argument("--max-rel-err", type=float, default=0.02, help="Tolerance on |rel_err| (default: 0.02)")
    p_oos_ew.add_argument(
        "--scales",
        default="mW,mH,1TeV,10TeV",
        help="Comma-separated scale labels to evaluate (default: mW,mH,1TeV,10TeV).",
    )
    p_oos_ew.set_defaults(func=cmd_oos_ew_mix)

    p_ew_sin2 = sub.add_parser(
        "ew-sin2",
        help="Predict sin^2thetaW(Q) from lattice-quantized alpha2^-1/alpha1_GUT^-1 anchors + SM/MSSM 1-loop running; optionally compare to supplied measurements",
    )
    p_ew_sin2.add_argument("--model", choices=["sm", "mssm"], default="sm", help="1-loop beta set for running (default: sm)")
    p_ew_sin2.add_argument(
        "--method",
        choices=["gauge_1loop", "gammaZ_1loop"],
        default="gammaZ_1loop",
        help="Prediction method (default: gammaZ_1loop).",
    )
    p_ew_sin2.add_argument("--base", type=float, default=360.0, help="Base constant to derive from (default: 360)")
    p_ew_sin2.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include.",
    )
    p_ew_sin2.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_ew_sin2.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_ew_sin2.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_ew_sin2.add_argument("--max-rel-err", type=float, default=0.02, help="Tolerance on |rel_err| when comparing to measurements (default: 0.02)")
    p_ew_sin2.add_argument(
        "--scales",
        default="mW,mH,1TeV,10TeV",
        help="Comma-separated scale labels to evaluate (default: mW,mH,1TeV,10TeV).",
    )
    p_ew_sin2.add_argument(
        "--measurement",
        action="append",
        nargs=3,
        metavar=("LABEL", "SIN2", "SIGMA"),
        default=[],
        help="Optional measurement triple: <label> <sin2> <sigma_or_0>. Label must match one of --scales entries.",
    )
    p_ew_sin2.set_defaults(func=cmd_ew_sin2)

    p_oos_ew_sin2 = sub.add_parser(
        "oos-ew-sin2",
        help="OOS: compare predicted sin^2thetaW(Q) (from lattice anchors + 1-loop running) to external registry-provided targets (tgt_sin2thetaW(...))",
    )
    p_oos_ew_sin2.add_argument(
        "--suite",
        choices=["registry-all", *sorted(ew_sin2_suites().keys())],
        default="registry-all",
        help="Which external sin2 target suite to run (default: registry-all). Use ew-independent-v1 for a frozen membership suite.",
    )
    p_oos_ew_sin2.add_argument("--model", choices=["sm", "mssm"], default="sm", help="1-loop beta set for running (default: sm)")
    p_oos_ew_sin2.add_argument(
        "--method",
        choices=["gauge_1loop", "gammaZ_1loop", "zpole_kappa_approx"],
        default="gammaZ_1loop",
        help="Prediction method (default: gammaZ_1loop). Use zpole_kappa_approx for Z-pole sin^2θ_eff^lept targets.",
    )
    p_oos_ew_sin2.add_argument("--base", type=float, default=360.0, help="Base constant to derive from (default: 360)")
    p_oos_ew_sin2.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include.",
    )
    p_oos_ew_sin2.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_oos_ew_sin2.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_oos_ew_sin2.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_oos_ew_sin2.add_argument("--max-rel-err", type=float, default=0.02, help="Tolerance on |rel_err| (default: 0.02)")
    p_oos_ew_sin2.add_argument(
        "--z-max",
        type=float,
        default=None,
        help="Optional sigma-based pass rule: require |z| <= z-max when sigma is available. "
        "Defaults: ew-independent-v1 -> 3; ew-independent-v2/v3 -> 2 (if not provided).",
    )
    p_oos_ew_sin2.add_argument(
        "--sigma-theory",
        type=float,
        default=None,
        help="Optional absolute theory/model uncertainty (added in quadrature to experimental sigma when computing z). "
        "Defaults: ew-independent-v* -> 0 (after the κ(Q) model upgrade).",
    )
    p_oos_ew_sin2.add_argument(
        "--scheme-prefix",
        type=str,
        default="",
        help="Optional scheme guardrail: require that each target's scheme string starts with this prefix. "
        "For --suite ew-independent-v*, defaults to 'sin2thetaW_eff:'.",
    )
    p_oos_ew_sin2.add_argument(
        "--report-theory-threshold",
        action="store_true",
        help="Report the minimum sigma_theory required (per target and overall) to satisfy the |z|<=z-max gate. "
        "Works when sigma is available (e.g. ew-independent-v1).",
    )
    p_oos_ew_sin2.add_argument(
        "--ci",
        action="store_true",
        help="Exit non-zero if any target FAILs (useful for CI). Note: ew-independent-v* suites are gates by default.",
    )
    p_oos_ew_sin2.set_defaults(func=cmd_oos_ew_sin2)

    p_oos_steps = sub.add_parser(
        "oos-steps",
        help="Step-signal OOS (C-independent): test whether anchor/target ratios land on integer Δm steps of φ",
    )
    p_oos_steps.add_argument("--suite", choices=["v1"], default="v1", help="Step-signal suite to run (default: v1)")
    p_oos_steps.add_argument(
        "--force",
        choices=["all", "em", "strong", "weak", "gravity"],
        default="all",
        help="Which force group(s) to run (default: all)",
    )
    p_oos_steps.add_argument(
        "--max-ratio-err",
        type=float,
        default=0.02,
        help="Tolerance on ratio error when snapping dm to nearest integer (default: 0.02)",
    )
    p_oos_steps.set_defaults(func=cmd_oos_steps)

    p_gbands = sub.add_parser("list-gravity-bands", help="List built-in gravity-wave frequency band presets")
    p_gbands.set_defaults(func=cmd_list_gravity_bands)

    p_fp = sub.add_parser("list-frequency-presets", help="List Option-2 frequency presets for EM/strong/weak")
    p_fp.set_defaults(func=cmd_list_frequency_presets)

    p_gc = sub.add_parser("list-gauge-Cs", help="List gauge-derived candidate C values from group invariants")
    p_gc.add_argument("--base", type=float, default=360.0, help="Base constant to derive from (default: 360)")
    p_gc.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include (see gauge_invariants.md).",
    )
    p_gc.set_defaults(func=cmd_list_gauge_candidates)

    p_sgc = sub.add_parser("scan-gauge-Cs", help="Scan gauge-derived C candidates across integer m for a target coupling")
    p_sgc.add_argument("--base", type=float, default=360.0, help="Base constant to derive from (default: 360)")
    p_sgc.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include (see gauge_invariants.md).",
    )
    p_sgc.add_argument("--target", default="1/alpha", help="Target key (see list-targets)")
    p_sgc.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_sgc.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_sgc.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_sgc.add_argument("--max-rel-err", type=float, default=0.05, help="Tolerance on |rel_err| (default: 0.05)")
    p_sgc.add_argument("--top", type=int, default=20, help="Show top N hits (default: 20)")
    p_sgc.set_defaults(func=cmd_scan_gauge_Cs)

    p_bases = sub.add_parser(
        "base-vs-alt-bases",
        help="Pre-registered-style scan: compare candidate base values by gauge-derived C menu size + strict-anchor hit quality",
    )
    p_bases.add_argument(
        "--bases",
        default="",
        help="Optional comma-separated list of integer bases to test (overrides the default frozen highly-composite list).",
    )
    p_bases.add_argument(
        "--no-hyper",
        dest="include_hyper",
        action="store_false",
        help="Disable the hypercharge anchor 1/alpha1_GUT(...) in the base-scan (default: include hyper).",
    )
    p_bases.add_argument(
        "--em-m-sign",
        choices=["any", "positive", "negative", "nonnegative", "nonpositive"],
        default="nonnegative",
        help="Constraint on the EM anchor's fitted m (default: nonnegative).",
    )
    p_bases.add_argument(
        "--strong-m-sign",
        choices=["any", "positive", "negative", "nonnegative", "nonpositive"],
        default="nonnegative",
        help="Constraint on the strong anchor's fitted m (default: nonnegative).",
    )
    p_bases.add_argument(
        "--weak-m-sign",
        choices=["any", "positive", "negative", "nonnegative", "nonpositive"],
        default="nonnegative",
        help="Constraint on the weak anchor's fitted m (default: nonnegative).",
    )
    p_bases.add_argument(
        "--hyper-m-sign",
        choices=["any", "positive", "negative", "nonnegative", "nonpositive"],
        default="nonnegative",
        help="Constraint on the hypercharge anchor's fitted m (default: nonnegative).",
    )
    p_bases.add_argument("--tol", type=float, default=0.05, help="Anchor fit tolerance on |rel_err| (default: 0.05)")
    p_bases.add_argument(
        "--max-nCs",
        dest="max_nCs",
        type=int,
        default=10,
        help="Shortness criterion: require deduplicated gauge-derived C menu size <= max-nCs (default: 10).",
    )
    p_bases.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include (kept fixed for comparability).",
    )
    p_bases.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_bases.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_bases.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_bases.set_defaults(func=cmd_base_vs_alt_bases)

    p_qg = sub.add_parser(
        "sweep-quantum-gravity",
        help="Sweep gravity mass scale (GeV) and find which scales allow gauge-C + integer-m fits under a GW band constraint",
    )
    p_qg.add_argument("--base", type=float, default=360.0, help="Base constant to derive Cs from (default: 360)")
    p_qg.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include.",
    )
    p_qg.add_argument("--max-rel-err", type=float, default=0.05, help="Tolerance on |rel_err| (default: 0.05)")
    p_qg.add_argument("--em-target", default="1/alpha", help="EM target key (default: 1/alpha)")
    p_qg.add_argument(
        "--strong-target",
        default="1/alpha_s_1loop_from_mZ(mH)",
        help="Strong target key (default: 1/alpha_s_1loop_from_mZ(mH))",
    )
    p_qg.add_argument(
        "--weak-target",
        default="1/alpha2(alpha(mZ),sin2_on_shell)",
        help="Weak target key (default: 1/alpha2(alpha(mZ),sin2_on_shell))",
    )
    p_qg.add_argument(
        "--gravity-mode",
        default="inverse",
        choices=["inverse", "direct"],
        help="Use 1/alpha_G(mass) or alpha_G(mass) as the gravity coupling target.",
    )
    p_qg.add_argument("--gravity-K", type=float, default=2.725, help="Gravity temperature K (default: CMB)")
    p_qg.add_argument("--gravity-band", default="cmb", choices=["cmb", "pta", "lisa", "ligo"], help="GW band constraint")
    p_qg.add_argument("--m-min", type=int, default=-256, help="Min m considered (default: -256)")
    p_qg.add_argument("--m-max", type=int, default=256, help="Max m considered (default: 256)")
    p_qg.add_argument("--scale-min-GeV", type=float, default=1e3, help="Min gravity mass scale in GeV (default: 1e3)")
    p_qg.add_argument("--scale-max-GeV", type=float, default=1e6, help="Max gravity mass scale in GeV (default: 1e6)")
    p_qg.add_argument("--n-scales", type=int, default=121, help="Number of log-spaced scales to test (default: 121)")
    p_qg.add_argument("--top", type=int, default=25, help="Show top N passing scales (default: 25)")
    p_qg.set_defaults(func=cmd_sweep_quantum_gravity)

    p_gut = sub.add_parser("gut-run", help="1-loop running test for alpha1_GUT, alpha2, alpha3 convergence")
    p_gut.add_argument("--model", choices=["sm", "mssm"], default="sm", help="1-loop beta set (default: sm)")
    p_gut.add_argument("--Q-min-GeV", type=float, default=1e2, help="Min Q in GeV (default: 1e2)")
    p_gut.add_argument("--Q-max-GeV", type=float, default=1e18, help="Max Q in GeV (default: 1e18)")
    p_gut.add_argument("--n", type=int, default=200, help="Number of log-spaced Q points (default: 200)")
    p_gut.add_argument("--top", type=int, default=12, help="Rows to print around best (default: 12)")
    p_gut.set_defaults(func=cmd_gut_run)

    p_gut_lat = sub.add_parser(
        "gut-run-lattice",
        help="1-loop convergence test like gut-run, but with lattice-quantized inputs for alpha1_GUT^-1, alpha2^-1, alpha3^-1",
    )
    p_gut_lat.add_argument("--model", choices=["sm", "mssm"], default="sm", help="1-loop beta set (default: sm)")
    p_gut_lat.add_argument("--Q-min-GeV", type=float, default=1e2, help="Min Q in GeV (default: 1e2)")
    p_gut_lat.add_argument("--Q-max-GeV", type=float, default=1e18, help="Max Q in GeV (default: 1e18)")
    p_gut_lat.add_argument("--n", type=int, default=200, help="Number of log-spaced Q points (default: 200)")
    p_gut_lat.add_argument("--base", type=float, default=360.0, help="Base constant to derive from (default: 360)")
    p_gut_lat.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include.",
    )
    p_gut_lat.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_gut_lat.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_gut_lat.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_gut_lat.set_defaults(func=cmd_gut_run_lattice)

    p_pg = sub.add_parser(
        "pair-forces-gaugeCs",
        help="Pair forces using only gauge-derived C candidates; Option-2 for EM/strong/weak; gravity fixed-K with optional GW band filter",
    )
    p_pg.add_argument("--base", type=float, default=360.0, help="Base constant to derive from (default: 360)")
    p_pg.add_argument(
        "--include",
        default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)",
        help="Comma-separated gauge C constructions to include.",
    )
    p_pg.add_argument("--C360-only", action="store_true", help="Restrict to C=360 only.")
    p_pg.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_pg.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_pg.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_pg.add_argument("--max-rel-err", type=float, default=0.05, help="Tolerance on coupling fit (default: 0.05)")
    p_pg.add_argument("--max-hits", type=int, default=30, help="Max hits per force target to consider (default: 30)")
    p_pg.add_argument("--max-results", type=int, default=20, help="Stop after this many results.")

    p_pg.add_argument("--em-target", default="1/alpha", help="EM target key (default: 1/alpha)")
    p_pg.add_argument("--strong-target", default="1/alpha_s_1loop_from_mZ(mH)", help="Strong target key (default: 1/alpha_s_1loop_from_mZ(mH))")
    p_pg.add_argument("--weak-target", default="1/alpha2(alpha(mZ),sin2_on_shell)", help="Weak target key (default: 1/alpha2(alpha(mZ),sin2_on_shell))")
    p_pg.add_argument(
        "--gravity-targets",
        default="1/alpha_G(p)",
        help="Comma-separated gravity target keys to try (default: 1/alpha_G(p))",
    )

    p_pg.add_argument("--em-F0", type=float, default=None, help="EM F0 (Hz)")
    p_pg.add_argument("--em-preset", default="em-lyman-alpha", help="EM frequency preset key (default: em-lyman-alpha)")
    p_pg.add_argument("--strong-F0", type=float, default=None, help="Strong F0 (Hz)")
    p_pg.add_argument("--strong-preset", default="strong-QCD-200MeV", help="Strong frequency preset key (default: strong-QCD-200MeV)")
    p_pg.add_argument("--weak-F0", type=float, default=None, help="Weak F0 (Hz)")
    p_pg.add_argument("--weak-preset", default="weak-W-80.379GeV", help="Weak frequency preset key (default: weak-W-80.379GeV)")

    p_pg.add_argument("--gravity-K", type=float, default=2.725, help="Gravity K (default: CMB)")
    p_pg.add_argument("--gravity-band", default="any", choices=["any", "ligo", "lisa", "pta", "cmb"], help="GW band filter")
    p_pg.add_argument("--gravity-f0-min", type=float, default=None, help="Override gravity F0 min (Hz)")
    p_pg.add_argument("--gravity-f0-max", type=float, default=None, help="Override gravity F0 max (Hz)")

    p_pg.add_argument("--em-m-sign", choices=["any", "positive", "negative"], default="any", help="EM m sign constraint")
    p_pg.add_argument("--strong-m-sign", choices=["any", "positive", "negative"], default="any", help="Strong m sign constraint")
    p_pg.add_argument("--weak-m-sign", choices=["any", "positive", "negative"], default="any", help="Weak m sign constraint")
    p_pg.add_argument("--gravity-m-sign", choices=["any", "positive", "negative"], default="any", help="Gravity m sign constraint")

    p_pg.add_argument("--em-K-min", type=float, default=None, help="EM K min")
    p_pg.add_argument("--em-K-max", type=float, default=None, help="EM K max")
    p_pg.add_argument("--strong-K-min", type=float, default=None, help="Strong K min")
    p_pg.add_argument("--strong-K-max", type=float, default=None, help="Strong K max")
    p_pg.add_argument("--weak-K-min", type=float, default=None, help="Weak K min")
    p_pg.add_argument("--weak-K-max", type=float, default=None, help="Weak K max")
    p_pg.add_argument("--gravity-K-min", type=float, default=None, help="Gravity K min")
    p_pg.add_argument("--gravity-K-max", type=float, default=None, help="Gravity K max")
    p_pg.set_defaults(func=cmd_pair_forces_gaugeCs)

    p_scan_all = sub.add_parser("scan-all", help="Scan the same (C,m) grid against every built-in target")
    p_scan_all.add_argument(
        "--set",
        dest="set_name",
        default="rotation-degrees",
        help="Candidate C set name (default: rotation-degrees). Use `list-sets` to see options.",
    )
    p_scan_all.add_argument(
        "--Cs",
        default="",
        help="Comma-separated override list of C values (e.g. 360,432,144). If provided, ignores --set.",
    )
    p_scan_all.add_argument("--m-min", type=float, default=-6.0, help="Min m (default: -6)")
    p_scan_all.add_argument("--m-max", type=float, default=6.0, help="Max m (default: 6)")
    p_scan_all.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_scan_all.add_argument("--m-integer", action="store_true", help="Force m values to be integers.")
    p_scan_all.add_argument(
        "--max-rel-err",
        type=float,
        default=0.05,
        help="Pass/fail threshold on best |rel_err| (default: 0.05 = 5%).",
    )
    p_scan_all.set_defaults(func=cmd_scan_all)

    p_sets = sub.add_parser("list-sets", help="List built-in candidate C sets")
    p_sets.set_defaults(func=cmd_list_sets)

    p_pair = sub.add_parser("pair-forces", help="Try pairing the 4 forces across both formulas (exploratory)")
    p_pair.add_argument(
        "--set",
        dest="set_name",
        default="rotation-degrees",
        help="Candidate C set name (default: rotation-degrees).",
    )
    p_pair.add_argument("--Cs", default="", help="Comma-separated override list of C values.")
    p_pair.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_pair.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_pair.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_pair.add_argument("--max-rel-err", type=float, default=0.05, help="Tolerance on G fit (default: 0.05)")

    # Default target choices (can change)
    p_pair.add_argument("--em-target", default="1/alpha", help="EM target key (default: 1/alpha)")
    p_pair.add_argument("--strong-target", default="1/alpha_s_1loop_from_mZ(mH)", help="Strong target key (default: 1/alpha_s_1loop_from_mZ(mH))")
    p_pair.add_argument("--weak-target", default="1/alpha2(alpha(mZ),sin2_on_shell)", help="Weak target key (default: 1/alpha2(alpha(mZ),sin2_on_shell))")
    p_pair.add_argument("--gravity-target", default="1/alpha_G(p)", help="Gravity target key (default: 1/alpha_G(p))")

    # Energy scale knobs for the quantum forces
    p_pair.add_argument("--em-eV", type=float, default=13.6, help="EM energy scale in eV (default: 13.6 eV)")
    p_pair.add_argument("--strong-MeV", type=float, default=200.0, help="Strong energy scale in MeV (default: 200 MeV)")
    p_pair.add_argument("--weak-GeV", type=float, default=80.379, help="Weak energy scale in GeV (default: mW)")

    # Macro temperature for gravity
    p_pair.add_argument("--gravity-K", type=float, default=2.725, help="Gravity temperature in K (default: CMB ~2.725K)")
    p_pair.set_defaults(func=cmd_pair_forces)

    p_pair_all = sub.add_parser(
        "pair-forces-all",
        help="Explore many target/energy options and print configurations where all 4 forces pass the G-fit tolerance",
    )
    p_pair_all.add_argument(
        "--set",
        dest="set_name",
        default="rotation-degrees",
        help="Candidate C set name (default: rotation-degrees).",
    )
    p_pair_all.add_argument("--Cs", default="", help="Comma-separated override list of C values.")
    p_pair_all.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_pair_all.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_pair_all.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_pair_all.add_argument("--max-rel-err", type=float, default=0.05, help="Tolerance on G fit (default: 0.05)")
    p_pair_all.add_argument("--max-results", type=int, default=20, help="Stop after this many PASS-ALL results.")
    p_pair_all.add_argument(
        "--max-hits-per-target",
        type=int,
        default=20,
        help="Limit number of (C,m) hits per target to consider (default: 20).",
    )

    p_pair_all.add_argument(
        "--em-targets",
        default="1/alpha",
        help="Comma-separated EM target keys to try.",
    )
    p_pair_all.add_argument(
        "--strong-targets",
        default="1/alpha_s_1loop_from_mZ(mH),1/alpha_s(mZ),1/alpha_s_1loop(10GeV)",
        help="Comma-separated strong target keys to try.",
    )
    p_pair_all.add_argument(
        "--weak-targets",
        default="1/alpha2(alpha(mZ),sin2_on_shell),1/alpha_w(mZ)",
        help="Comma-separated weak target keys to try.",
    )
    p_pair_all.add_argument(
        "--gravity-targets",
        default="alpha_G(p),1/alpha_G(p),alpha_G(e),1/alpha_G(e)",
        help="Comma-separated gravity target keys to try.",
    )

    p_pair_all.add_argument("--em-eVs", default="1,13.6,2.0", help="Comma-separated EM energy scales in eV.")
    p_pair_all.add_argument("--strong-MeVs", default="200,938", help="Comma-separated strong energy scales in MeV.")
    p_pair_all.add_argument("--weak-GeVs", default="80.379,91.1876,125", help="Comma-separated weak energy scales in GeV.")
    p_pair_all.add_argument("--gravity-K", type=float, default=2.725, help="Gravity temperature in K (default: CMB).")
    # Sanity filters (all optional; defaults are permissive)
    p_pair_all.add_argument(
        "--em-m-sign",
        default="any",
        choices=["any", "positive", "negative", "nonnegative", "nonpositive"],
        help="Require EM m sign (default: any).",
    )
    p_pair_all.add_argument(
        "--strong-m-sign",
        default="any",
        choices=["any", "positive", "negative", "nonnegative", "nonpositive"],
        help="Require strong m sign (default: any).",
    )
    p_pair_all.add_argument(
        "--weak-m-sign",
        default="any",
        choices=["any", "positive", "negative", "nonnegative", "nonpositive"],
        help="Require weak m sign (default: any).",
    )
    p_pair_all.add_argument(
        "--gravity-m-sign",
        default="any",
        choices=["any", "positive", "negative", "nonnegative", "nonpositive"],
        help="Require gravity m sign (default: any).",
    )
    p_pair_all.add_argument(
        "--require-order",
        default="any",
        choices=["any", "strong>weak>em", "strong>em>weak", "weak>strong>em"],
        help="Optional ordering constraint among positive-m forces (default: any).",
    )

    # F0 bounds (Hz). Use these to reject absurd macro/micro pairings.
    p_pair_all.add_argument("--em-f0-min", type=float, default=None, help="Min allowed EM F0 (Hz).")
    p_pair_all.add_argument("--em-f0-max", type=float, default=None, help="Max allowed EM F0 (Hz).")
    p_pair_all.add_argument("--strong-f0-min", type=float, default=None, help="Min allowed strong F0 (Hz).")
    p_pair_all.add_argument("--strong-f0-max", type=float, default=None, help="Max allowed strong F0 (Hz).")
    p_pair_all.add_argument("--weak-f0-min", type=float, default=None, help="Min allowed weak F0 (Hz).")
    p_pair_all.add_argument("--weak-f0-max", type=float, default=None, help="Max allowed weak F0 (Hz).")
    p_pair_all.add_argument("--gravity-f0-min", type=float, default=None, help="Min allowed gravity F0 (Hz).")
    p_pair_all.add_argument("--gravity-f0-max", type=float, default=None, help="Max allowed gravity F0 (Hz).")
    p_pair_all.add_argument(
        "--gravity-band",
        default="any",
        choices=["any", "ligo", "lisa", "pta", "cmb"],
        help="Shortcut preset for gravity-wave frequency band (sets gravity F0 bounds if none are provided).",
    )
    p_pair_all.set_defaults(func=cmd_pair_forces_all)

    p_spec = sub.add_parser("spectrum", help="Phi-lattice spectrum: all targets organized by m-value, revealing lattice topology")
    p_spec.add_argument("--base", type=float, default=360.0, help="Gauge base (default: 360)")
    p_spec.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_spec.add_argument("--m-min", type=int, default=-60, help="Min m (default: -60)")
    p_spec.add_argument("--m-max", type=int, default=60, help="Max m (default: 60)")
    p_spec.add_argument("--max-rel-err", type=float, default=0.05, help="Max |rel_err| to include (default: 0.05)")
    p_spec.add_argument("--filter", default="", help="Filter target names (case-insensitive substring)")
    p_spec.add_argument("--no-summary", action="store_true", help="Skip sector summary at end")
    p_spec.set_defaults(func=cmd_spectrum)

    p_gut_mssm = sub.add_parser("gut-compare", help="Compare SM vs MSSM 1-loop GUT convergence (with lattice-quantized inputs)")
    p_gut_mssm.add_argument("--base", type=float, default=360.0, help="Gauge base (default: 360)")
    p_gut_mssm.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_mssm.set_defaults(func=cmd_gut_run_mssm)

    p_gut_traj = sub.add_parser("gut-trajectory", help="GUT trajectory: show RG flow lattice addresses from mZ to M_Planck + GUT scale checks")
    p_gut_traj.add_argument("--base", type=float, default=360.0, help="Gauge base (default: 360)")
    p_gut_traj.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_traj.set_defaults(func=cmd_gut_trajectory)

    p_gut_val = sub.add_parser("gut-validate", help="GUT validation: independent predictions, 2-loop tightened search, threshold corrections, Fibonacci structure")
    p_gut_val.add_argument("--base", type=float, default=360.0, help="Gauge base (default: 360)")
    p_gut_val.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_val.set_defaults(func=cmd_gut_validate)

    p_gut_sig = sub.add_parser("gut-significance", help="Statistical significance: base permutation test, C-clustering, pre-registered predictions")
    p_gut_sig.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_sig.set_defaults(func=cmd_gut_significance)

    p_gut_deep = sub.add_parser("gut-deep", help="Deep analysis: address table, CKM/PMNS, lattice operations, 360-family, n_gap")
    p_gut_deep.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_deep.set_defaults(func=cmd_gut_deep)

    p_gut_predict = sub.add_parser("gut-predict", help="Extended predictions: RG lattice energies, m_W consistency, neutrinos, duality, vacancies")
    p_gut_predict.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_predict.set_defaults(func=cmd_gut_predict)

    p_gut_explore = sub.add_parser("gut-explore", help="Round 4: lattice arithmetic, duality axis, α_s data, thermodynamics, 360 number theory")
    p_gut_explore.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_explore.set_defaults(func=cmd_gut_explore)

    p_gut_spectrum = sub.add_parser("gut-spectrum", help="Round 5: fermion spectrum, symmetry breaking, new ratios, RG flow, outliers")
    p_gut_spectrum.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_spectrum.set_defaults(func=cmd_gut_spectrum)

    p_gut_pheno = sub.add_parser("gut-pheno", help="Round 6: proton decay, dark matter, muon g-2, information theory, neutrinos")
    p_gut_pheno.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_pheno.set_defaults(func=cmd_gut_pheno)

    p_gut_struct = sub.add_parser("gut-struct", help="Round 7: lattice renormalization, phase transitions, anomaly matching, hadron spectrum, topology")
    p_gut_struct.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_struct.set_defaults(func=cmd_gut_struct)

    p_gut_precision = sub.add_parser("gut-precision", help="Round 8: EW precision, action principle, CP violation, inflation")
    p_gut_precision.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_precision.set_defaults(func=cmd_gut_precision)

    p_gut_nuclear = sub.add_parser("gut-nuclear", help="Round 9: nuclear physics, astrophysics, bootstrap, coupling ratios, QCD")
    p_gut_nuclear.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_nuclear.set_defaults(func=cmd_gut_nuclear)

    p_gut_final = sub.add_parser("gut-final", help="Round 10: Yukawa hierarchy, vacuum energy, grid symmetries, scorecard")
    p_gut_final.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_final.set_defaults(func=cmd_gut_final)

    p_gut_deep2 = sub.add_parser("gut-deep2", help="Round 11: BH thermo, lepton universality, base-360 derivation, 91/588 analysis")
    p_gut_deep2.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_deep2.set_defaults(func=cmd_gut_deep2)

    p_gut_predict2 = sub.add_parser("gut-predict2", help="Round 12: math constants, condensed matter, sharp predictions")
    p_gut_predict2.add_argument("--include", default="base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter)", help="C constructions")
    p_gut_predict2.set_defaults(func=cmd_gut_predict2)

    p_opt2 = sub.add_parser("pair-forces-option2", help="Option-2 pairing: use F0 inputs for EM/strong/weak and solve K")
    p_opt2.add_argument("--set", dest="set_name", default="octave-union", help="Candidate C set name (default: octave-union).")
    p_opt2.add_argument("--Cs", default="", help="Comma-separated override list of C values.")
    p_opt2.add_argument("--m-min", type=float, default=-256.0, help="Min m (default: -256)")
    p_opt2.add_argument("--m-max", type=float, default=256.0, help="Max m (default: 256)")
    p_opt2.add_argument("--m-step", type=float, default=1.0, help="Step for m (default: 1)")
    p_opt2.add_argument("--max-rel-err", type=float, default=0.05, help="Tolerance on G fit (default: 0.05)")
    p_opt2.add_argument("--max-hits", type=int, default=30, help="Max hits per target to consider (default: 30)")
    p_opt2.add_argument("--max-results", type=int, default=20, help="Stop after this many solutions.")

    p_opt2.add_argument("--em-target", default="1/alpha", help="EM coupling target key (default: 1/alpha)")
    p_opt2.add_argument(
        "--strong-target",
        default="1/alpha_s_1loop_from_mZ(mH)",
        help="Strong coupling target key (default: 1/alpha_s_1loop_from_mZ(mH))",
    )
    p_opt2.add_argument(
        "--weak-target",
        default="1/alpha2(alpha(mZ),sin2_on_shell)",
        help="Weak coupling target key (default: 1/alpha2(alpha(mZ),sin2_on_shell))",
    )
    p_opt2.add_argument("--gravity-target", default="1/alpha_G(p)", help="Gravity coupling target key (default: 1/alpha_G(p))")

    # Frequencies: either provide numeric F0 or a preset key.
    p_opt2.add_argument("--em-F0", type=float, default=None, help="EM phenomenon frequency F0 (Hz)")
    p_opt2.add_argument("--em-preset", default="em-lyman-alpha", help="EM frequency preset key (default: em-lyman-alpha)")
    p_opt2.add_argument("--strong-F0", type=float, default=None, help="Strong phenomenon frequency F0 (Hz)")
    p_opt2.add_argument("--strong-preset", default="strong-QCD-200MeV", help="Strong frequency preset key (default: strong-QCD-200MeV)")
    p_opt2.add_argument("--weak-F0", type=float, default=None, help="Weak phenomenon frequency F0 (Hz)")
    p_opt2.add_argument("--weak-preset", default="weak-W-80.379GeV", help="Weak frequency preset key (default: weak-W-80.379GeV)")

    p_opt2.add_argument("--gravity-K", type=float, default=2.725, help="Gravity K if gravity F0 not provided (default: CMB)")
    p_opt2.add_argument("--gravity-F0", type=float, default=None, help="Optional gravity F0 (Hz) to solve K instead of using gravity-K")
    p_opt2.add_argument("--gravity-preset", default=None, help="Optional gravity frequency preset key (can reuse EM presets)")

    # Optional sanity constraints
    p_opt2.add_argument("--em-m-sign", choices=["any", "positive", "negative"], default="any", help="EM m sign constraint")
    p_opt2.add_argument("--strong-m-sign", choices=["any", "positive", "negative"], default="any", help="Strong m sign constraint")
    p_opt2.add_argument("--weak-m-sign", choices=["any", "positive", "negative"], default="any", help="Weak m sign constraint")
    p_opt2.add_argument("--gravity-m-sign", choices=["any", "positive", "negative"], default="any", help="Gravity m sign constraint")
    p_opt2.add_argument("--K-min", type=float, default=None, help="Optional minimum K allowed (applies to all forces)")
    p_opt2.add_argument("--K-max", type=float, default=None, help="Optional maximum K allowed (applies to all forces)")
    p_opt2.add_argument("--em-K-min", type=float, default=None, help="Optional EM K minimum")
    p_opt2.add_argument("--em-K-max", type=float, default=None, help="Optional EM K maximum")
    p_opt2.add_argument("--strong-K-min", type=float, default=None, help="Optional strong K minimum")
    p_opt2.add_argument("--strong-K-max", type=float, default=None, help="Optional strong K maximum")
    p_opt2.add_argument("--weak-K-min", type=float, default=None, help="Optional weak K minimum")
    p_opt2.add_argument("--weak-K-max", type=float, default=None, help="Optional weak K maximum")
    p_opt2.add_argument("--gravity-K-min", type=float, default=None, help="Optional gravity K minimum")
    p_opt2.add_argument("--gravity-K-max", type=float, default=None, help="Optional gravity K maximum")
    p_opt2.set_defaults(func=cmd_pair_forces_option2)

    return p


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())

