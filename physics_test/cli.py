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
from physics_test.oos import oos_suites, predictive_force_suites, resolve_oos_targets
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
from physics_test.gut import MSSM_1LOOP, SM_1LOOP, converge_score, find_best_convergence, run_alpha_inv
from physics_test.normalization import normalization_factor_for_force, normalization_families
from physics_test.steps import step_from_targets
from physics_test.rg_scales import lambda_qcd_from_alpha_s
from physics_test.oos_rg import rg_suites
from physics_test.rg_within_band import QCDRunSpec, alpha_s_from_ref, qcd_run_spec_from_key, scale_GeV_from_target_key
from physics_test.qed_running import alpha_inv_mZ_from_delta_alpha, qed_run_alpha_inv_1loop_from_ref
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
        inv2 = run_alpha_inv(inv_a2_0, Q0, Q, betas.b2)
        inv1 = run_alpha_inv(inv_a1_0, Q0, Q, betas.b1)
        invY = (5.0 / 3.0) * float(inv1)
        sin2_pred = float(inv2) / (float(inv2) + invY) if (float(inv2) + invY) != 0 else float("nan")

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

    # Collect registry-added external sin2 targets
    ext: list[tuple[str, float, float | None, float]] = []
    for t in all_targets.values():
        if not t.name.startswith("sin2thetaW("):
            continue
        if not str(t.note).startswith("Registry target:"):
            continue
        if t.Q_GeV is None:
            continue
        ext.append((t.name, float(t.value), t.sigma, float(t.Q_GeV)))

    if not ext:
        raise SystemExit(
            "No external sin2thetaW(Q) targets found. Add registry keys like "
            "'tgt_sin2thetaW(Qweak)' with fields value/sigma/Q_GeV/scheme/citation "
            "(you can use PHYSICS_TEST_TARGET_REGISTRY to point to a local registry file)."
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

    print("OOS: external sin2thetaW(Q) targets (registry-driven)")
    print(f"model = {betas.name}")
    print(f"tol(|rel_err|) = {args.max_rel_err}")
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
    for name, v, s, Q in sorted(ext, key=lambda x: x[3]):
        inv2 = run_alpha_inv(inv_a2_0, Q0, Q, betas.b2)
        inv1 = run_alpha_inv(inv_a1_0, Q0, Q, betas.b1)
        invY = (5.0 / 3.0) * float(inv1)
        pred = float(inv2) / (float(inv2) + invY) if (float(inv2) + invY) != 0 else float("nan")

        rel_err = (pred - v) / v if v != 0 else float("nan")
        ok = abs(rel_err) <= float(args.max_rel_err)
        status = "PASS" if ok else "FAIL"
        if ok:
            n_pass += 1
        n += 1

        if s is not None and s > 0:
            z = (pred - v) / float(s)
            n_sigma += 1
            chi2 += float(z * z)
            z_s = f"  z={z:.3g}"
        else:
            z_s = ""

        print(f"[{status}] {name:28s} Q={Q:g} GeV  meas={v:.12g}  pred={pred:.12g}  rel_err={rel_err:.3e}{z_s}")

    if n_sigma:
        print(f"\nSummary: {n_pass}/{n} PASS at tol={args.max_rel_err}  (sigma-annotated: n={n_sigma}, chi2={chi2:.6g})")
    else:
        print(f"\nSummary: {n_pass}/{n} PASS at tol={args.max_rel_err}")
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
    p_oos_ew_sin2.add_argument("--model", choices=["sm", "mssm"], default="sm", help="1-loop beta set for running (default: sm)")
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

