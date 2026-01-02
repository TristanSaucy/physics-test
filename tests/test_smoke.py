from __future__ import annotations


def test_registry_loads_alpha0() -> None:
    # Basic smoke test: registry loads and contains the canonical low-energy EM input.
    from physics_test.target_registry import get_measurement

    m = get_measurement("alpha0")
    assert m.value > 0
    assert m.sigma is not None
    assert m.scheme
    assert m.citation


def test_known_targets_include_registry_driven_sin2() -> None:
    # Smoke test: registry-driven targets (tgt_*) are surfaced via known_targets().
    from physics_test.targets import known_targets

    names = {t.name for t in known_targets()}

    # Core canonical targets
    assert "1/alpha" in names
    assert "sin2thetaW(on-shell)" in names

    # External low-Q targets added via data/targets.json as tgt_sin2thetaW(...)
    assert "sin2thetaW(Qweak)" in names
    assert "sin2thetaW(E158)" in names
    assert "sin2thetaW(nue_lowE)" in names
    assert "sin2thetaW(CHARMII)" in names
    assert "sin2thetaW(NuTeV)" in names
    assert "sin2thetaW(NuTeV_eq10)" in names
    assert "sin2thetaW(CsAPV)" in names
    assert "sin2thetaW(CsAPV_PDG)" in names
    assert "sin2thetaW(eDIS)" in names


def test_ew_independent_v1_suite_is_frozen_and_resolvable() -> None:
    from physics_test.oos import ew_sin2_suites
    from physics_test.targets import known_targets

    by_name = {t.name: t for t in known_targets()}
    suites = ew_sin2_suites()
    assert "ew-independent-v1" in suites
    assert "ew-exploratory-v1" in suites
    assert "ew-dis-exploratory-v1" in suites
    assert "ew-zpole-exploratory-v1" in suites

    # Ensure exploratory-only points do not leak into the CI-grade ew-independent suites.
    indep_v1 = {ot.key for ot in suites["ew-independent-v1"]}
    indep_v2 = {ot.key for ot in suites["ew-independent-v2"]}
    indep_v3 = {ot.key for ot in suites["ew-independent-v3"]}
    assert "sin2thetaW(CsAPV_PDG)" not in indep_v1
    assert "sin2thetaW(CsAPV_PDG)" not in indep_v2
    assert "sin2thetaW(CsAPV_PDG)" not in indep_v3
    assert "sin2thetaW(CHARMII)" not in indep_v1
    assert "sin2thetaW(CHARMII)" not in indep_v2
    assert "sin2thetaW(CHARMII)" not in indep_v3
    assert "sin2thetaW(NuTeV)" not in indep_v1
    assert "sin2thetaW(NuTeV)" not in indep_v2
    assert "sin2thetaW(NuTeV)" not in indep_v3
    assert "sin2thetaW(NuTeV_eq10)" not in indep_v1
    assert "sin2thetaW(NuTeV_eq10)" not in indep_v2
    assert "sin2thetaW(NuTeV_eq10)" not in indep_v3
    assert "sin2thetaW(LEP_SLC_Zpole)" not in indep_v1
    assert "sin2thetaW(LEP_SLC_Zpole)" not in indep_v2
    assert "sin2thetaW(LEP_SLC_Zpole)" not in indep_v3
    assert "sin2thetaW(Tevatron_Zpole)" not in indep_v1
    assert "sin2thetaW(Tevatron_Zpole)" not in indep_v2
    assert "sin2thetaW(Tevatron_Zpole)" not in indep_v3

    # Exploratory suite includes the extra points explicitly.
    exploratory = {ot.key for ot in suites["ew-exploratory-v1"]}
    assert "sin2thetaW(CsAPV_PDG)" in exploratory
    assert "sin2thetaW(CHARMII)" in exploratory
    assert "sin2thetaW(NuTeV)" not in exploratory
    assert "sin2thetaW(NuTeV_eq10)" not in exploratory
    assert "sin2thetaW(LEP_SLC_Zpole)" not in exploratory
    assert "sin2thetaW(Tevatron_Zpole)" not in exploratory

    ew_dis = {ot.key for ot in suites["ew-dis-exploratory-v1"]}
    assert "sin2thetaW(NuTeV)" in ew_dis
    assert "sin2thetaW(NuTeV_eq10)" in ew_dis

    ew_zpole = {ot.key for ot in suites["ew-zpole-exploratory-v1"]}
    assert "sin2thetaW(LEP_SLC_Zpole)" in ew_zpole
    assert "sin2thetaW(Tevatron_Zpole)" in ew_zpole

    for ot in suites["ew-independent-v1"]:
        assert ot.key in by_name
        t = by_name[ot.key]
        assert t.scheme is not None
        assert str(t.scheme).startswith("sin2thetaW_eff:")

    assert "ew-independent-v2" in suites
    for ot in suites["ew-independent-v2"]:
        assert ot.key in by_name
        t = by_name[ot.key]
        assert t.scheme is not None
        assert str(t.scheme).startswith("sin2thetaW_eff:")

    assert "ew-independent-v3" in suites
    for ot in suites["ew-independent-v3"]:
        assert ot.key in by_name
        t = by_name[ot.key]
        assert t.scheme is not None
        assert str(t.scheme).startswith("sin2thetaW_eff:")

    for ot in suites["ew-exploratory-v1"]:
        assert ot.key in by_name
        t = by_name[ot.key]
        assert t.scheme is not None
        assert str(t.scheme).startswith("sin2thetaW_eff:")

    for ot in suites["ew-dis-exploratory-v1"]:
        assert ot.key in by_name
        t = by_name[ot.key]
        assert t.scheme is not None
        assert str(t.scheme).startswith("sin2thetaW_on_shell:")

    for ot in suites["ew-zpole-exploratory-v1"]:
        assert ot.key in by_name
        t = by_name[ot.key]
        assert t.scheme is not None
        assert str(t.scheme).startswith("sin2thetaW_eff_lept:")


def test_ew_independent_v1_gate_exits_nonzero_on_fail() -> None:
    # Ensure the "one-command gate" semantics are real: suite ew-independent-v1 should
    # return non-zero when forced to fail.
    import subprocess
    import sys

    cmd = [
        sys.executable,
        "-m",
        "physics_test.cli",
        "oos-ew-sin2",
        "--suite",
        "ew-independent-v1",
        "--model",
        "sm",
        "--method",
        "gammaZ_1loop",
        "--z-max",
        "0.01",  # force failure
        "--m-min",
        "-6",
        "--m-max",
        "6",
        "--m-step",
        "1",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert res.returncode != 0


def test_ew_exploratory_v1_runs_and_is_nongating_by_default() -> None:
    # Ensure exploratory suite runs and does NOT act as a CI gate unless --ci is passed.
    import subprocess
    import sys

    cmd = [
        sys.executable,
        "-m",
        "physics_test.cli",
        "oos-ew-sin2",
        "--suite",
        "ew-exploratory-v1",
        "--model",
        "sm",
        "--method",
        "gammaZ_1loop",
        "--m-min",
        "-6",
        "--m-max",
        "6",
        "--m-step",
        "1",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert res.returncode == 0


def test_ew_dis_exploratory_v1_runs_and_is_nongating_by_default() -> None:
    # Ensure NuTeV DIS suite runs and does NOT act as a CI gate unless --ci is passed.
    import subprocess
    import sys

    cmd = [
        sys.executable,
        "-m",
        "physics_test.cli",
        "oos-ew-sin2",
        "--suite",
        "ew-dis-exploratory-v1",
        "--model",
        "sm",
        "--method",
        "gammaZ_1loop",
        "--m-min",
        "-6",
        "--m-max",
        "6",
        "--m-step",
        "1",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert res.returncode == 0
    assert "sin2thetaW(NuTeV)" in res.stdout


def test_ew_zpole_exploratory_v1_runs_and_is_nongating_by_default() -> None:
    # Ensure Z-pole effective-angle suite runs and does NOT act as a CI gate unless --ci is passed.
    import subprocess
    import sys

    cmd = [
        sys.executable,
        "-m",
        "physics_test.cli",
        "oos-ew-sin2",
        "--suite",
        "ew-zpole-exploratory-v1",
        "--model",
        "sm",
        "--method",
        "gammaZ_1loop",
        "--m-min",
        "-6",
        "--m-max",
        "6",
        "--m-step",
        "1",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert res.returncode == 0
    assert "sin2thetaW(LEP_SLC_Zpole)" in res.stdout


def test_registry_all_scheme_prefix_filters_targets() -> None:
    # Scheme isolation: registry-all should treat --scheme-prefix as a filter so that
    # eff-only and on-shell-only reports do not mix.
    import subprocess
    import sys

    cmd_eff = [
        sys.executable,
        "-m",
        "physics_test.cli",
        "oos-ew-sin2",
        "--suite",
        "registry-all",
        "--model",
        "sm",
        "--method",
        "gammaZ_1loop",
        "--scheme-prefix",
        "sin2thetaW_eff:",
        "--z-max",
        "2",
        "--m-min",
        "-6",
        "--m-max",
        "6",
        "--m-step",
        "1",
    ]
    res_eff = subprocess.run(cmd_eff, capture_output=True, text=True, timeout=60)
    assert res_eff.returncode == 0
    assert "sin2thetaW(Qweak)" in res_eff.stdout
    assert "sin2thetaW(NuTeV)" not in res_eff.stdout

    cmd_os = [
        sys.executable,
        "-m",
        "physics_test.cli",
        "oos-ew-sin2",
        "--suite",
        "registry-all",
        "--model",
        "sm",
        "--method",
        "gammaZ_1loop",
        "--scheme-prefix",
        "sin2thetaW_on_shell:",
        "--z-max",
        "2",
        "--m-min",
        "-6",
        "--m-max",
        "6",
        "--m-step",
        "1",
    ]
    res_os = subprocess.run(cmd_os, capture_output=True, text=True, timeout=60)
    assert res_os.returncode == 0
    assert "sin2thetaW(NuTeV)" in res_os.stdout
    assert "sin2thetaW(Qweak)" not in res_os.stdout


def test_base_vs_alt_bases_command_smoke() -> None:
    # Lightweight smoke test: ensure the pre-registered base scan command runs.
    import subprocess
    import sys

    cmd = [
        sys.executable,
        "-m",
        "physics_test.cli",
        "base-vs-alt-bases",
        "--bases",
        "360,420",
        "--tol",
        "0.05",
        "--max-nCs",
        "10",
        "--m-min",
        "-12",
        "--m-max",
        "12",
        "--m-step",
        "1",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert res.returncode == 0
    assert "Base-vs-alt-bases scan" in res.stdout

