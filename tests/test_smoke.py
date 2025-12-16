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
    assert "sin2thetaW(CsAPV)" in names
    assert "sin2thetaW(CsAPV_PDG)" in names
    assert "sin2thetaW(eDIS)" in names


def test_ew_independent_v1_suite_is_frozen_and_resolvable() -> None:
    from physics_test.oos import ew_sin2_suites
    from physics_test.targets import known_targets

    by_name = {t.name: t for t in known_targets()}
    suites = ew_sin2_suites()
    assert "ew-independent-v1" in suites
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

