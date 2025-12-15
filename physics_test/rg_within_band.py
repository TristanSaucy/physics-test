from __future__ import annotations

import re
from dataclasses import dataclass

from physics_test.forces import (
    alpha_s_run_1loop_from_ref,
    alpha_s_run_1loop_from_ref_nf_switch,
    alpha_s_run_2loop_from_ref,
    alpha_s_run_2loop_from_ref_nf_switch,
)
from physics_test.scales import M_T_GEV, scale_GeV


@dataclass(frozen=True)
class QCDRunSpec:
    """
    Deterministic QCD running configuration for "within-band" predictions.

    - loops: 1 or 2
    - nf_mode:
        - "const": use a fixed n_f
        - "nf56": switch from n_f=5 to n_f=6 at Q_switch (default ~mt)
    """

    loops: int = 1
    nf_mode: str = "const"  # "const" | "nf56"
    n_f: int = 5
    Q_switch_GeV: float = float(M_T_GEV)
    n_f_below: int = 5
    n_f_above: int = 6
    steps_per_unit_log: int = 500  # only used for 2-loop integrator


def extract_paren_label(key: str) -> str | None:
    """
    Extract the final "(...)" label from a target key.

    Examples:
      - "1/alpha_s_1loop_from_mZ(mH)" -> "mH"
      - "1/alpha_s(mZ)" -> "mZ"
      - "1/alpha" -> None
    """

    m = re.search(r"\(([^()]*)\)\s*$", key)
    return m.group(1) if m else None


def scale_GeV_from_target_key(key: str) -> float:
    lab = extract_paren_label(key)
    if lab is None:
        raise KeyError(f"Target key has no '(...)' scale label: {key!r}")
    return scale_GeV(lab)


def qcd_run_spec_from_key(key: str) -> QCDRunSpec:
    """
    Infer a deterministic QCD running spec from a target key string.

    This is used to avoid adding free knobs: if a target is labeled "2loop" or "nf56",
    we use the corresponding deterministic running prescription.
    """

    loops = 2 if "2loop" in key else 1
    nf_mode = "nf56" if "nf56" in key else "const"
    return QCDRunSpec(loops=loops, nf_mode=nf_mode)


def alpha_s_from_ref(Q_GeV: float, *, alpha_s_Q0: float, Q0_GeV: float, spec: QCDRunSpec) -> float:
    """
    Run alpha_s(Q) deterministically from a reference alpha_s(Q0).
    """

    if spec.loops == 1:
        if spec.nf_mode == "nf56":
            return alpha_s_run_1loop_from_ref_nf_switch(
                Q_GeV,
                alpha_s_Q0=alpha_s_Q0,
                Q0_GeV=Q0_GeV,
                Q_switch_GeV=spec.Q_switch_GeV,
                n_f_below=spec.n_f_below,
                n_f_above=spec.n_f_above,
            )
        return alpha_s_run_1loop_from_ref(Q_GeV, alpha_s_Q0=alpha_s_Q0, Q0_GeV=Q0_GeV, n_f=spec.n_f)

    if spec.loops == 2:
        if spec.nf_mode == "nf56":
            return alpha_s_run_2loop_from_ref_nf_switch(
                Q_GeV,
                alpha_s_Q0=alpha_s_Q0,
                Q0_GeV=Q0_GeV,
                Q_switch_GeV=spec.Q_switch_GeV,
                n_f_below=spec.n_f_below,
                n_f_above=spec.n_f_above,
                steps_per_unit_log=spec.steps_per_unit_log,
            )
        return alpha_s_run_2loop_from_ref(
            Q_GeV,
            alpha_s_Q0=alpha_s_Q0,
            Q0_GeV=Q0_GeV,
            n_f=spec.n_f,
            steps_per_unit_log=spec.steps_per_unit_log,
        )

    raise ValueError("spec.loops must be 1 or 2")


