from __future__ import annotations

import json
import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path


@dataclass(frozen=True)
class Measurement:
    """
    A small, auditable measured input (or pinned convention) with uncertainty metadata.

    - value: central value
    - sigma: 1σ uncertainty (None means unknown / not provided)
    - Q_GeV: reference scale if applicable (None if not applicable)
    - scheme: short scheme/context description
    - citation: human-readable citation hint (URL-free; keep bare URLs out of code/docs)
    """

    key: str
    value: float
    sigma: float | None = None
    Q_GeV: float | None = None
    scheme: str = ""
    citation: str = ""


def _default_registry_path() -> Path:
    # Allow overriding the built-in registry for local experiments without editing the repo.
    # Note: cache is per-process; set env var before running the CLI.
    env = os.getenv("PHYSICS_TEST_TARGET_REGISTRY", "").strip()
    if env:
        return Path(env)
    # repo_root / data / targets.json
    return Path(__file__).resolve().parents[1] / "data" / "targets.json"


@lru_cache(maxsize=1)
def load_registry(path: str | None = None) -> dict[str, Measurement]:
    """
    Load the measurement registry.

    By design, this registry is small. If the file is missing, return an empty dict.
    """

    p = Path(path) if path is not None else _default_registry_path()
    if not p.exists():
        return {}
    raw = json.loads(p.read_text(encoding="utf-8"))
    out: dict[str, Measurement] = {}
    for k, v in raw.items():
        if str(k).startswith("_"):
            continue
        if not isinstance(v, dict):
            continue
        out[str(k)] = Measurement(
            key=str(k),
            value=float(v["value"]),
            sigma=None if v.get("sigma", None) is None else float(v["sigma"]),
            Q_GeV=None if v.get("Q_GeV", None) is None else float(v["Q_GeV"]),
            scheme=str(v.get("scheme", "")),
            citation=str(v.get("citation", "")),
        )
    return out


def get_measurement(
    key: str,
    *,
    default_value: float | None = None,
    default_sigma: float | None = None,
    default_Q_GeV: float | None = None,
    default_scheme: str = "",
    default_citation: str = "",
) -> Measurement:
    """
    Fetch a measurement by key; if missing, return a Measurement populated with defaults.
    """

    reg = load_registry()
    if key in reg:
        return reg[key]
    if default_value is None:
        raise KeyError(f"Missing measurement key {key!r} and no default provided")
    return Measurement(
        key=str(key),
        value=float(default_value),
        sigma=default_sigma if default_sigma is None else float(default_sigma),
        Q_GeV=default_Q_GeV if default_Q_GeV is None else float(default_Q_GeV),
        scheme=str(default_scheme),
        citation=str(default_citation),
    )


