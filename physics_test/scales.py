from __future__ import annotations

"""
Small, shared numeric scales used across target definitions and RG-style reports.

We keep these in one place to avoid accidental drift between different parts of the codebase.
"""


# Pole-mass / benchmark scales (GeV) used in target keys like "(mW)", "(mZ)", etc.
M_W_GEV: float = 80.379
M_Z_GEV: float = 91.1876
M_H_GEV: float = 125.0
M_T_GEV: float = 172.76


def scale_GeV(label: str) -> float:
    """
    Convert a simple label used in target keys to an energy scale in GeV.

    Supported examples:
      - "mW", "mZ", "mH", "mt"
      - "1TeV", "10TeV", "40TeV", "100TeV"
      - numeric strings like "1000" (interpreted as GeV)
    """

    x = label.strip()
    if not x:
        raise KeyError("Empty scale label")

    # Common particle benchmark labels
    if x == "mW":
        return float(M_W_GEV)
    if x == "mZ":
        return float(M_Z_GEV)
    if x == "mH":
        return float(M_H_GEV)
    if x == "mt":
        return float(M_T_GEV)

    # Simple TeV parsing: "<number>TeV"
    if x.endswith("TeV"):
        num = x[: -len("TeV")]
        return float(num) * 1_000.0

    # Fallback: treat as a numeric GeV value
    return float(x)


