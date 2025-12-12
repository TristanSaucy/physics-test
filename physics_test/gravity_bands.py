from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GravityBand:
    name: str
    f_min_hz: float
    f_max_hz: float
    note: str = ""


def bands() -> list[GravityBand]:
    """
    Rough observational frequency bands for gravitational waves.
    These are order-of-magnitude ranges, not sharp cutoffs.
    """
    return [
        GravityBand("ligo", 10.0, 1000.0, "Ground-based interferometers (LIGO/Virgo/KAGRA)"),
        GravityBand("lisa", 1e-4, 1e-1, "Space-based interferometer (LISA)"),
        GravityBand("pta", 1e-9, 1e-7, "Pulsar Timing Arrays (nHz band)"),
        GravityBand("cmb", 1e-18, 1e-16, "Primordial/horizon-scale (CMB-sensitive)"),
    ]


