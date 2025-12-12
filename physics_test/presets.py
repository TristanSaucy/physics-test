from __future__ import annotations

from dataclasses import dataclass

from physics_test.units import energy_J_from_GeV, energy_J_from_MeV, energy_J_from_eV, frequency_Hz_from_energy_J
from physics_test import constants


@dataclass(frozen=True)
class FrequencyPreset:
    key: str
    F0_hz: float
    note: str = ""


def em_frequency_presets() -> list[FrequencyPreset]:
    """
    EM frequency presets (rough bands or common reference lines).
    These are not exact; they're just convenient anchors for Option-2 exploration.
    """
    return [
        FrequencyPreset("em-radio-1GHz", 1e9, "Radio ~1 GHz"),
        FrequencyPreset("em-microwave-100GHz", 1e11, "Microwave ~100 GHz"),
        FrequencyPreset("em-visible-500THz", 5e14, "Visible light ~500 THz (green-ish)"),
        FrequencyPreset("em-uv-1PHz", 1e15, "UV/near-UV ~1 PHz"),
        FrequencyPreset("em-lyman-alpha", 2.466e15, "Hydrogen Lyman-alpha line (~121.6 nm)"),
        FrequencyPreset("em-xray-1EHz", 1e18, "X-ray ~1 EHz (order-of-magnitude)"),
    ]


def particle_proxy_presets() -> list[FrequencyPreset]:
    """
    Particle-energy proxies expressed as frequencies f=E/h.
    This is a common way to attach a frequency scale to a particle mass/energy.
    """
    # These use simple energies, not precise particle masses.
    return [
        FrequencyPreset(
            "strong-timescale-1e-23s",
            1e23,
            "Phenomenon proxy: typical strong-interaction timescale ~1e-23 s -> ~1e23 Hz",
        ),
        FrequencyPreset(
            "weak-muon-decay",
            1.0 / 2.1969811e-6,
            "Phenomenon proxy: muon lifetime ~2.2 microseconds -> ~4.55e5 Hz",
        ),
        FrequencyPreset(
            "strong-QCD-200MeV",
            frequency_Hz_from_energy_J(energy_J_from_MeV(200.0)),
            "Proxy strong scale: ~200 MeV -> f=E/h",
        ),
        FrequencyPreset(
            "strong-proton-938MeV",
            frequency_Hz_from_energy_J(energy_J_from_MeV(938.0)),
            "Proxy strong scale: proton mass-energy -> f=E/h",
        ),
        FrequencyPreset(
            "weak-W-80.379GeV",
            frequency_Hz_from_energy_J(energy_J_from_GeV(80.379)),
            "Proxy weak scale: W mass-energy -> f=E/h",
        ),
        FrequencyPreset(
            "weak-Z-91.1876GeV",
            frequency_Hz_from_energy_J(energy_J_from_GeV(91.1876)),
            "Proxy weak scale: Z mass-energy -> f=E/h",
        ),
        FrequencyPreset(
            "em-hydrogen-13.6eV",
            frequency_Hz_from_energy_J(energy_J_from_eV(13.6)),
            "Proxy EM atomic scale: 13.6 eV -> f=E/h",
        ),
    ]


def thermal_presets() -> list[FrequencyPreset]:
    """
    Thermal frequency anchors using the universal scale kB*T/h.

    Note: in the model, F0 = phi^m * (kB*K/h), so these are most naturally interpreted
    as the m=0 baseline for a given temperature.
    """

    def _kbt_over_h(T_K: float) -> float:
        return (constants.BOLTZMANN * float(T_K)) / constants.PLANCK

    return [
        FrequencyPreset(
            "thermal-CMB-2.725K-kBT_over_h",
            _kbt_over_h(2.725),
            "Thermal baseline: kB*T/h at CMB temperature (m=0 baseline)",
        ),
        FrequencyPreset(
            "thermal-room-300K-kBT_over_h",
            _kbt_over_h(300.0),
            "Thermal baseline: kB*T/h at ~room temperature (m=0 baseline)",
        ),
        FrequencyPreset(
            "thermal-body-310K-kBT_over_h",
            _kbt_over_h(310.0),
            "Thermal baseline: kB*T/h at ~human body temperature (m=0 baseline)",
        ),
    ]


def get_preset(key: str) -> FrequencyPreset:
    for p in em_frequency_presets() + particle_proxy_presets() + thermal_presets():
        if p.key == key:
            return p
    raise KeyError(f"Unknown preset: {key}")


