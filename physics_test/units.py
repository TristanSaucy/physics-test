from __future__ import annotations

from physics_test import constants


EV_TO_J = 1.602176634e-19  # exact (J/eV)


def energy_J_from_eV(eV: float) -> float:
    return float(eV) * EV_TO_J


def energy_J_from_MeV(MeV: float) -> float:
    return energy_J_from_eV(float(MeV) * 1e6)


def energy_J_from_GeV(GeV: float) -> float:
    return energy_J_from_eV(float(GeV) * 1e9)


def temperature_K_from_energy_J(E_J: float) -> float:
    """K = E/kB."""
    return float(E_J) / constants.BOLTZMANN


def frequency_Hz_from_energy_J(E_J: float) -> float:
    """f = E/h."""
    return float(E_J) / constants.PLANCK


def mass_kg_from_energy_J(E_J: float) -> float:
    """m = E/c^2."""
    return float(E_J) / (constants.SPEED_OF_LIGHT**2)


def mass_kg_from_eV(eV: float) -> float:
    return mass_kg_from_energy_J(energy_J_from_eV(eV))


def mass_kg_from_MeV(MeV: float) -> float:
    return mass_kg_from_energy_J(energy_J_from_MeV(MeV))


def mass_kg_from_GeV(GeV: float) -> float:
    return mass_kg_from_energy_J(energy_J_from_GeV(GeV))


