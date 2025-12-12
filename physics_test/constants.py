from __future__ import annotations

import math

"""
Minimal constants module.

We avoid requiring SciPy so this repo runs out-of-the-box on a plain Python install.
Values are CODATA 2018 (commonly used; exact constants are now fixed by SI redefinition
for kB and h).
"""

# Exact by SI definition
BOLTZMANN = 1.380649e-23  # J/K
PLANCK = 6.62607015e-34  # J*s
HBAR = PLANCK / (2.0 * math.pi)  # J*s (derived exact)

# Exact by SI definition
SPEED_OF_LIGHT = 299_792_458.0  # m/s

# Dimensionless; CODATA 2018 value (commonly cited)
FINE_STRUCTURE = 7.2973525693e-3  # α

# Newtonian gravitational constant (measured, CODATA 2018-ish)
G_NEWTON = 6.67430e-11  # m^3 / (kg s^2)

# Masses (CODATA 2018-ish)
MASS_ELECTRON = 9.1093837015e-31  # kg
MASS_PROTON = 1.67262192369e-27  # kg

# Useful derived scale (dimensionless gravity coupling uses a chosen mass scale).
# Planck mass: sqrt(hbar*c/G) (approximately 2.176e-8 kg)
MASS_PLANCK = math.sqrt(HBAR * SPEED_OF_LIGHT / G_NEWTON)

PI = math.pi


