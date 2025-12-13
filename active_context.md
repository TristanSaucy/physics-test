## Active context (physics-test)

This file is meant to help you resume work in a fresh chat / new agent instance.

### Core idea being explored

We’re exploring a framework that links:

- a **dimensionless** “topological gauge” number:

\[
G = \\frac{C}{\\phi^m}
\]

to a **frequency**:

\[
F_0 = \\phi^m\\,\\frac{k_B K}{h}
\]

with integer \(m\\) interpreted as discrete harmonic “steps.”

Key invariant:

\[
G\\,F_0 = C\\,\\frac{k_B K}{h}
\]

Inverse (used for Option 2):

\[
K = \\frac{F_0 h}{k_B\\,\\phi^m}
\]

### Important modeling choices made so far

- **m is integer** everywhere (harmonic steps).
- **G is dimensionless**.
- **Targets** are dimensionless couplings (EM/strong/weak + gravity couplings via \\(\\alpha_G\\)).
- **K interpretation**: we’ve treated K as a scale parameter; for gravity it can be literal (CMB), for micro forces it often behaves better as an “energy-scale temperature.”
- **Thermal anchor**: the model’s natural baseline frequency at \(m=0\) is \(k_B T/h\). At \(T=310K\), \(k_B T/h \\approx 6.46\\,\\text{THz}\\); at \(m=1\), \(F_0 \\approx 10.45\\,\\text{THz}\\).

### “Non-arbitrary C” constraint strategy

We explored three tiers:

1) **Fixed C=360 only** (strictest)
   - Strong positive: \(360/\\phi^2\\) is close to \(1/\\alpha\\).
   - But 360 alone cannot match strong/weak/gravity couplings within 5%.

2) **Small discrete “topological families”** (degrees/432/fibonacci/etc.) and octave-scaled extensions (powers of two).

3) **Gauge-derived C candidates** (most promising non-arbitrary approach):
   - Generate C values from gauge invariants of U(1), SU(2), SU(3) using base 360.
   - This yields candidates like: 360, 180, 120, 60, 45, 15.

### Key empirical findings (high signal)

- **EM**: \(C=360, m=2\\) matches \(1/\\alpha\\) within ~0.34%.
- **Gauge-derived C scan** found <5% matches for:
  - strong coupling benchmark `alpha_s(mZ)` using C=60 (m=13) or C=15 (m=10)
  - weak coupling proxy `alpha_w(mZ)` using C=120 (m=17) or C=45 (m=15)
- **Gravity band constraint (CMB K)**:
  - If gravity uses K=CMB and you require gravity-wave frequencies in the **CMB/primordial band** (~1e-18..1e-16 Hz), the frequency equation implies gravity m roughly in **-137..-129**.
  - With gravity target `1/alpha_G(p)` and gauge-only C candidates, gravity coupling fit pushes to much more negative m (e.g. ~-175), so it doesn’t land in the CMB GW band.
  - Introducing “gravity type” targets `1/alpha_G(scale)` where scale is a **TeV mass** can reconcile both constraints.

### Commands you’ll likely run next

List targets / presets:

- `python -m physics_test.cli list-targets`
- `python -m physics_test.cli list-frequency-presets`
- `python -m physics_test.cli list-gravity-bands`

Check the EM example:

- `python -m physics_test.cli check-example`

Gauge-derived C exploration:

- `python -m physics_test.cli list-gauge-Cs`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha"`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha_s(mZ)"`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha_w(mZ)"`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha1_GUT(alpha(mZ),sin2)"`

Full pairing under gauge-derived C:

- `python -m physics_test.cli pair-forces-gaugeCs --gravity-band any`
- `python -m physics_test.cli pair-forces-gaugeCs --gravity-band cmb --gravity-targets "1/alpha_G(p),1/alpha_G(40TeV),1/alpha_G(100TeV)" ...`

Quantum-gravity mass sweep (CMB band):

- `python -m physics_test.cli sweep-quantum-gravity --gravity-band cmb --scale-min-GeV 1e3 --scale-max-GeV 1e6 --n-scales 121 --top 30`

### Open questions / next hypotheses

- Are the gauge-derived C choices (360/120/60/45/15) actually the “right” non-arbitrary set, or do we need more invariants (rank, Casimirs, etc.)?
- Gravity orientation is now frozen to inverse coupling targets (`1/alpha_G(mass)`).
- Ordinary-matter gravity is now frozen to the **proton** target `1/alpha_G(p)`, with **electron** `1/alpha_G(e)` as a required cross-check (not a free alternative).
- Remaining gravity open question: which mass anchors correspond to GW-band / primordial / quantum-gravity “types.”
- EM is now frozen to `1/alpha` (low-energy) for strict runs (minimal scheme baggage; matches the 360 anchor cleanly). `1/alpha(mZ)` remains exploratory.
- How to choose phenomenon frequencies F0 (Option 2) without introducing arbitrariness.


