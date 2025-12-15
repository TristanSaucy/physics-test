# Active context (physics-test)

This file is meant to help you resume work in a fresh chat / new agent instance.

## Core idea being explored

We’re exploring a framework that links:

- a **dimensionless** “topological gauge” number:

$$
G = \frac{C}{\phi^m}
$$

to a **frequency**:

$$
F_0 = \phi^m\,\frac{k_B K}{h}
$$

with integer $m$ interpreted as discrete harmonic “steps.”

Key invariant:

$$
G\,F_0 = C\,\frac{k_B K}{h}
$$

Inverse (used for Option 2):

$$
K = \frac{F_0 h}{k_B\,\phi^m}
$$

## Important modeling choices made so far

- **m is integer** everywhere (harmonic steps).
- **G is dimensionless**.
- **Targets** are dimensionless couplings (EM/strong/weak + gravity couplings via $\alpha_G$).
- **K interpretation**: we’ve treated K as a scale parameter; for gravity it can be literal (CMB), for micro forces it often behaves better as an “energy-scale temperature.”
- **Thermal anchor**: the model’s natural baseline frequency at $m=0$ is $k_B T/h$. At $T=310\,\mathrm{K}$, $k_B T/h \approx 6.46\,\mathrm{THz}$; at $m=1$, $F_0 \approx 10.45\,\mathrm{THz}$.

## “Non-arbitrary C” constraint strategy

We explored three tiers:

1) **Fixed C=360 only** (strictest)
   - Strong positive: $360/\phi^2$ is close to $1/\alpha$.
   - But 360 alone cannot match strong/weak/gravity couplings within 5%.

2) **Small discrete “topological families”** (degrees/432/fibonacci/etc.) and octave-scaled extensions (powers of two).

3) **Gauge-derived C candidates** (most promising non-arbitrary approach):
   - Generate C values from gauge invariants of U(1), SU(2), SU(3) using base 360.
   - This yields candidates like: 360, 180, 120, 60, 45, 15.

## Key empirical findings (high signal)

- **EM (strict)**: $C=360, m=2$ matches $1/\alpha$ within ~0.34%.
- **Strong/weak (strict inverse targets)**:
  - `1/alpha_s_1loop_from_mZ(mH)` best strict hit: $C=60, m=4$ (≈ −1.27%)
  - `1/alpha2(alpha(mZ),sin2_on_shell)` best strict hit: $C=120, m=3$ (≈ −0.73%)
  - `1/alpha1_GUT(alpha(mZ),sin2)` best strict hit: $C=60, m=0$ (≈ +1.66%)
- **Strong running “within-band” (new)**:
  - Step-only predictive misses at nearby scales (e.g. $m_W$, $m_t$, 1 TeV) suggest \(m\) is a *coarse band index*, not the whole story for running.
  - The strong-only test `oos-predictive-rg` (fit the lattice anchor once, then RG-run to other scales without re-fitting \(m\)) yields **passes at 2%** across the v2/v3/v4 strong-running OOS keys (typical errors \(\sim\)1–2%).
- **EM running “within-band” (new)**:
  - The EM OOS miss `1/alpha → 1/alpha(mZ)` is naturally explained by vacuum polarization (running of \(\alpha(Q)\)).
  - `oos-predictive-rg --force em` supports both a **PDG-style** Δα(mZ²) runner and a simple **1-loop threshold** runner; both yield a **pass at 2%** (typical error \(\sim\)1%).
- **Weak running “within-band” (new)**:
  - Step-only predictive cannot follow smooth EW running to high scales (fails at 1 TeV / 10 TeV in suite `v2`).
  - `oos-predictive-rg --suite v2 --force weak` uses deterministic **SM 1-loop** running for \(1/\alpha_2(Q)\) from the on-shell-defined \(\alpha_2(m_Z)\) anchor and yields **passes at 2%** across `mW`, `mH`, `1TeV`, `10TeV`.
- **Gravity band constraint (CMB K)**:
  - If gravity uses $K=2.725$ K and you require gravity-wave frequencies in a GW band, the frequency equation implies a narrow integer-$m$ window per band (CMB/PTA/LISA/LIGO).
  - Ordinary-matter inverse targets `1/alpha_G(p)` and `1/alpha_G(e)` fit at much more negative $m$ and therefore do **not** land in GW detector bands under CMB $K$.
  - Frozen GW-band “gravity types” (`1/alpha_G(GW_*)`) reconcile coupling-fit + band window under strict gauge-$C$; the end-to-end per-band results are tabulated in `gravity_types_report.md` and in `paper.md` §5.4.2.

## Commands you’ll likely run next

List targets / presets:

- `python -m physics_test.cli list-targets`
- `python -m physics_test.cli list-frequency-presets`
- `python -m physics_test.cli list-gravity-bands`
- See `F0_anchors.md` for the frozen Option‑2 frequency anchor menu (primary + cross-check presets).

Check the EM example:

- `python -m physics_test.cli check-example`

Gauge-derived C exploration:

- `python -m physics_test.cli list-gauge-Cs`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha"`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha_s_1loop_from_mZ(mH)"`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha2(alpha(mZ),sin2_on_shell)"`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha1_GUT(alpha(mZ),sin2)"`

Full pairing under gauge-derived C:

- `python -m physics_test.cli pair-forces-gaugeCs --gravity-band any`

RG-within-band predictive tests (no re-fitting \(m\) per target):

- `python -m physics_test.cli oos-predictive-rg --suite v1 --max-rel-err 0.02`
- `python -m physics_test.cli oos-predictive-rg --suite v2 --max-rel-err 0.02`

Strict all-forces per GW band (frozen Option‑2 presets; see `F0_anchors.md`):

- `python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band cmb  --gravity-targets "1/alpha_G(GW_CMB)"  --max-hits 10 --max-results 5`
- `python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band pta  --gravity-targets "1/alpha_G(GW_PTA)"  --max-hits 10 --max-results 5`
- `python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band lisa --gravity-targets "1/alpha_G(GW_LISA)" --max-hits 10 --max-results 5`
- `python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band ligo --gravity-targets "1/alpha_G(GW_LIGO)" --max-hits 10 --max-results 5`

Recommended strict Option‑2 anchor presets (frozen):

- `--em-preset em-lyman-alpha`
- `--strong-preset strong-QCD-200MeV`
- `--weak-preset weak-W-80.379GeV`

Quantum-gravity mass sweep (CMB band):

- `python -m physics_test.cli sweep-quantum-gravity --gravity-band cmb --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12`

## Open questions / next hypotheses

- Are the gauge-derived C choices (360/120/60/45/15) actually the “right” non-arbitrary set, or do we need more invariants (rank, Casimirs, etc.)?
- Gravity orientation is now frozen to inverse coupling targets (`1/alpha_G(mass)`).
- Ordinary-matter gravity is frozen to the **proton** target `1/alpha_G(p)`, with **electron** `1/alpha_G(e)` as a required cross-check (not a free alternative).
- GW-band gravity “types” are now frozen to these mass-anchor targets:
  - `1/alpha_G(GW_CMB)`, `1/alpha_G(GW_PTA)`, `1/alpha_G(GW_LISA)`, `1/alpha_G(GW_LIGO)`
- Planck/quantum gravity type is frozen to `1/alpha_G(mP)` (~1).
- EM is now frozen to `1/alpha` (low-energy) for strict runs (minimal scheme baggage; matches the 360 anchor cleanly). `1/alpha(mZ)` remains exploratory.
- Next: tighten tolerance (5% → 2% → 1%), add a simple score/ranking, and expand gauge-derived $C$ constructions only if the allowed set stays short.
