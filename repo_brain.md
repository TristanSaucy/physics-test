## Repo “brain” (durable context + how to think about the project)

If you want to keep iterating without re-deriving everything, read this first.

### 1) What we’re trying to test

We’re testing whether a minimal harmonic/topology model can connect:

- **dimensionless couplings** (pure numbers) via

\[
G = \\frac{C}{\\phi^m}
\]

- to **frequency/temperature scaling** via

\[
F_0 = \\phi^m\\,\\frac{k_B K}{h}
\]

with **integer m** as a discrete “trophic/harmonic step.”

The model’s built-in invariant:

\[
G\\,F_0 = C\\,\\frac{k_B K}{h}
\]

### 2) What has been most promising (signals)

#### EM + 360
- \(C=360, m=2\\) gives \(G\\approx 137.5078\\), close to **\(1/\\alpha\\)** (0.34%).
- This is the cleanest “360 is special” result so far.

#### Gauge-derived non-arbitrary C set
We derived candidate \(C\) values from simple gauge invariants (U(1), SU(2), SU(3)) using base=360:
- 360, 180, 120, 60, 45, 15 (and a couple duplicates)

Scanning only those C values, integer m, and 5% tolerance gives plausible fits for:
- EM (`1/alpha`) at C=360, m=2
- Strong (strict inverse benchmark `1/alpha_s(mZ)`) at C=60, m=4
- Weak (strict inverse benchmark `1/alpha_w(mZ)`) at C=120, m=3

This is valuable because it reduces “arbitrary C.”

#### Gravity: why “type” matters
With gravity K fixed to the CMB, the frequency equation forces gravity’s m into narrow windows for each GW band.
Reconciling gravity coupling fits with those m windows typically requires changing the **gravity coupling definition**, e.g. using `1/alpha_G(scale)` at TeV-scale “quantum gravity” mass scales.

### 3) Main risk / how to stay falsifiable

The risk is overfitting: if C is allowed to roam freely or if we keep changing target definitions ad hoc, we can always find matches.

To stay falsifiable:
- Keep C restricted to a **small, justified set** (e.g., gauge-derived).
- Fix target definitions and reference scales (strong/weak run with scale).
- Fix GW band constraints if you claim “gravity type X lives here.”
- Prefer experiments where a constraint *eliminates* most candidates.

### 4) The meaning of K

K can mean:
- literal thermodynamic temperature (works naturally for gravity+CMB), or
- an energy-scale temperature (E = kB*K) for micro physics.

Option 2 (“phenomenon-first”) is useful because it forces K to be implied by m and the chosen F0, rather than being tuned.

### 5) Where to work next (high-value)

1) Decide on a canonical mapping for “force → coupling target”
   - EM is now frozen to `1/alpha` (low-energy).
   - Strong: which inverse target and reference scale are frozen? (current strict benchmark: `1/alpha_s(mZ)`.)
   - Weak: which inverse target and reference scale are frozen? (current strict benchmark: `1/alpha_w(mZ)`.)
   - Gravity: which alpha_G definition corresponds to which band/type?

2) Decide how to anchor “phenomenon frequency” F0 for strong/weak in Option 2 so it’s not arbitrary.

3) Add per-force plausibility windows (K ranges, band ranges) only if they reflect a stated interpretation.

### 6) Commands to reproduce core results

- EM 360 example:
  - `python -m physics_test.cli check-example`

- Gauge-derived C list + scans:
  - `python -m physics_test.cli list-gauge-Cs`
  - `python -m physics_test.cli scan-gauge-Cs --target "1/alpha"`
  - `python -m physics_test.cli scan-gauge-Cs --target "1/alpha_s(mZ)"`
  - `python -m physics_test.cli scan-gauge-Cs --target "1/alpha_w(mZ)"`
  - `python -m physics_test.cli scan-gauge-Cs --target "1/alpha1_GUT(alpha(mZ),sin2)"`

- Full pairing under gauge-derived C:
  - `python -m physics_test.cli pair-forces-gaugeCs --gravity-band any`
  - `python -m physics_test.cli pair-forces-gaugeCs --gravity-band cmb --gravity-targets "1/alpha_G(p),1/alpha_G(GW_CMB)"`
  - `python -m physics_test.cli pair-forces-gaugeCs --gravity-band pta --gravity-targets "1/alpha_G(GW_PTA)"`
  - `python -m physics_test.cli pair-forces-gaugeCs --gravity-band lisa --gravity-targets "1/alpha_G(GW_LISA)"`
  - `python -m physics_test.cli pair-forces-gaugeCs --gravity-band ligo --gravity-targets "1/alpha_G(GW_LIGO)"`

- Quantum gravity mass sweep:
  - `python -m physics_test.cli sweep-quantum-gravity --gravity-band cmb --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12`

### 7) Visualization

Open `notebooks/topological_gauge_explorer.ipynb` and run the new plotting cells to see log10(G) vs m curves for gauge-derived C candidates and highlight “hit” points within tolerance for targets like:
- `sin2thetaW(mZ)`
- `alpha_w(mZ)`
- `alpha_s(mZ)`
- `1/alpha`


