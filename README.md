# physics-test

Tiny playground for exploring a proposed coupling:

$$
G = \frac{C}{\phi^m}
\qquad\text{and}\qquad
F_0 = \phi^m \frac{k_B K}{h}
$$

For the full manuscript-style writeup, see `paper.md`. For the frozen contract, see `freezing_the_rules.md`.

Where:

- **G**: dimensionless “topological gauge”
- **C**: dimensionless constant
- **φ**: golden ratio
- **m**: integer index (harmonic step)
- **F0**: frequency (Hz)
- **kB**: Boltzmann constant
- **K**: temperature (Kelvin)
- **h**: Planck constant

Derived relation:

$$
G\,F_0 = C\,\frac{k_B K}{h}
$$

## What this repo is doing (in one sentence)

We’re testing whether **non-arbitrary “topological constants”** $C$ together with integer “harmonic steps” $m$ can connect **dimensionless coupling strengths** (via $G$) to **frequency/temperature scales** (via $F_0$) in a coherent way.

## Key formulas (core model)

- **Topological gauge (dimensionless)**:

$$
G = \frac{C}{\phi^m}
$$

- **Frequency coupling**:

$$
F_0 = \phi^m\,\frac{k_B K}{h}
$$

- **Invariant relation (eliminates $\phi^m$)**:

$$
G\,F_0 = C\,\frac{k_B K}{h}
$$

- **Inverse frequency → temperature (used in “Option 2”)**:

$$
K = \frac{F_0 h}{k_B\,\phi^m}
$$

## Key findings so far

### 1) 360 looks special for EM (1/alpha alignment)

- Your example is reproduced exactly:
  - $G=360/\phi^2 \approx 137.507764$
  - This is close to **$1/\alpha\approx 137.035999$** (relative error $\sim 0.34\%$).

Run:

- `python -m physics_test.cli check-example`

### 2) Integer m behaves like a harmonic ladder

We enforced **integer $m$** throughout the CLI and scans, so each $m \to m \pm 1$ is a discrete step.

### 3) Non-arbitrary C values can come from gauge invariants (base=360)

We added a small gauge-invariant generator for **U(1), SU(2), SU(3)** and derived candidate $C$ values from rank/dimension/Coxeter-type invariants using a single base “360”. This produces familiar non-arbitrary values:

- **U(1)**: $C=360$
- **SU(2)**: $360/\dim=120$, $360/(\dim\cdot h)=60$, etc.
- **SU(3)**: $360/\dim=45$, $360/(\dim\cdot h)=15$, etc.

Run:

- `python -m physics_test.cli list-gauge-Cs`

### 4) With gauge-derived C only, EM/strong/weak coupling fits land within 5% (integer m)

Scanning only gauge-derived **C candidates** (`{360, 180, 120, 60, 45, 15}`) we found `(C, m)` pairs within 5% for:

- **EM**: target `1/alpha` via `C=360, m=2`
- **Strong** (benchmark at `mZ`): target `1/alpha_s(mZ)` via `C=60, m=4` (within 5%)
- **Weak** (benchmark at `mZ`): target `1/alpha_w(mZ)` via `C=120, m=3` (within 5%)

Note: `mZ` means $m_Z$, the **Z boson mass** (~91.1876 GeV). It’s a standard reference energy scale where strong/electroweak couplings are commonly quoted (because couplings “run” with energy).

Run:

- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha"`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha_s(mZ)"`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha_w(mZ)"`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha1_GUT(alpha(mZ),sin2)"`

### 5) Gravity “type” matters; CMB K implies a specific m window per GW band

If we fix gravity’s temperature $K$ to the **CMB** (default 2.725 K), then for a chosen GW frequency band the 2nd equation implies a narrow window of $m$ values (because $F_0 \propto \phi^m$).

We found:

- With CMB $K$, **CMB/primordial GW frequencies** correspond to roughly $m\approx -137$ to $-129$ (integerized).

### 6) Gravity targets: ordinary-matter anchor vs GW-band “types”

We introduced a family of gravity couplings:

$$
\alpha_G(m) = \frac{G_N m^2}{\hbar c}
\qquad\text{and}\qquad
\frac{1}{\alpha_G(m)}
$$

where $m$ is a chosen mass scale (proton/electron, or a TeV-scale mass interpreted as $m=E/c^2$).

**Frozen (strict) orientation:** use inverse gravity targets `1/alpha_G(mass)`.

**Frozen (ordinary matter gravity):**

- Canonical target: `1/alpha_G(p)` (proton mass scale)
- Cross-check (not a free alternative): `1/alpha_G(e)` (electron mass scale)

Under strict constraints:

- **C restricted to gauge-derived values only**
- **integer m**
- **coupling fits within 5%**
- **gravity K fixed to CMB**
- **gravity F0 constrained to CMB GW band**

we found viable solutions when using inverse gravity couplings at an **effective high-energy mass scale** for GW-band “gravity types” (not a single needle; see `paper.md` for the key identity and band tables).

Run:

- `python -m physics_test.cli sweep-quantum-gravity --gravity-band cmb --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12`
  - Try `--gravity-band pta`, `--gravity-band lisa`, or `--gravity-band ligo` to explore other GW “gravity types”.

## How to run the main workflows

### List targets and presets

- `python -m physics_test.cli list-targets`
- `python -m physics_test.cli list-frequency-presets`
- `python -m physics_test.cli list-gravity-bands`

### Option 1 (energy → K → F0)

- `python -m physics_test.cli pair-forces-all ...`

### Option 2 (phenomenon F0 → solve K)

- `python -m physics_test.cli pair-forces-option2 ...`

### Gauge-derived C only (non-arbitrary C set)

- `python -m physics_test.cli pair-forces-gaugeCs ...`

## Caveats (important)

- **Couplings run with scale**: strong/weak (and effective EM) depend on the energy scale; targets like `alpha_s(mZ)` (and therefore `1/alpha_s(mZ)`) are benchmark values at a reference scale.
- **What “K” means**: for micro physics, treating $K$ as an “energy-scale temperature” (i.e., $E=k_BK$) is often more coherent than interpreting it as literal thermodynamic temperature.
- **Avoiding overfitting**: as you allow larger families of C values, matches become easier. The most meaningful tests are the strict ones: discrete, justified C sets + fixed target definitions + fixed band constraints.

## Note on "different alpha" and the 4 forces (why we scan options)

- **EM**: The fine-structure constant **alpha in vacuum** is not different for visible vs X-ray photons. What *can* look different is an **effective coupling** due to:
  - **material response** (screening/dispersion in matter), and/or
  - **scale dependence** in high-energy processes (people write alpha(Q)).
- **Strong/weak**: their dimensionless couplings are **scale-dependent** (they “run” with energy), so there isn't one single number unless you specify a reference scale.
- **Gravity**: a common dimensionless coupling is $\alpha_G(m)=G_N m^2/(\hbar c)$, which depends on the chosen mass scale $m$.

This repo supports scanning **multiple target options** (different scales / definitions) to see what—if anything—lands on simple $(C,m)$ patterns.

## Quickstart

Install deps (pick one):

- `pip install -r requirements.txt`
- `pip install -e .`

Run:

- `python -m physics_test.cli check-example`
- `python -m physics_test.cli calc --C 360 --m 2 --K 300`
- `python -m physics_test.cli fits --m 2`
- `python -m physics_test.cli solve-K --m 2 --F0 1e13`
- `python -m physics_test.cli list-targets`
- `python -m physics_test.cli scan-all --set rotation-degrees --m-min -6 --m-max 6 --m-step 1 --m-integer`
