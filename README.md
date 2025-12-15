# physics-test

Tiny playground for exploring a proposed coupling:

$$
G = \frac{C}{\phi^m}
\qquad\text{and}\qquad
F_0 = \phi^m \frac{k_B K}{h}
$$

For the full manuscript-style writeup, see `paper.md`. For the frozen contract, see `freezing_the_rules.md`.
For the frozen Option‑2 anchor menu (F0 presets), see `F0_anchors.md`.
For the gravity-band/type report, see `gravity_types_report.md`.

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

Interpretive hypothesis (under the frozen inverse-coupling convention $G\sim 1/\alpha$): $m$ can be treated as a **discrete log-strength index**. In the strict best fits we see:

- strong: $m=4$
- weak: $m=3$
- EM: $m=2$
- ordinary-matter gravity (proton anchor): large negative $m$ (e.g. $\sim -175$)

See `paper.md` §6.3 for the full hypothesis and the simple $\phi$-step scaling argument.

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
- **Strong** (strict): target `1/alpha_s_1loop_from_mZ(mH)` via `C=60, m=4` (within 5%)
- **Weak** (strict): target `1/alpha2(alpha(mZ),sin2_on_shell)` via `C=120, m=3` (within 5%)

Note: `mZ` means $m_Z$, the **Z boson mass** (~91.1876 GeV). It’s a standard reference energy scale where strong/electroweak couplings are commonly quoted (because couplings “run” with energy). `mH` means $m_H$, the **Higgs boson mass** (~125 GeV), which we use as a fixed physical scale for the refined strong target.

Run:

- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha"`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha_s_1loop_from_mZ(mH)"`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha2(alpha(mZ),sin2_on_shell)"`
- `python -m physics_test.cli scan-gauge-Cs --target "1/alpha1_GUT(alpha(mZ),sin2)"`

### 5) Gravity “type” matters; CMB K implies a specific m window per GW band

If we fix gravity’s temperature $K$ to the **CMB** (default 2.725 K), then for a chosen GW frequency band the 2nd equation implies a narrow window of $m$ values (because $F_0 \propto \phi^m$).

We found:

- With CMB $K$, **CMB/primordial GW frequencies** correspond to roughly $m\approx -137$ to $-129$ (integerized).

### 6) Gravity targets: ordinary-matter anchor vs GW-band “types”

We introduced a family of (dimensionless) gravity couplings **parameterized by a mass scale** $M$:

$$
\alpha_G(M) = \frac{G_N M^2}{\hbar c}
\qquad\text{and}\qquad
\frac{1}{\alpha_G(M)} = \frac{\hbar c}{G_N M^2}
$$

This is why gravity gives a **family** of dimensionless coupling targets rather than one number: you must specify **which mass anchor** you mean (electron, proton, Planck, or an effective “GW type” anchor). Because it scales as $M^2$, changing the mass anchor changes the coupling by a fixed factor:

$$
\frac{\alpha_G(M_2)}{\alpha_G(M_1)}=\left(\frac{M_2}{M_1}\right)^2
\qquad\text{and}\qquad
\frac{1/\alpha_G(M_2)}{1/\alpha_G(M_1)}=\left(\frac{M_1}{M_2}\right)^2.
$$

In this repo, $M$ can be a particle mass (proton/electron), the Planck mass, or an effective mass scale inferred for a frozen GW-band “gravity type” (see `gravity_types_report.md`).

**Frozen (strict) orientation:** use inverse gravity targets `1/alpha_G(mass)`.

**Frozen (ordinary matter gravity):**

- Canonical target: `1/alpha_G(p)` (proton mass scale)
- Cross-check (not a free alternative): `1/alpha_G(e)` (electron mass scale)

**Frozen (GW-band gravity types; inverse coupling under CMB K):**

- `1/alpha_G(GW_CMB)` (CMB/primordial band anchor)
- `1/alpha_G(GW_PTA)` (PTA band anchor)
- `1/alpha_G(GW_LISA)` (LISA band anchor)
- `1/alpha_G(GW_LIGO)` (LIGO band anchor)
- `1/alpha_G(mP)` (Planck/quantum type; ~1)

Under strict constraints:

- **C restricted to gauge-derived values only**
- **integer m**
- **coupling fits within 5%**
- **gravity K fixed to CMB**
- **gravity F0 constrained to the chosen GW band window** (CMB/PTA/LISA/LIGO)

we found viable solutions when using inverse gravity couplings at an **effective high-energy mass scale** for GW-band “gravity types” (not a single needle; see `paper.md` for the key identity and band tables).
The strict end-to-end “all-forces per GW band” results are summarized in `gravity_types_report.md` and in `paper.md` §5.4.2.

Run:

- `python -m physics_test.cli sweep-quantum-gravity --gravity-band cmb --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12`
  - Try `--gravity-band pta`, `--gravity-band lisa`, or `--gravity-band ligo` to explore other GW “gravity types”.

## Reproduce strict core results (copy/paste)

```bash
# 360/phi^2 vs 1/alpha
python -m physics_test.cli check-example

# strict gauge-derived C set + core inverse-target scans
python -m physics_test.cli list-gauge-Cs
python -m physics_test.cli scan-gauge-Cs --target "1/alpha" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_s_1loop_from_mZ(mH)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha2(alpha(mZ),sin2_on_shell)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha1_GUT(alpha(mZ),sin2)" --max-rel-err 0.05

# strict all-forces per GW band (Option-2 anchors are frozen in F0_anchors.md)
python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band cmb  --gravity-targets "1/alpha_G(GW_CMB)"  --max-hits 10 --max-results 5
python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band pta  --gravity-targets "1/alpha_G(GW_PTA)"  --max-hits 10 --max-results 5
python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band lisa --gravity-targets "1/alpha_G(GW_LISA)" --max-hits 10 --max-results 5
python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band ligo --gravity-targets "1/alpha_G(GW_LIGO)" --max-hits 10 --max-results 5

# out-of-sample (frozen test suites)
python -m physics_test.cli oos-report --suite v1 --max-rel-err 0.02
python -m physics_test.cli oos-report --suite v2 --max-rel-err 0.02
python -m physics_test.cli oos-report --suite v3 --max-rel-err 0.02
python -m physics_test.cli oos-report --suite v4 --max-rel-err 0.02

# predictive OOS (fit one C per force from strict anchors, then hold C fixed)
python -m physics_test.cli oos-predictive --suite v1 --max-rel-err 0.02

# predictive OOS (strong + EM): within-band RG running from the lattice anchor (no re-fitting m per target)
python -m physics_test.cli oos-predictive-rg --suite v1 --max-rel-err 0.02

# step-signal OOS (C-independent): do ratios look like φ^integer?
python -m physics_test.cli oos-steps --suite v1 --max-ratio-err 0.02
python -m physics_test.cli oos-steps --suite v1 --max-ratio-err 0.05

# optional: apply a principled normalization family (see list-norm-families)
python -m physics_test.cli list-norm-families
python -m physics_test.cli oos-predictive --suite v1 --norm-family inv_C2_fund --max-rel-err 0.02

# RG / dimensional transmutation helper (shows where 'e' naturally enters via exp(-const/alpha))
python -m physics_test.cli rg-scales

# RG+phi OOS (compute Lambda_QCD implied by lattice-fit inverse strong couplings)
python -m physics_test.cli oos-rg --suite qcd-lambda-v1 --max-rel-err 0.06
python -m physics_test.cli oos-rg --suite qcd-lambda-v2 --max-rel-err 0.05
```

## How to run the main workflows

### List targets and presets

- `python -m physics_test.cli list-targets`
- `python -m physics_test.cli list-frequency-presets`
- `python -m physics_test.cli list-gravity-bands`

`list-targets` prints optional metadata (when available): **sigma (1σ)**, **reference scale Q**, and **scheme notes**.
The source of truth for these curated inputs is `data/targets.json` (see also `physics_test/target_registry.py`).

### Option 1 (energy → K → F0)

- Exploratory broad search:
  - `python -m physics_test.cli pair-forces-all --help`
  - `python -m physics_test.cli pair-forces-all`

### Option 2 (phenomenon F0 → solve K)

- Option‑2 pairing (general):
  - `python -m physics_test.cli pair-forces-option2 --help`

- Strict Option‑2 pairing (recommended; gauge-derived C only):
  - `python -m physics_test.cli pair-forces-gaugeCs --gravity-band cmb --gravity-targets "1/alpha_G(GW_CMB)"`

### Gauge-derived C only (non-arbitrary C set)

- `python -m physics_test.cli pair-forces-gaugeCs --help`

## Caveats (important)

- **Couplings run with scale / depend on scheme**: strong depends on the chosen energy scale; weak depends on both scale and $\sin^2\theta_W$ definition. Our strict targets reflect this by using a fixed 1-loop run from `alpha_s(mZ)` to `mH` and an on-shell $\sin^2\theta_W$-derived $\alpha_2$ (see `paper.md` §4.1.1).
- **Integer steps are not “all the running”**: the strongest evidence so far is that integer \(m\) behaves like a **coarse band index**, while **RG running happens within the band**. The command `oos-predictive-rg` (fit the lattice anchor once, then RG-run to other scales without re-fitting \(m\)) turns:
  - strong running OOS misses into **passes at 2%** across v2/v3/v4 strong cross-check keys (typical errors \(\sim\)1–2%), and
  - the EM OOS miss (`1/alpha → 1/alpha(mZ)`) into a **pass at 2%** under deterministic QED running (either a PDG-style Δα(mZ²) relation or a simple 1-loop threshold runner; both are available as `--runner` options).
  - weak within-band running at new scales (suite `v2`) into **passes at 2%** when using deterministic SM 1-loop running for \(1/\alpha_2(Q)\) from the on-shell-defined \(\alpha_2(m_Z)\) anchor.
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
