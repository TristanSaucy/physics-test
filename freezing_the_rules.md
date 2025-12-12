## Freezing the rules (so the framework stays predictive/falsifiable)

This doc is the “contract” for what we consider frozen versus exploratory.

### 1) Model (frozen)

We will treat these as fixed:

- \(\\phi = (1+\\sqrt{5})/2\)
- \(m \\in \\mathbb{Z}\) (integer harmonic steps)
- \(G\) is dimensionless and defined by:

\[
G = \\frac{C}{\\phi^m}
\]

- Frequency relation:

\[
F_0 = \\phi^m\\,\\frac{k_B K}{h}
\]

- Invariant:
\[
G\\,F_0 = C\\,\\frac{k_B K}{h}
\]

### 2) Allowed C (frozen for “strict mode”)

To avoid arbitrary fitting, strict mode uses **only gauge-derived C values** from base=360 and simple group invariants of U(1), SU(2), SU(3):

- \(C \\in \\{360, 180, 120, 60, 45, 15\\}\)

Generate/inspect:
- `python -m physics_test.cli list-gauge-Cs`

### 3) Target couplings (strong/weak candidates we consider “likely”)

We want a small menu of dimensionless targets that are common in physics literature.

#### 3.0) Frozen orientation (this is the big commitment)

We freeze the “gauge strength coordinate” as **inverse couplings**:

- Treat \(G\) as mapping to \(\\alpha^{-1}\)–type quantities (Rule A/B from the kinetic-term + RG-linearity logic).
- That means: in strict claims we use `1/alpha`-style targets for EM/strong/weak (and `1/alpha1_GUT` for hypercharge in GUT tests).

EM choice (frozen):
- **EM strict target = `1/alpha` (low-energy)**.
- `1/alpha(mZ)` is **exploratory** only (it’s scale/scheme dependent and meant for EW/GUT-context comparisons).

Direct couplings (`alpha`, `alpha_s`, `alpha_w`) are still available for exploratory runs, but **not** part of the frozen strict contract.

#### Strong (choose 1–2, but freeze which ones)

- **Benchmark**: `alpha_s(mZ)` (scale-dependent but standard reference)
- **Toy running probes** (exploratory but useful to see stability):
  - `alpha_s_1loop(2GeV)` (hadronic-ish scale)
  - `alpha_s_1loop(mZ)` (electroweak scale)

We also allow inverses as a hypothesis test (see §8 “inverse rule”):
- `1/alpha_s(mZ)`
- `1/alpha_s_1loop(Q)`

Pros/cons summary:
- **`alpha_s(mZ)`**
  - **Pros**: standard reference point; widely tabulated; easy to reproduce.
  - **Cons**: running + scheme dependence; “the strong coupling” is not a single number.
- **`1/alpha_s(mZ)`**
  - **Pros**: matches the common unification plotting convention (inverse couplings).
  - **Cons**: can become a “numerology escape hatch” if you don’t freeze orientation.

#### Weak / electroweak (choose 1–2, but freeze which ones)

- `alpha_w(mZ)` (our alpha_2 = g^2/(4π) proxy near EW scale)
- `sin2thetaW(mZ)` (dimensionless, scheme/scale dependent but widely quoted)
- `alpha2(alpha(mZ),sin2)` (derived)
- Optional hypercharge-like:
  - `alpha_Y(mZ)` and `alpha1(alpha(mZ),sin2)` (note normalization conventions)

Also available (useful for GUT-style running diagnostics):
- `alpha1_GUT(alpha(mZ),sin2)` (the \(5/3\) normalization) and its inverse

Pros/cons summary:
- **`alpha_w(mZ)` / `alpha2(alpha(mZ),sin2)`**
  - **Pros**: closely tied to \(SU(2)_L\); standard for unification-style comparisons.
  - **Cons**: depends on the input scheme (how you define/measure \(\sin^2\\theta_W\)).
- **`sin2thetaW(mZ)`**
  - **Pros**: widely quoted; clean dimensionless diagnostic of EW mixing.
  - **Cons**: not itself a coupling; interpreting it as “the weak strength” is indirect.
- **`alpha1_GUT(...)`**
  - **Pros**: puts hypercharge on the same footing as \(SU(2)\) and \(SU(3)\) for convergence tests.
  - **Cons**: normalization is a convention; can confuse interpretation unless explicitly stated.

### 4) Pass/fail criterion (frozen)

- Coupling fit threshold: **|relative error| ≤ 5%**

Command pattern:
- `python -m physics_test.cli scan-gauge-Cs --target "<target>" --max-rel-err 0.05`

### 5) Where we are on strong/weak scans (using gauge-derived C only)

These were found within 5% (examples; run scans to reproduce):

- Under the **frozen inverse** convention:
  - `1/alpha_s(mZ)`:
    - C=60, m=4 (≈ +3.21%)
  - `1/alpha_w(mZ)`:
    - C=120, m=3 (≈ -4.17%)
  - `1/alpha1_GUT(alpha(mZ),sin2)`:
    - C=60, m=0 (≈ +1.66%)

- `alpha_s(mZ)`:
  - C=60, m=13 (≈ -2.3%)
  - C=15, m=10 (≈ +3.4%)

- `alpha_w(mZ)` and `alpha2(alpha(mZ),sin2)`:
  - C=120, m=17 (≈ -0.6%)
  - C=45,  m=15 (≈ -2.4%)

- `sin2thetaW(mZ)`:
  - C=120, m=13 (≈ -0.39%)
  - C=45,  m=11 (≈ -2.2%)

- `alpha_Y(mZ)`:
  - C=60, m=18 (≈ +2.4%)

### 6) Gravity (not frozen yet)

Gravity spans broad frequency bands and multiple plausible dimensionless coupling definitions.

We will defer final gravity “rules” until:
- we decide which GW band is being targeted (CMB/PTA/LISA/LIGO), and
- we decide which gravity coupling definition corresponds to “that gravity type.”

### 7) Option 2 (phenomenon-first) vs Option 1 (energy-first)

- Option 2 (phenomenon-first): pick \(F_0\), solve \(K = F_0 h/(k_B\\phi^m)\).
- Option 1 (energy-first): pick \(K\), compute \(F_0\).

For falsifiability, Option 2 is preferred when \(F_0\) is genuinely observed (not chosen to fit).

### 8) “Principled inverse” rule (frozen once chosen per force)

Because allowing both a coupling and its inverse doubles the degrees of freedom, we treat “inverse targets” as **allowed but not free**.

Rules:
- **Per force, freeze ONE orientation**: either use the coupling (e.g. `alpha_s(mZ)`) OR its inverse (e.g. `1/alpha_s(mZ)`), but not both in strict claims.
- **A principled reason must be stated** for using an inverse target. Allowed reasons include:
  - **Unification convention**: GUT running tests are usually expressed in terms of \(\\alpha_i^{-1}\\) vs \(\\ln Q\\) (linear behavior at 1-loop).
  - **Topology-as-“resistance”** hypothesis: interpret \(G\\) as a “stiffness / action-like” measure, which more naturally maps to \(1/\\alpha\\) (large means weak interaction).
- If both orientations are explored, treat results as **exploratory** (not “frozen”), and penalize in any scoring.

### 9) Thermal anchors / “Brownian scale” (exploratory but useful)

There isn’t a single universal “Brownian frequency” (it depends on particle size and environment), but there *is* a universal **thermal frequency scale**:

- \(f_T = k_B T/h\\)

At \(T=310\\,K\):
- \(f_T \\approx 6.46\\,\\text{THz}\\) (this is exactly what `python -m physics_test.cli calc --m 0 --K 310` prints for \(F_0\\))
- at \(m=1\): \(F_0 = \\phi f_T \\approx 10.45\\,\\text{THz}\\)

These are now available as presets via:
- `python -m physics_test.cli list-frequency-presets`


