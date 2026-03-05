# Phi-Lattice Framework

**Status**: exploratory (not peer-reviewed)
**Repository**: `physics-test`
**Date**: 2026-01

---

## 1. What This Is (and What It Isn't)

**What it IS:**

- A minimal, auditable toy framework that asks: "Can a very simple discrete rule land near real, dimensionless physics numbers?"
- A codebase that freezes rules up front (so it cannot wiggle forever to fit anything).
- A way to separate a discrete hypothesis ("allowed values come from a short integer ladder") from the known continuous physics of running couplings with energy.

**What it is NOT:**

- Not a replacement for the Standard Model.
- Not a precision electroweak calculation.
- Not peer-reviewed.
- Not a claim that "phi explains physics." It is a hypothesis under stress-test.

**In one sentence:** this project tests whether a tiny discrete ladder $C/\phi^m$, with $C$ restricted by gauge-group integers and $m$ forced to be an integer, can serve as a "band selector," while ordinary renormalization-group running handles the smooth scale dependence within the band.

---

## 2. Core Model Definitions

### 2.1 The golden ratio

$$
\phi \equiv \frac{1+\sqrt{5}}{2} \approx 1.6180339887
$$

### 2.2 The phi-lattice

$$
G(C,m) = \frac{C}{\phi^m}, \qquad m \in \mathbb{Z}
$$

### 2.3 Frequency coupling

$$
F_0(m,K) = \phi^m \, \frac{k_B K}{h}
$$

### 2.4 Invariant relation

Eliminating $\phi^m$ from the two equations above yields:

$$
G \cdot F_0 = C \, \frac{k_B K}{h}
$$

### 2.5 Option-2 inversion (phenomenon-first)

Given an observed phenomenon frequency $F_0$, solve for temperature:

$$
K = \frac{F_0 \, h}{k_B \, \phi^m}
$$

### 2.6 Interpretation of the harmonic index m

From the frequency equation:

$$
\phi^m = \frac{h F_0}{k_B K}
$$

So $m$ is the base-$\phi$ logarithm of a physically meaningful ratio: quantum energy scale $hF_0$ relative to thermal scale $k_B K$.

- **$m > 0$**: $hF_0 > k_B K$ (phenomena "above thermal scale" at temperature $K$).
- **$m < 0$**: $hF_0 < k_B K$ (phenomena "below thermal scale" at temperature $K$).
- **$m = 0$**: baseline thermal frequency $f_T = k_B T / h$.

### 2.7 Parameter summary

| Symbol | Meaning |
|---|---|
| $G$ | Dimensionless "topological gauge" (coupling coordinate) |
| $C$ | Discrete "topological constant" from a restricted menu |
| $\phi$ | Golden ratio |
| $m$ | Integer band/harmonic index |
| $F_0$ | Frequency scale (Hz) |
| $k_B$ | Boltzmann constant |
| $K$ | Temperature (Kelvin) |
| $h$ | Planck constant |

---

## 3. Constants and Conventions

### 3.1 Coupling constants

$$
\alpha_i = \frac{g_i^2}{4\pi}
$$

$\alpha_i^{-1}$ is the inverse coupling coordinate used throughout this framework.

### 3.2 Dimensionless gravitational coupling

$$
\alpha_G(M) = \frac{G_N M^2}{\hbar c}, \qquad \frac{1}{\alpha_G(M)} = \frac{\hbar c}{G_N M^2}
$$

Unlike gauge couplings, $\alpha_G$ depends on the choice of mass scale $M$.

### 3.3 Why inverse couplings

- Gauge kinetic terms naturally involve $1/g^2$, and $\alpha \propto g^2$.
- At 1-loop order, inverse couplings run approximately linearly with $\ln Q$.
- Unification-style plots are traditionally drawn using $\alpha_i^{-1}$.
- Freezing one orientation (inverse) prevents doubling the hypothesis space.

---

## 4. Strict Contract (Anti-Overfitting Rules)

The biggest danger with pattern-hunting is numerology. This contract freezes key choices so failures actually mean something.

### 4.1 Frozen: integer stepping

$$
m \in \mathbb{Z} \text{ everywhere (no real-valued } m \text{ fitting)}
$$

### 4.2 Frozen: inverse-coupling orientation

Strict claims use **inverse couplings** ($1/\alpha$-style targets). Per force, one orientation is frozen. Allowed reasons for using inverses include the unification convention and the topology-as-"resistance" hypothesis.

### 4.3 Frozen: gauge-derived C menu

Generate candidates from a base value 360 and simple invariants of the Standard Model gauge groups U(1), SU(2), SU(3):

$$
C \in \{360, 180, 120, 60, 45, 15\}
$$

### 4.4 Frozen: base selection rule

**Rule**: choose the **smallest highly-composite integer base** such that the strict gauge-derived candidate set stays **short after de-duplication** and contains values that hit the strict EM/strong/weak anchor targets within tolerance using integer $m$, with the **micro-force anchor fits constrained to nonnegative $m$**.

Motivations:

- **Geometry/phase**: $360°$ is the conventional integer representation of a full $2\pi$ cycle.
- **Divisor richness**: 360's factorization makes it highly divisible, naturally yielding a small menu of rationally-related candidates.

### 4.5 Frozen: base-vs-alt-bases experiment protocol

Alternative bases (e.g., 420, 840, ...) belong to a separate, explicitly labeled experiment, not a rescue knob.

- **Protocol**: fix the invariant menu and strict target set; scan a pre-declared list of candidate bases and evaluate: size of the deduplicated gauge-derived $C$ set, strict-anchor hit rates at the frozen tolerance, and (optionally) survival on frozen OOS suites.
- **Decision rule**: the strict base is the smallest base meeting the criteria; changing it requires updating the frozen contract and re-running the pre-registered evaluation.

```bash
python -m physics_test.cli base-vs-alt-bases --tol 0.05 --max-nCs 10
```

### 4.6 Frozen: pass/fail thresholds

- Coupling fit threshold: $|\text{relative error}| \le 5\%$
- OOS tests typically reported at $2\%$
- Broader exploratory scan threshold: $5\%$

### 4.7 Frozen: target definitions

- **EM strict target**: `1/alpha` (low-energy fine-structure inverse). `1/alpha(mZ)` is exploratory only.
- **Strong strict target**: `1/alpha_s_1loop_from_mZ(mH)` (1-loop run from $\alpha_s(m_Z)$ to $m_H = 125$ GeV; no free $\Lambda_{\text{QCD}}$ knob).
- **Weak strict target**: `1/alpha2(alpha(mZ),sin2_on_shell)` (derived using on-shell $\sin^2\theta_W = 1 - m_W^2/m_Z^2$).
- **Gravity**: inverse coupling `1/alpha_G(\text{mass})` with frozen mass anchors per type.

### 4.8 Frozen: K interpretation

- **Gravity (GW types)**: $K$ is treated as **literal CMB temperature** (default 2.725 K).
- **Micro forces (EM/strong/weak)**: $K$ is treated as an **effective scale parameter** ($E = k_B K$), not as an ambient thermodynamic temperature.

### 4.9 Measurement registry

Key external inputs are stored in `data/targets.json`. Each entry includes (where known): **value**, **1$\sigma$ uncertainty**, optional **reference scale** $Q$, a short **scheme/notes** string, and a **citation hint**. The code reads this via `physics_test/target_registry.py`.

Any registry key of the form `tgt_<target_name>` is automatically surfaced by `list-targets` as a dimensionless target.

---

## 5. Where C Comes From: Gauge-Group Invariants

### 5.1 Purpose

Restrict $C$ to a **small, externally-motivated set** that comes from unambiguous properties of the gauge group, is reproducible from a few lines of math, and stays short after de-duplication.

### 5.2 Lie algebra invariants used

These are invariants of the **Lie algebra / root system** (representation-independent).

**Rank** $r$: dimension of a maximal torus. For $SU(N)$: $r = N-1$.

**Lie algebra dimension** $\dim(\mathfrak{g})$: number of generators (dimension of the adjoint representation). For $SU(N)$: $\dim = N^2 - 1$.

**Coxeter number** $h$: an integer invariant of the root system, defined as $h = 1 +$ (height of the highest root). For $SU(N)$: $h = N$.

**Dual Coxeter number** $h^\vee$: under the standard normalization (long roots have length-squared 2), equals the quadratic Casimir of the adjoint: $C_2(\text{adj}) = h^\vee$. For simply-laced groups (including $SU(N)$): $h = h^\vee$.

| Group | Rank $r$ | $\dim(\mathfrak{g})$ | $h$ | $h^\vee$ |
|---|---:|---:|---:|---:|
| U(1) | 1 | 1 | -- | -- |
| SU(2) | 1 | 3 | 2 | 2 |
| SU(3) | 2 | 8 | 3 | 3 |

### 5.3 Strict construction rules

The current strict generator takes base = 360 and produces:

- $C = \text{base}$
- $C = \text{base}/\dim$
- $C = \text{base}/h$
- $C = \text{base}/h^\vee$
- $C = \text{base}/(\dim \cdot h)$

Concrete values:

- **U(1)**: $360$
- **SU(2)**: $360/3 = 120$, $360/2 = 180$, $360/(3 \cdot 2) = 60$
- **SU(3)**: $360/8 = 45$, $360/3 = 120$, $360/(8 \cdot 3) = 15$

After de-duplication: $C \in \{360, 180, 120, 60, 45, 15\}$.

### 5.4 Additional root-system invariants (representation-independent)

**Number of roots**: $|\Phi| = \dim(\mathfrak{g}) - r$. For $SU(N)$: $|\Phi| = N(N-1)$.

**Number of positive roots**: $|\Phi^+| = (\dim(\mathfrak{g}) - r)/2$.

**Weyl group order** $|W|$: for $SU(N)$, $|W| = N!$.

**Center order** $|Z(G)|$: for $SU(N)$, $|Z| = N$.

### 5.5 Representation-dependent invariants (Casimirs and Dynkin indices)

**Quadratic Casimir** $C_2(R)$:

$$
\sum_a T_R^a T_R^a = C_2(R) \, \mathbf{1}
$$

**Dynkin index** $T(R)$:

$$
\text{Tr}_R(T^a T^b) = T(R) \, \delta^{ab}
$$

Related by:

$$
\dim(\mathfrak{g}) \, T(R) = \dim(R) \, C_2(R)
$$

For $SU(N)$ under standard normalization:

- Adjoint: $C_A = C_2(\text{adj}) = N$
- Fundamental: $C_F = C_2(F) = (N^2-1)/(2N)$
- $T(F) = 1/2$

These are dangerous for "$C$ mining" because you must pick a representation and a normalization convention. The extended `--include` menu (with `base/C2_fund`, `base/T_fund`, etc.) yields $C \in \{15, 45, 60, 120, 180, 270, 360, 480, 720\}$.

### 5.6 Principled expansion ruleset

If extending beyond the current menu, constraints to keep the list short:

1. **No free reps**: freeze one representation choice per gauge factor.
2. **No free normalizations**: freeze a convention (e.g., long roots have length-squared 2; $T(F) = 1/2$).
3. **Simple arithmetic only**: allow only base divided by a single invariant or by a product of at most two invariants.
4. **Integer/small-denominator filter**: require candidates to be integers or rationals with very small denominator.
5. **Deduplicate aggressively**: if multiple derivations yield the same $C$, keep the value once.
6. **Freeze the menu**: decide up front which invariants are allowed in strict mode; treat all others as exploratory.

### 5.7 CLI construction keys

Strict mode defaults: `base`, `base/dim`, `base/coxeter`, `base/dual_coxeter`, `base/(dim*coxeter)`.

Next-tier keys (exploratory): `base/rank`, `base/roots`, `base/positive_roots`, `base/weyl_order`, `base/center_order`, `base/C2_adj`, `base/C2_fund`, `base/T_fund`.

---

## 6. Target Definitions

### 6.1 EM (electromagnetism)

- **Strict anchor** (low energy / vacuum): `1/alpha`
- **OOS cross-check**: `1/alpha(mZ)` (effective coupling at the Z pole)

### 6.2 Strong (QCD)

Because $\alpha_s(Q)$ runs strongly, targets must specify scale and prescription.

- **Strict anchor**: `1/alpha_s_1loop_from_mZ(mH)` -- take benchmark $\alpha_s(m_Z)$ and run to $m_H = 125$ GeV using 1-loop QCD running with fixed $n_f = 5$ (no free $\Lambda_{\text{QCD}}$ knob).
- **OOS cross-checks**: `1/alpha_s_1loop_from_mZ(mW)`, `(...mt)`, `(...1TeV)`, `(...10TeV)`, plus `nf56` and `2loop` variants.

### 6.3 Weak / electroweak

Weak targets are scheme-sensitive (especially $\sin^2\theta_W$).

- **Strict anchor**: `1/alpha2(alpha(mZ),sin2_on_shell)` using the on-shell mixing angle:

$$
\sin^2\theta_W^{\text{OS}} \equiv 1 - \frac{m_W^2}{m_Z^2}, \qquad \alpha_2 \equiv \frac{\alpha(m_Z)}{\sin^2\theta_W^{\text{OS}}}
$$

- **Hypercharge (GUT-normalized)**: `1/alpha1_GUT(alpha(mZ),sin2_on_shell)` using the $5/3$ normalization factor.

### 6.4 Gravity

Frozen gravity targets and types:

| Target | Meaning | Mass anchor |
|---|---|---|
| `1/alpha_G(p)` | Ordinary matter (proton) | $m_p$ |
| `1/alpha_G(e)` | Cross-check (electron) | $m_e$ |
| `1/alpha_G(GW_CMB)` | CMB/primordial band | $2.93012 \times 10^4$ GeV |
| `1/alpha_G(GW_PTA)` | PTA band | $1.58009 \times 10^9$ GeV |
| `1/alpha_G(GW_LISA)` | LISA band | $2.15524 \times 10^{12}$ GeV |
| `1/alpha_G(GW_LIGO)` | LIGO band | $5.41086 \times 10^{13}$ GeV |
| `1/alpha_G(mP)` | Planck/quantum type | $\sim 10^{19}$ GeV |

The electron target is a mandatory cross-check (not a free alternative); it should reproduce an $m$-shift consistent with $(m_p/m_e)^2$.

### 6.5 Reference scale: mZ

$m_Z \approx 91.1876$ GeV (Z boson mass). Standard reference energy scale at which $\alpha_s$ and electroweak parameters are commonly quoted.

---

## 7. RG Running ("Within-Band") Formulas

The key insight: integer $m$ is a **coarse band index** (discrete anchor placement), while **RG flow supplies within-band motion** (continuous in $\ln Q$).

### 7.1 Generic 1-loop gauge running

$$
\frac{d(\alpha_i^{-1})}{d \ln \mu} = -\frac{b_i}{2\pi}
$$

$$
\alpha_i^{-1}(\mu) = \alpha_i^{-1}(\mu_0) - \frac{b_i}{2\pi} \ln\!\left(\frac{\mu}{\mu_0}\right)
$$

SM and MSSM 1-loop coefficients:

| Model | $b_1$ | $b_2$ | $b_3$ |
|---|---:|---:|---:|
| SM (1-loop) | $41/10$ | $-19/6$ | $-7$ |
| MSSM (1-loop) | $33/5$ | $1$ | $-3$ |

### 7.2 QCD running (strong sector)

**1-loop** (inverse coupling runs linearly in $\ln Q$):

$$
\alpha_s^{-1}(Q) = \alpha_s^{-1}(Q_0) + \frac{b_0}{2\pi} \ln\!\left(\frac{Q}{Q_0}\right), \qquad b_0 = 11 - \frac{2}{3} n_f
$$

**2-loop** (differential equation integrated numerically):

$$
\frac{d\alpha_s}{d\ln Q} = -\frac{\beta_0}{2\pi}\alpha_s^2 - \frac{\beta_1}{4\pi^2}\alpha_s^3
$$

$$
\beta_0 = 11 - \frac{2}{3} n_f, \qquad \beta_1 = 102 - \frac{38}{3} n_f
$$

**Threshold structure** (minimal model): $n_f$ changes when crossing heavy-quark masses (e.g. switch $5 \to 6$ at $Q \simeq m_t$). Treated as a deterministic option (no matching/fitting).

### 7.3 QED running (EM sector)

$$
\frac{d(\alpha^{-1})}{d\ln Q} = -\frac{2}{3\pi} \sum_f N_c \, Q_f^2
$$

The sum includes only fermions with mass $m_f < Q$ (sharp thresholds). This is not precision hadronic vacuum polarization, but it is deterministic and captures the correct direction/scale of running. For Z-pole work, registry-provided $\Delta\alpha(m_Z^2)$ pieces are also available.

### 7.4 Predictive modes (how running is applied)

Three tiers of testing:

- **OOS report** (`oos-report`): lets $C$ vary per target within the strict menu.
- **Predictive OOS** (`oos-predictive`): fits one best $C$ per force on an anchor, then holds $C$ fixed but allows integer $m$ to vary per target.
- **RG-predictive OOS** (`oos-predictive-rg`): fits the anchor once, then predicts other scales via RG **with no re-fitting of $m$**. This is the closest thing to "band + within-band running."

---

## 8. Electroweak Mixing and Z-Pole Mapping

### 8.1 On-shell weak mixing angle

$$
\sin^2\theta_W^{\text{OS}} = 1 - \frac{m_W^2}{m_Z^2}
$$

### 8.2 Coupling relations

$$
\alpha = \alpha_2 \, \sin^2\theta_W = \alpha_Y \, \cos^2\theta_W
$$

$$
\alpha_1^{\text{GUT}} = \frac{5}{3} \, \alpha_Y
$$

### 8.3 Inverse-coupling form for the mixing angle

$$
\sin^2\theta_W = \frac{\alpha_Y}{\alpha_2 + \alpha_Y} = \frac{\alpha_2^{-1}}{\alpha_2^{-1} + \alpha_Y^{-1}}
$$

$$
\alpha_Y^{-1} = \frac{5}{3} \, \alpha_{1,\text{GUT}}^{-1} \implies \sin^2\theta_W = \frac{\alpha_2^{-1}}{\alpha_2^{-1} + \frac{5}{3}\alpha_{1,\text{GUT}}^{-1}}
$$

### 8.4 Z-pole effective leptonic weak mixing angle

Z-pole measurements quote $\sin^2\theta_{\text{eff}}^{\text{lept}}$, a pseudo-observable defined from effective Z-lepton couplings:

$$
\sin^2\theta_{\text{eff}}^{\text{lept}} = \frac{1}{4}\left(1 - \frac{v_l}{a_l}\right)
$$

This folds in full EW radiative corrections by definition. Comparing a toy model directly to this number is not meaningful ("fails by construction").

### 8.5 Z-pole mapping method: zpole_kappa_approx

**Step A** (reference angle from lattice-anchored gauge couplings at $Q = m_Z$):

$$
\sin^2\theta_{\text{ref}}(m_Z) = \frac{\alpha_2^{-1}(m_Z)}{\alpha_2^{-1}(m_Z) + \alpha_Y^{-1}(m_Z)}
$$

**Step B** (map to the Z-pole pseudo-observable):

$$
\sin^2\theta_{\text{eff}}^{\text{lept}}(m_Z) \approx \kappa_Z \, \sin^2\theta_{\text{ref}}(m_Z)
$$

Toy but auditable approximation:

$$
\kappa_Z \approx 1 + \Delta\alpha_{\text{total}}(m_Z^2) - \Delta\rho_{\text{top}}
$$

$$
\Delta\alpha_{\text{total}}(m_Z^2) = \Delta\alpha_{\text{lept}} + \Delta\alpha_{\text{had}}^{(5)} + \Delta\alpha_{\text{top}}
$$

$\Delta\rho_{\text{top}}$ is the leading top-loop rho correction (computed/recorded as a target in the repo).

### 8.6 Frozen EW suites

**CI-grade (gating):**

- `ew-independent-v1`: Qweak, E158, $\nu_e$-lowE. Default: $|z| \le 3$.
- `ew-independent-v2`: extends v1 with Cs APV. Default: $|z| \le 2$.
- `ew-independent-v3`: extends v2 with eDIS. Default: $|z| \le 2$.

All use `--method gammaZ_1loop` and scheme prefix `sin2thetaW_eff:`.

**Exploratory (non-gating):**

- `ew-exploratory-v1`: extends v3 with interpretation-sensitive points (PDG Cs APV, CHARM II).
- `ew-dis-exploratory-v1`: NuTeV-style DIS extractions. Scheme prefix: `sin2thetaW_on_shell:`.
- `ew-zpole-exploratory-v1`: LEP/SLC and Tevatron Z-pole. Scheme prefix: `sin2thetaW_eff_lept:`. Frozen to `--method zpole_kappa_approx`.

### 8.7 NuTeV note

NuTeV extracted a weak mixing angle from neutrino scattering on a nuclear target and originally appeared "off" from SM expectations by a few sigma. The extracted number depends on nuclear effects, PDF uncertainties, strange/anti-strange quark distributions, charm threshold effects, and EW radiative corrections. Treated as exploratory-only.

---

## 9. Frequency and Temperature Anchors (Option-2)

### 9.1 Why $F_0 = E/h$ is a proxy

For **photons**, this is literally measured (spectral lines). For **particle masses / interaction scales**, $E/h$ is a timescale proxy rather than a directly observed oscillation. The "con" is that it introduces an extra modeling choice (which energy scale $E$ represents "the phenomenon"), which becomes an overfitting knob if allowed to float freely.

### 9.2 Principled selection rules

- **Rule P (pole/resonance)**: use a pole mass or sharp resonance energy.
- **Rule L (Lagrangian-parameter)**: use a fundamental scale parameter in the theory.
- **Rule M (measurement-context)**: use the energy scale at which the coupling target is defined.

Recommendation for falsifiability: use **Rule P** for EM/weak, and a hybrid of **Rule P** and **Rule L** for strong.

### 9.3 Frozen anchor menu

One primary anchor per force + 1-2 cross-checks (not free alternatives):

**EM** (F0 is directly measurable):

| Key | Description | Source |
|---|---|---|
| `em-lyman-alpha` (primary) | Hydrogen Lyman-$\alpha$ line ($\lambda \approx 121.6$ nm) | NIST ASD |
| `em-visible-500THz` | Broad human-visible band anchor | -- |
| `em-hydrogen-13.6eV` | Atomic energy scale $E/h$ | NIST |

**Strong** (energy/time proxies):

| Key | Description | Source |
|---|---|---|
| `strong-QCD-200MeV` (primary) | QCD/confinement-scale proxy ($\hbar c / 1\text{ fm} \approx 200$ MeV) | PDG |
| `strong-proton-938MeV` | Proton mass-energy scale | PDG |
| `strong-timescale-1e-23s` | Historical strong-interaction timescale | -- |

**Weak** (mediator mass proxies):

| Key | Description | Source |
|---|---|---|
| `weak-W-80.379GeV` (primary) | W boson pole mass | PDG |
| `weak-Z-91.1876GeV` | Z boson pole mass | PDG |
| `weak-muon-decay` | Muon lifetime rate proxy $F_0 \sim 1/\tau_\mu$ | PDG |

**Gravity**: not a single $F_0$ -- frozen by **bands/types**; $F_0$ constrained by the observational window (LIGO/LISA/PTA/CMB).

### 9.4 Blackbody peak factors (exploratory)

Two principled alternatives to $k_B T / h$:

- Peak in frequency-domain Planck spectrum: $\nu_{\text{peak}} = x \, k_B T / h$ with $x \approx 2.821439$.
- Peak in wavelength-domain spectrum (converted to frequency): $c/\lambda_{\text{peak}} = y \, k_B T / h$ with $y \approx 4.965114$.

---

## 10. Key Findings

### 10.1 Discrete anchor hits (strict gauge-derived C)

Under strict gauge-derived $C$ and integer $m$:

| Sector | Target | Target value | Best $(C,m)$ | $G(C,m)$ | Rel. err |
|---|---|---:|---:|---:|---:|
| EM | `1/alpha` | 137.035999 | (360, 2) | 137.507764 | +0.344% |
| Strong | `1/alpha_s_1loop_from_mZ(mH)` | 8.866605 | (60, 4) | 8.753882 | -1.27% |
| Weak | `1/alpha2(alpha(mZ),sin2_on_shell)` | 28.535657 | (120, 3) | 28.328157 | -0.73% |
| Hypercharge | `1/alpha1_GUT(alpha(mZ),sin2_on_shell)` | 59.651606 | (60, 0) | 60.000000 | +0.584% |

The EM hit is the original motivating example: $360/\phi^2$ lands near $1/\alpha$.

### 10.2 Integer-steps-only failure mode

When you try to explain scale-dependence using only integer $\Delta m$ steps, you predictably miss running-sensitive targets (strong near $m_W$/$m_t$/1 TeV, EM at $m_Z$). This is expected: couplings run continuously with $\ln Q$.

### 10.3 RG-within-band resolves those misses (no new fit knobs)

Using `oos-predictive-rg` (fit the anchor once, then RG-run to other scales with no re-fitting of $m$):

- **Strong**: running cross-checks become **passes at 2%** (typical errors ~1-2% across v2/v3/v4 strong running keys).
- **EM**: OOS cross-check `1/alpha -> 1/alpha(mZ)` becomes a **pass at 2%** under deterministic QED running.
- **Weak**: within-band running of $\alpha_2^{-1}(Q)$ becomes **passes at 2%** across $m_W$, $m_H$, 1 TeV, 10 TeV (suite `v2`).
- **Hypercharge**: within-band running of $\alpha_{1,\text{GUT}}^{-1}(Q)$ becomes **passes at 2%** across $m_W$, $m_H$, 1 TeV, 10 TeV (suite `v3`).
- **EW mixing**: derived $\sin^2\theta_W(Q)$ from $\alpha_2 + \alpha_{1,\text{GUT}}$ running passes at **2%** across $m_W$, $m_H$, 1 TeV, 10 TeV (typical offset ~1%).

Interpretation: $m$ behaves like a **coarse band index**, while **RG flow supplies within-band motion**.

### 10.4 m as a discrete strength coordinate (force hierarchy)

Under the frozen inverse-coupling convention $G \sim 1/\alpha$:

$$
G = \frac{C}{\phi^m} \implies \alpha \sim \frac{\phi^m}{C}
$$

$$
\log_\phi(\alpha) \sim m - \log_\phi(C)
$$

Holding $C$ fixed, each +1 step in $m$ multiplies $\alpha$ by $\phi \approx 1.618$. The strict best fits:

- Strong: $(C,m) = (60, 4)$
- Weak: $(C,m) = (120, 3)$
- EM: $(C,m) = (360, 2)$

This ordering ($\alpha_s > \alpha_w > \alpha$) is consistent with the usual hierarchy. For adjacent pairs:

$$
\frac{G_{\text{weak}}}{G_{\text{strong}}} = \frac{C_w}{C_s} \phi^{m_s - m_w} \approx 2\phi, \qquad \frac{G_{\text{EM}}}{G_{\text{weak}}} = \frac{C_{\text{em}}}{C_w} \phi^{m_w - m_{\text{em}}} \approx 3\phi
$$

A one-step $m$ shift plus a small integer $C$ ratio produces an $\mathcal{O}(1\text{-}10)$ hierarchy across the gauge forces.

For **ordinary-matter gravity** (proton anchor) at $(C,m) \approx (45, -175)$:

$$
\frac{G_{\text{grav}}}{G_{\text{EM}}} = \frac{C_g}{C_{\text{em}}} \phi^{m_{\text{em}} - m_g} \approx \frac{1}{8} \phi^{177} \sim 10^{36}
$$

Hypothesis: gauge forces occupy **nearby positive-$m$ levels** (few steps apart), while ordinary-matter gravity occupies a **far negative-$m$ level**, producing exponential weakness.

### 10.5 Delta-r lattice alignment (exploratory)

The on-shell electroweak radiative correction parameter $\Delta r$ from $\alpha(0), G_F, m_W, m_Z$:

- `1/delta_r(on-shell;alpha0,GF,mW,mZ)` $\approx 28.244$
- Best strict hit: $(C,m) = (120,3) \Rightarrow 120/\phi^3 \approx 28.328$ ($\approx +0.30\%$).

This does not mean the phi-lattice "explains $\Delta r$"; it is a notable numerical coincidence.

### 10.6 Radiative-piece structure (exploratory)

Suite `v6` probes inverse vacuum-polarization pieces:

- `1/delta_alpha_total(mZ2)` best strict hit: $(C,m) = (45, 2)$ ($\approx +1.53\%$).
- `1/delta_r(on-shell;alpha0,GF,mW,mZ)` best strict hit: $(C,m) = (120, 3)$ ($\approx +0.30\%$).

At 2% tolerance, 2/4 pass; at 5% tolerance, 4/4 pass.

Suite `v7`: `1/delta_rho_top(GF,mt)` lands within ~4% of strict lattice points (passes at 5% but not 2%).

### 10.7 External EW cross-check via G_F (exploratory)

Suite `v5` adds `1/alpha2_tree_from_GF(mW)` (tree-level extraction from $G_F, m_W$ with no $\Delta r$). Best strict hit is still $(C,m) = (120, 3)$ but off by ~3.9% (fails at 2%, passes at 5%). This is expected to be dominated by electroweak radiative corrections ($\Delta r$).

### 10.8 Thermal baseline sanity check

The model's universal thermal frequency scale at $m = 0$:

$$
f_T = \frac{k_B T}{h}
$$

At $T = 310$ K (human body temperature): $f_T \approx 6.46$ THz. At $m = 1$: $F_0 = \phi f_T \approx 10.45$ THz. This is a clean anchor (no new physics -- only a direct conversion between temperature and frequency).

### 10.9 Bio-adjacent probes (exploratory)

- "Brownian motion" does not define a single characteristic frequency without specifying a system. If a biology-adjacent anchor is needed, $k_B T / h$ is universal.
- Many biomolecules have broadband THz features (collective vibrational/rotational modes). A meaningful "hit" would require a reproducible narrow resonance or evidence for non-thermal coherence in biologically realistic conditions.
- If $F_0$ is in the ~6-12 THz range (thermal baseline at 300-310 K), the falsifiable question becomes whether strict $C$ candidates + frozen target choices yield an $m$ predicting a THz-band $F_0$ at realistic $K$ without tuning.

---

## 11. Gravity Extension

### 11.1 Dimensionless gravitational coupling

$$
\alpha_G(M) = \frac{G_N M^2}{\hbar c}, \qquad \frac{1}{\alpha_G(M)} = \frac{\hbar c}{G_N M^2}
$$

Unlike gauge couplings, $\alpha_G$ depends on choice of $M$, motivating different "gravity types" (macro, GW-band, quantum) corresponding to different frozen mass anchors.

### 11.2 Electron vs proton mass-ratio identity

Since $\alpha_G^{-1}(M) \propto 1/M^2$:

$$
\frac{1/\alpha_G(e)}{1/\alpha_G(p)} = \left(\frac{m_p}{m_e}\right)^2 \approx 3.37 \times 10^6, \qquad \log_\phi\!\left(\left(\frac{m_p}{m_e}\right)^2\right) \approx 31.24
$$

If $C$ were held fixed, moving from proton to electron as the mass anchor would shift the implied harmonic index by ~31 steps.

### 11.3 Predicted m-shift (including discrete C correction)

If $G_t \approx C/\phi^m$, then the best-fit $m$ is:

$$
m \approx \log_\phi\!\left(\frac{C}{G_t}\right) = \log_\phi(C) - \log_\phi(G_t)
$$

For two targets $G_{t,1}, G_{t,2}$ and two allowed constants $C_1, C_2$:

$$
\Delta m = m_2 - m_1 \approx \log_\phi\!\left(\frac{C_2}{C_1}\right) - \log_\phi\!\left(\frac{G_{t,2}}{G_{t,1}}\right)
$$

For electron vs proton inverse gravity: $G_{t,2}/G_{t,1} = (m_p/m_e)^2$. In the strict set, the best proton hit uses $C_1 = 45$ while the electron hit uses $C_2 = 360$, so $C_2/C_1 = 8$:

$$
\Delta m \approx \log_\phi(8) - \log_\phi\!\left(\left(\frac{m_p}{m_e}\right)^2\right) \approx 4.32 - 31.24 \approx -26.9
$$

This matches the observed strict hits ($m_p$ around -175, $m_e$ around -202), i.e. a shift of -27 integer steps, within rounding/discreteness.

### 11.4 Ordinary matter gravity fits

Under strict gauge-derived $C$ and 5%:

- `1/alpha_G(p)`: hits at $m \approx -175$ (e.g. $C = 45$) and $m \approx -173$ (e.g. $C = 120$).
- `1/alpha_G(e)`: hit at $m \approx -202$ (e.g. $C = 360$).

Freeze: "ordinary matter gravity" uses the proton anchor `1/alpha_G(p)`, with electron `1/alpha_G(e)` as a required cross-check.

### 11.5 GW-band definitions

Order-of-magnitude sensitivity windows of different observational techniques:

| Band | Frequency range | Platform |
|---|---|---|
| LIGO | 10-1000 Hz | Ground-based laser interferometers (LIGO Hanford/Livingston, Virgo, KAGRA) |
| LISA | $10^{-4}$-$10^{-1}$ Hz | Planned space-based interferometer (three spacecraft, heliocentric orbit) |
| PTA | $10^{-9}$-$10^{-7}$ Hz | Pulsar timing arrays (NANOGrav, EPTA, PPTA; combined under IPTA) |
| CMB/primordial | $\sim 10^{-18}$-$10^{-16}$ Hz | Indirect constraints from CMB anisotropy/polarization (not direct detection) |

### 11.6 Band-implied m windows under CMB K

With $K = 2.725$ K:

| Band | $m$ window |
|---|---|
| CMB | $\approx -137$ to $-129$ |
| PTA | $\approx -94$ to $-85$ |
| LISA | $\approx -70$ to $-57$ |
| LIGO | $\approx -46$ to $-38$ |

### 11.7 Band-constrained mass-scale solutions

Searching for mass scales $M$ where EM/strong/weak pass strict constraints AND gravity's implied $m$ lands in the selected GW band window (under CMB $K$):

| GW band | $m$ window | Best $M$ (GeV) | Best $m_g$ | $C$ label | $F_0$ (Hz) |
|---|---:|---:|---:|---|---:|
| CMB | $[-137, -129]$ | $2.93 \times 10^4$ | $-132$ | `SU(3):base/dim` | $1.47 \times 10^{-17}$ |
| PTA | $[-94, -85]$ | $1.59 \times 10^9$ | $-89$ | `SU(3):base/(dim*coxeter)` | $1.43 \times 10^{-8}$ |
| LISA | $[-70, -57]$ | $2.15 \times 10^{12}$ | $-59$ | `SU(3):base/(dim*coxeter)` | $2.65 \times 10^{-2}$ |
| LIGO | $[-46, -38]$ | $5.41 \times 10^{13}$ | $-39$ | `U(1):base` | $4.02 \times 10^2$ |

Caveat: these sweeps find many passing scales, not a single needle. The result is a discrete "compatibility band" rather than a unique prediction unless additional freezes are introduced.

### 11.8 Strict all-forces configurations per GW band

After freezing GW-band gravity targets and using frozen Option-2 phenomenon anchors (EM `em-lyman-alpha`, strong `strong-QCD-200MeV`, weak `weak-W-80.379GeV`):

| GW band | Gravity target | EM $(C,m,K)$ | Strong $(C,m,K)$ | Weak $(C,m,K)$ | Gravity $(C,m)$ | $F_0$ (Hz, $K = 2.725$ K) |
|---|---|---|---|---|---:|---:|
| CMB | `1/alpha_G(GW_CMB)` | (360, +2, $4.52 \times 10^4$ K) | (60, +4, $3.39 \times 10^{11}$ K) | (120, +3, $2.20 \times 10^{14}$ K) | (45, -132) | $1.47 \times 10^{-17}$ |
| PTA | `1/alpha_G(GW_PTA)` | (360, +2, $4.52 \times 10^4$ K) | (60, +4, $3.39 \times 10^{11}$ K) | (120, +3, $2.20 \times 10^{14}$ K) | (15, -89) | $1.43 \times 10^{-8}$ |
| LISA | `1/alpha_G(GW_LISA)` | (360, +2, $4.52 \times 10^4$ K) | (60, +4, $3.39 \times 10^{11}$ K) | (120, +3, $2.20 \times 10^{14}$ K) | (15, -59) | $2.65 \times 10^{-2}$ |
| LIGO | `1/alpha_G(GW_LIGO)` | (360, +2, $4.52 \times 10^4$ K) | (60, +4, $3.39 \times 10^{11}$ K) | (120, +3, $2.20 \times 10^{14}$ K) | (360, -39) | $4.02 \times 10^2$ |

**Key observation:** under frozen anchors, EM/strong/weak are invariant across GW bands; band selection only changes gravity's $(C,m)$.

Relative errors: EM $\approx +3.44 \times 10^{-3}$, strong $\approx -1.27 \times 10^{-2}$, weak $\approx -7.27 \times 10^{-3}$, gravity $\ll 1\%$ for GW-band types.

### 11.9 Strict gravity fit table

| Gravity target | Meaning | Best $(C,m)$ | Rel. err | $F_0$ ($K = 2.725$ K) | Band match |
|---|---|---:|---:|---:|---|
| `1/alpha_G(p)` | Ordinary matter (proton) | (45, -175) | $-6.07 \times 10^{-3}$ | $1.52 \times 10^{-26}$ Hz | none |
| `1/alpha_G(e)` | Cross-check (electron) | (360, -202) | $+3.58 \times 10^{-2}$ | $3.46 \times 10^{-32}$ Hz | none |
| `1/alpha_G(GW_CMB)` | GW CMB type | (45, -132) | $-3.05 \times 10^{-6}$ | $1.47 \times 10^{-17}$ Hz | CMB |
| `1/alpha_G(GW_PTA)` | GW PTA type | (15, -89) | $+3.97 \times 10^{-6}$ | $1.43 \times 10^{-8}$ Hz | PTA |
| `1/alpha_G(GW_LISA)` | GW LISA type | (15, -59) | $+3.38 \times 10^{-7}$ | $2.65 \times 10^{-2}$ Hz | LISA |
| `1/alpha_G(GW_LIGO)` | GW LIGO type | (360, -39) | $-2.30 \times 10^{-7}$ | $4.01 \times 10^2$ Hz | LIGO |
| `1/alpha_G(mP)` | Planck/quantum type | (120, 10) | $-2.43 \times 10^{-2}$ | very high | none |

Ordinary-matter targets do not land in detector bands under CMB $K$ -- consistent with them being a different "gravity type."

### 11.10 Tension with ordinary-mass targets under CMB K

If gravity's temperature is fixed to CMB ($K \approx 2.725$ K), for the proton inverse-gravity fit $m = -175$, the implied frequency is extremely low ($\sim 10^{-26}$ Hz), far below standard GW bands. This indicates either the "gravity temperature" is not the CMB for that target, or the relevant gravity coupling target is not `1/alpha_G(p)`, or "gravity type" must be treated differently.

### 11.11 Planck-strength coupling

Using the Planck mass $m_P$ as the mass anchor yields $\alpha_G(m_P) \approx 1$. Under strict gauge-$C$, there are 5% hits at positive $m$:

- $C = 45, m = 8 \Rightarrow G \approx 0.958$
- $C = 120, m = 10 \Rightarrow G \approx 0.976$

"Planck-strength gravity" sits naturally on the **positive-$m$** side of the ladder, while macro-weak gravity sits on the **negative-$m$** side.

---

## 12. Unification Diagnostic

### 12.1 Procedure

Choose inputs at $\mu_0 = m_Z$: $a_1^{-1}(\mu_0)$, $a_2^{-1}(\mu_0)$, $a_3^{-1}(\mu_0)$. Run each inverse coupling to $\mu$ across a log-spaced grid using 1-loop beta coefficients (SM or MSSM). Define:

$$
\text{score}(\mu) = \max \text{ pairwise difference among } \{a_1^{-1}(\mu), a_2^{-1}(\mu), a_3^{-1}(\mu)\}
$$

Report the $\mu$ where score is minimized.

### 12.2 Results

| Model | Best $Q$ | Max $\Delta\alpha^{-1}$ |
|---|---:|---:|
| SM baseline (`gut-run`) | $\sim 2.41 \times 10^{14}$ GeV | $\approx 3.65$ |
| SM lattice-quantized (`gut-run-lattice`) | $\sim 4.28 \times 10^{14}$ GeV | $\approx 2.16$ |
| MSSM lattice-quantized | $\sim 4.33 \times 10^{16}$ GeV | $\approx 1.55$ |

- SM 1-loop gives poor convergence (classic textbook result).
- MSSM 1-loop gives dramatically improved convergence near $\sim 10^{16}$ GeV (also classic).
- Lattice-quantized inputs appear to reduce the 1-loop mismatch in inverse-coupling convergence.

This is a diagnostic, not proof of unification.

---

## 13. Scoring, Uncertainties, and Falsifiability

### 13.1 Relative-error thresholds

The primary toy-model pass/fail criterion. Reason: toy models can be "close" at the percent level yet have enormous z-scores when experimental sigmas are tiny.

### 13.2 z-scores

$$
z = \frac{G_{\text{pred}} - G_{\text{target}}}{\sigma_{\text{eff}}}, \qquad \sigma_{\text{eff}} = \sqrt{\sigma_{\text{exp}}^2 + \sigma_{\text{theory}}^2}
$$

$\sigma_{\text{theory}}$ is an optional uncertainty floor used only when explicitly justified. For the Z pole, a dedicated mapping method ($\kappa_Z$) replaced the default sigma floor.

### 13.3 What counts as falsification

The project becomes falsifiable when these are frozen: discrete $C$ menu, integer $m$, fixed target definitions, fixed $F_0$ anchor menus, and a fixed tolerance threshold. Then:

- Once you fit an anchor, **RG-within-band predicts additional scales** with no further tuning.
- If those predictions fail broadly as you expand the OOS set (or tighten tolerance), that is a genuine failure mode.

### 13.4 Concrete falsification directions

- Tighten 2% to 1% on RG-predictive suites and see what survives.
- Add more EM scales (e.g. $m_\tau$, 10 GeV, 200 GeV) using external $\alpha(Q)$ references.
- Add more strong scales and threshold structure (still deterministic; no "fit $\Lambda$").
- Check whether the same discrete anchor $(C,m)$ remains stable under cross-check changes.

---

## 14. Limitations and Failure Modes

- **Scale dependence**: strong/weak couplings are not single numbers without specifying a reference scale and scheme.
- **Integer steps vs smooth running**: treating all scale dependence as integer $\Delta m$ steps fails where RG physics matters. The coherent interpretation is $m$ as a coarse band index with within-band motion from deterministic RG running.
- **Gravity ambiguity**: $\alpha_G$ depends on mass scale; without freezing "gravity type," one can sweep mass anchors and fit many outcomes.
- **Frequency anchoring**: Option 2 requires observed $F_0$. Post-hoc selection or expanding the menu remains an overfitting pathway.
- **Overfitting risk**: relaxing the strict $C$ set or allowing orientation flips rapidly expands the hypothesis space.

**Failure modes to watch:**

- If strict constraints cannot simultaneously fit EM/strong/weak under frozen targets at a tighter threshold (e.g. 1%), the "signal" may not be robust.
- If gravity requires sweeping broad mass ranges to find "some" fit for every band, the framework may be underconstrained without additional physical freezes.

### 14.1 Exploratory scan at 3% tolerance

With the extended `--include` menu, at **3%** tolerance:

- EM inverse and hypercharge (GUT norm.) inverse still fit.
- Strong/weak inverse do **not** fit at 3%.

Tightening from 5% to 3% would currently require a broader (but principled) $C$ menu, different strong/weak target choices, or revisiting the frozen orientation.

---

## 15. Work Remaining / Next Steps

### Done

- Frozen coupling orientation (inverse couplings), EM/strong/weak strict targets, gravity orientation and mass anchors (ordinary-matter + GW-band + Planck), strict $C$ candidates, base selection rule, Option-2 anchor menu, K interpretation, all-forces-per-band results, OOS test suites v1-v7, RG-within-band predictive tests, EW sin2 suites, GUT diagnostic, base-vs-alt-bases scaffold.

### Open

- Tighten tolerance and re-test survival rates (5% to 2% to 1%).
- Add a single "score" per configuration combining coupling fit errors, band constraint satisfaction, K plausibility windows, and penalties for free choices.
- Expand or justify gauge-derived $C$ constructions (Casimirs / reps) only if the allowed set stays short.
- Stress-test robustness under cross-check $F_0$ presets.
- Replace toy electroweak mappings with more standard precision-EW approximations (still frozen, still auditable).
- Add notebook sections for visualization (m vs target coupling, gravity band-implied m windows, sensitivity to target choice).
- Add CSV export for sweep results.

---

## 16. CLI Reference

### Setup

```bash
pip install -r requirements.txt
# or
pip install -e .
```

### Core checks

```bash
python -m physics_test.cli check-example
python -m physics_test.cli calc --C 360 --m 2 --K 300
python -m physics_test.cli fits --m 2
python -m physics_test.cli solve-K --m 2 --F0 1e13
```

### List targets and presets

```bash
python -m physics_test.cli list-targets
python -m physics_test.cli list-frequency-presets
python -m physics_test.cli list-gravity-bands
python -m physics_test.cli list-gauge-Cs
python -m physics_test.cli list-norm-families
```

### Strict gauge-derived C scans

```bash
python -m physics_test.cli scan-gauge-Cs --target "1/alpha" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_s_1loop_from_mZ(mH)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha2(alpha(mZ),sin2_on_shell)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha1_GUT(alpha(mZ),sin2_on_shell)" --max-rel-err 0.05
```

### Base-vs-alt-bases experiment

```bash
python -m physics_test.cli base-vs-alt-bases --tol 0.05 --max-nCs 10
```

### Out-of-sample reports (frozen test suites)

```bash
python -m physics_test.cli oos-report --suite v1 --max-rel-err 0.02
python -m physics_test.cli oos-report --suite v2 --max-rel-err 0.02
python -m physics_test.cli oos-report --suite v3 --max-rel-err 0.02
python -m physics_test.cli oos-report --suite v4 --max-rel-err 0.02
python -m physics_test.cli oos-report --suite v5 --max-rel-err 0.05
python -m physics_test.cli oos-report --suite v6 --max-rel-err 0.05
python -m physics_test.cli oos-report --suite v7 --max-rel-err 0.05
```

### Predictive OOS (fit one C per force, hold fixed)

```bash
python -m physics_test.cli oos-predictive --suite v1 --max-rel-err 0.02
python -m physics_test.cli oos-predictive --suite v1 --norm-family inv_C2_fund --max-rel-err 0.02
```

### RG-within-band predictive OOS

```bash
python -m physics_test.cli oos-predictive-rg --suite v1 --max-rel-err 0.02
python -m physics_test.cli oos-predictive-rg --suite v2 --max-rel-err 0.02
python -m physics_test.cli oos-predictive-rg --suite v3 --max-rel-err 0.02
```

### EW mixing and sin2thetaW

```bash
python -m physics_test.cli oos-ew-mix --max-rel-err 0.02
python -m physics_test.cli ew-sin2 --model sm --scales mW,1TeV
python -m physics_test.cli oos-ew-sin2 --model sm --max-rel-err 0.02
python -m physics_test.cli oos-ew-sin2 --suite ew-independent-v1
python -m physics_test.cli oos-ew-sin2 --suite ew-independent-v2
python -m physics_test.cli oos-ew-sin2 --suite ew-independent-v3
python -m physics_test.cli oos-ew-sin2 --suite ew-exploratory-v1
python -m physics_test.cli oos-ew-sin2 --suite ew-dis-exploratory-v1
python -m physics_test.cli oos-ew-sin2 --suite ew-zpole-exploratory-v1
```

### Step-signal OOS (C-independent)

```bash
python -m physics_test.cli oos-steps --suite v1 --max-ratio-err 0.02
python -m physics_test.cli oos-steps --suite v1 --max-ratio-err 0.05
```

### RG / Lambda_QCD consistency

```bash
python -m physics_test.cli rg-scales
python -m physics_test.cli oos-rg --suite qcd-lambda-v1 --max-rel-err 0.06
python -m physics_test.cli oos-rg --suite qcd-lambda-v2 --max-rel-err 0.05
```

### GUT convergence diagnostic

```bash
python -m physics_test.cli gut-run --model sm --Q-min-GeV 1e2 --Q-max-GeV 1e19 --n 2000
python -m physics_test.cli gut-run --model mssm --Q-min-GeV 1e2 --Q-max-GeV 1e19 --n 2000
python -m physics_test.cli gut-run-lattice --model sm --n 400
python -m physics_test.cli gut-run-lattice --model mssm --n 400
```

### Strict all-forces per GW band (Option-2)

```bash
python -m physics_test.cli pair-forces-gaugeCs \
  --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV \
  --weak-preset weak-W-80.379GeV \
  --gravity-band cmb --gravity-targets "1/alpha_G(GW_CMB)" \
  --max-hits 10 --max-results 5

python -m physics_test.cli pair-forces-gaugeCs \
  --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV \
  --weak-preset weak-W-80.379GeV \
  --gravity-band pta --gravity-targets "1/alpha_G(GW_PTA)" \
  --max-hits 10 --max-results 5

python -m physics_test.cli pair-forces-gaugeCs \
  --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV \
  --weak-preset weak-W-80.379GeV \
  --gravity-band lisa --gravity-targets "1/alpha_G(GW_LISA)" \
  --max-hits 10 --max-results 5

python -m physics_test.cli pair-forces-gaugeCs \
  --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV \
  --weak-preset weak-W-80.379GeV \
  --gravity-band ligo --gravity-targets "1/alpha_G(GW_LIGO)" \
  --max-hits 10 --max-results 5
```

### Gravity scans

```bash
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(p)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(e)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(GW_CMB)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(GW_PTA)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(GW_LISA)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(GW_LIGO)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(mP)" --max-rel-err 0.05
```

### Gravity mass sweeps

```bash
python -m physics_test.cli sweep-quantum-gravity --gravity-band cmb \
  --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12
python -m physics_test.cli sweep-quantum-gravity --gravity-band pta \
  --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12
python -m physics_test.cli sweep-quantum-gravity --gravity-band lisa \
  --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12
python -m physics_test.cli sweep-quantum-gravity --gravity-band ligo \
  --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12
```

### Broad exploratory pairing

```bash
python -m physics_test.cli pair-forces-all
python -m physics_test.cli pair-forces-option2 --help
python -m physics_test.cli scan-all --set rotation-degrees --m-min -6 --m-max 6 --m-step 1 --m-integer
```

### Extended gauge-invariant scans (exploratory)

```bash
python -m physics_test.cli list-gauge-Cs --include \
  "base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter),base/roots,base/positive_roots,base/weyl_order,base/center_order,base/C2_fund,base/T_fund"

python -m physics_test.cli scan-gauge-Cs --target "1/alpha" --include \
  "base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter),base/roots,base/positive_roots,base/weyl_order,base/center_order,base/C2_fund,base/T_fund,base/rank" \
  --max-rel-err 0.03
```

---

## 17. References

- [CODATA recommended values of the fundamental physical constants](https://physics.nist.gov/cuu/Constants/)
- [Review of Particle Physics (PDG)](https://pdg.lbl.gov/)
- [NIST Atomic Spectra Database](https://physics.nist.gov/PhysRefData/ASD/lines_form.html)
- Gravitational-wave detector frequency bands: see documentation for LIGO/Virgo/KAGRA, LISA, and PTA collaborations (band definitions in this repo are order-of-magnitude presets).
