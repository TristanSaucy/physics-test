# Draft manuscript — A topological-harmonic toy framework linking dimensionless couplings to frequency scales

**Status**: draft / exploratory (not peer-reviewed)  
**Repository**: `physics-test` (this repo)  
**Date**: 2025-12-13  

## Abstract

We explore a minimal “harmonic-step” toy framework that relates a dimensionless gauge-like quantity $G$ to physical frequency and temperature scales through two coupled equations:

$$
G=\frac{C}{\phi^m},
\qquad
F_0=\phi^m\,\frac{k_B K}{h},
\qquad
G\,F_0=C\,\frac{k_B K}{h}.
$$

Here $\phi$ is the golden ratio and $m\in\mathbb{Z}$ is an integer “harmonic index.”
$C$ is a discrete “topological constant” candidate, $k_B$ is Boltzmann’s constant,
$h$ is Planck’s constant, $K$ is a temperature scale, and $F_0$ is a frequency scale.
To control overfitting, we define a strict mode that (i) restricts $m$ to integers,
(ii) restricts $C$ to a small set of non-arbitrary values generated from Standard Model
gauge-group invariants from a base value 360, and (iii) freezes the canonical
“gauge strength coordinate” to inverse couplings ($1/\alpha$-type targets).
Under these constraints and a 5% relative-error threshold, we find simultaneous fits
for EM, strong, and weak sector targets using only $C\in\{360,180,120,60,45,15\}$.
We then investigate gravity by combining coupling fits with gravitational-wave band
constraints under a fixed CMB temperature; this motivates a “gravity type”
parameterization in terms of an effective mass scale $M$ in the dimensionless gravitational
coupling $\alpha_G(M)=G_N M^2/(\hbar c)$. Sweeps over $M$ reveal distinct mass-scale
bands compatible with CMB/PTA/LISA/LIGO frequency windows under strict constraints.
With frozen phenomenon-frequency anchors (Option 2), the strict end-to-end pairing yields a
single configuration per GW band: the EM/strong/weak hits are identical across bands
((360, 2), (60, 4), (120, 3)), while gravity’s $(C,m)$ shifts to land in the selected GW window.
We emphasize that this is exploratory and that falsifiability depends on maintaining a small
frozen menu of choices (especially gravity “types” and frequency anchors).

## 1. Motivation

The project goal is not to propose a completed theory, but to test whether a very small set of discrete numbers and an integer ladder can create a coherent bridge between:

- dimensionless couplings (pure numbers), and
- frequency/temperature scales associated with physical phenomena.

The main risk is numerology. The main countermeasure is to freeze the “rules of the game” up front (strict $C$, integer $m$, fixed target definitions and reference scales).

## 2. Model

### 2.1 Definitions

Let $\phi=(1+\sqrt{5})/2\approx 1.6180339887$. Define:

$$
G(C,m)=\frac{C}{\phi^m}, \quad m\in\mathbb{Z}.
$$

Define a coupled frequency scale:

$$
F_0(m,K)=\phi^m\,\frac{k_B K}{h}.
$$

Eliminating $\phi^m$ yields an invariant:

$$
G\,F_0=C\,\frac{k_B K}{h}.
$$

An inverse relation used frequently in “Option 2” analyses is:

$$
K=\frac{F_0 h}{k_B\,\phi^m}.
$$

### 2.2 Interpretation of the harmonic index m

From the frequency equation:

$$
\phi^m=\frac{hF_0}{k_B K}.
$$

So $m$ is the base-$\phi$ logarithm of a physically meaningful ratio: quantum energy scale $hF_0$ relative to thermal scale $k_BK$.

- **$m>0$** corresponds to $hF_0>k_BK$ (phenomena “above thermal scale” at temperature $K$).
- **$m<0$** corresponds to $hF_0<k_BK$ (phenomena “below thermal scale” at temperature $K$).

This ratio-based view supports (as a hypothesis) connecting the sign of $m$ to whether a phenomenon is thermally populated vs requires non-equilibrium driving, which in open systems can correlate with entropy export vs import. We treat this as interpretive, not proven.

## 3. Strict contract (anti-overfitting)

This repo uses a strict mode (see `freezing_the_rules.md`) intended to keep the framework falsifiable.

### 3.1 Frozen: integer stepping

- $m\in\mathbb{Z}$ everywhere (no real-valued $m$ fitting).

### 3.2 Frozen: non-arbitrary candidate $C$ set (“gauge-derived $C$”)

We generate candidates from a base value 360 and simple invariants of the Standard Model gauge groups U(1), SU(2), SU(3), using constructions such as:

- $C=\text{base}$
- $C=\text{base}/\dim$
- $C=\text{base}/(\dim\cdot h)$
- $C=\text{base}/h$

where $\dim$ is the Lie algebra dimension and $h$ is the Coxeter number (for $SU(N)$, $h=N$). De-duplicating yields the strict candidate set:

- $C\in\{360,180,120,60,45,15\}$.

Interpretation: these invariants are **representation-independent integers** attached to the Lie algebra / root system. For the SM factors:

- $SU(2)$: $\dim=3$, $h=h^\vee=2$
- $SU(3)$: $\dim=8$, $h=h^\vee=3$

This is intentionally “toy-simple”: divide a fixed base by 1–2 canonical integers and deduplicate. For a deeper discussion of these invariants (rank, Coxeter/dual Coxeter) and “next tier” candidates (Weyl group order, root counts, Casimirs/Dynkin indices and the representation/normalization pitfalls), see `gauge_invariants.md`.

### 3.3 Frozen: canonical target orientation (inverse couplings)

To avoid doubling the hypothesis space, we freeze the “gauge strength coordinate” to inverse couplings:

- EM strict target: `1/alpha` (low-energy fine-structure inverse)
- Strong strict target: `1/alpha_s_1loop_from_mZ(mH)` (1-loop run from `alpha_s(mZ)` to $m_H=125$ GeV; no free $\Lambda_{\mathrm{QCD}}$ knob)
- Weak strict target: `1/alpha2(alpha(mZ),sin2_on_shell)` (derived using on-shell $\sin^2\theta_W=1-m_W^2/m_Z^2$)

Motivation:

- Gauge kinetic terms are often written with coefficients $\propto 1/g^2$, and $\alpha\propto g^2$, so $1/\alpha$ is a natural variable.
- At 1-loop, $\alpha^{-1}(Q)$ runs linearly with $\ln Q$, which aligns better with “integer step” ladder intuition than $\alpha(Q)$ itself.

### 3.4 Pass/fail threshold

We define a hit if:

- $|\text{relative error}|\le 5\%$.

### 3.5 Frozen: Option 2 phenomenon-frequency anchors ($F_0$ menu)

When using Option 2 (fix $F_0$ and solve for $K$), $F_0$ is treated as an observed phenomenon frequency.
To limit degrees of freedom, this repo freezes a small menu of primary + cross-check anchors for EM/strong/weak
and records the citations/notes in `F0_anchors.md` (e.g., EM `em-lyman-alpha`, strong `strong-QCD-200MeV`, weak
`weak-W-80.379GeV`).

## 4. Results under strict constraints (gauge-derived C, integer m, inverse targets)

### 4.1 EM / strong / weak

Using strict gauge-derived $C$ candidates only, integer $m$, and the inverse targets above, we obtain the following hits (all reproducible from the CLI):

| Sector | Target | Target value | Best (C,m) | G(C,m) | rel. err |
|---|---:|---:|---:|---:|---:|
| EM | `1/alpha` | 137.035999 | (360, 2) | 137.507764 | +0.344% |
| Strong | `1/alpha_s_1loop_from_mZ(mH)` | 8.866605 | (60, 4) | 8.753882 | −1.27% |
| Weak | `1/alpha2(alpha(mZ),sin2_on_shell)` | 28.535657 | (120, 3) | 28.328157 | −0.73% |
| Hypercharge (GUT norm.) | `1/alpha1_GUT(alpha(mZ),sin2)` | 59.021547 | (60, 0) | 60.000000 | +1.66% |

Notes:

- The EM hit is the original motivating example: $360/\phi^2$ lands near $1/\alpha$.
- Strong and weak couplings “run” with energy; in strict mode, we freeze a reference scale at the Z boson mass.
- The hypercharge target uses the standard GUT normalization factor $5/3$.

#### 4.1.1 Target-definition sensitivity (and why we freeze these strict targets)

The strong/weak sector “targets” are not single immutable numbers: both are **scale-dependent**, and the electroweak quantities also depend on **scheme/definition** (especially $\sin^2\theta_W$). This matters because our strict $C$ set is very small, so a few-percent target shift can move a fit from “fail” to “pass” at tighter thresholds.

Two particularly relevant target-definition choices (implemented as target keys in `physics_test.targets`) are:

- **Weak (on-shell mixing angle)**: use the pole-mass definition
  $$
  \sin^2\theta_W^{\mathrm{OS}} \equiv 1-\frac{m_W^2}{m_Z^2},
  \qquad
  \alpha_2 \equiv \frac{\alpha(m_Z)}{\sin^2\theta_W},
  $$
  yielding the inverse target `1/alpha2(alpha(mZ),sin2_on_shell)`.

- **Strong (run to a fixed physical scale)**: run the benchmark $\alpha_s(m_Z)$ to $m_H$ with a fixed 1-loop prescription (no free $\Lambda_{\mathrm{QCD}}$), yielding `1/alpha_s_1loop_from_mZ(mH)`.

Under these definitions, the same $(C,m)$ patterns tighten from “few percent” to “$\sim$1%” without changing the strict $C$ set. We therefore freeze them as the strict strong/weak targets, and treat the older $m_Z$ / proxy targets as legacy comparisons.

| Sector | Refined target key | Best (C,m) | rel. err |
|---|---|---:|---:|
| Strong | `1/alpha_s_1loop_from_mZ(mH)` | (60, 4) | −1.27% |
| Weak | `1/alpha2(alpha(mZ),sin2_on_shell)` | (120, 3) | −0.73% |

### 4.2 What is “mZ”?

`mZ` means $m_Z$, the Z boson mass:

- $m_Z\approx 91.1876\,\text{GeV}$.

It is a standard reference energy scale at which $\alpha_s$ and electroweak parameters are commonly quoted, because it is experimentally clean and because couplings are scale-dependent.

### 4.3 1-loop GUT-style convergence diagnostic (cross-check)

This repo includes a toy 1-loop running test for $\alpha_1^{GUT},\alpha_2,\alpha_3$ in the conventional inverse-coupling coordinates. Using approximate $m_Z$ anchors, it reproduces the classic qualitative result:

- SM 1-loop: poor convergence ($\Delta \alpha^{-1}\sim\mathcal{O}(1)$)
- MSSM 1-loop: dramatically improved convergence near $Q\sim 10^{16}$ GeV

Example outputs (from `physics_test.cli gut-run`):

- **SM** best: $Q\approx 2.41\times 10^{14}$ GeV, max $\Delta\alpha^{-1}\approx 3.65$
- **MSSM** best: $Q\approx 2.09\times 10^{16}$ GeV, max $\Delta\alpha^{-1}\approx 0.052$

This does not validate the topology model, but it supports the “inverse coupling coordinate” choice as a principled convention for unification-style comparisons.

### 4.4 Thermal baseline sanity check (bio-adjacent anchor)

The model implies a universal thermal frequency scale at $m=0$:

$$
f_T=\frac{k_B T}{h}.
$$

At $T=310$ K (human body temperature):

- $f_T\approx 6.46$ THz
- $m=1\Rightarrow F_0=\phi f_T\approx 10.45$ THz

This is a clean anchor because it introduces no new physics—only a direct conversion between temperature and a frequency scale.

## 5. Gravity: inverse coupling targets, bands, and mass-scale sweeps

### 5.1 Dimensionless gravitational coupling

We use the standard dimensionless gravitational coupling at a chosen mass scale $M$:

$$
\alpha_G(M)=\frac{G_N M^2}{\hbar c},
\qquad
\frac{1}{\alpha_G(M)}=\frac{\hbar c}{G_N M^2}.
$$

Unlike gauge couplings, $\alpha_G$ depends on the choice of $M$. This motivates the concept of different “gravity types” (macro gravity, GW-band gravity, quantum gravity) corresponding to different frozen mass anchors.

### 5.2 Electron vs proton (inverse targets)

For inverse targets:

- `1/alpha_G(p)` $\approx 1.69\times 10^{38}$ (proton mass anchor)
- `1/alpha_G(e)` $\approx 5.71\times 10^{44}$ (electron mass anchor)

**Key identity (mass ratio $\Rightarrow$ coupling ratio):** since

$$
\alpha_G(M)=\frac{G_N M^2}{\hbar c},
\qquad
\frac{1}{\alpha_G(M)}=\frac{\hbar c}{G_N M^2},
$$

the inverse-coupling ratio is purely a squared mass ratio:

$$
\frac{1/\alpha_G(e)}{1/\alpha_G(p)}=\left(\frac{m_p}{m_e}\right)^2\approx 3.37\times 10^6,
\qquad
\log_\phi\left(\left(\frac{m_p}{m_e}\right)^2\right)\approx 31.24.
$$

So if $C$ were held fixed, moving from proton to electron as the mass anchor would shift the implied harmonic index by $\sim 31$ steps.

**How this predicts an $m$-shift in the topology model (including the discrete $C$ correction):** if we match a target value $G_t$ with

$$
G_t\approx \frac{C}{\phi^m},
$$

then (formally) the best-fit $m$ is:

$$
m \approx \log_\phi\!\left(\frac{C}{G_t}\right)=\log_\phi(C)-\log_\phi(G_t).
$$

For two targets $G_{t,1},G_{t,2}$ and two allowed constants $C_1,C_2$, the predicted shift is:

$$
\Delta m = m_2-m_1 \approx \log_\phi\!\left(\frac{C_2}{C_1}\right)-\log_\phi\!\left(\frac{G_{t,2}}{G_{t,1}}\right).
$$

For electron vs proton inverse gravity targets, $G_{t,2}/G_{t,1}=(m_p/m_e)^2$. In the strict gauge-derived $C$ set, the best proton hit uses $C_1=45$ while the electron hit uses $C_2=360$, so $C_2/C_1=8$ and:

$$
\Delta m \approx \log_\phi(8)-\log_\phi\!\left(\left(\frac{m_p}{m_e}\right)^2\right)
\approx 4.32-31.24
\approx -26.9.
$$

This matches the observed strict hits ($m_p$ around $-175$ and $m_e$ around $-202$), i.e. a shift of $-27$ integer steps, within rounding/discreteness.

Under strict gauge-derived $C$ and 5%:

- `1/alpha_G(p)` has hits at $m\approx -175$ (e.g. $C=45$) and $m\approx -173$ (e.g. $C=120$).
- `1/alpha_G(e)` has a hit at $m\approx -202$ (e.g. $C=360$).

**Freeze implied by this structure (falsifiability):** we freeze “ordinary matter gravity” to the proton anchor `1/alpha_G(p)` and require the electron `1/alpha_G(e)` to appear as a cross-check with the mass-ratio–predicted $m$-shift (up to allowed discrete $C$ ratios), rather than allowing “proton vs electron” to be chosen ad hoc.

### 5.3 Gravity + CMB temperature: a tension with ordinary-mass targets

If we fix gravity’s temperature to the CMB ($K\approx 2.725$ K), then the frequency equation predicts:

$$
F_0=\phi^m\frac{k_B K}{h}.
$$

For the proton inverse-gravity fit $m=-175$, the implied frequency is extremely low ($\sim 10^{-26}$ Hz), far below standard gravitational wave bands. This indicates that:

- either the “gravity temperature” $K$ is not the CMB for that target,
- or the relevant gravity coupling target is not `1/alpha_G(p)`,
- or “gravity type” must be treated differently.

### 5.4 Band-constrained mass-scale solutions (inverse gravity “types”)

This repo includes a sweep that searches for mass scales $M$ such that **all** of the following hold simultaneously:

- EM/strong/weak fits pass strict gauge-$C$ + integer $m$ + 5% tolerance
- gravity target is `1/alpha_G(M)` (inverse, consistent with frozen orientation)
- gravity uses $K=$ CMB
- gravity $m$ lies inside a chosen GW-band frequency window (CMB/PTA/LISA/LIGO)

Representative best hits found in wide sweeps:

| GW band | $m$ window (CMB K) | Best scale $M$ (GeV) | Best $m_g$ | Example $C$ label | $F_0$ (Hz) |
|---|---:|---:|---:|---|---:|
| CMB/primordial | $[-137,-129]$ | $2.93\times10^4$ | $-132$ | `SU(3):base/dim` | $1.47\times10^{-17}$ |
| PTA | $[-94,-85]$ | $1.59\times10^9$ | $-89$ | `SU(3):base/(dim*coxeter)` | $1.43\times10^{-8}$ |
| LISA | $[-70,-57]$ | $2.15\times10^{12}$ | $-59$ | `SU(3):base/(dim*coxeter)` | $2.65\times10^{-2}$ |
| LIGO | $[-46,-38]$ | $5.41\times10^{13}$ | $-39$ | `U(1):base` | $4.02\times10^{2}$ |

In this repo we now **freeze** these as named “gravity types” (targets `1/alpha_G(GW_CMB)`, `1/alpha_G(GW_PTA)`, `1/alpha_G(GW_LISA)`, `1/alpha_G(GW_LIGO)`) rather than continuing to sweep mass scales.

#### 5.4.1 Observational platforms and locations (what “bands” mean)

The frequency “bands” above are **order-of-magnitude sensitivity windows** of different observational techniques. They are not all the same kind of measurement:

- **LIGO band (10–1000 Hz)**: **ground-based laser interferometers** measuring strain on Earth.
  - LIGO: Hanford, Washington (USA) and Livingston, Louisiana (USA)
  - Virgo: Cascina (near Pisa), Italy
  - KAGRA: Kamioka, Japan (underground)
- **LISA band (0.1 mHz–0.1 Hz)**: **space-based laser interferometer** (planned; not operating yet).
  - Mission concept: three spacecraft in a heliocentric orbit trailing Earth, forming a triangular interferometer.
- **PTA band (nHz)**: **pulsar timing arrays**, i.e. radio telescopes on Earth timing millisecond pulsars and searching for correlated timing residuals.
  - Examples: NANOGrav, EPTA, PPTA; combined efforts under IPTA.
- **CMB/primordial band (~10⁻¹⁸–10⁻¹⁶ Hz)**: **indirect constraints**, not a direct interferometric detection.
  - These come from measurements of CMB temperature/polarization anisotropies (e.g., satellites at L2 and ground/balloon experiments), which constrain primordial tensor modes on horizon scales.

Connection to “celestial ordering” intuition: for a binary source, the GW frequency is typically $\sim 2\times$ the orbital frequency. Larger, slower systems (wide separations, higher total masses) sit at lower frequencies; smaller/compact systems sit at higher frequencies. However, detectability depends strongly on **compactness** and emitted strain—solar-system bodies radiate gravitational waves far too weak to detect with current techniques despite having orbital frequencies that can fall in these ranges.

Important caveat: these sweeps typically find **many** passing scales, not a single needle, so the result is better interpreted as a discrete “compatibility band” rather than a unique prediction unless additional freezes are introduced.

#### 5.4.2 Strict end-to-end “all-forces” configurations per GW band (frozen $F_0$ anchors)

After freezing GW-band gravity targets (see §5.4 and `gravity_types_report.md`), we can run a strict all-forces pairing search (`pair-forces-gaugeCs`) that simultaneously:

- matches EM/strong/weak inverse targets within 5% using strict gauge-derived $C$ and integer $m$,
- uses frozen phenomenon frequencies (Option 2) to solve for implied $K$ in EM/strong/weak, and
- requires gravity (with $K=2.725$ K) to land in the chosen GW band window while matching the chosen frozen gravity type target.

Using the primary frozen frequency anchors (EM `em-lyman-alpha`, strong `strong-QCD-200MeV`, weak `weak-W-80.379GeV`),
the best strict configuration per GW band is:

| GW band | Gravity target | EM (C,m,K) | Strong (C,m,K) | Weak (C,m,K) | Gravity (C,m) | $F_0$ (Hz, $K=2.725$ K) |
|---|---|---|---|---|---:|---:|
| CMB | `1/alpha_G(GW_CMB)` | (360, +2, $4.52\times 10^4$ K) | (60, +4, $3.39\times 10^{11}$ K) | (120, +3, $2.20\times 10^{14}$ K) | (45, −132) | $1.47\times 10^{-17}$ |
| PTA | `1/alpha_G(GW_PTA)` | (360, +2, $4.52\times 10^4$ K) | (60, +4, $3.39\times 10^{11}$ K) | (120, +3, $2.20\times 10^{14}$ K) | (15, −89) | $1.43\times 10^{-8}$ |
| LISA | `1/alpha_G(GW_LISA)` | (360, +2, $4.52\times 10^4$ K) | (60, +4, $3.39\times 10^{11}$ K) | (120, +3, $2.20\times 10^{14}$ K) | (15, −59) | $2.65\times 10^{-2}$ |
| LIGO | `1/alpha_G(GW_LIGO)` | (360, +2, $4.52\times 10^4$ K) | (60, +4, $3.39\times 10^{11}$ K) | (120, +3, $2.20\times 10^{14}$ K) | (360, −39) | $4.02\times 10^{2}$ |

**Key observation:** under these frozen anchors, EM/strong/weak are invariant across GW bands; band selection only changes gravity’s $(C,m)$.
This is a stronger “coherence” check than fitting each sector independently.

### 5.5 Planck-strength coupling and positive m

Using the Planck mass $m_P$ as the mass anchor yields $\alpha_G(m_P)\approx 1$, so both `alpha_G(mP)` and `1/alpha_G(mP)` are $\approx 1$. Under strict gauge-$C$, there are 5% hits at positive $m$, e.g.:

- $C=45, m=8 \Rightarrow G\approx 0.958$
- $C=120, m=10 \Rightarrow G\approx 0.976$

This supports the idea that “Planck-strength gravity” sits naturally on the **positive-$m$** side of the ladder, while macro-weak gravity sits on the **negative-$m$** side—though these correspond to very different implied frequencies at a fixed $K$.

## 6. Discussion

### 6.1 Why inverse couplings are the cleanest frozen choice

Allowing both $\alpha$ and $1/\alpha$ after the fact doubles the target menu and raises the risk of post-hoc fitting. Inverse couplings also have principled advantages:

- They appear naturally in gauge kinetic terms ($\propto 1/g^2$).
- They run linearly at 1-loop ($\alpha^{-1}(Q)$ is approximately affine in $\ln Q$).
- Unification diagnostics are conventionally expressed as line intersections in $\alpha_i^{-1}$.

### 6.2 A grounded view of “octaves,” sign of m, and entropy (hypothesis)

The identity $\phi^m=hF_0/(k_BK)$ suggests a concrete meaning for the sign of $m$:

- $m>0$: quantum scale exceeds thermal scale at $K$
- $m<0$: quantum scale lies below thermal scale at $K$

If a system sustains high-$m$ dynamics in a warm environment, it typically requires non-equilibrium structure and dissipation. This can be consistent with the intuition that “positive $m$” phenomena correlate with entropy export (maintained order), while “negative $m$” phenomena correlate with entropy import/bath-like behavior. This remains interpretive and should be treated as a hypothesis to test against concrete phenomenon choices.

### 6.3 Hypothesis: m as a discrete strength coordinate (why strong > weak > EM > gravity)

In strict mode we map the topology variable $G$ to **inverse couplings** ($G\sim 1/\alpha$). Under that convention, the model implies:

$$
G=\frac{C}{\phi^m}
\quad\Rightarrow\quad
\alpha \sim \frac{\phi^m}{C}.
$$

Taking a base-$\phi$ logarithm makes the role of $m$ explicit:

$$
\log_\phi(\alpha) \sim m-\log_\phi(C).
$$

So, **holding $C$ fixed**, each +1 step in $m$ multiplies $\alpha$ by a factor of $\phi\approx 1.618$ (and divides the inverse coupling by $\phi$). In that sense, $m$ can be interpreted as a **discrete log-strength index**.

In the strict best fits for the gauge forces (all at 5% tolerance, with frozen inverse targets at the $m_Z$ benchmark where applicable), the ordering is:

- strong: $(C,m)=(60,4)$
- weak: $(C,m)=(120,3)$
- EM: $(C,m)=(360,2)$

This is consistent with the usual hierarchy $\alpha_s>\alpha_w>\alpha$. Here, the **relative inverse-strengths** are explained by a small number of discrete factors coming from *both* the $m$ steps and the discrete $C$ ratios. For adjacent pairs:

$$
\frac{G_{\text{weak}}}{G_{\text{strong}}}
=\frac{C_w}{C_s}\,\phi^{m_s-m_w}
\approx 2\phi,
\qquad
\frac{G_{\text{EM}}}{G_{\text{weak}}}
=\frac{C_{em}}{C_w}\,\phi^{m_w-m_{em}}
\approx 3\phi,
$$

so a one-step $m$ shift plus a small integer $C$ ratio produces an $\mathcal{O}(1\text{–}10)$ hierarchy across the gauge forces.

For **ordinary-matter gravity**, using the canonical inverse target at the proton mass anchor gives a strict fit around $(C,m)\approx(45,-175)$ (see §5.2). Relative to the EM fit $(360,2)$, the inverse-coupling hierarchy is dominated by the large $m$ gap:

$$
\frac{G_{\text{grav}}}{G_{\text{EM}}}
=\frac{C_g}{C_{em}}\,\phi^{m_{em}-m_g}
\approx \frac{1}{8}\,\phi^{177}\sim 10^{36},
$$

which matches the observed enormous separation between $1/\alpha_G(p)$ and $1/\alpha$.

Interpretation (hypothesis): if $m$ is a real discrete “harmonic level,” then the observed hierarchy of interaction strengths could be reframed as: the gauge forces occupy **nearby positive-$m$ micro levels** (few steps apart), while ordinary-matter gravity occupies a **far negative-$m$ macro level**, producing an exponential weakness in the inverse-coupling coordinate. This is only meaningful under the strict contract; if $C$ or target definitions are allowed to drift freely, $m$ can be made to absorb almost anything.

### 6.4 Limitations and failure modes

- **Scale dependence**: strong/weak couplings are not single numbers without specifying a reference scale and scheme.
- **Gravity ambiguity**: $\alpha_G$ depends on mass scale; without freezing “gravity type” one can sweep mass anchors and fit many outcomes.
- **Frequency anchoring**: Option 2 requires observed $F_0$. This repo freezes a small anchor menu (`F0_anchors.md`), but post-hoc selection or expanding the menu remains an overfitting pathway.
- **Overfitting risk**: relaxing the strict $C$ set or allowing orientation flips rapidly expands the hypothesis space.

Failure modes worth watching:

- If strict constraints cannot simultaneously fit EM/strong/weak under frozen targets at a tighter threshold (e.g. 1%), the “signal” may not be robust.
- If gravity requires sweeping broad mass ranges to find “some” fit for every band, the framework may be underconstrained without additional physical freezes.

## 7. Next steps to increase falsifiability

1) **Done in this repo**: freeze a **small table of phenomenon $F_0$ anchors** (with notes/citations) for Option 2 (see `F0_anchors.md`).
2) **Done in this repo**: for gravity, freeze:
   - orientation (inverse),
   - temperature choice for GW-band gravity types (CMB),
   - a small set of mass/type anchors (ordinary matter + GW-band types; see §5.2–§5.4 and `gravity_types_report.md`).
3) Tighten tolerance and re-test survival rates (e.g., 5% → 2% → 1%).
4) Expand or justify $C$ derivations with additional group-invariant constructions only if they remain short and principled.

## 8. Reproducibility (CLI commands)

Core example:

```bash
python -m physics_test.cli check-example
```

Strict gauge-derived $C$ list:

```bash
python -m physics_test.cli list-gauge-Cs
```

Strict scans for the frozen inverse targets:

```bash
python -m physics_test.cli scan-gauge-Cs --target "1/alpha" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_s_1loop_from_mZ(mH)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha2(alpha(mZ),sin2_on_shell)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha1_GUT(alpha(mZ),sin2)" --max-rel-err 0.05
```

Strict all-forces configurations per GW band (Option 2; frozen $F_0$ presets):

```bash
python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band cmb  --gravity-targets "1/alpha_G(GW_CMB)"  --max-hits 10 --max-results 5
python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band pta  --gravity-targets "1/alpha_G(GW_PTA)"  --max-hits 10 --max-results 5
python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band lisa --gravity-targets "1/alpha_G(GW_LISA)" --max-hits 10 --max-results 5
python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band ligo --gravity-targets "1/alpha_G(GW_LIGO)" --max-hits 10 --max-results 5
```

Out-of-sample report (frozen test suites v1/v2):

```bash
python -m physics_test.cli oos-report --suite v1 --max-rel-err 0.02
python -m physics_test.cli oos-report --suite v2 --max-rel-err 0.02
```

Gravity band sweeps (inverse gravity targets, CMB K):

```bash
python -m physics_test.cli sweep-quantum-gravity --gravity-band cmb  --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12
python -m physics_test.cli sweep-quantum-gravity --gravity-band pta  --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12
python -m physics_test.cli sweep-quantum-gravity --gravity-band lisa --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12
python -m physics_test.cli sweep-quantum-gravity --gravity-band ligo --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12
```

GUT running diagnostic:

```bash
python -m physics_test.cli gut-run --model sm   --Q-min-GeV 1e2 --Q-max-GeV 1e19 --n 2000
python -m physics_test.cli gut-run --model mssm --Q-min-GeV 1e2 --Q-max-GeV 1e19 --n 2000
```

## References (informal starting points)

- [CODATA recommended values of the fundamental physical constants](https://physics.nist.gov/cuu/Constants/)
- [Review of Particle Physics (PDG)](https://pdg.lbl.gov/)
- Gravitational-wave detector frequency bands: see documentation for LIGO/Virgo/KAGRA, LISA, and PTA collaborations (band definitions in this repo are order-of-magnitude presets).
