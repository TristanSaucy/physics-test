# Concise model note (φ‑lattice + RG‑within‑band) — status: exploratory

This document is a **concise, falsifiability‑oriented** summary of the current state of the `physics-test` project:

- the **core formulas** and definitions,
- what is **frozen vs exploratory** (anti‑overfitting “contract”),
- how we handle **renormalization / running** (QCD + QED),
- the key **findings** and what is and is not predictive,
- and why this is exploratory (how it might be used to hunt new phenomena) without claiming it’s a finished theory.

For the longer narrative + details, see `paper.md` and `freezing_the_rules.md`.

---

## 1) Core model (frozen)

Define the golden ratio:

\[
\phi \equiv \frac{1+\sqrt{5}}{2}.
\]

Define a dimensionless “lattice” quantity:

\[
G(C,m) \equiv \frac{C}{\phi^m}, \qquad m\in\mathbb{Z}.
\]

Define a coupled frequency scale:

\[
F_0(m,K) \equiv \phi^m\,\frac{k_B K}{h}.
\]

Eliminate \(\phi^m\) to get the invariant:

\[
G\,F_0 = C\,\frac{k_B K}{h}.
\]

Option‑2 (“phenomenon‑first”) inversion used for frequency anchoring:

\[
K = \frac{F_0\,h}{k_B\,\phi^m}.
\]

---

## 2) Interpretation of parameters (working hypothesis)

- **\(m\)**: an integer **band / harmonic index**. It is *not* assumed to encode smooth RG running by itself.
- **\(C\)**: a discrete “topological constant candidate” restricted (in strict mode) to a small, justified menu.
- **\(G\)**: treated as a **dimensionless coupling coordinate**. In strict claims we use **inverse couplings** as targets (see below).
- **\(F_0, K\)**: used to connect the lattice band index to physical frequency/temperature scales (Option‑2 is preferred when \(F_0\) is genuinely observed).

---

## 3) Strict contract (anti‑overfitting / falsifiability)

Strict mode freezes:

- **Integer stepping**: \(m \in \mathbb{Z}\) only.
- **Target orientation**: strict claims use **inverse couplings** (`1/alpha`‑style).
- **Gauge‑derived \(C\) menu**: generate \(C\) from SM gauge group invariants using base \(=360\), then de‑duplicate.
  - Current strict candidate set:
    \[
    C\in\{360,180,120,60,45,15\}.
    \]
- **Reference scales / schemes** for running‑sensitive targets are frozen (examples below).
- **Pass/fail**: we typically report at **2%** for OOS tests and keep **5%** as a broader exploratory scan threshold.

This is what keeps the project from devolving into unconstrained numerology: the hypothesis space is deliberately small.

### 3.1 Measurement registry (toward rigor)

To make the framework auditable and uncertainty-aware, key external inputs are stored in:

- `data/targets.json`

Each entry includes (where known): **value**, **1σ uncertainty**, optional **reference scale** \(Q\), a short **scheme/notes** string, and a **citation hint**.
The code reads this via `physics_test/target_registry.py` and attaches metadata to `TargetConstant` objects (see `physics_test/targets.py`).

---

## 4) Target definitions (what we fit \(G\) to)

The framework compares \(G\) to dimensionless target constants from `physics_test/targets.py`.

### 4.1 EM (electromagnetism)

- **Strict EM anchor** (low energy / vacuum):
  - target: `1/alpha`
- **OOS EM cross‑check**:
  - target: `1/alpha(mZ)` (effective coupling at the Z pole; vacuum polarization)

### 4.2 Strong (QCD)

Because \(\alpha_s(Q)\) runs strongly, targets must specify scale and prescription.

- **Strict strong anchor**:
  - target: `1/alpha_s_1loop_from_mZ(mH)`
  - definition: take benchmark \(\alpha_s(m_Z)\) and run to \(m_H=125\) GeV using 1‑loop QCD running with fixed \(n_f=5\) (no free \(\Lambda_{\rm QCD}\) knob).

- **OOS strong cross‑checks**:
  - `1/alpha_s_1loop_from_mZ(mW)`, `(...mt)`, `(...1TeV)`, `(...10TeV)`
  - `nf56` and `2loop` variants (still deterministic, still no free parameters; see §5).

### 4.3 Weak / electroweak (EW)

Weak targets are scheme‑sensitive (especially \(\sin^2\theta_W\)).

- **Strict weak anchor**:
  - target: `1/alpha2(alpha(mZ),sin2_on_shell)`
  - definition: use on‑shell mixing angle:
    \[
    \sin^2\theta_W^{\rm OS} \equiv 1-\frac{m_W^2}{m_Z^2},
    \qquad
    \alpha_2 \equiv \frac{\alpha(m_Z)}{\sin^2\theta_W^{\rm OS}}.
    \]

EW “running within band” is not implemented yet (see §9).

---

## 5) Renormalization / “within‑band” running (the key upgrade)

Pure integer steps (varying only \(\Delta m\)) fail exactly where physics says they should: **couplings run continuously with \(\ln Q\)**.

So we introduced a stricter, more physical predictive mode:

- **Fit the anchor once** on the φ‑lattice (choose best \(C\) and integer \(m\)).
- **Hold \(C\) and the anchor \(m\) fixed** (band selection is discrete).
- Use a **deterministic RG prescription** to predict the coupling at new scales **without re‑fitting \(m\)**.

CLI:

```bash
python -m physics_test.cli oos-predictive-rg --suite v1 --max-rel-err 0.02
```

### 5.1 QCD running used (strong sector)

**1‑loop (inverse coupling runs linearly in \(\ln Q\))**:

\[
\alpha_s^{-1}(Q) = \alpha_s^{-1}(Q_0) + \frac{b_0}{2\pi}\ln\!\left(\frac{Q}{Q_0}\right),
\qquad
b_0 = 11-\frac{2}{3}n_f.
\]

**2‑loop** is implemented by integrating:

\[
\frac{d\alpha_s}{d\ln Q}
:=
-\frac{\beta_0}{2\pi}\alpha_s^2
-\frac{\beta_1}{4\pi^2}\alpha_s^3,
\quad
\beta_0 = 11-\frac{2}{3}n_f,\;
\beta_1 = 102-\frac{38}{3}n_f.
\]

Optional minimal threshold model: switch \(n_f=5\to 6\) at \(Q\simeq m_t\) (no matching).

### 5.2 QED running used (EM sector)

We implemented a deterministic **1‑loop QED threshold model** for the inverse coupling:

\[
\frac{d(\alpha^{-1})}{d\ln Q}
:=
-\frac{2}{3\pi}\sum_f N_c\,Q_f^2,
\]

integrated with **sharp thresholds** at fermion masses. In code:

- start from the lattice‑fit \(1/\alpha\) anchor,
- choose \(Q_0=m_e\) as a fixed reference scale for the toy model,
- run up to \(Q=m_Z\) by turning on each fermion contribution above its mass threshold.

This is **not precision electroweak vacuum polarization** (hadronic contributions are subtle), but it is deterministic and captures the correct direction/scale of running.

---

## 6) What we tested (predictability vs flexibility)

We distinguish:

- **OOS report** (`oos-report`): lets \(C\) vary per target within the strict menu (good for “does anything in the strict menu hit?”).
- **Predictive OOS** (`oos-predictive`): fits one best \(C\) per force on an anchor, then holds \(C\) fixed but still allows integer \(m\) to vary per target.
- **RG‑predictive OOS** (`oos-predictive-rg`): fits the anchor once, then predicts other scales via RG **with no re‑fitting of \(m\)**. This is the closest thing we have to “band + within‑band running.”

---

## 7) Key findings (high signal)

### 7.1 Discrete anchors (strict gauge‑derived \(C\))

Under strict gauge‑derived \(C\) and integer \(m\), the best anchor hits are:

- **EM**: `1/alpha` at \((C,m)=(360,2)\) (≈ +0.34%)
- **Strong**: `1/alpha_s_1loop_from_mZ(mH)` at \((C,m)=(60,4)\) (≈ −1.27%)
- **Weak**: `1/alpha2(alpha(mZ),sin2_on_shell)` at \((C,m)=(120,3)\) (≈ −0.73%)
- **Hypercharge (GUT‑normalized)**: `1/alpha1_GUT(alpha(mZ),sin2_on_shell)` at \((C,m)=(60,0)\) (≈ +0.58%)

### 7.2 The “integer‑steps only” failure mode

When you try to explain scale‑dependence using only integer \(\Delta m\) steps, you predictably miss running‑sensitive targets (strong near \(m_W\)/\(m_t\)/1 TeV, EM at \(m_Z\)).

### 7.3 RG‑within‑band resolves those misses (no new fit knobs)

Using `oos-predictive-rg`:

- **Strong** running cross‑checks that missed before become **passes at 2%** (typical errors ~1–2% across v2/v3/v4 strong running keys).
- **EM** OOS cross‑check `1/alpha → 1/alpha(mZ)` becomes a **pass at 2%** under deterministic QED running (either a PDG-style Δα(mZ²) relation or a simple 1-loop threshold runner; both are available as `--runner` options).
- **Weak** within‑band running of `alpha2^{-1}(Q)` becomes **passes at 2%** across `mW`, `mH`, `1TeV`, `10TeV` (suite `v2`).
- **Hypercharge** within‑band running of `alpha1_GUT^{-1}(Q)` becomes **passes at 2%** across `mW`, `mH`, `1TeV`, `10TeV` (suite `v3`).
- A derived **EW mixing** consistency check `oos-ew-mix` passes at **2%** across `mW`, `mH`, `1TeV`, `10TeV` (typical offset ~1%).
- An **external EW cross-check** using muon-decay `G_F` (suite `v5`) lands on the same strict weak lattice point \((C,m)=(120,3)\) but is **off by ~3.9%** at 1-loop/2% tolerance. This is expected to be dominated by **electroweak radiative corrections** (Δr) rather than φ‑stepping.

Interpretation: \(m\) behaves like a **coarse band index**, while **RG flow supplies within‑band motion**.

### 7.4 Δr shows a sharp lattice alignment (exploratory)

If you compute the on-shell electroweak radiative correction parameter \(\Delta r\) from \(\alpha(0),G_F,m_W,m_Z\),
its inverse lands very close to the strict weak lattice point:

- `1/delta_r(on-shell;alpha0,GF,mW,mZ)` ≈ 28.244
- best strict hit: \((C,m)=(120,3)\Rightarrow 120/\phi^3 \approx 28.328\) (≈ +0.30%).

This does **not** mean the φ‑lattice “explains Δr”; it is a notable numerical coincidence worth tracking as we add more independent EW inputs.

### 7.5 GUT‑style convergence improves with lattice‑quantized inputs (exploratory)

The command `gut-run-lattice` initializes \(\alpha_1^{-1},\alpha_2^{-1},\alpha_3^{-1}\) at \(m_Z\) using the φ‑lattice anchors (and runs the strong anchor from \(m_H\to m_Z\) deterministically), then scans for best 1‑loop convergence at high scales.

In current runs:

- Baseline `gut-run --model sm`: best score \(\approx 3.67\) at \(Q\sim 2.46\times 10^{14}\) GeV.
- `gut-run-lattice --model sm`: best score \(\approx 2.16\) at \(Q\sim 4.28\times 10^{14}\) GeV.
- `gut-run-lattice --model mssm`: best score \(\approx 1.55\) at \(Q\sim 4.33\times 10^{16}\) GeV.

This is not “proof of unification,” but it is a nontrivial diagnostic: lattice‑quantized inputs appear to reduce the 1‑loop mismatch in inverse‑coupling convergence.

### 7.6 Radiative-correction pieces show additional lattice structure (exploratory)

Suite `v6` probes inverse vacuum-polarization pieces and an inferred \(\Delta r\):

- `1/delta_alpha_total(mZ2)` best strict hit: \((C,m)=(45,2)\) (≈ +1.53%).
- `1/delta_r(on-shell;alpha0,GF,mW,mZ)` best strict hit: \((C,m)=(120,3)\) (≈ +0.30%).

At 2% tolerance, 2/4 pass; at 5% tolerance, 4/4 pass. These are **diagnostics** (not strict couplings), but they strengthen the pattern that “RG/log physics” quantities often land near the same small set of strict lattice points.

Suite `v7` extends this with a standard leading EW loop term:

- `1/delta_rho_top(GF,mt)` (inverse leading top-loop \(\Delta\rho\)) lands **within ~4%** of strict lattice points (e.g. \((C,m)=(180,1)\) or \((15,-4)\)), so it passes at 5% but not 2%.

---

## 8) Predictability and falsifiability (what would count as “wrong”)

The project becomes falsifiable when the following are frozen:

- discrete \(C\) menu,
- integer \(m\),
- fixed target definitions (including reference scales and schemes),
- fixed \(F_0\) anchor menus for Option‑2 (where used),
- and a fixed tolerance threshold.

Then:

- Once you fit an anchor, **RG‑within‑band predicts additional scales** with no further tuning.
- If those predictions fail broadly as you expand the OOS set (or tighten tolerance), that is a genuine failure mode of the framework.

### 8.1 Uncertainty-aware reporting (z-scores / χ²)

CLI reports now include:

- **z-score** \(z \equiv (G_{\text{pred}}-G_{\text{target}})/\sigma\) when a target has a registry-provided \(\sigma\),
- and simple **χ² sums** over sigma-annotated targets in a suite.

Important: these \(\sigma\) values are **measurement uncertainties only**. They do **not** include “model discrepancy” (the fact that the lattice is a toy ansatz).
So for ultra-precise quantities (e.g., \(\alpha\)), z-scores can be enormous even when the model is “close” at the percent level.
We therefore keep **relative-error thresholds** as the primary “toy-model” pass/fail criterion, while z/χ² is used to keep inputs honest and make residual structure obvious.

Concrete falsification directions:

- tighten 2% → 1% on the RG‑predictive suites and see what survives,
- add more EM scales (e.g., \(m_\tau\), 10 GeV, 200 GeV) using external \(\alpha(Q)\) references,
- add more strong scales and threshold structure (still deterministic; no “fit Λ”),
- check whether the same discrete anchor \((C,m)\) remains stable when you change only *cross‑check* choices (not free alternatives).

---

## 9) What’s still exploratory (and where “new phenomena” could hide)

This repo is exploratory in two ways:

1) **Band structure hypothesis**: integer \(m\) could correspond to a real discretization of “levels/bands” in how phenomena choose scales (resonance / harmonics / entropy‑flow heuristics). This is not established physics; it’s a hypothesis being stress‑tested under strict constraints.

2) **Search workflow**: once you freeze the model, you can use it to propose:
   - candidate **\(m\) windows** for phenomena constrained by frequency/temperature anchors (Option‑2),
   - candidate **cross‑check scales** where a coupling should land “near” a band boundary,
   - candidate **gravity “types”** tied to GW bands under fixed CMB \(K\) (see `gravity_types_report.md`).

Important: until the model makes a prediction that is *not* already “baked in” by a known RG formula (and survives strict OOS tests), this remains hypothesis‑generation, not discovery.

### 9.1 Next lever: EW/QED+EW unification layer

- Add a deterministic within‑band running mode for \(\alpha_1(Q)\), \(\alpha_2(Q)\), and \(\sin^2\theta_W(Q)\) under frozen scheme choices, then re‑run the EW OOS suite in the same spirit as strong+EM.
- **Status**: EW within‑band running is now implemented at SM 1‑loop for:
  - \(\alpha_2^{-1}(Q)\) (suite `v2`), and
  - \(\alpha_{1,\mathrm{GUT}}^{-1}(Q)\) (suite `v3`).
  A derived \(\sin^2\theta_W(Q)\) consistency check is available via `oos-ew-mix` (measured multi‑\(Q\) \(\sin^2\theta_W\) targets are still pending).

---

## 10) Reproduce (minimal)

```bash
# Baseline predictive OOS (step-only)
python -m physics_test.cli oos-predictive --suite v1 --max-rel-err 0.02

# RG-within-band predictive OOS (strong + EM + EW running)
python -m physics_test.cli oos-predictive-rg --suite v1 --max-rel-err 0.02
python -m physics_test.cli oos-predictive-rg --suite v2 --max-rel-err 0.02
python -m physics_test.cli oos-predictive-rg --suite v3 --max-rel-err 0.02

# EW mixing derived check (sin^2thetaW(Q) from alpha2 + alpha1_GUT running)
python -m physics_test.cli oos-ew-mix --max-rel-err 0.02

# GUT convergence diagnostic (baseline vs lattice-quantized initialization)
python -m physics_test.cli gut-run --model sm --n 400
python -m physics_test.cli gut-run-lattice --model sm --n 400
python -m physics_test.cli gut-run-lattice --model mssm --n 400

# Frozen OOS suites (C can vary per target, strict menu)
python -m physics_test.cli oos-report --suite v2 --max-rel-err 0.02
python -m physics_test.cli oos-report --suite v3 --max-rel-err 0.02
python -m physics_test.cli oos-report --suite v4 --max-rel-err 0.02
python -m physics_test.cli oos-report --suite v5 --max-rel-err 0.05
python -m physics_test.cli oos-report --suite v6 --max-rel-err 0.05
python -m physics_test.cli oos-report --suite v7 --max-rel-err 0.05
```
