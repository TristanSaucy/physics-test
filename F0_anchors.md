# Frequency anchors (F0) and temperature anchors (K)

This doc is about making Option‑2 (phenomenon‑first) and Option‑1 (energy‑first) **predictive** by freezing what counts as a legitimate $F_0$ and what $K$ means for each sector.

## 1) Why “$F_0 = E/h$” is called a proxy

When we write:

$$
F_0 \equiv \frac{E}{h},
$$

we are using Planck’s relation as a **unit conversion** from an energy scale to a frequency scale.

- For **photons**, this is literally measured: a monochromatic EM wave at frequency $f$ carries photon energy $E=hf$. So EM spectral lines are true $F_0$ measurements.
- For **particle masses / interaction scales**, $E/h$ is often a *timescale proxy* (or a quantum phase frequency) rather than a directly observed oscillation in time.
  - Example: $m_W c^2/h$ is a perfectly well-defined frequency scale, but there isn’t a lab instrument that directly sees a classical oscillation at that frequency the way it sees a radio tone.

So the “con” is not that $E/h$ is wrong—it’s that it introduces an extra modeling choice:

- which **energy scale $E$** represents “the phenomenon” for the force in question?

If we let $E$ float freely, it becomes an overfitting knob.

## 2) What a “principled choice” of energy scale can mean

A principled choice is one that is:

- **Externally fixed** (measured/defined elsewhere; not tuned to make a fit),
- **Standard in the literature**, and
- **Unique enough** that we can freeze it (one primary choice, at most a couple cross-checks).

Concrete principled rules you can choose from:

- **Rule P (pole/resonance rule)**: use a pole mass or a sharp resonance energy that defines a sector.
  - Weak: $m_W c^2$ or $m_Z c^2$.
  - Strong: a hadronic resonance scale (e.g., $m_p c^2$, $m_\pi c^2$) rather than an adjustable parameter.
  - EM: a specific atomic spectral line (Lyman‑α) or a standard atomic scale (Rydberg).

- **Rule L (Lagrangian-parameter rule)**: use a fundamental scale parameter in the theory.
  - Strong: $\Lambda_{QCD}$ (note: scheme dependence).
  - Weak: electroweak scale $v\approx 246$ GeV (from $G_F$).
  - Pro: conceptually fundamental; Con: scheme/convention baggage.

- **Rule M (measurement-context rule)**: use the energy scale at which the coupling target is defined.
  - Pro: reduces mismatch between “coupling at scale $Q$” and “phenomenon at scale $E$”.
  - Con: can erase sector differences if every coupling is quoted at the same $Q$ (e.g., $m_Z$).

Our recommendation for falsifiability: use **Rule P** for EM/weak, and a hybrid of **Rule P** and **Rule L** for strong (one primary + one cross-check).

## 3) Recommended “freeze menu” (primary + cross-checks)

The goal is to pick **one primary anchor per force** and keep 1–2 cross-check anchors that are not free alternatives but sanity checks.

These correspond to existing CLI presets in `physics_test/presets.py`.

### EM (F0 is directly measurable)

- **Primary**: `em-lyman-alpha` (Hydrogen Lyman‑α line; a clean spectral line)
- **Cross-checks**:
  - `em-visible-500THz` (broad human-visible band anchor)
  - `em-hydrogen-13.6eV` (atomic energy scale expressed as $E/h$; “series limit” style scale)

### Strong (no single “measured oscillation”; use energy/time proxies)

- **Primary**: `strong-QCD-200MeV` (QCD/confinement-scale proxy via $E/h$)
- **Cross-checks**:
  - `strong-proton-938MeV` (physical hadron mass scale via $E/h$)
  - `strong-timescale-1e-23s` (historical strong-interaction timescale proxy $F_0\sim 1/\tau$)

### Weak (mediator mass vs decay rate)

- **Primary**: `weak-W-80.379GeV` (pole mass energy proxy via $E/h$)
- **Cross-checks**:
  - `weak-Z-91.1876GeV` (same idea; nearby scale)
  - `weak-muon-decay` (true measured rate proxy $F_0\sim 1/\tau_\mu$; use as a *cross-check* because it’s a very different notion of “frequency”)

### Gravity (already frozen by type/band)

For gravity, we do not pick “one $F_0$”—we freeze **bands/types** and treat $F_0$ as constrained by the observational window (LIGO/LISA/PTA/CMB), while coupling targets use inverse gravity anchors:

- ordinary matter: `1/alpha_G(p)` (cross-check `1/alpha_G(e)`)
- GW types: `1/alpha_G(GW_CMB)`, `1/alpha_G(GW_PTA)`, `1/alpha_G(GW_LISA)`, `1/alpha_G(GW_LIGO)`
- Planck/quantum: `1/alpha_G(mP)`

## 4) How K fits in (what we freeze)

We currently treat:

- **Gravity**: $K$ literal (CMB) for GW‑band types.
- **Micro forces**: $K$ as an “energy-scale temperature” when using $E/h$ proxies:

$$
K \equiv \frac{E}{k_B}.
$$

When using truly observed frequencies (spectral lines or rates), Option‑2 implies:

$$
K = \frac{F_0 h}{k_B \phi^m},
$$

and then we interpret the implied $K$ as a scale parameter and check plausibility.

## 5) Reproduce with the CLI (Option‑2)

List presets:

```bash
python -m physics_test.cli list-frequency-presets
```

Use presets in the strict gauge‑C pairing workflow:

```bash
python -m physics_test.cli pair-forces-gaugeCs --em-preset em-lyman-alpha --strong-preset strong-QCD-200MeV --weak-preset weak-W-80.379GeV --gravity-band cmb --gravity-targets "1/alpha_G(p),1/alpha_G(GW_CMB)"
```

## 6) Citations and numeric source notes (to keep anchors defensible)

We separate:

- **the external quantity we cite** (wavelength, mass, lifetime, conventional scale), from
- **our derived frequency number** (usually via $f=c/\lambda$ or $f=E/h$).

The source of truth for the numeric presets used by the CLI is `physics_test/presets.py`.

### EM anchors

- **`em-lyman-alpha`**
  - **What it is**: Hydrogen Lyman‑α transition (vacuum wavelength near 121.6 nm).
  - **How we map to $F_0$**: $F_0=c/\lambda$ (we use an order-of-magnitude rounded wavelength).
  - **Citation**: NIST Atomic Spectra Database (lines): [`physics.nist.gov/PhysRefData/ASD/lines_form.html`](https://physics.nist.gov/PhysRefData/ASD/lines_form.html)
  - **Note**: our preset is intentionally rounded (it’s an anchor, not a high-precision spectroscopy database).

- **`em-hydrogen-13.6eV`**
  - **What it is**: hydrogenic atomic energy scale often quoted as 13.6 eV (Rydberg/ionization-scale proxy).
  - **How we map to $F_0$**: $F_0=E/h$.
  - **Citation**: NIST (fundamental constants + atomic physics references): [`physics.nist.gov/cuu/Constants/`](https://physics.nist.gov/cuu/Constants/)
  - **Note**: 13.6 eV is a conventional rounded anchor; if we ever need higher precision, we should pick a specific transition energy from NIST ASD instead.

### Weak anchors

- **`weak-W-80.379GeV`** and **`weak-Z-91.1876GeV`**
  - **What they are**: W and Z boson masses (pole mass benchmarks).
  - **How we map to $F_0$**: treat $E\approx mc^2$ and use $F_0=E/h$.
  - **Citation**: Particle Data Group (PDG) Review of Particle Physics: [`pdg.lbl.gov/`](https://pdg.lbl.gov/)
  - **Note**: PDG values can update over time; our presets are pinned to the numeric values in `physics_test/presets.py` for reproducibility.

- **`weak-muon-decay`**
  - **What it is**: muon lifetime $\tau_\mu$ (a directly measured decay constant).
  - **How we map to $F_0$**: $F_0\sim 1/\tau_\mu$.
  - **Citation**: PDG: [`pdg.lbl.gov/`](https://pdg.lbl.gov/)
  - **Note**: this is a *rate* anchor (not an $E/h$ proxy), so we keep it as a cross-check.

### Strong anchors

- **`strong-QCD-200MeV`**
  - **What it is**: a conventional confinement/hadronic scale anchor ($\sim$200 MeV).
  - **How we map to $F_0$**: $F_0=E/h$.
  - **Citations / justification**:
    - PDG QCD overview and scale dependence discussions: [`pdg.lbl.gov/`](https://pdg.lbl.gov/)
    - Rule-of-thumb physical justification: hadron size $\sim 1$ fm corresponds to momentum/energy scale $\hbar c/(1\text{ fm})\approx 200$ MeV.
  - **Note**: unlike W/Z masses, “the QCD scale” is not unique (scheme/threshold dependent). That’s why we keep at least one additional cross-check.

- **`strong-proton-938MeV`**
  - **What it is**: proton mass-energy scale (a concrete hadronic mass benchmark).
  - **How we map to $F_0$**: $F_0=E/h$ with $E\approx m_pc^2$.
  - **Citation**: PDG: [`pdg.lbl.gov/`](https://pdg.lbl.gov/)

- **`strong-timescale-1e-23s`**
  - **What it is**: historical order-of-magnitude “strong interaction timescale” proxy.
  - **How we map to $F_0$**: $F_0\sim 1/\tau$.
  - **Justification**: a characteristic length $\sim$1 fm implies a light-crossing time $\sim 1\text{ fm}/c \approx 3\times 10^{-24}$ s; $10^{-23}$ s is the same order of magnitude.
  - **Note**: we keep this as a cross-check only, because it is not a uniquely defined observable.
