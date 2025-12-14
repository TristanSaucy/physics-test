# Repo “brain” (durable context + how to think about the project)

If you want to keep iterating without re-deriving everything, read this first.

## 1) What we’re trying to test

We’re testing whether a minimal harmonic/topology model can connect:

- **dimensionless couplings** (pure numbers) via

$$
G = \frac{C}{\phi^m}
$$

- to **frequency/temperature scaling** via

$$
F_0 = \phi^m\,\frac{k_B K}{h}
$$

with **integer m** as a discrete “trophic/harmonic step.”

The model’s built-in invariant:

$$
G\,F_0 = C\,\frac{k_B K}{h}
$$

## 2) What has been most promising (signals)

### EM + 360

- $C=360, m=2$ gives $G\approx 137.5078$, close to **$1/\alpha$** (0.34%).
- This is the cleanest “360 is special” result so far.

### Gauge-derived non-arbitrary C set

We derived candidate $C$ values from simple gauge invariants (U(1), SU(2), SU(3)) using base=360:

- 360, 180, 120, 60, 45, 15 (and a couple duplicates)

Scanning only those C values, integer m, and 5% tolerance gives plausible fits for:

- EM (`1/alpha`) at C=360, m=2
- Strong (strict) `1/alpha_s_1loop_from_mZ(mH)` at C=60, m=4
- Weak (strict) `1/alpha2(alpha(mZ),sin2_on_shell)` at C=120, m=3

This is valuable because it reduces “arbitrary C.”

### Strict all-forces per GW band (new)

With frozen Option‑2 phenomenon anchors (see `F0_anchors.md`), strict gauge-$C$, and CMB $K$ for gravity:

- `pair-forces-gaugeCs` yields **one clean configuration per GW band** (CMB/PTA/LISA/LIGO).
- **EM/strong/weak are identical across bands**; only gravity’s $(C,m)$ shifts to land in the chosen band window.
- The compact tables + reproduce commands are in `gravity_types_report.md` and `paper.md` §5.4.2.

### Gravity: why “type” matters

With gravity K fixed to the CMB, the frequency equation forces gravity’s m into narrow windows for each GW band.
Reconciling gravity coupling fits with those m windows typically requires changing the **gravity coupling definition**, e.g. using `1/alpha_G(scale)` at TeV-scale “quantum gravity” mass scales.

## 3) Main risk / how to stay falsifiable

The risk is overfitting: if C is allowed to roam freely or if we keep changing target definitions ad hoc, we can always find matches.

To stay falsifiable:

- Keep C restricted to a **small, justified set** (e.g., gauge-derived).
- Fix target definitions and reference scales (strong/weak run with scale).
- Fix GW band constraints if you claim “gravity type X lives here.”
- Prefer experiments where a constraint *eliminates* most candidates.

## 4) The meaning of K

K can mean:

- literal thermodynamic temperature (works naturally for gravity+CMB), or
- an energy-scale temperature (E = kB*K) for micro physics.

Option 2 (“phenomenon-first”) is useful because it forces K to be implied by m and the chosen F0, rather than being tuned.

## 5) Where to work next (high-value)

Now that the strict contract, gravity types, and Option‑2 $F_0$ anchors are frozen, the highest value work is:

1) Tighten the tolerance and re-test survival rates (5% → 2% → 1%).
2) Add a simple **score** / ranking for configurations (penalize free choices; reward cross-check stability).
3) Expand or justify the gauge-derived $C$ construction (Casimirs / reps) only if the allowed set stays short.
4) Stress-test robustness under cross-check $F_0$ presets (do the same $(C,m)$ pattern persist?).

## 6) Reproduce strict core results (copy/paste)

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
```

Optional gravity mass sweep (inverse gravity; CMB $K$):

```bash
python -m physics_test.cli sweep-quantum-gravity --gravity-band cmb --scale-min-GeV 1e3 --scale-max-GeV 1e19 --n-scales 241 --top 12
```

## 7) Visualization

Open `notebooks/topological_gauge_explorer.ipynb` and run the new plotting cells to see log10(G) vs m curves for gauge-derived C candidates and highlight “hit” points within tolerance for targets like:

- `sin2thetaW(mZ)`
- `alpha_w(mZ)`
- `alpha_s(mZ)`
- `1/alpha`
