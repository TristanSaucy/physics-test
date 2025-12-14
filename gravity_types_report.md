# Gravity “types” and GW bands (strict report)

This report summarizes how the repo’s **frozen gravity coupling targets** relate to **gravitational-wave (GW) frequency bands** under the strict contract:

- Gauge-derived $C\in\{360,180,120,60,45,15\}$
- Integer $m$
- Inverse-coupling coordinate for gravity: `1/alpha_G(mass)`

It also clarifies what is “measured” in each band and what we mean by a “gravity type”.

## 1) Two different things: coupling target vs GW frequency band

- The **coupling target** is a *dimensionless number*:

$$
\alpha_G(M)=\frac{G_N M^2}{\hbar c},
\qquad
\frac{1}{\alpha_G(M)}=\frac{\hbar c}{G_N M^2}.
$$

This depends on a chosen mass scale $M$.

- The **GW band** is an *observational frequency window* for $F_0$ (Hz).

In this framework, the bridge is the frequency law (with $K$ interpreted per type):

$$
F_0=\phi^m\frac{k_BK}{h},
\qquad
G=\frac{C}{\phi^m}.
$$

So once you pick a target $G$ (coupling) and a $K$, you imply an $m$ and thus an $F_0$.

## 2) GW bands and where measurements are made

These are order-of-magnitude sensitivity windows (also summarized in `paper.md` §5.4.1):

- **LIGO band** (10–1000 Hz): ground-based laser interferometers (LIGO Hanford/Livingston; Virgo; KAGRA).
- **LISA band** (10⁻⁴–10⁻¹ Hz): planned space-based interferometer (three spacecraft; heliocentric orbit trailing Earth).
- **PTA band** (10⁻⁹–10⁻⁷ Hz): pulsar timing arrays (radio telescopes timing millisecond pulsars).
- **CMB/primordial** (~10⁻¹⁸–10⁻¹⁶ Hz): indirect constraints from CMB anisotropy/polarization (not a direct interferometric detection).

## 3) Band-implied $m$ windows under $K = K_{CMB}$

If we treat $K$ as literal CMB temperature ($K=2.725$ K), the frequency equation alone implies a narrow integer-$m$ window for each band (independent of $C$):

- **CMB band**: $m\approx -137..-129$
- **PTA band**: $m\approx -94..-85$
- **LISA band**: $m\approx -70..-57$
- **LIGO band**: $m\approx -46..-38$

## 4) Frozen gravity targets and what bands they land in

### 4.1 Ordinary matter gravity (not a GW-band type)

These are frozen as “ordinary matter” coupling targets:

- Canonical: `1/alpha_G(p)` (proton)
- Cross-check: `1/alpha_G(e)` (electron)

Under $K=2.725$ K, the strict best fits imply extremely small $F_0$ (far below GW detector bands), which is consistent with treating these as *not* “GW-band gravity types”.

### 4.2 Frozen GW-band gravity types (4 anchors)

These are frozen as named targets (mass anchors derived from strict band sweeps under CMB $K$):

- `1/alpha_G(GW_CMB)`
- `1/alpha_G(GW_PTA)`
- `1/alpha_G(GW_LISA)`
- `1/alpha_G(GW_LIGO)`

### 4.3 Frozen Planck/quantum type (+1)

- `1/alpha_G(mP)` (~1) (Planck mass scale)

This is frozen as a separate “quantum gravity” type and is not expected to land in the GW bands under CMB $K$.

## 5) Strict fit table (gauge-derived C only)

All rows below are strict gauge-derived $C$ hits within 5% using `scan-gauge-Cs` (best hit shown), with $F_0$ computed at $K=2.725$ K using `calc`.

| Gravity target | Frozen meaning | Best (C,m) | rel_err | $F_0(K=2.725K)$ | Band match |
|---|---|---:|---:|---:|---|
| `1/alpha_G(p)` | ordinary matter (proton) | (45, −175) | −6.07e−03 | 1.52e−26 Hz | none |
| `1/alpha_G(e)` | cross-check (electron) | (360, −202) | +3.58e−02 | 3.46e−32 Hz | none |
| `1/alpha_G(GW_CMB)` | GW_CMB type | (45, −132) | −3.05e−06 | 1.47e−17 Hz | CMB |
| `1/alpha_G(GW_PTA)` | GW_PTA type | (15, −89) | +3.97e−06 | 1.43e−08 Hz | PTA |
| `1/alpha_G(GW_LISA)` | GW_LISA type | (15, −59) | +3.38e−07 | 2.65e−02 Hz | LISA |
| `1/alpha_G(GW_LIGO)` | GW_LIGO type | (360, −39) | −2.30e−07 | 4.01e+02 Hz | LIGO |
| `1/alpha_G(mP)` | Planck/quantum type | (120, 10) | −2.43e−02 | (CMB $K$) very high | none |

Notes:

- The GW_* targets land cleanly in their intended GW bands by construction (and we now treat them as frozen “types” rather than continuing to sweep).
- The proton/electron ordinary-matter targets do not land in detector bands under CMB $K$; that’s a feature, not a bug, given they represent a different “gravity type.”

## 6) Strict “all-forces” configurations per GW band (Option-2 for EM/strong/weak)

Here we run the strict “pair all forces” search with:

- Gauge-derived $C\in\{360,180,120,60,45,15\}$
- Fixed phenomenon frequencies (frozen in `F0_anchors.md`)
- Targets: EM=`1/alpha`, strong=`1/alpha_s(mZ)`, weak=`1/alpha_w(mZ)`
- Gravity target chosen per band (one of the frozen GW_* targets)
- Gravity uses $K=2.725$ K and must land inside the named GW band window

**Key observation:** under the frozen EM/strong/weak anchors, the strict solution for EM/strong/weak is the same across bands; only the gravity $(C,m)$ changes to land in the chosen GW band.

Best configurations (first hit shown per band):

| GW band | Gravity target | EM (C,m,K) | Strong (C,m,K) | Weak (C,m,K) | Gravity (C,m) -> $F_0$ (Hz) |
|---|---|---|---|---|---|
| CMB | `1/alpha_G(GW_CMB)` | (360, +2, 4.5205e4 K) | (60, +4, 3.3862e11 K) | (120, +3, 2.2020e14 K) | (45, −132) -> 1.472e−17 |
| PTA | `1/alpha_G(GW_PTA)` | (360, +2, 4.5205e4 K) | (60, +4, 3.3862e11 K) | (120, +3, 2.2020e14 K) | (15, −89) -> 1.427e−08 |
| LISA | `1/alpha_G(GW_LISA)` | (360, +2, 4.5205e4 K) | (60, +4, 3.3862e11 K) | (120, +3, 2.2020e14 K) | (15, −59) -> 2.654e−02 |
| LIGO | `1/alpha_G(GW_LIGO)` | (360, +2, 4.5205e4 K) | (60, +4, 3.3862e11 K) | (120, +3, 2.2020e14 K) | (360, −39) -> 4.015e+02 |

Rel-error context (from the CLI run):

- EM rel_err ≈ +3.44e−03
- strong rel_err ≈ +3.21e−02
- weak rel_err ≈ −4.17e−02
- gravity rel_err is band-dependent but is $\ll 1\%$ for the PTA/LISA/LIGO rows above.

## 7) Reproduce

List gravity bands and targets:

```bash
python -m physics_test.cli list-gravity-bands
python -m physics_test.cli list-targets
```

Strict gauge-$C$ scans:

```bash
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(p)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(e)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(GW_CMB)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(GW_PTA)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(GW_LISA)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(GW_LIGO)" --max-rel-err 0.05
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_G(mP)" --max-rel-err 0.05
```

Compute the implied $F_0$ (CMB $K$) for a given strict hit:

```bash
python -m physics_test.cli calc --C 45 --m -132 --K 2.725
```

Strict all-forces configurations per band:

```bash
python -m physics_test.cli pair-forces-gaugeCs --gravity-band cmb  --gravity-targets "1/alpha_G(GW_CMB)"  --max-hits 10 --max-results 5
python -m physics_test.cli pair-forces-gaugeCs --gravity-band pta  --gravity-targets "1/alpha_G(GW_PTA)"  --max-hits 10 --max-results 5
python -m physics_test.cli pair-forces-gaugeCs --gravity-band lisa --gravity-targets "1/alpha_G(GW_LISA)" --max-hits 10 --max-results 5
python -m physics_test.cli pair-forces-gaugeCs --gravity-band ligo --gravity-targets "1/alpha_G(GW_LIGO)" --max-hits 10 --max-results 5
```
