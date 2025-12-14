# Bio sanity checks (exploratory)

This file is intentionally **not** part of strict/frozen claims. It’s a place to log “does this collide with biology?” probes.

## 1) Universal thermal frequency scale (anchor)

From the model:

- $F_0 = \phi^m\,\frac{k_B K}{h}$

So at a given temperature $K=T$, the **baseline** ($m=0$) is:

- $f_T = k_B T/h$

At human body temperature $T\approx 310\,\mathrm{K}$:

- $f_T \approx 6.46\,\mathrm{THz}$
- $m=1\Rightarrow F_0 = \phi f_T \approx 10.45\,\mathrm{THz}$

Reproduce:

- `python -m physics_test.cli calc --m 0 --K 310`
- `python -m physics_test.cli calc --m 1 --K 310`
- `python -m physics_test.cli list-frequency-presets` (thermal presets are listed)

## 2) “Brownian motion frequency” is not universal

“Brownian motion” does not define a single characteristic frequency without specifying a system (particle size, viscosity, trap stiffness, etc.).

If we want a *principled* biology-adjacent frequency anchor, $k_B T/h$ is universal; Brownian roll-off frequencies are system-dependent.

## 3) Microtubules in the THz band (status + how to test in this framework)

What we can say safely:

- **Many biomolecules (proteins, hydration shells)** have broadband dielectric/absorptive features in the sub-THz..THz range because THz probes collective vibrational/rotational modes.
- Showing **“there exists THz absorption”** is not the same as showing a **sharp resonance** or a **coherent, long-lived mode** that could plausibly play an organizing role at room/body temperature.
- In aqueous environments, **damping is typically strong**, so claims of long-lived coherence in the THz band are usually the hard part.

What would count as a meaningful “hit” (evidence tiers):

- **Tier A (weak)**: published THz time-domain spectroscopy / dielectric spectroscopy showing microtubule/tubulin-dependent features in the THz range (broad is fine).
- **Tier B (medium)**: a reproducible, relatively **narrow** resonance feature (or set of modes) attributable to microtubules/tubulin, not just bulk water/protein background.
- **Tier C (strong)**: evidence for **non-thermal coherence** (coherence times or phase correlations inconsistent with purely thermal noise) *in situ* or in a biologically realistic environment.

How this plugs into the model (if we pursue it):

- If we pick $F_0$ in the **~6–12 THz** neighborhood (thermal baseline at 300–310 K, and the $m=1$ step),
  then Option 2 converts that into an implied $m$ for any target coupling $G$ via the fitted $m$ from $G=C/\phi^m$.
- The falsifiable question becomes: **do strict C candidates + frozen target choices yield an $m$ that predicts a THz-band $F_0$ at realistic $K$** without tuning?

Open work (needs real citations):

- Curate a short list of experimental papers (THz-TDS/dielectric) on tubulin/microtubules.
- Add a small “observed THz features table” (frequency, linewidth, sample conditions).
