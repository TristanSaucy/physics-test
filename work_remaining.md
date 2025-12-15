# Work remaining / next steps

This is a living checklist for continuing the framework work in future sessions.

## 1) Freeze the “rules of the game” (reduce degrees of freedom)

- **Frozen (done): coupling orientation**
  - Canonical “gauge strength coordinate” is **inverse couplings** (`1/alpha`-style) for strict claims.
- **Frozen (done): EM target**
  - EM strict target is **`1/alpha` (low-energy)**. (`1/alpha(mZ)` remains exploratory only.)
- **Frozen (done): strong/weak default targets (strict)**
  - Strong strict target: **`1/alpha_s_1loop_from_mZ(mH)`** (1-loop run from `alpha_s(mZ)` to $m_H=125$ GeV; no free $\Lambda_{\mathrm{QCD}}$ knob).
  - Weak strict target: **`1/alpha2(alpha(mZ),sin2_on_shell)`** (derived using on-shell $\sin^2\theta_W=1-m_W^2/m_Z^2$).
- **Frozen (done): gravity orientation**
  - Gravity strict targets use **inverse coupling**: `1/alpha_G(mass)`.
- **Frozen (done): ordinary-matter gravity anchor + cross-check**
  - Canonical ordinary-matter gravity target: **`1/alpha_G(p)`** (proton mass scale).
  - Mandatory cross-check (not a free alternative): **`1/alpha_G(e)`** (electron mass scale), which should appear with an $m$-shift consistent with $(m_p/m_e)^2$ (see `paper.md` §5.2).
- **Frozen (done): other gravity “types” (GW-band / primordial / quantum)**
  - Frozen GW-band type anchors (inverse gravity, derived under CMB $K$; see `paper.md` §5.4):
    - **GW_CMB**: $M = 2.93012\times 10^4$ GeV (~29.3 TeV)
    - **GW_PTA**: $M = 1.58009\times 10^9$ GeV
    - **GW_LISA**: $M = 2.15524\times 10^{12}$ GeV
    - **GW_LIGO**: $M = 5.41086\times 10^{13}$ GeV
  - Frozen “Planck/quantum gravity” type (separate from GW bands):
    - **mP**: $M \sim 10^{19}$ GeV (Planck scale; $1/\alpha_G(m_P)\approx 1$, can land at positive $m$ under strict gauge-$C$)
- **Frozen (done): strict C candidates**
  - Strict gauge-derived set: **$C\in\{360,180,120,60,45,15\}$** (derived from SM gauge-group invariants from base=360).
  - Avoid large “octave-union” families unless explicitly justified (keep them exploratory).

## 2) Make K and F0 physically anchored (especially for micro forces)

- **Frozen (done): Option‑2 anchor menu + K interpretation**
  - See `F0_anchors.md` (rationale + frozen primary/cross-check presets).
  - Frozen interpretation: gravity GW‑types use literal CMB $K$; micro forces treat $K$ as an effective scale parameter.
- **Frozen (done): citations/notes for anchors**
  - Added a citations + numeric-notes section in `F0_anchors.md` tying each frozen preset back to PDG/NIST (and noting which are intentionally approximate anchors).

## 3) Gravity “types” / bands

- **Frozen (done): GW-band “gravity types”**
  - Frozen type anchors and targets are listed in section #1 and in `freezing_the_rules.md`.
- **Frozen (done): report**
  - Added `gravity_types_report.md` summarizing which frozen gravity targets land in which GW band under strict gauge-derived $C$ (including strict all-forces per-band configurations).

## 4) Elevate “gauge-derived C” from toy to principled

- Right now we derive C from base=360 using simple invariants (dim, Coxeter).
- Added background doc: `gauge_invariants.md` (dim/rank/Coxeter/dual Coxeter + “next tier” invariants like roots/Weyl group + Casimir/Dynkin index and the representation/normalization pitfalls).
- Next: test other invariant constructions (e.g., Casimir factors, representation choices) and document them.
- Goal: a *short* list of allowed C values with a clear derivation and minimal tuning.

## 5) Add falsifiable predictions / scoring

- Introduce a single “score” per configuration combining:
  - coupling fit error(s),
  - band constraint satisfaction,
  - any K plausibility windows,
  - penalties for free choices.
- Document what would count as a **failure mode** (no configurations survive under strict constraints).
- Added frozen out-of-sample test suites (v1/v2/v3/v4): `python -m physics_test.cli oos-report --suite v1|v2|v3|v4`
- Added predictive OOS (fit one C per force from strict anchors, then hold C fixed): `python -m physics_test.cli oos-predictive --suite v1`
- Added principled normalization-factor families for predictive OOS: `python -m physics_test.cli list-norm-families` and `python -m physics_test.cli oos-predictive --norm-family ...`
- Added step-signal OOS (C-independent): `python -m physics_test.cli oos-steps --suite v1 --max-ratio-err 0.02|0.05`
- Added optional discrete C-ratio step-signal diagnostic (exploratory): `python -m physics_test.cli oos-steps --allow-c-ratio`
- Added RG+phi OOS (Lambda_QCD consistency from lattice-fit inverse couplings): `python -m physics_test.cli oos-rg --suite qcd-lambda-v1|v2`

## 6) Make analysis easier to inspect

- Add a notebook section that:
  - plots m vs target coupling for fixed C sets,
  - plots gravity band-implied m windows,
  - shows sensitivity to target choice (alpha vs 1/alpha etc.).
- Add CSV export for sweep results.

## 7) Documentation hygiene

- **Done**: kept “living docs” in sync (`README.md`, `paper.md`, `freezing_the_rules.md`, `active_context.md`, `repo_brain.md`, `gravity_types_report.md`).
- **Done**: added a short, copy/paste “reproduce strict core results” command list in both `README.md` and `repo_brain.md`.
