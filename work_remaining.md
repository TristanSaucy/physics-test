## Work remaining / next steps

This is a living checklist for continuing the framework work in future sessions.

### 1) Freeze the “rules of the game” (reduce degrees of freedom)
- **Frozen (done): coupling orientation**
  - Canonical “gauge strength coordinate” is **inverse couplings** (`1/alpha`-style) for strict claims.
- **Frozen (done): EM target**
  - EM strict target is **`1/alpha` (low-energy)**. (`1/alpha(mZ)` remains exploratory only.)
- **Frozen (done): strong/weak default targets (strict)**
  - Strong strict target: **`1/alpha_s(mZ)`** (scale is explicit and frozen at \(m_Z\) for the benchmark).
  - Weak strict target: **`1/alpha_w(mZ)`** (same \(m_Z\) reference convention).
- **Frozen (done): gravity orientation**
  - Gravity strict targets use **inverse coupling**: `1/alpha_G(mass)`.
- **Frozen (done): ordinary-matter gravity anchor + cross-check**
  - Canonical ordinary-matter gravity target: **`1/alpha_G(p)`** (proton mass scale).
  - Mandatory cross-check (not a free alternative): **`1/alpha_G(e)`** (electron mass scale), which should appear with an \(m\)-shift consistent with \((m_p/m_e)^2\) (see `paper.md` §5.2).
- **Frozen (done): other gravity “types” (GW-band / primordial / quantum)**
  - Frozen GW-band type anchors (inverse gravity, derived under CMB K; see `paper.md` §5.4):
    - **GW_CMB**: \(M = 2.93012\times 10^4\) GeV (~29.3 TeV)
    - **GW_PTA**: \(M = 1.58009\times 10^9\) GeV
    - **GW_LISA**: \(M = 2.15524\times 10^{12}\) GeV
    - **GW_LIGO**: \(M = 5.41086\times 10^{13}\) GeV
  - Frozen “Planck/quantum gravity” type (separate from GW bands):
    - **mP**: \(M \sim 10^{19}\) GeV (Planck scale; \(1/\\alpha_G(m_P)\\approx 1\\), can land at positive \(m\) under strict gauge-\(C\))
- **Frozen (done): strict C candidates**
  - Strict gauge-derived set: **\(C\\in\\{360,180,120,60,45,15\\}\\)** (derived from SM gauge-group invariants from base=360).
  - Avoid large “octave-union” families unless explicitly justified (keep them exploratory).

### 2) Make K and F0 physically anchored (especially for micro forces)
- Decide whether **K is literal thermodynamic temperature** or **energy-scale temperature** (E = kB*K).
- For Option 2 (“phenomenon-first”), decide what counts as a legitimate phenomenon frequency:
  - EM: spectral line or band (e.g. Lyman-alpha, visible band).
  - Strong: characteristic hadronic/nuclear resonance frequency vs E/h proxy.
  - Weak: particle scale (W/Z) vs decay rates (muon lifetime) — these imply wildly different K.
- Add a **small curated table** of chosen F0 anchors with citations/notes so choices aren’t arbitrary.

### 3) Gravity “types” / bands
- Decide what “gravity type” corresponds to:
  - **CMB/primordial GW band** (m implied by K=CMB).
  - **PTA/LISA/LIGO bands** (likely require a different K and/or different coupling definition).
- Expand the gravity coupling menu systematically:
  - `alpha_G(mass)` and `1/alpha_G(mass)` for physically motivated mass scales.
  - Map “quantum gravity” mass scales to which GW band they align with under K=CMB.
- Add a report that summarizes **which gravity targets land in which band** under strict gauge-derived C.

### 4) Elevate “gauge-derived C” from toy to principled
- Right now we derive C from base=360 using simple invariants (dim, Coxeter).
- Next: test other invariant constructions (e.g., Casimir factors, representation choices) and document them.
- Goal: a *short* list of allowed C values with a clear derivation and minimal tuning.

### 5) Add falsifiable predictions / scoring
- Introduce a single “score” per configuration combining:
  - coupling fit error(s),
  - band constraint satisfaction,
  - any K plausibility windows,
  - penalties for free choices.
- Document what would count as a **failure mode** (no configurations survive under strict constraints).

### 6) Make analysis easier to inspect
- Add a notebook section that:
  - plots m vs target coupling for fixed C sets,
  - plots gravity band-implied m windows,
  - shows sensitivity to target choice (alpha vs 1/alpha etc.).
- Add CSV export for sweep results.

### 7) Documentation hygiene
- Keep these “living docs” in sync whenever rules/results change:
  - `freezing_the_rules.md` (the contract: what is frozen vs exploratory)
  - `active_context.md` (resume context + current defaults)
  - `repo_brain.md` (narrative explanation + reproduction commands)
  - `work_remaining.md` (next steps checklist)
  - `bio_sanity_checks.md` (exploratory bio-related probes)
- Keep a short, copy/paste-able “reproduce strict core results” command list in **both**:
  - `README.md` (user-facing quickstart)
  - `repo_brain.md` (durable context)


