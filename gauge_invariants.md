# Gauge-group invariants and “gauge-derived C” constructions

This note expands “#4” in `work_remaining.md`: how we derive non-arbitrary candidate values of the topological constant $C$ from **group invariants**, and what “more complex invariants” would mean without turning the framework into an overfitting machine.

## 1) Why use invariants at all?

The purpose of “gauge-derived $C$” is **not** to claim that a particular invariant *must* produce the right number. The purpose is to restrict $C$ to a **small, externally-motivated set** that:

- comes from unambiguous properties of the gauge group (not tuned to match data),
- is reproducible from a few lines of math, and
- stays short after de-duplication (anti-numerology constraint).

In this repo, strict mode currently freezes:

- **base** = 360
- groups = U(1), SU(2), SU(3)
- constructions = a small menu of “divide base by an invariant (or a simple product of invariants)”

The implementation lives in `physics_test/gauge_groups.py`.

## 2) “Simple” Lie-group invariants used today

These are invariants of the **Lie algebra / root system** (i.e., they don’t depend on a choice of representation).

### 2.1 Rank $r$

Informally: the number of independent commuting generators (dimension of a maximal torus). For classical groups:

- $SU(N)$: $r=N-1$

Rank is not yet used directly in the current strict $C$ generator, but it is a standard “first invariant” and becomes useful when counting roots (below).

### 2.2 Lie algebra dimension $\dim(\mathfrak{g})$

This is the number of generators; equivalently the dimension of the **adjoint representation**.

- $SU(N)$: $\dim = N^2-1$

Examples:

- $SU(2)$: $\dim=3$
- $SU(3)$: $\dim=8$

### 2.3 Coxeter number $h$

$h$ is an integer invariant of the root system. It can be defined multiple equivalent ways; one useful definition is:

- $h = 1 +$ (height of the highest root)

For $SU(N)$ (type $A_{N-1}$), the exponents are $1,2,\dots,N-1$ and:

- $h = N$

### 2.4 Dual Coxeter number $h^\vee$

$h^\vee$ is another root-system invariant. Under the standard normalization (long roots have length-squared 2), it equals the **quadratic Casimir of the adjoint**:

$$
C_2(\mathrm{adj}) = h^\vee.
$$

For simply-laced groups (including $SU(N)$), $h=h^\vee$.

For $SU(N)$:

- $h^\vee = N$

## 3) The current strict “gauge-derived $C$” menu (what it is and why it’s small)

The current strict generator takes a base value (default 360) and produces:

- $C=\text{base}$
- $C=\text{base}/\dim$
- $C=\text{base}/h$
- $C=\text{base}/h^\vee$
- $C=\text{base}/(\dim\cdot h)$

For U(1), only $\dim=1$ is used (no root system, so $h,h^\vee$ are undefined).

For SU(N), since $h=h^\vee=N$, the $h$ and $h^\vee$ constructions coincide.

### 3.1 Concrete values for the Standard Model factors

| Group | rank $r$ | $\dim(\mathfrak{g})$ | $h$ | $h^\vee$ |
|---|---:|---:|---:|---:|
| U(1) | 1 | 1 | — | — |
| SU(2) | 1 | 3 | 2 | 2 |
| SU(3) | 2 | 8 | 3 | 3 |

Using base = 360, the constructions yield:

- U(1): $360$
- SU(2): $360/3=120$, $360/2=180$, $360/(3\cdot 2)=60$
- SU(3): $360/8=45$, $360/3=120$, $360/(8\cdot 3)=15$

After de-duplication across U(1), SU(2), SU(3), strict mode freezes:

- $C\in\{360,180,120,60,45,15\}$

### 3.2 Construction keys implemented in the CLI

The CLI exposes “construction keys” through `--include` (comma-separated). Strict mode defaults to:

- `base`
- `base/dim`
- `base/coxeter`
- `base/dual_coxeter`
- `base/(dim*coxeter)`

The next-tier keys implemented (exploratory unless you freeze them) are:

- `base/rank`
- `base/roots`
- `base/positive_roots`
- `base/weyl_order`
- `base/center_order`
- `base/C2_adj` (alias for `base/dual_coxeter`)
- `base/C2_fund` (SU(N) fundamental Casimir, standard normalization)
- `base/T_fund` (SU(N) fundamental Dynkin index, standard normalization)

Example:

```bash
python -m physics_test.cli list-gauge-Cs --include "base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter),base/roots,base/positive_roots,base/weyl_order,base/center_order,base/C2_fund,base/T_fund"
```

### 3.3 Exploratory scan snapshot (3% tolerance)

We ran an exploratory check at **3%** tolerance to see whether the next-tier keys materially improve the strict inverse-coupling fits.

With the extended `--include` menu above, the unique candidate set becomes:

- $C\in\{15,45,60,120,180,270,360,480,720\}$ (the new values come from `base/C2_fund` and `base/T_fund`)

Results at **3%**:

- **EM inverse** still fits:
  - `1/alpha`: $(C,m)=(360,2)$ with rel_err ≈ +0.344%
- **Hypercharge (GUT norm.) inverse** still fits:
  - `1/alpha1_GUT(alpha(mZ),sin2)`: $(C,m)=(60,0)$ with rel_err ≈ +1.66%
- **Strong/weak inverse** do **not** fit at 3% under this candidate set:
  - `1/alpha_s(mZ)`: **no hits**; best remains $(60,4)$ with rel_err ≈ +3.21%
  - `1/alpha_w(mZ)`: **no hits**; best remains $(120,3)$ with rel_err ≈ −4.17%

This implies that tightening the strict inverse-coupling contract from 5% → 3% would currently require either:

- a broader (but still principled) $C$ menu, or
- different strong/weak target choices (scale/scheme), or
- revisiting the frozen orientation (exploratory only).

Reproduce:

```bash
python -m physics_test.cli scan-gauge-Cs --target "1/alpha" --include "base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter),base/roots,base/positive_roots,base/weyl_order,base/center_order,base/C2_fund,base/T_fund,base/rank" --max-rel-err 0.03
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_s(mZ)" --include "base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter),base/roots,base/positive_roots,base/weyl_order,base/center_order,base/C2_fund,base/T_fund,base/rank" --max-rel-err 0.03
python -m physics_test.cli scan-gauge-Cs --target "1/alpha_w(mZ)" --include "base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter),base/roots,base/positive_roots,base/weyl_order,base/center_order,base/C2_fund,base/T_fund,base/rank" --max-rel-err 0.03
python -m physics_test.cli scan-gauge-Cs --target "1/alpha1_GUT(alpha(mZ),sin2)" --include "base,base/dim,base/coxeter,base/dual_coxeter,base/(dim*coxeter),base/roots,base/positive_roots,base/weyl_order,base/center_order,base/C2_fund,base/T_fund,base/rank" --max-rel-err 0.03
```

## 4) “More complex” invariants we can consider next (and what they buy us)

There are two broad classes:

- **(A) Representation-independent invariants** (safer; fewer degrees of freedom)
- **(B) Representation-dependent invariants** (dangerous unless representation + normalization are frozen)

### 4.1 (A) More root-system invariants (representation-independent)

These can provide **alternative derivations** of the same integers we already see (good) or generate small new candidates (potentially useful).

- **Number of roots**:

$$
|\Phi| = \dim(\mathfrak{g}) - r.
$$

For SU(N): $|\Phi| = N(N-1)$.

- **Number of positive roots**:

$$
|\Phi^+| = \frac{\dim(\mathfrak{g}) - r}{2}.
$$

- **Weyl group order** $|W|$:
  - For SU(N): $|W| = N!$ (the symmetric group $S_N$).

- **Center order** $|Z(G)|$:
  - For SU(N): $|Z|=N$.

For $SU(N)$, these invariants are related: the product of $(e_i+1)$ over exponents equals $|W|$; and for $A_{N-1}$ the exponents are $1,2,\dots,N-1$, so $|W|=N!$.

**Why these are attractive:** they are canonical integers of the group/root system, no representation choice needed.

**Why they can be risky:** $|W|$ grows very fast, and naive constructions like base/$|W|$ can become tiny or non-integer, so we need strict rules about what constructions are allowed.

### 4.2 (B) Casimirs and Dynkin indices (representation-dependent)

Many physically meaningful quantities (beta functions, anomaly coefficients, interaction strengths in specific channels) involve **representation-dependent** invariants.

- **Quadratic Casimir** $C_2(R)$ (depends on representation $R$):

$$
\sum_a T_R^a T_R^a = C_2(R)\,\mathbf{1}.
$$

- **Dynkin index** $T(R)$:

$$
\mathrm{Tr}_R(T^a T^b) = T(R)\,\delta^{ab}.
$$

These are related by:

$$
\dim(\mathfrak{g})\,T(R)=\dim(R)\,C_2(R).
$$

For $SU(N)$ under standard normalization:

- adjoint: $C_A = C_2(\mathrm{adj}) = N$
- fundamental: $C_F = C_2(F) = \frac{N^2-1}{2N}$
- $T(F)=1/2$

**Why these are attractive:** they are central in real gauge theory calculations.

**Why they are dangerous for “$C$ mining”:**

- you must pick a representation ($F$ vs adjoint vs higher reps),
- you must pick a normalization convention,
- many outcomes are rational numbers, which can explode the candidate set unless strongly constrained.

## 5) A “principled ruleset” for expanding $C$ without exploding the hypothesis space

If we extend beyond the current toy menu, the main goal is to keep the allowed list short. A workable set of constraints is:

- **Rule 1 (no free reps)**: if a construction uses representation-dependent invariants ($C_2(R)$, $T(R)$, $\dim(R)$), then **freeze one representation choice** (e.g., adjoint only; or fundamental only) *per gauge factor*, and document it.
- **Rule 2 (no free normalizations)**: freeze a convention (e.g., long roots have length-squared 2; $T(F)=1/2$ for SU(N)).
- **Rule 3 (simple arithmetic only)**: allow only base divided by a single invariant, or by a product of at most two invariants. Avoid “bespoke” linear combinations.
- **Rule 4 (integer/small-denominator filter)**: require candidates to be integers (or rationals with a very small denominator) so the set does not explode.
- **Rule 5 (deduplicate aggressively)**: if multiple derivations yield the same $C$, that is a *feature*; keep the value once.
- **Rule 6 (freeze the menu)**: decide up front which invariants are allowed in strict mode, and treat all others as exploratory.

## 6) Reproduce (current strict behavior)

List the currently frozen strict candidates:

```bash
python -m physics_test.cli list-gauge-Cs
```

See:

- `physics_test/gauge_groups.py` for the current invariant fields and constructions
- `paper.md` §3.2 for the strict contract summary
- `freezing_the_rules.md` §2 for the frozen strict $C$ set
