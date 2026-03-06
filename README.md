# physics-test

A minimal, auditable toy framework that tests whether a small discrete rule can approximate dimensionless coupling constants across EM, strong, weak, and gravitational sectors.

## Core hypothesis

$$
G = \frac{C}{\phi^m}, \qquad F_0 = \phi^m \frac{k_B K}{h}, \qquad G \cdot F_0 = C \frac{k_B K}{h}
$$

With $\phi$ = golden ratio, $m \in \mathbb{Z}$, and $C$ restricted to a small gauge-derived menu: $C \in \{360, 180, 120, 60, 45, 15\}$.

## Documentation

All formulas, findings, anti-overfitting rules, results tables, and CLI commands are consolidated in:

- **[PHI_LATTICE_FRAMEWORK.md](PHI_LATTICE_FRAMEWORK.md)** -- the single reference document for the entire project.

Curated measurement inputs (values, uncertainties, schemes) live in:

- **[data/targets.json](data/targets.json)**

## Quickstart

```bash
pip install -r requirements.txt

python -m physics_test.cli check-example
python -m physics_test.cli list-gauge-Cs
python -m physics_test.cli scan-gauge-Cs --target "1/alpha" --max-rel-err 0.05
python -m physics_test.cli oos-predictive-rg --suite v1 --max-rel-err 0.02
python -m physics_test.cli list-targets
```

## Key results at a glance

| Sector | Target | Best $(C,m)$ | Rel. err |
|---|---|---:|---:|
| EM | `1/alpha` | (360, 2) | +0.34% |
| Strong | `1/alpha_s_1loop_from_mZ(mH)` | (60, 4) | -1.27% |
| Weak | `1/alpha2(alpha(mZ),sin2_on_shell)` | (120, 3) | -0.73% |
| Hypercharge | `1/alpha1_GUT(alpha(mZ),sin2_on_shell)` | (60, 0) | +0.58% |

RG-within-band running (no re-fitting of $m$) resolves all strong/EM/weak/hypercharge OOS misses to passes at 2%.

See [PHI_LATTICE_FRAMEWORK.md](PHI_LATTICE_FRAMEWORK.md) for the full picture: gravity extension, GW-band configurations, GUT diagnostics, EW mixing, and all CLI commands.
