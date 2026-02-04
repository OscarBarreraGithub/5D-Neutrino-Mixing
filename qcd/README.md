# QCD Running Coupling Module

Computes the strong coupling constant $\alpha_s(\mu)$ at arbitrary energy scales via numerical integration of the $\overline{\text{MS}}$ beta function. Primarily used to evaluate $\alpha_s$ at 1–10 TeV for neutron EDM calculations in the Randall-Sundrum model.

## Quick Start

```python
from qcd import alpha_s

# alpha_s at 3 TeV (3-loop, default)
a_s = alpha_s(3000.0)
print(f"alpha_s(3 TeV) = {a_s:.4f}")   # 0.0796

# alpha_s at multiple scales
from qcd import alpha_s_array
import numpy as np
scales = np.array([1000, 3000, 5000, 10000])
print(alpha_s_array(scales))
```

## Physics

### The beta function

The QCD coupling constant $\alpha_s = g^2/(4\pi)$ runs with the energy scale $\mu$ according to the renormalization group equation. Starting from the gauge coupling beta function

$$
\mu\frac{dg}{d\mu} = -\frac{\beta_0\, g^3}{16\pi^2} - \frac{\beta_1\, g^5}{(16\pi^2)^2} - \cdots
$$

and using $\alpha_s = g^2/(4\pi)$, one derives

$$
\mu^2 \frac{d\alpha_s}{d\mu^2} = -\frac{\alpha_s^2}{4\pi}\left[\beta_0 + \beta_1 \frac{\alpha_s}{4\pi} + \beta_2 \left(\frac{\alpha_s}{4\pi}\right)^{\!2} + \beta_3 \left(\frac{\alpha_s}{4\pi}\right)^{\!3} + \cdots\right]
$$

The key point is that the expansion parameter is $\alpha_s/(4\pi)$, _not_ $\alpha_s/(2\pi)$. This factor of 2 arises from $\alpha_s = g^2/(4\pi)$ together with the change of variable from $\ln\mu$ to $\ln\mu^2$.

### Beta function coefficients

The $\overline{\text{MS}}$ coefficients depend on $n_f$, the number of active quark flavors:

| Loop order | Coefficient | Formula |
|------------|-------------|---------|
| 1-loop | $\beta_0$ | $11 - \frac{2}{3}n_f$ |
| 2-loop | $\beta_1$ | $102 - \frac{38}{3}n_f$ |
| 3-loop | $\beta_2$ | $\frac{2857}{2} - \frac{5033}{18}n_f + \frac{325}{54}n_f^2$ |
| 4-loop | $\beta_3$ | van Ritbergen–Vermaseren–Larin (1997), involves $\zeta(3)$ |

For $n_f = 5$ (between $m_b$ and $m_t$): $\beta_0 = 23/3 \approx 7.67$.

For $n_f = 6$ (above $m_t$): $\beta_0 = 7$.

The positivity of $\beta_0$ for $n_f \le 16$ guarantees **asymptotic freedom**: $\alpha_s(\mu) \to 0$ as $\mu \to \infty$.

### Numerical integration

We define the evolution variable $t = \ln(\mu^2)$ so the RG equation becomes

$$
\frac{d\alpha_s}{dt} = \beta(\alpha_s)
$$

with no explicit $\mu$ on the right-hand side. This ODE is integrated numerically using `scipy.integrate.solve_ivp` (Runge–Kutta 4/5) from a reference point $\alpha_s(M_Z) = 0.1180$ (PDG 2024).

### Flavor thresholds

As $\mu$ crosses a quark mass threshold, the number of active flavors changes and the beta function coefficients shift discontinuously. For running from $M_Z$ to TeV scales, the only crossing is the top quark at $m_t = 172.69$ GeV:

$$
M_Z \;\xrightarrow{n_f = 5}\; m_t \;\xrightarrow{n_f = 6}\; \mu_{\text{target}}
$$

The integration is split into segments at each threshold. At leading order, $\alpha_s$ is continuous across thresholds (no matching correction). Higher-order decoupling corrections (Chetyrkin–Kuhn–Sturm 2006) are $\lesssim 0.1\%$ at the top threshold — well below the $\pm 0.8\%$ uncertainty on $\alpha_s(M_Z)$ itself.

### Reference values

| Scale | $\alpha_s$ (3-loop) |
|-------|:-------------------:|
| $M_Z = 91.19$ GeV | 0.1180 (input) |
| $m_t = 172.7$ GeV | 0.1076 |
| 1 TeV | 0.0885 |
| 3 TeV | 0.0796 |
| 5 TeV | 0.0761 |
| 10 TeV | 0.0718 |

At $\mu = 3$ TeV (a typical KK scale), $\alpha_s/\pi \approx 0.025$, so perturbative QCD corrections are well-controlled.

## API Reference

### `alpha_s(mu, n_loops=3, alpha_s_ref=0.1180, mu_ref=91.1876, thresholds=None)`

Compute $\alpha_s(\mu)$ by integrating the beta function from the reference scale.

- **`mu`** — target energy scale in GeV
- **`n_loops`** — loop order (1–4), default 3 (NNLO)
- **`alpha_s_ref`** — reference coupling value, default $\alpha_s(M_Z)$
- **`mu_ref`** — reference scale in GeV, default $M_Z$
- **`thresholds`** — list of `(mass, n_f_below, n_f_above)` tuples; pass `[]` to disable

### `alpha_s_array(mu_values, ...)`

Convenience wrapper that evaluates `alpha_s` at each element of an array.

### `beta_coefficients(n_f, n_loops)`

Returns the array $[\beta_0, \beta_1, \ldots]$ of length `n_loops` for a given $n_f$.

### `beta_rhs(alpha_s, n_f, n_loops)`

Evaluates the beta function RHS for use in ODE integration.

### `beta_0(n_f)` through `beta_3(n_f)`

Individual coefficient functions, useful for analytic cross-checks.

## Module Structure

```
qcd/
├── __init__.py        # Package exports
├── constants.py       # PDG 2024 values (alpha_s(M_Z), quark masses, thresholds)
├── beta_function.py   # Beta coefficients and ODE right-hand side
├── running.py         # Numerical RG integration (main logic)
├── alphaS.ipynb       # Demonstration notebook with plots
└── README.md          # This file
```

## References

- PDG 2024, "Quantum Chromodynamics" review — $\alpha_s(M_Z)$, quark masses, beta function formulas
- van Ritbergen, Vermaseren, Larin, Phys. Lett. B400 (1997) 379 [[hep-ph/9701390](https://arxiv.org/abs/hep-ph/9701390)] — 4-loop $\beta_3$
- Chetyrkin, Kuhn, Sturm, Eur. Phys. J. C48 (2006) 107 — threshold matching corrections
