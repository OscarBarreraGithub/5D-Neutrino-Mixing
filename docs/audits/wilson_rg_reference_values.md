# Wilson-Coefficient RG Reference Values

Audit date: 2026-05-15  
Code target: `quarkConstraints/deltaf2.py::_evolve_wilsons` and
`quarkConstraints/qcd_running.py` on branch `audit/wilson-rg`

## External references

- Buras, Jager, Urban, "Master Formulae for Delta F=2 NLO-QCD Factors in the
  Standard Model and Beyond", arXiv:hep-ph/0102316, Nucl. Phys. B605 (2001)
  600-624: https://arxiv.org/abs/hep-ph/0102316.  The arXiv abstract states
  that the formulae relate Wilson coefficients between a high scale and lower
  hadronic scales in the NDR scheme, and notes that the strongest RG effects
  occur in LR and scalar sectors.
- Buras, Misiak, Urban, "Two-Loop QCD Anomalous Dimensions of
  Flavour-Changing Four-Quark Operators Within and Beyond the Standard Model",
  arXiv:hep-ph/0005183, Nucl. Phys. B586 (2000) 397-426:
  https://arxiv.org/abs/hep-ph/0005183.  This is the two-loop ADM source
  behind the NLO formulae.
- FLAG Review 2024, arXiv:2411.04268:
  https://arxiv.org/abs/2411.04268.  FLAG reviews `B_K` and the additional BSM
  bag parameters for neutral-meson mixing and provides the lattice scale and
  scheme context used by the Phase 2 hole #5 bag audit.

## Operator basis extracted for this code path

BMU LR basis:

```text
Q1_LR^BMU = (bar h gamma_mu P_L q)(bar h gamma^mu P_R q)
Q2_LR^BMU = (bar h P_L q)(bar h P_R q)
```

Code basis:

```text
O1_VLL, O1_VRR, O4_LR, O5_LR
Q1_LR^BMU = -2 O5_LR
Q2_LR^BMU = O4_LR
C_BMU = (-C5_LR / 2, C4_LR)
```

The VLL/VRR entries are the BMU current-current operators.  The LR entries are
the conventional scalar `O4/O5` basis used with the code's `B4/B5` matrix
elements.

## LO anomalous dimensions used as reference

The Wilson coefficient evolution convention used here is

```text
C(mu_low) = U(mu_low, mu_high) C(mu_high)
eta = alpha_s(mu_high) / alpha_s(mu_low)
U_LO = V diag(eta ** (gamma_i^(0) / (2 beta0))) V^-1
```

Reference one-loop entries:

```text
gamma_VLL^(0) = 4

gamma_LR,BMU^(0) for [C1_LR^BMU, C2_LR^BMU] =
[[  2,   0],
 [ 12, -16]]

gamma_LR,code^(0) for [C4_LR, C5_LR] =
[[-16, -6],
 [  0,  2]]
```

The scalar-basis off-diagonal is therefore `gamma_45 = -6`, which maps
high-scale `C5_LR` into low-scale `C4_LR`.  The zero lower-left entry means high-scale
`C4_LR` does not generate `C5_LR` at LO in this basis.

## Threshold and alpha_s reference path

For the audit benchmark `mu_high = 3000 GeV` and `mu_low = 2 GeV`, the fixed
LO path is:

| Segment | Active flavours | alpha_s upper | alpha_s lower |
|---|---:|---:|---:|
| `3000 -> 163.5 GeV` | 6 | `0.0804135531797` | `0.108763928953` |
| `163.5 -> 4.18 GeV` | 5 | `0.108763928953` | `0.211846290829` |
| `4.18 -> 2 GeV` | 4 | `0.211846290829` | `0.267186117745` |

The charm threshold at `1.27 GeV` is below the kaon endpoint and is not crossed
for the default `2 GeV` run.

## Audited unit-vector evolution

The standalone audit script `scripts/audit_wilson_rg.py` compares the in-code
evolution with the closed-form LO reference above.  Its reference values for
`3 TeV -> 2 GeV` are:

| High-scale unit coefficient | Low-scale vector `(C1_VLL, C1_VRR, C4_LR, C5_LR)` |
|---|---|
| `C1_VLL = 1` | `(0.729130912171, 0, 0, 0)` |
| `C1_VRR = 1` | `(0, 0.729130912171, 0, 0)` |
| `C4_LR = 1` | `(0, 0, 3.53816397486, 0)` |
| `C5_LR = 1` | `(0, 0, 0.894757448992, 0.853891627884)` |

The maximum relative discrepancy between the code and the closed-form reference
is `1.300e-16`.

## Bag-parameter endpoint context

The Phase 2 hole #5 audit fixed the bag constants without changing RG
endpoints:

| System | Bag-parameter convention from hole #5 | Current code endpoint |
|---|---|---|
| Kaon `B_K` | MSbar/NDR at `2 GeV` | `2 GeV` |
| Kaon `B4`, `B5` | FLAG BSM averages quoted at `3 GeV` | `2 GeV` |
| `B_d`, `B_s` bags | MSbar/NDR at `mu = m_b` | `2 GeV` default |
| `D0` bags | comparable lattice inputs at `3 GeV` | `2 GeV` default |

The Wilson-RG audit does not change these hadronic endpoints.  Per-system
endpoint alignment should be a separate orchestrated change because it changes
all Delta F=2 scan outputs.
