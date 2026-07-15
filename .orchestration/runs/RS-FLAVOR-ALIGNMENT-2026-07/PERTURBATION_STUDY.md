# Yukawa perturbation study — results and how to reproduce

**Run:** RS-FLAVOR-ALIGNMENT-2026-07. **Date:** 2026-07-15.
**Question (from the PI):** treat the failure of an epsilon_K survivor under a Yukawa
perturbation as data. Add noise to the entries, see how far you can move before it
fails and by how much; treat the response as a *linear transformation* on Yukawa
space that varies point to point; then find something interesting about that field of
linear transformations.

## What was done

Four surviving point classes at M_KK = 3 TeV:
- `flat_typical` — flat-anarchic survivor, magnitude-suppressed (small |C4|).
- `flat_tuned`   — flat-anarchic survivor, phase-aligned (large |C4|, theta_K~0). The
  "fine-tuned off-diagonal" point.
- `nelson_barr`  — real Y_d (CP in the up sector).
- `u2`           — rank-one/U(2) light-family-symmetric ensemble.

Three analyses (`scripts/yukawa_perturbation_study.py`, plots
`scripts/plot_perturbation_study.py`, data `perturb_study.npz`):

### A. Naive multiplicative noise  `Y_ij -> Y_ij*(1 + sigma * z_ij)`, z ~ N(0,1) per entry
Scan sigma; record eps_K pass-fraction AND median eps_K/bound (how badly it fails).
`fig_perturb_noise.png`.

| sigma | flat_tuned pass | flat_typical | u2 | nelson_barr |
|------:|---:|---:|---:|---:|
| 1e-3 | 88% | 100% | 100% | 100% |
| 3e-3 | 67% | 100% | 100% | 100% |
| 1e-2 | 29% | 99% | 100% | 100% |
| 3e-2 | 13% | 80% | 100% | 100% |
| 1e-1 |  2% | 55% | 85% | 100% |
| 3e-1 |  3% | 26% | 63% | 100% |

- **`flat_tuned` breaks at sub-percent noise** (drops below 50% by sigma~5e-3) and the
  failure is violent (median eps_K/bound reaches 36x). This confirms the PI's guess:
  the phase-cancellation survivor is destroyed by a ~0.1-1% perturbation.
- **`nelson_barr` never breaks** (pass = 100% at every sigma, median eps_K = 0 exactly).
  Reason: real multiplicative noise (z real) keeps Y_d REAL, so theta_K stays 0 by the
  symmetry. Real perturbations are exactly harmless.
- **`u2` is robust to ~10-30%**, `flat_typical` to ~3-10% (both magnitude-suppressed).

### B. Linear response: gradient  g = d(Im C4)/d(Y_d)  (18 real components)
The local linear operator. |g| = inverse tuning radius; direction = dangerous
perturbation. `fig_perturb_gradient.png`.

| class | \|grad\| | tuning radius delta_min = S_bound/\|grad\| |
|---|---:|---:|
| flat_typical | 1.4e-17 | 0.15 |
| flat_tuned   | 3.5e-16 | **7.9e-4** |
| u2           | 2.6e-17 | 0.024 |
| nelson_barr  | 3.5e-16 | (real dirs: infinite) |

- `flat_tuned` has the largest gradient and the smallest tuning radius (~0.08%): an
  isolated fragile point.
- Magnitude-suppressed points (`flat_typical`, `u2`) have small gradients / large radii.
- `nelson_barr`'s gradient is entirely in imaginary directions (top entries all `imY..`).

### C. The field of linear transformations: gradient direction + CP-purity
Collect the unit gradient over an ensemble per class; PCA it. `fig_perturb_pca.png`.
Plus the single scalar **CP-purity** = |grad_imag|^2 / |grad|^2 (fraction of the
sensitivity in CP-breaking/imaginary directions).

| class | PCA eff. dim (of 18) | CP-purity |
|---|---:|---:|
| flat_typical | 13.4 | 0.61 |
| nelson_barr  | 7.7 (**sharp cliff at 9**) | **1.000** |
| u2           | 10.7 | 0.874 |
| flat_tuned (point) | -- | 0.46 |

- **Nelson-Barr: the gradient field lives EXACTLY in a 9-dimensional subspace** (the 9
  imaginary components); PCA variance crashes to 1e-32 at component 10. So the real-Y_d
  half of Yukawa space (9 real dims) is an *exactly flat* protected manifold, and only
  the CP-breaking half is dangerous. CP-purity = 1.000 confirms it.
- **Anarchy is isotropic** (CP-purity ~0.5): real and imaginary perturbations are
  equally dangerous, so a survivor is an isolated point with no protected directions.
- **U(2)** sits between: CP-purity 0.87 (mostly protected; residual from third-family
  leakage), and its robustness comes from the small gradient magnitude (Part B).

## The interpretation (the "something interesting")

The epsilon_K constraint defines, at each Yukawa point, a **local linear operator**
g(Y) = grad(Im C4) — a covector on the 18-real-dim Y_d tangent space. Two invariants
of this field diagnose naturalness:
1. **|g(Y)|** = the inverse tuning radius. Large (fragile) at phase-tuned points, small
   (robust) where the magnitude is suppressed.
2. **the DIRECTION of g(Y)**, collected over an ensemble, is a fingerprint of the
   protecting symmetry: isotropic (fills Yukawa space) for anarchy, confined to the
   symmetry-breaking coset for a protected point. For Nelson-Barr it is *exactly* the
   9-dim CP-breaking subspace; the orthogonal 9-dim real subspace is an exactly flat
   protected manifold.

So "the survivors form a submanifold" becomes quantitative: the protected flat is the
kernel of the linear-response operator, and its dimension (9 for NB) and orientation
(CP-breaking vs real; light-family vs third-family) tell you exactly which symmetry to
gauge.

## How to reproduce
```
source ~/.bashrc && conda activate ising_bootstrap
export LD_LIBRARY_PATH=$HOME/.conda/envs/ising_bootstrap/lib:$LD_LIBRARY_PATH
export MPLBACKEND=Agg PYTHONPATH=$(pwd)
python scripts/yukawa_perturbation_study.py --part all --ndir 150 --nens 50 --seed 7
python scripts/plot_perturbation_study.py
```
Outputs: `perturb_study.npz`, `fig_perturb_noise.png`, `fig_perturb_gradient.png`,
`fig_perturb_pca.png` (all in this run dir).

## Continuation / next steps (what to do next)
- **Additive per-entry map (finer):** current gradient is `Im C4` only. Extend to the
  FULL constraint vector (masses, CKM, eps_K, dm_K, B_d, B_s, D) -> a Jacobian matrix
  per point; its singular-value spectrum gives the stiff/soft directions of the *whole*
  allowed region, not just eps_K. The protected flat is the intersection of the kernels.
- **U(2) exact cliff:** u2 did not show an exact PCA cliff (its protection is magnitude
  + leak, not an exact null at the Im-C4 level). Test whether the FULL Jacobian null
  space for u2 is the light-family coset (expected dimension countable from c_Q1=c_Q2).
- **Curvature:** the tuning radius uses the linear term; for the phase sheets the second
  order matters (budget-scaling exponent p in codex_C_tuning_spec.md). Add the Hessian.
- **Calibrated d_n:** replace the Im-C4 proxy with the APS EDM invariant field.
```
