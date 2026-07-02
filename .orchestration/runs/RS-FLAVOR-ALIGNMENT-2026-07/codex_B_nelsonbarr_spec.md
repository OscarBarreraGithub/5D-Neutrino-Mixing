**Nelson-Barr Draw Spec**

Use an improved two-parameter NB draw. The CP spurion strength in the up sector is `rho_cp = O(1)`. The small parameter `eta_leak` controls symmetry-breaking leakage of the same CP spurion into the down sector.

Parametrization:

```text
Y_d = R_d + i eta_leak rho_cp S_d
Y_u = R_u + i rho_cp S_u
```

where `R_u,R_d` are real anarchic matrices with Bauer-like entry moduli, and `S_u,S_d` are rank-one real spurion textures. The strict down-real limit is `eta_leak = 0`; then `theta_K` is zero at tree level while CKM CP remains in `Y_u`.

```python
import numpy as np

def _draw_real_bauer_matrix(rng, y_min=0.1, y_max=3.0):
    mag = rng.uniform(y_min, y_max, size=(3, 3))
    sign = rng.choice(np.array([-1.0, 1.0]), size=(3, 3))
    return sign * mag

def _draw_real_bauer_vector(rng, y_min=0.1, y_max=3.0):
    mag = rng.uniform(y_min, y_max, size=3)
    sign = rng.choice(np.array([-1.0, 1.0]), size=3)
    return sign * mag

def _rms(x):
    return float(np.sqrt(np.mean(np.abs(x) ** 2)))

def _rank_one_like(a, b, target_rms):
    raw = np.outer(a, b)
    return raw * (target_rms / max(_rms(raw), 1e-12))

def draw_nelson_barr_yukawas(
    rng,
    *,
    y_min=0.1,
    y_max=3.0,
    rho_cp=1.0,
    eta_leak=0.0,
):
    R_u = _draw_real_bauer_matrix(rng, y_min, y_max)
    R_d = _draw_real_bauer_matrix(rng, y_min, y_max)

    # One CP-odd spurion direction shared through the left doublet index.
    a_Q = _draw_real_bauer_vector(rng, y_min, y_max)
    b_u = _draw_real_bauer_vector(rng, y_min, y_max)
    b_d = _draw_real_bauer_vector(rng, y_min, y_max)

    S_u = _rank_one_like(a_Q, b_u, _rms(R_u))
    S_d = _rank_one_like(a_Q, b_d, _rms(R_d))

    Y_u = R_u.astype(np.complex128) + 1j * rho_cp * S_u
    Y_d = R_d.astype(np.complex128) + 1j * eta_leak * rho_cp * S_d

    meta = {
        "rho_cp": float(rho_cp),
        "eta_leak": float(eta_leak),
        "Su_rms": _rms(S_u),
        "Sd_rms": _rms(S_d),
    }
    return Y_u, Y_d, meta
```

Do not use the literal minimal mode `Y_u = R_u + i eta a b^T` as the main NB test if `eta < 1e-2`: then generically `J(eta) = eta J_1 + O(eta^3)`, so CKM CP is also turned off. In the NB mode above,

```text
J(eta_leak, rho_cp) = J_0(rho_cp) + O(eta_leak)
```

so `rho_cp ~ 1` gives nonzero CKM Jarlskog while `eta_leak` controls only down-sector CP leakage.

**Predictions**

With `eta_leak = 0`, `M_d` is real, so `U_L^d,U_R^d` can be chosen real and `G_L^d,G_R^d` are real. Therefore `sin Phi_12 = 0` up to numerical noise, but `|G_L12 G_R12|` is anarchic.

For small leakage,

```text
theta_K,CP = asin(|sin Phi_12|) = O(eta_leak)
|C4_12| = |G_L12 G_R12| / M_KK^2 = anarchic * (1 + O(eta_leak^2))
R_K = 2 |M12_K^LR| / Delta m_K = anarchic * (1 + O(eta_leak^2))
```

Up-sector phases are not suppressed: `PhiD_12 = arg(G_L^u12 G_R^u12)` should remain broad/O(1). EDM-like down-sector CP invariants vanish at `eta_leak=0` and scale as `O(eta_leak)`.

**Test Protocol**

Modify [scripts/instrument_epsK_phase.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scripts/instrument_epsK_phase.py:1):

1. Add CLI args:
   `--draw-mode {bauer,nelson_barr}`, `--rho-cp`, `--eta-leak`.

2. In the draw loop:
   for `draw_mode="nelson_barr"`, replace `_draw_bauer_matrix` calls with `draw_nelson_barr_yukawas(...)`. Keep `_fn_c_values(... common_cd=False ...)` so magnitudes stay S1-like, not S2-like.

3. Add the missing J gate already present in `scripts/anarchic_bauer_s1.py`:
   record `J`, `abs_J_over_target`, `j_log`, and include `j_log <= log(j_factor)` in `passes_pdg`.

4. Add observables:
   `sinPhi12 = sin(Phi_12)`, `thetaK_folded = 0.5*angle(exp(2j*Phi_12))`;
   `GUL12_abs/arg`, `GUR12_abs/arg`, `C4Dabs_12`, `PhiD_12`, `sinPhiD12`;
   preferably also `R_K_LR = 2*abs(M12_K_LR)/DELTA_M_K`, not only full `ratio_dm_K`.

5. Scan at `M_KK=3 TeV`:
   `eta_leak = 0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2`, fixed `rho_cp=1`. If J acceptance is poor, scan `rho_cp = 0.3, 1, 3` and choose the value with stable PDG+J yield.

PASS criterion confirming the mechanism:

```text
On PDG+J-passing points at M_KK = 3 TeV:
eps_K pass fraction > 50% for eta_leak <= 1e-2
median R_K_LR within 20% of S1 complex-anarchy baseline
median |G_L12| and |G_R12| each within 20% of S1 baseline
median |sin Phi_12| scales linearly with eta_leak
median |sin PhiD_12| remains O(1), e.g. > 0.3
median |J|/J_PDG remains stable as eta_leak -> 0
```

**Distinct From S2/FPR**

Observable separation:

```text
NB:  |G_L12| ~ S1, |G_R12| ~ S1, R_K large, |sin Phi_12| << 1
S2:  |G_L12| ~ S1, |G_R12| suppressed, R_K small, phase not the main effect
FPR: |G_L12| suppressed, |G_R12| ~ S1, R_K small, phase not the main effect
```

So NB predicts epsilon_K survival with Delta m_K still large. S2/FPR survive by reducing the LR magnitude itself. This is the clean discriminator: plot `R_K_LR` versus `|sin Phi_12|` and color by `|G_L12|,|G_R12|`. NB occupies large-`R_K`, small-phase space; S2/FPR move to small-`R_K` magnitude-suppressed space.