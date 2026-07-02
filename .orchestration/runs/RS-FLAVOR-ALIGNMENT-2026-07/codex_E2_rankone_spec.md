**ANALYSIS/SPEC: Wrong-Ensemble Axis**

I did not write repo files or run long jobs.

**1. Implementable Ensemble**

Do not make this two independent `_draw_bauer_matrix(...)` calls. The U(2) selection rule needs a paired draw of `(Y_u,Y_d)` with a shared light-family structure. Also, a pure Yukawa replacement is not enough if `_fn_c_values` then regenerates split `c_{Q1},c_{Q2}` and `c_{d1},c_{d2}`: that silently breaks the U(2) protection. The correct scan mode is “paired Yukawa draw + U(2)-aware c fit”.

Use the bi-rank-one generalization
\[
Y_f=\sum_{a=1}^3 s_{fa}\, n^Q_{fa}\, n^{f\dagger}_{fa},\quad
s_f=(\epsilon'_f y_{f1},\epsilon_f y_{f2},y_{f3}).
\]
The user’s literal `n n^\dag` is the Hermitian special case `n^Q=n^f`; for quarks I would keep separate left/right vectors so RH-down alignment can be explicit.

```python
import math
import numpy as np

LAM, A = 0.225, 0.826

def _rot(i, j, theta, phi=0.0):
    r = np.eye(3, dtype=complex)
    c, s = math.cos(theta), math.sin(theta) * np.exp(-1j * phi)
    r[i, i] = r[j, j] = c
    r[i, j] = s
    r[j, i] = -s.conjugate()
    return r

def _haar_u2(rng):
    z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    q, rr = np.linalg.qr(z)
    q *= np.exp(-1j * np.angle(np.diag(rr)))
    u = np.eye(3, dtype=complex)
    u[:2, :2] = q
    return u

def _ckm_like(rng, cp="sequestered"):
    lam = LAM * rng.lognormal(0.0, 0.04)
    a = A * rng.lognormal(0.0, 0.08)
    eta = 0.35 if cp == "sequestered" else rng.uniform(-1.0, 1.0)
    rho = rng.normal(0.14, 0.08)
    th12, th23, th13 = lam, a * lam**2, a * lam**3 * math.hypot(rho, eta)
    delta = math.atan2(eta, rho)
    return _rot(1, 2, th23) @ _rot(0, 2, th13, delta) @ _rot(0, 1, th12)

def _rankone_yukawa(left_frame, right_frame, s):
    y = np.zeros((3, 3), dtype=complex)
    for a in range(3):
        y += s[a] * np.outer(left_frame[:, a], right_frame[:, a].conjugate())
    return y

def draw_rankone_u2_yukawas(rng, y_min=0.5, y_max=3.0,
                            cp="sequestered", rh_down_leak=0.0):
    # Eigenvalue hierarchy. Widths are intentional: masses still pass through SVD gate.
    yu3, yd3 = rng.uniform(y_min, y_max, 2)
    su = np.array([
        yu3 * (2.0e-5 * rng.lognormal(0, 0.5)),   # first: loop + mass suppressed
        yu3 * (7.0e-3 * rng.lognormal(0, 0.4)),   # second
        yu3,
    ])
    sd = np.array([
        yd3 * (1.2e-3 * rng.lognormal(0, 0.5)),
        yd3 * (2.2e-2 * rng.lognormal(0, 0.4)),
        yd3,
    ])

    # Shared U(2)_Q light plane; CKM comes from relative up/down left frames.
    q0 = _haar_u2(rng)
    q_u = q0
    q_d = q0 @ _ckm_like(rng, cp=cp)

    # RH frames: arbitrary light U(2) rotations are harmless if c_{f1}=c_{f2}.
    r_u = _haar_u2(rng) @ _rot(1, 2, LAM**2 * rng.uniform(0, 1)) @ _rot(0, 2, LAM**3 * rng.uniform(0, 1))
    r_d = _haar_u2(rng)
    if rh_down_leak:
        r_d = r_d @ _rot(1, 2, rh_down_leak * LAM**2) @ _rot(0, 2, rh_down_leak * LAM**3)

    Yu = _rankone_yukawa(q_u, r_u, su)
    Yd = _rankone_yukawa(q_d, r_d, sd)
    return Yu, Yd, {"s_u": su, "s_d": sd, "rh_down_leak": rh_down_leak}
```

Parameter count: flat Bauer anarchy samples 18 complex entries = 36 real. A literal Hermitian rank-one/projector pair has about `2*(3 eigenvalues + 6 orientation)=18` real. The U(2)-paired version above has roughly 6 singular values + 4 CKM-like left parameters + 4 up-RH leakage + 0-2 down-RH leakage + one sequestered CP phase, so about 14-17 active real parameters. Light-plane U(2) rotations can still be drawn as nuisance parameters, but exact light degeneracy makes them unphysical for KK-gluon `12`.

Masses are reproduced by the singular hierarchy plus profiles. CKM is reproduced by `q_u^\dag q_d`, not by random anarchic minors. The c-fit should enforce
`c_Q1=c_Q2`, `c_d1=c_d2` at leading order; first/second mass splittings then live in `s_{f1}/s_{f2}`, not in light-family profile splitting.

**2. Physical Prediction**

Flat anarchy gives the standard RS estimate
\[
(G_L^d)_{12}\sim g_s f_{Q1}f_{Q2},\quad
(G_R^d)_{12}\sim g_s f_{d1}f_{d2},
\]
so
\[
C_4^K\sim g_s^2\,{m_dm_s\over v^2Y_*^2M_{KK}^2}
\]
with an O(1) phase. This is the CFW/Bauer/Blanke epsilon_K wall: generic KK-gluon floors around 20 TeV in fully anarchic RS are repeatedly found in the literature. ([arxiv.org](https://arxiv.org/abs/0804.1954)) ([arxiv.org](https://arxiv.org/abs/0809.1073))

The U(2) rank-one ensemble changes the selection rule:
\[
F_Q^2=\mathrm{diag}(F_{Q\ell}^2,F_{Q\ell}^2,F_{Q3}^2),\quad
F_d^2=\mathrm{diag}(F_{d\ell}^2,F_{d\ell}^2,F_{d3}^2).
\]
Pure Cabibbo rotation inside the light doublet then gives no `12` KK-gluon coupling. The leading term needs two third-family spurion insertions:
\[
(G_L^d)_{12}\sim g_s(F_{Q3}^2-F_{Q\ell}^2)V_{td}^*V_{ts}\sim g_s F_{Q3}^2\lambda^5,
\]
and RH-down alignment gives
\[
(G_R^d)_{12}=0+\mathcal{O}(\delta_{d12},\theta^R_{13}\theta^R_{23}).
\]
If the only RH leakage is CKM-sized, the ratio to anarchic `GR12` is roughly
\[
{F_{d3}^2\lambda^5\over f_{d1}f_{d2}}
\sim {m_b^2\over m_dm_s}\lambda^{10}\sim 10^{-2}.
\]
So `C4` drops by about `10^-2` and the epsilon_K mass floor drops by `sqrt(10^-2)`, i.e. from ~20 TeV to ~2-3 TeV. With exact Santiago/FPR-style RH-down protection the leading LR operator is zero and the residual floor is set by U(2)-breaking leakage, radiative misalignment, or EW/collider limits, not by random phase luck. FPR and Santiago are the relevant alignment precedents. ([arxiv.org](https://arxiv.org/abs/0710.1869)) ([arxiv.org](https://arxiv.org/abs/0806.1230))

**3. Interaction With F1**

F1 says the flat-anarchy survivors at `M_KK<=3 TeV` are mostly rare small-`|C4|` points, with a smaller phase-aligned population. The rank-one/U(2) ensemble should move the whole cloud, not just select a tail.

Expected scatter:

- Flat anarchy: broad radial cloud in `(Delta m_K, epsilon_K)`; epsilon_K survivors either near the origin or in a phase strip with large `Delta m_K` but small imaginary part.
- Rank-one+U(2): both `epsilon_K` and `Delta m_K` contract toward the origin because `|C4|` is structurally small.
- Add CP sequester: `Phi_12` is no longer uniform; non-CKM phases collapse, but the important diagnostic is still the `C4abs_12` distribution shifting down before cuts.

So yes: it should convert F1’s “accidentally small magnitude” survivors into a “structurally small magnitude” ensemble. The check is whether `median(C4abs_12)` at 3 TeV is already below the flat-anarchy epsilon_K bound, not whether `|sin Phi_12|` is unusually small.

**4. Decisive Numerical Comparison**

Run one controlled comparison after implementing the ensemble:

\[
M_{50}^{\rm all}=\min M_{KK}\quad\text{such that}\quad
P(\text{PDG masses+CKM and all } \Delta F=2 \text{ cuts pass})\ge 50\%.
\]

Compare four lanes with identical scan budget and diagnostics:

`flat Bauer S1` vs `RH-down alignment only` vs `rank-one+U(2) only` vs `rank-one+U(2)+RH-down alignment+CP sequester`.

The synthesis wins only if the final lane has `M50_all ~ 2-3 TeV` while also showing:

- `median |C4_12|` suppressed by `>=10^-2` vs flat anarchy at fixed `M_KK`,
- phase-aligned fraction not needed, e.g. passes do not require `|sin Phi_12| < 0.1`,
- D-mixing and B-mixing do not become the new floor.

This beats anarchy by lowering the median floor, and beats pure alignment by showing the low floor survives realistic low-rank mass generation plus controlled CP rather than an exact down-sector limit.

**5. Novelty Statement**

Published ingredients exist: Greljo-Thomsen rank-one sequential flavor generation, U(2) light-family flavor as a middle ground between anarchy and MFV, FPR/Santiago/CPSW alignment in RS, Vecchi flavor branes, and Cheung-Fitzpatrick-Randall CP sequestering. ([arxiv.org](https://arxiv.org/abs/2309.11547)) ([arxiv.org](https://arxiv.org/abs/1105.2296)) ([arxiv.org](https://arxiv.org/abs/0907.0474)) ([arxiv.org](https://arxiv.org/abs/1206.4701)) ([arxiv.org](https://arxiv.org/abs/0711.4421))

What appears unpublished, and should be flagged as such, is the complete calculable package: sequential rank-one Yukawa spurions implemented as the RS flavor ensemble, an accidental U(2)^5 light-family bulk symmetry, RH-down minimal flavor protection, and CP sequestering, all tested numerically against the KK-gluon LR `C4` epsilon_K floor in the same generator. Low confidence only on the exact 2-3 TeV number until the U(2)-aware c-fit and radiative misalignment terms are included.