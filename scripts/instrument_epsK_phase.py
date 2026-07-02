#!/usr/bin/env python
"""Instrumented anarchic RS draw: retain the complex KK-gluon down couplings so we
can measure the epsilon_K survival MECHANISM (magnitude suppression vs phase
cancellation) with the TRUE phase Phi = arg[(G_L)_12 (G_R)_12], not a proxy.

Reuses the exact draw + fit + coupling machinery of scripts/anarchic_bauer_s1.py
(so numbers match the production S1/S2 ensembles), but captures per draw:
  - GL_ij, GR_ij (complex) for ij in {12,13,23}   [down-sector KK-gluon couplings]
  - Phi12 = arg(GL12*GR12), |C4_12| = |GL12*GR12|/M_KK^2   [the eps_K driver]
  - ratio_eps_K, ratio_dm_K, ratio_B_{d,s}, ratio_D        [observables]
The point: eps_K ~ |C4| sin(Phi),  Delta m_K ~ |C4| cos(Phi).  A survivor is either
small |C4| (magnitude) or small sin(Phi) (phase).  Now we can tell which, exactly.

RESEARCH TOOL (run RS-FLAVOR-ALIGNMENT-2026-07). Not wired into production.
"""
from __future__ import annotations
import argparse, math
import numpy as np
import pandas as pd

from warpConfig.wavefuncs import f_IR
from scripts.anarchic_bauer_s1 import (
    _draw_bauer_matrix, _fn_c_values, SCENARIOS,
    DEFAULT_K_GEV, DEFAULT_V_GEV,
)
from scripts.run_rs_anarchy import (
    _ordered_svd, _build_kk_gluon_couplings, _load_pdg_targets,
    jarlskog_invariant,
)


# ---------------------------------------------------------------------------
# Draw ensembles
# ---------------------------------------------------------------------------
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


def draw_nelson_barr_yukawas(rng, *, y_min=0.1, y_max=3.0, rho_cp=1.0, eta_leak=0.0):
    """Nelson-Barr-like draw (Codex-B spec): real anarchic down magnitudes, CP
    routed through the up sector via a shared rank-one spurion on the left-doublet
    index.  eta_leak=0 => Y_d real => sin(Phi_12)=0 while |C_4| stays anarchic and
    CKM CP survives in Y_u (rho_cp~1)."""
    R_u = _draw_real_bauer_matrix(rng, y_min, y_max)
    R_d = _draw_real_bauer_matrix(rng, y_min, y_max)
    a_Q = _draw_real_bauer_vector(rng, y_min, y_max)   # shared CP-odd left direction
    b_u = _draw_real_bauer_vector(rng, y_min, y_max)
    b_d = _draw_real_bauer_vector(rng, y_min, y_max)
    S_u = _rank_one_like(a_Q, b_u, _rms(R_u))
    S_d = _rank_one_like(a_Q, b_d, _rms(R_d))
    Y_u = R_u.astype(np.complex128) + 1j * rho_cp * S_u
    Y_d = R_d.astype(np.complex128) + 1j * eta_leak * rho_cp * S_d
    return Y_u, Y_d
from quarkConstraints.deltaf2 import (
    evaluate_delta_f2_constraints, compute_delta_f2_wilsons,
    evaluate_delta_mk_with_running,
)


def _instrument_draw(Y_u, Y_d, f_Q, f_u, f_d, M_KK_GeV, xi_KK, targets,
                     mass_factor, ckm_factor, j_factor):
    """One draw -> dict with observables AND the complex down couplings."""
    D_Q, D_u, D_d = np.diag(f_Q), np.diag(f_u), np.diag(f_d)
    M_u = DEFAULT_V_GEV * D_Q @ Y_u @ D_u
    M_d = DEFAULT_V_GEV * D_Q @ Y_d @ D_d
    U_L_u, m_up, U_R_u = _ordered_svd(M_u)
    U_L_d, m_dn, U_R_d = _ordered_svd(M_d)
    ckm = U_L_u.conj().T @ U_L_d
    J = float(jarlskog_invariant(ckm))

    # PDG gate (same tolerances as production _eval_draw), now incl. Jarlskog J
    # so CP is actually reproduced (needed for the Nelson-Barr comparison).
    up_log = np.log(np.maximum(m_up, 1e-30) / targets["up_masses_GeV"])
    dn_log = np.log(np.maximum(m_dn, 1e-30) / targets["down_masses_GeV"])
    ckm_log = np.array([
        abs(math.log(max(abs(ckm[0, 1]), 1e-30) / targets["abs_V_us"])),
        abs(math.log(max(abs(ckm[1, 2]), 1e-30) / targets["abs_V_cb"])),
        abs(math.log(max(abs(ckm[0, 2]), 1e-30) / targets["abs_V_ub"])),
    ])
    j_log = abs(math.log(max(abs(J), 1e-30) / abs(targets["J"])))
    passes_pdg = bool(
        (np.abs(up_log) <= math.log(mass_factor)).all()
        and (np.abs(dn_log) <= math.log(mass_factor)).all()
        and (ckm_log <= math.log(ckm_factor)).all()
        and j_log <= math.log(j_factor)
    )

    couplings = _build_kk_gluon_couplings(
        M_KK=M_KK_GeV, xi_KK=xi_KK, f_Q=f_Q, f_u=f_u, f_d=f_d,
        U_L_u=U_L_u, U_L_d=U_L_d, U_R_u=U_R_u, U_R_d=U_R_d,
    )
    GL = np.asarray(couplings.left_down)    # 3x3 complex, g_s * U_Ld^dag F_Q^2 U_Ld
    GR = np.asarray(couplings.right_down)   # 3x3 complex, g_s * U_Rd^dag F_d^2 U_Rd
    GLu = np.asarray(couplings.left_up)     # up-sector, for D0/up CP (CP displacement)
    GRu = np.asarray(couplings.right_up)

    df2 = evaluate_delta_f2_constraints(couplings, M_KK=M_KK_GeV, xi_KK=xi_KK)
    bs = df2.by_system
    ratio_eps_K = float(bs["K"].ratio_to_bound)
    try:
        unev = compute_delta_f2_wilsons(couplings, M_KK=M_KK_GeV, xi_KK=xi_KK)
        w_k = next((w for w in unev if w.input.key == "epsilon_k"), None)
        ratio_dm_K = float(evaluate_delta_mk_with_running(w_k, mu_had=2.0).ratio_to_exp) if w_k else float("nan")
    except Exception:
        ratio_dm_K = float("nan")

    out = dict(passes_pdg=passes_pdg, J=J,
               ratio_eps_K=ratio_eps_K, ratio_dm_K=ratio_dm_K,
               ratio_B_d=float(bs["B_d"].ratio_to_bound),
               ratio_B_s=float(bs["B_s"].ratio_to_bound),
               ratio_D=float(bs["D"].ratio_to_bound))
    # complex down couplings for the 3 off-diagonal sectors
    for (i, j), tag in [((0, 1), "12"), ((0, 2), "13"), ((1, 2), "23")]:
        gl, gr = GL[i, j], GR[i, j]
        out[f"GL{tag}_abs"] = float(abs(gl)); out[f"GL{tag}_arg"] = float(np.angle(gl))
        out[f"GR{tag}_abs"] = float(abs(gr)); out[f"GR{tag}_arg"] = float(np.angle(gr))
        prod = gl * gr
        out[f"C4abs_{tag}"] = float(abs(prod) / M_KK_GeV**2)
        out[f"Phi_{tag}"] = float(np.angle(prod))           # total phase; eps ~ sin(Phi)
    # up-sector 1-2 (D0 / where CP is displaced under down-sector CP alignment)
    glu, gru = GLu[0, 1], GRu[0, 1]
    out["GLu12_abs"] = float(abs(glu)); out["GRu12_abs"] = float(abs(gru))
    out["C4Dabs_12"] = float(abs(glu * gru) / M_KK_GeV**2)
    out["PhiD_12"] = float(np.angle(glu * gru))

    # EDM CP invariant (Codex-A / Agashe-Perez-Soni): tree KK-gluon EDM vanishes
    # (diagonal couplings real); the RS neutron EDM comes from the one-loop
    # Higgs/KK-fermion dipole ~ Im[B_11/A_11] with the MISALIGNED spurion
    # B = F_Q (Y_u Y_u^d + Y_d Y_d^d) Y_d F_d.  I_d is rephasing-invariant; the
    # absolute d_n normalization is O(3)-parametric and is NOT applied here.
    FQ, Fd, Fu = np.diag(f_Q), np.diag(f_d), np.diag(f_u)
    YuYu = Y_u @ Y_u.conj().T
    YdYd = Y_d @ Y_d.conj().T
    A_d = U_L_d.conj().T @ FQ @ Y_d @ Fd @ U_R_d
    B_d = U_L_d.conj().T @ FQ @ (YuYu + YdYd) @ Y_d @ Fd @ U_R_d
    S_d = B_d[0, 0] / A_d[0, 0] if abs(A_d[0, 0]) > 0 else 0.0 + 0.0j
    out["edm_Id"] = float(np.imag(S_d))          # down-quark EDM CP invariant
    out["edm_Sd_abs"] = float(abs(S_d))
    A_u = U_L_u.conj().T @ FQ @ Y_u @ Fu @ U_R_u
    B_u = U_L_u.conj().T @ FQ @ (YuYu + YdYd) @ Y_u @ Fu @ U_R_u
    S_u = B_u[0, 0] / A_u[0, 0] if abs(A_u[0, 0]) > 0 else 0.0 + 0.0j
    out["edm_Iu"] = float(np.imag(S_u))          # up-quark EDM CP invariant
    return out


def run(scenario, M_KK_TeV_list, n_draws, seed, xi_KK=1.0,
        draw_mode="bauer", rho_cp=1.0, eta_leak=0.0):
    targets = _load_pdg_targets()
    sc = SCENARIOS[scenario]
    y_max = sc["y_max"]; c_max = sc["c_max"]
    common_cd = (scenario == "S2")
    rng = np.random.default_rng(seed)
    mass_factor, ckm_factor, j_factor = 3.0, 3.0, 10.0
    rows = []
    for M_KK_TeV in M_KK_TeV_list:
        M_KK_GeV = M_KK_TeV * 1000.0
        Lambda_IR = M_KK_GeV / xi_KK
        epsilon = (math.exp(-sc["L"]) if sc["L"] is not None else Lambda_IR / DEFAULT_K_GEV)
        c_lo_scan, c_hi_scan = -c_max, 0.5
        n_ok = 0
        for _ in range(n_draws):
            c_u3 = float(rng.uniform(c_lo_scan, c_hi_scan))
            try:
                if draw_mode == "nelson_barr":
                    Y_u, Y_d = draw_nelson_barr_yukawas(
                        rng, y_min=0.1, y_max=y_max, rho_cp=rho_cp, eta_leak=eta_leak)
                else:
                    Y_u = _draw_bauer_matrix(rng, 0.1, y_max)
                    Y_d = _draw_bauer_matrix(rng, 0.1, y_max)
                c_Q, c_u, c_d = _fn_c_values(c_u3, epsilon, targets, Y_u=Y_u, Y_d=Y_d,
                                             common_cd=common_cd, rng=rng, c_jitter=0.0)
                f_Q, f_u, f_d = f_IR(c_Q, epsilon), f_IR(c_u, epsilon), f_IR(c_d, epsilon)
                r = _instrument_draw(Y_u, Y_d, f_Q, f_u, f_d, M_KK_GeV, xi_KK, targets,
                                     mass_factor, ckm_factor, j_factor)
            except (ValueError, np.linalg.LinAlgError, FloatingPointError):
                continue
            except Exception:
                continue
            r["M_KK_TeV"] = M_KK_TeV; r["scenario"] = scenario
            r["draw_mode"] = draw_mode; r["eta_leak"] = float(eta_leak)
            rows.append(r); n_ok += 1
        print(f"  M_KK={M_KK_TeV} TeV: {n_ok}/{n_draws} evaluated", flush=True)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", default="S1")
    ap.add_argument("--mkk", default="1.5,2.0,3.0", help="comma TeV list")
    ap.add_argument("--ndraws", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--draw-mode", default="bauer", choices=["bauer", "nelson_barr"])
    ap.add_argument("--rho-cp", type=float, default=1.0)
    ap.add_argument("--eta-leak", type=float, default=0.0)
    ap.add_argument("--out", default=".orchestration/runs/RS-FLAVOR-ALIGNMENT-2026-07/instrument_S1.parquet")
    a = ap.parse_args()
    mkk = [float(x) for x in a.mkk.split(",")]
    print(f"instrumenting scenario={a.scenario} mode={a.draw_mode} "
          f"rho_cp={a.rho_cp} eta_leak={a.eta_leak} M_KK={mkk} n={a.ndraws}")
    df = run(a.scenario, mkk, a.ndraws, a.seed,
             draw_mode=a.draw_mode, rho_cp=a.rho_cp, eta_leak=a.eta_leak)
    df.to_parquet(a.out)
    print(f"wrote {len(df)} rows -> {a.out}")
