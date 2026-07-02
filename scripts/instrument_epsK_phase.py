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
)
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

    # PDG gate (same tolerances as production _eval_draw)
    up_log = np.log(np.maximum(m_up, 1e-30) / targets["up_masses_GeV"])
    dn_log = np.log(np.maximum(m_dn, 1e-30) / targets["down_masses_GeV"])
    ckm_log = np.array([
        abs(math.log(max(abs(ckm[0, 1]), 1e-30) / targets["abs_V_us"])),
        abs(math.log(max(abs(ckm[1, 2]), 1e-30) / targets["abs_V_cb"])),
        abs(math.log(max(abs(ckm[0, 2]), 1e-30) / targets["abs_V_ub"])),
    ])
    passes_pdg = bool(
        (np.abs(up_log) <= math.log(mass_factor)).all()
        and (np.abs(dn_log) <= math.log(mass_factor)).all()
        and (ckm_log <= math.log(ckm_factor)).all()
    )

    couplings = _build_kk_gluon_couplings(
        M_KK=M_KK_GeV, xi_KK=xi_KK, f_Q=f_Q, f_u=f_u, f_d=f_d,
        U_L_u=U_L_u, U_L_d=U_L_d, U_R_u=U_R_u, U_R_d=U_R_d,
    )
    GL = np.asarray(couplings.left_down)    # 3x3 complex, g_s * U_Ld^dag F_Q^2 U_Ld
    GR = np.asarray(couplings.right_down)   # 3x3 complex, g_s * U_Rd^dag F_d^2 U_Rd

    df2 = evaluate_delta_f2_constraints(couplings, M_KK=M_KK_GeV, xi_KK=xi_KK)
    bs = df2.by_system
    ratio_eps_K = float(bs["K"].ratio_to_bound)
    try:
        unev = compute_delta_f2_wilsons(couplings, M_KK=M_KK_GeV, xi_KK=xi_KK)
        w_k = next((w for w in unev if w.input.key == "epsilon_k"), None)
        ratio_dm_K = float(evaluate_delta_mk_with_running(w_k, mu_had=2.0).ratio_to_exp) if w_k else float("nan")
    except Exception:
        ratio_dm_K = float("nan")

    out = dict(passes_pdg=passes_pdg,
               ratio_eps_K=ratio_eps_K, ratio_dm_K=ratio_dm_K,
               ratio_B_d=float(bs["B_d"].ratio_to_bound),
               ratio_B_s=float(bs["B_s"].ratio_to_bound),
               ratio_D=float(bs["D"].ratio_to_bound))
    # complex down couplings for the 3 off-diagonal sectors
    for (i, j), tag in [((0, 1), "12"), ((0, 2), "13"), ((1, 2), "23")]:
        gl, gr = GL[i, j], GR[i, j]
        out[f"GL{tag}_abs"] = float(abs(gl)); out[f"GL{tag}_arg"] = float(np.angle(gl))
        out[f"GR{tag}_abs"] = float(abs(gr)); out[f"GR{tag}_arg"] = float(np.angle(gr))
        # C4-like driver for this sector
        prod = gl * gr
        out[f"C4abs_{tag}"] = float(abs(prod) / M_KK_GeV**2)
        out[f"Phi_{tag}"] = float(np.angle(prod))           # total phase; eps ~ sin(Phi)
    return out


def run(scenario, M_KK_TeV_list, n_draws, seed, xi_KK=1.0):
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
            rows.append(r); n_ok += 1
        print(f"  M_KK={M_KK_TeV} TeV: {n_ok}/{n_draws} evaluated", flush=True)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", default="S1")
    ap.add_argument("--mkk", default="1.5,2.0,3.0", help="comma TeV list")
    ap.add_argument("--ndraws", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out", default=".orchestration/runs/RS-FLAVOR-ALIGNMENT-2026-07/instrument_S1.parquet")
    a = ap.parse_args()
    mkk = [float(x) for x in a.mkk.split(",")]
    print(f"instrumenting scenario={a.scenario} M_KK={mkk} n={a.ndraws}")
    df = run(a.scenario, mkk, a.ndraws, a.seed)
    df.to_parquet(a.out)
    print(f"wrote {len(df)} rows -> {a.out}")
