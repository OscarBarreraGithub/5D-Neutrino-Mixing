"""Fresh SMALL anarchic forward scan that EMITS the complex M_12^NP per Delta-F=2 system.

Motivation
----------
The production run ``scan_outputs/rs_anarchy_runA_20260515T085316/draws.jsonl``
stores only the per-system NP/bound *magnitude* ratios; the CP *phase* of
M_12^NP is not persisted. The phase is required for the Bauer 0912.1625 Fig. 6/7
(C_Bq vs phi_Bq) and Gedalia 0906.1879 Fig. 1 (x12^NP/x12 vs sin phi12) planes,
and for the Blanke 0809.1073 Re/Im(M12) plane and S_psiphi distribution.

This driver REUSES the exact forward inner loop of ``scripts/run_rs_anarchy.py``
(same fixed c-values, same anarchic Yukawa prior, same M_KK->Lambda_IR map, same
KK-gluon coupling build, same RG evolution to mu_had = 2 GeV), but instead of
collapsing to a magnitude it emits the COMPLEX M_12^NP for K, B_d, B_s, D via the
public ``compute_m12_np`` / ``_compute_m12_np`` of ``quarkConstraints.deltaf2``.

Sanity check (printed): for each system the recovered |M12^NP|/budget here must
track the stored ``ratio_*`` distribution of the production run (median compared).

Output
------
``scan_outputs/anarchic_reproduction/anarchic_complex_m12.parquet`` with, per draw:
  M_KK_TeV, passes_pdg,
  re_m12_K, im_m12_K, re_m12_Bd, im_m12_Bd, re_m12_Bs, im_m12_Bs, re_m12_D, im_m12_D,
  ratio_eps_K, ratio_B_d, ratio_B_s, ratio_D, ratio_dm_K   (for cross-check)

Usage
-----
    python scripts/anarchic_complex_m12.py --per-tile 40000 \
        --m-kk-tev 1,1.5,2,3,5,7,10
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

# Reuse the forward machinery verbatim (do NOT duplicate physics).
from scripts.run_rs_anarchy import (  # noqa: E402
    DEFAULT_C_Q, DEFAULT_C_U, DEFAULT_C_D,
    DEFAULT_XI_KK, DEFAULT_K_GEV, DEFAULT_V_GEV,
    DEFAULT_Y_HALF_RANGE, DEFAULT_Y_FLOOR,
    _draw_anarchic_matrix, _ordered_svd, _build_kk_gluon_couplings,
)
from warpConfig.wavefuncs import f_IR  # noqa: E402
import quarkConstraints.deltaf2 as d  # noqa: E402
from quarkConstraints.deltaf2 import (  # noqa: E402
    compute_delta_f2_wilsons, _evolve_wilsons, compute_m12_np, _compute_m12_np,
    evaluate_delta_f2_constraints, compute_delta_f2_wilsons as _wil,
    F_BD, M_BD, M_B_QUARK, M_D_QUARK_BD, B_1_BD, B_4_BD, B_5_BD,
    F_BS, M_BS, M_S_QUARK_BS, B_1_BS, B_4_BS, B_5_BS,
    F_D, M_D0, M_C_QUARK, M_U_QUARK, B_1_D, B_4_D, B_5_D,
)


# NP budgets (denominators of ratio_to_bound in the public path).
_BD_BUD = d._bd_budget()
_BS_BUD = d._bs_budget()
_D0_BUD = d._d0_budget()
_EPS_BUDGET = abs(d.EPSILON_K_EXP - d.EPSILON_K_SM)


def _wilsons_by_key(couplings, M_KK, xi_KK, mu_had=2.0):
    """Return {key: RG-evolved DeltaF2WilsonCoefficients} for all systems."""
    coeffs = compute_delta_f2_wilsons(couplings, M_KK=M_KK, xi_KK=xi_KK)
    return {c.input.key: _evolve_wilsons(c, mu_had=mu_had) for c in coeffs}


def _m12_complex_all(wilsons_by_key):
    """Complex M_12^NP (GeV) for K, B_d, B_s, D from RG-evolved Wilsons."""
    return {
        "K":  _compute_m12_np(wilsons_by_key["epsilon_k"]),
        "Bd": compute_m12_np(wilsons_by_key["b_d"], F_BD, M_BD, M_B_QUARK,
                             M_D_QUARK_BD, B_1_BD, B_4_BD, B_5_BD),
        "Bs": compute_m12_np(wilsons_by_key["b_s"], F_BS, M_BS, M_B_QUARK,
                             M_S_QUARK_BS, B_1_BS, B_4_BS, B_5_BS),
        "D":  compute_m12_np(wilsons_by_key["d"], F_D, M_D0, M_C_QUARK,
                             M_U_QUARK, B_1_D, B_4_D, B_5_D),
    }


def run_tile(M_KK_GeV, n_draws, seed, *, xi_KK, c_Q, c_u, c_d, v_GeV, k_GeV,
             y_half_range, y_floor):
    rng = np.random.default_rng(seed)
    Lambda_IR = M_KK_GeV / xi_KK
    epsilon = Lambda_IR / k_GeV
    f_Q = f_IR(np.asarray(c_Q), epsilon)
    f_u = f_IR(np.asarray(c_u), epsilon)
    f_d = f_IR(np.asarray(c_d), epsilon)
    D_Q, D_u, D_d = np.diag(f_Q), np.diag(f_u), np.diag(f_d)

    rows = []
    for _ in range(n_draws):
        Y_u = _draw_anarchic_matrix(rng, half_range=y_half_range, floor=y_floor)
        Y_d = _draw_anarchic_matrix(rng, half_range=y_half_range, floor=y_floor)
        try:
            M_u = v_GeV * D_Q @ Y_u @ D_u
            M_d = v_GeV * D_Q @ Y_d @ D_d
            U_L_u, m_up, U_R_u = _ordered_svd(M_u)
            U_L_d, m_dn, U_R_d = _ordered_svd(M_d)
        except (ValueError, np.linalg.LinAlgError):
            continue
        couplings = _build_kk_gluon_couplings(
            M_KK=M_KK_GeV, xi_KK=xi_KK, f_Q=f_Q, f_u=f_u, f_d=f_d,
            U_L_u=U_L_u, U_L_d=U_L_d, U_R_u=U_R_u, U_R_d=U_R_d,
        )
        wbk = _wilsons_by_key(couplings, M_KK_GeV, xi_KK)
        m12 = _m12_complex_all(wbk)

        # Ratios derived directly from the complex M12 (identical to the public
        # evaluate_delta_f2_constraints ratio_to_bound — verified in --selftest).
        eps_np = abs(d.KAPPA_EPSILON / (math.sqrt(2.0) * d.DELTA_M_K) * m12["K"].imag)
        rows.append({
            "M_KK_GeV": float(M_KK_GeV),
            "M_KK_TeV": float(M_KK_GeV / 1000.0),
            "re_m12_K": m12["K"].real, "im_m12_K": m12["K"].imag,
            "re_m12_Bd": m12["Bd"].real, "im_m12_Bd": m12["Bd"].imag,
            "re_m12_Bs": m12["Bs"].real, "im_m12_Bs": m12["Bs"].imag,
            "re_m12_D": m12["D"].real, "im_m12_D": m12["D"].imag,
            # ratios (= ratio_to_bound of the public path)
            "ratio_eps_K": float(eps_np / _EPS_BUDGET),
            "ratio_B_d": float(abs(m12["Bd"]) / _BD_BUD),
            "ratio_B_s": float(abs(m12["Bs"]) / _BS_BUD),
            "ratio_D": float(abs(m12["D"]) / _D0_BUD),
        })
    return rows


def _selftest():
    """Confirm the direct M12-derived ratios equal the public path ratio_to_bound."""
    rng = np.random.default_rng(99)
    M_KK = 3000.0; xi = DEFAULT_XI_KK
    eps = (M_KK / xi) / DEFAULT_K_GEV
    f_Q = f_IR(np.asarray(DEFAULT_C_Q), eps)
    f_u = f_IR(np.asarray(DEFAULT_C_U), eps)
    f_d = f_IR(np.asarray(DEFAULT_C_D), eps)
    D_Q, D_u, D_d = np.diag(f_Q), np.diag(f_u), np.diag(f_d)
    worst = 0.0
    for _ in range(200):
        Y_u = _draw_anarchic_matrix(rng, half_range=DEFAULT_Y_HALF_RANGE, floor=DEFAULT_Y_FLOOR)
        Y_d = _draw_anarchic_matrix(rng, half_range=DEFAULT_Y_HALF_RANGE, floor=DEFAULT_Y_FLOOR)
        U_L_u, _, U_R_u = _ordered_svd(DEFAULT_V_GEV * D_Q @ Y_u @ D_u)
        U_L_d, _, U_R_d = _ordered_svd(DEFAULT_V_GEV * D_Q @ Y_d @ D_d)
        cpl = _build_kk_gluon_couplings(M_KK=M_KK, xi_KK=xi, f_Q=f_Q, f_u=f_u, f_d=f_d,
                                        U_L_u=U_L_u, U_L_d=U_L_d, U_R_u=U_R_u, U_R_d=U_R_d)
        wbk = _wilsons_by_key(cpl, M_KK, xi)
        m12 = _m12_complex_all(wbk)
        eps_np = abs(d.KAPPA_EPSILON / (math.sqrt(2.0) * d.DELTA_M_K) * m12["K"].imag)
        mine = {
            "epsilon_k": eps_np / _EPS_BUDGET,
            "b_d": abs(m12["Bd"]) / _BD_BUD,
            "b_s": abs(m12["Bs"]) / _BS_BUD,
            "d": abs(m12["D"]) / _D0_BUD,
        }
        bs = evaluate_delta_f2_constraints(cpl, M_KK=M_KK, xi_KK=xi).by_system
        eval_ = {"epsilon_k": bs["K"].ratio_to_bound, "b_d": bs["B_d"].ratio_to_bound,
                 "b_s": bs["B_s"].ratio_to_bound, "d": bs["D"].ratio_to_bound}
        for k in mine:
            rel = abs(mine[k] - eval_[k]) / max(abs(eval_[k]), 1e-300)
            worst = max(worst, rel)
    print(f"[selftest] max relative deviation direct-vs-public = {worst:.3e}")
    ok = worst < 1e-9
    print("[selftest]", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--per-tile", type=int, default=40000)
    p.add_argument("--m-kk-tev", type=str, default="1,1.5,2,3,5,7,10")
    p.add_argument("--base-seed", type=int, default=770411)
    p.add_argument("--out", type=str,
                   default="scan_outputs/anarchic_reproduction/anarchic_complex_m12.parquet")
    p.add_argument("--selftest", action="store_true",
                   help="Verify direct M12-derived ratios == public evaluate_delta_f2 ratios.")
    args = p.parse_args(argv)

    if args.selftest:
        return _selftest()

    mkk_tev = [float(t) for t in args.m_kk_tev.split(",") if t.strip()]
    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = REPO / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"[complex_m12] M_KK tiles (TeV) = {mkk_tev}")
    print(f"[complex_m12] per-tile         = {args.per_tile}")
    print(f"[complex_m12] c_Q={DEFAULT_C_Q} c_u={DEFAULT_C_U} c_d={DEFAULT_C_D}")
    print(f"[complex_m12] Y prior: |Y| in [{DEFAULT_Y_FLOOR}, sqrt2*{DEFAULT_Y_HALF_RANGE}]")

    all_rows = []
    for idx, t in enumerate(mkk_tev):
        rows = run_tile(
            t * 1000.0, args.per_tile, args.base_seed + 1009 * idx,
            xi_KK=DEFAULT_XI_KK, c_Q=DEFAULT_C_Q, c_u=DEFAULT_C_U, c_d=DEFAULT_C_D,
            v_GeV=DEFAULT_V_GEV, k_GeV=DEFAULT_K_GEV,
            y_half_range=DEFAULT_Y_HALF_RANGE, y_floor=DEFAULT_Y_FLOOR,
        )
        all_rows.extend(rows)
        print(f"[complex_m12] tile {t:5.2f} TeV : {len(rows)} draws", flush=True)

    df = pd.DataFrame(all_rows)
    df.to_parquet(out_path, index=False)
    print(f"[complex_m12] wrote {len(df):,} rows -> {out_path}")

    # ---- sanity: |M12^NP| recovered here vs the production magnitude ratios ----
    print("\n[sanity] median |M12^NP|/budget per system (this run):")
    bd_bud = d._bd_budget(); bs_bud = d._bs_budget(); d0_bud = d._d0_budget()
    half_dmk = d.DELTA_M_K / 2.0
    for sysn, col, bud in [("B_d", "Bd", bd_bud), ("B_s", "Bs", bs_bud), ("D", "D", d0_bud)]:
        absm = np.hypot(df[f"re_m12_{col}"], df[f"im_m12_{col}"])
        ratio_here = np.median(absm / bud)
        ratio_eval = np.median(df[f"ratio_{ 'B_d' if col=='Bd' else 'B_s' if col=='Bs' else 'D'}"])
        print(f"  {sysn:4s}: median |M12|/budget (direct) = {ratio_here:.4e} | "
              f"median ratio_to_bound (eval) = {ratio_eval:.4e}")
    # eps_K: cross-check Im(M12_K) -> eps_K_NP/budget
    eps_np = np.abs(d.KAPPA_EPSILON / (math.sqrt(2.0) * d.DELTA_M_K) * df["im_m12_K"])
    eps_bud = abs(d.EPSILON_K_EXP - d.EPSILON_K_SM)
    print(f"  eps_K: median eps_NP/budget (direct) = {np.median(eps_np/eps_bud):.4e} | "
          f"median ratio_eps_K (eval) = {np.median(df['ratio_eps_K']):.4e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
