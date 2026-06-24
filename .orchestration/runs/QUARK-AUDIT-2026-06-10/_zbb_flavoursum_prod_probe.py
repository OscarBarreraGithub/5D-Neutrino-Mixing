"""INDEPENDENT production-floor impact probe for the Zbb flavour-sum (1-2c)->(1+2c) fix.

Replays REAL fitted scan draws through the full
  fit_quark_sector -> build_from_rs_ew_inputs
path and, at each point, compares the Casagrande fermion-KK Zbb shift and the
TOTAL z_delta_g_L/R_d[2,2] (what T010 consumes) computed with:
  - AFTER  : the live code  (denominator 1 + 2 c)   [agent's fix]
  - BEFORE : the old code   (denominator 1 - 2 c)   [pre-fix bug]
by monkeypatching build_rs_zbb_fermion_kk_mixing's denominator.

It also evaluates the T010 np_shift_ratio both ways at each point.

Read-only w.r.t. the repo source (monkeypatch is in-process only).
"""
from __future__ import annotations
import json, glob, sys
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from flavor_catalog_constraints import point_builder, registry
from quarkConstraints.benchmarks import default_quark_targets
from quarkConstraints.couplings import compute_quark_kk_gluon_couplings
from quarkConstraints.fit import QuarkFitSeed, fit_quark_sector
from quarkConstraints.rs_ew_spectrum import RSEWSpectrum, RSEWOverlapSplineCache
import scripts.run_full_catalog_scan as scan
import quarkConstraints.rs_ew_couplings as rsew

K_GEV = 1.2209e19
XI_KK = scan.DEFAULT_XI_KK
N_GAUGE = scan.DEFAULT_N_GAUGE_MODES
QUAD = scan.DEFAULT_QUADRATURE_ORDER
MINM, MAXM = 16, 512
RTOL = scan.DEFAULT_OVERLAP_RTOL
SPLG = scan.DEFAULT_SPLINE_GRID_SIZE
SPLV = scan.DEFAULT_SPLINE_VERIFY_POINTS
CMIN, CMAX = 0.3, 0.9
QFR = scan.DEFAULT_QUARK_FIT_R
QFMAX = scan.DEFAULT_QUARK_FIT_MAX_NFEV
EWM = "minimal_rs"


def seed_from_row(params):
    s = params["quark_yukawa_seed"]
    y_u = np.array(s["Y_u_re"]) + 1j*np.array(s["Y_u_im"])
    y_d = np.array(s["Y_d_re"]) + 1j*np.array(s["Y_d_im"])
    us, ul, ur = scan._svd_seed_parts(y_u)
    ds, dl, dr = scan._svd_seed_parts(y_d)
    return QuarkFitSeed(up_singular_values=us, down_singular_values=ds, overall_scale=1.0,
                        up_left=ul, up_right=ur, down_left=dl, down_right=dr)


def build_cache(lam):
    sp = RSEWSpectrum.build(lambda_ir_gev=lam, k_gev=K_GEV, n_gauge_modes=N_GAUGE,
                            quadrature_order=QUAD, model_label=EWM)
    ca = RSEWOverlapSplineCache.build(sp, c_min=CMIN, c_max=CMAX, grid_size=SPLG,
                                      include_omega=True, verify_points=SPLV, rel_tol=RTOL,
                                      min_modes=MINM, max_modes=MAXM)
    return sp, ca


def replay(params, mkk, lam, sp, ca):
    seed = seed_from_row(params)
    sol = fit_quark_sector(default_quark_targets(), r=QFR, seed=seed, Lambda_IR=lam,
                           k=K_GEV, max_nfev=QFMAX, fit_orientation=True)
    fit_result = sol.result
    rs_point = point_builder.build_from_rs_ew_inputs(
        fit_result, Lambda_IR=lam, k=K_GEV, n_gauge_modes=N_GAUGE, quadrature_order=QUAD,
        min_overlap_modes=MINM, max_overlap_modes=MAXM, overlap_rel_tol=RTOL, ew_model=EWM,
        spectrum=sp, rs_ew_cache=ca, **scan.QUARK_ONLY_BUILD_INCLUDE_FLAGS,
        lepton_yukawa_result=None, raw={"tile_id": 0, "draw_id": 0, "params": params})
    qc = compute_quark_kk_gluon_couplings(fit_result, M_KK=mkk, xi_KK=XI_KK, g_s_star=None)
    point = point_builder.make_point(raw=rs_point.raw, **dict(rs_point.extras),
                                      quark_mass_basis_couplings=qc, kk_gluon_mass_gev=float(mkk))
    res = scan._evaluate_constraint_ids(point, ("T010",))
    rsc = rs_point.extras["rs_ew_couplings"]
    meta = rsc.metadata.get("zbb_fermion_kk_mixing")
    return {
        "dgL_total": float(np.real(rsc.z_delta_g_L_d[2, 2])),
        "dgR_total": float(np.real(rsc.z_delta_g_R_d[2, 2])),
        "dgL_ferm": float(np.real(meta["delta_g_L_b"])) if meta else np.nan,
        "dgR_ferm": float(np.real(meta["delta_g_R_b"])) if meta else np.nan,
        "T010_ratio": float(res["T010"].ratio) if res["T010"].ratio is not None else np.nan,
        "T010_pass": bool(res["T010"].passes),
        # fitted Yukawa ratios driving the flavour sum
        "row_ratio": None, "col_ratio": None,
        "fit_result": fit_result,
    }


# ---- the monkeypatch machinery: force denominator sign ----
_ORIG = rsew.build_rs_zbb_fermion_kk_mixing
import math


def make_patched(sign):
    """sign=+1 -> 1+2c (live); sign=-1 -> 1-2c (old bug)."""
    def patched(quark_fit_result, *, spectrum):
        bs = rsew._required_attr(quark_fit_result, "bulk_state")
        c_q = rsew._real_triplet_required_attr(bs, "c_Q")
        c_d = rsew._real_triplet_required_attr(bs, "c_d")
        f_q = rsew._positive_real_triplet_required_attr(bs, "F_Q")
        f_d = rsew._positive_real_triplet_required_attr(bs, "F_d")
        y_d = rsew._complex_matrix_required_attr(bs, "Y_d_bulk_basis")
        masses_down = rsew._real_triplet_required_attr(quark_fit_result, "masses_down")
        m_b = float(masses_down[2])
        lambda_ir = rsew._positive_float(rsew._required_attr(spectrum, "lambda_ir_gev"), "spectrum.lambda_ir_gev")
        y33 = float(abs(y_d[2, 2])**2)
        profile_b_q = rsew._casagrande_zbb_B_profile_triplet(c_q, f_q, name="B_Q")
        profile_b_d = rsew._casagrande_zbb_B_profile_triplet(c_d, f_d, name="B_d")
        row_ratio = np.array(np.abs(y_d[2, :])**2 / y33, dtype=float)
        column_ratio = np.array(np.abs(y_d[:, 2])**2 / y33, dtype=float)
        f_d3_sq = float(f_d[2])**2
        f_q3_sq = float(f_q[2])**2
        ls_d = ls_q = 0.0
        for i in (0, 1):
            dd = 1.0 + sign*2.0*float(c_d[i])
            dq = 1.0 + sign*2.0*float(c_q[i])
            ls_d += float(row_ratio[i])/dd
            ls_q += float(column_ratio[i])/dq
        B_d = float(profile_b_d[2] + (1.0/(2.0*f_d3_sq))*ls_d)
        B_Q = float(profile_b_q[2] + (1.0/(2.0*f_q3_sq))*ls_q)
        pref = float(m_b*m_b/(2.0*lambda_ir*lambda_ir))
        obj = _ORIG(quark_fit_result, spectrum=spectrum)  # build the proper dataclass
        # overwrite the four fields via dataclasses.replace
        import dataclasses
        return dataclasses.replace(obj, B_Q=B_Q, B_d=B_d,
                                   delta_g_L_b=float(pref*B_d), delta_g_R_b=float(-pref*B_Q))
    return patched


def main():
    registry.discover()
    # representative FITTED scan points across M_KK tiles
    tile_glob = sorted(glob.glob(str(REPO/"scan_outputs/wq_quarkonly_20260604T195252/tile-*.jsonl")))
    per_tile = 6
    # pick a spread of tiles (low/mid/high M_KK)
    pick = [tile_glob[0], tile_glob[len(tile_glob)//3], tile_glob[2*len(tile_glob)//3], tile_glob[-1]]
    print(f"{'M_KK':>6} {'draw':>5} | {'dgL_f AFT':>11} {'dgL_f BEF':>11} {'ratioAB':>8} | "
          f"{'dgL_tot AFT':>12} {'dgL_tot BEF':>12} {'tot %chg':>9} | "
          f"{'T010 AFT':>10} {'T010 BEF':>10} {'pass A/B':>8}")
    print("-"*135)
    summary = []
    for tp in pick:
        with open(tp) as fh:
            draws = [json.loads(l) for l in fh][:60]
        draws = [d for d in draws if not d.get("skipped")][:per_tile]
        if not draws:
            continue
        lam = float(draws[0]["params"]["Lambda_IR"]); mkk = float(draws[0]["params"]["M_KK"])
        sp, ca = build_cache(lam)
        for d in draws:
            params = d["params"]
            try:
                rsew.build_rs_zbb_fermion_kk_mixing = make_patched(+1)
                aft = replay(params, mkk, lam, sp, ca)
                rsew.build_rs_zbb_fermion_kk_mixing = make_patched(-1)
                bef = replay(params, mkk, lam, sp, ca)
            except Exception as exc:
                print(f"{mkk/1000:6.1f} {d['draw_id']:5} | replay failed: {exc}")
                continue
            finally:
                rsew.build_rs_zbb_fermion_kk_mixing = _ORIG
            rab = (bef["dgL_ferm"]/aft["dgL_ferm"]) if aft["dgL_ferm"] else float("nan")
            totchg = 100.0*(bef["dgL_total"]-aft["dgL_total"])/abs(aft["dgL_total"]) if aft["dgL_total"] else float("nan")
            print(f"{mkk/1000:6.1f} {d['draw_id']:5} | {aft['dgL_ferm']:11.3e} {bef['dgL_ferm']:11.3e} {rab:8.2f} | "
                  f"{aft['dgL_total']:12.4e} {bef['dgL_total']:12.4e} {totchg:8.2f}% | "
                  f"{aft['T010_ratio']:10.4f} {bef['T010_ratio']:10.4f} {str(aft['T010_pass'])[0]}/{str(bef['T010_pass'])[0]:>5}")
            summary.append((mkk, aft, bef, totchg))
    # aggregate
    print("\n=== AGGREGATE ===")
    if summary:
        totchgs = np.array([s[3] for s in summary if np.isfinite(s[3])])
        t010_aft = np.array([s[1]["T010_ratio"] for s in summary if np.isfinite(s[1]["T010_ratio"])])
        t010_bef = np.array([s[2]["T010_ratio"] for s in summary if np.isfinite(s[2]["T010_ratio"])])
        print(f"total z_delta_g_L_d[2,2] % change BEFORE vs AFTER: "
              f"median {np.median(np.abs(totchgs)):.3f}%  max {np.max(np.abs(totchgs)):.3f}%")
        print(f"T010 ratio AFTER : median {np.median(t010_aft):.4f}  max {np.max(t010_aft):.4f}")
        print(f"T010 ratio BEFORE: median {np.median(t010_bef):.4f}  max {np.max(t010_bef):.4f}")
        flips = sum(1 for s in summary if s[1]["T010_pass"] != s[2]["T010_pass"])
        print(f"T010 pass/fail FLIPS between BEFORE and AFTER: {flips} / {len(summary)}")


if __name__ == "__main__":
    main()
