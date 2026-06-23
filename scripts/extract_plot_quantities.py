#!/usr/bin/env python3
"""Replay stored scan draws through the FIXED working-tree code to extract the
per-draw physical quantities needed for the literature-comparison plots, and to
VALIDATE the replay against the stored constraint ratios.

For each sampled draw we reconstruct the exact draw (the stored anarchic seed
matrices Y_u, Y_d), rebuild the QuarkFitSeed via the scan's own SVD helper, run
the scan's own deterministic pipeline:

    fit_quark_sector -> build_from_rs_ew_inputs -> compute_quark_kk_gluon_couplings
    -> make_point -> evaluate(T010, T011, EW001, K001, C001, C002, B001, B003)

and read off:
  * delta_g_L_b, delta_g_R_b  (TOTAL bb Z-coupling shift = gauge-KK + fermion-KK
    admixture, the mass-basis [2,2] entry of z_delta_g_L/R_d that T010/R_b
    consumes; 0807.4937 Fig. 8 plane)
  * S, T                      (oblique proxy, 0807.4937 Fig. 4 plane)
  * |eps_K|, M12 D-mixing pieces, ratios for every constraint

VALIDATION GATE: the recomputed ratio_T010 and ratio_EW001 (from the replayed
adapters) are compared to the STORED ratios per draw.  If they do not match to
~1e-6 the script reports the mismatch loudly.

Usage:
    python scripts/extract_plot_quantities.py \
        --scan-dir scan_outputs/wq_quarkonly_20260622T090807 \
        --per-tile 250 --out scan_outputs/plot_quantities.parquet
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from flavor_catalog_constraints import point_builder, registry  # noqa: E402
from quarkConstraints.benchmarks import default_quark_targets  # noqa: E402
from quarkConstraints.couplings import compute_quark_kk_gluon_couplings  # noqa: E402
from quarkConstraints.fit import QuarkFitSeed, fit_quark_sector  # noqa: E402
from quarkConstraints.rs_ew_spectrum import RSEWSpectrum, RSEWOverlapSplineCache  # noqa: E402

# Reuse the scan's own helpers so the replay is byte-identical.
import scripts.run_full_catalog_scan as scan  # noqa: E402

# Constraint IDs whose ratio we recompute / extract.
EXTRACT_IDS = ("T010", "T011", "EW001", "K001", "C001", "C002", "B001", "B003")

# SM reference Zbb couplings used in Casagrande et al. 0807.4937 Fig. 8.
SM_GL_B = -0.42114
SM_GR_B = 0.077420

# Scan config knobs (must match the stored run_report config).
K_GEV = 1.2209e19
XI_KK = scan.DEFAULT_XI_KK
N_GAUGE_MODES = scan.DEFAULT_N_GAUGE_MODES
QUADRATURE_ORDER = scan.DEFAULT_QUADRATURE_ORDER
MIN_OVERLAP_MODES = 16
MAX_OVERLAP_MODES = 512
OVERLAP_RTOL = scan.DEFAULT_OVERLAP_RTOL
SPLINE_GRID = scan.DEFAULT_SPLINE_GRID_SIZE
SPLINE_VERIFY = scan.DEFAULT_SPLINE_VERIFY_POINTS
C_MIN, C_MAX = 0.3, 0.9
QUARK_FIT_R = scan.DEFAULT_QUARK_FIT_R
QUARK_FIT_MAX_NFEV = scan.DEFAULT_QUARK_FIT_MAX_NFEV
EW_MODEL = "minimal_rs"


def _seed_from_row(params: dict) -> QuarkFitSeed:
    """Rebuild the exact QuarkFitSeed from the stored anarchic seed matrices."""
    s = params["quark_yukawa_seed"]
    y_u = np.array(s["Y_u_re"], dtype=float) + 1j * np.array(s["Y_u_im"], dtype=float)
    y_d = np.array(s["Y_d_re"], dtype=float) + 1j * np.array(s["Y_d_im"], dtype=float)
    up_s, up_left, up_right = scan._svd_seed_parts(y_u)
    down_s, down_left, down_right = scan._svd_seed_parts(y_d)
    return QuarkFitSeed(
        up_singular_values=up_s,
        down_singular_values=down_s,
        overall_scale=1.0,
        up_left=up_left,
        up_right=up_right,
        down_left=down_left,
        down_right=down_right,
    )


def _build_spectrum_cache(lambda_ir_gev: float):
    spectrum = RSEWSpectrum.build(
        lambda_ir_gev=lambda_ir_gev,
        k_gev=K_GEV,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        model_label=EW_MODEL,
    )
    cache = RSEWOverlapSplineCache.build(
        spectrum,
        c_min=C_MIN,
        c_max=C_MAX,
        grid_size=SPLINE_GRID,
        include_omega=True,
        verify_points=SPLINE_VERIFY,
        rel_tol=OVERLAP_RTOL,
        min_modes=MIN_OVERLAP_MODES,
        max_modes=MAX_OVERLAP_MODES,
    )
    return spectrum, cache


def _replay_draw(params: dict, mkk_gev: float, lambda_ir_gev: float, spectrum, cache):
    """Replay one draw; return (extras_dict, results_dict)."""
    seed = _seed_from_row(params)
    sol = fit_quark_sector(
        default_quark_targets(),
        r=QUARK_FIT_R,
        seed=seed,
        Lambda_IR=lambda_ir_gev,
        k=K_GEV,
        max_nfev=QUARK_FIT_MAX_NFEV,
        fit_orientation=True,
    )
    fit_result = sol.result
    rs_point = point_builder.build_from_rs_ew_inputs(
        fit_result,
        Lambda_IR=lambda_ir_gev,
        k=K_GEV,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        min_overlap_modes=MIN_OVERLAP_MODES,
        max_overlap_modes=MAX_OVERLAP_MODES,
        overlap_rel_tol=OVERLAP_RTOL,
        ew_model=EW_MODEL,
        spectrum=spectrum,
        rs_ew_cache=cache,
        **scan.QUARK_ONLY_BUILD_INCLUDE_FLAGS,
        lepton_yukawa_result=None,
        raw={"tile_id": 0, "draw_id": 0, "params": params},
    )
    quark_couplings = compute_quark_kk_gluon_couplings(
        fit_result, M_KK=mkk_gev, xi_KK=XI_KK, g_s_star=None
    )
    point = point_builder.make_point(
        raw=rs_point.raw,
        **dict(rs_point.extras),
        quark_mass_basis_couplings=quark_couplings,
        kk_gluon_mass_gev=float(mkk_gev),
    )
    results = scan._evaluate_constraint_ids(point, EXTRACT_IDS)
    return rs_point.extras, results, fit_result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-dir", default="scan_outputs/wq_quarkonly_20260622T090807")
    ap.add_argument("--per-tile", type=int, default=250,
                    help="draws sampled per M_KK tile")
    ap.add_argument("--out", default="scan_outputs/plot_quantities.parquet")
    ap.add_argument("--max-tiles", type=int, default=None)
    args = ap.parse_args()

    registry.discover()

    scan_dir = (REPO / args.scan_dir) if not Path(args.scan_dir).is_absolute() else Path(args.scan_dir)
    tiles = sorted(scan_dir.glob("tile-*.jsonl"))
    if args.max_tiles is not None:
        tiles = tiles[: args.max_tiles]

    rows = []
    val_T010 = []
    val_EW001 = []
    for tile_path in tiles:
        # Read the first `per-tile` rows of this tile.
        sampled = []
        with tile_path.open() as fh:
            for i, line in enumerate(fh):
                if len(sampled) >= args.per_tile:
                    break
                r = json.loads(line)
                if r.get("skipped"):
                    continue
                sampled.append(r)
        if not sampled:
            continue
        lambda_ir = float(sampled[0]["params"]["Lambda_IR"])
        mkk = float(sampled[0]["params"]["M_KK"])
        spectrum, cache = _build_spectrum_cache(lambda_ir)
        print(f"[{tile_path.name}] M_KK={mkk/1000:.1f} TeV  Lambda_IR={lambda_ir:.1f}  "
              f"replaying {len(sampled)} draws ...", flush=True)
        for r in sampled:
            params = r["params"]
            try:
                extras, results, fit_result = _replay_draw(params, mkk, lambda_ir, spectrum, cache)
            except Exception as exc:  # noqa: BLE001
                print(f"  draw {r['draw_id']} replay failed: {exc}", flush=True)
                continue
            rsc = extras["rs_ew_couplings"]
            # TOTAL down-type bb Z-coupling shift that R_b / T010 actually sees:
            # the mass-basis [2,2] entry of z_delta_g_L_d / z_delta_g_R_d, which
            # already contains the universal/non-universal GAUGE-KK piece
            # (rs_ew_couplings.py:739) PLUS the Casagrande fermion-KK admixture
            # added on top (rs_ew_couplings.py:749). This is the exact object the
            # T010/Z-pole adapter consumes (T010.py:664-668,
            # _coupling_entry(rsc, "z_delta_g_L_d", 2, 2)). The previous code read
            # rsc.metadata["zbb_fermion_kk_mixing"], i.e. the fermion-only piece,
            # which omitted the dominant gauge contribution and so put the Zbb
            # stripe in the wrong place / wrong sign.
            dgL = float(np.real(rsc.z_delta_g_L_d[2, 2]))
            dgR = float(np.real(rsc.z_delta_g_R_d[2, 2]))
            # Fermion-only sub-piece kept as a diagnostic column (NOT the plotted
            # quantity) for the convention-note cross-check.
            _zbb_meta = rsc.metadata.get("zbb_fermion_kk_mixing")
            dgL_fermion = (
                float(np.real(_zbb_meta["delta_g_L_b"]))
                if _zbb_meta is not None else np.nan
            )
            dgR_fermion = (
                float(np.real(_zbb_meta["delta_g_R_b"]))
                if _zbb_meta is not None else np.nan
            )
            # results is already a dict[str, ConstraintResult].
            res = results
            row = {
                "tile_id": r["tile_id"],
                "draw_id": r["draw_id"],
                "M_KK_GeV": mkk,
                "Lambda_IR_GeV": lambda_ir,
                "delta_g_L_b": dgL,
                "delta_g_R_b": dgR,
                # Fermion-only sub-piece (diagnostic; NOT what R_b/T010 sees).
                "delta_g_L_b_fermion_only": dgL_fermion,
                "delta_g_R_b_fermion_only": dgR_fermion,
                # Absolute couplings on the 0807.4937 Fig. 8 SM reference point
                # (g_L^b_SM=-0.42114, g_R^b_SM=+0.077420).
                "g_L_b": SM_GL_B + dgL,
                "g_R_b": SM_GR_B + dgR,
            }
            # Stored vs recomputed ratios for the validation gate + the plots.
            for cid in EXTRACT_IDS:
                stored = r["constraints"].get(cid, {}).get("ratio")
                recomputed = float(res[cid].ratio) if cid in res and res[cid].ratio is not None else np.nan
                row[f"ratio_{cid}_stored"] = stored
                row[f"ratio_{cid}_recomp"] = recomputed
            # EW001 predicted S,T (from diagnostics) for the S-T plane.
            ew = res.get("EW001")
            if ew is not None and ew.diagnostics:
                diag = ew.diagnostics
                row["S_pred"] = float(diag.get("s_prediction", np.nan))
                row["T_pred"] = float(diag.get("t_prediction", np.nan))
                row["EW001_m_kk_gev"] = float(diag.get("m_kk_gev", np.nan))
                row["EW001_predicted_chi2"] = float(ew.predicted) if ew.predicted is not None else np.nan
            # K001 absolute |eps_K^NP|  (predicted is the absolute NP value).
            k001 = res.get("K001")
            if k001 is not None:
                row["epsK_NP_pred"] = float(k001.predicted) if k001.predicted is not None else np.nan
                kd = k001.diagnostics or {}
                if "im_m12_np_gev" in kd:
                    row["K001_im_m12_np_gev"] = float(np.real(kd["im_m12_np_gev"]))
            # C-mixing diagnostics for the D0 (Gedalia et al.) plot.
            # C001 carries |M12^NP| and the D-mixing inputs (x_D, y_D, Delta m_D).
            c001 = res.get("C001")
            if c001 is not None and c001.diagnostics:
                cd = c001.diagnostics
                row["C001_abs_m12_np_gev"] = float(np.real(cd.get("abs_m12_np_gev", np.nan)))
                row["C001_delta_m_d_exp_gev"] = float(cd.get("delta_m_d_experimental_gev", np.nan))
                row["C001_x_d_percent"] = float(cd.get("x_d_percent", np.nan))
                row["C001_y_d_percent"] = float(cd.get("y_d_percent", np.nan))
            # C002 carries the NP phase (sin 2 sigma_D proxy) and CP-odd fraction.
            c002 = res.get("C002")
            if c002 is not None and getattr(c002, "diagnostics", None):
                c2d = c002.diagnostics
                row["C002_abs_m12_np_gev"] = float(np.real(c2d.get("abs_m12_np_gev", np.nan)))
                row["C002_m12_np_phase_deg"] = float(np.real(c2d.get("m12_np_phase_degrees", np.nan)))
                row["C002_m12_np_phase_sine_abs"] = float(np.real(c2d.get("m12_np_phase_sine_abs", np.nan)))
                row["C002_cp_odd_fraction"] = float(np.real(c2d.get("cp_odd_fraction", np.nan)))
                row["C002_amp_ratio_to_dmd_half"] = float(
                    np.real(c2d.get("amplitude_ratio_to_delta_m_d_half", np.nan))
                )
            row["passes_T010"] = bool(res["T010"].passes) if "T010" in res else None
            row["passes_EW001"] = bool(res["EW001"].passes) if "EW001" in res else None
            row["survives_all_HARD_strict"] = bool(r.get("survives_all_HARD_strict"))
            row["survives_all_HARD_inclusive"] = bool(r.get("survives_all_HARD_inclusive"))
            rows.append(row)

            if row["ratio_T010_stored"] is not None and not np.isnan(row["ratio_T010_recomp"]):
                val_T010.append((row["ratio_T010_stored"], row["ratio_T010_recomp"]))
            if row["ratio_EW001_stored"] is not None and not np.isnan(row["ratio_EW001_recomp"]):
                val_EW001.append((row["ratio_EW001_stored"], row["ratio_EW001_recomp"]))

    df = pd.DataFrame(rows)
    out_path = (REPO / args.out) if not Path(args.out).is_absolute() else Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)
    print(f"\nWrote {len(df)} rows -> {out_path}")

    # ---- VALIDATION GATE ----
    def _report(name, pairs):
        if not pairs:
            print(f"VALIDATION {name}: NO PAIRS (cannot validate)")
            return None
        arr = np.array(pairs, dtype=float)
        stored, recomp = arr[:, 0], arr[:, 1]
        rel = np.abs(recomp - stored) / np.maximum(np.abs(stored), 1e-300)
        absd = np.abs(recomp - stored)
        worst = float(np.max(rel))
        med = float(np.median(rel))
        n_bad = int(np.sum(rel > 1e-6))
        print(f"VALIDATION {name}: n={len(pairs)}  median_rel={med:.2e}  "
              f"max_rel={worst:.2e}  max_abs={float(np.max(absd)):.2e}  "
              f"n(rel>1e-6)={n_bad}")
        return worst

    print("\n==== VALIDATION GATE (recomputed vs stored ratios) ====")
    w_t = _report("ratio_T010", val_T010)
    w_e = _report("ratio_EW001", val_EW001)
    ok = (w_t is not None and w_t < 1e-6) and (w_e is not None and w_e < 1e-6)
    print(f"\nGATE {'PASS' if ok else 'FAIL'} (threshold 1e-6)")
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
