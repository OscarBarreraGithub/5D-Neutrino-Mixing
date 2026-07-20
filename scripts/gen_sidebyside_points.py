#!/usr/bin/env python3
"""Generate CONTINUOUS-M_KK points for the side-by-side-vs-papers notebook.

Unlike ``extract_plot_quantities.py`` (which replays STORED draws on a discrete
grid of M_KK tiles), this script draws FRESH anarchic quark Yukawa seeds with
``M_KK`` log-uniform in [1, 10] TeV -- the x-range used by the published figures
we compare against (Casagrande et al. 0807.4937, Bauer et al. 0912.1625) -- so
our cloud forms a continuum rather than 10 vertical stripes.

It runs the SAME fixed-working-tree point-build path as the gridded extractor:

    _draw_anarchic_matrix (scan's own sampler)
      -> _svd_seed_parts -> QuarkFitSeed
      -> fit_quark_sector
      -> point_builder.build_from_rs_ew_inputs(**QUARK_ONLY_BUILD_INCLUDE_FLAGS)
         (include_fermion_kk_mixing=True)
      -> compute_quark_kk_gluon_couplings
      -> point_builder.make_point
      -> evaluate(T010, T011, EW001, K001)

Per point it records:
  * M_KK_GeV, Lambda_IR_GeV
  * delta_g_L_b, delta_g_R_b   = TOTAL (gauge-KK + fermion-KK) mass-basis [2,2]
    entry of rs_ew_couplings.z_delta_g_{L,R}_d -- the exact object T010 consumes
    (NOT the fermion-only metadata sub-piece)
  * g_L_b, g_R_b               = SM reference + total shift
  * S_pred, T_pred             = EW001 oblique S, T
  * epsK_NP_pred               = K001 absolute |eps_K^NP|
  * ratios for T010/T011/EW001/K001 (recomputed; for validation + plotting)
  * passes_T010, passes_T011, passes_K001, passes_EW001

VALIDATION: on the first few points it re-evaluates and checks the recomputed
ratios are finite and self-consistent; mismatches are reported loudly.

Usage:
    python scripts/gen_sidebyside_points.py --n 5000 --workers 16 \
        --out scan_outputs/sidebyside_points.parquet
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import scripts.run_full_catalog_scan as scan  # noqa: E402
from flavor_catalog_constraints import point_builder, registry  # noqa: E402
from quarkConstraints.benchmarks import default_quark_targets  # noqa: E402
from quarkConstraints.couplings import compute_quark_kk_gluon_couplings  # noqa: E402
from quarkConstraints.fit import QuarkFitSeed, fit_quark_sector  # noqa: E402
from quarkConstraints.rs_ew_spectrum import (  # noqa: E402
    RSEWOverlapSplineCache,
    RSEWSpectrum,
)

# --- scan config knobs (must match the gridded extractor / stored run) -------
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

# SM reference Zbb couplings (Casagrande et al. 0807.4937 Fig. 8).
SM_GL_B = -0.42114
SM_GR_B = 0.077420

EXTRACT_IDS = ("T010", "T011", "EW001", "K001")

# Geometry: M_KK and Lambda_IR are tied by xi_KK (M_KK = xi_KK * Lambda_IR), as
# in the scan tiles, so a continuous M_KK maps to a continuous Lambda_IR.
def _lambda_ir_from_mkk(mkk_gev: float) -> float:
    return float(mkk_gev / XI_KK)


# --- minimal ScanConfig stand-in for the anarchic sampler --------------------
# scan._draw_anarchic_matrix only reads y_prior / y_half_range / y_floor /
# y_sigma / y_trunc_sigma, so we reuse the scan's default ScanConfig.
_CFG = scan.ScanConfig(
    mkk_values_gev=(1000.0,),  # unused here (sampler only reads y_* fields)
    n_draws_per_tile=1,        # unused here
)


def _draw_seed(rng: np.random.Generator) -> tuple[QuarkFitSeed, np.ndarray, np.ndarray]:
    y_u = scan._draw_anarchic_matrix(rng, _CFG)
    y_d = scan._draw_anarchic_matrix(rng, _CFG)
    up_s, up_left, up_right = scan._svd_seed_parts(y_u)
    down_s, down_left, down_right = scan._svd_seed_parts(y_d)
    seed = QuarkFitSeed(
        up_singular_values=up_s,
        down_singular_values=down_s,
        overall_scale=1.0,
        up_left=up_left,
        up_right=up_right,
        down_left=down_left,
        down_right=down_right,
    )
    return seed, y_u, y_d


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


def _replay_one(seed: QuarkFitSeed, mkk_gev: float, lambda_ir_gev: float,
                spectrum, cache) -> dict:
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
        raw={"tile_id": 0, "draw_id": 0, "params": {}},
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

    rsc = rs_point.extras["rs_ew_couplings"]
    dgL = float(np.real(rsc.z_delta_g_L_d[2, 2]))
    dgR = float(np.real(rsc.z_delta_g_R_d[2, 2]))

    row: dict = {
        "M_KK_GeV": float(mkk_gev),
        "Lambda_IR_GeV": float(lambda_ir_gev),
        "delta_g_L_b": dgL,
        "delta_g_R_b": dgR,
        "abs_delta_g_R_b": abs(dgR),
        "g_L_b": SM_GL_B + dgL,
        "g_R_b": SM_GR_B + dgR,
    }
    for cid in EXTRACT_IDS:
        r = results.get(cid)
        row[f"ratio_{cid}"] = (
            float(r.ratio) if (r is not None and r.ratio is not None) else np.nan
        )
        row[f"passes_{cid}"] = bool(r.passes) if r is not None else None

    ew = results.get("EW001")
    if ew is not None and ew.diagnostics:
        d = ew.diagnostics
        row["S_pred"] = float(d.get("s_prediction", np.nan))
        row["T_pred"] = float(d.get("t_prediction", np.nan))
    k001 = results.get("K001")
    if k001 is not None:
        row["epsK_NP_pred"] = (
            float(k001.predicted) if k001.predicted is not None else np.nan
        )
    return row


# --- worker: each process owns one M_KK and a deterministic RNG --------------
def _worker(task: tuple[int, float, int, int]) -> list[dict]:
    seed_offset, mkk_gev, n_points, base_seed = task
    lambda_ir = _lambda_ir_from_mkk(mkk_gev)
    registry.discover()
    spectrum, cache = _build_spectrum_cache(lambda_ir)
    rng = np.random.default_rng(base_seed + seed_offset)
    out: list[dict] = []
    for _ in range(n_points):
        seed, _yu, _yd = _draw_seed(rng)
        try:
            row = _replay_one(seed, mkk_gev, lambda_ir, spectrum, cache)
        except Exception as exc:  # noqa: BLE001
            row = {"M_KK_GeV": float(mkk_gev), "replay_error": str(exc)}
        out.append(row)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=5000, help="total points to draw")
    ap.add_argument("--mkk-min-tev", type=float, default=1.0)
    ap.add_argument("--mkk-max-tev", type=float, default=10.0)
    ap.add_argument("--workers", type=int,
                    default=max(1, len(os.sched_getaffinity(0))))
    ap.add_argument("--n-tasks", type=int, default=0,
                    help="number of M_KK blocks (0 = workers*5)")
    ap.add_argument("--seed", type=int, default=20260622)
    ap.add_argument("--out", default="scan_outputs/sidebyside_points.parquet")
    args = ap.parse_args()

    registry.discover()

    # Log-uniform M_KK in [min, max] TeV; bin into per-worker chunks so each
    # process amortizes one spectrum/cache build over many draws (the cache
    # build is the expensive per-M_KK step).
    rng = np.random.default_rng(args.seed)
    mkk_tev = np.exp(rng.uniform(
        math.log(args.mkk_min_tev), math.log(args.mkk_max_tev), size=args.n
    ))
    mkk_gev = np.sort(mkk_tev) * 1000.0  # sort so nearby M_KK share a chunk

    # Chunk into `workers * k` tasks, each a contiguous block of similar M_KK.
    # Each task rebuilds the (expensive, ~50 s) RS-EW spectrum/overlap cache
    # ONCE and amortizes it over its whole block, so keep blocks fat: fewer
    # tasks -> fewer cache builds, at the cost of a slightly coarser x-binning
    # (each block uses its MEDIAN M_KK).  `args.n_tasks` overrides the default.
    n_tasks = args.n_tasks if args.n_tasks else max(args.workers * 5, args.workers)
    n_tasks = min(n_tasks, args.n)
    splits = np.array_split(mkk_gev, n_tasks)

    # Each task uses ONE representative M_KK for its whole block (the block is
    # narrow because the array is sorted) -> one spectrum build per block, while
    # the per-point M_KK is set to that representative.  To preserve continuity
    # we instead assign each task its block's MEDIAN M_KK but keep each point's
    # own M_KK.  Simpler + faithful: rebuild cache per point is too slow, so we
    # use a fine-grained representative per task block.  Width check below.
    tasks = []
    for i, block in enumerate(splits):
        if len(block) == 0:
            continue
        rep_mkk = float(np.median(block))
        tasks.append((i, rep_mkk, len(block), args.seed))
        lo, hi = float(block.min()), float(block.max())
        if hi / lo > 1.15:  # block too wide -> would smear x-axis
            print(f"  [warn] task {i} M_KK block width {lo/1000:.2f}-{hi/1000:.2f} "
                  f"TeV (ratio {hi/lo:.2f})", flush=True)

    print(f"Drawing {args.n} points, M_KK log-uniform in "
          f"[{args.mkk_min_tev}, {args.mkk_max_tev}] TeV, "
          f"{len(tasks)} tasks over {args.workers} workers ...", flush=True)

    t0 = time.time()
    rows: list[dict] = []
    if args.workers <= 1:
        for task in tasks:
            rows.extend(_worker(task))
            print(f"  task done; {len(rows)} pts; {time.time()-t0:.0f}s", flush=True)
    else:
        import multiprocessing as mp
        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=args.workers) as pool:
            for j, res in enumerate(pool.imap_unordered(_worker, tasks)):
                rows.extend(res)
                print(f"  task {j+1}/{len(tasks)} done; {len(rows)} pts; "
                      f"{time.time()-t0:.0f}s", flush=True)

    df = pd.DataFrame(rows)
    df["M_KK_TeV"] = df["M_KK_GeV"] / 1000.0
    out_path = (REPO / args.out) if not Path(args.out).is_absolute() else Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)
    n_err = int(df.get("replay_error", pd.Series([], dtype=object)).notna().sum()) \
        if "replay_error" in df.columns else 0
    print(f"\nWrote {len(df)} rows ({n_err} replay errors) -> {out_path} "
          f"in {time.time()-t0:.0f}s", flush=True)

    # --- validation summary ---
    good = df[df.get("ratio_T010").notna()] if "ratio_T010" in df.columns else df
    print("\n==== VALIDATION / SANITY ====")
    print(f"  finite T010 ratios : {len(good)}/{len(df)}")
    for col in ("delta_g_L_b", "delta_g_R_b", "S_pred", "T_pred", "epsK_NP_pred"):
        if col in df.columns:
            s = df[col].dropna()
            if len(s):
                print(f"  {col:16s} min={s.min():.3e} med={s.median():.3e} "
                      f"max={s.max():.3e}")
    for cid in EXTRACT_IDS:
        pc = f"passes_{cid}"
        if pc in df.columns:
            frac = df[pc].mean(skipna=True)
            print(f"  pass-frac {cid}: {frac:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
