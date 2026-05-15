"""Recompute per-quark mass + per-CKM-element residuals for an accepted-points CSV.

Re-runs the deterministic quark fit for each accepted point and writes a
sibling NPZ file containing per-row residual vectors. This isolates the
expensive refit step from notebook execution so the latter can run on a
low-resource node.

Usage
-----
    python scripts/compute_per_quark_residuals.py \
        --run scan_outputs/dense_20260506T141321 \
        --workers 32 \
        --out scan_outputs/dense_20260506T141321/derived/residuals_per_point.npz
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(var, "1")

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _refit_one(args):
    r, overall_scale, lambda_ir, k = args
    from quarkConstraints.benchmarks import default_quark_targets, default_spurion_seed
    from quarkConstraints.fit import fit_quark_sector

    sol = fit_quark_sector(
        default_quark_targets(),
        r=float(r), overall_scale=float(overall_scale), seed=default_spurion_seed(),
        k=float(k), Lambda_IR=float(lambda_ir), max_nfev=200, fit_orientation=True,
    )
    res = sol.result
    return (
        np.asarray(res.mass_residuals_up, dtype=float),
        np.asarray(res.mass_residuals_down, dtype=float),
        np.asarray(res.ckm_residuals, dtype=float),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path, help="run dir")
    parser.add_argument("--accepted-csv", type=Path, default=None,
                        help="defaults to <run>/derived/accepted_points.csv")
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--out", type=Path, default=None,
                        help="output NPZ path (default: <run>/derived/residuals_per_point.npz)")
    parser.add_argument("--xi-kk", type=float, default=2.4487, help="publication ξ_KK")
    parser.add_argument("--k", type=float, default=1.2209e19, help="AdS curvature k (GeV)")
    args = parser.parse_args()

    accepted_csv = args.accepted_csv or (args.run / "derived" / "accepted_points.csv")
    out_path = args.out or (args.run / "derived" / "residuals_per_point.npz")

    df = pd.read_csv(accepted_csv)
    print(f"loaded {len(df):,} accepted points from {accepted_csv}", flush=True)

    lambda_ir = (df["m_gkk_TeV"].astype(float).values * 1000.0) / float(args.xi_kk)
    items = list(zip(
        df["r"].astype(float).values,
        df["overall_scale"].astype(float).values,
        lambda_ir,
        [float(args.k)] * len(df),
    ))

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        out = list(pool.map(_refit_one, items, chunksize=64))
    dt = time.time() - t0
    print(f"re-fit {len(out):,} points in {dt:.1f} s using {args.workers} workers", flush=True)

    mres_up = np.array([r[0] for r in out])
    mres_dn = np.array([r[1] for r in out])
    ckm_res = np.array([r[2] for r in out])

    np.savez(
        out_path,
        point_id=df["point_id"].astype(str).values,
        r=df["r"].astype(float).values,
        overall_scale=df["overall_scale"].astype(float).values,
        m_gkk_TeV=df["m_gkk_TeV"].astype(float).values,
        mass_residuals_up=mres_up,
        mass_residuals_down=mres_dn,
        ckm_residuals=ckm_res,
    )
    print(f"wrote {out_path} ({out_path.stat().st_size / 1e6:.2f} MB)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
