"""Extract compact per-draw reproduction quantities from an RS-anarchy forward run.

The forward anarchy ensemble (``scripts/run_rs_anarchy.py``) draws COMPLEX O(1)
anarchic Yukawas + FIXED bulk c-values, computes the quark mass matrices, SVDs
them, forms the CKM, and evaluates all five Delta-F=2 systems FORWARD (no fit).
It keeps every draw and tags ``passes_pdg`` (masses+CKM within a factor of the
PDG targets). This is exactly the ACPS/Bauer/Casagrande "anarchic forward"
procedure that produces the multi-decade |eps_K| cloud.

This script STREAMS that run's ``draws.jsonl`` and writes a compact downsampled
table with the physical observables the literature figures need on their native
axes, plus BOTH a paper-era and a current experimental-input column where the
toggle matters (eps_K and the B-mixing windows). No re-scan is required: the
per-draw NP amplitudes are recovered from the stored ratios, which are linear in
the (fixed) experimental NP budget.

Recovered per-draw quantities
-----------------------------
* ``eps_k_np``    = ratio_eps_K * BUDGET_CENTRAL          [absolute |eps_K^NP|]
* ``dm_bd_np``    = ratio_B_d   * BOUND_Bd                 [GeV, |M12^NP|_d-like]
* ``dm_bs_np``    = ratio_B_s   * BOUND_Bs                 [GeV]
* ``dm_d_np``     = ratio_D     * BOUND_D                  [GeV]
* ``dm_k_np``     = ratio_Delta_m_K * (DELTA_M_K/2)        [GeV]

The stored ratio is ``NP_amplitude / bound``; multiplying back by the SAME bound
the run used recovers the bound-independent physical NP amplitude, which we then
re-compare against paper-era OR current windows at analysis time.

Usage
-----
    python scripts/anarchic_reproduction_extract.py \
        --draws scan_outputs/rs_anarchy_runA_20260515T085316/draws.jsonl \
        --per-tile 40000 \
        --out scan_outputs/anarchic_reproduction/anarchic_draws.parquet
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

# Bounds / budgets the forward run used (so ratio * bound recovers the NP amp).
# These are read straight from quarkConstraints.deltaf2 so the recovery is exact.
from quarkConstraints.deltaf2 import (  # noqa: E402
    EPSILON_K_EXP,
    EPSILON_K_SM,
    DELTA_M_K,
    DEFAULT_DELTA_F2_INPUTS_V1,
)

BUDGET_CENTRAL = abs(EPSILON_K_EXP - EPSILON_K_SM)  # ~6.7e-5
_BOUNDS = {i.key: i.bound for i in DEFAULT_DELTA_F2_INPUTS_V1}
BOUND_EPS_K = _BOUNDS["epsilon_k"]
BOUND_BD = _BOUNDS["b_d"]
BOUND_BS = _BOUNDS["b_s"]
BOUND_D = _BOUNDS["d"]
HALF_DM_K = DELTA_M_K / 2.0


def stream_extract(draws_path: Path, per_tile: int, seed: int = 12345) -> pd.DataFrame:
    """Reservoir-sample ``per_tile`` rows per M_KK tile and extract quantities."""
    rng = np.random.default_rng(seed)
    # Reservoir per tile to keep memory bounded over the 8M-row file.
    reservoirs: dict[float, list] = defaultdict(list)
    counts: dict[float, int] = defaultdict(int)
    n_total = 0
    with draws_path.open() as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not r.get("ok") or "deltaf2_ratios" not in r:
                continue
            n_total += 1
            mkk = float(r["M_KK_GeV"])
            counts[mkk] += 1
            res = reservoirs[mkk]
            if len(res) < per_tile:
                res.append(r)
            else:
                j = rng.integers(0, counts[mkk])
                if j < per_tile:
                    res[j] = r
            if n_total % 1_000_000 == 0:
                print(f"  streamed {n_total:,} rows...", flush=True)

    rows = []
    for mkk, res in reservoirs.items():
        for r in res:
            ratios = r["deltaf2_ratios"]
            passes = r.get("deltaf2_passes", {})
            rk = ratios.get("epsilon_K", np.nan)
            eps_k_np = rk * BUDGET_CENTRAL if np.isfinite(rk) else np.nan
            rows.append(
                {
                    "M_KK_GeV": mkk,
                    "M_KK_TeV": mkk / 1000.0,
                    "passes_pdg": bool(r.get("passes_pdg", False)),
                    # recovered absolute NP amplitudes
                    "eps_k_np": eps_k_np,
                    "ratio_eps_K": rk,
                    "ratio_dm_K": ratios.get("Delta_m_K", np.nan),
                    "ratio_B_d": ratios.get("Delta_m_Bd", np.nan),
                    "ratio_B_s": ratios.get("Delta_m_Bs", np.nan),
                    "ratio_D": ratios.get("Delta_m_D0", np.nan),
                    "dm_bd_np": ratios.get("Delta_m_Bd", np.nan) * BOUND_BD,
                    "dm_bs_np": ratios.get("Delta_m_Bs", np.nan) * BOUND_BS,
                    "dm_d_np": ratios.get("Delta_m_D0", np.nan) * BOUND_D,
                    "dm_k_np": ratios.get("Delta_m_K", np.nan) * HALF_DM_K,
                    # per-draw DeltaF2 pass flags (current bounds, central budget)
                    "pass_eps_K_current": bool(passes.get("epsilon_K", False)),
                    "pass_B_d": bool(passes.get("Delta_m_Bd", False)),
                    "pass_B_s": bool(passes.get("Delta_m_Bs", False)),
                    "pass_D": bool(passes.get("Delta_m_D0", False)),
                    # masses / CKM for diagnostics
                    "abs_V_us": r.get("abs_V_us", np.nan),
                    "abs_V_cb": r.get("abs_V_cb", np.nan),
                    "abs_V_ub": r.get("abs_V_ub", np.nan),
                    "J": r.get("J", np.nan),
                }
            )
    df = pd.DataFrame(rows)
    print(f"  total usable rows streamed: {n_total:,}; tiles: "
          f"{sorted(df['M_KK_TeV'].unique())}")
    print(f"  per-tile sampled counts:\n{df.groupby('M_KK_TeV').size()}")
    return df


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--draws", type=str,
                   default="scan_outputs/rs_anarchy_runA_20260515T085316/draws.jsonl")
    p.add_argument("--per-tile", type=int, default=40000)
    p.add_argument("--out", type=str,
                   default="scan_outputs/anarchic_reproduction/anarchic_draws.parquet")
    p.add_argument("--seed", type=int, default=12345)
    args = p.parse_args()

    draws_path = Path(args.draws)
    if not draws_path.is_absolute():
        draws_path = REPO / draws_path
    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = REPO / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"[extract] draws    = {draws_path}")
    print(f"[extract] per-tile = {args.per_tile}")
    print(f"[extract] BUDGET_CENTRAL (|eps_K^exp - eps_K^SM|) = {BUDGET_CENTRAL:.4e}")
    print(f"[extract] EPSILON_K_SM = {EPSILON_K_SM:.4e}, EXP = {EPSILON_K_EXP:.4e}")

    df = stream_extract(draws_path, args.per_tile, seed=args.seed)
    df.to_parquet(out_path, index=False)
    print(f"[extract] wrote {len(df):,} rows -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
