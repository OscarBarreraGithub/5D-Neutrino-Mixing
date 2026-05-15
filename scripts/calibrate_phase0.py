"""Phase-0 calibration sweep against the PDG-2024 quark target bundle.

Per plan v3 §7 / §8, Phase-0 is a deterministic dry scan that runs the new
PDG-2024 target objective with **acceptance gating disabled** so we can
record the per-quark log-residual + per-CKM relative-residual
distribution at every converged point. The 50/95/max percentiles inform
the production tolerance floors.

This script is the *machinery only*: the actual production scan covers
several million grid points and takes hours. Run with the user's preferred
grid, then commit the resulting JSON to ``data/phase0_residual_calibration.json``.

Example
-------
    python scripts/calibrate_phase0.py \
        --r-grid 0.05,0.1,0.2,0.3,0.5,0.8,1.2 \
        --overall-scale-grid 1.5,2.0,2.8,4.0,6.0 \
        --lambda-ir-grid 1000,3000,10000 \
        --output data/phase0_residual_calibration.json

Outputs JSON with the layout

    {
      "edition": "PDG 2024",
      "scale_GeV": 163.5,
      "n_points": <int>,
      "per_quark_residual_log": {
          "u": {"p50": ..., "p95": ..., "max": ...},
          ...
      },
      "per_ckm_residual_relative": {
          "Vus": {"p50": ..., "p95": ..., "max": ...},
          ...
      },
      "recommended_floor_per_quark": {"u": ..., ...},
      "raw": {"up": [[...], ...], "down": [[...], ...], "ckm": [[...], ...]}
    }
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List, Tuple

import numpy as np

from quarkConstraints.benchmarks import default_quark_targets
from quarkConstraints.pdg_quark_masses import (
    PDG_QUARK_MASSES_EDITION,
    pdg_2sigma_relative_at_scale,
)
from quarkConstraints.scales import DEFAULT_QUARK_FIT_SCALE_GEV
from quarkConstraints.scan import (
    MASS_TOLERANCE_FLOOR,
    QuarkScanConfig,
    run_quark_scan,
)


_UP_FLAVORS = ("u", "c", "t")
_DOWN_FLAVORS = ("d", "s", "b")
_CKM_LABELS = ("Vus", "Vcb", "Vub", "J")


def _percentiles(values: np.ndarray) -> dict[str, float]:
    if values.size == 0:
        return {"p50": float("nan"), "p95": float("nan"), "max": float("nan")}
    return {
        "p50": float(np.percentile(values, 50)),
        "p95": float(np.percentile(values, 95)),
        "max": float(np.max(values)),
    }


def _parse_grid(arg: str) -> np.ndarray:
    return np.array([float(x) for x in arg.split(",") if x.strip()], dtype=float)


def run(
    *,
    r_grid: np.ndarray,
    overall_scale_grid: np.ndarray,
    lambda_ir_grid: np.ndarray,
    output_path: Path,
    progress_every: int = 100,
) -> dict:
    config = QuarkScanConfig(
        r_values=r_grid,
        overall_scale_values=overall_scale_grid,
        Lambda_IR_values=lambda_ir_grid,
        record_git_metadata=False,
        targets=default_quark_targets(),
        apply_acceptance_gate=False,
    )
    rows = run_quark_scan(config, progress_every=progress_every)

    up_residuals: List[List[float]] = []
    down_residuals: List[List[float]] = []
    ckm_residuals: List[List[float]] = []
    for row in rows:
        if not bool(row.get("fit_success", False)):
            continue
        # rows store serialized arrays as ";"-joined strings; the in-memory
        # row also keeps the residuals via the QuarkFitResult, but here we
        # extract from the serialized columns to keep this script free of
        # internal-API coupling.
        up = np.array([float(x) for x in str(row["masses_up"]).split(";")])
        down = np.array([float(x) for x in str(row["masses_down"]).split(";")])
        targets = config.targets
        up_log = np.abs(np.log(np.maximum(up, 1e-30) / targets.up_masses))
        down_log = np.abs(np.log(np.maximum(down, 1e-30) / targets.down_masses))
        ckm = np.array([float(x) for x in str(row["ckm_observables"]).split(";")])
        ckm_rel = np.abs((ckm - targets.ckm_observables) / np.abs(targets.ckm_observables))
        up_residuals.append(up_log.tolist())
        down_residuals.append(down_log.tolist())
        ckm_residuals.append(ckm_rel.tolist())

    up_arr = np.array(up_residuals, dtype=float) if up_residuals else np.empty((0, 3))
    down_arr = np.array(down_residuals, dtype=float) if down_residuals else np.empty((0, 3))
    ckm_arr = np.array(ckm_residuals, dtype=float) if ckm_residuals else np.empty((0, 4))

    pdg_2sigma = pdg_2sigma_relative_at_scale(DEFAULT_QUARK_FIT_SCALE_GEV)

    per_quark: dict[str, dict[str, float]] = {}
    recommended_floor: dict[str, float] = {}
    for i, f in enumerate(_UP_FLAVORS):
        stats = _percentiles(up_arr[:, i] if up_arr.size else np.array([]))
        per_quark[f] = stats
        recommended_floor[f] = max(MASS_TOLERANCE_FLOOR, stats["p95"] if up_arr.size else 0.0)
    for i, f in enumerate(_DOWN_FLAVORS):
        stats = _percentiles(down_arr[:, i] if down_arr.size else np.array([]))
        per_quark[f] = stats
        recommended_floor[f] = max(MASS_TOLERANCE_FLOOR, stats["p95"] if down_arr.size else 0.0)

    per_ckm: dict[str, dict[str, float]] = {}
    for i, label in enumerate(_CKM_LABELS):
        per_ckm[label] = _percentiles(ckm_arr[:, i] if ckm_arr.size else np.array([]))

    payload = {
        "edition": PDG_QUARK_MASSES_EDITION,
        "scale_GeV": float(DEFAULT_QUARK_FIT_SCALE_GEV),
        "n_points": int(up_arr.shape[0]),
        "per_quark_residual_log": per_quark,
        "per_ckm_residual_relative": per_ckm,
        "pdg_2sigma_relative": {f: pdg_2sigma[f] for f in _UP_FLAVORS + _DOWN_FLAVORS},
        "recommended_floor_per_quark": recommended_floor,
        "deterministic_floor": MASS_TOLERANCE_FLOOR,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    return payload


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r-grid", type=str, default="0.05,0.1,0.2,0.4")
    parser.add_argument("--overall-scale-grid", type=str, default="2.8")
    parser.add_argument("--lambda-ir-grid", type=str, default="3000")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/phase0_residual_calibration.json"),
    )
    parser.add_argument("--progress-every", type=int, default=100)
    args = parser.parse_args(argv)
    payload = run(
        r_grid=_parse_grid(args.r_grid),
        overall_scale_grid=_parse_grid(args.overall_scale_grid),
        lambda_ir_grid=_parse_grid(args.lambda_ir_grid),
        output_path=args.output,
        progress_every=args.progress_every,
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover - script entry point
    raise SystemExit(main())
