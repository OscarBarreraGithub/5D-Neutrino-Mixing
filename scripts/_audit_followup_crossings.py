"""Audit M_KK^min crossings per follow-up run.

For each run dir: load draws.jsonl, restrict to passes_pdg, compute per-draw
M_KK_min_TeV = (M_KK_t / 1000) * sqrt(max_ratio), and report 50% and 95%
crossings in:
  - perturbative convention (g_s ~ 1.05) -- raw values
  - g_s* = 3 strong-coupling rescaling   -- multiply by sqrt((3/1.05)^2) = 2.857

Writes JSON to scan_outputs/followup_crossings_summary.json and prints a
plain-text table to stdout.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
GS_PERT = 1.05
GS_STAR = 3.0
GS_RESCALE = GS_STAR / GS_PERT  # = 2.857

RUNS = {
    "RUNA":              "scan_outputs/rs_anarchy_runA_20260507T065040",
    "RUN3_BASELINE":     "scan_outputs/rs_anarchy_run3_baseline_20260507T065529",
    "RUN3_QTOP_SHIFTED": "scan_outputs/rs_anarchy_run3_qtop_shifted_20260507T065529",
    "RUN3_MOREUV":       "scan_outputs/rs_anarchy_run3_moreUV_20260507T070052",
    "RUN3_MOREIR":       "scan_outputs/rs_anarchy_run3_moreIR_20260507T065523",
    "RUNB_NARROW":       "scan_outputs/rs_anarchy_runB_narrow_uniform_20260507T065416",
    "RUNB_WIDE":         "scan_outputs/rs_anarchy_runB_wide_uniform_20260507T065351",
    "RUNB_GAUSS":        "scan_outputs/rs_anarchy_runB_gaussian_3sigma_20260507T065141",
    "RUNC":              "scan_outputs/rs_anarchy_runC_20260507T065322",
}


def _load_pdg_passing_mkk_min(draws_path: Path) -> np.ndarray:
    out = []
    with draws_path.open() as fh:
        for line in fh:
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if not row.get("ok") or not row.get("passes_pdg"):
                continue
            mr = row.get("max_ratio")
            if mr is None or not np.isfinite(mr) or mr <= 0:
                continue
            mkk_t_tev = float(row["M_KK_GeV"]) / 1000.0
            out.append(mkk_t_tev * math.sqrt(float(mr)))
    return np.asarray(out)


def main():
    summary = {
        "convention": {
            "g_s_pert": GS_PERT,
            "g_s_star": GS_STAR,
            "rescale_factor": GS_RESCALE,
            "M_KK_units": "TeV",
            "M_KK_min_per_draw": "(M_KK_tile_GeV / 1000) * sqrt(max_ratio_pdg)",
        },
        "runs": {},
    }
    print("=== M_KK^min crossings per run (TeV) ===")
    print(f"{'run':22s}  {'n_pdg':>10s}  "
          f"{'p50_pert':>9s} {'p95_pert':>9s}  "
          f"{'p50_gs3':>9s} {'p95_gs3':>9s}")
    for tag, rd in RUNS.items():
        draws = REPO / rd / "draws.jsonl"
        if not draws.exists():
            print(f"{tag:22s}  MISSING: {draws}")
            summary["runs"][tag] = {"error": "draws.jsonl not found",
                                     "path": str(draws)}
            continue
        arr = _load_pdg_passing_mkk_min(draws)
        if arr.size == 0:
            print(f"{tag:22s}  {0:>10}     "
                  f"{'--':>9s} {'--':>9s}      {'--':>9s} {'--':>9s}")
            summary["runs"][tag] = {
                "n_pdg_pass": 0,
                "p50_pert_TeV": None,
                "p95_pert_TeV": None,
                "p50_gs_star_3_TeV": None,
                "p95_gs_star_3_TeV": None,
                "path": str(draws.parent),
            }
            continue
        p50p = float(np.percentile(arr, 50))
        p95p = float(np.percentile(arr, 95))
        p50s = p50p * GS_RESCALE
        p95s = p95p * GS_RESCALE
        print(f"{tag:22s}  {arr.size:>10,}    "
              f"{p50p:>9.2f} {p95p:>9.2f}    "
              f"{p50s:>9.2f} {p95s:>9.2f}")
        summary["runs"][tag] = {
            "n_pdg_pass": int(arr.size),
            "p50_pert_TeV": p50p,
            "p95_pert_TeV": p95p,
            "p50_gs_star_3_TeV": p50s,
            "p95_gs_star_3_TeV": p95s,
            "path": str(draws.parent),
        }

    out_json = REPO / "scan_outputs/followup_crossings_summary.json"
    out_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nwrote {out_json}")


if __name__ == "__main__":
    main()
