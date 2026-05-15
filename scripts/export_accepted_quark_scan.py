"""Export all accepted points from a merged quark scan to a flat CSV.

Reads `scan_outputs/<run>/merged/results.jsonl`, applies the publication
convention (physical KK-gluon mass `m_g^(1) = xi_kk * Lambda_IR` with
`xi_kk = GAUGE_KK_ROOT_NN ≈ 2.449`), keeps points that pass all five
Δ F = 2 systems under that convention, and writes a CSV with one row per
accepted point plus a sibling provenance JSON.

Because the ΔF=2 tree-level amplitudes scale as 1/M^2, the rescaling from the
scan's `xi_scan = 1` bookkeeping to `xi_kk ≈ 2.449` is exact: M_KK grows by
`xi_kk` and every ratio shrinks by `1/xi_kk^2`. The resulting CSV matches the
data underlying the publication figures in `results/figures/publication/`.

Usage
-----
python scripts/export_accepted_quark_scan.py \
    --run scan_outputs/dense_20260414T213617 \
    [--out-dir <run>/derived]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from quarkConstraints.scales import GAUGE_KK_ROOT_NN  # noqa: E402


SYSTEMS = ["epsilon_K", "K", "B_d", "B_s", "D0"]

CSV_COLUMNS = [
    "r",
    "m_gkk_TeV",
    "overall_scale",
    "binding_system",
    "max_ratio_to_bound",
    "ratio_epsilon_K",
    "ratio_K",
    "ratio_B_d",
    "ratio_B_s",
    "ratio_D0",
    "point_id",
]


def _binding_system(ratios: dict) -> str:
    return max(ratios.items(), key=lambda kv: kv[1])[0]


def _rescaled_row(point: dict, xi_kk: float) -> dict | None:
    """Rescale one JSONL row into publication convention and return a CSV row.

    Returns None if the point either has no ratio dict or fails the
    max-ratio ≤ 1 cut after rescaling.
    """
    raw_ratios = point.get("ratio_to_bound_by_system", {})
    if not raw_ratios:
        return None

    xi_scan = float(point.get("xi_KK", 1.0))
    axis_scale = xi_kk / xi_scan
    ratio_scale = 1.0 / (axis_scale ** 2)

    rescaled = {s: float(raw_ratios.get(s, 0.0)) * ratio_scale for s in SYSTEMS}
    max_ratio = max(rescaled.values())
    if max_ratio > 1.0:
        return None

    m_gkk_gev = float(point["M_KK"]) * axis_scale
    out = {
        "r": float(point["r"]),
        "m_gkk_TeV": m_gkk_gev / 1e3,
        "overall_scale": float(point["overall_scale"]),
        "binding_system": _binding_system(rescaled),
        "max_ratio_to_bound": max_ratio,
        "point_id": point["point_id"],
    }
    for sys_id in SYSTEMS:
        out[f"ratio_{sys_id}"] = rescaled[sys_id]
    return out


def _file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path,
                        help="Scan run directory (containing merged/results.jsonl)")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory (default: <run>/derived)")
    parser.add_argument("--xi-kk", type=float, default=GAUGE_KK_ROOT_NN,
                        help=f"Publication convention xi_kk (default: {GAUGE_KK_ROOT_NN:.6f})")
    args = parser.parse_args()

    jsonl_path = args.run / "merged" / "results.jsonl"
    manifest_path = args.run / "merged" / "manifest.json"
    config_path = args.run / "config.json"
    out_dir = args.out_dir or (args.run / "derived")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    total = 0
    converged = 0
    fit_rejected = 0
    with jsonl_path.open() as fh:
        for line in fh:
            total += 1
            pt = json.loads(line)
            if not pt.get("fit_success", pt.get("fit_converged", True)):
                continue
            if pt.get("fit_score", 0.0) > 0.1:
                fit_rejected += 1
                continue
            converged += 1
            row = _rescaled_row(pt, args.xi_kk)
            if row is not None:
                rows.append(row)

    rows.sort(key=lambda row: (row["r"], row["m_gkk_TeV"], row["overall_scale"]))

    csv_path = out_dir / "accepted_points.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    unique_r = sorted({row["r"] for row in rows})
    manifest = json.load(manifest_path.open()) if manifest_path.exists() else {}
    provenance = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_run_dir": str(args.run.resolve()),
        "source_results_jsonl": str(jsonl_path.resolve()),
        "source_results_sha256": _file_sha256(jsonl_path),
        "source_manifest_config_hash": manifest.get("config_hash"),
        "source_config_json": str(config_path.resolve()) if config_path.exists() else None,
        "total_points_in_jsonl": total,
        "points_passing_fit_filter": converged,
        "points_rejected_by_fit_score": fit_rejected,
        "accepted_points_publication_convention": len(rows),
        "unique_r_values": len(unique_r),
        "xi_kk": args.xi_kk,
        "xi_kk_source": "quarkConstraints.scales.GAUGE_KK_ROOT_NN"
                         if args.xi_kk == GAUGE_KK_ROOT_NN else "user-supplied",
        "rescaling": {
            "m_gkk_GeV": "Lambda_IR_scan * xi_kk",
            "ratio_<sys>": "ratio_scan_<sys> / xi_kk**2",
            "acceptance": "max(rescaled ratios over all 5 systems) <= 1",
        },
        "csv_columns": CSV_COLUMNS,
        "systems": SYSTEMS,
        "binding_system_definition": "argmax over systems of rescaled ratio_to_bound",
        "conventions_note": (
            "CSV values are in the publication convention: m_gkk is the physical "
            "first KK-gluon mode (not the scan's Lambda_IR bookkeeping), and all "
            "ratios are rescaled by 1/xi_kk**2 to match the modern DeltaF=2 "
            "figures in results/figures/publication/. Fit-quality filter "
            "(fit_score <= 0.1) matches the filter used in plot_publication_figures.py."
        ),
        "generator_script": "scripts/export_accepted_quark_scan.py",
    }
    (out_dir / "accepted_points.provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )

    print(f"Wrote {csv_path} ({len(rows)} rows, {len(unique_r)} unique r; "
          f"{converged} converged input points, {fit_rejected} fit-rejected)")


if __name__ == "__main__":
    main()
