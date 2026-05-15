"""Export accepted quark-scan points with fitted bulk-basis Yukawas.

Reads the canonical derived ``accepted_points.csv`` for a scan run, preserves
its 11 publication-convention columns unchanged, re-runs the deterministic
quark-sector fit for each accepted point, and appends:

- row-major real/imaginary entries of ``Y_u_bulk_basis`` and
  ``Y_d_bulk_basis`` from ``quarkConstraints.fit.evaluate_quark_fit``;
- the diagonal bulk-mass parameters ``c_Q``, ``c_u``, and ``c_d`` in the
  bulk-mass eigenbasis.

The scan cache is inspected by convention but not used for the reconstruction:
the cache stores rounded reported warm-start seeds, not the raw optimizer point
whose ``QuarkFitResult`` was used in the original scan rows. Re-fitting from
the deterministic default seed reproduces the stored fit scores.

Usage
-----
python scripts/export_accepted_quark_scan_with_yukawas.py \
    --run scan_outputs/dense_20260414T213617 \
    --workers 32
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# Keep each worker to one BLAS/OpenMP thread; parallelism is across points.
for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from quarkConstraints.scales import GAUGE_KK_ROOT_NN  # noqa: E402


BASE_CSV_COLUMNS = [
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
SYSTEMS = ["epsilon_K", "K", "B_d", "B_s", "D0"]
MATRIX_PREFIXES = ("Y_u", "Y_d")
YUKAWA_COLUMNS = [
    f"{prefix}_{i}{j}_{part}"
    for prefix in MATRIX_PREFIXES
    for i in range(1, 4)
    for j in range(1, 4)
    for part in ("re", "im")
]
C_COLUMNS = [
    *(f"c_Q{i}" for i in range(1, 4)),
    *(f"c_u{i}" for i in range(1, 4)),
    *(f"c_d{i}" for i in range(1, 4)),
]
CSV_COLUMNS = [*BASE_CSV_COLUMNS, *YUKAWA_COLUMNS, *C_COLUMNS]


@dataclass(frozen=True)
class WorkItem:
    point_id: str
    r: float
    overall_scale: float
    lambda_ir: float
    k: float
    expected_fit_score: float
    max_nfev: int
    fit_orientation: bool


@dataclass(frozen=True)
class WorkerResult:
    point_id: str
    values: dict[str, str]
    fit_success: bool
    fit_score: float
    expected_fit_score: float
    fit_score_abs_delta: float
    nfev: int
    max_abs_log_mass_residual: float
    max_abs_ckm_observable_residual: float
    c_min: float
    c_max: float


def _file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _format_float(value: float) -> str:
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ValueError(f"non-finite numeric value {value!r}")
    return format(numeric, ".17g")


def _flat_complex_values(prefix: str, matrix: np.ndarray) -> dict[str, str]:
    arr = np.asarray(matrix, dtype=np.complex128)
    if arr.shape != (3, 3):
        raise ValueError(f"{prefix} matrix has shape {arr.shape}, expected (3, 3)")
    out: dict[str, str] = {}
    for i in range(3):
        for j in range(3):
            out[f"{prefix}_{i + 1}{j + 1}_re"] = _format_float(arr[i, j].real)
            out[f"{prefix}_{i + 1}{j + 1}_im"] = _format_float(arr[i, j].imag)
    return out


def _solve_one(item: WorkItem) -> WorkerResult:
    from quarkConstraints.benchmarks import default_quark_targets, default_spurion_seed
    from quarkConstraints.fit import fit_quark_sector

    solution = fit_quark_sector(
        default_quark_targets(),
        r=item.r,
        overall_scale=item.overall_scale,
        seed=default_spurion_seed(),
        k=item.k,
        Lambda_IR=item.lambda_ir,
        max_nfev=item.max_nfev,
        fit_orientation=item.fit_orientation,
    )
    result = solution.result
    state = result.bulk_state

    values: dict[str, str] = {}
    values.update(_flat_complex_values("Y_u", state.Y_u_bulk_basis))
    values.update(_flat_complex_values("Y_d", state.Y_d_bulk_basis))
    for prefix, arr in (
        ("c_Q", state.c_Q),
        ("c_u", state.c_u),
        ("c_d", state.c_d),
    ):
        c_values = np.asarray(arr, dtype=float)
        if c_values.shape != (3,):
            raise ValueError(f"{prefix} has shape {c_values.shape}, expected (3,)")
        for idx, value in enumerate(c_values, start=1):
            values[f"{prefix}{idx}"] = _format_float(value)

    mass_residuals = np.concatenate(
        [result.mass_residuals_up, result.mass_residuals_down]
    )
    ckm_residuals = result.ckm_residuals
    c_all = np.concatenate([state.c_Q, state.c_u, state.c_d])
    fit_score = float(result.score)
    return WorkerResult(
        point_id=item.point_id,
        values=values,
        fit_success=bool(solution.success),
        fit_score=fit_score,
        expected_fit_score=float(item.expected_fit_score),
        fit_score_abs_delta=abs(fit_score - float(item.expected_fit_score)),
        nfev=int(solution.nfev),
        max_abs_log_mass_residual=float(np.max(np.abs(mass_residuals))),
        max_abs_ckm_observable_residual=float(np.max(np.abs(ckm_residuals))),
        c_min=float(np.min(c_all)),
        c_max=float(np.max(c_all)),
    )


def _read_accepted_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != BASE_CSV_COLUMNS:
            raise ValueError(
                f"{path} schema mismatch: expected {BASE_CSV_COLUMNS!r}, "
                f"found {reader.fieldnames!r}"
            )
        return list(reader)


def _load_jsonl_metadata(
    jsonl_path: Path,
    accepted_point_ids: set[str],
) -> tuple[dict[str, dict[str, Any]], int, int, int]:
    metadata: dict[str, dict[str, Any]] = {}
    total = 0
    points_passing_fit_filter = 0
    points_rejected_by_fit_score = 0
    with jsonl_path.open() as fh:
        for line in fh:
            total += 1
            point = json.loads(line)
            if point.get("fit_success", point.get("fit_converged", True)):
                if float(point.get("fit_score", 0.0)) <= 0.1:
                    points_passing_fit_filter += 1
                else:
                    points_rejected_by_fit_score += 1
            if point["point_id"] in accepted_point_ids:
                metadata[point["point_id"]] = point
    missing = sorted(accepted_point_ids - metadata.keys())
    if missing:
        raise ValueError(f"{len(missing)} accepted point_ids missing from JSONL; first={missing[0]}")
    return metadata, total, points_passing_fit_filter, points_rejected_by_fit_score


def _read_scan_options(config_path: Path) -> tuple[int, bool]:
    if not config_path.exists():
        return 120, True
    payload = json.loads(config_path.read_text())
    config = payload.get("config", payload)
    return int(config.get("max_nfev", 120)), bool(config.get("fit_orientation", True))


def _build_work_items(
    rows: list[dict[str, str]],
    metadata: dict[str, dict[str, Any]],
    *,
    max_nfev: int,
    fit_orientation: bool,
) -> list[WorkItem]:
    items: list[WorkItem] = []
    for row in rows:
        point = metadata[row["point_id"]]
        items.append(
            WorkItem(
                point_id=row["point_id"],
                r=float(point["r"]),
                overall_scale=float(point["overall_scale"]),
                lambda_ir=float(point["Lambda_IR"]),
                k=float(point["k"]),
                expected_fit_score=float(point["fit_score"]),
                max_nfev=max_nfev,
                fit_orientation=fit_orientation,
            )
        )
    return items


def _cache_summary(run_dir: Path) -> dict[str, Any]:
    cache_files = list(run_dir.glob("shards/shard-*-of-*/cache/*.seed.json"))
    return {
        "cache_seed_json_count": len(cache_files),
        "cache_seed_glob": "shards/shard-*-of-*/cache/*.seed.json",
        "used_for_export": False,
        "reason_not_used": (
            "These files store rounded reported warm-start seeds; direct "
            "evaluation does not reproduce the stored raw-fit score. The "
            "exporter re-runs the deterministic fit from default_spurion_seed()."
        ),
    }


def _write_csv(path: Path, rows: list[dict[str, str]], results: list[WorkerResult]) -> None:
    result_by_id = {result.point_id: result for result in results}
    if len(result_by_id) != len(results):
        raise ValueError("duplicate point_id in worker results")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            point_id = row["point_id"]
            if point_id not in result_by_id:
                raise ValueError(f"missing fitted values for {point_id}")
            out = dict(row)
            out.update(result_by_id[point_id].values)
            writer.writerow(out)


def _verification_summary(
    rows: list[dict[str, str]],
    results: list[WorkerResult],
    output_csv: Path,
) -> dict[str, Any]:
    result_by_id = {result.point_id: result for result in results}
    byte_match_ok = True
    checked = 0
    with output_csv.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != CSV_COLUMNS:
            raise ValueError(f"{output_csv} schema mismatch after write")
        for original, exported in zip(rows, reader, strict=True):
            checked += 1
            for column in BASE_CSV_COLUMNS:
                if exported[column] != original[column]:
                    byte_match_ok = False
                    break
            if not byte_match_ok:
                break
    return {
        "base_column_string_match": byte_match_ok,
        "base_column_rows_checked": checked,
        "max_fit_score": max(result.fit_score for result in results),
        "max_fit_score_abs_delta_vs_jsonl": max(
            result.fit_score_abs_delta for result in results
        ),
        "max_abs_log_mass_residual": max(
            result.max_abs_log_mass_residual for result in results
        ),
        "max_abs_ckm_observable_residual": max(
            result.max_abs_ckm_observable_residual for result in results
        ),
        "c_min": min(result.c_min for result in results),
        "c_max": max(result.c_max for result in results),
        "fit_success_count": sum(1 for result in results if result.fit_success),
        "max_nfev_observed": max(result.nfev for result in results),
        "sample_point_ids": [
            rows[idx]["point_id"]
            for idx in sorted(
                set(
                    [
                        0,
                        max(0, len(rows) // 4),
                        max(0, len(rows) // 2),
                        max(0, (3 * len(rows)) // 4),
                        len(rows) - 1,
                    ]
                )
            )
        ],
        "sample_fit_scores": {
            point_id: result_by_id[point_id].fit_score
            for point_id in [
                rows[idx]["point_id"]
                for idx in sorted(
                    set(
                        [
                            0,
                            max(0, len(rows) // 4),
                            max(0, len(rows) // 2),
                            max(0, (3 * len(rows)) // 4),
                            len(rows) - 1,
                        ]
                    )
                )
            ]
        },
    }


def _write_provenance(
    path: Path,
    *,
    base_provenance_path: Path,
    accepted_csv_path: Path,
    jsonl_path: Path,
    manifest_path: Path,
    config_path: Path,
    run_dir: Path,
    row_count: int,
    total_points: int,
    points_passing_fit_filter: int,
    points_rejected_by_fit_score: int,
    workers: int,
    wall_seconds: float,
    max_nfev: int,
    fit_orientation: bool,
    verification: dict[str, Any],
    limited_to: int | None,
    cache: dict[str, Any],
) -> None:
    base_provenance = (
        json.loads(base_provenance_path.read_text())
        if base_provenance_path.exists()
        else {}
    )
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    provenance = dict(base_provenance)
    provenance.update(
        {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "source_run_dir": str(run_dir.resolve()),
            "source_results_jsonl": str(jsonl_path.resolve()),
            "source_results_sha256": _file_sha256(jsonl_path),
            "source_manifest_config_hash": manifest.get("config_hash"),
            "source_config_json": str(config_path.resolve()) if config_path.exists() else None,
            "source_accepted_points_csv": str(accepted_csv_path.resolve()),
            "source_accepted_points_sha256": _file_sha256(accepted_csv_path),
            "total_points_in_jsonl": total_points,
            "points_passing_fit_filter": points_passing_fit_filter,
            "points_rejected_by_fit_score": points_rejected_by_fit_score,
            "accepted_points_publication_convention": row_count,
            "accepted_points_with_yukawas_publication_convention": row_count,
            "xi_kk": GAUGE_KK_ROOT_NN,
            "xi_kk_source": "quarkConstraints.scales.GAUGE_KK_ROOT_NN",
            "csv_columns": CSV_COLUMNS,
            "base_csv_columns": BASE_CSV_COLUMNS,
            "yukawa_columns": YUKAWA_COLUMNS,
            "c_columns": C_COLUMNS,
            "total_csv_columns": len(CSV_COLUMNS),
            "systems": SYSTEMS,
            "generator_script": "scripts/export_accepted_quark_scan_with_yukawas.py",
            "limited_to_first_n_rows": limited_to,
            "reconstruction": {
                "source": "deterministic re-fit",
                "fit_function": "quarkConstraints.fit.fit_quark_sector",
                "seed_function": "quarkConstraints.benchmarks.default_spurion_seed",
                "target_function": "quarkConstraints.benchmarks.default_quark_targets",
                "max_nfev": max_nfev,
                "fit_orientation": fit_orientation,
                "worker_count": workers,
                "wall_seconds": wall_seconds,
                "cache": cache,
            },
            "yukawa_schema": {
                "basis": "bulk-mass eigenbasis",
                "source_attributes": [
                    "QuarkFitResult.bulk_state.Y_u_bulk_basis",
                    "QuarkFitResult.bulk_state.Y_d_bulk_basis",
                ],
                "column_order": "Y_u then Y_d; each matrix row-major with 1-based indices; real part then imaginary part",
                "columns": YUKAWA_COLUMNS,
                "complex_encoding": "<matrix>_<row><col>_re and <matrix>_<row><col>_im",
                "includes_overall_scale": True,
                "overall_scale_note": (
                    "The fit absorbs the scan-grid overall_scale into the fitted "
                    "singular values before optimization; the exported matrices "
                    "are the actual dimensionless Yukawas used to build the mass "
                    "matrices, not unit-norm spurions."
                ),
                "post_fit_rotation": (
                    "No mass-basis post-fit rotation is applied to these matrices; "
                    "they are the bulk-basis Yukawas used before SVD extraction of "
                    "the physical quark masses and CKM."
                ),
            },
            "c_schema": {
                "basis": "bulk-mass eigenbasis",
                "source_attributes": [
                    "QuarkFitResult.bulk_state.c_Q",
                    "QuarkFitResult.bulk_state.c_u",
                    "QuarkFitResult.bulk_state.c_d",
                ],
                "columns": C_COLUMNS,
                "order": (
                    "For each species, entries 1..3 follow the ordered Hermitian "
                    "spectrum used by quarkConstraints.model._ordered_hermitian_spectrum."
                ),
                "species": {
                    "c_Q": "left-handed quark doublet bulk-mass diagonal entries",
                    "c_u": "right-handed up-type singlet bulk-mass diagonal entries",
                    "c_d": "right-handed down-type singlet bulk-mass diagonal entries",
                },
                "map": "BulkMassMap: c = 0.72 - 0.42 * lambda / (lambda + 1)",
            },
            "verification": verification,
        }
    )
    path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")


def _default_workers() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except Exception:
        return os.cpu_count() or 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run",
        required=True,
        type=Path,
        help="Scan run directory (containing derived/accepted_points.csv and merged/results.jsonl)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: <run>/derived)",
    )
    parser.add_argument(
        "--accepted-csv",
        type=Path,
        default=None,
        help="Base accepted-points CSV (default: <out-dir>/accepted_points.csv)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=_default_workers(),
        help="Parallel worker processes (default: CPU affinity count)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process the first N accepted rows; writes a sample-named CSV by default",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=None,
        help="Output CSV path (default: accepted_points_with_yukawas.csv, or .sample<N>.csv with --limit)",
    )
    parser.add_argument(
        "--provenance-json",
        type=Path,
        default=None,
        help="Output provenance path (default: sibling .provenance.json)",
    )
    args = parser.parse_args()

    if args.workers <= 0:
        raise ValueError("--workers must be positive")
    if args.limit is not None and args.limit <= 0:
        raise ValueError("--limit must be positive when supplied")

    run_dir = args.run
    out_dir = args.out_dir or (run_dir / "derived")
    accepted_csv_path = args.accepted_csv or (out_dir / "accepted_points.csv")
    jsonl_path = run_dir / "merged" / "results.jsonl"
    manifest_path = run_dir / "merged" / "manifest.json"
    config_path = run_dir / "config.json"
    base_provenance_path = out_dir / "accepted_points.provenance.json"
    suffix = f".sample{args.limit}" if args.limit is not None and args.output_csv is None else ""
    output_csv = args.output_csv or (out_dir / f"accepted_points_with_yukawas{suffix}.csv")
    provenance_json = args.provenance_json or output_csv.with_suffix(".provenance.json")

    rows = _read_accepted_rows(accepted_csv_path)
    if args.limit is not None:
        rows = rows[: args.limit]
    accepted_point_ids = {row["point_id"] for row in rows}
    metadata, total_points, points_passing_fit_filter, points_rejected_by_fit_score = (
        _load_jsonl_metadata(jsonl_path, accepted_point_ids)
    )
    max_nfev, fit_orientation = _read_scan_options(config_path)
    work_items = _build_work_items(
        rows,
        metadata,
        max_nfev=max_nfev,
        fit_orientation=fit_orientation,
    )

    start = time.perf_counter()
    if args.workers == 1:
        results = [_solve_one(item) for item in work_items]
    else:
        chunksize = max(1, min(64, len(work_items) // (args.workers * 4) or 1))
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            results = list(executor.map(_solve_one, work_items, chunksize=chunksize))
    wall_seconds = time.perf_counter() - start

    _write_csv(output_csv, rows, results)
    verification = _verification_summary(rows, results, output_csv)
    _write_provenance(
        provenance_json,
        base_provenance_path=base_provenance_path,
        accepted_csv_path=accepted_csv_path,
        jsonl_path=jsonl_path,
        manifest_path=manifest_path,
        config_path=config_path,
        run_dir=run_dir,
        row_count=len(rows),
        total_points=total_points,
        points_passing_fit_filter=points_passing_fit_filter,
        points_rejected_by_fit_score=points_rejected_by_fit_score,
        workers=args.workers,
        wall_seconds=wall_seconds,
        max_nfev=max_nfev,
        fit_orientation=fit_orientation,
        verification=verification,
        limited_to=args.limit,
        cache=_cache_summary(run_dir),
    )

    print(
        f"Wrote {output_csv} ({len(rows)} rows, {len(CSV_COLUMNS)} columns; "
        f"{wall_seconds:.2f}s with {args.workers} worker(s))"
    )


if __name__ == "__main__":
    main()
