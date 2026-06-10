#!/usr/bin/env python3
"""Build paired WQ quark-only minimal-vs-custodial comparison artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence


SCHEMA_ID = "wq_quarkonly_minimal_vs_custodial_manifest_v1"
MINIMAL_RS_EW_MODEL = "minimal_rs"
CUSTODIAL_RS_PLR_EW_MODEL = "custodial_rs_plr"
PAIRING_KEY = ["r", "mkk_tev", "draw_seed"]
PHYSICS_TAGS = {
    "custodial_tree_status": "included",
    "custodial_top_partner_loop_status": "deferred",
    "custodial_fcnc_mode": "pr1_minimal_offdiag",
    "oblique_proxy_status": "tree_level_custodial_proxy_top_partner_loop_deferred",
}
RUN_IDENTITY_FIELDS = {
    "array_job_id",
    "created_utc",
    "dirty",
    "ew_model",
    "git_sha",
    "job_id",
    "output_dir",
    "output_root",
    "path",
    "run_id",
    "slurm_job_id",
}


class ComparisonValidationError(RuntimeError):
    """Raised when two run roots cannot be paired one-to-one."""


@dataclass(frozen=True)
class RunRows:
    root: Path
    ew_model: str
    rows: dict[tuple[float, float, int], dict[str, Any]]
    config_hashes: tuple[str, ...]
    git_shas: tuple[str, ...]
    scan_plan: dict[str, Any]


def build_comparison(minimal_root: Path, custodial_root: Path, output_dir: Path) -> dict[str, Any]:
    minimal_root = minimal_root.resolve()
    custodial_root = custodial_root.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    minimal = _load_run(minimal_root, ew_model=MINIMAL_RS_EW_MODEL)
    custodial = _load_run(custodial_root, ew_model=CUSTODIAL_RS_PLR_EW_MODEL)
    scan_plan_equivalent = _scan_plans_equivalent(minimal.scan_plan, custodial.scan_plan)
    if not scan_plan_equivalent:
        raise ComparisonValidationError("scan_plan.json differs beyond run identity fields")

    validation = _validate_pairing(minimal.rows, custodial.rows)
    paired_rows, veto_rows, survival_rows, constraint_rows = _build_artifact_rows(
        minimal,
        custodial,
    )

    artifacts = {
        "readme": "README.md",
        "manifest": "manifest.json",
        "schema": "schema.json",
        "paired_draws": "paired_draws.parquet",
        "paired_vetoes": "paired_vetoes.parquet",
        "survival_by_r_mkk": "survival_by_r_mkk.csv",
        "constraint_veto_by_r_mkk": "constraint_veto_by_r_mkk.csv",
        "run_index": "run_index.json",
    }

    _write_paired_draws(output_dir / artifacts["paired_draws"], paired_rows)
    _write_paired_vetoes(output_dir / artifacts["paired_vetoes"], veto_rows)
    _write_csv(output_dir / artifacts["survival_by_r_mkk"], SURVIVAL_COLUMNS, survival_rows)
    _write_csv(
        output_dir / artifacts["constraint_veto_by_r_mkk"],
        CONSTRAINT_VETO_COLUMNS,
        constraint_rows,
    )
    _write_json(output_dir / artifacts["schema"], _schema_payload())
    _write_text(
        output_dir / artifacts["readme"],
        _readme_text(minimal.root, custodial.root),
    )

    manifest = _manifest_payload(
        minimal,
        custodial,
        artifacts=artifacts,
        validation={
            **validation,
            "scan_plan_equivalent": scan_plan_equivalent,
            "passed": True,
        },
    )
    run_index = _run_index_payload(minimal.root, custodial.root, output_dir, artifacts)
    _write_json(output_dir / artifacts["run_index"], run_index)
    _write_json(output_dir / artifacts["manifest"], manifest)
    return manifest


def _load_run(root: Path, *, ew_model: str) -> RunRows:
    scan_plan = _load_scan_plan(root)
    rows: dict[tuple[float, float, int], dict[str, Any]] = {}
    duplicate_count = 0
    for path in _jsonl_paths(root):
        with path.open("r", encoding="utf-8") as fh:
            for line_no, line in enumerate(fh, start=1):
                if not line.strip():
                    continue
                try:
                    raw = json.loads(line)
                except json.JSONDecodeError as exc:
                    raise ComparisonValidationError(
                        f"{path}:{line_no}: invalid JSONL row: {exc}"
                    ) from exc
                row = _normalize_row(raw, path=path, line_no=line_no)
                key = (row["r"], row["mkk_tev"], row["draw_seed"])
                if key in rows:
                    duplicate_count += 1
                    continue
                rows[key] = row
    if duplicate_count:
        raise ComparisonValidationError(f"{root}: duplicate pairing keys: {duplicate_count}")
    if not rows:
        raise ComparisonValidationError(f"{root}: no tile-*.jsonl rows found")

    config_hashes = {
        str(row["config_hash"])
        for row in rows.values()
        if row.get("config_hash")
    }
    git_shas = {
        str(row["git_sha"])
        for row in rows.values()
        if row.get("git_sha")
    }
    for summary in _summary_payloads(root):
        if summary.get("config_hash"):
            config_hashes.add(str(summary["config_hash"]))
        provenance = dict(summary.get("provenance") or {})
        if provenance.get("git_sha"):
            git_shas.add(str(provenance["git_sha"]))
    return RunRows(
        root=root,
        ew_model=ew_model,
        rows=rows,
        config_hashes=tuple(sorted(config_hashes)),
        git_shas=tuple(sorted(git_shas)),
        scan_plan=scan_plan,
    )


def _normalize_row(raw: Mapping[str, Any], *, path: Path, line_no: int) -> dict[str, Any]:
    try:
        params = dict(raw["params"])
        draw_seed = int(raw["seed"])
        mkk_gev = float(params["M_KK"])
        r_source = raw.get("quark_fit_r", params.get("quark_fit_r"))
        r_value = float(r_source)
        tile_id = int(raw["tile_id"])
        draw_id = int(raw["draw_id"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ComparisonValidationError(
            f"{path}:{line_no}: row lacks required WQ pairing fields"
        ) from exc

    constraints = dict(raw.get("constraints") or {})
    provenance = dict(raw.get("provenance") or {})
    return {
        "r": r_value,
        "mkk_tev": mkk_gev / 1000.0,
        "mkk_gev": mkk_gev,
        "draw_seed": draw_seed,
        "tile_id": tile_id,
        "draw_id": draw_id,
        "skipped": bool(raw.get("skipped", False)),
        "skip_reason": _nullable_str(raw.get("skip_reason")),
        "survives_strict": bool(raw.get("survives_all_HARD_strict", False)),
        "survives_inclusive": bool(raw.get("survives_all_HARD_inclusive", False)),
        "excluded_by_rigorous": _string_list(raw.get("excluded_by_rigorous", [])),
        "excluded_by_proxy": _string_list(raw.get("excluded_by_proxy", [])),
        "hard_not_evaluated": _string_list(raw.get("hard_not_evaluated", [])),
        "constraints": constraints,
        "config_hash": _nullable_str(provenance.get("config_hash")),
        "git_sha": _nullable_str(provenance.get("git_sha")),
    }


def _validate_pairing(
    minimal_rows: Mapping[tuple[float, float, int], Mapping[str, Any]],
    custodial_rows: Mapping[tuple[float, float, int], Mapping[str, Any]],
) -> dict[str, Any]:
    minimal_keys = set(minimal_rows)
    custodial_keys = set(custodial_rows)
    minimal_only = minimal_keys - custodial_keys
    custodial_only = custodial_keys - minimal_keys
    minimal_seed_groups = {
        int(row["draw_seed"]): (float(row["r"]), float(row["mkk_tev"]))
        for row in minimal_rows.values()
    }
    custodial_seed_groups = {
        int(row["draw_seed"]): (float(row["r"]), float(row["mkk_tev"]))
        for row in custodial_rows.values()
    }
    seed_group_mismatches = sum(
        1
        for seed in set(minimal_seed_groups) & set(custodial_seed_groups)
        if minimal_seed_groups[seed] != custodial_seed_groups[seed]
    )
    validation = {
        "minimal_rows": len(minimal_rows),
        "custodial_rows": len(custodial_rows),
        "paired_rows": len(minimal_keys & custodial_keys),
        "duplicate_key_counts": {
            "minimal": 0,
            "custodial": 0,
        },
        "missing_pair_counts": {
            "minimal_only": len(minimal_only),
            "custodial_only": len(custodial_only),
            "draw_seed_group_mismatches": seed_group_mismatches,
        },
    }
    if minimal_only or custodial_only or seed_group_mismatches:
        raise ComparisonValidationError(
            "run roots do not pair one-to-one on (r, mkk_tev, draw_seed): "
            f"{validation['missing_pair_counts']}"
        )
    return validation


def _build_artifact_rows(
    minimal: RunRows,
    custodial: RunRows,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    paired_rows: list[dict[str, Any]] = []
    veto_rows: list[dict[str, Any]] = []
    survival = defaultdict(_survival_counter)
    constraint_veto = defaultdict(_constraint_counter)

    for key in sorted(minimal.rows):
        min_row = minimal.rows[key]
        cust_row = custodial.rows[key]
        paired_rows.append(
            {
                "r": min_row["r"],
                "mkk_tev": min_row["mkk_tev"],
                "mkk_gev": min_row["mkk_gev"],
                "draw_seed": min_row["draw_seed"],
                "tile_id_minimal": min_row["tile_id"],
                "tile_id_custodial": cust_row["tile_id"],
                "draw_id_minimal": min_row["draw_id"],
                "draw_id_custodial": cust_row["draw_id"],
                "minimal_skipped": min_row["skipped"],
                "custodial_skipped": cust_row["skipped"],
                "minimal_skip_reason": min_row["skip_reason"],
                "custodial_skip_reason": cust_row["skip_reason"],
                "minimal_survives_strict": min_row["survives_strict"],
                "minimal_survives_inclusive": min_row["survives_inclusive"],
                "custodial_survives_strict": cust_row["survives_strict"],
                "custodial_survives_inclusive": cust_row["survives_inclusive"],
                "minimal_excluded_by_rigorous": min_row["excluded_by_rigorous"],
                "minimal_excluded_by_proxy": min_row["excluded_by_proxy"],
                "custodial_excluded_by_rigorous": cust_row["excluded_by_rigorous"],
                "custodial_excluded_by_proxy": cust_row["excluded_by_proxy"],
                "minimal_hard_not_evaluated": min_row["hard_not_evaluated"],
                "custodial_hard_not_evaluated": cust_row["hard_not_evaluated"],
                "minimal_config_hash": min_row["config_hash"],
                "custodial_config_hash": cust_row["config_hash"],
                "minimal_git_sha": min_row["git_sha"],
                "custodial_git_sha": cust_row["git_sha"],
            }
        )
        _accumulate_survival(survival[(min_row["r"], min_row["mkk_tev"])], "minimal", min_row)
        _accumulate_survival(survival[(cust_row["r"], cust_row["mkk_tev"])], "custodial", cust_row)
        for model, row in (
            (MINIMAL_RS_EW_MODEL, min_row),
            (CUSTODIAL_RS_PLR_EW_MODEL, cust_row),
        ):
            _accumulate_constraint_vetoes(constraint_veto, model, row)
            veto_rows.extend(_draw_veto_rows(model, row))

    survival_rows = [_finalize_survival_row(key, values) for key, values in sorted(survival.items())]
    constraint_rows = [
        _finalize_constraint_row(key, values)
        for key, values in sorted(constraint_veto.items())
    ]
    return paired_rows, veto_rows, survival_rows, constraint_rows


def _draw_veto_rows(ew_model: str, row: Mapping[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    constraints = dict(row.get("constraints") or {})
    for pid, raw_item in sorted(constraints.items()):
        item = dict(raw_item)
        tag = _nullable_str(item.get("tag"))
        severity = _nullable_str(item.get("severity"))
        if (
            severity == "HARD"
            and bool(item.get("evaluated"))
            and item.get("passes") is False
            and tag in {"rigorous", "proxy"}
        ):
            out.append(_veto_row(ew_model, row, pid, item, veto_class=str(tag)))
    for pid in row.get("hard_not_evaluated", []):
        item = dict(constraints.get(pid) or {})
        out.append(
            {
                "r": row["r"],
                "mkk_tev": row["mkk_tev"],
                "draw_seed": row["draw_seed"],
                "ew_model": ew_model,
                "constraint_id": str(pid),
                "veto_class": "not_evaluated",
                "tag": None,
                "severity": _nullable_str(item.get("severity", "HARD")),
                "ratio": None,
                "passes": None,
                "evaluated": False,
                "active": _nullable_bool(item.get("active")),
            }
        )
    return out


def _veto_row(
    ew_model: str,
    row: Mapping[str, Any],
    pid: str,
    item: Mapping[str, Any],
    *,
    veto_class: str,
) -> dict[str, Any]:
    return {
        "r": row["r"],
        "mkk_tev": row["mkk_tev"],
        "draw_seed": row["draw_seed"],
        "ew_model": ew_model,
        "constraint_id": str(pid),
        "veto_class": veto_class,
        "tag": _nullable_str(item.get("tag")),
        "severity": _nullable_str(item.get("severity")),
        "ratio": _nullable_float(item.get("ratio")),
        "passes": _nullable_bool(item.get("passes")),
        "evaluated": bool(item.get("evaluated")),
        "active": _nullable_bool(item.get("active")),
    }


def _accumulate_survival(counter: dict[str, int], prefix: str, row: Mapping[str, Any]) -> None:
    counter[f"{prefix}_rows"] += 1
    if row["skipped"]:
        counter[f"{prefix}_skipped"] += 1
        return
    counter[f"{prefix}_evaluated"] += 1
    if row["survives_strict"]:
        counter[f"{prefix}_survives_strict"] += 1
    if row["survives_inclusive"]:
        counter[f"{prefix}_survives_inclusive"] += 1


def _accumulate_constraint_vetoes(
    aggregates: dict[tuple[str, float, float, str, str | None, str | None], dict[str, int]],
    ew_model: str,
    row: Mapping[str, Any],
) -> None:
    for pid, raw_item in dict(row.get("constraints") or {}).items():
        item = dict(raw_item)
        tag = _nullable_str(item.get("tag"))
        severity = _nullable_str(item.get("severity"))
        key = (ew_model, row["r"], row["mkk_tev"], str(pid), tag, severity)
        counter = aggregates[key]
        counter["points"] += 1
        if item.get("evaluated"):
            counter["evaluated"] += 1
        if item.get("active"):
            counter["active"] += 1
        failed = item.get("passes") is False
        if failed:
            counter["failed"] += 1
        if failed and severity == "HARD" and item.get("evaluated") and tag in {"rigorous", "proxy"}:
            counter["vetoed"] += 1


def _finalize_survival_row(
    key: tuple[float, float],
    values: Mapping[str, int],
) -> dict[str, Any]:
    r_value, mkk_tev = key
    row: dict[str, Any] = {"r": r_value, "mkk_tev": mkk_tev}
    for prefix in ("minimal", "custodial"):
        rows = int(values.get(f"{prefix}_rows", 0))
        evaluated = int(values.get(f"{prefix}_evaluated", 0))
        skipped = int(values.get(f"{prefix}_skipped", 0))
        strict = int(values.get(f"{prefix}_survives_strict", 0))
        inclusive = int(values.get(f"{prefix}_survives_inclusive", 0))
        row.update(
            {
                f"{prefix}_rows": rows,
                f"{prefix}_evaluated": evaluated,
                f"{prefix}_skipped": skipped,
                f"{prefix}_survives_strict": strict,
                f"{prefix}_survives_inclusive": inclusive,
                f"{prefix}_strict_survival_frac": _safe_div(strict, evaluated),
                f"{prefix}_inclusive_survival_frac": _safe_div(inclusive, evaluated),
            }
        )
    row["delta_strict_survival_frac"] = _nullable_delta(
        row["custodial_strict_survival_frac"],
        row["minimal_strict_survival_frac"],
    )
    row["delta_inclusive_survival_frac"] = _nullable_delta(
        row["custodial_inclusive_survival_frac"],
        row["minimal_inclusive_survival_frac"],
    )
    return row


def _finalize_constraint_row(
    key: tuple[str, float, float, str, str | None, str | None],
    values: Mapping[str, int],
) -> dict[str, Any]:
    ew_model, r_value, mkk_tev, constraint_id, tag, severity = key
    points = int(values.get("points", 0))
    vetoed = int(values.get("vetoed", 0))
    return {
        "ew_model": ew_model,
        "r": r_value,
        "mkk_tev": mkk_tev,
        "constraint_id": constraint_id,
        "tag": tag,
        "severity": severity,
        "points": points,
        "evaluated": int(values.get("evaluated", 0)),
        "active": int(values.get("active", 0)),
        "failed": int(values.get("failed", 0)),
        "vetoed": vetoed,
        "veto_fraction": _safe_div(vetoed, points),
    }


def _manifest_payload(
    minimal: RunRows,
    custodial: RunRows,
    *,
    artifacts: Mapping[str, str],
    validation: Mapping[str, Any],
) -> dict[str, Any]:
    plan = minimal.scan_plan
    return {
        "schema": SCHEMA_ID,
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "minimal_run": {
            "path": str(minimal.root),
            "run_id": minimal.root.name,
            "git_sha": _single_or_list(minimal.git_shas),
            "config_hashes": list(minimal.config_hashes),
            "ew_model": MINIMAL_RS_EW_MODEL,
        },
        "custodial_run": {
            "path": str(custodial.root),
            "run_id": custodial.root.name,
            "git_sha": _single_or_list(custodial.git_shas),
            "config_hashes": list(custodial.config_hashes),
            "ew_model": CUSTODIAL_RS_PLR_EW_MODEL,
        },
        "grid": {
            "r_grid": list(plan.get("r_grid", _sorted_unique(row["r"] for row in minimal.rows.values()))),
            "mkk_tev": list(
                plan.get("mkk_tev", _sorted_unique(row["mkk_tev"] for row in minimal.rows.values()))
            ),
        },
        "seed_contract": {
            "base_seed": plan.get("base_seed"),
            "tile_seed_stride": plan.get("tile_seed_stride"),
            "shard_seed_block": plan.get("shard_seed_block"),
            "draw_seed_formula": dict(plan.get("seed_disjointness") or {}).get("draw_seed_formula"),
            "tile_seed_formula": dict(plan.get("seed_disjointness") or {}).get("tile_seed_formula"),
            "per_shard_base_formula": dict(plan.get("seed_disjointness") or {}).get(
                "per_shard_base_formula"
            ),
        },
        "pairing_key": list(PAIRING_KEY),
        "artifacts": dict(artifacts),
        "validation": dict(validation),
        "physics_tags": dict(PHYSICS_TAGS),
    }


def _run_index_payload(
    minimal_root: Path,
    custodial_root: Path,
    output_dir: Path,
    artifacts: Mapping[str, str],
) -> dict[str, Any]:
    return {
        "schema": "wq_quarkonly_comparison_run_index_v1",
        "minimal": {
            "raw_root": str(minimal_root),
            "cache": _existing_or_default(minimal_root / "wq_quarkonly_cache.parquet"),
            "analysis": _existing_or_default(minimal_root / "analysis"),
        },
        "custodial": {
            "raw_root": str(custodial_root),
            "cache": _existing_or_default(custodial_root / "wq_quarkonly_cache.parquet"),
            "analysis": _existing_or_default(custodial_root / "analysis"),
        },
        "comparison": {
            key: str(output_dir / rel_path)
            for key, rel_path in artifacts.items()
        },
    }


def _schema_payload() -> dict[str, Any]:
    return {
        "schema": "wq_quarkonly_minimal_vs_custodial_schema_v1",
        "artifacts": {
            "paired_draws.parquet": {
                "description": "One row per normalized paired draw key.",
                "columns": _column_schema(PAIRED_DRAWS_COLUMNS),
            },
            "paired_vetoes.parquet": {
                "description": "Long format constraint vetoes and explicit HARD coverage gaps.",
                "columns": _column_schema(PAIRED_VETOES_COLUMNS),
            },
            "survival_by_r_mkk.csv": {
                "description": "Survival aggregates per (r, M_KK) for both models.",
                "columns": _column_schema(SURVIVAL_COLUMN_META),
            },
            "constraint_veto_by_r_mkk.csv": {
                "description": "Constraint aggregate counts per model, (r, M_KK), id, tag, severity.",
                "columns": _column_schema(CONSTRAINT_VETO_COLUMN_META),
            },
        },
    }


def _column_schema(columns: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    return {str(item["name"]): {k: v for k, v in item.items() if k != "name"} for item in columns}


def _readme_text(minimal_root: Path, custodial_root: Path) -> str:
    return "\n".join(
        [
            "# WQ Quark-Only Minimal vs Custodial Comparison",
            "",
            f"This compares baseline non-custodial `{minimal_root}` against custodial `{custodial_root}`.",
            "",
            "Rows are paired by normalized `(r, mkk_tev, draw_seed)`.",
            "Raw scan rows use `seed`, not `draw_seed`; this builder normalizes `row[\"seed\"]`.",
            "Raw scan rows use `params.M_KK` in GeV; this builder writes `mkk_tev` in TeV.",
            "Raw quark-only rows use top-level `quark_fit_r`; `params.quark_fit_r` is only a fallback.",
            "",
            "Survival uses the existing row semantics: `survives_all_HARD_strict` and "
            "`survives_all_HARD_inclusive` from each raw row.",
            "Skipped fit rows remain in `paired_draws.parquet` but are not counted as evaluated "
            "survival in `survival_by_r_mkk.csv`.",
            "",
            "W8 top-partner loop status for this W9 scan: deferred.",
            "",
        ]
    )


def _jsonl_paths(root: Path) -> list[Path]:
    return sorted(
        path
        for path in root.rglob("tile-*.jsonl")
        if "comparison" not in path.relative_to(root).parts and ".tmp." not in path.name
    )


def _summary_payloads(root: Path) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for path in sorted(root.rglob("tile-*.summary.json")):
        if "comparison" in path.relative_to(root).parts:
            continue
        try:
            out.append(json.loads(path.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError):
            continue
    return out


def _load_scan_plan(root: Path) -> dict[str, Any]:
    path = root / "scan_plan.json"
    if not path.is_file():
        raise ComparisonValidationError(f"{root}: missing scan_plan.json")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ComparisonValidationError(f"{path}: invalid scan_plan.json: {exc}") from exc


def _scan_plans_equivalent(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    return _scrub_run_identity(left) == _scrub_run_identity(right)


def _scrub_run_identity(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {
            str(k): _scrub_run_identity(v)
            for k, v in sorted(value.items())
            if str(k) not in RUN_IDENTITY_FIELDS
        }
    if isinstance(value, list):
        return [_scrub_run_identity(item) for item in value]
    return value


def _write_paired_draws(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    schema = pa.schema(
        [
            ("r", pa.float64()),
            ("mkk_tev", pa.float64()),
            ("mkk_gev", pa.float64()),
            ("draw_seed", pa.int64()),
            ("tile_id_minimal", pa.int32()),
            ("tile_id_custodial", pa.int32()),
            ("draw_id_minimal", pa.int32()),
            ("draw_id_custodial", pa.int32()),
            ("minimal_skipped", pa.bool_()),
            ("custodial_skipped", pa.bool_()),
            ("minimal_skip_reason", pa.string()),
            ("custodial_skip_reason", pa.string()),
            ("minimal_survives_strict", pa.bool_()),
            ("minimal_survives_inclusive", pa.bool_()),
            ("custodial_survives_strict", pa.bool_()),
            ("custodial_survives_inclusive", pa.bool_()),
            ("minimal_excluded_by_rigorous", pa.list_(pa.string())),
            ("minimal_excluded_by_proxy", pa.list_(pa.string())),
            ("custodial_excluded_by_rigorous", pa.list_(pa.string())),
            ("custodial_excluded_by_proxy", pa.list_(pa.string())),
            ("minimal_hard_not_evaluated", pa.list_(pa.string())),
            ("custodial_hard_not_evaluated", pa.list_(pa.string())),
            ("minimal_config_hash", pa.string()),
            ("custodial_config_hash", pa.string()),
            ("minimal_git_sha", pa.string()),
            ("custodial_git_sha", pa.string()),
        ]
    )
    table = pa.Table.from_pylist([dict(row) for row in rows], schema=schema)
    _atomic_write(path, lambda tmp: pq.write_table(table, tmp, compression="snappy"))


def _write_paired_vetoes(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    schema = pa.schema(
        [
            ("r", pa.float64()),
            ("mkk_tev", pa.float64()),
            ("draw_seed", pa.int64()),
            ("ew_model", pa.string()),
            ("constraint_id", pa.string()),
            ("veto_class", pa.string()),
            ("tag", pa.string()),
            ("severity", pa.string()),
            ("ratio", pa.float64()),
            ("passes", pa.bool_()),
            ("evaluated", pa.bool_()),
            ("active", pa.bool_()),
        ]
    )
    table = pa.Table.from_pylist([dict(row) for row in rows], schema=schema)
    _atomic_write(path, lambda tmp: pq.write_table(table, tmp, compression="snappy"))


def _write_csv(path: Path, columns: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    def writer(tmp: Path) -> None:
        with tmp.open("w", encoding="utf-8", newline="") as fh:
            csv_writer = csv.DictWriter(fh, fieldnames=list(columns), extrasaction="ignore")
            csv_writer.writeheader()
            for row in rows:
                csv_writer.writerow({name: _csv_cell(row.get(name)) for name in columns})

    _atomic_write(path, writer)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    _write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _write_text(path: Path, text: str) -> None:
    def writer(tmp: Path) -> None:
        tmp.write_text(text, encoding="utf-8")

    _atomic_write(path, writer)


def _atomic_write(path: Path, writer: Callable[[Path], None]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    try:
        writer(tmp)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def _survival_counter() -> dict[str, int]:
    return defaultdict(int)


def _constraint_counter() -> dict[str, int]:
    return defaultdict(int)


def _string_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return [str(item) for item in value]


def _nullable_str(value: Any) -> str | None:
    if value is None:
        return None
    return str(value)


def _nullable_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _nullable_bool(value: Any) -> bool | None:
    if value is None:
        return None
    return bool(value)


def _safe_div(numerator: int | float, denominator: int | float) -> float | None:
    denominator = float(denominator)
    if denominator == 0.0:
        return None
    return float(numerator) / denominator


def _nullable_delta(left: float | None, right: float | None) -> float | None:
    if left is None or right is None:
        return None
    return float(left) - float(right)


def _csv_cell(value: Any) -> Any:
    if value is None:
        return ""
    return value


def _sorted_unique(values: Sequence[float] | Any) -> list[float]:
    return sorted({float(value) for value in values})


def _single_or_list(values: Sequence[str]) -> str | list[str] | None:
    if not values:
        return None
    if len(values) == 1:
        return values[0]
    return list(values)


def _existing_or_default(path: Path) -> str | None:
    return str(path) if path.exists() else None


def _build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build paired WQ quark-only minimal-vs-custodial comparison artifacts.",
    )
    parser.add_argument("--minimal", required=True, type=Path)
    parser.add_argument("--custodial", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_argparser().parse_args(argv)
    build_comparison(args.minimal, args.custodial, args.output)
    return 0


PAIRED_DRAWS_COLUMNS = [
    {"name": "r", "type": "float64", "unit": "dimensionless", "nullable": False, "description": "Quark-fit r grid value."},
    {"name": "mkk_tev", "type": "float64", "unit": "TeV", "nullable": False, "description": "KK scale in TeV."},
    {"name": "mkk_gev", "type": "float64", "unit": "GeV", "nullable": False, "description": "KK scale in GeV."},
    {"name": "draw_seed", "type": "int64", "unit": None, "nullable": False, "description": "Normalized raw row seed."},
    {"name": "tile_id_minimal", "type": "int32", "unit": None, "nullable": False, "description": "Minimal-run tile id."},
    {"name": "tile_id_custodial", "type": "int32", "unit": None, "nullable": False, "description": "Custodial-run tile id."},
    {"name": "draw_id_minimal", "type": "int32", "unit": None, "nullable": False, "description": "Minimal-run draw id."},
    {"name": "draw_id_custodial", "type": "int32", "unit": None, "nullable": False, "description": "Custodial-run draw id."},
    {"name": "minimal_skipped", "type": "bool", "unit": None, "nullable": False, "description": "Whether the minimal draw skipped before evaluation."},
    {"name": "custodial_skipped", "type": "bool", "unit": None, "nullable": False, "description": "Whether the custodial draw skipped before evaluation."},
    {"name": "minimal_skip_reason", "type": "string", "unit": None, "nullable": True, "description": "Minimal skip reason."},
    {"name": "custodial_skip_reason", "type": "string", "unit": None, "nullable": True, "description": "Custodial skip reason."},
    {"name": "minimal_survives_strict", "type": "bool", "unit": None, "nullable": False, "description": "Minimal strict HARD survival."},
    {"name": "minimal_survives_inclusive", "type": "bool", "unit": None, "nullable": False, "description": "Minimal inclusive HARD survival."},
    {"name": "custodial_survives_strict", "type": "bool", "unit": None, "nullable": False, "description": "Custodial strict HARD survival."},
    {"name": "custodial_survives_inclusive", "type": "bool", "unit": None, "nullable": False, "description": "Custodial inclusive HARD survival."},
    {"name": "minimal_excluded_by_rigorous", "type": "list<string>", "unit": None, "nullable": False, "description": "Minimal rigorous HARD veto ids."},
    {"name": "minimal_excluded_by_proxy", "type": "list<string>", "unit": None, "nullable": False, "description": "Minimal proxy HARD veto ids."},
    {"name": "custodial_excluded_by_rigorous", "type": "list<string>", "unit": None, "nullable": False, "description": "Custodial rigorous HARD veto ids."},
    {"name": "custodial_excluded_by_proxy", "type": "list<string>", "unit": None, "nullable": False, "description": "Custodial proxy HARD veto ids."},
    {"name": "minimal_hard_not_evaluated", "type": "list<string>", "unit": None, "nullable": False, "description": "Minimal HARD coverage gaps."},
    {"name": "custodial_hard_not_evaluated", "type": "list<string>", "unit": None, "nullable": False, "description": "Custodial HARD coverage gaps."},
    {"name": "minimal_config_hash", "type": "string", "unit": None, "nullable": True, "description": "Minimal config hash."},
    {"name": "custodial_config_hash", "type": "string", "unit": None, "nullable": True, "description": "Custodial config hash."},
    {"name": "minimal_git_sha", "type": "string", "unit": None, "nullable": True, "description": "Minimal git SHA."},
    {"name": "custodial_git_sha", "type": "string", "unit": None, "nullable": True, "description": "Custodial git SHA."},
]

PAIRED_VETOES_COLUMNS = [
    {"name": "r", "type": "float64", "unit": "dimensionless", "nullable": False, "description": "Quark-fit r grid value."},
    {"name": "mkk_tev", "type": "float64", "unit": "TeV", "nullable": False, "description": "KK scale in TeV."},
    {"name": "draw_seed", "type": "int64", "unit": None, "nullable": False, "description": "Normalized raw row seed."},
    {"name": "ew_model", "type": "enum<string>", "unit": None, "nullable": False, "description": "minimal_rs or custodial_rs_plr."},
    {"name": "constraint_id", "type": "string", "unit": None, "nullable": False, "description": "Constraint process id."},
    {"name": "veto_class", "type": "enum<string>", "unit": None, "nullable": False, "description": "rigorous, proxy, or not_evaluated."},
    {"name": "tag", "type": "string", "unit": None, "nullable": True, "description": "Harness classification tag."},
    {"name": "severity", "type": "string", "unit": None, "nullable": True, "description": "Constraint severity."},
    {"name": "ratio", "type": "float64", "unit": None, "nullable": True, "description": "Constraint ratio when available."},
    {"name": "passes", "type": "bool", "unit": None, "nullable": True, "description": "Constraint pass flag when evaluated."},
    {"name": "evaluated", "type": "bool", "unit": None, "nullable": False, "description": "Whether the constraint was evaluated."},
    {"name": "active", "type": "bool", "unit": None, "nullable": True, "description": "Whether the constraint was active."},
]

SURVIVAL_COLUMNS = [
    "r",
    "mkk_tev",
    "minimal_rows",
    "minimal_evaluated",
    "minimal_skipped",
    "minimal_survives_strict",
    "minimal_survives_inclusive",
    "minimal_strict_survival_frac",
    "minimal_inclusive_survival_frac",
    "custodial_rows",
    "custodial_evaluated",
    "custodial_skipped",
    "custodial_survives_strict",
    "custodial_survives_inclusive",
    "custodial_strict_survival_frac",
    "custodial_inclusive_survival_frac",
    "delta_strict_survival_frac",
    "delta_inclusive_survival_frac",
]

SURVIVAL_COLUMN_META = [
    {"name": name, "type": "float64" if name in {"r", "mkk_tev"} or name.endswith("_frac") else "int64", "unit": "TeV" if name == "mkk_tev" else None, "nullable": name.endswith("_frac"), "description": name.replace("_", " ")}
    for name in SURVIVAL_COLUMNS
]

CONSTRAINT_VETO_COLUMNS = [
    "ew_model",
    "r",
    "mkk_tev",
    "constraint_id",
    "tag",
    "severity",
    "points",
    "evaluated",
    "active",
    "failed",
    "vetoed",
    "veto_fraction",
]

CONSTRAINT_VETO_COLUMN_META = [
    {"name": "ew_model", "type": "enum<string>", "unit": None, "nullable": False, "description": "minimal_rs or custodial_rs_plr"},
    {"name": "r", "type": "float64", "unit": "dimensionless", "nullable": False, "description": "Quark-fit r grid value"},
    {"name": "mkk_tev", "type": "float64", "unit": "TeV", "nullable": False, "description": "KK scale in TeV"},
    {"name": "constraint_id", "type": "string", "unit": None, "nullable": False, "description": "Constraint process id"},
    {"name": "tag", "type": "string", "unit": None, "nullable": True, "description": "Harness tag"},
    {"name": "severity", "type": "string", "unit": None, "nullable": True, "description": "Constraint severity"},
    {"name": "points", "type": "int64", "unit": None, "nullable": False, "description": "Constraint records"},
    {"name": "evaluated", "type": "int64", "unit": None, "nullable": False, "description": "Evaluated records"},
    {"name": "active", "type": "int64", "unit": None, "nullable": False, "description": "Active records"},
    {"name": "failed", "type": "int64", "unit": None, "nullable": False, "description": "Failed records"},
    {"name": "vetoed", "type": "int64", "unit": None, "nullable": False, "description": "Evaluated HARD rigorous/proxy failures"},
    {"name": "veto_fraction", "type": "float64", "unit": None, "nullable": True, "description": "vetoed / points"},
]


if __name__ == "__main__":
    raise SystemExit(main())
