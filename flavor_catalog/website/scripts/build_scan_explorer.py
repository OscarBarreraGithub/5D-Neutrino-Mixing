#!/usr/bin/env python3
"""Build the compact Scan Explorer data artifact from quark-only scan JSONL.

The raw scans are large, so this script streams rows and keeps only aggregate
counts plus compact per-cell numeric arrays for exact percentile summaries.
"""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from array import array
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np

try:
    import orjson
except ImportError:  # pragma: no cover - exercised only where orjson is absent.
    orjson = None


REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from quarkConstraints.benchmarks import default_quark_targets
from quarkConstraints.fit import QuarkFitSeed, fit_quark_sector
from scripts.run_full_catalog_scan import (
    DEFAULT_QUARK_FIT_MAX_NFEV,
    _svd_seed_parts,
)


MINIMAL_ROOT = REPO_ROOT / "scan_outputs/wq_quarkonly_1M_20128400"
CUSTODIAL_ROOT = REPO_ROOT / "scan_outputs/wq_quarkonly_1M_custodial_20675555"
OUTPUT_PATH = REPO_ROOT / "flavor_catalog/website/src/content/scan_explorer.json"
SUMMARY_PATH = REPO_ROOT / ".orchestration/runs/WEB-EXPLORER/data_impl_summary.md"
ENTRIES_DIR = REPO_ROOT / "flavor_catalog/website/src/content/entries"

FLOOR_THRESHOLD = 0.5
R_GRID = [0.05, 0.1, 0.25, 0.5, 1.0]
MKK_GRID_TEV = [1, 2, 3, 5, 7, 10, 15, 20, 30, 50]
MODEL_ROOTS = {
    "minimal": MINIMAL_ROOT,
    "custodial": CUSTODIAL_ROOT,
}
RAW_EW_MODEL = {
    "minimal": "minimal_rs",
    "custodial": "custodial_rs_plr",
}

REP_TARGETS = [
    ("low r (down-localized)", "custodial", 0.05),
    ("r = 0.25", "custodial", 0.25),
    ("r = 1.0 (up-dominated)", "custodial", 1.0),
]
REP_MKK_PREFERENCE = [3, 2, 5, 7, 10, 1, 15, 20, 30, 50]


@dataclass
class SingularCell:
    up: list[array] = field(default_factory=lambda: [array("d") for _ in range(3)])
    down: list[array] = field(default_factory=lambda: [array("d") for _ in range(3)])

    def append(self, up_values: Iterable[float], down_values: Iterable[float]) -> None:
        for idx, value in enumerate(up_values):
            self.up[idx].append(float(value))
        for idx, value in enumerate(down_values):
            self.down[idx].append(float(value))

    @property
    def n(self) -> int:
        return len(self.up[0])


@dataclass
class RepCandidate:
    label: str
    ew_model: str
    r: float
    mkk_tev: int
    rank: int
    row: dict[str, Any]


def _json_loads(line: bytes) -> dict[str, Any]:
    if orjson is not None:
        return orjson.loads(line)
    return json.loads(line.decode("utf-8"))


def _r_key(value: float) -> str:
    return str(float(value))


def _mkk_key(value: int | float) -> str:
    value_f = float(value)
    if value_f.is_integer():
        return str(int(value_f))
    return str(value_f)


def _mkk_tev_from_row(row: dict[str, Any]) -> int:
    mkk_gev = float(row["params"]["M_KK"])
    mkk_tev = mkk_gev / 1000.0
    rounded = int(round(mkk_tev))
    if not math.isclose(mkk_tev, rounded, rel_tol=0.0, abs_tol=1.0e-9):
        raise ValueError(f"unexpected non-integer M_KK grid value: {mkk_gev}")
    return rounded


def _round_fixed(value: float, digits: int = 6) -> float:
    return round(float(value), digits)


def _round_sig(value: float, sig: int = 6) -> float:
    value = float(value)
    if value == 0.0 or not math.isfinite(value):
        return value
    return float(f"{value:.{sig}g}")


def _round_matrix_abs(matrix: np.ndarray, sig: int = 4) -> list[list[float]]:
    return [[_round_sig(abs(value), sig=sig) for value in row] for row in matrix]


def _iter_jsonl_paths(root: Path) -> list[Path]:
    return sorted(root.glob("r*/shard-*/tile-*.jsonl"))


def _git_sha() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _load_entry(cid: str) -> dict[str, Any]:
    path = ENTRIES_DIR / f"{cid}.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _constraint_label(cid: str, entry: dict[str, Any]) -> str:
    return (
        entry.get("title")
        or entry.get("standard_notation")
        or entry.get("process_name")
        or cid
    )


def _constraint_group(cid: str, entry: dict[str, Any]) -> str:
    family = str(entry.get("family") or "")
    if cid.startswith("CR") or family == "collider_rs":
        return "collider"
    if cid.startswith(("K", "B", "C")) or family in {"kaon", "beauty", "charm"}:
        return "meson_mixing"
    if cid.startswith(("EW", "T")) or family == "top_higgs_ew":
        return "electroweak"
    if cid.startswith("E") or family == "edm_neutrino":
        return "edm"
    return family or "other"


def _tag_for(cid: str, tag_counts: dict[str, Counter[str]]) -> str:
    counts = tag_counts.get(cid)
    if not counts:
        return "unknown"
    priority = {"rigorous": 0, "proxy": 1, "partial": 2, "unknown": 3}
    return sorted(counts.items(), key=lambda item: (-item[1], priority.get(item[0], 9), item[0]))[0][0]


def _empty_veto_tree() -> dict[str, dict[str, dict[str, list[float]]]]:
    return {model: {_r_key(r): {} for r in R_GRID} for model in MODEL_ROOTS}


def _is_strict_vetoing_result(result: dict[str, Any]) -> bool:
    tag = str(result.get("tag") or "unknown")
    return (
        bool(result.get("active"))
        and bool(result.get("evaluated"))
        and str(result.get("severity")) == "HARD"
        and result.get("passes") is False
        and tag in ("rigorous", "proxy")
    )


def _compact_rep_row(row: dict[str, Any]) -> dict[str, Any]:
    params = row["params"]
    return {
        "seed": row["seed"],
        "quark_fit_r": float(row["quark_fit_r"]),
        "params": {
            "Lambda_IR": float(params["Lambda_IR"]),
            "M_KK": float(params["M_KK"]),
            "k": float(params["k"]),
            "quark_yukawa_seed": params["quark_yukawa_seed"],
        },
        "fit_diagnostics": {
            "fitted_up_yukawa_singular_values": row["fit_diagnostics"][
                "fitted_up_yukawa_singular_values"
            ],
            "fitted_down_yukawa_singular_values": row["fit_diagnostics"][
                "fitted_down_yukawa_singular_values"
            ],
        },
    }


def _maybe_take_rep(
    candidates: dict[tuple[str, float], RepCandidate],
    *,
    model: str,
    row: dict[str, Any],
    r_value: float,
    mkk_tev: int,
) -> None:
    if not bool(row.get("survives_all_HARD_strict")):
        return
    for label, target_model, target_r in REP_TARGETS:
        if model != target_model or not math.isclose(r_value, target_r, abs_tol=1.0e-12):
            continue
        try:
            rank = REP_MKK_PREFERENCE.index(mkk_tev)
        except ValueError:
            rank = len(REP_MKK_PREFERENCE)
        key = (target_model, target_r)
        current = candidates.get(key)
        if current is None or rank < current.rank:
            candidates[key] = RepCandidate(
                label=label,
                ew_model=model,
                r=target_r,
                mkk_tev=mkk_tev,
                rank=rank,
                row=_compact_rep_row(row),
            )


def _percentiles(values: array) -> list[float]:
    if not values:
        return []
    data = np.frombuffer(values, dtype=np.float64)
    return [_round_sig(x, sig=6) for x in np.quantile(data, [0.25, 0.5, 0.75])]


def _stream_scan_roots() -> tuple[
    dict[tuple[str, str, int], int],
    dict[tuple[str, str, int], int],
    dict[tuple[str, str, int], Counter[str]],
    dict[tuple[str, str, int], int],
    dict[tuple[str, str, int], SingularCell],
    dict[str, Counter[str]],
    dict[tuple[str, float], RepCandidate],
    dict[str, int],
]:
    raw_counts: dict[tuple[str, str, int], int] = defaultdict(int)
    evaluated_counts: dict[tuple[str, str, int], int] = defaultdict(int)
    veto_counts: dict[tuple[str, str, int], Counter[str]] = defaultdict(Counter)
    joint_veto_counts: dict[tuple[str, str, int], int] = defaultdict(int)
    singular_cells: dict[tuple[str, str, int], SingularCell] = defaultdict(SingularCell)
    tag_counts: dict[str, Counter[str]] = defaultdict(Counter)
    rep_candidates: dict[tuple[str, float], RepCandidate] = {}
    rows_by_model: dict[str, int] = {}

    for model, root in MODEL_ROOTS.items():
        paths = _iter_jsonl_paths(root)
        if not paths:
            raise FileNotFoundError(f"no JSONL scan tiles found under {root}")

        model_rows = 0
        model_evaluated = 0
        for path_idx, path in enumerate(paths, start=1):
            with path.open("rb") as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    row = _json_loads(line)
                    model_rows += 1
                    r_value = float(row.get("quark_fit_r", row["params"].get("quark_fit_r")))
                    r_key = _r_key(r_value)
                    mkk_tev = _mkk_tev_from_row(row)
                    key = (model, r_key, mkk_tev)
                    raw_counts[key] += 1

                    if row.get("skipped"):
                        continue

                    model_evaluated += 1
                    evaluated_counts[key] += 1

                    diagnostics = row.get("fit_diagnostics") or {}
                    up_singular = diagnostics.get("fitted_up_yukawa_singular_values")
                    down_singular = diagnostics.get("fitted_down_yukawa_singular_values")
                    if up_singular and down_singular:
                        singular_cells[key].append(up_singular, down_singular)

                    constraints = row.get("constraints") or {}
                    row_joint_vetoed = False
                    for cid, result in constraints.items():
                        tag = str(result.get("tag") or "unknown")
                        tag_counts[cid][tag] += 1
                        # Count a veto only for a tag that can actually veto in
                        # the harness/comparison semantics (tag in {rigorous,
                        # proxy}); a partial/stub HARD failure never vetoes
                        # strict OR inclusive floors, so it must not appear in
                        # the explorer's veto fractions either (slice-5 F7).
                        if _is_strict_vetoing_result(result):
                            veto_counts[key][cid] += 1
                            row_joint_vetoed = True
                    if row_joint_vetoed:
                        joint_veto_counts[key] += 1

                    _maybe_take_rep(
                        rep_candidates,
                        model=model,
                        row=row,
                        r_value=r_value,
                        mkk_tev=mkk_tev,
                    )

            if path_idx % 100 == 0:
                print(
                    f"streamed {model}: {path_idx}/{len(paths)} files, "
                    f"{model_rows:,} rows",
                    file=sys.stderr,
                    flush=True,
                )

        rows_by_model[model] = model_rows
        print(
            f"streamed {model}: {model_rows:,} rows, {model_evaluated:,} evaluated",
            file=sys.stderr,
            flush=True,
        )

    return (
        raw_counts,
        evaluated_counts,
        veto_counts,
        joint_veto_counts,
        singular_cells,
        tag_counts,
        rep_candidates,
        rows_by_model,
    )


def _build_constraints(
    *,
    evaluated_counts: dict[tuple[str, str, int], int],
    veto_counts: dict[tuple[str, str, int], Counter[str]],
    tag_counts: dict[str, Counter[str]],
) -> list[dict[str, Any]]:
    contributing_ids: set[str] = set()
    for cid in tag_counts:
        for model in MODEL_ROOTS:
            for r_value in R_GRID:
                r_key = _r_key(r_value)
                for mkk_tev in MKK_GRID_TEV:
                    key = (model, r_key, mkk_tev)
                    denom = evaluated_counts.get(key, 0)
                    if denom <= 0:
                        continue
                    if veto_counts.get(key, Counter()).get(cid, 0) > 0:
                        contributing_ids.add(cid)
                        break
                if cid in contributing_ids:
                    break
            if cid in contributing_ids:
                break

    constraints = []
    for cid in sorted(contributing_ids):
        entry = _load_entry(cid)
        constraints.append(
            {
                "id": cid,
                "label": _constraint_label(cid, entry),
                "tag": _tag_for(cid, tag_counts),
                "group": _constraint_group(cid, entry),
                "default_on": True,
            }
        )

    group_order = {"electroweak": 0, "meson_mixing": 1, "collider": 2, "edm": 3, "other": 4}
    constraints.sort(key=lambda item: (group_order.get(item["group"], 9), item["id"]))
    return constraints


def _build_joint_veto_tree(
    *,
    evaluated_counts: dict[tuple[str, str, int], int],
    joint_veto_counts: dict[tuple[str, str, int], int],
) -> dict[str, dict[str, list[float]]]:
    joint_veto: dict[str, dict[str, list[float]]] = {
        model: {_r_key(r): [] for r in R_GRID} for model in MODEL_ROOTS
    }
    for model in MODEL_ROOTS:
        for r_value in R_GRID:
            r_key = _r_key(r_value)
            values = []
            for mkk_tev in MKK_GRID_TEV:
                key = (model, r_key, mkk_tev)
                denom = evaluated_counts.get(key, 0)
                fraction = 0.0 if denom <= 0 else joint_veto_counts.get(key, 0) / denom
                values.append(_round_fixed(fraction, digits=6))
            joint_veto[model][r_key] = values
    return joint_veto


def _build_veto_tree(
    constraints: list[dict[str, Any]],
    *,
    evaluated_counts: dict[tuple[str, str, int], int],
    veto_counts: dict[tuple[str, str, int], Counter[str]],
) -> dict[str, dict[str, dict[str, list[float]]]]:
    veto = _empty_veto_tree()
    for model in MODEL_ROOTS:
        for r_value in R_GRID:
            r_key = _r_key(r_value)
            for constraint in constraints:
                cid = constraint["id"]
                values = []
                for mkk_tev in MKK_GRID_TEV:
                    key = (model, r_key, mkk_tev)
                    denom = evaluated_counts.get(key, 0)
                    fraction = 0.0 if denom <= 0 else veto_counts.get(key, Counter()).get(cid, 0) / denom
                    values.append(_round_fixed(fraction, digits=6))
                veto[model][r_key][cid] = values
    return veto


def _build_bare_floor(evaluated_counts: dict[tuple[str, str, int], int]) -> dict[str, dict[str, int | None]]:
    bare_floor: dict[str, dict[str, int | None]] = {model: {} for model in MODEL_ROOTS}
    for model in MODEL_ROOTS:
        for r_value in R_GRID:
            r_key = _r_key(r_value)
            floor = None
            for mkk_tev in MKK_GRID_TEV:
                if evaluated_counts.get((model, r_key, mkk_tev), 0) > 0:
                    floor = mkk_tev
                    break
            bare_floor[model][r_key] = floor
    return bare_floor


def _build_yukawa_tree(
    singular_cells: dict[tuple[str, str, int], SingularCell]
) -> dict[str, dict[str, dict[str, dict[str, Any]]]]:
    yukawa: dict[str, dict[str, dict[str, dict[str, Any]]]] = {
        model: {_r_key(r): {} for r in R_GRID} for model in MODEL_ROOTS
    }
    for model in MODEL_ROOTS:
        for r_value in R_GRID:
            r_key = _r_key(r_value)
            for mkk_tev in MKK_GRID_TEV:
                key = (model, r_key, mkk_tev)
                cell = singular_cells.get(key)
                if cell is None or cell.n == 0:
                    continue
                yukawa[model][r_key][_mkk_key(mkk_tev)] = {
                    "up": [_percentiles(values) for values in cell.up],
                    "down": [_percentiles(values) for values in cell.down],
                    "n": cell.n,
                }
    return yukawa


def _matrix_from_seed(seed_payload: dict[str, Any], prefix: str) -> np.ndarray:
    return np.asarray(seed_payload[f"{prefix}_re"], dtype=np.float64) + 1j * np.asarray(
        seed_payload[f"{prefix}_im"], dtype=np.float64
    )


def _reconstruct_rep_matrices(
    rep_candidates: dict[tuple[str, float], RepCandidate]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    reps: list[dict[str, Any]] = []
    residuals: list[dict[str, Any]] = []
    targets = default_quark_targets()

    for label, model, r_value in REP_TARGETS:
        candidate = rep_candidates.get((model, r_value))
        if candidate is None:
            raise RuntimeError(f"no surviving representative found for {model} r={r_value}")

        row = candidate.row
        seed_payload = row["params"]["quark_yukawa_seed"]
        y_u_seed = _matrix_from_seed(seed_payload, "Y_u")
        y_d_seed = _matrix_from_seed(seed_payload, "Y_d")
        up_s, up_left, up_right = _svd_seed_parts(y_u_seed)
        down_s, down_left, down_right = _svd_seed_parts(y_d_seed)
        seed = QuarkFitSeed(
            up_singular_values=up_s,
            down_singular_values=down_s,
            overall_scale=1.0,
            up_left=up_left,
            up_right=up_right,
            down_left=down_left,
            down_right=down_right,
        )
        solution = fit_quark_sector(
            targets,
            r=float(row["quark_fit_r"]),
            seed=seed,
            Lambda_IR=float(row["params"]["Lambda_IR"]),
            k=float(row["params"]["k"]),
            max_nfev=DEFAULT_QUARK_FIT_MAX_NFEV,
            fit_orientation=True,
        )

        y_u = np.asarray(solution.result.point.Y_u, dtype=np.complex128)
        y_d = np.asarray(solution.result.point.Y_d, dtype=np.complex128)
        up_svd = np.sort(np.linalg.svd(y_u, compute_uv=False))
        down_svd = np.sort(np.linalg.svd(y_d, compute_uv=False))
        stored_up = np.asarray(
            row["fit_diagnostics"]["fitted_up_yukawa_singular_values"], dtype=np.float64
        )
        stored_down = np.asarray(
            row["fit_diagnostics"]["fitted_down_yukawa_singular_values"], dtype=np.float64
        )
        up_residual = float(np.max(np.abs(up_svd - stored_up)))
        down_residual = float(np.max(np.abs(down_svd - stored_down)))
        if up_residual >= 1.0e-6 or down_residual >= 1.0e-6:
            raise AssertionError(
                f"representative SVD mismatch for {label}: "
                f"up={up_residual:.3e}, down={down_residual:.3e}"
            )

        reps.append(
            {
                "label": label,
                "ew_model": candidate.ew_model,
                "r": candidate.r,
                "mkk_tev": candidate.mkk_tev,
                "seed": row["seed"],
                "Yu_abs": _round_matrix_abs(y_u, sig=4),
                "Yd_abs": _round_matrix_abs(y_d, sig=4),
                "up_singular": [_round_sig(x, sig=6) for x in stored_up],
                "down_singular": [_round_sig(x, sig=6) for x in stored_down],
            }
        )
        residuals.append(
            {
                "label": label,
                "up": up_residual,
                "down": down_residual,
                "seed": row["seed"],
                "mkk_tev": candidate.mkk_tev,
            }
        )

    return reps, residuals


def _load_comparison_rows() -> dict[tuple[str, str, str, str], float]:
    path = CUSTODIAL_ROOT / "comparison/constraint_veto_by_r_mkk.csv"
    rows: dict[tuple[str, str, str, str], float] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows[
                (
                    row["ew_model"],
                    str(float(row["r"])),
                    _mkk_key(float(row["mkk_tev"])),
                    row["constraint_id"],
                )
            ] = float(row["veto_fraction"])
    return rows


def _cross_checks(
    constraints: list[dict[str, Any]],
    veto: dict[str, dict[str, dict[str, list[float]]]],
) -> list[dict[str, Any]]:
    comparison = _load_comparison_rows()
    included = {item["id"] for item in constraints}
    checks: list[dict[str, Any]] = []

    preferred = [
        ("custodial", 0.05, 1, "EW001"),
        ("custodial", 0.1, 1, "B003"),
        ("custodial", 0.25, 1, "B013"),
        ("custodial", 1.0, 2, "B004"),
        ("minimal", 0.1, 1, "B003"),
    ]
    for model, r_value, mkk_tev, cid in preferred:
        if cid not in included:
            continue
        r_key = _r_key(r_value)
        try:
            idx = MKK_GRID_TEV.index(mkk_tev)
        except ValueError:
            continue
        csv_value = comparison.get((RAW_EW_MODEL[model], r_key, _mkk_key(mkk_tev), cid))
        if csv_value is None:
            continue
        generated = veto[model][r_key][cid][idx]
        checks.append(
            {
                "ew_model": model,
                "r": r_value,
                "mkk_tev": mkk_tev,
                "constraint_id": cid,
                "generated": generated,
                "comparison_csv": _round_fixed(csv_value, digits=6),
                "abs_diff": _round_sig(abs(generated - csv_value), sig=3),
            }
        )
    return checks


def _write_summary(
    *,
    size_bytes: int,
    constraints: list[dict[str, Any]],
    cross_checks: list[dict[str, Any]],
    residuals: list[dict[str, Any]],
    rows_by_model: dict[str, int],
) -> None:
    lines = [
        "# Scan Explorer Data Implementation Summary",
        "",
        f"- Built `{OUTPUT_PATH.relative_to(REPO_ROOT)}` ({size_bytes:,} bytes).",
        f"- Streamed rows: "
        + ", ".join(f"{model}={count:,}" for model, count in rows_by_model.items())
        + ".",
        f"- Included {len(constraints)} constraints: "
        + ", ".join(item["id"] for item in constraints)
        + ".",
        "- Cross-checks against `comparison/constraint_veto_by_r_mkk.csv`:",
    ]
    for check in cross_checks:
        lines.append(
            "  - {constraint_id} {ew_model} r={r} M_KK={mkk_tev} TeV: "
            "json={generated}, csv={comparison_csv}, abs_diff={abs_diff}".format(**check)
        )
    lines.append("- Representative matrix SVD residuals:")
    for residual in residuals:
        lines.append(
            "  - {label} seed={seed} M_KK={mkk_tev} TeV: "
            "up={up:.3e}, down={down:.3e}".format(**residual)
        )
    lines.extend(["", "DATA-READY"])
    SUMMARY_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build() -> dict[str, Any]:
    (
        raw_counts,
        evaluated_counts,
        veto_counts,
        joint_veto_counts,
        singular_cells,
        tag_counts,
        rep_candidates,
        rows_by_model,
    ) = _stream_scan_roots()

    constraints = _build_constraints(
        evaluated_counts=evaluated_counts,
        veto_counts=veto_counts,
        tag_counts=tag_counts,
    )
    veto = _build_veto_tree(
        constraints,
        evaluated_counts=evaluated_counts,
        veto_counts=veto_counts,
    )
    joint_veto = _build_joint_veto_tree(
        evaluated_counts=evaluated_counts,
        joint_veto_counts=joint_veto_counts,
    )
    yukawa = _build_yukawa_tree(singular_cells)
    rep_matrices, residuals = _reconstruct_rep_matrices(rep_candidates)

    raw_cell_sizes = set(raw_counts.values())
    n_draws_per_cell = max(raw_cell_sizes) if raw_cell_sizes else None
    artifact = {
        "meta": {
            "r_grid": R_GRID,
            "mkk_grid_tev": MKK_GRID_TEV,
            "minimal_root": str(MINIMAL_ROOT.relative_to(REPO_ROOT)),
            "custodial_root": str(CUSTODIAL_ROOT.relative_to(REPO_ROOT)),
            "floor_threshold": FLOOR_THRESHOLD,
            "envelope_floor_policy": (
                "default_all_constraints_uses_true_joint_veto_fraction_any_active_"
                "hard_rigorous_proxy_failure; custom browser subsets use stored "
                "per-constraint curves as an optimistic lower-bound proxy"
            ),
            "mkk_convention": "physical first KK mass m1 = 2.45 * Lambda_IR",
            "n_draws_per_cell": n_draws_per_cell,
            "generated_from_git": _git_sha(),
        },
        "constraints": constraints,
        "joint_veto": joint_veto,
        "veto": veto,
        "bare_floor_tev": _build_bare_floor(evaluated_counts),
        "yukawa": yukawa,
        "rep_matrices": rep_matrices,
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(
        json.dumps(artifact, ensure_ascii=True, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )

    size_bytes = OUTPUT_PATH.stat().st_size
    cross_checks = _cross_checks(constraints, veto)
    _write_summary(
        size_bytes=size_bytes,
        constraints=constraints,
        cross_checks=cross_checks,
        residuals=residuals,
        rows_by_model=rows_by_model,
    )

    print(f"JSON size: {size_bytes:,} bytes")
    print("constraints: " + ", ".join(item["id"] for item in constraints))
    print("cross-checks:")
    for check in cross_checks:
        print(
            "  {constraint_id} {ew_model} r={r} M_KK={mkk_tev} TeV: "
            "json={generated} csv={comparison_csv} abs_diff={abs_diff}".format(**check)
        )
    print("rep-matrix SVD residuals:")
    for residual in residuals:
        print(
            "  {label}: up={up:.3e} down={down:.3e} "
            "seed={seed} M_KK={mkk_tev} TeV".format(**residual)
        )
    print(f"summary: {SUMMARY_PATH.relative_to(REPO_ROOT)}")
    return artifact


def main() -> None:
    build()


if __name__ == "__main__":
    main()
