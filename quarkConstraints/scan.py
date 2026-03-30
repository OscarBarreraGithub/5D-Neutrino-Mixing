"""Scan wrapper for the quark-sector MFV benchmark workflow."""

from __future__ import annotations

import csv
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .benchmarks import QuarkTargets, default_quark_targets, default_spurion_seed
from .fit import fit_quark_sector
from .model import BulkMassMap
from .proxies import summarize_flavor_diagnostics


CSV_COLUMNS = [
    "sample_index",
    "git_commit",
    "dirty_tree",
    "rng_seed_global",
    "r",
    "overall_scale",
    "Lambda_IR",
    "k",
    "fit_success",
    "fit_score",
    "residual_norm",
    "mass_residual_max",
    "ckm_residual_max",
    "proxy_h_rs",
    "down_alignment",
    "up_alignment",
    "alignment_ratio",
    "passes_all",
    "reject_reason",
    "c_Q",
    "c_u",
    "c_d",
    "F_Q",
    "F_u",
    "F_d",
    "masses_up",
    "masses_down",
    "ckm_observables",
]


def _as_1d_float_array(name: str, values: Sequence[float] | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"{name} must be a non-empty 1D array")
    return arr


def _run_git(args: List[str], repo_root: Path) -> str:
    proc = subprocess.run(
        ["git", *args],
        cwd=str(repo_root),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return ""
    return proc.stdout.strip()


def _resolve_git_metadata(enabled: bool) -> Tuple[str, Optional[bool]]:
    if not enabled:
        return "disabled", None
    repo_root = Path(__file__).resolve().parents[1]
    git_commit = _run_git(["rev-parse", "HEAD"], repo_root)
    if not git_commit:
        return "unknown", None
    status = _run_git(["status", "--porcelain"], repo_root)
    return git_commit, len(status) > 0


def _serialize_array(values: np.ndarray) -> str:
    return ";".join(f"{float(x):.12g}" for x in np.asarray(values, dtype=float))


@dataclass
class QuarkScanConfig:
    """Configuration for a deterministic quark-sector MFV scan."""

    r_values: np.ndarray = field(default_factory=lambda: np.array([0.10, 0.25, 0.40], dtype=float))
    overall_scale_values: np.ndarray = field(default_factory=lambda: np.array([3.0], dtype=float))
    Lambda_IR_values: np.ndarray = field(default_factory=lambda: np.array([3000.0], dtype=float))
    k: float = 1.2209e19
    targets: QuarkTargets = field(default_factory=default_quark_targets)
    bulk_mass_map: BulkMassMap = field(default_factory=BulkMassMap)
    max_nfev: int = 120
    fit_orientation: bool = True
    max_mass_log_residual: float = 0.10
    max_ckm_relative_residual: float = 0.10
    max_alignment_ratio: float = 5.0
    record_git_metadata: bool = True
    rng_seed_global: Optional[int] = None

    def __post_init__(self) -> None:
        self.r_values = _as_1d_float_array("r_values", self.r_values)
        self.overall_scale_values = _as_1d_float_array(
            "overall_scale_values", self.overall_scale_values
        )
        self.Lambda_IR_values = _as_1d_float_array("Lambda_IR_values", self.Lambda_IR_values)
        if self.k <= 0.0:
            raise ValueError("k must be positive")
        if self.max_nfev <= 0:
            raise ValueError("max_nfev must be positive")


def _classify_solution(
    mass_log_residual: float,
    ckm_relative_residual: float,
    alignment_ratio: float,
    fit_success: bool,
    config: QuarkScanConfig,
) -> tuple[bool, str]:
    reasons = []
    if not fit_success:
        reasons.append("fit_failed")
    if mass_log_residual > config.max_mass_log_residual:
        reasons.append("mass_fit")
    if ckm_relative_residual > config.max_ckm_relative_residual:
        reasons.append("ckm_fit")
    if alignment_ratio > config.max_alignment_ratio:
        reasons.append("alignment")
    return len(reasons) == 0, ",".join(reasons) if reasons else "accepted"


def run_quark_scan(
    config: QuarkScanConfig,
    *,
    output_csv: str | None = None,
    progress_every: int = 100,
) -> List[Dict[str, object]]:
    """Run the quark-sector scan and optionally write a CSV file."""
    git_commit, dirty_tree = _resolve_git_metadata(config.record_git_metadata)
    rows: List[Dict[str, object]] = []
    handle = None
    writer = None
    if output_csv is not None:
        handle = open(output_csv, "w", newline="", encoding="utf-8")
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()

    try:
        sample_index = 0
        for Lambda_IR in config.Lambda_IR_values:
            for overall_scale in config.overall_scale_values:
                current_seed = default_spurion_seed()
                for r_val in config.r_values:
                    solution = fit_quark_sector(
                        config.targets,
                        r=float(r_val),
                        overall_scale=float(overall_scale),
                        seed=current_seed,
                        k=config.k,
                        Lambda_IR=float(Lambda_IR),
                        bulk_map=config.bulk_mass_map,
                        max_nfev=config.max_nfev,
                        fit_orientation=config.fit_orientation,
                    )
                    current_seed = solution.seed
                    result = solution.result
                    diagnostics = summarize_flavor_diagnostics(result, m_kk=float(Lambda_IR))
                    mass_log_residual = float(
                        np.max(
                            np.abs(
                                np.concatenate(
                                    [result.mass_residuals_up, result.mass_residuals_down]
                                )
                            )
                        )
                    )
                    ckm_relative_residual = float(np.max(np.abs(result.ckm_residuals)))
                    alignment_ratio = float(
                        diagnostics.diagnostics.down_offdiag_ratio_in_q_basis
                        / max(diagnostics.diagnostics.up_offdiag_ratio_in_q_basis, 1e-30)
                    )
                    passes_all, reject_reason = _classify_solution(
                        mass_log_residual,
                        ckm_relative_residual,
                        alignment_ratio,
                        solution.success,
                        config,
                    )
                    row: Dict[str, object] = {
                        "sample_index": sample_index,
                        "git_commit": git_commit,
                        "dirty_tree": dirty_tree,
                        "rng_seed_global": config.rng_seed_global,
                        "r": float(r_val),
                        "overall_scale": float(overall_scale),
                        "Lambda_IR": float(Lambda_IR),
                        "k": config.k,
                        "fit_success": solution.success,
                        "fit_score": result.score,
                        "residual_norm": result.residual_norm,
                        "mass_residual_max": mass_log_residual,
                        "ckm_residual_max": ckm_relative_residual,
                        "proxy_h_rs": diagnostics.h_rs_proxy,
                        "down_alignment": diagnostics.diagnostics.down_offdiag_ratio_in_q_basis,
                        "up_alignment": diagnostics.diagnostics.up_offdiag_ratio_in_q_basis,
                        "alignment_ratio": alignment_ratio,
                        "passes_all": passes_all,
                        "reject_reason": reject_reason,
                        "c_Q": _serialize_array(result.state.c_Q),
                        "c_u": _serialize_array(result.state.c_u),
                        "c_d": _serialize_array(result.state.c_d),
                        "F_Q": _serialize_array(result.state.F_Q),
                        "F_u": _serialize_array(result.state.F_u),
                        "F_d": _serialize_array(result.state.F_d),
                        "masses_up": _serialize_array(result.masses_up),
                        "masses_down": _serialize_array(result.masses_down),
                        "ckm_observables": _serialize_array(result.ckm_observables),
                    }
                    rows.append(row)
                    if writer is not None:
                        writer.writerow(row)
                    sample_index += 1
                    if progress_every and sample_index % progress_every == 0:
                        print(f"[quark scan] completed {sample_index} points")
    finally:
        if handle is not None:
            handle.close()

    return rows
