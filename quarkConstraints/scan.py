"""Scan wrapper for the quark-sector MFV benchmark workflow."""

from __future__ import annotations

import csv
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .benchmarks import (
    _FIXED_SCALE_TARGETS_PDG2024_MT_V1,
    QuarkTargets,
    default_quark_targets,
    default_spurion_seed,
)
from .couplings import compute_quark_kk_gluon_couplings
from .deltaf2 import evaluate_delta_f2_constraints
from .fit import fit_quark_sector
from .model import BulkMassMap
from .proxies import summarize_flavor_diagnostics
from .scales import (
    DEFAULT_QUARK_BENCHMARK_H_RS_MAX,
    DEFAULT_QUARK_XI_KK,
    default_quark_m_kk_from_lambda_ir,
)

# Per-flavor mass tolerances at mu_common = m_t(m_t).
#
# Initial production policy: use PDG 2024 2sigma relative uncertainties with
# a deterministic floor of 0.003 (per-mille minimum). After the user runs
# ``scripts/calibrate_phase0.py``, the floor will be replaced by the
# 95th-percentile log-residual from the Phase-0 dry scan (gating off).
#
# TODO(Phase-0 calibration): replace ``MASS_TOLERANCE_FLOOR`` with the
# per-flavor calibration loaded from
# ``data/phase0_residual_calibration.json`` once the scan has been run.
MASS_TOLERANCE_FLOOR = 0.003

# Per-element CKM 2sigma relative tolerances, derived from PDG 2024 §12
# global-fit averages of (|V_us|, |V_cb|, |V_ub|, J):
#
#   |V_us| = 0.22501 ± 0.00050   → 2sigma rel = 2 * 0.00050 / 0.22501 ≈ 0.0044
#   |V_cb| = 0.04183 ± 0.00150   → 2sigma rel = 2 * 0.00150 / 0.04183 ≈ 0.0717
#   |V_ub| = 0.00382 ± 0.00020   → 2sigma rel = 2 * 0.00020 / 0.00382 ≈ 0.1047
#   J      = 3.08e-5 ± 0.13e-5   → 2sigma rel = 2 * 0.13   / 3.08    ≈ 0.0844
#
# |V_ub| spans the inclusive/exclusive PDG-average envelope; 10% covers the
# full 2sigma band so neither tree-level extraction is silently rejected.
# Order matches ``QuarkTargets.ckm_observables``: (|V_us|, |V_cb|, |V_ub|, J).
_CKM_2SIGMA_RELATIVE_DEFAULT = np.array([0.0044, 0.072, 0.10, 0.085], dtype=float)


def _default_per_flavor_mass_tolerances() -> tuple[np.ndarray, np.ndarray]:
    """Return the per-flavor 2sigma relative tolerance vectors.

    Layout:
        up   = (m_u, m_c, m_t)
        down = (m_d, m_s, m_b)

    The values combine PDG 2024 2sigma relative uncertainties (from the
    PDG-derived target bundle) with a uniform floor (`MASS_TOLERANCE_FLOOR`).
    """
    bundle = _FIXED_SCALE_TARGETS_PDG2024_MT_V1
    up_pdg = np.asarray(bundle["up_2sigma_relative"], dtype=float)
    down_pdg = np.asarray(bundle["down_2sigma_relative"], dtype=float)
    up = np.maximum(up_pdg, MASS_TOLERANCE_FLOOR)
    down = np.maximum(down_pdg, MASS_TOLERANCE_FLOOR)
    return up, down

CSV_COLUMNS = [
    "sample_index",
    "git_commit",
    "dirty_tree",
    "rng_seed_global",
    "r",
    "overall_scale",
    "xi_KK",
    "Lambda_IR",
    "M_KK",
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
    "deltaf2_model",
    "deltaf2_input_bundle",
    "deltaf2_passes",
    "deltaf2_max_ratio",
    "epsilon_k_ratio",
    "epsilon_k_passes",
    "b_d_mix_ratio",
    "b_d_mix_passes",
    "b_s_mix_ratio",
    "b_s_mix_passes",
    "d_mix_ratio",
    "d_mix_passes",
    "passes_all",
    "reject_reason",
    "bulk_mass_map_model",
    "bulk_mass_map_c_uv",
    "bulk_mass_map_c_ir",
    "bulk_mass_map_eigen_scale",
    "fit_parameterization",
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
    """Configuration for a deterministic quark-sector MFV scan.

    Acceptance gating uses two layers of mass / CKM tolerances:

    * Per-quark vectors (``mass_tolerance_up``, ``mass_tolerance_down``)
      and per-element CKM tolerance (``ckm_tolerance``) — these are the
      production gates from the plan v3 implementation.
    * Legacy scalars (``max_mass_log_residual``, ``max_ckm_relative_residual``)
      — kept for backwards-compatibility with the legacy 0.10 gate. A
      candidate must clear *both* the per-quark and the legacy thresholds.

    Set ``apply_acceptance_gate=False`` to disable mass/CKM/proxy/alignment
    gating entirely (used by the Phase-0 calibration sweep).
    """

    r_values: np.ndarray = field(default_factory=lambda: np.array([0.10, 0.25, 0.40], dtype=float))
    overall_scale_values: np.ndarray = field(default_factory=lambda: np.array([3.0], dtype=float))
    Lambda_IR_values: np.ndarray = field(default_factory=lambda: np.array([3000.0], dtype=float))
    xi_KK: float = DEFAULT_QUARK_XI_KK
    k: float = 1.2209e19
    targets: QuarkTargets = field(default_factory=default_quark_targets)
    bulk_mass_map: BulkMassMap = field(default_factory=BulkMassMap)
    max_nfev: int = 120
    fit_orientation: bool = True
    max_mass_log_residual: float = 0.10
    max_ckm_relative_residual: float = 0.10
    mass_tolerance_up: np.ndarray = field(
        default_factory=lambda: _default_per_flavor_mass_tolerances()[0]
    )
    mass_tolerance_down: np.ndarray = field(
        default_factory=lambda: _default_per_flavor_mass_tolerances()[1]
    )
    ckm_tolerance: np.ndarray = field(
        default_factory=lambda: _CKM_2SIGMA_RELATIVE_DEFAULT.copy()
    )
    max_proxy_h_rs: float = DEFAULT_QUARK_BENCHMARK_H_RS_MAX
    max_alignment_ratio: float = 6.0
    record_git_metadata: bool = True
    rng_seed_global: Optional[int] = None
    apply_acceptance_gate: bool = True

    def __post_init__(self) -> None:
        self.r_values = _as_1d_float_array("r_values", self.r_values)
        self.overall_scale_values = _as_1d_float_array(
            "overall_scale_values", self.overall_scale_values
        )
        self.Lambda_IR_values = _as_1d_float_array("Lambda_IR_values", self.Lambda_IR_values)
        if self.k <= 0.0:
            raise ValueError("k must be positive")
        if self.xi_KK <= 0.0:
            raise ValueError("xi_KK must be positive")
        if self.max_nfev <= 0:
            raise ValueError("max_nfev must be positive")
        if self.max_proxy_h_rs <= 0.0:
            raise ValueError("max_proxy_h_rs must be positive")
        if np.any(self.Lambda_IR_values <= 0.0):
            raise ValueError("Lambda_IR_values must be positive")
        self.mass_tolerance_up = _as_1d_float_array("mass_tolerance_up", self.mass_tolerance_up)
        self.mass_tolerance_down = _as_1d_float_array(
            "mass_tolerance_down", self.mass_tolerance_down
        )
        self.ckm_tolerance = _as_1d_float_array("ckm_tolerance", self.ckm_tolerance)
        if self.mass_tolerance_up.size != 3:
            raise ValueError("mass_tolerance_up must have shape (3,)")
        if self.mass_tolerance_down.size != 3:
            raise ValueError("mass_tolerance_down must have shape (3,)")
        if self.ckm_tolerance.size != 4:
            raise ValueError("ckm_tolerance must have shape (4,)")
        if np.any(self.mass_tolerance_up <= 0.0):
            raise ValueError("mass_tolerance_up must be strictly positive")
        if np.any(self.mass_tolerance_down <= 0.0):
            raise ValueError("mass_tolerance_down must be strictly positive")
        if np.any(self.ckm_tolerance <= 0.0):
            raise ValueError("ckm_tolerance must be strictly positive")
        if self.rng_seed_global is not None:
            raise ValueError(
                "rng_seed_global is not yet supported for stochastic seeding; "
                "leave it unset until randomized scan initialization is implemented"
            )


_UP_FLAVOR_LABELS = ("u", "c", "t")
_DOWN_FLAVOR_LABELS = ("d", "s", "b")
_CKM_OBSERVABLE_LABELS = ("Vus", "Vcb", "Vub", "J")


def _per_quark_mass_failures(
    up_log_residuals: np.ndarray,
    down_log_residuals: np.ndarray,
    *,
    tol_up: np.ndarray,
    tol_down: np.ndarray,
) -> list[str]:
    failures: list[str] = []
    for label, residual, tol in zip(_UP_FLAVOR_LABELS, np.abs(up_log_residuals), tol_up):
        if residual > tol:
            failures.append(f"mass_{label}")
    for label, residual, tol in zip(_DOWN_FLAVOR_LABELS, np.abs(down_log_residuals), tol_down):
        if residual > tol:
            failures.append(f"mass_{label}")
    return failures


def _per_ckm_failures(
    ckm_residuals: np.ndarray,
    *,
    tol_ckm: np.ndarray,
) -> list[str]:
    failures: list[str] = []
    for label, residual, tol in zip(_CKM_OBSERVABLE_LABELS, np.abs(ckm_residuals), tol_ckm):
        if residual > tol:
            failures.append(f"ckm_{label}")
    return failures


def _classify_solution(
    mass_log_residual: float,
    ckm_relative_residual: float,
    proxy_h_rs: float,
    alignment_ratio: float,
    deltaf2_summary,
    fit_success: bool,
    config: QuarkScanConfig,
    *,
    up_log_residuals: Optional[np.ndarray] = None,
    down_log_residuals: Optional[np.ndarray] = None,
    ckm_residuals: Optional[np.ndarray] = None,
) -> tuple[bool, str]:
    """Apply the production acceptance gate.

    The function records gating failures in a fixed order; if
    ``config.apply_acceptance_gate`` is False the function still returns the
    full list of failure reasons (so Phase-0 calibration sweeps can log
    everything) but always reports ``passes = True``.
    """
    reasons: list[str] = []
    if not fit_success:
        reasons.append("fit_failed")
    if mass_log_residual > config.max_mass_log_residual:
        reasons.append("mass_fit")
    if ckm_relative_residual > config.max_ckm_relative_residual:
        reasons.append("ckm_fit")
    if up_log_residuals is not None and down_log_residuals is not None:
        reasons.extend(
            _per_quark_mass_failures(
                up_log_residuals,
                down_log_residuals,
                tol_up=config.mass_tolerance_up,
                tol_down=config.mass_tolerance_down,
            )
        )
    if ckm_residuals is not None:
        reasons.extend(
            _per_ckm_failures(ckm_residuals, tol_ckm=config.ckm_tolerance)
        )
    if proxy_h_rs > config.max_proxy_h_rs:
        reasons.append("proxy_h_rs")
    if alignment_ratio > config.max_alignment_ratio:
        reasons.append("alignment")
    system_labels = {
        "K": "epsilon_k",
        "B_d": "b_d_mix",
        "B_s": "b_s_mix",
        "D": "d_mix",
    }
    for system, label in system_labels.items():
        if not deltaf2_summary.by_system[system].passes:
            reasons.append(label)
    if not config.apply_acceptance_gate:
        # Phase-0 calibration: surface every reason but never gate.
        return True, "accepted_gate_disabled" if not reasons else (
            "accepted_gate_disabled:" + ",".join(reasons)
        )
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
            M_KK = default_quark_m_kk_from_lambda_ir(float(Lambda_IR), xi_KK=config.xi_KK)
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
                    diagnostics = summarize_flavor_diagnostics(result, m_kk=M_KK)
                    deltaf2 = evaluate_delta_f2_constraints(
                        # perturbative g_s (legacy repo_v1 behavior)
                        compute_quark_kk_gluon_couplings(result, M_KK=M_KK, xi_KK=config.xi_KK, g_s_star=None),
                        M_KK=M_KK,
                    )
                    deltaf2_by_system = deltaf2.by_system
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
                        diagnostics.h_rs_proxy,
                        alignment_ratio,
                        deltaf2,
                        solution.success,
                        config,
                        up_log_residuals=result.mass_residuals_up,
                        down_log_residuals=result.mass_residuals_down,
                        ckm_residuals=result.ckm_residuals,
                    )
                    row: Dict[str, object] = {
                        "sample_index": sample_index,
                        "git_commit": git_commit,
                        "dirty_tree": dirty_tree,
                        "rng_seed_global": config.rng_seed_global,
                        "r": float(r_val),
                        "overall_scale": float(overall_scale),
                        "xi_KK": config.xi_KK,
                        "Lambda_IR": float(Lambda_IR),
                        "M_KK": M_KK,
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
                        "deltaf2_model": deltaf2.operator_convention,
                        "deltaf2_input_bundle": deltaf2.input_bundle,
                        "deltaf2_passes": deltaf2.passes_all,
                        "deltaf2_max_ratio": deltaf2.max_ratio_to_bound,
                        "epsilon_k_ratio": deltaf2_by_system["K"].ratio_to_bound,
                        "epsilon_k_passes": deltaf2_by_system["K"].passes,
                        "b_d_mix_ratio": deltaf2_by_system["B_d"].ratio_to_bound,
                        "b_d_mix_passes": deltaf2_by_system["B_d"].passes,
                        "b_s_mix_ratio": deltaf2_by_system["B_s"].ratio_to_bound,
                        "b_s_mix_passes": deltaf2_by_system["B_s"].passes,
                        "d_mix_ratio": deltaf2_by_system["D"].ratio_to_bound,
                        "d_mix_passes": deltaf2_by_system["D"].passes,
                        "passes_all": passes_all,
                        "reject_reason": reject_reason,
                        "bulk_mass_map_model": "saturating_eigenvalue_window_v1",
                        "bulk_mass_map_c_uv": config.bulk_mass_map.c_uv,
                        "bulk_mass_map_c_ir": config.bulk_mass_map.c_ir,
                        "bulk_mass_map_eigen_scale": config.bulk_mass_map.eigen_scale,
                        "fit_parameterization": (
                            "singular_values_plus_left_rotations"
                            if config.fit_orientation
                            else "singular_values_only_fixed_orientations"
                        ),
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
