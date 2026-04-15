"""Operational scan machinery for the modern quark lane.

This module builds deterministic, shardable, resumable scan plumbing around the
existing modern evaluation/artifact/verifier surface. It is intentionally
honest about scope: this is operational scan machinery, not a final modern
allowed-region claim while later modern-physics milestones remain open.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from ..benchmarks import default_quark_targets, default_spurion_seed
from ..fit import QuarkFitSeed, QuarkFitSolution, QuarkTargets, fit_quark_sector
from ..model import RotationParameters
from ..scales import DEFAULT_QUARK_XI_KK, default_quark_m_kk_from_lambda_ir
from .artifacts import (
    build_modern_point_artifact,
    write_modern_point_artifact,
)
from .bridge_artifacts import (
    build_modern_point_bridge_artifact,
    write_modern_point_bridge_artifact,
)
from .conventions import MODERN_LANE_ID
from .evaluation import evaluate_modern_point
from .inputs import MODERN_DEFAULT_INPUT_BUNDLE_ID, MODERN_DEFAULT_INPUT_PROVENANCE_ID
from .phenomenology import (
    MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS,
    MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID,
    build_modern_point_phenomenology_artifact,
    write_modern_point_phenomenology_artifact,
)

MODERN_SCAN_CONFIG_SCHEMA_ID = "quarkConstraints.modern.scan.config.v1"
MODERN_SCAN_POINT_SCHEMA_ID = "quarkConstraints.modern.scan.point.v1"
MODERN_SCAN_RESULT_SCHEMA_ID = "quarkConstraints.modern.scan.result.v1"
MODERN_SCAN_SHARD_MANIFEST_SCHEMA_ID = "quarkConstraints.modern.scan.shard_manifest.v1"
MODERN_SCAN_MERGE_MANIFEST_SCHEMA_ID = "quarkConstraints.modern.scan.merge_manifest.v1"
MODERN_SCAN_RUN_CONFIG_SCHEMA_ID = "quarkConstraints.modern.scan.run_config.v1"
MODERN_SCAN_VERIFICATION_SCHEMA_ID = "quarkConstraints.modern.scan.verification.v1"
MODERN_SCAN_SCHEMA_VERSION = 1
MODERN_SCAN_CLAIM_LEVEL = "operational_scan_only"
MODERN_SCAN_CLAIM_LEVEL_ID = MODERN_SCAN_CLAIM_LEVEL
MODERN_SCAN_ROW_SCHEMA_ID = MODERN_SCAN_RESULT_SCHEMA_ID
MODERN_SCAN_PHYSICS_BACKEND_ID = "modern_point_eval_over_repo_v1_deltaf2"
MODERN_SCAN_NOTES = (
    "Operational scan machinery only. This module enumerates points, runs fits, "
    "exports artifacts, and executes verifier subprocesses, but it does not by "
    "itself upgrade the repo to a final modern allowed-region claim while "
    "broader claim-bearing physics milestones remain open."
)

_REPO_ROOT = Path(__file__).resolve().parents[2]


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n"


def _stable_hash(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _format_float(value: float) -> str:
    return format(float(value), ".12g")


def _require_positive_float(name: str, value: float) -> float:
    numeric = float(value)
    if not np.isfinite(numeric) or numeric <= 0.0:
        raise ValueError(f"{name} must be a positive finite float")
    return numeric


def _require_nonnegative_int(name: str, value: int) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a non-negative integer")
    return value


def _require_positive_int(name: str, value: int) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _as_sorted_unique_1d_float_array(
    name: str,
    values: Sequence[float] | np.ndarray,
) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"{name} must be a non-empty 1D array")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite floats")
    arr = np.unique(arr.astype(float))
    arr.sort()
    return arr


def _run_git(args: list[str], repo_root: Path) -> str:
    completed = subprocess.run(
        ["git", *args],
        cwd=str(repo_root),
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        return ""
    return completed.stdout.strip()


def _resolve_git_metadata(enabled: bool) -> tuple[str, bool | None]:
    if not enabled:
        return "disabled", None
    git_commit = _run_git(["rev-parse", "HEAD"], _REPO_ROOT)
    if not git_commit:
        return "unknown", None
    status = _run_git(["status", "--porcelain"], _REPO_ROOT)
    return git_commit, bool(status)


def _targets_payload(targets: QuarkTargets) -> dict[str, Any]:
    return {
        "label": targets.label,
        "up_masses": [float(value) for value in targets.up_masses],
        "down_masses": [float(value) for value in targets.down_masses],
        "ckm_real": [[float(value) for value in row] for row in np.real(targets.ckm)],
        "ckm_imag": [[float(value) for value in row] for row in np.imag(targets.ckm)],
    }


def _write_json(path: Path, payload: Mapping[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_canonical_json(payload), encoding="utf-8")
    return path


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> Path:
    """Write JSON atomically via temp file + rename to avoid race conditions."""
    import tempfile
    path.parent.mkdir(parents=True, exist_ok=True)
    content = _canonical_json(payload)
    fd, tmp = tempfile.mkstemp(dir=path.parent, suffix=".tmp", prefix=".config_")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(content)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
    return path


def _append_jsonl(path: Path, payload: Mapping[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True, allow_nan=False) + "\n")
    return path


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _parse_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.exists():
        return rows
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            try:
                row = json.loads(stripped)
            except json.JSONDecodeError as exc:
                raise ValueError(f"invalid JSON on line {line_number} of {path}") from exc
            if not isinstance(row, dict):
                raise ValueError(f"line {line_number} of {path} must decode to a JSON object")
            rows.append(row)
    return rows


def _display_path(path: Path, *, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path.resolve())


def _rotation_to_dict(rotation: RotationParameters) -> dict[str, float]:
    return {
        "theta12": float(rotation.theta12),
        "theta13": float(rotation.theta13),
        "theta23": float(rotation.theta23),
        "delta": float(rotation.delta),
    }


def _rotation_from_dict(payload: Mapping[str, Any]) -> RotationParameters:
    return RotationParameters(
        theta12=float(payload["theta12"]),
        theta13=float(payload["theta13"]),
        theta23=float(payload["theta23"]),
        delta=float(payload["delta"]),
    )


def seed_to_dict(seed: Any) -> dict[str, Any]:
    return {
        "up_singular_values": [float(value) for value in seed.up_singular_values],
        "down_singular_values": [float(value) for value in seed.down_singular_values],
        "overall_scale": float(seed.overall_scale),
        "up_left": _rotation_to_dict(seed.up_left),
        "up_right": _rotation_to_dict(seed.up_right),
        "down_left": _rotation_to_dict(seed.down_left),
        "down_right": _rotation_to_dict(seed.down_right),
    }


def seed_from_dict(payload: Mapping[str, Any]) -> QuarkFitSeed:
    return QuarkFitSeed(
        up_singular_values=np.asarray(payload["up_singular_values"], dtype=float),
        down_singular_values=np.asarray(payload["down_singular_values"], dtype=float),
        overall_scale=float(payload["overall_scale"]),
        up_left=_rotation_from_dict(payload["up_left"]),
        up_right=_rotation_from_dict(payload["up_right"]),
        down_left=_rotation_from_dict(payload["down_left"]),
        down_right=_rotation_from_dict(payload["down_right"]),
    )


@dataclass(frozen=True)
class ModernScanConfig:
    """Frozen deterministic config for the operational modern scan layer."""

    schema_id: str = MODERN_SCAN_CONFIG_SCHEMA_ID
    claim_level: str = MODERN_SCAN_CLAIM_LEVEL
    r_values: np.ndarray = field(
        default_factory=lambda: np.array([0.10, 0.25, 0.40], dtype=float)
    )
    overall_scale_values: np.ndarray = field(default_factory=lambda: np.array([3.0], dtype=float))
    Lambda_IR_values: np.ndarray = field(default_factory=lambda: np.array([3000.0], dtype=float))
    xi_KK: float = DEFAULT_QUARK_XI_KK
    k: float = 1.2209e19
    g_s_star: float | None = 3.0
    targets: QuarkTargets = field(default_factory=default_quark_targets)
    max_nfev: int = 120
    fit_orientation: bool = True
    record_git_metadata: bool = True
    notes: str = MODERN_SCAN_NOTES

    def __post_init__(self) -> None:
        if self.schema_id != MODERN_SCAN_CONFIG_SCHEMA_ID:
            raise ValueError(f"schema_id must be exactly {MODERN_SCAN_CONFIG_SCHEMA_ID!r}")
        if self.claim_level != MODERN_SCAN_CLAIM_LEVEL:
            raise ValueError(f"claim_level must be exactly {MODERN_SCAN_CLAIM_LEVEL!r}")
        object.__setattr__(
            self,
            "r_values",
            _as_sorted_unique_1d_float_array("r_values", self.r_values),
        )
        object.__setattr__(
            self,
            "overall_scale_values",
            _as_sorted_unique_1d_float_array("overall_scale_values", self.overall_scale_values),
        )
        object.__setattr__(
            self,
            "Lambda_IR_values",
            _as_sorted_unique_1d_float_array("Lambda_IR_values", self.Lambda_IR_values),
        )
        if np.any(self.r_values <= 0.0):
            raise ValueError("r_values must be strictly positive")
        if np.any(self.overall_scale_values <= 0.0):
            raise ValueError("overall_scale_values must be strictly positive")
        if np.any(self.Lambda_IR_values <= 0.0):
            raise ValueError("Lambda_IR_values must be strictly positive")
        object.__setattr__(self, "xi_KK", _require_positive_float("xi_KK", self.xi_KK))
        object.__setattr__(self, "k", _require_positive_float("k", self.k))
        if self.g_s_star is not None:
            g_s_star_val = float(self.g_s_star)
            if g_s_star_val <= 0.0:
                raise ValueError("g_s_star must be positive when provided")
            object.__setattr__(self, "g_s_star", g_s_star_val)
        if not isinstance(self.targets, QuarkTargets):
            raise ValueError("targets must be a QuarkTargets instance")
        object.__setattr__(self, "max_nfev", _require_positive_int("max_nfev", self.max_nfev))
        if not isinstance(self.fit_orientation, bool):
            raise ValueError("fit_orientation must be a bool")
        if not isinstance(self.record_git_metadata, bool):
            raise ValueError("record_git_metadata must be a bool")
        if not isinstance(self.notes, str) or not self.notes.strip():
            raise ValueError("notes must be a non-empty string")
        object.__setattr__(self, "notes", self.notes.strip())

    @property
    def total_points(self) -> int:
        return int(
            len(self.r_values) * len(self.overall_scale_values) * len(self.Lambda_IR_values)
        )

    @property
    def config_hash(self) -> str:
        return _stable_hash(self.as_dict())

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "claim_level": self.claim_level,
            "r_values": [float(value) for value in self.r_values],
            "overall_scale_values": [float(value) for value in self.overall_scale_values],
            "Lambda_IR_values": [float(value) for value in self.Lambda_IR_values],
            "xi_KK": self.xi_KK,
            "k": self.k,
            "g_s_star": self.g_s_star,
            "targets": _targets_payload(self.targets),
            "max_nfev": self.max_nfev,
            "fit_orientation": self.fit_orientation,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ModernScanConfig":
        targets_payload = dict(payload["targets"])
        ckm_real = np.asarray(targets_payload["ckm_real"], dtype=float)
        ckm_imag = np.asarray(targets_payload["ckm_imag"], dtype=float)
        # Preserve historical meaning: configs written before g_s* was
        # introduced used perturbative g_s (i.e. g_s_star=None).  Default
        # to None so we don't silently inject g_s*=3 into old configs.
        raw_g_s_star = payload.get("g_s_star", None)
        g_s_star = None if raw_g_s_star is None else float(raw_g_s_star)
        return cls(
            schema_id=str(payload["schema_id"]),
            claim_level=str(payload["claim_level"]),
            r_values=np.asarray(payload["r_values"], dtype=float),
            overall_scale_values=np.asarray(payload["overall_scale_values"], dtype=float),
            Lambda_IR_values=np.asarray(payload["Lambda_IR_values"], dtype=float),
            xi_KK=float(payload["xi_KK"]),
            k=float(payload["k"]),
            g_s_star=g_s_star,
            targets=QuarkTargets(
                up_masses=np.asarray(targets_payload["up_masses"], dtype=float),
                down_masses=np.asarray(targets_payload["down_masses"], dtype=float),
                ckm=ckm_real + 1j * ckm_imag,
                label=str(targets_payload["label"]),
            ),
            max_nfev=int(payload["max_nfev"]),
            fit_orientation=bool(payload["fit_orientation"]),
            notes=str(payload["notes"]),
        )


@dataclass(frozen=True)
class ModernScanPoint:
    """One deterministic scan point in the operational modern scan layer."""

    schema_id: str = MODERN_SCAN_POINT_SCHEMA_ID
    config_hash: str = ""
    point_index: int = 0
    point_id: str = ""
    point_label: str = ""
    r: float = 0.0
    overall_scale: float = 0.0
    Lambda_IR: float = 0.0
    M_KK: float = 0.0
    xi_KK: float = 0.0
    k: float = 0.0

    def __post_init__(self) -> None:
        if self.schema_id != MODERN_SCAN_POINT_SCHEMA_ID:
            raise ValueError(f"schema_id must be exactly {MODERN_SCAN_POINT_SCHEMA_ID!r}")
        if not isinstance(self.config_hash, str) or not self.config_hash:
            raise ValueError("config_hash must be a non-empty string")
        object.__setattr__(self, "point_index", _require_nonnegative_int("point_index", self.point_index))
        if not isinstance(self.point_id, str) or not self.point_id:
            raise ValueError("point_id must be a non-empty string")
        if not isinstance(self.point_label, str) or not self.point_label:
            raise ValueError("point_label must be a non-empty string")
        object.__setattr__(self, "r", float(self.r))
        object.__setattr__(self, "overall_scale", _require_positive_float("overall_scale", self.overall_scale))
        object.__setattr__(self, "Lambda_IR", _require_positive_float("Lambda_IR", self.Lambda_IR))
        object.__setattr__(self, "M_KK", _require_positive_float("M_KK", self.M_KK))
        object.__setattr__(self, "xi_KK", _require_positive_float("xi_KK", self.xi_KK))
        object.__setattr__(self, "k", _require_positive_float("k", self.k))

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "config_hash": self.config_hash,
            "point_index": self.point_index,
            "point_id": self.point_id,
            "point_label": self.point_label,
            "r": self.r,
            "overall_scale": self.overall_scale,
            "Lambda_IR": self.Lambda_IR,
            "M_KK": self.M_KK,
            "xi_KK": self.xi_KK,
            "k": self.k,
        }


@dataclass(frozen=True)
class ModernScanResultRecord:
    """One exported per-point scan record."""

    schema_id: str = MODERN_SCAN_RESULT_SCHEMA_ID
    claim_level: str = MODERN_SCAN_CLAIM_LEVEL
    lane_id: str = MODERN_LANE_ID
    config_hash: str = ""
    input_bundle_id: str = MODERN_DEFAULT_INPUT_BUNDLE_ID
    input_provenance_id: str = MODERN_DEFAULT_INPUT_PROVENANCE_ID
    git_commit: str = "disabled"
    dirty_tree: bool | None = None
    point_id: str = ""
    point_index: int = 0
    point_label: str = ""
    shard_index: int = 0
    shard_count: int = 1
    r: float = 0.0
    overall_scale: float = 0.0
    Lambda_IR: float = 0.0
    M_KK: float = 0.0
    xi_KK: float = 0.0
    k: float = 0.0
    fit_success: bool = False
    fit_message: str = ""
    fit_nfev: int = 0
    initial_score: float = 0.0
    fit_score: float = 0.0
    residual_norm: float = 0.0
    phenomenology_passes: bool = False
    verifier_ok: bool = False
    bridge_verifier_ok: bool = False
    phenomenology_verifier_ok: bool = False
    accepted: bool = False
    phenomenology_release_scope_id: str = MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID
    non_cp_acceptance_system_ids: tuple[str, ...] = (
        MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS
    )
    failing_non_cp_system_ids: tuple[str, ...] = ()
    non_cp_ratio_to_bound_by_system: Mapping[str, float] = field(default_factory=dict)
    max_non_cp_ratio_to_bound: float = 0.0
    diagnostic_only_system_ids: tuple[str, ...] = ()
    diagnostic_failing_system_ids: tuple[str, ...] = ()
    diagnostic_ratio_to_bound_by_system: Mapping[str, float] = field(default_factory=dict)
    blocked_system_ids: tuple[str, ...] = ()
    max_ratio_to_bound: float = 0.0
    failing_system_ids: tuple[str, ...] = ()
    ratio_to_bound_by_system: Mapping[str, float] = field(default_factory=dict)
    artifact_path: str = ""
    bridge_artifact_path: str = ""
    phenomenology_artifact_path: str = ""
    verifier_issue_codes: tuple[str, ...] = ()
    bridge_verifier_issue_codes: tuple[str, ...] = ()
    phenomenology_verifier_issue_codes: tuple[str, ...] = ()
    notes: str = MODERN_SCAN_NOTES

    def __post_init__(self) -> None:
        if self.schema_id != MODERN_SCAN_RESULT_SCHEMA_ID:
            raise ValueError(f"schema_id must be exactly {MODERN_SCAN_RESULT_SCHEMA_ID!r}")
        if self.claim_level != MODERN_SCAN_CLAIM_LEVEL:
            raise ValueError(f"claim_level must be exactly {MODERN_SCAN_CLAIM_LEVEL!r}")
        if self.lane_id != MODERN_LANE_ID:
            raise ValueError(f"lane_id must be exactly {MODERN_LANE_ID!r}")
        if not isinstance(self.config_hash, str) or not self.config_hash:
            raise ValueError("config_hash must be a non-empty string")
        if not isinstance(self.input_bundle_id, str) or not self.input_bundle_id:
            raise ValueError("input_bundle_id must be a non-empty string")
        if not isinstance(self.input_provenance_id, str) or not self.input_provenance_id:
            raise ValueError("input_provenance_id must be a non-empty string")
        if not isinstance(self.git_commit, str) or not self.git_commit:
            raise ValueError("git_commit must be a non-empty string")
        if not isinstance(self.point_id, str) or not self.point_id:
            raise ValueError("point_id must be a non-empty string")
        object.__setattr__(self, "point_index", _require_nonnegative_int("point_index", self.point_index))
        object.__setattr__(self, "shard_index", _require_nonnegative_int("shard_index", self.shard_index))
        object.__setattr__(self, "shard_count", _require_positive_int("shard_count", self.shard_count))
        if not isinstance(self.point_label, str) or not self.point_label:
            raise ValueError("point_label must be a non-empty string")
        if not isinstance(self.fit_message, str):
            raise ValueError("fit_message must be a string")
        object.__setattr__(self, "fit_nfev", _require_nonnegative_int("fit_nfev", self.fit_nfev))
        if not isinstance(self.artifact_path, str) or not self.artifact_path:
            raise ValueError("artifact_path must be a non-empty string")
        if not isinstance(self.bridge_artifact_path, str) or not self.bridge_artifact_path:
            raise ValueError("bridge_artifact_path must be a non-empty string")
        if (
            not isinstance(self.phenomenology_artifact_path, str)
            or not self.phenomenology_artifact_path
        ):
            raise ValueError("phenomenology_artifact_path must be a non-empty string")
        object.__setattr__(self, "verifier_issue_codes", tuple(str(code) for code in self.verifier_issue_codes))
        object.__setattr__(
            self,
            "bridge_verifier_issue_codes",
            tuple(str(code) for code in self.bridge_verifier_issue_codes),
        )
        object.__setattr__(
            self,
            "phenomenology_verifier_issue_codes",
            tuple(str(code) for code in self.phenomenology_verifier_issue_codes),
        )
        object.__setattr__(
            self,
            "phenomenology_release_scope_id",
            str(self.phenomenology_release_scope_id),
        )
        if (
            self.phenomenology_release_scope_id
            != MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID
        ):
            raise ValueError(
                "phenomenology_release_scope_id must match the frozen modern "
                "non-CP release scope"
            )
        object.__setattr__(
            self,
            "non_cp_acceptance_system_ids",
            tuple(str(system) for system in self.non_cp_acceptance_system_ids),
        )
        if (
            self.non_cp_acceptance_system_ids
            != MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS
        ):
            raise ValueError(
                "non_cp_acceptance_system_ids must be exactly "
                f"{MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS!r}"
            )
        object.__setattr__(
            self,
            "failing_non_cp_system_ids",
            tuple(str(system) for system in self.failing_non_cp_system_ids),
        )
        object.__setattr__(
            self,
            "diagnostic_only_system_ids",
            tuple(str(system) for system in self.diagnostic_only_system_ids),
        )
        object.__setattr__(
            self,
            "diagnostic_failing_system_ids",
            tuple(str(system) for system in self.diagnostic_failing_system_ids),
        )
        object.__setattr__(
            self,
            "blocked_system_ids",
            tuple(str(system) for system in self.blocked_system_ids),
        )
        object.__setattr__(self, "failing_system_ids", tuple(str(system) for system in self.failing_system_ids))
        if not isinstance(self.non_cp_ratio_to_bound_by_system, Mapping):
            raise ValueError("non_cp_ratio_to_bound_by_system must be a mapping")
        canonical_non_cp_ratios = {
            str(key): float(value)
            for key, value in sorted(self.non_cp_ratio_to_bound_by_system.items())
        }
        if tuple(canonical_non_cp_ratios) != tuple(
            system_id
            for system_id in self.non_cp_acceptance_system_ids
            if system_id in canonical_non_cp_ratios
        ):
            unexpected = tuple(
                key
                for key in canonical_non_cp_ratios
                if key not in self.non_cp_acceptance_system_ids
            )
            if unexpected:
                raise ValueError(
                    "non_cp_ratio_to_bound_by_system contains systems outside the "
                    "acceptance-bearing scope"
                )
        object.__setattr__(
            self,
            "non_cp_ratio_to_bound_by_system",
            canonical_non_cp_ratios,
        )
        expected_max_non_cp_ratio = max(canonical_non_cp_ratios.values(), default=0.0)
        object.__setattr__(
            self,
            "max_non_cp_ratio_to_bound",
            float(self.max_non_cp_ratio_to_bound),
        )
        if not np.isclose(
            self.max_non_cp_ratio_to_bound,
            expected_max_non_cp_ratio,
            rtol=0.0,
            atol=0.0,
        ):
            raise ValueError(
                "max_non_cp_ratio_to_bound must equal the maximum acceptance-bearing ratio"
            )
        if not isinstance(self.diagnostic_ratio_to_bound_by_system, Mapping):
            raise ValueError("diagnostic_ratio_to_bound_by_system must be a mapping")
        canonical_diagnostic_ratios = {
            str(key): float(value)
            for key, value in sorted(self.diagnostic_ratio_to_bound_by_system.items())
        }
        unexpected_diagnostics = tuple(
            key
            for key in canonical_diagnostic_ratios
            if key not in self.diagnostic_only_system_ids
        )
        if unexpected_diagnostics:
            raise ValueError(
                "diagnostic_ratio_to_bound_by_system contains systems outside the "
                "diagnostic-only scope"
            )
        object.__setattr__(
            self,
            "diagnostic_ratio_to_bound_by_system",
            canonical_diagnostic_ratios,
        )
        if not isinstance(self.ratio_to_bound_by_system, Mapping):
            raise ValueError("ratio_to_bound_by_system must be a mapping")
        canonical_ratios = {
            str(key): float(value)
            for key, value in sorted(self.ratio_to_bound_by_system.items())
        }
        object.__setattr__(self, "ratio_to_bound_by_system", canonical_ratios)
        expected_legacy_keys = tuple(
            sorted(
                set(self.non_cp_ratio_to_bound_by_system)
                | set(self.diagnostic_ratio_to_bound_by_system)
            )
        )
        if tuple(canonical_ratios) != expected_legacy_keys:
            raise ValueError(
                "ratio_to_bound_by_system must mirror the union of non-CP and "
                "diagnostic evaluated systems"
            )
        # Validate that failing lists are consistent with ratios (set
        # comparison to avoid ordering sensitivity).
        expected_failing_non_cp = frozenset(
            sid for sid, ratio in self.non_cp_ratio_to_bound_by_system.items()
            if ratio > 1.0
        )
        if frozenset(self.failing_non_cp_system_ids) != expected_failing_non_cp:
            raise ValueError(
                "failing_non_cp_system_ids must match the failing acceptance-bearing systems"
            )
        expected_failing_diag = frozenset(
            sid for sid, ratio in self.diagnostic_ratio_to_bound_by_system.items()
            if ratio > 1.0
        )
        if frozenset(self.diagnostic_failing_system_ids) != expected_failing_diag:
            raise ValueError(
                "diagnostic_failing_system_ids must match the failing diagnostic systems"
            )
        expected_failing_system_ids = tuple(
            system_id
            for system_id in ("epsilon_K", "K", "B_d", "B_s", "D0")
            if (
                system_id in self.failing_non_cp_system_ids
                or system_id in self.diagnostic_failing_system_ids
            )
        )
        if tuple(self.failing_system_ids) != expected_failing_system_ids:
            raise ValueError(
                "failing_system_ids must mirror the union of non-CP and diagnostic failures"
            )
        # No systems are blocked in the full Delta F = 2 release
        if self.blocked_system_ids != ():
            raise ValueError(
                "blocked_system_ids must be empty in the full Delta F = 2 release"
            )
        expected_max_ratio = max(canonical_ratios.values(), default=0.0)
        object.__setattr__(self, "max_ratio_to_bound", float(self.max_ratio_to_bound))
        if not np.isclose(
            self.max_ratio_to_bound,
            expected_max_ratio,
            rtol=0.0,
            atol=0.0,
        ):
            raise ValueError(
                "max_ratio_to_bound must equal the maximum evaluated ratio across "
                "acceptance-bearing and diagnostic systems"
            )
        if not isinstance(self.notes, str) or not self.notes.strip():
            raise ValueError("notes must be a non-empty string")
        object.__setattr__(self, "notes", self.notes.strip())

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "claim_level": self.claim_level,
            "lane_id": self.lane_id,
            "config_hash": self.config_hash,
            "input_bundle_id": self.input_bundle_id,
            "input_provenance_id": self.input_provenance_id,
            "git_commit": self.git_commit,
            "dirty_tree": self.dirty_tree,
            "point_id": self.point_id,
            "point_index": self.point_index,
            "point_label": self.point_label,
            "shard_index": self.shard_index,
            "shard_count": self.shard_count,
            "r": self.r,
            "overall_scale": self.overall_scale,
            "Lambda_IR": self.Lambda_IR,
            "M_KK": self.M_KK,
            "xi_KK": self.xi_KK,
            "k": self.k,
            "fit_success": self.fit_success,
            "fit_message": self.fit_message,
            "fit_nfev": self.fit_nfev,
            "initial_score": self.initial_score,
            "fit_score": self.fit_score,
            "residual_norm": self.residual_norm,
            "phenomenology_passes": self.phenomenology_passes,
            "verifier_ok": self.verifier_ok,
            "bridge_verifier_ok": self.bridge_verifier_ok,
            "phenomenology_verifier_ok": self.phenomenology_verifier_ok,
            "accepted": self.accepted,
            "phenomenology_release_scope_id": self.phenomenology_release_scope_id,
            "non_cp_acceptance_system_ids": list(self.non_cp_acceptance_system_ids),
            "failing_non_cp_system_ids": list(self.failing_non_cp_system_ids),
            "non_cp_ratio_to_bound_by_system": dict(self.non_cp_ratio_to_bound_by_system),
            "max_non_cp_ratio_to_bound": self.max_non_cp_ratio_to_bound,
            "diagnostic_only_system_ids": list(self.diagnostic_only_system_ids),
            "diagnostic_failing_system_ids": list(
                self.diagnostic_failing_system_ids
            ),
            "diagnostic_ratio_to_bound_by_system": dict(
                self.diagnostic_ratio_to_bound_by_system
            ),
            "blocked_system_ids": list(self.blocked_system_ids),
            "max_ratio_to_bound": self.max_ratio_to_bound,
            "failing_system_ids": list(self.failing_system_ids),
            "ratio_to_bound_by_system": dict(self.ratio_to_bound_by_system),
            "artifact_path": self.artifact_path,
            "bridge_artifact_path": self.bridge_artifact_path,
            "phenomenology_artifact_path": self.phenomenology_artifact_path,
            "verifier_issue_codes": list(self.verifier_issue_codes),
            "bridge_verifier_issue_codes": list(self.bridge_verifier_issue_codes),
            "phenomenology_verifier_issue_codes": list(
                self.phenomenology_verifier_issue_codes
            ),
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ModernScanResultRecord":
        return cls(
            schema_id=str(payload["schema_id"]),
            claim_level=str(payload["claim_level"]),
            lane_id=str(payload.get("lane_id", MODERN_LANE_ID)),
            config_hash=str(payload["config_hash"]),
            input_bundle_id=str(
                payload.get("input_bundle_id", MODERN_DEFAULT_INPUT_BUNDLE_ID)
            ),
            input_provenance_id=str(
                payload.get("input_provenance_id", MODERN_DEFAULT_INPUT_PROVENANCE_ID)
            ),
            git_commit=str(payload.get("git_commit", "disabled")),
            dirty_tree=payload.get("dirty_tree"),
            point_id=str(payload["point_id"]),
            point_index=int(payload["point_index"]),
            point_label=str(payload["point_label"]),
            shard_index=int(payload["shard_index"]),
            shard_count=int(payload["shard_count"]),
            r=float(payload["r"]),
            overall_scale=float(payload["overall_scale"]),
            Lambda_IR=float(payload["Lambda_IR"]),
            M_KK=float(payload["M_KK"]),
            xi_KK=float(payload["xi_KK"]),
            k=float(payload["k"]),
            fit_success=bool(payload["fit_success"]),
            fit_message=str(payload["fit_message"]),
            fit_nfev=int(payload["fit_nfev"]),
            initial_score=float(payload["initial_score"]),
            fit_score=float(payload["fit_score"]),
            residual_norm=float(payload["residual_norm"]),
            phenomenology_passes=bool(payload["phenomenology_passes"]),
            verifier_ok=bool(payload["verifier_ok"]),
            bridge_verifier_ok=bool(payload.get("bridge_verifier_ok", False)),
            phenomenology_verifier_ok=bool(
                payload.get("phenomenology_verifier_ok", False)
            ),
            accepted=bool(payload["accepted"]),
            phenomenology_release_scope_id=str(
                payload.get(
                    "phenomenology_release_scope_id",
                    MODERN_POINT_PHENOMENOLOGY_RELEASE_SCOPE_ID,
                )
            ),
            non_cp_acceptance_system_ids=tuple(
                str(item)
                for item in payload.get(
                    "non_cp_acceptance_system_ids",
                    MODERN_POINT_PHENOMENOLOGY_NON_CP_ACCEPTANCE_SYSTEM_IDS,
                )
            ),
            failing_non_cp_system_ids=tuple(
                str(item) for item in payload.get("failing_non_cp_system_ids", ())
            ),
            non_cp_ratio_to_bound_by_system={
                str(key): float(value)
                for key, value in dict(
                    payload.get("non_cp_ratio_to_bound_by_system", {})
                ).items()
            },
            max_non_cp_ratio_to_bound=float(
                payload.get("max_non_cp_ratio_to_bound", 0.0)
            ),
            diagnostic_only_system_ids=tuple(
                str(item) for item in payload.get("diagnostic_only_system_ids", ())
            ),
            diagnostic_failing_system_ids=tuple(
                str(item)
                for item in payload.get("diagnostic_failing_system_ids", ())
            ),
            diagnostic_ratio_to_bound_by_system={
                str(key): float(value)
                for key, value in dict(
                    payload.get("diagnostic_ratio_to_bound_by_system", {})
                ).items()
            },
            blocked_system_ids=tuple(
                str(item) for item in payload.get("blocked_system_ids", ())
            ),
            max_ratio_to_bound=float(payload["max_ratio_to_bound"]),
            failing_system_ids=tuple(str(item) for item in payload["failing_system_ids"]),
            ratio_to_bound_by_system={
                str(key): float(value)
                for key, value in dict(payload["ratio_to_bound_by_system"]).items()
            },
            artifact_path=str(payload["artifact_path"]),
            bridge_artifact_path=str(payload.get("bridge_artifact_path", "")),
            phenomenology_artifact_path=str(
                payload.get("phenomenology_artifact_path", "")
            ),
            verifier_issue_codes=tuple(str(item) for item in payload["verifier_issue_codes"]),
            bridge_verifier_issue_codes=tuple(
                str(item) for item in payload.get("bridge_verifier_issue_codes", ())
            ),
            phenomenology_verifier_issue_codes=tuple(
                str(item)
                for item in payload.get("phenomenology_verifier_issue_codes", ())
            ),
            notes=str(payload["notes"]),
        )


@dataclass(frozen=True)
class ModernScanShardManifest:
    """One shard-level manifest for a deterministic modern scan run."""

    schema_id: str = MODERN_SCAN_SHARD_MANIFEST_SCHEMA_ID
    claim_level: str = MODERN_SCAN_CLAIM_LEVEL
    lane_id: str = MODERN_LANE_ID
    config_hash: str = ""
    input_bundle_id: str = MODERN_DEFAULT_INPUT_BUNDLE_ID
    input_provenance_id: str = MODERN_DEFAULT_INPUT_PROVENANCE_ID
    shard_index: int = 0
    shard_count: int = 1
    assigned_point_count: int = 0
    completed_point_count: int = 0
    skipped_existing_point_count: int = 0
    accepted_point_count: int = 0
    verifier_failed_point_count: int = 0
    bridge_verifier_failed_point_count: int = 0
    phenomenology_verifier_failed_point_count: int = 0
    complete: bool = False
    results_path: str = ""
    artifact_dir: str = ""
    bridge_artifact_dir: str = ""
    phenomenology_artifact_dir: str = ""
    cache_dir: str = ""
    run_config_path: str = ""
    git_commit: str = "disabled"
    dirty_tree: bool | None = None
    notes: str = MODERN_SCAN_NOTES

    def __post_init__(self) -> None:
        if self.schema_id != MODERN_SCAN_SHARD_MANIFEST_SCHEMA_ID:
            raise ValueError(
                f"schema_id must be exactly {MODERN_SCAN_SHARD_MANIFEST_SCHEMA_ID!r}"
            )
        if self.claim_level != MODERN_SCAN_CLAIM_LEVEL:
            raise ValueError(f"claim_level must be exactly {MODERN_SCAN_CLAIM_LEVEL!r}")
        if self.lane_id != MODERN_LANE_ID:
            raise ValueError(f"lane_id must be exactly {MODERN_LANE_ID!r}")
        if not isinstance(self.config_hash, str) or not self.config_hash:
            raise ValueError("config_hash must be a non-empty string")
        if not isinstance(self.input_bundle_id, str) or not self.input_bundle_id:
            raise ValueError("input_bundle_id must be a non-empty string")
        if not isinstance(self.input_provenance_id, str) or not self.input_provenance_id:
            raise ValueError("input_provenance_id must be a non-empty string")
        object.__setattr__(self, "shard_index", _require_nonnegative_int("shard_index", self.shard_index))
        object.__setattr__(self, "shard_count", _require_positive_int("shard_count", self.shard_count))
        object.__setattr__(
            self,
            "assigned_point_count",
            _require_nonnegative_int("assigned_point_count", self.assigned_point_count),
        )
        object.__setattr__(
            self,
            "completed_point_count",
            _require_nonnegative_int("completed_point_count", self.completed_point_count),
        )
        object.__setattr__(
            self,
            "skipped_existing_point_count",
            _require_nonnegative_int(
                "skipped_existing_point_count",
                self.skipped_existing_point_count,
            ),
        )
        object.__setattr__(
            self,
            "accepted_point_count",
            _require_nonnegative_int("accepted_point_count", self.accepted_point_count),
        )
        object.__setattr__(
            self,
            "verifier_failed_point_count",
            _require_nonnegative_int(
                "verifier_failed_point_count",
                self.verifier_failed_point_count,
            ),
        )
        object.__setattr__(
            self,
            "bridge_verifier_failed_point_count",
            _require_nonnegative_int(
                "bridge_verifier_failed_point_count",
                self.bridge_verifier_failed_point_count,
            ),
        )
        object.__setattr__(
            self,
            "phenomenology_verifier_failed_point_count",
            _require_nonnegative_int(
                "phenomenology_verifier_failed_point_count",
                self.phenomenology_verifier_failed_point_count,
            ),
        )
        if not isinstance(self.complete, bool):
            raise ValueError("complete must be a bool")
        if not isinstance(self.results_path, str) or not self.results_path:
            raise ValueError("results_path must be a non-empty string")
        if not isinstance(self.artifact_dir, str) or not self.artifact_dir:
            raise ValueError("artifact_dir must be a non-empty string")
        if not isinstance(self.bridge_artifact_dir, str) or not self.bridge_artifact_dir:
            raise ValueError("bridge_artifact_dir must be a non-empty string")
        if (
            not isinstance(self.phenomenology_artifact_dir, str)
            or not self.phenomenology_artifact_dir
        ):
            raise ValueError("phenomenology_artifact_dir must be a non-empty string")
        if not isinstance(self.cache_dir, str) or not self.cache_dir:
            raise ValueError("cache_dir must be a non-empty string")
        if not isinstance(self.run_config_path, str) or not self.run_config_path:
            raise ValueError("run_config_path must be a non-empty string")
        if not isinstance(self.git_commit, str) or not self.git_commit:
            raise ValueError("git_commit must be a non-empty string")
        if self.completed_point_count > self.assigned_point_count:
            raise ValueError("completed_point_count cannot exceed assigned_point_count")
        if self.accepted_point_count > self.completed_point_count:
            raise ValueError("accepted_point_count cannot exceed completed_point_count")
        if self.verifier_failed_point_count > self.completed_point_count:
            raise ValueError("verifier_failed_point_count cannot exceed completed_point_count")
        if self.bridge_verifier_failed_point_count > self.completed_point_count:
            raise ValueError("bridge_verifier_failed_point_count cannot exceed completed_point_count")
        if self.phenomenology_verifier_failed_point_count > self.completed_point_count:
            raise ValueError(
                "phenomenology_verifier_failed_point_count cannot exceed completed_point_count"
            )
        if not isinstance(self.notes, str) or not self.notes.strip():
            raise ValueError("notes must be a non-empty string")
        object.__setattr__(self, "notes", self.notes.strip())

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "claim_level": self.claim_level,
            "lane_id": self.lane_id,
            "config_hash": self.config_hash,
            "input_bundle_id": self.input_bundle_id,
            "input_provenance_id": self.input_provenance_id,
            "shard_index": self.shard_index,
            "shard_count": self.shard_count,
            "assigned_point_count": self.assigned_point_count,
            "completed_point_count": self.completed_point_count,
            "skipped_existing_point_count": self.skipped_existing_point_count,
            "accepted_point_count": self.accepted_point_count,
            "verifier_failed_point_count": self.verifier_failed_point_count,
            "bridge_verifier_failed_point_count": self.bridge_verifier_failed_point_count,
            "phenomenology_verifier_failed_point_count": self.phenomenology_verifier_failed_point_count,
            "complete": self.complete,
            "results_path": self.results_path,
            "artifact_dir": self.artifact_dir,
            "bridge_artifact_dir": self.bridge_artifact_dir,
            "phenomenology_artifact_dir": self.phenomenology_artifact_dir,
            "cache_dir": self.cache_dir,
            "run_config_path": self.run_config_path,
            "git_commit": self.git_commit,
            "dirty_tree": self.dirty_tree,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ModernScanShardManifest":
        return cls(
            schema_id=str(payload["schema_id"]),
            claim_level=str(payload["claim_level"]),
            lane_id=str(payload["lane_id"]),
            config_hash=str(payload["config_hash"]),
            input_bundle_id=str(payload["input_bundle_id"]),
            input_provenance_id=str(payload["input_provenance_id"]),
            shard_index=int(payload["shard_index"]),
            shard_count=int(payload["shard_count"]),
            assigned_point_count=int(payload["assigned_point_count"]),
            completed_point_count=int(payload["completed_point_count"]),
            skipped_existing_point_count=int(payload["skipped_existing_point_count"]),
            accepted_point_count=int(payload["accepted_point_count"]),
            verifier_failed_point_count=int(payload["verifier_failed_point_count"]),
            bridge_verifier_failed_point_count=int(
                payload["bridge_verifier_failed_point_count"]
            ),
            phenomenology_verifier_failed_point_count=int(
                payload.get("phenomenology_verifier_failed_point_count", 0)
            ),
            complete=bool(payload["complete"]),
            results_path=str(payload["results_path"]),
            artifact_dir=str(payload["artifact_dir"]),
            bridge_artifact_dir=str(payload["bridge_artifact_dir"]),
            phenomenology_artifact_dir=str(
                payload.get("phenomenology_artifact_dir", "")
            ),
            cache_dir=str(payload["cache_dir"]),
            run_config_path=str(payload["run_config_path"]),
            git_commit=str(payload["git_commit"]),
            dirty_tree=payload.get("dirty_tree"),
            notes=str(payload["notes"]),
        )


@dataclass(frozen=True)
class ModernMergedScanManifest:
    """Merged reduction manifest across shard outputs."""

    schema_id: str = MODERN_SCAN_MERGE_MANIFEST_SCHEMA_ID
    claim_level: str = MODERN_SCAN_CLAIM_LEVEL
    lane_id: str = MODERN_LANE_ID
    config_hash: str = ""
    input_bundle_id: str = MODERN_DEFAULT_INPUT_BUNDLE_ID
    input_provenance_id: str = MODERN_DEFAULT_INPUT_PROVENANCE_ID
    shard_count: int = 1
    shard_manifest_count: int = 0
    total_point_count: int = 0
    accepted_point_count: int = 0
    verifier_failed_point_count: int = 0
    bridge_verifier_failed_point_count: int = 0
    phenomenology_verifier_failed_point_count: int = 0
    complete: bool = False
    merged_results_path: str = ""
    source_manifest_paths: tuple[str, ...] = ()
    notes: str = MODERN_SCAN_NOTES

    def __post_init__(self) -> None:
        if self.schema_id != MODERN_SCAN_MERGE_MANIFEST_SCHEMA_ID:
            raise ValueError(
                f"schema_id must be exactly {MODERN_SCAN_MERGE_MANIFEST_SCHEMA_ID!r}"
            )
        if self.claim_level != MODERN_SCAN_CLAIM_LEVEL:
            raise ValueError(f"claim_level must be exactly {MODERN_SCAN_CLAIM_LEVEL!r}")
        if self.lane_id != MODERN_LANE_ID:
            raise ValueError(f"lane_id must be exactly {MODERN_LANE_ID!r}")
        if not isinstance(self.config_hash, str) or not self.config_hash:
            raise ValueError("config_hash must be a non-empty string")
        if not isinstance(self.input_bundle_id, str) or not self.input_bundle_id:
            raise ValueError("input_bundle_id must be a non-empty string")
        if not isinstance(self.input_provenance_id, str) or not self.input_provenance_id:
            raise ValueError("input_provenance_id must be a non-empty string")
        object.__setattr__(self, "shard_count", _require_positive_int("shard_count", self.shard_count))
        object.__setattr__(
            self,
            "shard_manifest_count",
            _require_nonnegative_int("shard_manifest_count", self.shard_manifest_count),
        )
        object.__setattr__(
            self,
            "total_point_count",
            _require_nonnegative_int("total_point_count", self.total_point_count),
        )
        object.__setattr__(
            self,
            "accepted_point_count",
            _require_nonnegative_int("accepted_point_count", self.accepted_point_count),
        )
        object.__setattr__(
            self,
            "verifier_failed_point_count",
            _require_nonnegative_int(
                "verifier_failed_point_count",
                self.verifier_failed_point_count,
            ),
        )
        object.__setattr__(
            self,
            "bridge_verifier_failed_point_count",
            _require_nonnegative_int(
                "bridge_verifier_failed_point_count",
                self.bridge_verifier_failed_point_count,
            ),
        )
        object.__setattr__(
            self,
            "phenomenology_verifier_failed_point_count",
            _require_nonnegative_int(
                "phenomenology_verifier_failed_point_count",
                self.phenomenology_verifier_failed_point_count,
            ),
        )
        if not isinstance(self.complete, bool):
            raise ValueError("complete must be a bool")
        if not isinstance(self.merged_results_path, str) or not self.merged_results_path:
            raise ValueError("merged_results_path must be a non-empty string")
        if self.accepted_point_count > self.total_point_count:
            raise ValueError("accepted_point_count cannot exceed total_point_count")
        if self.verifier_failed_point_count > self.total_point_count:
            raise ValueError("verifier_failed_point_count cannot exceed total_point_count")
        if self.bridge_verifier_failed_point_count > self.total_point_count:
            raise ValueError("bridge_verifier_failed_point_count cannot exceed total_point_count")
        if self.phenomenology_verifier_failed_point_count > self.total_point_count:
            raise ValueError(
                "phenomenology_verifier_failed_point_count cannot exceed total_point_count"
            )
        object.__setattr__(
            self,
            "source_manifest_paths",
            tuple(str(path) for path in self.source_manifest_paths),
        )
        if not isinstance(self.notes, str) or not self.notes.strip():
            raise ValueError("notes must be a non-empty string")
        object.__setattr__(self, "notes", self.notes.strip())

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "claim_level": self.claim_level,
            "lane_id": self.lane_id,
            "config_hash": self.config_hash,
            "input_bundle_id": self.input_bundle_id,
            "input_provenance_id": self.input_provenance_id,
            "shard_count": self.shard_count,
            "shard_manifest_count": self.shard_manifest_count,
            "total_point_count": self.total_point_count,
            "accepted_point_count": self.accepted_point_count,
            "verifier_failed_point_count": self.verifier_failed_point_count,
            "bridge_verifier_failed_point_count": self.bridge_verifier_failed_point_count,
            "phenomenology_verifier_failed_point_count": self.phenomenology_verifier_failed_point_count,
            "complete": self.complete,
            "merged_results_path": self.merged_results_path,
            "source_manifest_paths": list(self.source_manifest_paths),
            "notes": self.notes,
        }


@dataclass(frozen=True)
class _ArtifactSource:
    schema_id: str
    point_id: str
    point_label: str
    policy: Any
    M_KK: float
    xi_KK: float
    verdicts: tuple[Any, ...]


@dataclass(frozen=True)
class VerifierSubprocessResult:
    ok: bool
    issue_codes: tuple[str, ...]
    message: str = ""


def _point_id_for_values(
    *,
    config_hash: str,
    r: float,
    overall_scale: float,
    Lambda_IR: float,
    point_index: int,
) -> str:
    digest = _stable_hash(
        {
            "config_hash": config_hash,
            "r": float(r),
            "overall_scale": float(overall_scale),
            "Lambda_IR": float(Lambda_IR),
        }
    )[:16]
    return f"mqs-{point_index:08d}-{digest}"


def _point_label(point_id: str, *, r: float, overall_scale: float, Lambda_IR: float) -> str:
    return (
        f"{point_id}:r={_format_float(r)}:"
        f"overall_scale={_format_float(overall_scale)}:"
        f"Lambda_IR={_format_float(Lambda_IR)}"
    )


def enumerate_modern_scan_points(config: ModernScanConfig) -> tuple[ModernScanPoint, ...]:
    """Enumerate the deterministic modern scan grid in canonical order."""

    points: list[ModernScanPoint] = []
    point_index = 0
    for Lambda_IR in config.Lambda_IR_values:
        M_KK = default_quark_m_kk_from_lambda_ir(float(Lambda_IR), xi_KK=config.xi_KK)
        for overall_scale in config.overall_scale_values:
            for r_value in config.r_values:
                point_id = _point_id_for_values(
                    config_hash=config.config_hash,
                    r=float(r_value),
                    overall_scale=float(overall_scale),
                    Lambda_IR=float(Lambda_IR),
                    point_index=point_index,
                )
                points.append(
                    ModernScanPoint(
                        config_hash=config.config_hash,
                        point_index=point_index,
                        point_id=point_id,
                        point_label=_point_label(
                            point_id,
                            r=float(r_value),
                            overall_scale=float(overall_scale),
                            Lambda_IR=float(Lambda_IR),
                        ),
                        r=float(r_value),
                        overall_scale=float(overall_scale),
                        Lambda_IR=float(Lambda_IR),
                        M_KK=M_KK,
                        xi_KK=config.xi_KK,
                        k=config.k,
                    )
                )
                point_index += 1
    return tuple(points)


def point_in_shard(point: ModernScanPoint, *, shard_index: int, shard_count: int) -> bool:
    """Return whether a deterministic point belongs to a given shard."""

    _require_nonnegative_int("shard_index", shard_index)
    _require_positive_int("shard_count", shard_count)
    if shard_index >= shard_count:
        raise ValueError("shard_index must be strictly smaller than shard_count")
    return point.point_index % shard_count == shard_index


def enumerate_modern_scan_shard_points(
    config: ModernScanConfig,
    *,
    shard_index: int = 0,
    shard_count: int = 1,
    max_points: int | None = None,
) -> tuple[ModernScanPoint, ...]:
    """Return the canonical point list assigned to one shard."""

    _require_nonnegative_int("shard_index", shard_index)
    _require_positive_int("shard_count", shard_count)
    if shard_index >= shard_count:
        raise ValueError("shard_index must be strictly smaller than shard_count")
    if max_points is not None:
        _require_positive_int("max_points", max_points)
    points = tuple(
        point
        for point in enumerate_modern_scan_points(config)
        if point_in_shard(point, shard_index=shard_index, shard_count=shard_count)
    )
    if max_points is not None:
        return points[:max_points]
    return points


def _write_run_config_snapshot(run_dir: Path, config: ModernScanConfig) -> Path:
    run_config_path = run_dir / "config.json"
    payload = {
        "schema_id": MODERN_SCAN_RUN_CONFIG_SCHEMA_ID,
        "claim_level": MODERN_SCAN_CLAIM_LEVEL,
        "config_hash": config.config_hash,
        "config": config.as_dict(),
        "record_git_metadata": config.record_git_metadata,
        "notes": MODERN_SCAN_NOTES,
    }
    if run_config_path.exists():
        # Retry reads to handle concurrent shard writes (atomic rename
        # ensures we never see a partial file, but a brief window exists
        # where the file doesn't exist between unlink and rename).
        for _attempt in range(5):
            try:
                existing = _read_json(run_config_path)
                break
            except (json.JSONDecodeError, FileNotFoundError):
                import time
                time.sleep(0.2)
        else:
            existing = _read_json(run_config_path)
        if existing.get("schema_id") == MODERN_SCAN_CONFIG_SCHEMA_ID:
            existing_config = ModernScanConfig.from_dict(existing)
            if existing_config.config_hash != config.config_hash:
                raise ValueError(
                    f"run config at {run_config_path} has hash {existing_config.config_hash!r}, "
                    f"expected {config.config_hash!r}"
                )
            _atomic_write_json(run_config_path, payload)
            return run_config_path
        existing_hash = existing.get("config_hash")
        if existing_hash != config.config_hash:
            raise ValueError(
                f"run config at {run_config_path} has hash {existing_hash!r}, "
                f"expected {config.config_hash!r}"
            )
        return run_config_path
    _atomic_write_json(run_config_path, payload)
    return run_config_path


def _load_existing_result_records(
    *,
    run_dir: Path,
    results_path: Path,
    expected_config_hash: str,
) -> dict[str, ModernScanResultRecord]:
    records: dict[str, ModernScanResultRecord] = {}
    for row in _parse_jsonl(results_path):
        record = ModernScanResultRecord.from_dict(row)
        if record.config_hash != expected_config_hash:
            raise ValueError(
                f"existing record {record.point_id!r} in {results_path} carries "
                f"config_hash {record.config_hash!r}, expected {expected_config_hash!r}"
            )
        artifact_path = run_dir / record.artifact_path
        if not artifact_path.exists():
            raise FileNotFoundError(
                f"resume cannot skip {record.point_id!r}: missing artifact {artifact_path}"
            )
        bridge_artifact_path = run_dir / record.bridge_artifact_path
        if not bridge_artifact_path.exists():
            raise FileNotFoundError(
                "resume cannot skip "
                f"{record.point_id!r}: missing bridge artifact {bridge_artifact_path}"
            )
        phenomenology_artifact_path = run_dir / record.phenomenology_artifact_path
        if not phenomenology_artifact_path.exists():
            raise FileNotFoundError(
                "resume cannot skip "
                f"{record.point_id!r}: missing phenomenology artifact "
                f"{phenomenology_artifact_path}"
            )
        if record.point_id in records:
            raise ValueError(f"duplicate point_id {record.point_id!r} found in {results_path}")
        records[record.point_id] = record
    return records


def _seed_cache_path(cache_dir: Path, point_id: str) -> Path:
    return cache_dir / f"{point_id}.seed.json"


def _read_seed_cache(cache_dir: Path, point_id: str) -> QuarkFitSeed | None:
    path = _seed_cache_path(cache_dir, point_id)
    if not path.exists():
        return None
    return seed_from_dict(_read_json(path))


def _write_seed_cache(cache_dir: Path, point_id: str, seed: QuarkFitSeed) -> Path:
    return _write_json(_seed_cache_path(cache_dir, point_id), seed_to_dict(seed))


def _verifier_subprocess(artifact_path: Path, *, verifier_name: str) -> VerifierSubprocessResult:
    script = """
import json
import sys
from pathlib import Path

from quarkConstraints.modern.verifier import (
    verify_artifact_path,
    verify_bridge_artifact_path,
    verify_phenomenology_artifact_path,
)

artifact_path = Path(sys.argv[1])
verifier_name = sys.argv[2]
verifier = {
    "point": verify_artifact_path,
    "bridge": verify_bridge_artifact_path,
    "phenomenology": verify_phenomenology_artifact_path,
}[verifier_name]
try:
    report = verifier(artifact_path)
except Exception as exc:
    print(json.dumps({
        "ok": False,
        "issue_codes": [type(exc).__name__],
        "message": str(exc),
    }))
else:
    print(json.dumps({
        "ok": report.ok,
        "issue_codes": list(report.issue_codes),
        "message": "",
    }))
"""
    completed = subprocess.run(
        [sys.executable, "-c", script, str(artifact_path), verifier_name],
        cwd=str(_REPO_ROOT),
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        return VerifierSubprocessResult(
            ok=False,
            issue_codes=("verifier_subprocess_failed",),
            message=(completed.stderr or completed.stdout).strip(),
        )
    payload = json.loads(completed.stdout)
    return VerifierSubprocessResult(
        ok=bool(payload["ok"]),
        issue_codes=tuple(str(code) for code in payload["issue_codes"]),
        message=str(payload.get("message", "")),
    )


def _result_record_from_run(
    *,
    point: ModernScanPoint,
    shard_index: int,
    shard_count: int,
    git_commit: str,
    dirty_tree: bool | None,
    solution: QuarkFitSolution,
    evaluation: Any,
    artifact_path: Path,
    bridge_artifact_path: Path,
    phenomenology_artifact: Any,
    phenomenology_artifact_path: Path,
    run_dir: Path,
    verifier_result: VerifierSubprocessResult,
    bridge_verifier_result: VerifierSubprocessResult,
    phenomenology_verifier_result: VerifierSubprocessResult,
) -> ModernScanResultRecord:
    ratio_to_bound_by_system = {
        system_result.system_id: float(system_result.ratio_to_bound)
        for system_result in phenomenology_artifact.system_results
        if system_result.ratio_to_bound is not None
    }
    non_cp_ratio_to_bound_by_system = {
        system_result.system_id: float(system_result.ratio_to_bound)
        for system_result in phenomenology_artifact.system_results
        if (
            system_result.included_in_non_cp_acceptance
            and system_result.ratio_to_bound is not None
        )
    }
    diagnostic_ratio_to_bound_by_system = {
        system_result.system_id: float(system_result.ratio_to_bound)
        for system_result in phenomenology_artifact.system_results
        if (
            system_result.evaluated_from_bridge
            and not system_result.included_in_non_cp_acceptance
            and system_result.ratio_to_bound is not None
        )
    }
    # Derive failing lists from the ratios directly (not from the
    # phenomenology sidecar's pass/fail booleans) to guarantee
    # self-consistency in the scan row.
    failing_non_cp_from_ratios = tuple(
        sid for sid in ("epsilon_K", "K", "B_d", "B_s", "D0")
        if sid in non_cp_ratio_to_bound_by_system
        and non_cp_ratio_to_bound_by_system[sid] > 1.0
    )
    failing_diag_from_ratios = tuple(
        sid for sid in ("epsilon_K", "K", "B_d", "B_s", "D0")
        if sid in diagnostic_ratio_to_bound_by_system
        and diagnostic_ratio_to_bound_by_system[sid] > 1.0
    )
    failing_system_ids = tuple(
        sid for sid in ("epsilon_K", "K", "B_d", "B_s", "D0")
        if sid in failing_non_cp_from_ratios or sid in failing_diag_from_ratios
    )
    diagnostic_only_system_ids = tuple(
        system_result.system_id
        for system_result in phenomenology_artifact.system_results
        if (
            system_result.evaluated_from_bridge
            and not system_result.included_in_non_cp_acceptance
        )
    )
    diagnostic_failing_system_ids = failing_diag_from_ratios
    blocked_system_ids = tuple(
        system_result.system_id
        for system_result in phenomenology_artifact.system_results
        if not system_result.evaluated_from_bridge
    )
    accepted = bool(
        solution.success
        and len(failing_non_cp_from_ratios) == 0
        and verifier_result.ok
        and bridge_verifier_result.ok
        and phenomenology_verifier_result.ok
    )
    return ModernScanResultRecord(
        config_hash=point.config_hash,
        input_bundle_id=MODERN_DEFAULT_INPUT_BUNDLE_ID,
        input_provenance_id=MODERN_DEFAULT_INPUT_PROVENANCE_ID,
        git_commit=git_commit,
        dirty_tree=dirty_tree,
        point_id=point.point_id,
        point_index=point.point_index,
        point_label=point.point_label,
        shard_index=shard_index,
        shard_count=shard_count,
        r=point.r,
        overall_scale=point.overall_scale,
        Lambda_IR=point.Lambda_IR,
        M_KK=point.M_KK,
        xi_KK=point.xi_KK,
        k=point.k,
        fit_success=bool(solution.success),
        fit_message=str(solution.message),
        fit_nfev=int(solution.nfev),
        initial_score=float(solution.initial_score),
        fit_score=float(solution.result.score),
        residual_norm=float(solution.result.residual_norm),
        phenomenology_passes=bool(phenomenology_artifact.non_cp_passes),
        verifier_ok=verifier_result.ok,
        bridge_verifier_ok=bridge_verifier_result.ok,
        phenomenology_verifier_ok=phenomenology_verifier_result.ok,
        accepted=accepted,
        phenomenology_release_scope_id=phenomenology_artifact.release_scope_id,
        non_cp_acceptance_system_ids=phenomenology_artifact.non_cp_acceptance_system_ids,
        failing_non_cp_system_ids=failing_non_cp_from_ratios,
        non_cp_ratio_to_bound_by_system=non_cp_ratio_to_bound_by_system,
        max_non_cp_ratio_to_bound=max(
            non_cp_ratio_to_bound_by_system.values(),
            default=0.0,
        ),
        diagnostic_only_system_ids=diagnostic_only_system_ids,
        diagnostic_failing_system_ids=diagnostic_failing_system_ids,
        diagnostic_ratio_to_bound_by_system=diagnostic_ratio_to_bound_by_system,
        blocked_system_ids=blocked_system_ids,
        max_ratio_to_bound=max(ratio_to_bound_by_system.values(), default=0.0),
        failing_system_ids=failing_system_ids,
        ratio_to_bound_by_system=ratio_to_bound_by_system,
        artifact_path=str(artifact_path.relative_to(run_dir)),
        bridge_artifact_path=str(bridge_artifact_path.relative_to(run_dir)),
        phenomenology_artifact_path=str(
            phenomenology_artifact_path.relative_to(run_dir)
        ),
        verifier_issue_codes=verifier_result.issue_codes,
        bridge_verifier_issue_codes=bridge_verifier_result.issue_codes,
        phenomenology_verifier_issue_codes=phenomenology_verifier_result.issue_codes,
        notes=MODERN_SCAN_NOTES,
    )


def _shard_manifest_from_records(
    *,
    config_hash: str,
    shard_index: int,
    shard_count: int,
    assigned_point_count: int,
    skipped_existing_point_count: int,
    records: Sequence[ModernScanResultRecord],
    run_dir: Path,
    results_path: Path,
    artifacts_dir: Path,
    bridge_artifacts_dir: Path,
    phenomenology_artifacts_dir: Path,
    cache_dir: Path,
    run_config_path: Path,
    git_commit: str,
    dirty_tree: bool | None,
) -> ModernScanShardManifest:
    completed_point_count = len(records)
    accepted_point_count = sum(1 for record in records if record.accepted)
    verifier_failed_point_count = sum(1 for record in records if not record.verifier_ok)
    bridge_verifier_failed_point_count = sum(
        1 for record in records if not record.bridge_verifier_ok
    )
    phenomenology_verifier_failed_point_count = sum(
        1 for record in records if not record.phenomenology_verifier_ok
    )
    return ModernScanShardManifest(
        config_hash=config_hash,
        input_bundle_id=MODERN_DEFAULT_INPUT_BUNDLE_ID,
        input_provenance_id=MODERN_DEFAULT_INPUT_PROVENANCE_ID,
        shard_index=shard_index,
        shard_count=shard_count,
        assigned_point_count=assigned_point_count,
        completed_point_count=completed_point_count,
        skipped_existing_point_count=skipped_existing_point_count,
        accepted_point_count=accepted_point_count,
        verifier_failed_point_count=verifier_failed_point_count,
        bridge_verifier_failed_point_count=bridge_verifier_failed_point_count,
        phenomenology_verifier_failed_point_count=phenomenology_verifier_failed_point_count,
        complete=completed_point_count == assigned_point_count,
        results_path=_display_path(results_path, root=run_dir),
        artifact_dir=_display_path(artifacts_dir, root=run_dir),
        bridge_artifact_dir=_display_path(bridge_artifacts_dir, root=run_dir),
        phenomenology_artifact_dir=_display_path(
            phenomenology_artifacts_dir,
            root=run_dir,
        ),
        cache_dir=_display_path(cache_dir, root=run_dir),
        run_config_path=_display_path(run_config_path, root=run_dir),
        git_commit=git_commit,
        dirty_tree=dirty_tree,
        notes=MODERN_SCAN_NOTES,
    )


def run_modern_scan_shard(
    config: ModernScanConfig,
    *,
    output_dir: str | Path,
    shard_index: int = 0,
    shard_count: int = 1,
    resume: bool = True,
    max_points: int | None = None,
    progress_every: int = 0,
) -> ModernScanShardManifest:
    """Run one deterministic modern scan shard and write resumable outputs."""

    _require_nonnegative_int("shard_index", shard_index)
    _require_positive_int("shard_count", shard_count)
    if shard_index >= shard_count:
        raise ValueError("shard_index must be strictly smaller than shard_count")
    if progress_every < 0:
        raise ValueError("progress_every must be non-negative")
    if max_points is not None:
        _require_positive_int("max_points", max_points)

    run_dir = Path(output_dir)
    shard_dir = run_dir / "shards" / f"shard-{shard_index:05d}-of-{shard_count:05d}"
    artifacts_dir = shard_dir / "artifacts"
    bridge_artifacts_dir = shard_dir / "bridge_artifacts"
    phenomenology_artifacts_dir = shard_dir / "phenomenology_artifacts"
    cache_dir = shard_dir / "cache"
    results_path = shard_dir / "results.jsonl"
    manifest_path = shard_dir / "manifest.json"
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    bridge_artifacts_dir.mkdir(parents=True, exist_ok=True)
    phenomenology_artifacts_dir.mkdir(parents=True, exist_ok=True)
    run_config_path = _write_run_config_snapshot(run_dir, config)
    git_commit, dirty_tree = _resolve_git_metadata(config.record_git_metadata)

    if not resume and results_path.exists():
        raise FileExistsError(f"{results_path} already exists; rerun with resume=True or a new output_dir")

    points = enumerate_modern_scan_shard_points(
        config,
        shard_index=shard_index,
        shard_count=shard_count,
        max_points=max_points,
    )
    existing_records = _load_existing_result_records(
        run_dir=run_dir,
        results_path=results_path,
        expected_config_hash=config.config_hash,
    )
    skipped_existing_point_count = 0

    records = dict(existing_records)
    initial_manifest = _shard_manifest_from_records(
        config_hash=config.config_hash,
        shard_index=shard_index,
        shard_count=shard_count,
        assigned_point_count=len(points),
        skipped_existing_point_count=0,
        records=tuple(sorted(records.values(), key=lambda record: record.point_index)),
        run_dir=run_dir,
        results_path=results_path,
        artifacts_dir=artifacts_dir,
        bridge_artifacts_dir=bridge_artifacts_dir,
        phenomenology_artifacts_dir=phenomenology_artifacts_dir,
        cache_dir=cache_dir,
        run_config_path=run_config_path,
        git_commit=git_commit,
        dirty_tree=dirty_tree,
    )
    _write_json(manifest_path, initial_manifest.as_dict())

    for point in points:
        if point.point_id in records:
            skipped_existing_point_count += 1
            continue
        cached_seed = _read_seed_cache(cache_dir, point.point_id)
        solution = fit_quark_sector(
            config.targets,
            r=point.r,
            overall_scale=point.overall_scale,
            seed=default_spurion_seed() if cached_seed is None else cached_seed,
            k=point.k,
            Lambda_IR=point.Lambda_IR,
            max_nfev=config.max_nfev,
            fit_orientation=config.fit_orientation,
        )
        evaluation = evaluate_modern_point(
            solution.result,
            M_KK=point.M_KK,
            xi_KK=point.xi_KK,
            point_id=point.point_id,
            point_label=point.point_label,
            g_s_star=config.g_s_star,
        )
        artifact = build_modern_point_artifact(
            _ArtifactSource(
                schema_id=evaluation.schema_id,
                point_id=point.point_id,
                point_label=point.point_label,
                policy=evaluation.policy,
                M_KK=evaluation.M_KK,
                xi_KK=evaluation.xi_KK,
                verdicts=evaluation.verdicts,
            )
        )
        bridge_artifact = build_modern_point_bridge_artifact(evaluation)
        phenomenology_artifact = build_modern_point_phenomenology_artifact(
            bridge_artifact,
            policy=evaluation.policy,
        )
        artifact_path = artifacts_dir / f"{point.point_id}.json"
        bridge_artifact_path = bridge_artifacts_dir / f"{point.point_id}.bridge.json"
        phenomenology_artifact_path = (
            phenomenology_artifacts_dir / f"{point.point_id}.phenomenology.json"
        )
        write_modern_point_artifact(artifact, artifact_path)
        write_modern_point_bridge_artifact(bridge_artifact, bridge_artifact_path)
        write_modern_point_phenomenology_artifact(
            phenomenology_artifact,
            phenomenology_artifact_path,
        )
        verifier_result = _verifier_subprocess(artifact_path, verifier_name="point")
        bridge_verifier_result = _verifier_subprocess(
            bridge_artifact_path,
            verifier_name="bridge",
        )
        phenomenology_verifier_result = _verifier_subprocess(
            phenomenology_artifact_path,
            verifier_name="phenomenology",
        )
        record = _result_record_from_run(
            point=point,
            shard_index=shard_index,
            shard_count=shard_count,
            git_commit=git_commit,
            dirty_tree=dirty_tree,
            solution=solution,
            evaluation=evaluation,
            artifact_path=artifact_path,
            bridge_artifact_path=bridge_artifact_path,
            phenomenology_artifact=phenomenology_artifact,
            phenomenology_artifact_path=phenomenology_artifact_path,
            run_dir=run_dir,
            verifier_result=verifier_result,
            bridge_verifier_result=bridge_verifier_result,
            phenomenology_verifier_result=phenomenology_verifier_result,
        )
        records[record.point_id] = record
        _append_jsonl(results_path, record.as_dict())
        _write_seed_cache(cache_dir, point.point_id, solution.seed)

        current_manifest = _shard_manifest_from_records(
            config_hash=config.config_hash,
            shard_index=shard_index,
            shard_count=shard_count,
            assigned_point_count=len(points),
            skipped_existing_point_count=skipped_existing_point_count,
            records=tuple(sorted(records.values(), key=lambda item: item.point_index)),
            run_dir=run_dir,
            results_path=results_path,
            artifacts_dir=artifacts_dir,
            bridge_artifacts_dir=bridge_artifacts_dir,
            phenomenology_artifacts_dir=phenomenology_artifacts_dir,
            cache_dir=cache_dir,
            run_config_path=run_config_path,
            git_commit=git_commit,
            dirty_tree=dirty_tree,
        )
        _write_json(manifest_path, current_manifest.as_dict())
        if progress_every and len(records) % progress_every == 0:
            print(
                f"[modern.scan] shard {shard_index}/{shard_count} "
                f"completed {len(records)}/{len(points)} assigned points",
                file=sys.stderr,
            )

    final_manifest = _shard_manifest_from_records(
        config_hash=config.config_hash,
        shard_index=shard_index,
        shard_count=shard_count,
        assigned_point_count=len(points),
        skipped_existing_point_count=skipped_existing_point_count,
        records=tuple(sorted(records.values(), key=lambda record: record.point_index)),
        run_dir=run_dir,
        results_path=results_path,
        artifacts_dir=artifacts_dir,
        bridge_artifacts_dir=bridge_artifacts_dir,
        phenomenology_artifacts_dir=phenomenology_artifacts_dir,
        cache_dir=cache_dir,
        run_config_path=run_config_path,
        git_commit=git_commit,
        dirty_tree=dirty_tree,
    )
    _write_json(manifest_path, final_manifest.as_dict())
    return final_manifest


def merge_modern_scan_shards(
    run_dir: str | Path,
    *,
    output_dir: str | Path | None = None,
    require_complete: bool = True,
) -> ModernMergedScanManifest:
    """Merge deterministic shard outputs into one reduction file."""

    root = Path(run_dir)
    merged_dir = root / "merged" if output_dir is None else Path(output_dir)
    manifest_paths = sorted(root.glob("shards/shard-*-of-*/manifest.json"))
    if not manifest_paths:
        raise FileNotFoundError(f"no shard manifests found under {root / 'shards'}")

    shard_manifests = tuple(
        ModernScanShardManifest.from_dict(_read_json(path)) for path in manifest_paths
    )
    config_hashes = {manifest.config_hash for manifest in shard_manifests}
    if len(config_hashes) != 1:
        raise ValueError(f"shard manifests disagree on config_hash: {sorted(config_hashes)!r}")
    shard_counts = {manifest.shard_count for manifest in shard_manifests}
    if len(shard_counts) != 1:
        raise ValueError(f"shard manifests disagree on shard_count: {sorted(shard_counts)!r}")
    shard_count = shard_manifests[0].shard_count
    shard_indexes = {manifest.shard_index for manifest in shard_manifests}
    if require_complete:
        expected_indexes = set(range(shard_count))
        missing = sorted(expected_indexes - shard_indexes)
        if missing:
            raise ValueError(f"missing shard manifests for shard indexes {missing!r}")
        incomplete = sorted(
            manifest.shard_index for manifest in shard_manifests if not manifest.complete
        )
        if incomplete:
            raise ValueError(f"incomplete shard manifests present for shard indexes {incomplete!r}")

    all_records: list[ModernScanResultRecord] = []
    seen_point_ids: set[str] = set()
    for manifest in sorted(shard_manifests, key=lambda item: item.shard_index):
        results_path = root / manifest.results_path
        for row in _parse_jsonl(results_path):
            record = ModernScanResultRecord.from_dict(row)
            if record.point_id in seen_point_ids:
                raise ValueError(f"duplicate point_id {record.point_id!r} encountered while merging")
            seen_point_ids.add(record.point_id)
            all_records.append(record)

    all_records.sort(key=lambda record: record.point_index)
    merged_results_path = merged_dir / "results.jsonl"
    merged_dir.mkdir(parents=True, exist_ok=True)
    with merged_results_path.open("w", encoding="utf-8") as handle:
        for record in all_records:
            handle.write(json.dumps(record.as_dict(), sort_keys=True, allow_nan=False) + "\n")

    merged_manifest = ModernMergedScanManifest(
        config_hash=next(iter(config_hashes)),
        shard_count=shard_count,
        shard_manifest_count=len(shard_manifests),
        total_point_count=len(all_records),
        accepted_point_count=sum(1 for record in all_records if record.accepted),
        verifier_failed_point_count=sum(1 for record in all_records if not record.verifier_ok),
        bridge_verifier_failed_point_count=sum(
            1 for record in all_records if not record.bridge_verifier_ok
        ),
        phenomenology_verifier_failed_point_count=sum(
            1 for record in all_records if not record.phenomenology_verifier_ok
        ),
        complete=len(shard_manifests) == shard_count and all(
            manifest.complete for manifest in shard_manifests
        ),
        merged_results_path=_display_path(merged_results_path, root=root),
        source_manifest_paths=tuple(str(path.relative_to(root)) for path in manifest_paths),
        notes=MODERN_SCAN_NOTES,
    )
    _write_json(merged_dir / "manifest.json", merged_manifest.as_dict())
    return merged_manifest


def verify_merged_scan(
    config: ModernScanConfig,
    *,
    output_root: str | Path,
    total_shards: int,
) -> dict[str, Any]:
    """Verify merged scan outputs without importing the standalone verifier."""

    root = Path(output_root)
    merged_manifest_path = root / "merged" / "manifest.json"
    merged_results_path = root / "merged" / "results.jsonl"
    issues: list[str] = []
    if not merged_manifest_path.exists():
        issues.append("missing_merged_manifest")
    if not merged_results_path.exists():
        issues.append("missing_merged_results")
    rows = _parse_jsonl(merged_results_path) if merged_results_path.exists() else []
    records = [ModernScanResultRecord.from_dict(row) for row in rows]
    if len(records) != len(enumerate_modern_scan_points(config)):
        issues.append("point_count_mismatch")
    if len({record.point_id for record in records}) != len(records):
        issues.append("duplicate_point_ids")
    if [record.point_index for record in records] != sorted(record.point_index for record in records):
        issues.append("merged_records_not_sorted")
    for shard_index in range(total_shards):
        if not (root / "shards" / f"shard-{shard_index:05d}-of-{total_shards:05d}" / "manifest.json").exists():
            issues.append("missing_shard_manifest")
            break
    for record in records:
        artifact_path = root / record.artifact_path
        bridge_artifact_path = root / record.bridge_artifact_path
        phenomenology_artifact_path = root / record.phenomenology_artifact_path
        if not artifact_path.exists():
            issues.append("missing_artifact")
            break
        if not bridge_artifact_path.exists():
            issues.append("missing_bridge_artifact")
            break
        if not phenomenology_artifact_path.exists():
            issues.append("missing_phenomenology_artifact")
            break
        if not record.verifier_ok:
            issues.append("artifact_verifier_failed")
            break
        if not record.bridge_verifier_ok:
            issues.append("bridge_artifact_verifier_failed")
            break
        if not record.phenomenology_verifier_ok:
            issues.append("phenomenology_artifact_verifier_failed")
            break
    payload = {
        "schema_id": MODERN_SCAN_VERIFICATION_SCHEMA_ID,
        "lane_id": MODERN_LANE_ID,
        "claim_level": MODERN_SCAN_CLAIM_LEVEL,
        "config_hash": config.config_hash,
        "input_bundle_id": MODERN_DEFAULT_INPUT_BUNDLE_ID,
        "input_provenance_id": MODERN_DEFAULT_INPUT_PROVENANCE_ID,
        "expected_point_count": len(enumerate_modern_scan_points(config)),
        "observed_point_count": len(records),
        "accepted_point_count": sum(1 for record in records if record.accepted),
        "verifier_failed_point_count": sum(1 for record in records if not record.verifier_ok),
        "bridge_verifier_failed_point_count": sum(
            1 for record in records if not record.bridge_verifier_ok
        ),
        "phenomenology_verifier_failed_point_count": sum(
            1 for record in records if not record.phenomenology_verifier_ok
        ),
        "total_shards": int(total_shards),
        "ok": len(issues) == 0,
        "issues": issues,
        "notes": MODERN_SCAN_NOTES,
    }
    _write_json(root / "merged" / "verification.json", payload)
    return payload


def smoke_scan_config() -> ModernScanConfig:
    return ModernScanConfig(
        r_values=np.array([0.25], dtype=float),
        overall_scale_values=np.array([2.8], dtype=float),
        Lambda_IR_values=np.array([3000.0], dtype=float),
        max_nfev=6,
    )


def pilot_scan_config() -> ModernScanConfig:
    return ModernScanConfig(
        r_values=np.array([0.10, 0.25], dtype=float),
        overall_scale_values=np.array([2.8, 3.0], dtype=float),
        Lambda_IR_values=np.array([3000.0, 4500.0], dtype=float),
        max_nfev=20,
    )


def exploratory_scan_config() -> ModernScanConfig:
    return ModernScanConfig(
        r_values=np.logspace(np.log10(0.05), np.log10(0.8), 8),
        overall_scale_values=np.array([3.0], dtype=float),
        Lambda_IR_values=np.logspace(np.log10(1500.0), np.log10(10000.0), 8),
        max_nfev=80,
    )


def production_scan_config() -> ModernScanConfig:
    return ModernScanConfig(
        r_values=np.logspace(np.log10(0.03), np.log10(1.0), 20),
        overall_scale_values=np.array([2.0, 3.0, 4.0, 5.0], dtype=float),
        Lambda_IR_values=np.logspace(np.log10(1000.0), np.log10(15000.0), 15),
        max_nfev=120,
    )


def dense_scan_config() -> ModernScanConfig:
    """Dense scan preset: 100 r × 20 overall_scale × 50 Lambda_IR = 100,000 points."""
    return ModernScanConfig(
        r_values=np.logspace(np.log10(0.02), np.log10(2.0), 100),
        overall_scale_values=np.linspace(1.5, 6.0, 20),
        Lambda_IR_values=np.logspace(np.log10(500.0), np.log10(20000.0), 50),
        max_nfev=120,
    )


def dense_wide_yukawa_scan_config() -> ModernScanConfig:
    """Wide Yukawa diagnostic: 50 r × 40 overall_scale × 50 Lambda_IR = 100,000 points."""
    return ModernScanConfig(
        r_values=np.logspace(np.log10(0.02), np.log10(2.0), 50),
        overall_scale_values=np.logspace(np.log10(0.001), np.log10(100.0), 40),
        Lambda_IR_values=np.logspace(np.log10(500.0), np.log10(20000.0), 50),
        max_nfev=120,
    )


def write_scan_config(config: ModernScanConfig, path: str | Path) -> Path:
    return _write_json(Path(path), config.as_dict())


def read_scan_config(path: str | Path) -> ModernScanConfig:
    payload = _read_json(Path(path))
    if payload.get("schema_id") == MODERN_SCAN_RUN_CONFIG_SCHEMA_ID:
        payload = dict(payload["config"])
    return ModernScanConfig.from_dict(payload)


def build_point_id(
    config: ModernScanConfig,
    *,
    r: float,
    overall_scale: float,
    Lambda_IR: float,
) -> str:
    for point in enumerate_modern_scan_points(config):
        if (
            point.r == float(r)
            and point.overall_scale == float(overall_scale)
            and point.Lambda_IR == float(Lambda_IR)
        ):
            return point.point_id
    raise ValueError("requested point is not present in config")


def enumerate_scan_points(config: ModernScanConfig) -> tuple[ModernScanPoint, ...]:
    return enumerate_modern_scan_points(config)


def scan_points_for_shard(
    config: ModernScanConfig,
    *,
    shard_id: int,
    total_shards: int,
) -> tuple[ModernScanPoint, ...]:
    return enumerate_modern_scan_shard_points(
        config,
        shard_index=shard_id,
        shard_count=total_shards,
    )


def run_scan_shard(
    config: ModernScanConfig,
    *,
    output_root: str | Path,
    shard_id: int,
    total_shards: int,
    resume: bool = True,
) -> dict[str, Any]:
    manifest = run_modern_scan_shard(
        config,
        output_dir=output_root,
        shard_index=shard_id,
        shard_count=total_shards,
        resume=resume,
    )
    return manifest.as_dict()


def merge_scan_shards(
    config: ModernScanConfig,
    *,
    output_root: str | Path,
    total_shards: int,
) -> dict[str, Any]:
    del config, total_shards
    manifest = merge_modern_scan_shards(output_root)
    return manifest.as_dict()


def _parse_float_csv(values: str) -> list[float]:
    items = [item.strip() for item in values.split(",") if item.strip()]
    if not items:
        raise ValueError("expected at least one float")
    return [float(item) for item in items]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m quarkConstraints.modern.scan",
        description=(
            "Operational modern quark scan machinery. This CLI is shardable and "
            "resumable, but it is not by itself a final allowed-region claim."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    preset_parser = subparsers.add_parser(
        "write-preset",
        help="write a smoke or pilot config snapshot",
    )
    preset_parser.add_argument("preset", choices=("smoke", "pilot", "exploratory", "production", "dense", "dense_wide_yukawa"))
    preset_parser.add_argument("output")

    run_parser = subparsers.add_parser(
        "run-shard",
        help="run one deterministic scan shard and write artifacts/results",
    )
    run_parser.add_argument("--output-dir", required=True)
    run_parser.add_argument("--config")
    run_parser.add_argument("--r-values", default="0.1,0.25,0.4")
    run_parser.add_argument("--overall-scale-values", default="3.0")
    run_parser.add_argument("--lambda-ir-values", default="3000.0")
    run_parser.add_argument("--xi-kk", type=float, default=DEFAULT_QUARK_XI_KK)
    run_parser.add_argument("--k", type=float, default=1.2209e19)
    run_parser.add_argument("--max-nfev", type=int, default=120)
    run_parser.add_argument("--shard-index", type=int, default=0)
    run_parser.add_argument("--shard-count", type=int, default=1)
    run_parser.add_argument("--max-points", type=int)
    run_parser.add_argument("--progress-every", type=int, default=0)
    run_parser.add_argument(
        "--fit-orientation",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    run_parser.add_argument(
        "--record-git-metadata",
        action=argparse.BooleanOptionalAction,
        default=True,
    )

    merge_parser = subparsers.add_parser(
        "merge",
        help="merge shard outputs into one reduction file",
    )
    merge_parser.add_argument("--run-dir", required=True)
    merge_parser.add_argument("--output-dir")
    merge_parser.add_argument(
        "--require-complete",
        action=argparse.BooleanOptionalAction,
        default=True,
    )

    verify_parser = subparsers.add_parser(
        "verify",
        help="verify merged outputs without re-running point evaluation",
    )
    verify_parser.add_argument("--run-dir", required=True)
    verify_parser.add_argument("--config")
    verify_parser.add_argument("--total-shards", type=int, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point for `python -m quarkConstraints.modern.scan`."""

    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command == "write-preset":
        _preset_configs = {
            "smoke": smoke_scan_config,
            "pilot": pilot_scan_config,
            "exploratory": exploratory_scan_config,
            "production": production_scan_config,
            "dense": dense_scan_config,
            "dense_wide_yukawa": dense_wide_yukawa_scan_config,
        }
        config = _preset_configs[args.preset]()
        write_scan_config(config, args.output)
        return 0
    if args.command == "run-shard":
        config = (
            read_scan_config(args.config)
            if args.config
            else ModernScanConfig(
                r_values=_parse_float_csv(args.r_values),
                overall_scale_values=_parse_float_csv(args.overall_scale_values),
                Lambda_IR_values=_parse_float_csv(args.lambda_ir_values),
                xi_KK=args.xi_kk,
                k=args.k,
                max_nfev=args.max_nfev,
                fit_orientation=args.fit_orientation,
                record_git_metadata=args.record_git_metadata,
            )
        )
        manifest = run_modern_scan_shard(
            config,
            output_dir=args.output_dir,
            shard_index=args.shard_index,
            shard_count=args.shard_count,
            max_points=args.max_points,
            progress_every=args.progress_every,
        )
        print(_canonical_json(manifest.as_dict()), end="")
        return 0
    if args.command == "merge":
        manifest = merge_modern_scan_shards(
            args.run_dir,
            output_dir=args.output_dir,
            require_complete=args.require_complete,
        )
        print(_canonical_json(manifest.as_dict()), end="")
        return 0
    if args.command == "verify":
        config = (
            read_scan_config(args.config)
            if args.config
            else read_scan_config(Path(args.run_dir) / "config.json")
        )
        payload = verify_merged_scan(
            config,
            output_root=args.run_dir,
            total_shards=args.total_shards,
        )
        print(_canonical_json(payload), end="")
        return 0
    parser.error(f"unknown command {args.command!r}")
    return 2


__all__ = [
    "MODERN_SCAN_CLAIM_LEVEL",
    "MODERN_SCAN_CLAIM_LEVEL_ID",
    "MODERN_SCAN_CONFIG_SCHEMA_ID",
    "MODERN_SCAN_MERGE_MANIFEST_SCHEMA_ID",
    "MODERN_SCAN_NOTES",
    "MODERN_SCAN_PHYSICS_BACKEND_ID",
    "MODERN_SCAN_POINT_SCHEMA_ID",
    "MODERN_SCAN_ROW_SCHEMA_ID",
    "MODERN_SCAN_RESULT_SCHEMA_ID",
    "MODERN_SCAN_RUN_CONFIG_SCHEMA_ID",
    "MODERN_SCAN_SCHEMA_VERSION",
    "MODERN_SCAN_SHARD_MANIFEST_SCHEMA_ID",
    "MODERN_SCAN_VERIFICATION_SCHEMA_ID",
    "ModernMergedScanManifest",
    "ModernScanConfig",
    "ModernScanPoint",
    "ModernScanResultRecord",
    "ModernScanShardManifest",
    "VerifierSubprocessResult",
    "build_point_id",
    "enumerate_scan_points",
    "enumerate_modern_scan_points",
    "enumerate_modern_scan_shard_points",
    "dense_scan_config",
    "exploratory_scan_config",
    "main",
    "merge_scan_shards",
    "merge_modern_scan_shards",
    "pilot_scan_config",
    "production_scan_config",
    "point_in_shard",
    "read_scan_config",
    "run_modern_scan_shard",
    "run_scan_shard",
    "scan_points_for_shard",
    "seed_from_dict",
    "seed_to_dict",
    "smoke_scan_config",
    "verify_merged_scan",
    "write_scan_config",
]


if __name__ == "__main__":
    raise SystemExit(main())
