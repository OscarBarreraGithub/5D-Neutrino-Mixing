"""Artifact-only sidecar export for the modern QS2 bridge."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .conventions import MODERN_LANE_ID

MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID = "quarkConstraints.modern.artifacts.bridge.v1"
MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION = 1
MODERN_POINT_COUPLINGS_SCHEMA_ID = "quarkConstraints.modern.couplings.point.v1"
MODERN_POINT_MATCHING_SCHEMA_ID = "quarkConstraints.modern.matching.point.v1"
MODERN_POINT_MATCHING_SYSTEM_SCHEMA_ID = "quarkConstraints.modern.matching.system.v1"
MODERN_DEFAULT_INPUTS_SCHEMA_ID = "quarkConstraints.modern.inputs.modern_default_inputs.v1"
MODERN_POINT_EVALUATION_SCHEMA_ID = "quarkConstraints.modern.evaluation.v1"
MODERN_MATCHING_SYSTEM_IDS = ("K", "B_d", "B_s", "D0")

_COUPLINGS_KEYS = (
    "schema_id",
    "lane_id",
    "input_bundle_schema_id",
    "input_bundle_id",
    "input_provenance_id",
    "input_resolution_policy_id",
    "qcd_metadata_id",
    "alpha_s_policy_id",
    "operator_convention_id",
    "ckm_target_id",
    "quark_mass_target_id",
    "point_id",
    "point_label",
    "Lambda_IR",
    "M_KK",
    "xi_KK",
    "alpha_s",
    "g_s",
    "left_overlap",
    "right_up_overlap",
    "right_down_overlap",
    "left_up",
    "left_down",
    "right_up",
    "right_down",
    "notes",
)
_MATCHING_SYSTEM_KEYS = (
    "schema_id",
    "system_id",
    "observable_id",
    "backend_system_id",
    "backend_key",
    "display_name",
    "sector_id",
    "generations",
    "operator_basis_id",
    "matching_scale_GeV",
    "weight_policy_id",
    "left_coupling",
    "right_coupling",
    "c1_vll",
    "c1_vrr",
    "c4_lr",
    "c5_lr",
    "note",
)
_MATCHING_KEYS = (
    "schema_id",
    "lane_id",
    "couplings_schema_id",
    "input_bundle_schema_id",
    "input_bundle_id",
    "input_provenance_id",
    "input_resolution_policy_id",
    "qcd_metadata_id",
    "alpha_s_policy_id",
    "operator_basis_id",
    "weight_policy_id",
    "point_id",
    "point_label",
    "Lambda_IR",
    "M_KK",
    "xi_KK",
    "alpha_s",
    "g_s",
    "system_matches",
    "notes",
)


class ArtifactSchemaError(ValueError):
    """Raised when a modern bridge artifact violates the frozen JSON contract."""


def _require_text(name: str, value: Any) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ArtifactSchemaError(f"{name} must be a non-empty string")
    return value.strip()


def _require_int(name: str, value: Any) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise ArtifactSchemaError(f"{name} must be an integer")
    return value


def _require_positive_float(name: str, value: Any) -> float:
    numeric = float(value)
    if not math.isfinite(numeric) or numeric <= 0.0:
        raise ArtifactSchemaError(f"{name} must be a positive finite float")
    return numeric


def _require_finite_float(name: str, value: Any) -> float:
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ArtifactSchemaError(f"{name} must be finite")
    return numeric


def _require_mapping(value: Any, *, context: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ArtifactSchemaError(f"{context} must be a JSON object")
    return value


def _require_sequence(value: Any, *, context: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise ArtifactSchemaError(f"{context} must be a JSON array")
    return value


def _require_exact(name: str, value: Any, *, expected: str | int) -> str | int:
    if value != expected:
        raise ArtifactSchemaError(f"{name} must be exactly {expected!r}")
    return value


def _require_exact_keys(
    mapping: Mapping[str, Any],
    *,
    expected: tuple[str, ...],
    context: str,
) -> None:
    keys = tuple(mapping.keys())
    if set(keys) != set(expected):
        missing = tuple(key for key in expected if key not in mapping)
        unexpected = tuple(key for key in keys if key not in expected)
        parts: list[str] = []
        if missing:
            parts.append(f"missing keys {missing!r}")
        if unexpected:
            parts.append(f"unexpected keys {unexpected!r}")
        raise ArtifactSchemaError(f"{context} must contain exactly {expected!r} ({'; '.join(parts)})")


def _deep_clone_json(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True, allow_nan=False))


def _canonical_complex_payload(name: str, value: Any) -> dict[str, float]:
    payload = _require_mapping(value, context=name)
    _require_exact_keys(payload, expected=("real", "imag"), context=name)
    return {
        "real": _require_finite_float(f"{name}.real", payload["real"]),
        "imag": _require_finite_float(f"{name}.imag", payload["imag"]),
    }


def _canonical_matrix_payload(name: str, value: Any) -> dict[str, list[list[float]]]:
    payload = _require_mapping(value, context=name)
    _require_exact_keys(payload, expected=("real", "imag"), context=name)
    matrices: dict[str, list[list[float]]] = {}
    for component in ("real", "imag"):
        rows = _require_sequence(payload[component], context=f"{name}.{component}")
        if len(rows) != 3:
            raise ArtifactSchemaError(f"{name}.{component} must have exactly 3 rows")
        canonical_rows: list[list[float]] = []
        for row_index, row in enumerate(rows):
            entries = _require_sequence(row, context=f"{name}.{component}[{row_index}]")
            if len(entries) != 3:
                raise ArtifactSchemaError(
                    f"{name}.{component}[{row_index}] must have exactly 3 columns"
                )
            canonical_rows.append(
                [
                    _require_finite_float(
                        f"{name}.{component}[{row_index}][{column_index}]",
                        entry,
                    )
                    for column_index, entry in enumerate(entries)
                ]
            )
        matrices[component] = canonical_rows
    return matrices


def _observable_id_for_system(system_id: str) -> str:
    if system_id == "K":
        return "epsilon_K"
    if system_id in {"B_d", "B_s", "D0"}:
        return system_id
    raise ArtifactSchemaError(f"system_id must be one of {MODERN_MATCHING_SYSTEM_IDS!r}")


def _backend_system_id_for_system(system_id: str) -> str:
    if system_id in {"K", "B_d", "B_s"}:
        return system_id
    if system_id == "D0":
        return "D"
    raise ArtifactSchemaError(f"system_id must be one of {MODERN_MATCHING_SYSTEM_IDS!r}")


def _canonical_couplings_payload(value: Any) -> dict[str, Any]:
    payload = _require_mapping(value, context="bridge_artifact.couplings")
    _require_exact_keys(payload, expected=_COUPLINGS_KEYS, context="bridge_artifact.couplings")
    return {
        "schema_id": str(
            _require_exact(
                "bridge_artifact.couplings.schema_id",
                _require_text("bridge_artifact.couplings.schema_id", payload["schema_id"]),
                expected=MODERN_POINT_COUPLINGS_SCHEMA_ID,
            )
        ),
        "lane_id": str(
            _require_exact(
                "bridge_artifact.couplings.lane_id",
                _require_text("bridge_artifact.couplings.lane_id", payload["lane_id"]),
                expected=MODERN_LANE_ID,
            )
        ),
        "input_bundle_schema_id": str(
            _require_exact(
                "bridge_artifact.couplings.input_bundle_schema_id",
                _require_text(
                    "bridge_artifact.couplings.input_bundle_schema_id",
                    payload["input_bundle_schema_id"],
                ),
                expected=MODERN_DEFAULT_INPUTS_SCHEMA_ID,
            )
        ),
        "input_bundle_id": _require_text(
            "bridge_artifact.couplings.input_bundle_id", payload["input_bundle_id"]
        ),
        "input_provenance_id": _require_text(
            "bridge_artifact.couplings.input_provenance_id", payload["input_provenance_id"]
        ),
        "input_resolution_policy_id": _require_text(
            "bridge_artifact.couplings.input_resolution_policy_id",
            payload["input_resolution_policy_id"],
        ),
        "qcd_metadata_id": _require_text(
            "bridge_artifact.couplings.qcd_metadata_id", payload["qcd_metadata_id"]
        ),
        "alpha_s_policy_id": _require_text(
            "bridge_artifact.couplings.alpha_s_policy_id", payload["alpha_s_policy_id"]
        ),
        "operator_convention_id": _require_text(
            "bridge_artifact.couplings.operator_convention_id",
            payload["operator_convention_id"],
        ),
        "ckm_target_id": _require_text(
            "bridge_artifact.couplings.ckm_target_id", payload["ckm_target_id"]
        ),
        "quark_mass_target_id": _require_text(
            "bridge_artifact.couplings.quark_mass_target_id",
            payload["quark_mass_target_id"],
        ),
        "point_id": _require_text("bridge_artifact.couplings.point_id", payload["point_id"]),
        "point_label": _require_text(
            "bridge_artifact.couplings.point_label", payload["point_label"]
        ),
        "Lambda_IR": _require_positive_float(
            "bridge_artifact.couplings.Lambda_IR", payload["Lambda_IR"]
        ),
        "M_KK": _require_positive_float("bridge_artifact.couplings.M_KK", payload["M_KK"]),
        "xi_KK": _require_positive_float(
            "bridge_artifact.couplings.xi_KK", payload["xi_KK"]
        ),
        "alpha_s": _require_positive_float(
            "bridge_artifact.couplings.alpha_s", payload["alpha_s"]
        ),
        "g_s": _require_positive_float("bridge_artifact.couplings.g_s", payload["g_s"]),
        "left_overlap": _canonical_matrix_payload(
            "bridge_artifact.couplings.left_overlap", payload["left_overlap"]
        ),
        "right_up_overlap": _canonical_matrix_payload(
            "bridge_artifact.couplings.right_up_overlap", payload["right_up_overlap"]
        ),
        "right_down_overlap": _canonical_matrix_payload(
            "bridge_artifact.couplings.right_down_overlap", payload["right_down_overlap"]
        ),
        "left_up": _canonical_matrix_payload(
            "bridge_artifact.couplings.left_up", payload["left_up"]
        ),
        "left_down": _canonical_matrix_payload(
            "bridge_artifact.couplings.left_down", payload["left_down"]
        ),
        "right_up": _canonical_matrix_payload(
            "bridge_artifact.couplings.right_up", payload["right_up"]
        ),
        "right_down": _canonical_matrix_payload(
            "bridge_artifact.couplings.right_down", payload["right_down"]
        ),
        "notes": _require_text("bridge_artifact.couplings.notes", payload["notes"]),
    }


def _canonical_matching_system_payload(value: Any) -> dict[str, Any]:
    payload = _require_mapping(value, context="bridge_artifact.matching.system")
    _require_exact_keys(
        payload,
        expected=_MATCHING_SYSTEM_KEYS,
        context="bridge_artifact.matching.system",
    )
    system_id = _require_text("bridge_artifact.matching.system.system_id", payload["system_id"])
    if system_id not in MODERN_MATCHING_SYSTEM_IDS:
        raise ArtifactSchemaError(
            f"bridge_artifact.matching.system.system_id must be one of {MODERN_MATCHING_SYSTEM_IDS!r}"
        )
    observable_id = _require_text(
        "bridge_artifact.matching.system.observable_id", payload["observable_id"]
    )
    if observable_id != _observable_id_for_system(system_id):
        raise ArtifactSchemaError(
            f"bridge_artifact.matching.system.observable_id for {system_id!r} must be exactly "
            f"{_observable_id_for_system(system_id)!r}"
        )
    backend_system_id = _require_text(
        "bridge_artifact.matching.system.backend_system_id",
        payload["backend_system_id"],
    )
    if backend_system_id != _backend_system_id_for_system(system_id):
        raise ArtifactSchemaError(
            f"bridge_artifact.matching.system.backend_system_id for {system_id!r} must be exactly "
            f"{_backend_system_id_for_system(system_id)!r}"
        )
    generations = _require_sequence(
        payload["generations"], context="bridge_artifact.matching.system.generations"
    )
    if len(generations) != 2:
        raise ArtifactSchemaError(
            "bridge_artifact.matching.system.generations must have length 2"
        )
    i = _require_int("bridge_artifact.matching.system.generations[0]", generations[0])
    j = _require_int("bridge_artifact.matching.system.generations[1]", generations[1])
    if i == j or i not in {0, 1, 2} or j not in {0, 1, 2}:
        raise ArtifactSchemaError(
            "bridge_artifact.matching.system.generations must be distinct integers in {0, 1, 2}"
        )
    return {
        "schema_id": str(
            _require_exact(
                "bridge_artifact.matching.system.schema_id",
                _require_text("bridge_artifact.matching.system.schema_id", payload["schema_id"]),
                expected=MODERN_POINT_MATCHING_SYSTEM_SCHEMA_ID,
            )
        ),
        "system_id": system_id,
        "observable_id": observable_id,
        "backend_system_id": backend_system_id,
        "backend_key": _require_text(
            "bridge_artifact.matching.system.backend_key", payload["backend_key"]
        ),
        "display_name": _require_text(
            "bridge_artifact.matching.system.display_name", payload["display_name"]
        ),
        "sector_id": _require_text(
            "bridge_artifact.matching.system.sector_id", payload["sector_id"]
        ),
        "generations": [i, j],
        "operator_basis_id": _require_text(
            "bridge_artifact.matching.system.operator_basis_id",
            payload["operator_basis_id"],
        ),
        "matching_scale_GeV": _require_positive_float(
            "bridge_artifact.matching.system.matching_scale_GeV",
            payload["matching_scale_GeV"],
        ),
        "weight_policy_id": _require_text(
            "bridge_artifact.matching.system.weight_policy_id",
            payload["weight_policy_id"],
        ),
        "left_coupling": _canonical_complex_payload(
            "bridge_artifact.matching.system.left_coupling",
            payload["left_coupling"],
        ),
        "right_coupling": _canonical_complex_payload(
            "bridge_artifact.matching.system.right_coupling",
            payload["right_coupling"],
        ),
        "c1_vll": _canonical_complex_payload(
            "bridge_artifact.matching.system.c1_vll", payload["c1_vll"]
        ),
        "c1_vrr": _canonical_complex_payload(
            "bridge_artifact.matching.system.c1_vrr", payload["c1_vrr"]
        ),
        "c4_lr": _canonical_complex_payload(
            "bridge_artifact.matching.system.c4_lr", payload["c4_lr"]
        ),
        "c5_lr": _canonical_complex_payload(
            "bridge_artifact.matching.system.c5_lr", payload["c5_lr"]
        ),
        "note": _require_text("bridge_artifact.matching.system.note", payload["note"]),
    }


def _canonical_matching_payload(value: Any) -> dict[str, Any]:
    payload = _require_mapping(value, context="bridge_artifact.matching")
    _require_exact_keys(payload, expected=_MATCHING_KEYS, context="bridge_artifact.matching")
    system_matches = [
        _canonical_matching_system_payload(item)
        for item in _require_sequence(
            payload["system_matches"],
            context="bridge_artifact.matching.system_matches",
        )
    ]
    if len(system_matches) != len(MODERN_MATCHING_SYSTEM_IDS):
        raise ArtifactSchemaError(
            "bridge_artifact.matching.system_matches must contain the frozen modern system set"
        )
    if tuple(item["system_id"] for item in system_matches) != MODERN_MATCHING_SYSTEM_IDS:
        raise ArtifactSchemaError(
            "bridge_artifact.matching.system_matches must preserve the frozen modern system order"
        )
    return {
        "schema_id": str(
            _require_exact(
                "bridge_artifact.matching.schema_id",
                _require_text("bridge_artifact.matching.schema_id", payload["schema_id"]),
                expected=MODERN_POINT_MATCHING_SCHEMA_ID,
            )
        ),
        "lane_id": str(
            _require_exact(
                "bridge_artifact.matching.lane_id",
                _require_text("bridge_artifact.matching.lane_id", payload["lane_id"]),
                expected=MODERN_LANE_ID,
            )
        ),
        "couplings_schema_id": str(
            _require_exact(
                "bridge_artifact.matching.couplings_schema_id",
                _require_text(
                    "bridge_artifact.matching.couplings_schema_id",
                    payload["couplings_schema_id"],
                ),
                expected=MODERN_POINT_COUPLINGS_SCHEMA_ID,
            )
        ),
        "input_bundle_schema_id": str(
            _require_exact(
                "bridge_artifact.matching.input_bundle_schema_id",
                _require_text(
                    "bridge_artifact.matching.input_bundle_schema_id",
                    payload["input_bundle_schema_id"],
                ),
                expected=MODERN_DEFAULT_INPUTS_SCHEMA_ID,
            )
        ),
        "input_bundle_id": _require_text(
            "bridge_artifact.matching.input_bundle_id", payload["input_bundle_id"]
        ),
        "input_provenance_id": _require_text(
            "bridge_artifact.matching.input_provenance_id", payload["input_provenance_id"]
        ),
        "input_resolution_policy_id": _require_text(
            "bridge_artifact.matching.input_resolution_policy_id",
            payload["input_resolution_policy_id"],
        ),
        "qcd_metadata_id": _require_text(
            "bridge_artifact.matching.qcd_metadata_id", payload["qcd_metadata_id"]
        ),
        "alpha_s_policy_id": _require_text(
            "bridge_artifact.matching.alpha_s_policy_id",
            payload["alpha_s_policy_id"],
        ),
        "operator_basis_id": _require_text(
            "bridge_artifact.matching.operator_basis_id",
            payload["operator_basis_id"],
        ),
        "weight_policy_id": _require_text(
            "bridge_artifact.matching.weight_policy_id",
            payload["weight_policy_id"],
        ),
        "point_id": _require_text("bridge_artifact.matching.point_id", payload["point_id"]),
        "point_label": _require_text(
            "bridge_artifact.matching.point_label", payload["point_label"]
        ),
        "Lambda_IR": _require_positive_float(
            "bridge_artifact.matching.Lambda_IR", payload["Lambda_IR"]
        ),
        "M_KK": _require_positive_float("bridge_artifact.matching.M_KK", payload["M_KK"]),
        "xi_KK": _require_positive_float(
            "bridge_artifact.matching.xi_KK", payload["xi_KK"]
        ),
        "alpha_s": _require_positive_float(
            "bridge_artifact.matching.alpha_s", payload["alpha_s"]
        ),
        "g_s": _require_positive_float("bridge_artifact.matching.g_s", payload["g_s"]),
        "system_matches": system_matches,
        "notes": _require_text("bridge_artifact.matching.notes", payload["notes"]),
    }


def _source_attr(source: Any, name: str) -> Any:
    if not hasattr(source, name):
        raise ArtifactSchemaError(f"source must provide {name}")
    return getattr(source, name)


@dataclass(frozen=True)
class ModernPointBridgeArtifactV1:
    """Frozen sidecar artifact carrying the explicit modern QS2 bridge."""

    schema_id: str = MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID
    schema_version: int = MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION
    lane_id: str = MODERN_LANE_ID
    input_bundle_schema_id: str = MODERN_DEFAULT_INPUTS_SCHEMA_ID
    input_bundle_id: str = ""
    input_provenance_id: str = ""
    input_resolution_policy_id: str = ""
    qcd_metadata_id: str = ""
    alpha_s_policy_id: str = ""
    point_id: str = ""
    point_label: str = ""
    coupling_schema_id: str = MODERN_POINT_COUPLINGS_SCHEMA_ID
    matching_schema_id: str = MODERN_POINT_MATCHING_SCHEMA_ID
    couplings: Mapping[str, Any] = field(default_factory=dict)
    matching: Mapping[str, Any] = field(default_factory=dict)
    notes: str = (
        "Explicit modern QS2 bridge sidecar only. This artifact exports the "
        "pointwise coupling and matching bridge with schema/lineage checks, "
        "but it does not claim a full EFT or RG package."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            str(
                _require_exact(
                    "schema_id",
                    _require_text("schema_id", self.schema_id),
                    expected=MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID,
                )
            ),
        )
        object.__setattr__(self, "schema_version", _require_int("schema_version", self.schema_version))
        if self.schema_version != MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION:
            raise ArtifactSchemaError(
                f"schema_version must be exactly {MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION}"
            )
        object.__setattr__(
            self,
            "lane_id",
            str(_require_exact("lane_id", _require_text("lane_id", self.lane_id), expected=MODERN_LANE_ID)),
        )
        object.__setattr__(
            self,
            "input_bundle_schema_id",
            str(
                _require_exact(
                    "input_bundle_schema_id",
                    _require_text("input_bundle_schema_id", self.input_bundle_schema_id),
                    expected=MODERN_DEFAULT_INPUTS_SCHEMA_ID,
                )
            ),
        )
        object.__setattr__(self, "input_bundle_id", _require_text("input_bundle_id", self.input_bundle_id))
        object.__setattr__(
            self,
            "input_provenance_id",
            _require_text("input_provenance_id", self.input_provenance_id),
        )
        object.__setattr__(
            self,
            "input_resolution_policy_id",
            _require_text("input_resolution_policy_id", self.input_resolution_policy_id),
        )
        object.__setattr__(self, "qcd_metadata_id", _require_text("qcd_metadata_id", self.qcd_metadata_id))
        object.__setattr__(
            self,
            "alpha_s_policy_id",
            _require_text("alpha_s_policy_id", self.alpha_s_policy_id),
        )
        object.__setattr__(self, "point_id", _require_text("point_id", self.point_id))
        object.__setattr__(self, "point_label", _require_text("point_label", self.point_label))
        object.__setattr__(
            self,
            "coupling_schema_id",
            str(
                _require_exact(
                    "coupling_schema_id",
                    _require_text("coupling_schema_id", self.coupling_schema_id),
                    expected=MODERN_POINT_COUPLINGS_SCHEMA_ID,
                )
            ),
        )
        object.__setattr__(
            self,
            "matching_schema_id",
            str(
                _require_exact(
                    "matching_schema_id",
                    _require_text("matching_schema_id", self.matching_schema_id),
                    expected=MODERN_POINT_MATCHING_SCHEMA_ID,
                )
            ),
        )
        canonical_couplings = _canonical_couplings_payload(self.couplings)
        canonical_matching = _canonical_matching_payload(self.matching)
        if canonical_couplings["point_id"] != self.point_id or canonical_matching["point_id"] != self.point_id:
            raise ArtifactSchemaError("bridge point_id must match nested couplings and matching")
        if canonical_couplings["point_label"] != self.point_label or canonical_matching["point_label"] != self.point_label:
            raise ArtifactSchemaError("bridge point_label must match nested couplings and matching")
        if canonical_couplings["input_bundle_id"] != self.input_bundle_id or canonical_matching["input_bundle_id"] != self.input_bundle_id:
            raise ArtifactSchemaError("bridge input_bundle_id must match nested couplings and matching")
        if canonical_couplings["input_provenance_id"] != self.input_provenance_id or canonical_matching["input_provenance_id"] != self.input_provenance_id:
            raise ArtifactSchemaError("bridge input_provenance_id must match nested couplings and matching")
        if canonical_couplings["input_resolution_policy_id"] != self.input_resolution_policy_id or canonical_matching["input_resolution_policy_id"] != self.input_resolution_policy_id:
            raise ArtifactSchemaError("bridge input_resolution_policy_id must match nested couplings and matching")
        if canonical_couplings["input_bundle_schema_id"] != self.input_bundle_schema_id or canonical_matching["input_bundle_schema_id"] != self.input_bundle_schema_id:
            raise ArtifactSchemaError("bridge input_bundle_schema_id must match nested couplings and matching")
        if canonical_couplings["qcd_metadata_id"] != self.qcd_metadata_id or canonical_matching["qcd_metadata_id"] != self.qcd_metadata_id:
            raise ArtifactSchemaError("bridge qcd_metadata_id must match nested couplings and matching")
        if canonical_couplings["alpha_s_policy_id"] != self.alpha_s_policy_id or canonical_matching["alpha_s_policy_id"] != self.alpha_s_policy_id:
            raise ArtifactSchemaError("bridge alpha_s_policy_id must match nested couplings and matching")
        if canonical_matching["couplings_schema_id"] != self.coupling_schema_id:
            raise ArtifactSchemaError("matching.couplings_schema_id must match coupling_schema_id")
        if not math.isclose(
            float(canonical_couplings["M_KK"]),
            float(canonical_matching["M_KK"]),
            rel_tol=0.0,
            abs_tol=0.0,
        ):
            raise ArtifactSchemaError("couplings.M_KK and matching.M_KK must agree")
        if not math.isclose(
            float(canonical_couplings["xi_KK"]),
            float(canonical_matching["xi_KK"]),
            rel_tol=0.0,
            abs_tol=0.0,
        ):
            raise ArtifactSchemaError("couplings.xi_KK and matching.xi_KK must agree")
        object.__setattr__(self, "couplings", canonical_couplings)
        object.__setattr__(self, "matching", canonical_matching)
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    @property
    def M_KK(self) -> float:
        return float(self.couplings["M_KK"])

    @property
    def xi_KK(self) -> float:
        return float(self.couplings["xi_KK"])

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_id": self.schema_id,
            "schema_version": self.schema_version,
            "lane_id": self.lane_id,
            "input_bundle_schema_id": self.input_bundle_schema_id,
            "input_bundle_id": self.input_bundle_id,
            "input_provenance_id": self.input_provenance_id,
            "input_resolution_policy_id": self.input_resolution_policy_id,
            "qcd_metadata_id": self.qcd_metadata_id,
            "alpha_s_policy_id": self.alpha_s_policy_id,
            "point_id": self.point_id,
            "point_label": self.point_label,
            "coupling_schema_id": self.coupling_schema_id,
            "matching_schema_id": self.matching_schema_id,
            "couplings": _deep_clone_json(self.couplings),
            "matching": _deep_clone_json(self.matching),
            "notes": self.notes,
        }

    def to_json(self) -> str:
        return json.dumps(self.as_dict(), indent=2, sort_keys=True, allow_nan=False) + "\n"

    def write_json(self, path: str | Path) -> None:
        Path(path).write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ModernPointBridgeArtifactV1":
        mapping = _require_mapping(data, context="bridge_artifact")
        _require_exact_keys(
            mapping,
            expected=(
                "schema_id",
                "schema_version",
                "lane_id",
                "input_bundle_schema_id",
                "input_bundle_id",
                "input_provenance_id",
                "input_resolution_policy_id",
                "qcd_metadata_id",
                "alpha_s_policy_id",
                "point_id",
                "point_label",
                "coupling_schema_id",
                "matching_schema_id",
                "couplings",
                "matching",
                "notes",
            ),
            context="bridge_artifact",
        )
        return cls(
            schema_id=_require_text("bridge_artifact.schema_id", mapping["schema_id"]),
            schema_version=_require_int(
                "bridge_artifact.schema_version",
                mapping["schema_version"],
            ),
            lane_id=_require_text("bridge_artifact.lane_id", mapping["lane_id"]),
            input_bundle_schema_id=_require_text(
                "bridge_artifact.input_bundle_schema_id",
                mapping["input_bundle_schema_id"],
            ),
            input_bundle_id=_require_text(
                "bridge_artifact.input_bundle_id", mapping["input_bundle_id"]
            ),
            input_provenance_id=_require_text(
                "bridge_artifact.input_provenance_id",
                mapping["input_provenance_id"],
            ),
            input_resolution_policy_id=_require_text(
                "bridge_artifact.input_resolution_policy_id",
                mapping["input_resolution_policy_id"],
            ),
            qcd_metadata_id=_require_text(
                "bridge_artifact.qcd_metadata_id", mapping["qcd_metadata_id"]
            ),
            alpha_s_policy_id=_require_text(
                "bridge_artifact.alpha_s_policy_id",
                mapping["alpha_s_policy_id"],
            ),
            point_id=_require_text("bridge_artifact.point_id", mapping["point_id"]),
            point_label=_require_text("bridge_artifact.point_label", mapping["point_label"]),
            coupling_schema_id=_require_text(
                "bridge_artifact.coupling_schema_id",
                mapping["coupling_schema_id"],
            ),
            matching_schema_id=_require_text(
                "bridge_artifact.matching_schema_id",
                mapping["matching_schema_id"],
            ),
            couplings=_require_mapping(mapping["couplings"], context="bridge_artifact.couplings"),
            matching=_require_mapping(mapping["matching"], context="bridge_artifact.matching"),
            notes=_require_text("bridge_artifact.notes", mapping["notes"]),
        )

    @classmethod
    def from_json(cls, payload: str | bytes) -> "ModernPointBridgeArtifactV1":
        return cls.from_dict(json.loads(payload))

    @classmethod
    def from_source(cls, source: Any) -> "ModernPointBridgeArtifactV1":
        source_schema_id = _require_text("source.schema_id", _source_attr(source, "schema_id"))
        if source_schema_id != MODERN_POINT_EVALUATION_SCHEMA_ID:
            raise ArtifactSchemaError(
                f"source.schema_id must be exactly {MODERN_POINT_EVALUATION_SCHEMA_ID!r}"
            )
        couplings = _source_attr(source, "couplings")
        matching = _source_attr(source, "matching")
        if hasattr(couplings, "as_dict"):
            couplings = couplings.as_dict()
        if hasattr(matching, "as_dict"):
            matching = matching.as_dict()
        return cls(
            input_bundle_schema_id=_require_text(
                "source.input_bundle_schema_id",
                _source_attr(source, "input_bundle_schema_id"),
            ),
            input_bundle_id=_require_text(
                "source.input_bundle_id",
                _source_attr(source, "input_bundle_id"),
            ),
            input_provenance_id=_require_text(
                "source.input_provenance_id",
                _source_attr(source, "input_provenance_id"),
            ),
            input_resolution_policy_id=_require_text(
                "source.input_resolution_policy_id",
                _source_attr(source, "input_resolution_policy_id"),
            ),
            qcd_metadata_id=_require_text(
                "source.qcd_metadata_id",
                _require_mapping(couplings, context="source.couplings")["qcd_metadata_id"],
            ),
            alpha_s_policy_id=_require_text(
                "source.alpha_s_policy_id",
                _require_mapping(couplings, context="source.couplings")["alpha_s_policy_id"],
            ),
            point_id=_require_text("source.point_id", _source_attr(source, "point_id")),
            point_label=_require_text("source.point_label", _source_attr(source, "point_label")),
            couplings=_require_mapping(couplings, context="source.couplings"),
            matching=_require_mapping(matching, context="source.matching"),
        )


def build_modern_point_bridge_artifact(source: Any) -> ModernPointBridgeArtifactV1:
    """Build one frozen modern bridge sidecar artifact from a point evaluation."""

    return ModernPointBridgeArtifactV1.from_source(source)


def default_modern_point_bridge_artifact(source: Any) -> ModernPointBridgeArtifactV1:
    """Return the canonical frozen modern bridge sidecar artifact."""

    return build_modern_point_bridge_artifact(source)


def bridge_artifact_from_dict(data: Mapping[str, Any]) -> ModernPointBridgeArtifactV1:
    """Build a modern bridge sidecar artifact from JSON-compatible data."""

    return ModernPointBridgeArtifactV1.from_dict(data)


def read_modern_point_bridge_artifact(path: str | Path) -> ModernPointBridgeArtifactV1:
    """Read one modern bridge sidecar artifact from disk."""

    return bridge_artifact_from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


def write_modern_point_bridge_artifact(
    artifact: ModernPointBridgeArtifactV1,
    path: str | Path,
) -> Path:
    """Write one modern bridge sidecar artifact to disk."""

    destination = Path(path)
    destination.write_text(artifact.to_json(), encoding="utf-8")
    return destination


__all__ = [
    "ArtifactSchemaError",
    "MODERN_DEFAULT_INPUTS_SCHEMA_ID",
    "MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID",
    "MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION",
    "MODERN_POINT_COUPLINGS_SCHEMA_ID",
    "MODERN_POINT_MATCHING_SCHEMA_ID",
    "MODERN_POINT_MATCHING_SYSTEM_SCHEMA_ID",
    "ModernPointBridgeArtifactV1",
    "bridge_artifact_from_dict",
    "build_modern_point_bridge_artifact",
    "default_modern_point_bridge_artifact",
    "read_modern_point_bridge_artifact",
    "write_modern_point_bridge_artifact",
]
