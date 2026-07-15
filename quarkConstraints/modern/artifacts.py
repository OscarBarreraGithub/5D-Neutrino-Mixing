"""Frozen per-point artifact export for the modern quark lane."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any

from .conventions import MODERN_LANE_ID
from .phenomenology import (
    MODERN_PHENOMENOLOGY_POLICY_ID,
    MODERN_PHENOMENOLOGY_SCHEMA_ID,
    MODERN_PHENOMENOLOGY_SYSTEM_IDS,
    ModernPhenomenologyPolicy,
    ModernPhenomenologySystemPolicy,
    default_modern_phenomenology_policy,
)

MODERN_POINT_ARTIFACT_SCHEMA_ID = "quarkConstraints.modern.artifacts.point.v1"
MODERN_POINT_ARTIFACT_SCHEMA_VERSION = 1
MODERN_POINT_ARTIFACT_HEADER_SCHEMA_ID = "quarkConstraints.modern.artifacts.header.v1"
MODERN_POINT_ARTIFACT_VERDICT_SCHEMA_ID = "quarkConstraints.modern.evaluation.verdict.v1"
MODERN_POINT_ARTIFACT_EVALUATION_SCHEMA_ID = "quarkConstraints.modern.evaluation.v1"
MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS = ("K", "B_d", "B_s", "D0")
MODERN_POINT_ARTIFACT_REQUIRED_OBSERVABLE_IDS = ("epsilon_K", "B_d", "B_s", "D0")
MODERN_POINT_ARTIFACT_REQUIRED_POLICY_SYSTEM_IDS = MODERN_PHENOMENOLOGY_SYSTEM_IDS
MODERN_POINT_ARTIFACT_BACKEND_SYSTEM_IDS = ("K", "B_d", "B_s", "D")
MODERN_POINT_ARTIFACT_BACKEND_KEYS = ("epsilon_k", "b_d", "b_s", "d")
MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID = "quarkConstraints.modern.artifacts.bridge.v1"
MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION = 1
MODERN_POINT_BRIDGE_COUPLINGS_SCHEMA_ID = "quarkConstraints.modern.couplings.point.v1"
MODERN_POINT_BRIDGE_MATCHING_SCHEMA_ID = "quarkConstraints.modern.matching.point.v1"
MODERN_POINT_BRIDGE_MATCHING_SYSTEM_SCHEMA_ID = (
    "quarkConstraints.modern.matching.system.v1"
)
MODERN_POINT_BRIDGE_MATCHING_SYSTEM_IDS = ("K", "B_d", "B_s", "D0")


class ArtifactSchemaError(ValueError):
    """Raised when a modern artifact violates the frozen JSON contract."""


def _require_text(name: str, value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ArtifactSchemaError(f"{name} must be a non-empty string")
    return value.strip()


def _require_text_value(value: Any, *, context: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ArtifactSchemaError(f"{context} must be a non-empty string")
    return value.strip()


def _require_exact(name: str, value: str | int, *, expected: str | int) -> str | int:
    if value != expected:
        raise ArtifactSchemaError(f"{name} must be exactly {expected!r}")
    return value


def _require_int(name: str, value: int) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise ArtifactSchemaError(f"{name} must be an integer")
    return value


def _require_positive_float(name: str, value: float) -> float:
    numeric = float(value)
    if not math.isfinite(numeric) or numeric <= 0.0:
        raise ArtifactSchemaError(f"{name} must be positive")
    return numeric


def _require_nonnegative_float(name: str, value: float) -> float:
    numeric = float(value)
    if not math.isfinite(numeric) or numeric < 0.0:
        raise ArtifactSchemaError(f"{name} must be non-negative")
    return numeric


def _require_finite_float(name: str, value: float) -> float:
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ArtifactSchemaError(f"{name} must be a finite float")
    return numeric


def _require_mapping(value: Any, *, context: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ArtifactSchemaError(f"{context} must be a JSON object")
    return value


def _require_sequence(value: Any, *, context: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise ArtifactSchemaError(f"{context} must be a JSON array")
    return value


def _require_exact_keys(mapping: Mapping[str, Any], *, expected: tuple[str, ...], context: str) -> None:
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


def _dump_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _canonical_operator_sizes(weighted_operator_sizes: Mapping[str, float]) -> MappingProxyType[str, float]:
    canonical = dict(sorted((str(key), float(value)) for key, value in weighted_operator_sizes.items()))
    return MappingProxyType(canonical)


def _require_complex(name: str, value: Any) -> complex:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ArtifactSchemaError(f"{name} must be finite")
    return number


def _complex_payload(value: complex) -> dict[str, float]:
    return {"real": float(value.real), "imag": float(value.imag)}


def _complex_from_payload(name: str, payload: Any) -> complex:
    mapping = _require_mapping(payload, context=name)
    _require_exact_keys(mapping, expected=("real", "imag"), context=name)
    return _require_complex(name, float(mapping["real"]) + 1j * float(mapping["imag"]))


def _matrix_payload(matrix: tuple[tuple[complex, ...], ...]) -> dict[str, object]:
    return {
        "real": [[float(entry.real) for entry in row] for row in matrix],
        "imag": [[float(entry.imag) for entry in row] for row in matrix],
    }


def _matrix_from_payload(name: str, payload: Any) -> tuple[tuple[complex, ...], ...]:
    mapping = _require_mapping(payload, context=name)
    _require_exact_keys(mapping, expected=("real", "imag"), context=name)
    real_rows = tuple(_require_sequence(mapping["real"], context=f"{name}.real"))
    imag_rows = tuple(_require_sequence(mapping["imag"], context=f"{name}.imag"))
    if len(real_rows) != 3 or len(imag_rows) != 3:
        raise ArtifactSchemaError(f"{name} payload must contain 3x3 real and imag arrays")
    matrix: list[tuple[complex, ...]] = []
    for row_index, (real_row, imag_row) in enumerate(zip(real_rows, imag_rows)):
        real_entries = tuple(_require_sequence(real_row, context=f"{name}.real[{row_index}]"))
        imag_entries = tuple(_require_sequence(imag_row, context=f"{name}.imag[{row_index}]"))
        if len(real_entries) != 3 or len(imag_entries) != 3:
            raise ArtifactSchemaError(f"{name} payload rows must each contain three entries")
        matrix.append(
            tuple(
                _require_complex(
                    f"{name}[{row_index}][{column_index}]",
                    float(real_entries[column_index]) + 1j * float(imag_entries[column_index]),
                )
                for column_index in range(3)
            )
        )
    return tuple(matrix)


def _canonical_complex_matrix(name: str, value: Any) -> tuple[tuple[complex, ...], ...]:
    if isinstance(value, Mapping):
        return _matrix_from_payload(name, value)
    rows = tuple(_require_sequence(value, context=name))
    if len(rows) != 3:
        raise ArtifactSchemaError(f"{name} must contain exactly 3 rows")
    matrix: list[tuple[complex, ...]] = []
    for row_index, row in enumerate(rows):
        entries = tuple(_require_sequence(row, context=f"{name}[{row_index}]"))
        if len(entries) != 3:
            raise ArtifactSchemaError(
                f"{name}[{row_index}] must contain exactly 3 columns"
            )
        matrix.append(
            tuple(
                _require_complex(f"{name}[{row_index}][{column_index}]", entry)
                for column_index, entry in enumerate(entries)
            )
        )
    return tuple(matrix)


def _require_system_id(name: str, value: str, expected: str) -> str:
    value = _require_text(name, value)
    if value != expected:
        raise ValueError(f"{name} must be exactly {expected!r}")
    return value


def _policy_system_id_for_system(system_id: str) -> str:
    if system_id == "K":
        return "epsilon_K"
    if system_id in {"B_d", "B_s", "D0"}:
        return system_id
    raise ValueError(f"system_id must be one of {MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS!r}")


def _backend_system_id_for_system(system_id: str) -> str:
    try:
        index = MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS.index(system_id)
    except ValueError as exc:
        raise ValueError(f"system_id must be one of {MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS!r}") from exc
    return MODERN_POINT_ARTIFACT_BACKEND_SYSTEM_IDS[index]


def _backend_key_for_system(system_id: str) -> str:
    try:
        index = MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS.index(system_id)
    except ValueError as exc:
        raise ValueError(f"system_id must be one of {MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS!r}") from exc
    return MODERN_POINT_ARTIFACT_BACKEND_KEYS[index]


def _observable_id_for_system(system_id: str) -> str:
    try:
        index = MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS.index(system_id)
    except ValueError as exc:
        raise ValueError(f"system_id must be one of {MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS!r}") from exc
    return MODERN_POINT_ARTIFACT_REQUIRED_OBSERVABLE_IDS[index]


def _policy_system_for_system(
    policy: ModernPhenomenologyPolicy,
    system_id: str,
) -> ModernPhenomenologySystemPolicy:
    requested = _policy_system_id_for_system(system_id)
    return policy.system_policy(requested)


def _modern_system_policy_from_dict(data: Mapping[str, Any]) -> ModernPhenomenologySystemPolicy:
    mapping = _require_mapping(data, context="policy.system")
    _require_exact_keys(
        mapping,
        expected=("schema_id", "system_id", "policy_id", "display_name", "notes"),
        context="policy.system",
    )
    return ModernPhenomenologySystemPolicy(
        schema_id=_require_text_value(mapping["schema_id"], context="policy.system.schema_id"),
        system_id=_require_text_value(mapping["system_id"], context="policy.system.system_id"),
        policy_id=_require_text_value(mapping["policy_id"], context="policy.system.policy_id"),
        display_name=_require_text_value(
            mapping["display_name"], context="policy.system.display_name"
        ),
        notes=_require_text_value(mapping["notes"], context="policy.system.notes"),
    )


def _modern_policy_from_dict(data: Mapping[str, Any]) -> ModernPhenomenologyPolicy:
    mapping = _require_mapping(data, context="policy")
    _require_exact_keys(
        mapping,
        expected=("schema_id", "policy_id", "lane_id", "registry_schema_id", "systems", "notes"),
        context="policy",
    )
    systems = tuple(
        _modern_system_policy_from_dict(item)
        for item in _require_sequence(mapping["systems"], context="policy.systems")
    )
    return ModernPhenomenologyPolicy(
        schema_id=_require_text_value(mapping["schema_id"], context="policy.schema_id"),
        policy_id=_require_text_value(mapping["policy_id"], context="policy.policy_id"),
        lane_id=_require_text_value(mapping["lane_id"], context="policy.lane_id"),
        registry_schema_id=_require_text_value(
            mapping["registry_schema_id"], context="policy.registry_schema_id"
        ),
        systems=systems,
        notes=_require_text_value(mapping["notes"], context="policy.notes"),
    )


def _source_attr(source: Any, name: str) -> Any:
    if not hasattr(source, name):
        raise ArtifactSchemaError(f"source must provide {name}")
    return getattr(source, name)


def _bridge_payload_from_source(source: Any, *, context: str) -> Mapping[str, Any]:
    if isinstance(source, Mapping):
        return _require_mapping(source, context=context)
    if hasattr(source, "as_dict"):
        return _require_mapping(source.as_dict(), context=context)
    raise ArtifactSchemaError(f"{context} must be a mapping or provide as_dict()")


@dataclass(frozen=True)
class ModernPointArtifactHeader:
    """Frozen header for one exported modern point artifact."""

    schema_id: str = MODERN_POINT_ARTIFACT_HEADER_SCHEMA_ID
    lane_id: str = MODERN_LANE_ID
    point_id: str = ""
    point_label: str = ""
    policy_schema_id: str = MODERN_PHENOMENOLOGY_SCHEMA_ID
    policy_id: str = MODERN_PHENOMENOLOGY_POLICY_ID
    policy_system_ids: tuple[str, ...] = MODERN_POINT_ARTIFACT_REQUIRED_POLICY_SYSTEM_IDS
    evaluation_schema_id: str = MODERN_POINT_ARTIFACT_EVALUATION_SCHEMA_ID
    verdict_schema_id: str = MODERN_POINT_ARTIFACT_VERDICT_SCHEMA_ID
    verdict_system_ids: tuple[str, ...] = MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS
    verdict_observable_ids: tuple[str, ...] = MODERN_POINT_ARTIFACT_REQUIRED_OBSERVABLE_IDS
    M_KK: float = 0.0
    xi_KK: float = 0.0
    verdict_count: int = 4
    notes: str = (
        "Per-point modern artifact header only. No scan loop, caching, warm-starts, "
        "manifests, or SLURM plumbing are represented here."
    )

    def __post_init__(self) -> None:
        _require_system_id("schema_id", self.schema_id, MODERN_POINT_ARTIFACT_HEADER_SCHEMA_ID)
        _require_system_id("lane_id", self.lane_id, MODERN_LANE_ID)
        object.__setattr__(self, "point_id", _require_text("point_id", self.point_id))
        object.__setattr__(self, "point_label", _require_text("point_label", self.point_label))
        _require_system_id("policy_schema_id", self.policy_schema_id, MODERN_PHENOMENOLOGY_SCHEMA_ID)
        _require_system_id("policy_id", self.policy_id, MODERN_PHENOMENOLOGY_POLICY_ID)
        _require_system_id(
            "evaluation_schema_id",
            self.evaluation_schema_id,
            MODERN_POINT_ARTIFACT_EVALUATION_SCHEMA_ID,
        )
        _require_system_id(
            "verdict_schema_id",
            self.verdict_schema_id,
            MODERN_POINT_ARTIFACT_VERDICT_SCHEMA_ID,
        )
        if tuple(self.policy_system_ids) != MODERN_POINT_ARTIFACT_REQUIRED_POLICY_SYSTEM_IDS:
            raise ValueError(
                "policy_system_ids must explicitly cover epsilon_K, K, B_d, B_s, and D0"
            )
        if tuple(self.verdict_system_ids) != MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS:
            raise ValueError("verdict_system_ids must explicitly cover K, B_d, B_s, and D0")
        if tuple(self.verdict_observable_ids) != MODERN_POINT_ARTIFACT_REQUIRED_OBSERVABLE_IDS:
            raise ValueError(
                "verdict_observable_ids must explicitly cover epsilon_K, B_d, B_s, and D0"
            )
        object.__setattr__(self, "M_KK", _require_positive_float("M_KK", self.M_KK))
        object.__setattr__(self, "xi_KK", _require_positive_float("xi_KK", self.xi_KK))
        object.__setattr__(self, "verdict_count", _require_int("verdict_count", self.verdict_count))
        if self.verdict_count != len(MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS):
            raise ValueError("verdict_count must be 4 for the frozen modern point artifact")
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "lane_id": self.lane_id,
            "point_id": self.point_id,
            "point_label": self.point_label,
            "policy_schema_id": self.policy_schema_id,
            "policy_id": self.policy_id,
            "policy_system_ids": list(self.policy_system_ids),
            "evaluation_schema_id": self.evaluation_schema_id,
            "verdict_schema_id": self.verdict_schema_id,
            "verdict_system_ids": list(self.verdict_system_ids),
            "verdict_observable_ids": list(self.verdict_observable_ids),
            "M_KK": self.M_KK,
            "xi_KK": self.xi_KK,
            "verdict_count": self.verdict_count,
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ModernPointArtifactHeader":
        mapping = _require_mapping(data, context="header")
        _require_exact_keys(
            mapping,
            expected=(
                "schema_id",
                "lane_id",
                "point_id",
                "point_label",
                "policy_schema_id",
                "policy_id",
                "policy_system_ids",
                "evaluation_schema_id",
                "verdict_schema_id",
                "verdict_system_ids",
                "verdict_observable_ids",
                "M_KK",
                "xi_KK",
                "verdict_count",
                "notes",
            ),
            context="header",
        )
        return cls(
            schema_id=_require_text_value(mapping["schema_id"], context="header.schema_id"),
            lane_id=_require_text_value(mapping["lane_id"], context="header.lane_id"),
            point_id=_require_text_value(mapping["point_id"], context="header.point_id"),
            point_label=_require_text_value(mapping["point_label"], context="header.point_label"),
            policy_schema_id=_require_text_value(
                mapping["policy_schema_id"], context="header.policy_schema_id"
            ),
            policy_id=_require_text_value(mapping["policy_id"], context="header.policy_id"),
            policy_system_ids=tuple(
                _require_text_value(item, context="header.policy_system_ids")
                for item in _require_sequence(mapping["policy_system_ids"], context="header.policy_system_ids")
            ),
            evaluation_schema_id=_require_text_value(
                mapping["evaluation_schema_id"], context="header.evaluation_schema_id"
            ),
            verdict_schema_id=_require_text_value(
                mapping["verdict_schema_id"], context="header.verdict_schema_id"
            ),
            verdict_system_ids=tuple(
                _require_text_value(item, context="header.verdict_system_ids")
                for item in _require_sequence(mapping["verdict_system_ids"], context="header.verdict_system_ids")
            ),
            verdict_observable_ids=tuple(
                _require_text_value(item, context="header.verdict_observable_ids")
                for item in _require_sequence(
                    mapping["verdict_observable_ids"], context="header.verdict_observable_ids"
                )
            ),
            M_KK=_require_positive_float("header.M_KK", mapping["M_KK"]),
            xi_KK=_require_positive_float("header.xi_KK", mapping["xi_KK"]),
            verdict_count=_require_int("header.verdict_count", mapping["verdict_count"]),
            notes=_require_text_value(mapping["notes"], context="header.notes"),
        )


@dataclass(frozen=True)
class ModernPointArtifactVerdict:
    """Frozen verdict record for one explicit modern neutral-meson system."""

    schema_id: str = MODERN_POINT_ARTIFACT_VERDICT_SCHEMA_ID
    system_id: str = ""
    observable_id: str = ""
    policy_system_id: str = ""
    backend_system_id: str = ""
    backend_key: str = ""
    policy_id: str = ""
    policy_display_name: str = ""
    policy_notes: str = ""
    passes: bool = False
    ratio_to_bound: float = 0.0
    effective_amplitude: float = 0.0
    coherent_amplitude: float = 0.0
    bound: float = 0.0
    dominant_operator: str = ""
    dominant_operator_size: float = 0.0
    weighted_operator_sizes: Mapping[str, float] = field(default_factory=dict)
    note: str = ""

    def __post_init__(self) -> None:
        _require_system_id("schema_id", self.schema_id, MODERN_POINT_ARTIFACT_VERDICT_SCHEMA_ID)
        object.__setattr__(self, "system_id", _require_text("system_id", self.system_id))
        object.__setattr__(self, "observable_id", _require_text("observable_id", self.observable_id))
        object.__setattr__(self, "policy_system_id", _require_text("policy_system_id", self.policy_system_id))
        object.__setattr__(self, "backend_system_id", _require_text("backend_system_id", self.backend_system_id))
        object.__setattr__(self, "backend_key", _require_text("backend_key", self.backend_key))
        object.__setattr__(self, "policy_id", _require_text("policy_id", self.policy_id))
        object.__setattr__(self, "policy_display_name", _require_text("policy_display_name", self.policy_display_name))
        object.__setattr__(self, "policy_notes", _require_text("policy_notes", self.policy_notes))
        if not isinstance(self.passes, bool):
            raise ValueError("passes must be a bool")
        object.__setattr__(self, "ratio_to_bound", _require_nonnegative_float("ratio_to_bound", self.ratio_to_bound))
        object.__setattr__(
            self,
            "effective_amplitude",
            _require_finite_float("effective_amplitude", self.effective_amplitude),
        )
        object.__setattr__(
            self,
            "coherent_amplitude",
            _require_finite_float("coherent_amplitude", self.coherent_amplitude),
        )
        object.__setattr__(self, "bound", _require_positive_float("bound", self.bound))
        object.__setattr__(self, "dominant_operator", _require_text("dominant_operator", self.dominant_operator))
        object.__setattr__(
            self,
            "dominant_operator_size",
            _require_nonnegative_float("dominant_operator_size", self.dominant_operator_size),
        )
        object.__setattr__(self, "note", _require_text("note", self.note))
        if self.system_id not in MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS:
            raise ValueError(
                f"system_id must be one of {MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS!r}"
            )
        expected_observable = _observable_id_for_system(self.system_id)
        expected_policy_system = _policy_system_id_for_system(self.system_id)
        expected_backend_system = _backend_system_id_for_system(self.system_id)
        expected_backend_key = _backend_key_for_system(self.system_id)
        if self.observable_id != expected_observable:
            raise ValueError(
                f"observable_id for {self.system_id!r} must be exactly {expected_observable!r}"
            )
        if self.policy_system_id != expected_policy_system:
            raise ValueError(
                f"policy_system_id for {self.system_id!r} must be exactly {expected_policy_system!r}"
            )
        if self.backend_system_id != expected_backend_system:
            raise ValueError(
                f"backend_system_id for {self.system_id!r} must be exactly {expected_backend_system!r}"
            )
        if self.backend_key != expected_backend_key:
            raise ValueError(
                f"backend_key for {self.system_id!r} must be exactly {expected_backend_key!r}"
            )
        if not isinstance(self.weighted_operator_sizes, Mapping):
            raise ValueError("weighted_operator_sizes must be a mapping")
        canonical_operator_sizes = dict(
            sorted(
                (
                    str(key),
                    _require_nonnegative_float("weighted_operator_sizes", value),
                )
                for key, value in self.weighted_operator_sizes.items()
            )
        )
        if not canonical_operator_sizes:
            raise ValueError("weighted_operator_sizes must contain at least one operator")
        if self.dominant_operator not in canonical_operator_sizes:
            raise ValueError("dominant_operator must be present in weighted_operator_sizes")
        dominant_size = canonical_operator_sizes[self.dominant_operator]
        max_size = max(canonical_operator_sizes.values())
        if not math.isclose(
            self.dominant_operator_size,
            dominant_size,
            rel_tol=0.0,
            abs_tol=0.0,
        ):
            raise ValueError("dominant_operator_size must match weighted_operator_sizes[dominant_operator]")
        if not math.isclose(dominant_size, max_size, rel_tol=0.0, abs_tol=0.0):
            raise ValueError("dominant_operator must name an operator with the maximum size")
        if not math.isclose(
            self.dominant_operator_size,
            max_size,
            rel_tol=0.0,
            abs_tol=0.0,
        ):
            raise ValueError("dominant_operator_size must equal the maximum operator size")
        object.__setattr__(
            self,
            "weighted_operator_sizes",
            MappingProxyType(canonical_operator_sizes),
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "system_id": self.system_id,
            "observable_id": self.observable_id,
            "policy_system_id": self.policy_system_id,
            "backend_system_id": self.backend_system_id,
            "backend_key": self.backend_key,
            "policy_id": self.policy_id,
            "policy_display_name": self.policy_display_name,
            "policy_notes": self.policy_notes,
            "passes": self.passes,
            "ratio_to_bound": self.ratio_to_bound,
            "effective_amplitude": self.effective_amplitude,
            "coherent_amplitude": self.coherent_amplitude,
            "bound": self.bound,
            "dominant_operator": self.dominant_operator,
            "dominant_operator_size": self.dominant_operator_size,
            "weighted_operator_sizes": dict(self.weighted_operator_sizes),
            "note": self.note,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ModernPointArtifactVerdict":
        mapping = _require_mapping(data, context="verdict")
        _require_exact_keys(
            mapping,
            expected=(
                "schema_id",
                "system_id",
                "observable_id",
                "policy_system_id",
                "backend_system_id",
                "backend_key",
                "policy_id",
                "policy_display_name",
                "policy_notes",
                "passes",
                "ratio_to_bound",
                "effective_amplitude",
                "coherent_amplitude",
                "bound",
                "dominant_operator",
                "dominant_operator_size",
                "weighted_operator_sizes",
                "note",
            ),
            context="verdict",
        )
        return cls(
            schema_id=_require_text_value(mapping["schema_id"], context="verdict.schema_id"),
            system_id=_require_text_value(mapping["system_id"], context="verdict.system_id"),
            observable_id=_require_text_value(mapping["observable_id"], context="verdict.observable_id"),
            policy_system_id=_require_text_value(
                mapping["policy_system_id"], context="verdict.policy_system_id"
            ),
            backend_system_id=_require_text_value(
                mapping["backend_system_id"], context="verdict.backend_system_id"
            ),
            backend_key=_require_text_value(mapping["backend_key"], context="verdict.backend_key"),
            policy_id=_require_text_value(mapping["policy_id"], context="verdict.policy_id"),
            policy_display_name=_require_text_value(
                mapping["policy_display_name"], context="verdict.policy_display_name"
            ),
            policy_notes=_require_text_value(mapping["policy_notes"], context="verdict.policy_notes"),
            passes=mapping["passes"],
            ratio_to_bound=_require_nonnegative_float(
                "verdict.ratio_to_bound", mapping["ratio_to_bound"]
            ),
            effective_amplitude=_require_finite_float(
                "verdict.effective_amplitude", mapping["effective_amplitude"]
            ),
            coherent_amplitude=_require_finite_float(
                "verdict.coherent_amplitude", mapping["coherent_amplitude"]
            ),
            bound=_require_positive_float("verdict.bound", mapping["bound"]),
            dominant_operator=_require_text_value(
                mapping["dominant_operator"], context="verdict.dominant_operator"
            ),
            dominant_operator_size=_require_nonnegative_float(
                "verdict.dominant_operator_size", mapping["dominant_operator_size"]
            ),
            weighted_operator_sizes={
                _require_text_value(item, context="verdict.weighted_operator_sizes.key"): _require_nonnegative_float(
                    "verdict.weighted_operator_sizes.value", value
                )
                for item, value in _require_mapping(
                    mapping["weighted_operator_sizes"], context="verdict.weighted_operator_sizes"
                ).items()
            },
            note=_require_text_value(mapping["note"], context="verdict.note"),
        )


@dataclass(frozen=True)
class ModernPointArtifactV1:
    """Frozen per-point modern artifact export."""

    schema_id: str = MODERN_POINT_ARTIFACT_SCHEMA_ID
    schema_version: int = MODERN_POINT_ARTIFACT_SCHEMA_VERSION
    header: ModernPointArtifactHeader = field(default_factory=ModernPointArtifactHeader)
    policy: ModernPhenomenologyPolicy = field(default_factory=default_modern_phenomenology_policy)
    verdicts: tuple[ModernPointArtifactVerdict, ...] = ()
    notes: str = (
        "Deterministic per-point modern artifact export. It carries no scan loop, "
        "caching, warm-starts, manifests, or SLURM plumbing."
    )

    def __post_init__(self) -> None:
        _require_system_id("schema_id", self.schema_id, MODERN_POINT_ARTIFACT_SCHEMA_ID)
        _require_int("schema_version", self.schema_version)
        if self.schema_version != MODERN_POINT_ARTIFACT_SCHEMA_VERSION:
            raise ValueError(
                f"schema_version must be exactly {MODERN_POINT_ARTIFACT_SCHEMA_VERSION}"
            )
        if not isinstance(self.header, ModernPointArtifactHeader):
            raise ValueError("header must be a ModernPointArtifactHeader instance")
        if not isinstance(self.policy, ModernPhenomenologyPolicy):
            raise ValueError("policy must be a ModernPhenomenologyPolicy instance")
        object.__setattr__(self, "notes", _require_text("notes", self.notes))
        if tuple(self.policy.system_ids) != MODERN_PHENOMENOLOGY_SYSTEM_IDS:
            raise ValueError(
                "policy must explicitly cover epsilon_K, K, B_d, B_s, and D0"
            )
        if self.header.policy_schema_id != self.policy.schema_id:
            raise ValueError("header.policy_schema_id must match the exported policy")
        if self.header.policy_id != self.policy.policy_id:
            raise ValueError("header.policy_id must match the exported policy")
        if tuple(self.header.policy_system_ids) != self.policy.system_ids:
            raise ValueError("header.policy_system_ids must match the exported policy")
        if self.header.verdict_count != len(self.verdicts):
            raise ValueError("header.verdict_count must match the number of verdicts")
        if not isinstance(self.verdicts, tuple):
            raise ValueError("verdicts must be a tuple of ModernPointArtifactVerdict instances")
        normalized_verdicts = tuple(self.verdicts)
        if len(normalized_verdicts) != len(MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS):
            raise ValueError("verdicts must contain exactly K, B_d, B_s, and D0")
        for verdict in normalized_verdicts:
            if not isinstance(verdict, ModernPointArtifactVerdict):
                raise ValueError("verdicts must contain only ModernPointArtifactVerdict instances")
        normalized_system_ids = tuple(verdict.system_id for verdict in normalized_verdicts)
        if normalized_system_ids != MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS:
            unknown = tuple(
                system_id
                for system_id in normalized_system_ids
                if system_id not in MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS
            )
            missing = tuple(
                system_id
                for system_id in MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS
                if system_id not in normalized_system_ids
            )
            if unknown:
                raise ValueError(f"verdicts contains unknown systems: {unknown!r}")
            if missing:
                raise ValueError(f"verdicts is missing systems: {missing!r}")
            raise ValueError("verdicts must preserve the frozen modern artifact order")
        if tuple(verdict.observable_id for verdict in normalized_verdicts) != MODERN_POINT_ARTIFACT_REQUIRED_OBSERVABLE_IDS:
            raise ValueError(
                "verdicts must explicitly cover epsilon_K, B_d, B_s, and D0 in order"
            )
        for verdict in normalized_verdicts:
            policy_system = self.policy.system_policy(verdict.policy_system_id)
            if verdict.policy_id != policy_system.policy_id:
                raise ValueError(
                    f"verdict policy_id for {verdict.system_id!r} must match the supplied policy entry"
                )
            if verdict.policy_display_name != policy_system.display_name:
                raise ValueError(
                    f"verdict policy_display_name for {verdict.system_id!r} must match the supplied policy entry"
                )
            if verdict.policy_notes != policy_system.notes:
                raise ValueError(
                    f"verdict policy_notes for {verdict.system_id!r} must match the supplied policy entry"
                )
        object.__setattr__(self, "verdicts", normalized_verdicts)

    @property
    def system_ids(self) -> tuple[str, ...]:
        return tuple(verdict.system_id for verdict in self.verdicts)

    @property
    def verdict_by_system(self) -> dict[str, ModernPointArtifactVerdict]:
        return {verdict.system_id: verdict for verdict in self.verdicts}

    @property
    def point_id(self) -> str:
        return self.header.point_id

    @property
    def point_label(self) -> str:
        return self.header.point_label

    @property
    def M_KK(self) -> float:
        return self.header.M_KK

    @property
    def xi_KK(self) -> float:
        return self.header.xi_KK

    @property
    def passes_all(self) -> bool:
        return all(verdict.passes for verdict in self.verdicts)

    @property
    def max_ratio_to_bound(self) -> float:
        if not self.verdicts:
            return 0.0
        return max(verdict.ratio_to_bound for verdict in self.verdicts)

    def verdict_for(self, system_id: str) -> ModernPointArtifactVerdict:
        requested = _require_text("system_id", system_id)
        for verdict in self.verdicts:
            if verdict.system_id == requested:
                return verdict
        raise ValueError(
            f"system_id must be one of {MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS!r}"
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "schema_version": self.schema_version,
            "header": self.header.as_dict(),
            "policy": self.policy.as_dict(),
            "verdicts": [verdict.as_dict() for verdict in self.verdicts],
            "notes": self.notes,
        }

    def to_json(self) -> str:
        return _dump_json(self.as_dict())

    def write_json(self, path: str | Path) -> None:
        Path(path).write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ModernPointArtifactV1":
        mapping = _require_mapping(data, context="artifact")
        _require_exact_keys(
            mapping,
            expected=("schema_id", "schema_version", "header", "policy", "verdicts", "notes"),
            context="artifact",
        )
        header = ModernPointArtifactHeader.from_dict(
            _require_mapping(mapping["header"], context="artifact.header")
        )
        policy = _modern_policy_from_dict(
            _require_mapping(mapping["policy"], context="artifact.policy")
        )
        verdicts = tuple(
            ModernPointArtifactVerdict.from_dict(_require_mapping(item, context="artifact.verdicts"))
            for item in _require_sequence(mapping["verdicts"], context="artifact.verdicts")
        )
        return cls(
            schema_id=_require_text_value(mapping["schema_id"], context="artifact.schema_id"),
            schema_version=_require_int("artifact.schema_version", mapping["schema_version"]),
            header=header,
            policy=policy,
            verdicts=verdicts,
            notes=_require_text_value(mapping["notes"], context="artifact.notes"),
        )

    @classmethod
    def from_json(cls, payload: str | bytes) -> "ModernPointArtifactV1":
        return cls.from_dict(json.loads(payload))

    @classmethod
    def from_source(cls, source: Any, *, policy: ModernPhenomenologyPolicy | None = None) -> "ModernPointArtifactV1":
        """Build an artifact from one point evaluation-like source object."""
        source_schema_id = _require_text("source.schema_id", _source_attr(source, "schema_id"))
        if source_schema_id != MODERN_POINT_ARTIFACT_EVALUATION_SCHEMA_ID:
            raise ValueError(
                f"source.schema_id must be exactly {MODERN_POINT_ARTIFACT_EVALUATION_SCHEMA_ID!r}"
            )
        resolved_policy = policy if policy is not None else _source_attr(source, "policy")
        if not isinstance(resolved_policy, ModernPhenomenologyPolicy):
            raise ValueError("policy must be a ModernPhenomenologyPolicy instance")
        point_label = _require_text("source.point_label", _source_attr(source, "point_label"))
        M_KK = _require_positive_float("source.M_KK", _source_attr(source, "M_KK"))
        xi_KK = _require_positive_float("source.xi_KK", _source_attr(source, "xi_KK"))
        verdicts = tuple(
            ModernPointArtifactVerdict(
                schema_id=_require_text_value(_source_attr(verdict, "schema_id"), context="verdict.schema_id"),
                system_id=_require_text_value(_source_attr(verdict, "system_id"), context="verdict.system_id"),
                observable_id=_require_text_value(
                    _source_attr(verdict, "observable_id"), context="verdict.observable_id"
                ),
                policy_system_id=_require_text_value(
                    _source_attr(verdict, "policy_system_id"), context="verdict.policy_system_id"
                ),
                backend_system_id=_require_text_value(
                    _source_attr(verdict, "backend_system_id"), context="verdict.backend_system_id"
                ),
                backend_key=_require_text_value(_source_attr(verdict, "backend_key"), context="verdict.backend_key"),
                policy_id=_require_text_value(_source_attr(verdict, "policy_id"), context="verdict.policy_id"),
                policy_display_name=_require_text_value(
                    _source_attr(verdict, "policy_display_name"), context="verdict.policy_display_name"
                ),
                policy_notes=_require_text_value(_source_attr(verdict, "policy_notes"), context="verdict.policy_notes"),
                passes=_source_attr(verdict, "passes"),
                ratio_to_bound=_require_nonnegative_float(
                    "verdict.ratio_to_bound", _source_attr(verdict, "ratio_to_bound")
                ),
                effective_amplitude=_require_finite_float(
                    "verdict.effective_amplitude", _source_attr(verdict, "effective_amplitude")
                ),
                coherent_amplitude=_require_finite_float(
                    "verdict.coherent_amplitude", _source_attr(verdict, "coherent_amplitude")
                ),
                bound=_require_positive_float("verdict.bound", _source_attr(verdict, "bound")),
                dominant_operator=_require_text_value(
                    _source_attr(verdict, "dominant_operator"), context="verdict.dominant_operator"
                ),
                dominant_operator_size=_require_nonnegative_float(
                    "verdict.dominant_operator_size", _source_attr(verdict, "dominant_operator_size")
                ),
                weighted_operator_sizes=dict(_source_attr(verdict, "weighted_operator_sizes")),
                note=_require_text_value(_source_attr(verdict, "note"), context="verdict.note"),
            )
            for verdict in _require_sequence(_source_attr(source, "verdicts"), context="source.verdicts")
        )
        source_point_id = _source_attr(source, "point_id") if hasattr(source, "point_id") else point_label
        header = ModernPointArtifactHeader(
            point_id=_require_text_value(source_point_id, context="source.point_id"),
            point_label=point_label,
            policy_system_ids=tuple(resolved_policy.system_ids),
            M_KK=M_KK,
            xi_KK=xi_KK,
            verdict_count=len(verdicts),
        )
        return cls(header=header, policy=resolved_policy, verdicts=verdicts)


def build_modern_point_artifact(
    source: Any,
    *,
    policy: ModernPhenomenologyPolicy | None = None,
) -> ModernPointArtifactV1:
    """Build one frozen modern per-point artifact from a point-evaluation source."""

    return ModernPointArtifactV1.from_source(source, policy=policy)


def default_modern_point_artifact(
    source: Any,
    *,
    policy: ModernPhenomenologyPolicy | None = None,
) -> ModernPointArtifactV1:
    """Return the canonical frozen modern point artifact."""

    return build_modern_point_artifact(source, policy=policy)


def artifact_from_dict(data: Mapping[str, Any]) -> ModernPointArtifactV1:
    """Build a modern point artifact from JSON-compatible data."""

    return ModernPointArtifactV1.from_dict(data)


def read_modern_point_artifact(path: str | Path) -> ModernPointArtifactV1:
    """Read one modern per-point artifact from disk."""

    return artifact_from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


def write_modern_point_artifact(artifact: ModernPointArtifactV1, path: str | Path) -> Path:
    """Write one modern per-point artifact to disk."""

    destination = Path(path)
    destination.write_text(artifact.to_json(), encoding="utf-8")
    return destination


@dataclass(frozen=True)
class ModernPointBridgeCouplingsRecord:
    """Artifact-side serialized modern coupling bridge."""

    schema_id: str = MODERN_POINT_BRIDGE_COUPLINGS_SCHEMA_ID
    lane_id: str = MODERN_LANE_ID
    input_bundle_schema_id: str = ""
    input_bundle_id: str = ""
    input_provenance_id: str = ""
    input_resolution_policy_id: str = ""
    qcd_metadata_id: str = ""
    alpha_s_policy_id: str = ""
    operator_convention_id: str = ""
    ckm_target_id: str = ""
    quark_mass_target_id: str = ""
    point_id: str = ""
    point_label: str = ""
    Lambda_IR: float = 0.0
    M_KK: float = 0.0
    xi_KK: float = 0.0
    alpha_s: float = 0.0
    g_s: float = 0.0
    g_s_4d: float = 0.0
    g_eff: float = 0.0
    g_s_multiplier: float = 0.0
    coupling_policy_id: str = ""
    left_overlap: tuple[tuple[complex, ...], ...] = ()
    right_up_overlap: tuple[tuple[complex, ...], ...] = ()
    right_down_overlap: tuple[tuple[complex, ...], ...] = ()
    left_up: tuple[tuple[complex, ...], ...] = ()
    left_down: tuple[tuple[complex, ...], ...] = ()
    right_up: tuple[tuple[complex, ...], ...] = ()
    right_down: tuple[tuple[complex, ...], ...] = ()
    notes: str = (
        "Artifact-side serialization of the modern pointwise KK-gluon couplings."
    )

    def __post_init__(self) -> None:
        _require_system_id("schema_id", self.schema_id, MODERN_POINT_BRIDGE_COUPLINGS_SCHEMA_ID)
        _require_system_id("lane_id", self.lane_id, MODERN_LANE_ID)
        object.__setattr__(self, "input_bundle_schema_id", _require_text("input_bundle_schema_id", self.input_bundle_schema_id))
        object.__setattr__(self, "input_bundle_id", _require_text("input_bundle_id", self.input_bundle_id))
        object.__setattr__(self, "input_provenance_id", _require_text("input_provenance_id", self.input_provenance_id))
        object.__setattr__(self, "input_resolution_policy_id", _require_text("input_resolution_policy_id", self.input_resolution_policy_id))
        object.__setattr__(self, "qcd_metadata_id", _require_text("qcd_metadata_id", self.qcd_metadata_id))
        object.__setattr__(self, "alpha_s_policy_id", _require_text("alpha_s_policy_id", self.alpha_s_policy_id))
        object.__setattr__(self, "operator_convention_id", _require_text("operator_convention_id", self.operator_convention_id))
        object.__setattr__(self, "ckm_target_id", _require_text("ckm_target_id", self.ckm_target_id))
        object.__setattr__(self, "quark_mass_target_id", _require_text("quark_mass_target_id", self.quark_mass_target_id))
        object.__setattr__(self, "point_id", _require_text("point_id", self.point_id))
        object.__setattr__(self, "point_label", _require_text("point_label", self.point_label))
        object.__setattr__(self, "Lambda_IR", _require_positive_float("Lambda_IR", self.Lambda_IR))
        object.__setattr__(self, "M_KK", _require_positive_float("M_KK", self.M_KK))
        object.__setattr__(self, "xi_KK", _require_positive_float("xi_KK", self.xi_KK))
        object.__setattr__(self, "alpha_s", _require_positive_float("alpha_s", self.alpha_s))
        object.__setattr__(self, "g_s", _require_positive_float("g_s", self.g_s))
        object.__setattr__(self, "g_s_4d", _require_positive_float("g_s_4d", self.g_s_4d))
        object.__setattr__(self, "g_eff", _require_positive_float("g_eff", self.g_eff))
        object.__setattr__(
            self,
            "g_s_multiplier",
            _require_positive_float("g_s_multiplier", self.g_s_multiplier),
        )
        object.__setattr__(
            self,
            "coupling_policy_id",
            _require_text("coupling_policy_id", self.coupling_policy_id),
        )
        object.__setattr__(self, "left_overlap", _canonical_complex_matrix("left_overlap", self.left_overlap))
        object.__setattr__(self, "right_up_overlap", _canonical_complex_matrix("right_up_overlap", self.right_up_overlap))
        object.__setattr__(self, "right_down_overlap", _canonical_complex_matrix("right_down_overlap", self.right_down_overlap))
        object.__setattr__(self, "left_up", _canonical_complex_matrix("left_up", self.left_up))
        object.__setattr__(self, "left_down", _canonical_complex_matrix("left_down", self.left_down))
        object.__setattr__(self, "right_up", _canonical_complex_matrix("right_up", self.right_up))
        object.__setattr__(self, "right_down", _canonical_complex_matrix("right_down", self.right_down))
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "lane_id": self.lane_id,
            "input_bundle_schema_id": self.input_bundle_schema_id,
            "input_bundle_id": self.input_bundle_id,
            "input_provenance_id": self.input_provenance_id,
            "input_resolution_policy_id": self.input_resolution_policy_id,
            "qcd_metadata_id": self.qcd_metadata_id,
            "alpha_s_policy_id": self.alpha_s_policy_id,
            "operator_convention_id": self.operator_convention_id,
            "ckm_target_id": self.ckm_target_id,
            "quark_mass_target_id": self.quark_mass_target_id,
            "point_id": self.point_id,
            "point_label": self.point_label,
            "Lambda_IR": self.Lambda_IR,
            "M_KK": self.M_KK,
            "xi_KK": self.xi_KK,
            "alpha_s": self.alpha_s,
            "g_s": self.g_s,
            "g_s_4d": self.g_s_4d,
            "g_eff": self.g_eff,
            "g_s_multiplier": self.g_s_multiplier,
            "coupling_policy_id": self.coupling_policy_id,
            "left_overlap": _matrix_payload(self.left_overlap),
            "right_up_overlap": _matrix_payload(self.right_up_overlap),
            "right_down_overlap": _matrix_payload(self.right_down_overlap),
            "left_up": _matrix_payload(self.left_up),
            "left_down": _matrix_payload(self.left_down),
            "right_up": _matrix_payload(self.right_up),
            "right_down": _matrix_payload(self.right_down),
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ModernPointBridgeCouplingsRecord":
        mapping = _require_mapping(data, context="bridge.couplings")
        _require_exact_keys(
            mapping,
            expected=(
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
                "g_s_4d",
                "g_eff",
                "g_s_multiplier",
                "coupling_policy_id",
                "left_overlap",
                "right_up_overlap",
                "right_down_overlap",
                "left_up",
                "left_down",
                "right_up",
                "right_down",
                "notes",
            ),
            context="bridge.couplings",
        )
        return cls(
            schema_id=_require_text_value(mapping["schema_id"], context="bridge.couplings.schema_id"),
            lane_id=_require_text_value(mapping["lane_id"], context="bridge.couplings.lane_id"),
            input_bundle_schema_id=_require_text_value(mapping["input_bundle_schema_id"], context="bridge.couplings.input_bundle_schema_id"),
            input_bundle_id=_require_text_value(mapping["input_bundle_id"], context="bridge.couplings.input_bundle_id"),
            input_provenance_id=_require_text_value(mapping["input_provenance_id"], context="bridge.couplings.input_provenance_id"),
            input_resolution_policy_id=_require_text_value(mapping["input_resolution_policy_id"], context="bridge.couplings.input_resolution_policy_id"),
            qcd_metadata_id=_require_text_value(mapping["qcd_metadata_id"], context="bridge.couplings.qcd_metadata_id"),
            alpha_s_policy_id=_require_text_value(mapping["alpha_s_policy_id"], context="bridge.couplings.alpha_s_policy_id"),
            operator_convention_id=_require_text_value(mapping["operator_convention_id"], context="bridge.couplings.operator_convention_id"),
            ckm_target_id=_require_text_value(mapping["ckm_target_id"], context="bridge.couplings.ckm_target_id"),
            quark_mass_target_id=_require_text_value(mapping["quark_mass_target_id"], context="bridge.couplings.quark_mass_target_id"),
            point_id=_require_text_value(mapping["point_id"], context="bridge.couplings.point_id"),
            point_label=_require_text_value(mapping["point_label"], context="bridge.couplings.point_label"),
            Lambda_IR=_require_positive_float("bridge.couplings.Lambda_IR", mapping["Lambda_IR"]),
            M_KK=_require_positive_float("bridge.couplings.M_KK", mapping["M_KK"]),
            xi_KK=_require_positive_float("bridge.couplings.xi_KK", mapping["xi_KK"]),
            alpha_s=_require_positive_float("bridge.couplings.alpha_s", mapping["alpha_s"]),
            g_s=_require_positive_float("bridge.couplings.g_s", mapping["g_s"]),
            g_s_4d=_require_positive_float("bridge.couplings.g_s_4d", mapping["g_s_4d"]),
            g_eff=_require_positive_float("bridge.couplings.g_eff", mapping["g_eff"]),
            g_s_multiplier=_require_positive_float(
                "bridge.couplings.g_s_multiplier",
                mapping["g_s_multiplier"],
            ),
            coupling_policy_id=_require_text_value(
                mapping["coupling_policy_id"],
                context="bridge.couplings.coupling_policy_id",
            ),
            left_overlap=_matrix_from_payload("bridge.couplings.left_overlap", mapping["left_overlap"]),
            right_up_overlap=_matrix_from_payload("bridge.couplings.right_up_overlap", mapping["right_up_overlap"]),
            right_down_overlap=_matrix_from_payload("bridge.couplings.right_down_overlap", mapping["right_down_overlap"]),
            left_up=_matrix_from_payload("bridge.couplings.left_up", mapping["left_up"]),
            left_down=_matrix_from_payload("bridge.couplings.left_down", mapping["left_down"]),
            right_up=_matrix_from_payload("bridge.couplings.right_up", mapping["right_up"]),
            right_down=_matrix_from_payload("bridge.couplings.right_down", mapping["right_down"]),
            notes=_require_text_value(mapping["notes"], context="bridge.couplings.notes"),
        )

    @classmethod
    def from_source(cls, source: Any) -> "ModernPointBridgeCouplingsRecord":
        return cls.from_dict(_bridge_payload_from_source(source, context="bridge.couplings"))


@dataclass(frozen=True)
class ModernPointBridgeMatchingSystemRecord:
    """Artifact-side serialized matching record for one neutral-meson system."""

    schema_id: str = MODERN_POINT_BRIDGE_MATCHING_SYSTEM_SCHEMA_ID
    system_id: str = ""
    observable_id: str = ""
    backend_system_id: str = ""
    backend_key: str = ""
    display_name: str = ""
    sector_id: str = ""
    generations: tuple[int, int] = (0, 1)
    operator_basis_id: str = ""
    matching_scale_GeV: float = 0.0
    weight_policy_id: str = ""
    left_coupling: complex = 0.0j
    right_coupling: complex = 0.0j
    c1_vll: complex = 0.0j
    c1_vrr: complex = 0.0j
    c4_lr: complex = 0.0j
    c5_lr: complex = 0.0j
    note: str = ""

    def __post_init__(self) -> None:
        _require_system_id("schema_id", self.schema_id, MODERN_POINT_BRIDGE_MATCHING_SYSTEM_SCHEMA_ID)
        object.__setattr__(self, "system_id", _require_text("system_id", self.system_id))
        object.__setattr__(self, "observable_id", _require_text("observable_id", self.observable_id))
        object.__setattr__(self, "backend_system_id", _require_text("backend_system_id", self.backend_system_id))
        object.__setattr__(self, "backend_key", _require_text("backend_key", self.backend_key))
        object.__setattr__(self, "display_name", _require_text("display_name", self.display_name))
        object.__setattr__(self, "sector_id", _require_text("sector_id", self.sector_id))
        if self.system_id not in MODERN_POINT_BRIDGE_MATCHING_SYSTEM_IDS:
            raise ValueError(f"system_id must be one of {MODERN_POINT_BRIDGE_MATCHING_SYSTEM_IDS!r}")
        if not isinstance(self.generations, tuple) or len(self.generations) != 2:
            raise ValueError("generations must be a tuple of two indices")
        i, j = self.generations
        if not isinstance(i, int) or not isinstance(j, int) or i == j or i not in {0, 1, 2} or j not in {0, 1, 2}:
            raise ValueError("generations must be distinct integers in {0, 1, 2}")
        object.__setattr__(self, "operator_basis_id", _require_text("operator_basis_id", self.operator_basis_id))
        object.__setattr__(self, "matching_scale_GeV", _require_positive_float("matching_scale_GeV", self.matching_scale_GeV))
        object.__setattr__(self, "weight_policy_id", _require_text("weight_policy_id", self.weight_policy_id))
        object.__setattr__(self, "left_coupling", _require_complex("left_coupling", self.left_coupling))
        object.__setattr__(self, "right_coupling", _require_complex("right_coupling", self.right_coupling))
        object.__setattr__(self, "c1_vll", _require_complex("c1_vll", self.c1_vll))
        object.__setattr__(self, "c1_vrr", _require_complex("c1_vrr", self.c1_vrr))
        object.__setattr__(self, "c4_lr", _require_complex("c4_lr", self.c4_lr))
        object.__setattr__(self, "c5_lr", _require_complex("c5_lr", self.c5_lr))
        object.__setattr__(self, "note", _require_text("note", self.note))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "system_id": self.system_id,
            "observable_id": self.observable_id,
            "backend_system_id": self.backend_system_id,
            "backend_key": self.backend_key,
            "display_name": self.display_name,
            "sector_id": self.sector_id,
            "generations": list(self.generations),
            "operator_basis_id": self.operator_basis_id,
            "matching_scale_GeV": self.matching_scale_GeV,
            "weight_policy_id": self.weight_policy_id,
            "left_coupling": _complex_payload(self.left_coupling),
            "right_coupling": _complex_payload(self.right_coupling),
            "c1_vll": _complex_payload(self.c1_vll),
            "c1_vrr": _complex_payload(self.c1_vrr),
            "c4_lr": _complex_payload(self.c4_lr),
            "c5_lr": _complex_payload(self.c5_lr),
            "note": self.note,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ModernPointBridgeMatchingSystemRecord":
        mapping = _require_mapping(data, context="bridge.matching.system")
        _require_exact_keys(
            mapping,
            expected=(
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
            ),
            context="bridge.matching.system",
        )
        return cls(
            schema_id=_require_text_value(mapping["schema_id"], context="bridge.matching.system.schema_id"),
            system_id=_require_text_value(mapping["system_id"], context="bridge.matching.system.system_id"),
            observable_id=_require_text_value(mapping["observable_id"], context="bridge.matching.system.observable_id"),
            backend_system_id=_require_text_value(mapping["backend_system_id"], context="bridge.matching.system.backend_system_id"),
            backend_key=_require_text_value(mapping["backend_key"], context="bridge.matching.system.backend_key"),
            display_name=_require_text_value(mapping["display_name"], context="bridge.matching.system.display_name"),
            sector_id=_require_text_value(mapping["sector_id"], context="bridge.matching.system.sector_id"),
            generations=tuple(int(item) for item in _require_sequence(mapping["generations"], context="bridge.matching.system.generations")),
            operator_basis_id=_require_text_value(mapping["operator_basis_id"], context="bridge.matching.system.operator_basis_id"),
            matching_scale_GeV=_require_positive_float("bridge.matching.system.matching_scale_GeV", mapping["matching_scale_GeV"]),
            weight_policy_id=_require_text_value(mapping["weight_policy_id"], context="bridge.matching.system.weight_policy_id"),
            left_coupling=_complex_from_payload("bridge.matching.system.left_coupling", mapping["left_coupling"]),
            right_coupling=_complex_from_payload("bridge.matching.system.right_coupling", mapping["right_coupling"]),
            c1_vll=_complex_from_payload("bridge.matching.system.c1_vll", mapping["c1_vll"]),
            c1_vrr=_complex_from_payload("bridge.matching.system.c1_vrr", mapping["c1_vrr"]),
            c4_lr=_complex_from_payload("bridge.matching.system.c4_lr", mapping["c4_lr"]),
            c5_lr=_complex_from_payload("bridge.matching.system.c5_lr", mapping["c5_lr"]),
            note=_require_text_value(mapping["note"], context="bridge.matching.system.note"),
        )


@dataclass(frozen=True)
class ModernPointBridgeMatchingRecord:
    """Artifact-side serialized matching bridge."""

    schema_id: str = MODERN_POINT_BRIDGE_MATCHING_SCHEMA_ID
    lane_id: str = MODERN_LANE_ID
    couplings_schema_id: str = MODERN_POINT_BRIDGE_COUPLINGS_SCHEMA_ID
    input_bundle_schema_id: str = ""
    input_bundle_id: str = ""
    input_provenance_id: str = ""
    input_resolution_policy_id: str = ""
    qcd_metadata_id: str = ""
    alpha_s_policy_id: str = ""
    operator_basis_id: str = ""
    weight_policy_id: str = ""
    point_id: str = ""
    point_label: str = ""
    Lambda_IR: float = 0.0
    M_KK: float = 0.0
    xi_KK: float = 0.0
    alpha_s: float = 0.0
    g_s: float = 0.0
    system_matches: tuple[ModernPointBridgeMatchingSystemRecord, ...] = field(default_factory=tuple)
    notes: str = (
        "Artifact-side serialization of the modern pointwise matching bridge."
    )

    def __post_init__(self) -> None:
        _require_system_id("schema_id", self.schema_id, MODERN_POINT_BRIDGE_MATCHING_SCHEMA_ID)
        _require_system_id("lane_id", self.lane_id, MODERN_LANE_ID)
        _require_system_id("couplings_schema_id", self.couplings_schema_id, MODERN_POINT_BRIDGE_COUPLINGS_SCHEMA_ID)
        object.__setattr__(self, "input_bundle_schema_id", _require_text("input_bundle_schema_id", self.input_bundle_schema_id))
        object.__setattr__(self, "input_bundle_id", _require_text("input_bundle_id", self.input_bundle_id))
        object.__setattr__(self, "input_provenance_id", _require_text("input_provenance_id", self.input_provenance_id))
        object.__setattr__(self, "input_resolution_policy_id", _require_text("input_resolution_policy_id", self.input_resolution_policy_id))
        object.__setattr__(self, "qcd_metadata_id", _require_text("qcd_metadata_id", self.qcd_metadata_id))
        object.__setattr__(self, "alpha_s_policy_id", _require_text("alpha_s_policy_id", self.alpha_s_policy_id))
        object.__setattr__(self, "operator_basis_id", _require_text("operator_basis_id", self.operator_basis_id))
        object.__setattr__(self, "weight_policy_id", _require_text("weight_policy_id", self.weight_policy_id))
        object.__setattr__(self, "point_id", _require_text("point_id", self.point_id))
        object.__setattr__(self, "point_label", _require_text("point_label", self.point_label))
        object.__setattr__(self, "Lambda_IR", _require_positive_float("Lambda_IR", self.Lambda_IR))
        object.__setattr__(self, "M_KK", _require_positive_float("M_KK", self.M_KK))
        object.__setattr__(self, "xi_KK", _require_positive_float("xi_KK", self.xi_KK))
        object.__setattr__(self, "alpha_s", _require_positive_float("alpha_s", self.alpha_s))
        object.__setattr__(self, "g_s", _require_positive_float("g_s", self.g_s))
        if not isinstance(self.system_matches, tuple):
            raise ValueError("system_matches must be a tuple")
        normalized = tuple(self.system_matches)
        if len(normalized) != len(MODERN_POINT_BRIDGE_MATCHING_SYSTEM_IDS):
            raise ValueError("system_matches must contain the frozen modern system set")
        for item in normalized:
            if not isinstance(item, ModernPointBridgeMatchingSystemRecord):
                raise ValueError("system_matches must contain only ModernPointBridgeMatchingSystemRecord instances")
        if tuple(item.system_id for item in normalized) != MODERN_POINT_BRIDGE_MATCHING_SYSTEM_IDS:
            raise ValueError("system_matches must preserve the frozen modern system order")
        object.__setattr__(self, "system_matches", normalized)
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    @property
    def system_ids(self) -> tuple[str, ...]:
        return tuple(item.system_id for item in self.system_matches)

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "lane_id": self.lane_id,
            "couplings_schema_id": self.couplings_schema_id,
            "input_bundle_schema_id": self.input_bundle_schema_id,
            "input_bundle_id": self.input_bundle_id,
            "input_provenance_id": self.input_provenance_id,
            "input_resolution_policy_id": self.input_resolution_policy_id,
            "qcd_metadata_id": self.qcd_metadata_id,
            "alpha_s_policy_id": self.alpha_s_policy_id,
            "operator_basis_id": self.operator_basis_id,
            "weight_policy_id": self.weight_policy_id,
            "point_id": self.point_id,
            "point_label": self.point_label,
            "Lambda_IR": self.Lambda_IR,
            "M_KK": self.M_KK,
            "xi_KK": self.xi_KK,
            "alpha_s": self.alpha_s,
            "g_s": self.g_s,
            "system_matches": [item.as_dict() for item in self.system_matches],
            "notes": self.notes,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "ModernPointBridgeMatchingRecord":
        mapping = _require_mapping(data, context="bridge.matching")
        _require_exact_keys(
            mapping,
            expected=(
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
            ),
            context="bridge.matching",
        )
        return cls(
            schema_id=_require_text_value(mapping["schema_id"], context="bridge.matching.schema_id"),
            lane_id=_require_text_value(mapping["lane_id"], context="bridge.matching.lane_id"),
            couplings_schema_id=_require_text_value(mapping["couplings_schema_id"], context="bridge.matching.couplings_schema_id"),
            input_bundle_schema_id=_require_text_value(mapping["input_bundle_schema_id"], context="bridge.matching.input_bundle_schema_id"),
            input_bundle_id=_require_text_value(mapping["input_bundle_id"], context="bridge.matching.input_bundle_id"),
            input_provenance_id=_require_text_value(mapping["input_provenance_id"], context="bridge.matching.input_provenance_id"),
            input_resolution_policy_id=_require_text_value(mapping["input_resolution_policy_id"], context="bridge.matching.input_resolution_policy_id"),
            qcd_metadata_id=_require_text_value(mapping["qcd_metadata_id"], context="bridge.matching.qcd_metadata_id"),
            alpha_s_policy_id=_require_text_value(mapping["alpha_s_policy_id"], context="bridge.matching.alpha_s_policy_id"),
            operator_basis_id=_require_text_value(mapping["operator_basis_id"], context="bridge.matching.operator_basis_id"),
            weight_policy_id=_require_text_value(mapping["weight_policy_id"], context="bridge.matching.weight_policy_id"),
            point_id=_require_text_value(mapping["point_id"], context="bridge.matching.point_id"),
            point_label=_require_text_value(mapping["point_label"], context="bridge.matching.point_label"),
            Lambda_IR=_require_positive_float("bridge.matching.Lambda_IR", mapping["Lambda_IR"]),
            M_KK=_require_positive_float("bridge.matching.M_KK", mapping["M_KK"]),
            xi_KK=_require_positive_float("bridge.matching.xi_KK", mapping["xi_KK"]),
            alpha_s=_require_positive_float("bridge.matching.alpha_s", mapping["alpha_s"]),
            g_s=_require_positive_float("bridge.matching.g_s", mapping["g_s"]),
            system_matches=tuple(
                ModernPointBridgeMatchingSystemRecord.from_dict(
                    _require_mapping(item, context="bridge.matching.system_matches")
                )
                for item in _require_sequence(mapping["system_matches"], context="bridge.matching.system_matches")
            ),
            notes=_require_text_value(mapping["notes"], context="bridge.matching.notes"),
        )

    @classmethod
    def from_source(cls, source: Any) -> "ModernPointBridgeMatchingRecord":
        return cls.from_dict(_bridge_payload_from_source(source, context="bridge.matching"))


@dataclass(frozen=True)
class ModernPointBridgeArtifactV1:
    """Frozen sidecar artifact carrying the explicit modern QS2 bridge."""

    schema_id: str = MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID
    schema_version: int = MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION
    lane_id: str = MODERN_LANE_ID
    input_bundle_schema_id: str = ""
    input_bundle_id: str = ""
    input_provenance_id: str = ""
    input_resolution_policy_id: str = ""
    qcd_metadata_id: str = ""
    alpha_s_policy_id: str = ""
    point_id: str = ""
    point_label: str = ""
    coupling_schema_id: str = MODERN_POINT_BRIDGE_COUPLINGS_SCHEMA_ID
    matching_schema_id: str = MODERN_POINT_BRIDGE_MATCHING_SCHEMA_ID
    couplings: ModernPointBridgeCouplingsRecord = field(default_factory=ModernPointBridgeCouplingsRecord)
    matching: ModernPointBridgeMatchingRecord = field(default_factory=ModernPointBridgeMatchingRecord)
    notes: str = (
        "Explicit modern QS2 bridge sidecar only. This artifact exports the "
        "pointwise coupling and matching bridge with schema and lineage checks, "
        "but it does not claim a full EFT or RG package."
    )

    def __post_init__(self) -> None:
        _require_system_id("schema_id", self.schema_id, MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID)
        _require_int("schema_version", self.schema_version)
        if self.schema_version != MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION:
            raise ValueError(
                f"schema_version must be exactly {MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION}"
            )
        _require_system_id("lane_id", self.lane_id, MODERN_LANE_ID)
        object.__setattr__(self, "input_bundle_schema_id", _require_text("input_bundle_schema_id", self.input_bundle_schema_id))
        object.__setattr__(self, "input_bundle_id", _require_text("input_bundle_id", self.input_bundle_id))
        object.__setattr__(self, "input_provenance_id", _require_text("input_provenance_id", self.input_provenance_id))
        object.__setattr__(self, "input_resolution_policy_id", _require_text("input_resolution_policy_id", self.input_resolution_policy_id))
        object.__setattr__(self, "qcd_metadata_id", _require_text("qcd_metadata_id", self.qcd_metadata_id))
        object.__setattr__(self, "alpha_s_policy_id", _require_text("alpha_s_policy_id", self.alpha_s_policy_id))
        object.__setattr__(self, "point_id", _require_text("point_id", self.point_id))
        object.__setattr__(self, "point_label", _require_text("point_label", self.point_label))
        _require_system_id("coupling_schema_id", self.coupling_schema_id, MODERN_POINT_BRIDGE_COUPLINGS_SCHEMA_ID)
        _require_system_id("matching_schema_id", self.matching_schema_id, MODERN_POINT_BRIDGE_MATCHING_SCHEMA_ID)
        if not isinstance(self.couplings, ModernPointBridgeCouplingsRecord):
            raise ValueError("couplings must be a ModernPointBridgeCouplingsRecord instance")
        if not isinstance(self.matching, ModernPointBridgeMatchingRecord):
            raise ValueError("matching must be a ModernPointBridgeMatchingRecord instance")
        if self.couplings.schema_id != self.coupling_schema_id:
            raise ValueError("couplings.schema_id must match coupling_schema_id")
        if self.matching.schema_id != self.matching_schema_id:
            raise ValueError("matching.schema_id must match matching_schema_id")
        if self.couplings.point_id != self.point_id or self.matching.point_id != self.point_id:
            raise ValueError("bridge point_id must match nested couplings and matching")
        if self.couplings.point_label != self.point_label or self.matching.point_label != self.point_label:
            raise ValueError("bridge point_label must match nested couplings and matching")
        if self.couplings.input_bundle_schema_id != self.input_bundle_schema_id or self.matching.input_bundle_schema_id != self.input_bundle_schema_id:
            raise ValueError("bridge input_bundle_schema_id must match nested couplings and matching")
        if self.couplings.input_bundle_id != self.input_bundle_id or self.matching.input_bundle_id != self.input_bundle_id:
            raise ValueError("bridge input_bundle_id must match nested couplings and matching")
        if self.couplings.input_provenance_id != self.input_provenance_id or self.matching.input_provenance_id != self.input_provenance_id:
            raise ValueError("bridge input_provenance_id must match nested couplings and matching")
        if self.couplings.input_resolution_policy_id != self.input_resolution_policy_id or self.matching.input_resolution_policy_id != self.input_resolution_policy_id:
            raise ValueError("bridge input_resolution_policy_id must match nested couplings and matching")
        if self.couplings.qcd_metadata_id != self.qcd_metadata_id or self.matching.qcd_metadata_id != self.qcd_metadata_id:
            raise ValueError("bridge qcd_metadata_id must match nested couplings and matching")
        if self.couplings.alpha_s_policy_id != self.alpha_s_policy_id or self.matching.alpha_s_policy_id != self.alpha_s_policy_id:
            raise ValueError("bridge alpha_s_policy_id must match nested couplings and matching")
        if not math.isclose(self.couplings.M_KK, self.matching.M_KK, rel_tol=0.0, abs_tol=0.0):
            raise ValueError("couplings.M_KK and matching.M_KK must agree")
        if not math.isclose(self.couplings.xi_KK, self.matching.xi_KK, rel_tol=0.0, abs_tol=0.0):
            raise ValueError("couplings.xi_KK and matching.xi_KK must agree")
        if self.matching.couplings_schema_id != self.coupling_schema_id:
            raise ValueError("matching.couplings_schema_id must match coupling_schema_id")
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    @property
    def M_KK(self) -> float:
        return self.couplings.M_KK

    @property
    def xi_KK(self) -> float:
        return self.couplings.xi_KK

    def as_dict(self) -> dict[str, object]:
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
            "couplings": self.couplings.as_dict(),
            "matching": self.matching.as_dict(),
            "notes": self.notes,
        }

    def to_json(self) -> str:
        return _dump_json(self.as_dict())

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
            schema_id=_require_text_value(mapping["schema_id"], context="bridge_artifact.schema_id"),
            schema_version=_require_int("bridge_artifact.schema_version", mapping["schema_version"]),
            lane_id=_require_text_value(mapping["lane_id"], context="bridge_artifact.lane_id"),
            input_bundle_schema_id=_require_text_value(mapping["input_bundle_schema_id"], context="bridge_artifact.input_bundle_schema_id"),
            input_bundle_id=_require_text_value(mapping["input_bundle_id"], context="bridge_artifact.input_bundle_id"),
            input_provenance_id=_require_text_value(mapping["input_provenance_id"], context="bridge_artifact.input_provenance_id"),
            input_resolution_policy_id=_require_text_value(mapping["input_resolution_policy_id"], context="bridge_artifact.input_resolution_policy_id"),
            qcd_metadata_id=_require_text_value(mapping["qcd_metadata_id"], context="bridge_artifact.qcd_metadata_id"),
            alpha_s_policy_id=_require_text_value(mapping["alpha_s_policy_id"], context="bridge_artifact.alpha_s_policy_id"),
            point_id=_require_text_value(mapping["point_id"], context="bridge_artifact.point_id"),
            point_label=_require_text_value(mapping["point_label"], context="bridge_artifact.point_label"),
            coupling_schema_id=_require_text_value(mapping["coupling_schema_id"], context="bridge_artifact.coupling_schema_id"),
            matching_schema_id=_require_text_value(mapping["matching_schema_id"], context="bridge_artifact.matching_schema_id"),
            couplings=ModernPointBridgeCouplingsRecord.from_dict(
                _require_mapping(mapping["couplings"], context="bridge_artifact.couplings")
            ),
            matching=ModernPointBridgeMatchingRecord.from_dict(
                _require_mapping(mapping["matching"], context="bridge_artifact.matching")
            ),
            notes=_require_text_value(mapping["notes"], context="bridge_artifact.notes"),
        )

    @classmethod
    def from_json(cls, payload: str | bytes) -> "ModernPointBridgeArtifactV1":
        return cls.from_dict(json.loads(payload))

    @classmethod
    def from_source(cls, source: Any) -> "ModernPointBridgeArtifactV1":
        couplings = ModernPointBridgeCouplingsRecord.from_source(_source_attr(source, "couplings"))
        matching = ModernPointBridgeMatchingRecord.from_source(_source_attr(source, "matching"))
        point_id = _require_text_value(_source_attr(source, "point_id"), context="source.point_id")
        point_label = _require_text_value(_source_attr(source, "point_label"), context="source.point_label")
        input_bundle_schema_id = (
            _require_text_value(_source_attr(source, "input_bundle_schema_id"), context="source.input_bundle_schema_id")
            if hasattr(source, "input_bundle_schema_id")
            else couplings.input_bundle_schema_id
        )
        input_bundle_id = (
            _require_text_value(_source_attr(source, "input_bundle_id"), context="source.input_bundle_id")
            if hasattr(source, "input_bundle_id")
            else couplings.input_bundle_id
        )
        input_provenance_id = (
            _require_text_value(_source_attr(source, "input_provenance_id"), context="source.input_provenance_id")
            if hasattr(source, "input_provenance_id")
            else couplings.input_provenance_id
        )
        input_resolution_policy_id = (
            _require_text_value(_source_attr(source, "input_resolution_policy_id"), context="source.input_resolution_policy_id")
            if hasattr(source, "input_resolution_policy_id")
            else couplings.input_resolution_policy_id
        )
        return cls(
            input_bundle_schema_id=input_bundle_schema_id,
            input_bundle_id=input_bundle_id,
            input_provenance_id=input_provenance_id,
            input_resolution_policy_id=input_resolution_policy_id,
            qcd_metadata_id=couplings.qcd_metadata_id,
            alpha_s_policy_id=couplings.alpha_s_policy_id,
            point_id=point_id,
            point_label=point_label,
            couplings=couplings,
            matching=matching,
        )


def build_modern_point_bridge_artifact(source: Any) -> ModernPointBridgeArtifactV1:
    """Build one frozen modern bridge sidecar from a point-evaluation source."""

    return ModernPointBridgeArtifactV1.from_source(source)


def default_modern_point_bridge_artifact(source: Any) -> ModernPointBridgeArtifactV1:
    """Return the canonical frozen modern bridge sidecar."""

    return build_modern_point_bridge_artifact(source)


def bridge_artifact_from_dict(data: Mapping[str, Any]) -> ModernPointBridgeArtifactV1:
    """Build a modern bridge artifact from JSON-compatible data."""

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
    "MODERN_POINT_ARTIFACT_BACKEND_KEYS",
    "MODERN_POINT_ARTIFACT_BACKEND_SYSTEM_IDS",
    "MODERN_POINT_ARTIFACT_EVALUATION_SCHEMA_ID",
    "MODERN_POINT_ARTIFACT_HEADER_SCHEMA_ID",
    "MODERN_POINT_ARTIFACT_REQUIRED_OBSERVABLE_IDS",
    "MODERN_POINT_ARTIFACT_REQUIRED_POLICY_SYSTEM_IDS",
    "MODERN_POINT_ARTIFACT_REQUIRED_SYSTEM_IDS",
    "MODERN_POINT_ARTIFACT_SCHEMA_ID",
    "MODERN_POINT_ARTIFACT_SCHEMA_VERSION",
    "MODERN_POINT_ARTIFACT_VERDICT_SCHEMA_ID",
    "MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_ID",
    "MODERN_POINT_BRIDGE_ARTIFACT_SCHEMA_VERSION",
    "MODERN_POINT_BRIDGE_COUPLINGS_SCHEMA_ID",
    "MODERN_POINT_BRIDGE_MATCHING_SCHEMA_ID",
    "MODERN_POINT_BRIDGE_MATCHING_SYSTEM_IDS",
    "MODERN_POINT_BRIDGE_MATCHING_SYSTEM_SCHEMA_ID",
    "ModernPointArtifactHeader",
    "ModernPointArtifactV1",
    "ModernPointArtifactVerdict",
    "ModernPointBridgeArtifactV1",
    "ModernPointBridgeCouplingsRecord",
    "ModernPointBridgeMatchingRecord",
    "ModernPointBridgeMatchingSystemRecord",
    "artifact_from_dict",
    "bridge_artifact_from_dict",
    "build_modern_point_artifact",
    "build_modern_point_bridge_artifact",
    "default_modern_point_artifact",
    "default_modern_point_bridge_artifact",
    "read_modern_point_bridge_artifact",
    "read_modern_point_artifact",
    "write_modern_point_bridge_artifact",
    "write_modern_point_artifact",
]
