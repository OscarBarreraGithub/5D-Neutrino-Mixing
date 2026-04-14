"""Paper-owned structural benchmarks and custom seeded reference builders."""

from __future__ import annotations

import importlib
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .conventions import PAPER_0710_1869_MODE_ID, PAPER_0710_1869_PAPER_ID
from .inputs import (
    Paper07101869Eq3Example,
    Paper07101869PhysicalSeedToProfileContract,
    Paper07101869TableIInputs,
    default_paper_0710_1869_eq3_example,
    default_paper_0710_1869_physical_seed_to_profile_contract,
    default_paper_0710_1869_table_i_inputs,
)
from .validation import (
    require_known_schema_id,
    require_member,
    require_nonempty_identifier,
    require_positive_finite,
)

PAPER_0710_1869_BENCHMARK_SCHEMA_ID = "quarkConstraints.paper_0710_1869.benchmarks.v1"
_PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID = (
    "profile_inputs.table_i_reference.diagnostics_only.v1"
)


def _as_positive_triplet(
    name: str, values: tuple[float, float, float] | list[float] | np.ndarray
) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must contain only positive values")
    return arr.astype(float, copy=True)


def _paper_model_module() -> Any:
    return importlib.import_module(".model", __package__)


def _rotation_parameters_type() -> type[Any]:
    return getattr(_paper_model_module(), "Paper07101869RotationParameters")


def _default_rotation_parameters() -> Any:
    return _rotation_parameters_type()()


def _merge_metadata_with_frozen_provenance(
    *,
    canonical_metadata: dict[str, Any],
    metadata: Mapping[str, Any] | None,
) -> dict[str, Any]:
    if metadata is None:
        return canonical_metadata
    if not isinstance(metadata, Mapping):
        raise ValueError("metadata must be a mapping")

    resolved_metadata = dict(canonical_metadata)
    for key, value in dict(metadata).items():
        if key in canonical_metadata and value != canonical_metadata[key]:
            raise ValueError(
                f"metadata[{key!r}] must not override the frozen canonical provenance value"
            )
        resolved_metadata[key] = canonical_metadata.get(key, value)
    return resolved_metadata


@dataclass(frozen=True)
class Paper07101869Benchmark:
    """Closed structural benchmark definition sourced from the paper text."""

    schema_id: str = PAPER_0710_1869_BENCHMARK_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    benchmark_id: str = "pr1.table_i_eq3_example.v1"
    label: str = "table_i_eq3_example"
    extraction_version: int = 1
    table_i_inputs: Paper07101869TableIInputs = field(
        default_factory=default_paper_0710_1869_table_i_inputs
    )
    eq3_example: Paper07101869Eq3Example = field(
        default_factory=default_paper_0710_1869_eq3_example
    )
    status: str = "sourced_structural_only"
    notes: str = (
        "Frozen paper-owned benchmark inputs only. "
        "No fit, EFT matching, RG running, or observable logic is implemented here."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_BENCHMARK_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        object.__setattr__(self, "paper_id", require_nonempty_identifier("paper_id", self.paper_id))
        object.__setattr__(
            self,
            "benchmark_id",
            require_nonempty_identifier("benchmark_id", self.benchmark_id),
        )
        object.__setattr__(self, "label", require_nonempty_identifier("label", self.label))
        if not isinstance(self.extraction_version, int) or isinstance(
            self.extraction_version, bool
        ):
            raise ValueError("extraction_version must be an integer")
        if self.extraction_version < 1:
            raise ValueError("extraction_version must be positive")
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        require_member("status", self.status, ("sourced_structural_only",))
        if not isinstance(self.table_i_inputs, Paper07101869TableIInputs):
            raise ValueError("table_i_inputs must be a Paper07101869TableIInputs")
        if not isinstance(self.eq3_example, Paper07101869Eq3Example):
            raise ValueError("eq3_example must be a Paper07101869Eq3Example")
        object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "benchmark_id": self.benchmark_id,
            "label": self.label,
            "extraction_version": self.extraction_version,
            "status": self.status,
            "notes": self.notes,
            "table_i_inputs": self.table_i_inputs.as_dict(),
            "eq3_example": self.eq3_example.as_dict(),
        }


@dataclass(frozen=True)
class Paper07101869BenchmarkSpurionSeed:
    """Explicit Yukawa seed required to build a custom structural-reference point."""

    up_singular_values: tuple[float, float, float]
    down_singular_values: tuple[float, float, float]
    overall_scale: float
    up_left: Any = field(default_factory=_default_rotation_parameters)
    up_right: Any = field(default_factory=_default_rotation_parameters)
    down_left: Any = field(default_factory=_default_rotation_parameters)
    down_right: Any = field(default_factory=_default_rotation_parameters)
    notes: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "up_singular_values",
            tuple(_as_positive_triplet("up_singular_values", self.up_singular_values)),
        )
        object.__setattr__(
            self,
            "down_singular_values",
            tuple(_as_positive_triplet("down_singular_values", self.down_singular_values)),
        )
        object.__setattr__(
            self,
            "overall_scale",
            require_positive_finite("overall_scale", self.overall_scale),
        )
        rotation_type = _rotation_parameters_type()
        for name in ("up_left", "up_right", "down_left", "down_right"):
            if not isinstance(getattr(self, name), rotation_type):
                raise ValueError(f"{name} must be a Paper07101869RotationParameters")
        if self.notes is not None:
            object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))


def default_paper_0710_1869_pr1_benchmark() -> Paper07101869Benchmark:
    """Return the PR1 structural benchmark assembled from frozen paper inputs."""
    return Paper07101869Benchmark()


def paper_0710_1869_pr1_benchmarks() -> tuple[Paper07101869Benchmark, ...]:
    """Return the closed set of sourced PR1 benchmark definitions."""
    return (default_paper_0710_1869_pr1_benchmark(),)


def build_paper_0710_1869_benchmark_point(
    benchmark: Paper07101869Benchmark,
    seed: Paper07101869BenchmarkSpurionSeed,
    *,
    profile_input_policy_id: str = _PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID,
    label: str | None = None,
    metadata: Mapping[str, Any] | None = None,
) -> Any:
    """Compatibility wrapper for a custom seeded structural-reference point."""
    return build_paper_0710_1869_seeded_structural_reference_point(
        benchmark,
        seed,
        profile_input_policy_id=profile_input_policy_id,
        label=label,
        metadata=metadata,
    )


def build_paper_0710_1869_seeded_structural_reference_point(
    benchmark: Paper07101869Benchmark,
    seed: Paper07101869BenchmarkSpurionSeed,
    *,
    profile_input_policy_id: str = _PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID,
    label: str | None = None,
    metadata: Mapping[str, Any] | None = None,
) -> Any:
    """Build a custom seeded point from a structural benchmark plus explicit spurion seed."""
    if not isinstance(benchmark, Paper07101869Benchmark):
        raise ValueError("benchmark must be a Paper07101869Benchmark")
    if not isinstance(seed, Paper07101869BenchmarkSpurionSeed):
        raise ValueError("seed must be a Paper07101869BenchmarkSpurionSeed")

    model = _paper_model_module()
    if profile_input_policy_id != getattr(
        model,
        "PAPER_0710_1869_REFERENCE_ONLY_PROFILE_INPUT_POLICY_ID",
    ):
        raise ValueError(
            "seeded structural-reference points only support the diagnostics-only "
            "Table I reference profile policy"
        )
    table_benchmark = model.build_table_i_benchmark_from_inputs(benchmark.table_i_inputs)
    resolved_metadata = _merge_metadata_with_frozen_provenance(
        canonical_metadata={
        "benchmark_id": benchmark.benchmark_id,
        "benchmark_status": benchmark.status,
        "benchmark_label": benchmark.label,
        "benchmark_reference_kind": "structural_only",
        "benchmark_eq3_example": {
            "a": benchmark.eq3_example.a,
            "r": benchmark.eq3_example.r,
            "theta12_deg": benchmark.eq3_example.theta12_deg,
            "theta23_deg": benchmark.eq3_example.theta23_deg,
            "theta13_deg": benchmark.eq3_example.theta13_deg,
            "delta": benchmark.eq3_example.delta,
        },
        },
        metadata=metadata,
    )

    note_parts = [
        (
            "Structural benchmark inputs supply quoted Table I sectors and the Eq. (3) "
            "reference bundle only. The returned object is a custom seeded structural "
            "reference point, not a sourced exact paper benchmark. Yukawa singular values "
            "and both left- and right-handed flavor rotations come from the explicit spurion seed."
        )
    ]
    if seed.notes is not None:
        note_parts.append(seed.notes)

    return model.build_paper_0710_1869_point_from_singular_values(
        up_singular_values=seed.up_singular_values,
        down_singular_values=seed.down_singular_values,
        overall_scale=seed.overall_scale,
        r=benchmark.eq3_example.r,
        up_left=seed.up_left,
        up_right=seed.up_right,
        down_left=seed.down_left,
        down_right=seed.down_right,
        q_sector=table_benchmark.q_sector,
        u_sector=table_benchmark.u_sector,
        d_sector=table_benchmark.d_sector,
        profile_input_policy_id=profile_input_policy_id,
        label=label or f"{benchmark.label}_seeded_structural_reference",
        construction_id=model.PAPER_0710_1869_BENCHMARK_SEED_CONSTRUCTION_ID,
        metadata=resolved_metadata,
        notes=" ".join(note_parts),
    )


def build_paper_0710_1869_seeded_physical_point(
    benchmark: Paper07101869Benchmark,
    seed: Paper07101869BenchmarkSpurionSeed,
    *,
    physical_contract: Paper07101869PhysicalSeedToProfileContract | None = None,
    label: str | None = None,
    metadata: Mapping[str, Any] | None = None,
) -> Any:
    """Build a hard-gated QS1 physical point from a structural benchmark plus seed."""
    if not isinstance(benchmark, Paper07101869Benchmark):
        raise ValueError("benchmark must be a Paper07101869Benchmark")
    if not isinstance(seed, Paper07101869BenchmarkSpurionSeed):
        raise ValueError("seed must be a Paper07101869BenchmarkSpurionSeed")

    resolved_contract = (
        default_paper_0710_1869_physical_seed_to_profile_contract()
        if physical_contract is None
        else physical_contract
    )
    if not isinstance(resolved_contract, Paper07101869PhysicalSeedToProfileContract):
        raise ValueError(
            "physical_contract must be a Paper07101869PhysicalSeedToProfileContract"
        )

    model = _paper_model_module()
    resolved_metadata = _merge_metadata_with_frozen_provenance(
        canonical_metadata={
        "benchmark_id": benchmark.benchmark_id,
        "benchmark_status": benchmark.status,
        "benchmark_label": benchmark.label,
        "benchmark_reference_kind": "physical_qs1_seed_to_profile",
        "physical_contract_schema_id": resolved_contract.schema_id,
        "physical_mapping_policy_id": resolved_contract.mapping_policy.policy_id,
        "physical_universal_term_policy_id": resolved_contract.universal_term_policy.policy_id,
        "benchmark_eq3_example": {
            "a": benchmark.eq3_example.a,
            "r": benchmark.eq3_example.r,
            "theta12_deg": benchmark.eq3_example.theta12_deg,
            "theta23_deg": benchmark.eq3_example.theta23_deg,
            "theta13_deg": benchmark.eq3_example.theta13_deg,
            "delta": benchmark.eq3_example.delta,
        },
        },
        metadata=metadata,
    )

    note_parts = [
        (
            "Structural benchmark inputs supply the paper-owned Eq. (3) reference bundle "
            "and benchmark provenance only. The returned object is a custom seeded QS1 "
            "physical seed-to-profile point. Its c and F values must be derived from "
            "the seeded MFV eigensystems under the exact frozen affine-per-sector "
            "contract, never from quoted Table I alias attachment."
        )
    ]
    if seed.notes is not None:
        note_parts.append(seed.notes)

    return model.build_paper_0710_1869_physical_point_from_singular_values(
        up_singular_values=seed.up_singular_values,
        down_singular_values=seed.down_singular_values,
        physical_contract=resolved_contract,
        overall_scale=seed.overall_scale,
        r=benchmark.eq3_example.r,
        up_left=seed.up_left,
        up_right=seed.up_right,
        down_left=seed.down_left,
        down_right=seed.down_right,
        label=label or f"{benchmark.label}_seeded_physical_qs1",
        construction_id=model.PAPER_0710_1869_PHYSICAL_BENCHMARK_SEED_CONSTRUCTION_ID,
        metadata=resolved_metadata,
        notes=" ".join(note_parts),
    )


__all__ = [
    "PAPER_0710_1869_BENCHMARK_SCHEMA_ID",
    "Paper07101869Benchmark",
    "Paper07101869BenchmarkSpurionSeed",
    "build_paper_0710_1869_benchmark_point",
    "build_paper_0710_1869_seeded_physical_point",
    "build_paper_0710_1869_seeded_structural_reference_point",
    "default_paper_0710_1869_pr1_benchmark",
    "paper_0710_1869_pr1_benchmarks",
]
