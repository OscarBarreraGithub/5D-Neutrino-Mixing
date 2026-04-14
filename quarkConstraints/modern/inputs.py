"""Versioned input families for the modern quark lane."""

from __future__ import annotations

from dataclasses import dataclass, field

from .conventions import (
    MODERN_DEFAULT_BUNDLE_FAMILY_ID,
    MODERN_DEFAULT_BUNDLE_SLUG,
    MODERN_DEFAULT_INPUT_BUNDLE_ID,
    MODERN_DEFAULT_INPUT_PROVENANCE_ID,
    MODERN_INPUT_REGISTRY_SCHEMA_ID,
    MODERN_LANE_ID,
    MODERN_STRICT_PAPER_RESOLUTION_POLICY_ID,
    MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID,
    MODERN_STRICT_PAPER_BUNDLE_SLUG,
    MODERN_STRICT_PAPER_INPUT_BUNDLE_ID,
    MODERN_STRICT_PAPER_INPUT_PROVENANCE_ID,
    MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_ID,
    MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_NAMES,
    MODERN_STRICT_PAPER_SOURCE_SCHEMA_IDS,
    build_modern_bundle_id,
    build_modern_provenance_id,
    ModernLaneConventions,
    default_modern_lane_conventions,
)

MODERN_STRICT_PAPER_INPUTS_SCHEMA_ID = "quarkConstraints.modern.inputs.strict_paper_inputs.v1"
MODERN_DEFAULT_INPUTS_SCHEMA_ID = "quarkConstraints.modern.inputs.modern_default_inputs.v1"
MODERN_DEFAULT_RESOLUTION_POLICY_ID = (
    "quarkConstraints.modern.inputs.modern_default.explicit_repo_payloads.v1"
)
MODERN_DEFAULT_PROVENANCE_RECORD_SCHEMA_ID = (
    "quarkConstraints.modern.inputs.modern_default.provenance_record.v1"
)
MODERN_DEFAULT_NEUTRAL_MESON_INPUT_SCHEMA_ID = (
    "quarkConstraints.modern.inputs.modern_default.neutral_meson_input.v1"
)
MODERN_DEFAULT_OPERATOR_WEIGHT_POLICY_SCHEMA_ID = (
    "quarkConstraints.modern.inputs.modern_default.operator_weight_policy.v1"
)
MODERN_DEFAULT_CKM_TARGET_SCHEMA_ID = (
    "quarkConstraints.modern.inputs.modern_default.ckm_target.v1"
)
MODERN_DEFAULT_QUARK_MASS_TARGET_SCHEMA_ID = (
    "quarkConstraints.modern.inputs.modern_default.quark_mass_target.v1"
)
MODERN_DEFAULT_QCD_METADATA_SCHEMA_ID = (
    "quarkConstraints.modern.inputs.modern_default.qcd_metadata.v1"
)
MODERN_DEFAULT_NEUTRAL_MESON_SYSTEM_IDS = ("epsilon_K", "B_d", "B_s", "D0")
MODERN_DEFAULT_WEIGHT_POLICY_ID = (
    "quarkConstraints.modern.inputs.modern_default.operator_weights.repo_deltaf2_v1"
)
MODERN_DEFAULT_FIXED_SCALE_TARGET_ID = (
    "quarkConstraints.modern.inputs.modern_default.fixed_scale_targets.mu_3tev.v1"
)
MODERN_DEFAULT_QCD_METADATA_ID = (
    "quarkConstraints.modern.inputs.modern_default.qcd_reference.mu_3tev.v1"
)
MODERN_DEFAULT_ALPHA_S_POLICY_ID = (
    "qcd.alpha_s.high_precision.threshold_matched.reference_only.v1"
)
MODERN_DEFAULT_TARGET_SCALE_GEV = 3000.0
MODERN_DEFAULT_REFERENCE_ALPHA_S_3TEV = 0.0797
MODERN_DEFAULT_PROVENANCE_RECORD_IDS = (
    "quarkConstraints.modern.inputs.modern_default.source.fixed_scale_targets.v1",
    "quarkConstraints.modern.inputs.modern_default.source.deltaf2_inputs.v1",
    "quarkConstraints.modern.inputs.modern_default.source.scale_conventions.v1",
    "quarkConstraints.modern.inputs.modern_default.source.alpha_s_reference.v1",
)


def _require_text(name: str, value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _require_exact(name: str, value: str, *, expected: str) -> str:
    if value != expected:
        raise ValueError(f"{name} must be exactly {expected!r}")
    return value


def _require_member(name: str, value: str, allowed: tuple[str, ...]) -> str:
    if value not in allowed:
        raise ValueError(f"{name} must be one of {allowed!r}")
    return value


def _require_positive_float(name: str, value: float) -> float:
    numeric = float(value)
    if numeric <= 0.0:
        raise ValueError(f"{name} must be positive")
    return numeric


def _require_flavor_triplet(name: str, value: tuple[str, ...]) -> tuple[str, str, str]:
    if not isinstance(value, tuple):
        raise ValueError(f"{name} must be a tuple of three flavor ids")
    if len(value) != 3:
        raise ValueError(f"{name} must contain exactly three flavor ids")
    normalized = tuple(_require_text(name, item) for item in value)
    if len(set(normalized)) != 3:
        raise ValueError(f"{name} must not repeat flavor ids")
    return normalized


def _require_mass_triplet(name: str, value: tuple[float, ...]) -> tuple[float, float, float]:
    if not isinstance(value, tuple):
        raise ValueError(f"{name} must be a tuple of three positive masses")
    if len(value) != 3:
        raise ValueError(f"{name} must contain exactly three masses")
    return tuple(_require_positive_float(name, item) for item in value)


def _require_generations(name: str, value: tuple[int, int]) -> tuple[int, int]:
    if not isinstance(value, tuple):
        raise ValueError(f"{name} must be a tuple of two generation indices")
    if len(value) != 2:
        raise ValueError(f"{name} must contain exactly two generation indices")
    first, second = value
    if not isinstance(first, int) or not isinstance(second, int):
        raise ValueError(f"{name} must contain integer generation indices")
    if first == second or not (0 <= first < 3) or not (0 <= second < 3):
        raise ValueError(f"{name} must contain distinct generation indices in {{0, 1, 2}}")
    return first, second


def _signature_items(mapping: dict[str, object]) -> tuple[tuple[str, object], ...]:
    return tuple(sorted(mapping.items()))


def _signature_dict(items: tuple[tuple[str, object], ...]) -> dict[str, object]:
    return dict(items)


MODERN_STRICT_PAPER_TABLE_I_PAYLOAD_SIGNATURE = _signature_items(
    {
        "schema_id": "quarkConstraints.paper_0710_1869.table_i_inputs.v1",
        "mode_id": "paper_0710_1869",
        "sectors": (
            (
                "Q",
                (
                    (1, 0.64, 0.002),
                    (2, 0.59, 0.01),
                    (3, 0.46, 0.2),
                ),
            ),
            (
                "u",
                (
                    (1, 0.68, 0.0007),
                    (2, 0.53, 0.06),
                    (3, -0.06, 0.8),
                ),
            ),
            (
                "d",
                (
                    (1, 0.65, 0.002),
                    (2, 0.60, 0.008),
                    (3, 0.58, 0.02),
                ),
            ),
        ),
    }
)
MODERN_STRICT_PAPER_EQ3_PAYLOAD_SIGNATURE = _signature_items(
    {
        "schema_id": "quarkConstraints.paper_0710_1869.eq3_example.v1",
        "mode_id": "paper_0710_1869",
        "source_kind": "equation",
        "locator_label": "Eq. (3)",
        "detail": "example parameters (a, r, theta12, theta23, theta13, delta)",
        "citation": "arXiv:0710.1869, Eq. (3)",
        "a": 0.8,
        "r": 0.3,
        "theta12_deg": 115.0,
        "theta23_deg": 65.0,
        "theta13_deg": 70.0,
        "delta": 0.6,
    }
)
MODERN_STRICT_PAPER_CANONICAL_SNAPSHOTS = (
    MODERN_STRICT_PAPER_TABLE_I_PAYLOAD_SIGNATURE,
    MODERN_STRICT_PAPER_EQ3_PAYLOAD_SIGNATURE,
)


@dataclass(frozen=True)
class ModernStrictPaperSnapshot:
    """Frozen canonical snapshot for one strict-paper input block."""

    schema_id: str
    mode_id: str
    source_kind: str
    locator_label: str
    detail: str
    citation: str
    payload_signature: tuple[tuple[str, object], ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "schema_id", _require_text("schema_id", self.schema_id))
        object.__setattr__(self, "mode_id", _require_text("mode_id", self.mode_id))
        object.__setattr__(self, "source_kind", _require_text("source_kind", self.source_kind))
        object.__setattr__(self, "locator_label", _require_text("locator_label", self.locator_label))
        object.__setattr__(self, "detail", _require_text("detail", self.detail))
        object.__setattr__(self, "citation", _require_text("citation", self.citation))
        object.__setattr__(self, "payload_signature", tuple(self.payload_signature))
        if not self.payload_signature:
            raise ValueError("payload_signature must be non-empty")

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "source_kind": self.source_kind,
            "locator_label": self.locator_label,
            "detail": self.detail,
            "citation": self.citation,
            "payload_signature": _signature_dict(self.payload_signature),
        }


def _build_canonical_strict_paper_inputs() -> tuple[ModernStrictPaperSnapshot, ModernStrictPaperSnapshot]:
    return (
        ModernStrictPaperSnapshot(
            schema_id=MODERN_STRICT_PAPER_SOURCE_SCHEMA_IDS[0],
            mode_id="paper_0710_1869",
            source_kind="table",
            locator_label="Table I",
            detail="c_Q/f_Q eigenvalues",
            citation="arXiv:0710.1869, Table I",
            payload_signature=MODERN_STRICT_PAPER_TABLE_I_PAYLOAD_SIGNATURE,
        ),
        ModernStrictPaperSnapshot(
            schema_id=MODERN_STRICT_PAPER_SOURCE_SCHEMA_IDS[1],
            mode_id="paper_0710_1869",
            source_kind="equation",
            locator_label="Eq. (3)",
            detail="example parameters (a, r, theta12, theta23, theta13, delta)",
            citation="arXiv:0710.1869, Eq. (3)",
            payload_signature=MODERN_STRICT_PAPER_EQ3_PAYLOAD_SIGNATURE,
        ),
    )


@dataclass(frozen=True)
class ModernStrictPaperInputs:
    """Thin modern adapter over the existing paper-lane sourced inputs."""

    schema_id: str = MODERN_STRICT_PAPER_INPUTS_SCHEMA_ID
    lane_id: str = MODERN_LANE_ID
    family_id: str = MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID
    bundle_slug: str = MODERN_STRICT_PAPER_BUNDLE_SLUG
    provenance_slug: str = MODERN_STRICT_PAPER_BUNDLE_SLUG
    bundle_id: str = MODERN_STRICT_PAPER_INPUT_BUNDLE_ID
    provenance_id: str = MODERN_STRICT_PAPER_INPUT_PROVENANCE_ID
    source_lane_id: str = "paper_0710_1869"
    source_snapshot_id: str = MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_ID
    source_snapshot_names: tuple[str, str] = MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_NAMES
    source_schema_ids: tuple[str, str] = MODERN_STRICT_PAPER_SOURCE_SCHEMA_IDS
    source_resolution_policy_id: str = MODERN_STRICT_PAPER_RESOLUTION_POLICY_ID
    notes: str = (
        "Registry adapter only. This family stores resolver metadata, not eager "
        "paper payloads, and it does not define any modern physics."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact("schema_id", self.schema_id, expected=MODERN_STRICT_PAPER_INPUTS_SCHEMA_ID),
        )
        object.__setattr__(self, "lane_id", _require_exact("lane_id", self.lane_id, expected=MODERN_LANE_ID))
        object.__setattr__(
            self,
            "family_id",
            _require_exact("family_id", self.family_id, expected=MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID),
        )
        object.__setattr__(
            self,
            "bundle_slug",
            _require_exact("bundle_slug", self.bundle_slug, expected=MODERN_STRICT_PAPER_BUNDLE_SLUG),
        )
        object.__setattr__(
            self,
            "provenance_slug",
            _require_exact(
                "provenance_slug",
                self.provenance_slug,
                expected=MODERN_STRICT_PAPER_BUNDLE_SLUG,
            ),
        )
        object.__setattr__(
            self,
            "bundle_id",
            _require_exact("bundle_id", self.bundle_id, expected=MODERN_STRICT_PAPER_INPUT_BUNDLE_ID),
        )
        object.__setattr__(
            self,
            "provenance_id",
            _require_exact(
                "provenance_id",
                self.provenance_id,
                expected=MODERN_STRICT_PAPER_INPUT_PROVENANCE_ID,
            ),
        )
        object.__setattr__(
            self,
            "source_lane_id",
            _require_exact("source_lane_id", self.source_lane_id, expected="paper_0710_1869"),
        )
        object.__setattr__(
            self,
            "source_snapshot_id",
            _require_exact(
                "source_snapshot_id",
                self.source_snapshot_id,
                expected=MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_ID,
            ),
        )
        object.__setattr__(
            self,
            "source_snapshot_names",
            tuple(self.source_snapshot_names),
        )
        if self.source_snapshot_names != MODERN_STRICT_PAPER_SOURCE_SNAPSHOT_NAMES:
            raise ValueError("source_snapshot_names must remain the frozen paper snapshot ids")
        object.__setattr__(
            self,
            "source_schema_ids",
            tuple(self.source_schema_ids),
        )
        if self.source_schema_ids != MODERN_STRICT_PAPER_SOURCE_SCHEMA_IDS:
            raise ValueError("source_schema_ids must remain the frozen paper adapter schemas")
        object.__setattr__(
            self,
            "source_resolution_policy_id",
            _require_exact(
                "source_resolution_policy_id",
                self.source_resolution_policy_id,
                expected=MODERN_STRICT_PAPER_RESOLUTION_POLICY_ID,
            ),
        )
        object.__setattr__(self, "notes", _require_text("notes", self.notes))
        _require_member("lane_id", self.lane_id, (MODERN_LANE_ID,))

    def resolve_paper_inputs(self) -> tuple[object, object]:
        """Return the canonical strict-paper snapshot without importing paper runtime."""
        resolved = _build_canonical_strict_paper_inputs()
        self.validate_resolved_paper_inputs(resolved)
        return resolved

    def validate_resolved_paper_inputs(self, paper_inputs: tuple[object, ...]) -> tuple[object, ...]:
        """Reject aliased or mutated paper snapshots under the frozen IDs."""
        if not isinstance(paper_inputs, tuple):
            raise ValueError("paper_inputs must be a tuple of canonical paper snapshots")
        if len(paper_inputs) != 2:
            raise ValueError("paper_inputs must contain the canonical Table I and Eq. (3) snapshots")
        first, second = paper_inputs
        if not isinstance(first, ModernStrictPaperSnapshot):
            raise ValueError("paper_inputs[0] must be a ModernStrictPaperSnapshot instance")
        if not isinstance(second, ModernStrictPaperSnapshot):
            raise ValueError("paper_inputs[1] must be a ModernStrictPaperSnapshot instance")
        if first.schema_id != self.source_schema_ids[0]:
            raise ValueError("paper_inputs[0] must keep the canonical Table I schema_id")
        if second.schema_id != self.source_schema_ids[1]:
            raise ValueError("paper_inputs[1] must keep the canonical Eq. (3) schema_id")
        if first.as_dict() != _build_canonical_strict_paper_inputs()[0].as_dict():
            raise ValueError("paper_inputs[0] must match the canonical strict-paper Table I snapshot")
        if second.as_dict() != _build_canonical_strict_paper_inputs()[1].as_dict():
            raise ValueError("paper_inputs[1] must match the canonical strict-paper Eq. (3) snapshot")
        return paper_inputs

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "lane_id": self.lane_id,
            "family_id": self.family_id,
            "bundle_slug": self.bundle_slug,
            "provenance_slug": self.provenance_slug,
            "bundle_id": self.bundle_id,
            "provenance_id": self.provenance_id,
            "source_lane_id": self.source_lane_id,
            "source_snapshot_id": self.source_snapshot_id,
            "source_snapshot_names": list(self.source_snapshot_names),
            "source_schema_ids": list(self.source_schema_ids),
            "source_resolution_policy_id": self.source_resolution_policy_id,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class ModernDefaultProvenanceRecord:
    """Frozen provenance record for one explicit modern default payload source."""

    schema_id: str = MODERN_DEFAULT_PROVENANCE_RECORD_SCHEMA_ID
    record_id: str = ""
    source_kind: str = ""
    source_path: str = ""
    source_symbol: str = ""
    description: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact(
                "schema_id",
                self.schema_id,
                expected=MODERN_DEFAULT_PROVENANCE_RECORD_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "record_id", _require_text("record_id", self.record_id))
        object.__setattr__(self, "source_kind", _require_text("source_kind", self.source_kind))
        object.__setattr__(self, "source_path", _require_text("source_path", self.source_path))
        object.__setattr__(self, "source_symbol", _require_text("source_symbol", self.source_symbol))
        object.__setattr__(self, "description", _require_text("description", self.description))

    def as_dict(self) -> dict[str, str]:
        return {
            "schema_id": self.schema_id,
            "record_id": self.record_id,
            "source_kind": self.source_kind,
            "source_path": self.source_path,
            "source_symbol": self.source_symbol,
            "description": self.description,
        }


@dataclass(frozen=True)
class ModernDefaultNeutralMesonInput:
    """Explicit repo-owned neutral-meson target in the modern default family."""

    schema_id: str = MODERN_DEFAULT_NEUTRAL_MESON_INPUT_SCHEMA_ID
    system_id: str = ""
    observable_kind: str = ""
    backend_key: str = ""
    display_name: str = ""
    column_name: str = ""
    reject_reason: str = ""
    sector_id: str = ""
    generations: tuple[int, int] = (0, 1)
    bound: float = 0.0
    weight_policy_id: str = MODERN_DEFAULT_WEIGHT_POLICY_ID
    provenance_record_id: str = MODERN_DEFAULT_PROVENANCE_RECORD_IDS[1]
    note: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact(
                "schema_id",
                self.schema_id,
                expected=MODERN_DEFAULT_NEUTRAL_MESON_INPUT_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "system_id",
            _require_member(
                "system_id",
                _require_text("system_id", self.system_id),
                MODERN_DEFAULT_NEUTRAL_MESON_SYSTEM_IDS,
            ),
        )
        object.__setattr__(self, "observable_kind", _require_text("observable_kind", self.observable_kind))
        object.__setattr__(self, "backend_key", _require_text("backend_key", self.backend_key))
        object.__setattr__(self, "display_name", _require_text("display_name", self.display_name))
        object.__setattr__(self, "column_name", _require_text("column_name", self.column_name))
        object.__setattr__(self, "reject_reason", _require_text("reject_reason", self.reject_reason))
        object.__setattr__(
            self,
            "sector_id",
            _require_member("sector_id", _require_text("sector_id", self.sector_id), ("down", "up")),
        )
        object.__setattr__(self, "generations", _require_generations("generations", self.generations))
        object.__setattr__(self, "bound", _require_positive_float("bound", self.bound))
        object.__setattr__(
            self,
            "weight_policy_id",
            _require_exact(
                "weight_policy_id",
                self.weight_policy_id,
                expected=MODERN_DEFAULT_WEIGHT_POLICY_ID,
            ),
        )
        object.__setattr__(
            self,
            "provenance_record_id",
            _require_exact(
                "provenance_record_id",
                self.provenance_record_id,
                expected=MODERN_DEFAULT_PROVENANCE_RECORD_IDS[1],
            ),
        )
        object.__setattr__(self, "note", _require_text("note", self.note))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "system_id": self.system_id,
            "observable_kind": self.observable_kind,
            "backend_key": self.backend_key,
            "display_name": self.display_name,
            "column_name": self.column_name,
            "reject_reason": self.reject_reason,
            "sector_id": self.sector_id,
            "generations": list(self.generations),
            "bound": self.bound,
            "weight_policy_id": self.weight_policy_id,
            "provenance_record_id": self.provenance_record_id,
            "note": self.note,
        }


@dataclass(frozen=True)
class ModernDefaultOperatorWeightPolicy:
    """Shared operator-weight surrogate used by the current repo-owned bounds."""

    schema_id: str = MODERN_DEFAULT_OPERATOR_WEIGHT_POLICY_SCHEMA_ID
    policy_id: str = MODERN_DEFAULT_WEIGHT_POLICY_ID
    operator_convention_id: str = "kk_gluon_tree_v1"
    acceptance_semantics_id: str = "largest_weighted_operator_amplitude.v1"
    reference_scale_GeV: float = MODERN_DEFAULT_TARGET_SCALE_GEV
    ll_weight: float = 1.0
    rr_weight: float = 1.0
    lr1_weight: float = 7.0
    lr2_weight: float = 2.0
    provenance_record_id: str = MODERN_DEFAULT_PROVENANCE_RECORD_IDS[1]
    notes: str = (
        "Shared conservative operator-weight surrogate copied from the current "
        "repo-owned Delta F = 2 benchmark bundle at mu = 3 TeV."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact(
                "schema_id",
                self.schema_id,
                expected=MODERN_DEFAULT_OPERATOR_WEIGHT_POLICY_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "policy_id",
            _require_exact("policy_id", self.policy_id, expected=MODERN_DEFAULT_WEIGHT_POLICY_ID),
        )
        object.__setattr__(
            self,
            "operator_convention_id",
            _require_text("operator_convention_id", self.operator_convention_id),
        )
        object.__setattr__(
            self,
            "acceptance_semantics_id",
            _require_text("acceptance_semantics_id", self.acceptance_semantics_id),
        )
        object.__setattr__(
            self,
            "reference_scale_GeV",
            _require_positive_float("reference_scale_GeV", self.reference_scale_GeV),
        )
        object.__setattr__(self, "ll_weight", _require_positive_float("ll_weight", self.ll_weight))
        object.__setattr__(self, "rr_weight", _require_positive_float("rr_weight", self.rr_weight))
        object.__setattr__(self, "lr1_weight", _require_positive_float("lr1_weight", self.lr1_weight))
        object.__setattr__(self, "lr2_weight", _require_positive_float("lr2_weight", self.lr2_weight))
        object.__setattr__(
            self,
            "provenance_record_id",
            _require_exact(
                "provenance_record_id",
                self.provenance_record_id,
                expected=MODERN_DEFAULT_PROVENANCE_RECORD_IDS[1],
            ),
        )
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "policy_id": self.policy_id,
            "operator_convention_id": self.operator_convention_id,
            "acceptance_semantics_id": self.acceptance_semantics_id,
            "reference_scale_GeV": self.reference_scale_GeV,
            "ll_weight": self.ll_weight,
            "rr_weight": self.rr_weight,
            "lr1_weight": self.lr1_weight,
            "lr2_weight": self.lr2_weight,
            "provenance_record_id": self.provenance_record_id,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class ModernDefaultCKMTarget:
    """Frozen CKM target metadata for the repo-owned fixed-scale benchmark."""

    schema_id: str = MODERN_DEFAULT_CKM_TARGET_SCHEMA_ID
    target_id: str = MODERN_DEFAULT_FIXED_SCALE_TARGET_ID
    scale_GeV: float = MODERN_DEFAULT_TARGET_SCALE_GEV
    parameterization_id: str = "ckm_like_unitary.rotation_parameters.v1"
    theta12: float = 0.2274
    theta13: float = 0.00368
    theta23: float = 0.0415
    delta: float = 1.196
    provenance_record_id: str = MODERN_DEFAULT_PROVENANCE_RECORD_IDS[0]
    notes: str = (
        "Repo-owned CKM target metadata for the fixed-scale mu = 3 TeV fit "
        "bundle. Stores only the explicit angles and phase used to build the "
        "target unitary."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact(
                "schema_id",
                self.schema_id,
                expected=MODERN_DEFAULT_CKM_TARGET_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "target_id",
            _require_exact("target_id", self.target_id, expected=MODERN_DEFAULT_FIXED_SCALE_TARGET_ID),
        )
        object.__setattr__(self, "scale_GeV", _require_positive_float("scale_GeV", self.scale_GeV))
        object.__setattr__(self, "parameterization_id", _require_text("parameterization_id", self.parameterization_id))
        object.__setattr__(self, "theta12", float(self.theta12))
        object.__setattr__(self, "theta13", float(self.theta13))
        object.__setattr__(self, "theta23", float(self.theta23))
        object.__setattr__(self, "delta", float(self.delta))
        object.__setattr__(
            self,
            "provenance_record_id",
            _require_exact(
                "provenance_record_id",
                self.provenance_record_id,
                expected=MODERN_DEFAULT_PROVENANCE_RECORD_IDS[0],
            ),
        )
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "target_id": self.target_id,
            "scale_GeV": self.scale_GeV,
            "parameterization_id": self.parameterization_id,
            "theta12": self.theta12,
            "theta13": self.theta13,
            "theta23": self.theta23,
            "delta": self.delta,
            "provenance_record_id": self.provenance_record_id,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class ModernDefaultQuarkMassTarget:
    """Frozen quark-mass target metadata for the repo-owned fixed-scale fit."""

    schema_id: str = MODERN_DEFAULT_QUARK_MASS_TARGET_SCHEMA_ID
    target_id: str = MODERN_DEFAULT_FIXED_SCALE_TARGET_ID
    scale_GeV: float = MODERN_DEFAULT_TARGET_SCALE_GEV
    scheme_id: str = "repo_fixed_scale_running_masses.v1"
    up_flavor_ids: tuple[str, str, str] = ("u", "c", "t")
    up_masses_GeV: tuple[float, float, float] = (0.0013, 0.62, 172.0)
    down_flavor_ids: tuple[str, str, str] = ("d", "s", "b")
    down_masses_GeV: tuple[float, float, float] = (0.0028, 0.057, 2.86)
    provenance_record_id: str = MODERN_DEFAULT_PROVENANCE_RECORD_IDS[0]
    notes: str = (
        "Repo-owned fixed-scale quark mass targets at mu = 3 TeV. This block "
        "freezes only the explicit benchmark masses used by the current fit "
        "layer and does not imply a wider running-mass package."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact(
                "schema_id",
                self.schema_id,
                expected=MODERN_DEFAULT_QUARK_MASS_TARGET_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "target_id",
            _require_exact("target_id", self.target_id, expected=MODERN_DEFAULT_FIXED_SCALE_TARGET_ID),
        )
        object.__setattr__(self, "scale_GeV", _require_positive_float("scale_GeV", self.scale_GeV))
        object.__setattr__(self, "scheme_id", _require_text("scheme_id", self.scheme_id))
        object.__setattr__(
            self,
            "up_flavor_ids",
            _require_flavor_triplet("up_flavor_ids", tuple(self.up_flavor_ids)),
        )
        object.__setattr__(
            self,
            "up_masses_GeV",
            _require_mass_triplet("up_masses_GeV", tuple(self.up_masses_GeV)),
        )
        object.__setattr__(
            self,
            "down_flavor_ids",
            _require_flavor_triplet("down_flavor_ids", tuple(self.down_flavor_ids)),
        )
        object.__setattr__(
            self,
            "down_masses_GeV",
            _require_mass_triplet("down_masses_GeV", tuple(self.down_masses_GeV)),
        )
        object.__setattr__(
            self,
            "provenance_record_id",
            _require_exact(
                "provenance_record_id",
                self.provenance_record_id,
                expected=MODERN_DEFAULT_PROVENANCE_RECORD_IDS[0],
            ),
        )
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "target_id": self.target_id,
            "scale_GeV": self.scale_GeV,
            "scheme_id": self.scheme_id,
            "up_flavor_ids": list(self.up_flavor_ids),
            "up_masses_GeV": list(self.up_masses_GeV),
            "down_flavor_ids": list(self.down_flavor_ids),
            "down_masses_GeV": list(self.down_masses_GeV),
            "provenance_record_id": self.provenance_record_id,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class ModernDefaultQCDMetadata:
    """Reference-only QCD metadata for the current repo-owned benchmark family."""

    schema_id: str = MODERN_DEFAULT_QCD_METADATA_SCHEMA_ID
    metadata_id: str = MODERN_DEFAULT_QCD_METADATA_ID
    target_scale_GeV: float = MODERN_DEFAULT_TARGET_SCALE_GEV
    matching_scale_GeV: float = MODERN_DEFAULT_TARGET_SCALE_GEV
    xi_KK: float = 1.0
    alpha_s_reference_scale_GeV: float = MODERN_DEFAULT_TARGET_SCALE_GEV
    alpha_s_reference_value: float = MODERN_DEFAULT_REFERENCE_ALPHA_S_3TEV
    alpha_s_policy_id: str = MODERN_DEFAULT_ALPHA_S_POLICY_ID
    alpha_s_precision: str = "high"
    scale_convention_provenance_record_id: str = MODERN_DEFAULT_PROVENANCE_RECORD_IDS[2]
    alpha_s_provenance_record_id: str = MODERN_DEFAULT_PROVENANCE_RECORD_IDS[3]
    notes: str = (
        "Reference-only QCD metadata for the repo-owned mu = 3 TeV benchmark. "
        "The alpha_s value is a frozen cross-check number, not a claim of a "
        "full RG or EFT package."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact(
                "schema_id",
                self.schema_id,
                expected=MODERN_DEFAULT_QCD_METADATA_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "metadata_id",
            _require_exact("metadata_id", self.metadata_id, expected=MODERN_DEFAULT_QCD_METADATA_ID),
        )
        object.__setattr__(
            self,
            "target_scale_GeV",
            _require_positive_float("target_scale_GeV", self.target_scale_GeV),
        )
        object.__setattr__(
            self,
            "matching_scale_GeV",
            _require_positive_float("matching_scale_GeV", self.matching_scale_GeV),
        )
        object.__setattr__(self, "xi_KK", _require_positive_float("xi_KK", self.xi_KK))
        object.__setattr__(
            self,
            "alpha_s_reference_scale_GeV",
            _require_positive_float("alpha_s_reference_scale_GeV", self.alpha_s_reference_scale_GeV),
        )
        object.__setattr__(
            self,
            "alpha_s_reference_value",
            _require_positive_float("alpha_s_reference_value", self.alpha_s_reference_value),
        )
        object.__setattr__(self, "alpha_s_policy_id", _require_text("alpha_s_policy_id", self.alpha_s_policy_id))
        object.__setattr__(self, "alpha_s_precision", _require_text("alpha_s_precision", self.alpha_s_precision))
        object.__setattr__(
            self,
            "scale_convention_provenance_record_id",
            _require_exact(
                "scale_convention_provenance_record_id",
                self.scale_convention_provenance_record_id,
                expected=MODERN_DEFAULT_PROVENANCE_RECORD_IDS[2],
            ),
        )
        object.__setattr__(
            self,
            "alpha_s_provenance_record_id",
            _require_exact(
                "alpha_s_provenance_record_id",
                self.alpha_s_provenance_record_id,
                expected=MODERN_DEFAULT_PROVENANCE_RECORD_IDS[3],
            ),
        )
        object.__setattr__(self, "notes", _require_text("notes", self.notes))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "metadata_id": self.metadata_id,
            "target_scale_GeV": self.target_scale_GeV,
            "matching_scale_GeV": self.matching_scale_GeV,
            "xi_KK": self.xi_KK,
            "alpha_s_reference_scale_GeV": self.alpha_s_reference_scale_GeV,
            "alpha_s_reference_value": self.alpha_s_reference_value,
            "alpha_s_policy_id": self.alpha_s_policy_id,
            "alpha_s_precision": self.alpha_s_precision,
            "scale_convention_provenance_record_id": self.scale_convention_provenance_record_id,
            "alpha_s_provenance_record_id": self.alpha_s_provenance_record_id,
            "notes": self.notes,
        }


def default_modern_default_provenance_records() -> tuple[ModernDefaultProvenanceRecord, ...]:
    """Return the frozen provenance records for the modern default bundle."""
    return (
        ModernDefaultProvenanceRecord(
            record_id=MODERN_DEFAULT_PROVENANCE_RECORD_IDS[0],
            source_kind="repo_owned_numeric_bundle",
            source_path="quarkConstraints/benchmarks.py",
            source_symbol="_FIXED_SCALE_TARGETS_MU_3TEV_V1",
            description=(
                "Repo-owned fixed-scale mu = 3 TeV target bundle used for the "
                "modern default CKM and quark-mass metadata."
            ),
        ),
        ModernDefaultProvenanceRecord(
            record_id=MODERN_DEFAULT_PROVENANCE_RECORD_IDS[1],
            source_kind="repo_owned_numeric_bundle",
            source_path="quarkConstraints/deltaf2.py",
            source_symbol="DEFAULT_DELTA_F2_INPUTS_V1",
            description=(
                "Repo-owned Delta F = 2 v1 surrogate bundle used here only as "
                "explicit modern bound and operator-weight payload data."
            ),
        ),
        ModernDefaultProvenanceRecord(
            record_id=MODERN_DEFAULT_PROVENANCE_RECORD_IDS[2],
            source_kind="repo_owned_scale_convention",
            source_path="quarkConstraints/scales.py",
            source_symbol="DEFAULT_QUARK_TARGET_SCALE_GEV, DEFAULT_QUARK_XI_KK",
            description=(
                "Repo-owned target-scale and KK-mass-convention metadata used by "
                "the current quark fit and scan machinery."
            ),
        ),
        ModernDefaultProvenanceRecord(
            record_id=MODERN_DEFAULT_PROVENANCE_RECORD_IDS[3],
            source_kind="repo_owned_reference_cross_check",
            source_path="tests/test_alpha_s.py",
            source_symbol="alpha_s(3000 GeV) ~= 0.0797",
            description=(
                "Frozen alpha_s reference value at mu = 3 TeV used only as "
                "modern QCD metadata, not as a general-purpose RG contract."
            ),
        ),
    )


def default_modern_default_neutral_meson_inputs() -> tuple[ModernDefaultNeutralMesonInput, ...]:
    """Return the explicit repo-owned neutral-meson targets for the default bundle."""
    return (
        ModernDefaultNeutralMesonInput(
            system_id="epsilon_K",
            observable_kind="cp_violation_surrogate",
            backend_key="epsilon_k",
            display_name="epsilon_K",
            column_name="epsilon_k_ratio",
            reject_reason="epsilon_k",
            sector_id="down",
            generations=(0, 1),
            bound=2.0e-8,
            note="Repo-owned kaon CP-violating surrogate bound from the v1 Delta F = 2 slice.",
        ),
        ModernDefaultNeutralMesonInput(
            system_id="B_d",
            observable_kind="mixing_amplitude_surrogate",
            backend_key="b_d",
            display_name="B_d mixing",
            column_name="b_d_ratio",
            reject_reason="b_d_mix",
            sector_id="down",
            generations=(0, 2),
            bound=4.0e-7,
            note="Repo-owned B_d neutral-meson surrogate bound at the 3 TeV matching scale.",
        ),
        ModernDefaultNeutralMesonInput(
            system_id="B_s",
            observable_kind="mixing_amplitude_surrogate",
            backend_key="b_s",
            display_name="B_s mixing",
            column_name="b_s_ratio",
            reject_reason="b_s_mix",
            sector_id="down",
            generations=(1, 2),
            bound=5.5e-6,
            note="Repo-owned B_s neutral-meson surrogate bound at the 3 TeV matching scale.",
        ),
        ModernDefaultNeutralMesonInput(
            system_id="D0",
            observable_kind="mixing_amplitude_surrogate",
            backend_key="d",
            display_name="D mixing",
            column_name="d_ratio",
            reject_reason="d_mix",
            sector_id="up",
            generations=(0, 1),
            bound=8.5e-9,
            note="Repo-owned conservative D0 surrogate bound from the v1 Delta F = 2 slice.",
        ),
    )


@dataclass(frozen=True)
class ModernDefaultInputs:
    """Explicit modern default payload family for the current repo-owned benchmark slice."""

    schema_id: str = MODERN_DEFAULT_INPUTS_SCHEMA_ID
    lane_id: str = MODERN_LANE_ID
    family_id: str = MODERN_DEFAULT_BUNDLE_FAMILY_ID
    bundle_slug: str = MODERN_DEFAULT_BUNDLE_SLUG
    provenance_slug: str = MODERN_DEFAULT_BUNDLE_SLUG
    bundle_id: str = MODERN_DEFAULT_INPUT_BUNDLE_ID
    provenance_id: str = MODERN_DEFAULT_INPUT_PROVENANCE_ID
    source_lane_id: str = MODERN_LANE_ID
    source_resolution_policy_id: str = MODERN_DEFAULT_RESOLUTION_POLICY_ID
    neutral_meson_inputs: tuple[ModernDefaultNeutralMesonInput, ...] = field(
        default_factory=default_modern_default_neutral_meson_inputs
    )
    operator_weight_policy: ModernDefaultOperatorWeightPolicy = field(
        default_factory=ModernDefaultOperatorWeightPolicy
    )
    ckm_target: ModernDefaultCKMTarget = field(default_factory=ModernDefaultCKMTarget)
    quark_mass_target: ModernDefaultQuarkMassTarget = field(default_factory=ModernDefaultQuarkMassTarget)
    qcd_metadata: ModernDefaultQCDMetadata = field(default_factory=ModernDefaultQCDMetadata)
    provenance_records: tuple[ModernDefaultProvenanceRecord, ...] = field(
        default_factory=default_modern_default_provenance_records
    )
    paper_inputs: tuple[object, ...] = ()
    notes: str = (
        "Explicit modern numeric payloads copied from the current repo-owned "
        "benchmark layer. This family remains disjoint from the paper lane and "
        "does not claim a full EFT, RG, or hadronic package."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact("schema_id", self.schema_id, expected=MODERN_DEFAULT_INPUTS_SCHEMA_ID),
        )
        object.__setattr__(self, "lane_id", _require_exact("lane_id", self.lane_id, expected=MODERN_LANE_ID))
        object.__setattr__(
            self,
            "family_id",
            _require_exact("family_id", self.family_id, expected=MODERN_DEFAULT_BUNDLE_FAMILY_ID),
        )
        object.__setattr__(
            self,
            "bundle_slug",
            _require_exact("bundle_slug", self.bundle_slug, expected=MODERN_DEFAULT_BUNDLE_SLUG),
        )
        object.__setattr__(
            self,
            "provenance_slug",
            _require_exact("provenance_slug", self.provenance_slug, expected=MODERN_DEFAULT_BUNDLE_SLUG),
        )
        object.__setattr__(
            self,
            "bundle_id",
            _require_exact("bundle_id", self.bundle_id, expected=MODERN_DEFAULT_INPUT_BUNDLE_ID),
        )
        object.__setattr__(
            self,
            "provenance_id",
            _require_exact("provenance_id", self.provenance_id, expected=MODERN_DEFAULT_INPUT_PROVENANCE_ID),
        )
        object.__setattr__(
            self,
            "source_lane_id",
            _require_exact("source_lane_id", self.source_lane_id, expected=MODERN_LANE_ID),
        )
        object.__setattr__(
            self,
            "source_resolution_policy_id",
            _require_exact(
                "source_resolution_policy_id",
                self.source_resolution_policy_id,
                expected=MODERN_DEFAULT_RESOLUTION_POLICY_ID,
            ),
        )
        object.__setattr__(self, "notes", _require_text("notes", self.notes))
        if not isinstance(self.neutral_meson_inputs, tuple):
            raise ValueError("neutral_meson_inputs must be a tuple")
        normalized_system_inputs = tuple(self.neutral_meson_inputs)
        if len(normalized_system_inputs) != len(MODERN_DEFAULT_NEUTRAL_MESON_SYSTEM_IDS):
            raise ValueError("neutral_meson_inputs must contain the frozen default system set")
        for item in normalized_system_inputs:
            if not isinstance(item, ModernDefaultNeutralMesonInput):
                raise ValueError(
                    "neutral_meson_inputs must contain only ModernDefaultNeutralMesonInput instances"
                )
        if tuple(item.system_id for item in normalized_system_inputs) != MODERN_DEFAULT_NEUTRAL_MESON_SYSTEM_IDS:
            raise ValueError("neutral_meson_inputs must preserve the frozen system order")
        canonical_system_inputs = default_modern_default_neutral_meson_inputs()
        if tuple(item.as_dict() for item in normalized_system_inputs) != tuple(
            item.as_dict() for item in canonical_system_inputs
        ):
            raise ValueError("neutral_meson_inputs must match the frozen default payload")
        object.__setattr__(self, "neutral_meson_inputs", normalized_system_inputs)
        if not isinstance(self.operator_weight_policy, ModernDefaultOperatorWeightPolicy):
            raise ValueError(
                "operator_weight_policy must be a ModernDefaultOperatorWeightPolicy instance"
            )
        if self.operator_weight_policy.as_dict() != ModernDefaultOperatorWeightPolicy().as_dict():
            raise ValueError("operator_weight_policy must match the frozen default payload")
        if not isinstance(self.ckm_target, ModernDefaultCKMTarget):
            raise ValueError("ckm_target must be a ModernDefaultCKMTarget instance")
        if self.ckm_target.as_dict() != ModernDefaultCKMTarget().as_dict():
            raise ValueError("ckm_target must match the frozen default payload")
        if not isinstance(self.quark_mass_target, ModernDefaultQuarkMassTarget):
            raise ValueError("quark_mass_target must be a ModernDefaultQuarkMassTarget instance")
        if self.quark_mass_target.as_dict() != ModernDefaultQuarkMassTarget().as_dict():
            raise ValueError("quark_mass_target must match the frozen default payload")
        if not isinstance(self.qcd_metadata, ModernDefaultQCDMetadata):
            raise ValueError("qcd_metadata must be a ModernDefaultQCDMetadata instance")
        if self.qcd_metadata.as_dict() != ModernDefaultQCDMetadata().as_dict():
            raise ValueError("qcd_metadata must match the frozen default payload")
        if not isinstance(self.provenance_records, tuple):
            raise ValueError("provenance_records must be a tuple")
        normalized_provenance_records = tuple(self.provenance_records)
        if len(normalized_provenance_records) != len(MODERN_DEFAULT_PROVENANCE_RECORD_IDS):
            raise ValueError("provenance_records must contain the frozen provenance set")
        for record in normalized_provenance_records:
            if not isinstance(record, ModernDefaultProvenanceRecord):
                raise ValueError(
                    "provenance_records must contain only ModernDefaultProvenanceRecord instances"
                )
        if tuple(record.record_id for record in normalized_provenance_records) != MODERN_DEFAULT_PROVENANCE_RECORD_IDS:
            raise ValueError("provenance_records must preserve the frozen provenance order")
        if tuple(record.as_dict() for record in normalized_provenance_records) != tuple(
            record.as_dict() for record in default_modern_default_provenance_records()
        ):
            raise ValueError("provenance_records must match the frozen default payload")
        object.__setattr__(self, "provenance_records", normalized_provenance_records)
        if self.operator_weight_policy.reference_scale_GeV != self.qcd_metadata.matching_scale_GeV:
            raise ValueError(
                "operator_weight_policy.reference_scale_GeV must match qcd_metadata.matching_scale_GeV"
            )
        if self.ckm_target.scale_GeV != self.quark_mass_target.scale_GeV:
            raise ValueError("ckm_target.scale_GeV must match quark_mass_target.scale_GeV")
        if self.ckm_target.scale_GeV != self.qcd_metadata.target_scale_GeV:
            raise ValueError("ckm_target.scale_GeV must match qcd_metadata.target_scale_GeV")
        provenance_record_ids = {record.record_id for record in normalized_provenance_records}
        if self.operator_weight_policy.provenance_record_id not in provenance_record_ids:
            raise ValueError("operator_weight_policy.provenance_record_id must resolve inside provenance_records")
        if self.ckm_target.provenance_record_id not in provenance_record_ids:
            raise ValueError("ckm_target.provenance_record_id must resolve inside provenance_records")
        if self.quark_mass_target.provenance_record_id not in provenance_record_ids:
            raise ValueError("quark_mass_target.provenance_record_id must resolve inside provenance_records")
        if self.qcd_metadata.scale_convention_provenance_record_id not in provenance_record_ids:
            raise ValueError(
                "qcd_metadata.scale_convention_provenance_record_id must resolve inside provenance_records"
            )
        if self.qcd_metadata.alpha_s_provenance_record_id not in provenance_record_ids:
            raise ValueError(
                "qcd_metadata.alpha_s_provenance_record_id must resolve inside provenance_records"
            )
        for item in normalized_system_inputs:
            if item.weight_policy_id != self.operator_weight_policy.policy_id:
                raise ValueError(
                    "neutral_meson_inputs weight_policy_id must match operator_weight_policy.policy_id"
                )
            if item.provenance_record_id not in provenance_record_ids:
                raise ValueError(
                    "neutral_meson_inputs provenance_record_id must resolve inside provenance_records"
                )
        if not isinstance(self.paper_inputs, tuple):
            raise ValueError("paper_inputs must be a tuple")
        object.__setattr__(self, "paper_inputs", tuple(self.paper_inputs))
        if self.paper_inputs:
            raise ValueError("paper_inputs must remain empty for modern_default_inputs")
        _require_member("lane_id", self.lane_id, (MODERN_LANE_ID,))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "lane_id": self.lane_id,
            "family_id": self.family_id,
            "bundle_slug": self.bundle_slug,
            "provenance_slug": self.provenance_slug,
            "bundle_id": self.bundle_id,
            "provenance_id": self.provenance_id,
            "source_lane_id": self.source_lane_id,
            "source_resolution_policy_id": self.source_resolution_policy_id,
            "neutral_meson_inputs": [item.as_dict() for item in self.neutral_meson_inputs],
            "operator_weight_policy": self.operator_weight_policy.as_dict(),
            "ckm_target": self.ckm_target.as_dict(),
            "quark_mass_target": self.quark_mass_target.as_dict(),
            "qcd_metadata": self.qcd_metadata.as_dict(),
            "provenance_records": [record.as_dict() for record in self.provenance_records],
            "paper_inputs": [input_bundle.as_dict() for input_bundle in self.paper_inputs],
            "notes": self.notes,
        }


@dataclass(frozen=True)
class ModernInputRegistry:
    """Closed modern registry containing the two disjoint bundle families."""

    schema_id: str = MODERN_INPUT_REGISTRY_SCHEMA_ID
    lane_id: str = MODERN_LANE_ID
    conventions: ModernLaneConventions = field(default_factory=default_modern_lane_conventions)
    strict_paper_inputs: ModernStrictPaperInputs = field(default_factory=ModernStrictPaperInputs)
    modern_default_inputs: ModernDefaultInputs = field(default_factory=ModernDefaultInputs)
    notes: str = (
        "Modern registry fence only. It freezes explicit versioned input "
        "payloads and provenance for the current repo-owned benchmark slice, "
        "but it still does not define couplings, EFT matching, RG, observables, "
        "scan engine, caching, warm-starts, artifacts, verifier, or parameter "
        "search."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            _require_exact("schema_id", self.schema_id, expected=MODERN_INPUT_REGISTRY_SCHEMA_ID),
        )
        object.__setattr__(self, "lane_id", _require_exact("lane_id", self.lane_id, expected=MODERN_LANE_ID))
        object.__setattr__(self, "notes", _require_text("notes", self.notes))
        if not isinstance(self.conventions, ModernLaneConventions):
            raise ValueError("conventions must be a ModernLaneConventions instance")
        if not isinstance(self.strict_paper_inputs, ModernStrictPaperInputs):
            raise ValueError("strict_paper_inputs must be a ModernStrictPaperInputs instance")
        if not isinstance(self.modern_default_inputs, ModernDefaultInputs):
            raise ValueError("modern_default_inputs must be a ModernDefaultInputs instance")
        if self.conventions.lane_id != self.lane_id:
            raise ValueError("conventions.lane_id must match registry lane_id")
        if self.strict_paper_inputs.lane_id != self.lane_id:
            raise ValueError("strict_paper_inputs.lane_id must match registry lane_id")
        if self.modern_default_inputs.lane_id != self.lane_id:
            raise ValueError("modern_default_inputs.lane_id must match registry lane_id")
        if self.strict_paper_inputs.family_id == self.modern_default_inputs.family_id:
            raise ValueError("bundle families must remain disjoint")
        if self.strict_paper_inputs.bundle_id == self.modern_default_inputs.bundle_id:
            raise ValueError("bundle ids must remain disjoint")
        if self.strict_paper_inputs.provenance_id == self.modern_default_inputs.provenance_id:
            raise ValueError("provenance ids must remain disjoint")
        if self.strict_paper_inputs.source_lane_id != "paper_0710_1869":
            raise ValueError("strict_paper_inputs must adapt the paper lane")
        if self.modern_default_inputs.source_lane_id != MODERN_LANE_ID:
            raise ValueError("modern_default_inputs must remain a modern namespace")
        if self.strict_paper_inputs.bundle_id != build_modern_bundle_id(
            MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID,
            MODERN_STRICT_PAPER_BUNDLE_SLUG,
        ):
            raise ValueError("strict_paper_inputs.bundle_id must follow the modern naming rule")
        if self.strict_paper_inputs.provenance_id != build_modern_provenance_id(
            MODERN_STRICT_PAPER_BUNDLE_FAMILY_ID,
            MODERN_STRICT_PAPER_BUNDLE_SLUG,
        ):
            raise ValueError(
                "strict_paper_inputs.provenance_id must follow the modern naming rule"
            )
        if self.modern_default_inputs.bundle_id != build_modern_bundle_id(
            MODERN_DEFAULT_BUNDLE_FAMILY_ID,
            MODERN_DEFAULT_BUNDLE_SLUG,
        ):
            raise ValueError("modern_default_inputs.bundle_id must follow the modern naming rule")
        if self.modern_default_inputs.provenance_id != build_modern_provenance_id(
            MODERN_DEFAULT_BUNDLE_FAMILY_ID,
            MODERN_DEFAULT_BUNDLE_SLUG,
        ):
            raise ValueError(
                "modern_default_inputs.provenance_id must follow the modern naming rule"
            )

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "lane_id": self.lane_id,
            "conventions": self.conventions.as_dict(),
            "strict_paper_inputs": self.strict_paper_inputs.as_dict(),
            "modern_default_inputs": self.modern_default_inputs.as_dict(),
            "notes": self.notes,
        }


def default_modern_strict_paper_inputs() -> ModernStrictPaperInputs:
    """Return the strict-paper modern adapter bundle."""
    return ModernStrictPaperInputs()


def default_modern_default_inputs() -> ModernDefaultInputs:
    """Return the explicit modern default numeric input bundle."""
    return ModernDefaultInputs()


def default_modern_input_registry() -> ModernInputRegistry:
    """Return the frozen modern input registry."""
    return ModernInputRegistry()
