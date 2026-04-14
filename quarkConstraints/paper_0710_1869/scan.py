"""Structural scan contracts for the dedicated paper-facing 0710.1869 mode."""

from __future__ import annotations

from dataclasses import dataclass, field

from .conventions import (
    PAPER_0710_1869_CONVENTIONS_SCHEMA_ID,
    PAPER_0710_1869_MODE_ID,
    PAPER_0710_1869_PAPER_ID,
    Paper07101869Conventions,
    default_paper_0710_1869_conventions,
)
from .scales import (
    PAPER_0710_1869_SCALES_SCHEMA_ID,
    Paper07101869ScalePoint,
    default_paper_0710_1869_scales,
)
from .validation import require_known_schema_id, require_member, require_nonempty_identifier

PAPER_0710_1869_SCAN_SCHEMA_ID = "quarkConstraints.paper_0710_1869.scan.v1"
PAPER_0710_1869_STRICT_PAPER_CLAIM_LEVEL_ID = "strict_paper"
PAPER_0710_1869_STRICT_PAPER_SCAN_STATUS = "strict_paper"
PAPER_0710_1869_STRICT_PAPER_INPUT_BUNDLE_ID = (
    "pr1.table_i_eq3_example.strict_paper.input_bundle.v1"
)
PAPER_0710_1869_STRICT_PAPER_INPUT_PROVENANCE_ID = (
    "pr1.table_i_eq3_example.strict_paper.input_provenance.v1"
)
PAPER_0710_1869_STRICT_PAPER_BENCHMARK_ID = "pr1.table_i_eq3_example.v1"


@dataclass(frozen=True)
class Paper07101869ScanRequest:
    """Canonical scan request for the paper-only package."""

    schema_id: str = PAPER_0710_1869_SCAN_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    conventions: Paper07101869Conventions = field(
        default_factory=default_paper_0710_1869_conventions
    )
    scale_points: tuple[Paper07101869ScalePoint, ...] = field(
        default_factory=lambda: (default_paper_0710_1869_scales(),)
    )

    def __post_init__(self) -> None:
        require_known_schema_id(
            "schema_id",
            self.schema_id,
            expected=PAPER_0710_1869_SCAN_SCHEMA_ID,
        )
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        if not isinstance(self.conventions, Paper07101869Conventions):
            raise ValueError("conventions must be a Paper07101869Conventions instance")

        points = tuple(self.scale_points)
        if not points:
            raise ValueError("scale_points must contain at least one scale point")

        labels = [point.label for point in points]
        if len(labels) != len(set(labels)):
            raise ValueError("scale_points labels must be unique within one scan request")

        for point in points:
            if not isinstance(point, Paper07101869ScalePoint):
                raise ValueError("scale_points must contain only Paper07101869ScalePoint items")
            if point.mode_id != self.mode_id:
                raise ValueError("scale point mode_id must match the scan request mode_id")

        if self.conventions.mode_id != self.mode_id:
            raise ValueError("conventions mode_id must match the scan request mode_id")

        object.__setattr__(self, "scale_points", points)


@dataclass(frozen=True)
class Paper07101869ScanRow:
    """Deterministic structural row emitted by the paper-only scan contract."""

    point_id: str
    mode_id: str
    conventions_schema_id: str
    scales_schema_id: str
    Lambda_IR_GeV: float
    m_g1_GeV: float
    mu_match_GeV: float
    mu_gs_GeV: float
    m_KK_eff_GeV: float | None
    propagator_mass_GeV: float
    status: str = "structural_only"
    note: str = (
        "Paper-mode structural row only. Observable evaluation is intentionally "
        "not implemented in this package yet."
    )


@dataclass(frozen=True)
class Paper07101869StrictPaperScanRow:
    """Frozen strict-paper scan row emitted beside the structural compatibility path."""

    point_id: str
    mode_id: str
    paper_id: str
    claim_level_id: str
    conventions_schema_id: str
    scales_schema_id: str
    provenance_policy_id: str
    verifier_policy_id: str
    seed_to_profile_mapping_policy_id: str
    universal_term_coefficient_policy_id: str
    profile_derivation_policy_id: str
    benchmark_id: str
    input_bundle_id: str
    input_provenance_id: str
    Lambda_IR_GeV: float
    m_g1_GeV: float
    mu_match_GeV: float
    mu_gs_GeV: float
    m_KK_eff_GeV: float | None
    propagator_mass_GeV: float
    status: str = PAPER_0710_1869_STRICT_PAPER_SCAN_STATUS
    note: str = (
        "Strict-paper reproduction scan contract only. This row freezes the published "
        "benchmark/input provenance boundary for the dedicated paper-mode scan surface."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "point_id",
            require_nonempty_identifier("point_id", self.point_id),
        )
        object.__setattr__(
            self,
            "mode_id",
            require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,)),
        )
        object.__setattr__(
            self,
            "paper_id",
            require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,)),
        )
        object.__setattr__(
            self,
            "claim_level_id",
            require_member(
                "claim_level_id",
                self.claim_level_id,
                (PAPER_0710_1869_STRICT_PAPER_CLAIM_LEVEL_ID,),
            ),
        )
        object.__setattr__(
            self,
            "conventions_schema_id",
            require_known_schema_id(
                "conventions_schema_id",
                self.conventions_schema_id,
                expected=PAPER_0710_1869_CONVENTIONS_SCHEMA_ID,
            ),
        )
        object.__setattr__(
            self,
            "scales_schema_id",
            require_known_schema_id(
                "scales_schema_id",
                self.scales_schema_id,
                expected=PAPER_0710_1869_SCALES_SCHEMA_ID,
            ),
        )
        for field_name in (
            "provenance_policy_id",
            "verifier_policy_id",
            "seed_to_profile_mapping_policy_id",
            "universal_term_coefficient_policy_id",
            "profile_derivation_policy_id",
            "note",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        object.__setattr__(
            self,
            "benchmark_id",
            require_member(
                "benchmark_id",
                self.benchmark_id,
                (PAPER_0710_1869_STRICT_PAPER_BENCHMARK_ID,),
            ),
        )
        object.__setattr__(
            self,
            "input_bundle_id",
            require_member(
                "input_bundle_id",
                self.input_bundle_id,
                (PAPER_0710_1869_STRICT_PAPER_INPUT_BUNDLE_ID,),
            ),
        )
        object.__setattr__(
            self,
            "input_provenance_id",
            require_member(
                "input_provenance_id",
                self.input_provenance_id,
                (PAPER_0710_1869_STRICT_PAPER_INPUT_PROVENANCE_ID,),
            ),
        )
        object.__setattr__(
            self,
            "status",
            require_member(
                "status",
                self.status,
                (PAPER_0710_1869_STRICT_PAPER_SCAN_STATUS,),
            ),
        )


def build_structural_scan_rows(
    request: Paper07101869ScanRequest,
) -> tuple[Paper07101869ScanRow, ...]:
    """Emit deterministic scan rows without routing through repo-v1 observables."""
    rows = []
    for point in request.scale_points:
        rows.append(
            Paper07101869ScanRow(
                point_id=point.label,
                mode_id=request.mode_id,
                conventions_schema_id=request.conventions.schema_id,
                scales_schema_id=point.schema_id,
                Lambda_IR_GeV=point.Lambda_IR_GeV,
                m_g1_GeV=point.m_g1_GeV,
                mu_match_GeV=point.mu_match_GeV,
                mu_gs_GeV=point.mu_gs_GeV,
                m_KK_eff_GeV=point.m_KK_eff_GeV,
                propagator_mass_GeV=point.propagator_mass_GeV,
            )
        )
    return tuple(rows)


def build_strict_paper_scan_rows(
    request: Paper07101869ScanRequest,
    *,
    benchmark_id: str = PAPER_0710_1869_STRICT_PAPER_BENCHMARK_ID,
    input_bundle_id: str = PAPER_0710_1869_STRICT_PAPER_INPUT_BUNDLE_ID,
    input_provenance_id: str = PAPER_0710_1869_STRICT_PAPER_INPUT_PROVENANCE_ID,
) -> tuple[Paper07101869StrictPaperScanRow, ...]:
    """Emit frozen strict-paper scan rows beside the structural compatibility surface."""
    resolved_benchmark_id = require_member(
        "benchmark_id",
        benchmark_id,
        (PAPER_0710_1869_STRICT_PAPER_BENCHMARK_ID,),
    )
    resolved_input_bundle_id = require_member(
        "input_bundle_id",
        input_bundle_id,
        (PAPER_0710_1869_STRICT_PAPER_INPUT_BUNDLE_ID,),
    )
    resolved_input_provenance_id = require_member(
        "input_provenance_id",
        input_provenance_id,
        (PAPER_0710_1869_STRICT_PAPER_INPUT_PROVENANCE_ID,),
    )
    rows = []
    for point in request.scale_points:
        rows.append(
            Paper07101869StrictPaperScanRow(
                point_id=point.label,
                mode_id=request.mode_id,
                paper_id=request.conventions.paper_id,
                claim_level_id=PAPER_0710_1869_STRICT_PAPER_CLAIM_LEVEL_ID,
                conventions_schema_id=request.conventions.schema_id,
                scales_schema_id=point.schema_id,
                provenance_policy_id=request.conventions.provenance_policy_id,
                verifier_policy_id=request.conventions.verifier_policy_id,
                seed_to_profile_mapping_policy_id=(
                    request.conventions.seed_to_profile_mapping_policy_id
                ),
                universal_term_coefficient_policy_id=(
                    request.conventions.universal_term_coefficient_policy_id
                ),
                profile_derivation_policy_id=request.conventions.profile_derivation_policy_id,
                benchmark_id=resolved_benchmark_id,
                input_bundle_id=resolved_input_bundle_id,
                input_provenance_id=resolved_input_provenance_id,
                Lambda_IR_GeV=point.Lambda_IR_GeV,
                m_g1_GeV=point.m_g1_GeV,
                mu_match_GeV=point.mu_match_GeV,
                mu_gs_GeV=point.mu_gs_GeV,
                m_KK_eff_GeV=point.m_KK_eff_GeV,
                propagator_mass_GeV=point.propagator_mass_GeV,
            )
        )
    return tuple(rows)
