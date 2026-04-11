"""Structural scan contracts for the dedicated paper-facing 0710.1869 mode."""

from __future__ import annotations

from dataclasses import dataclass, field

from .conventions import (
    PAPER_0710_1869_MODE_ID,
    Paper07101869Conventions,
    default_paper_0710_1869_conventions,
)
from .scales import Paper07101869ScalePoint, default_paper_0710_1869_scales
from .validation import require_known_schema_id, require_member

PAPER_0710_1869_SCAN_SCHEMA_ID = "quarkConstraints.paper_0710_1869.scan.v1"


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
