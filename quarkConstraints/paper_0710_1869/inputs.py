"""Sourced benchmark inputs owned by the dedicated ``paper_0710_1869`` mode."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .conventions import PAPER_0710_1869_MODE_ID, PAPER_0710_1869_PAPER_ID
from .validation import (
    require_known_schema_id,
    require_member,
    require_nonempty_identifier,
    require_positive_finite,
)

PAPER_0710_1869_SOURCE_SCHEMA_ID = "quarkConstraints.paper_0710_1869.source_ref.v1"
PAPER_0710_1869_TABLE_I_SCHEMA_ID = "quarkConstraints.paper_0710_1869.table_i_inputs.v1"
PAPER_0710_1869_EQ3_SCHEMA_ID = "quarkConstraints.paper_0710_1869.eq3_example.v1"


def _require_finite_float(name: str, value: float) -> float:
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ValueError(f"{name} must be a finite float")
    return numeric


@dataclass(frozen=True)
class Paper07101869SourceRef:
    """Frozen citation metadata for one sourced paper value block."""

    schema_id: str = PAPER_0710_1869_SOURCE_SCHEMA_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    source_kind: str = "table"
    locator_label: str = "Table I"
    detail: str = "c_Q/f_Q eigenvalues"
    citation: str = "arXiv:0710.1869, Table I"
    notes: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_SOURCE_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "paper_id", require_nonempty_identifier("paper_id", self.paper_id))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        object.__setattr__(
            self,
            "source_kind",
            require_member("source_kind", self.source_kind, ("table", "equation")),
        )
        object.__setattr__(
            self,
            "locator_label",
            require_nonempty_identifier("locator_label", self.locator_label),
        )
        object.__setattr__(self, "detail", require_nonempty_identifier("detail", self.detail))
        object.__setattr__(self, "citation", require_nonempty_identifier("citation", self.citation))
        if self.notes is not None:
            object.__setattr__(self, "notes", require_nonempty_identifier("notes", self.notes))

    def as_dict(self) -> dict[str, str | None]:
        return {
            "schema_id": self.schema_id,
            "paper_id": self.paper_id,
            "source_kind": self.source_kind,
            "locator_label": self.locator_label,
            "detail": self.detail,
            "citation": self.citation,
            "notes": self.notes,
        }


@dataclass(frozen=True)
class Paper07101869EigenvalueEntry:
    """One generation entry for the sourced Table I profile eigenvalues."""

    generation: int
    c_value: float
    f_value: float

    def __post_init__(self) -> None:
        if not isinstance(self.generation, int) or isinstance(self.generation, bool):
            raise ValueError("generation must be an integer")
        require_member("generation", str(self.generation), ("1", "2", "3"))
        object.__setattr__(self, "c_value", _require_finite_float("c_value", self.c_value))
        object.__setattr__(self, "f_value", require_positive_finite("f_value", self.f_value))

    def as_dict(self) -> dict[str, int | float]:
        return {
            "generation": self.generation,
            "c_value": self.c_value,
            "f_value": self.f_value,
        }


@dataclass(frozen=True)
class Paper07101869SectorTableInputs:
    """One sourced Table I sector block."""

    sector_id: str
    source: Paper07101869SourceRef
    entries: tuple[Paper07101869EigenvalueEntry, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "sector_id",
            require_member("sector_id", self.sector_id, ("Q", "u", "d")),
        )
        if not isinstance(self.source, Paper07101869SourceRef):
            raise ValueError("source must be a Paper07101869SourceRef")
        if len(self.entries) != 3:
            raise ValueError("entries must contain exactly three generations")
        if not all(isinstance(entry, Paper07101869EigenvalueEntry) for entry in self.entries):
            raise ValueError("entries must contain only Paper07101869EigenvalueEntry items")
        generations = tuple(entry.generation for entry in self.entries)
        if generations != (1, 2, 3):
            raise ValueError("entries must be ordered by generation 1, 2, 3")

    def as_dict(self) -> dict[str, object]:
        return {
            "sector_id": self.sector_id,
            "source": self.source.as_dict(),
            "entries": [entry.as_dict() for entry in self.entries],
        }

    @classmethod
    def from_pairs(
        cls,
        *,
        sector_id: str,
        source: Paper07101869SourceRef,
        pairs: tuple[tuple[float, float], tuple[float, float], tuple[float, float]],
    ) -> "Paper07101869SectorTableInputs":
        entries = tuple(
            Paper07101869EigenvalueEntry(
                generation=index,
                c_value=values[0],
                f_value=values[1],
            )
            for index, values in enumerate(pairs, start=1)
        )
        return cls(sector_id=sector_id, source=source, entries=entries)


@dataclass(frozen=True)
class Paper07101869TableIInputs:
    """Closed, sourced copy of the Table I eigenvalue inputs."""

    schema_id: str = PAPER_0710_1869_TABLE_I_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    sectors: tuple[Paper07101869SectorTableInputs, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_TABLE_I_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        if len(self.sectors) != 3:
            raise ValueError("sectors must contain exactly Q, u, and d")
        if not all(isinstance(sector, Paper07101869SectorTableInputs) for sector in self.sectors):
            raise ValueError("sectors must contain only Paper07101869SectorTableInputs")
        sector_ids = tuple(sector.sector_id for sector in self.sectors)
        if sector_ids != ("Q", "u", "d"):
            raise ValueError("sectors must be ordered as Q, u, d")

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "sectors": [sector.as_dict() for sector in self.sectors],
        }


@dataclass(frozen=True)
class Paper07101869Eq3Example:
    """Sourced structural copy of the Eq. (3) example parameters."""

    schema_id: str = PAPER_0710_1869_EQ3_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    source: Paper07101869SourceRef = field(
        default_factory=lambda: Paper07101869SourceRef(
            source_kind="equation",
            locator_label="Eq. (3)",
            detail="example parameters (a, r, theta12, theta23, theta13, delta)",
            citation="arXiv:0710.1869, Eq. (3)",
            notes=(
                "Angles are carried in degrees exactly as extracted in the local paper notes. "
                "delta is carried as the quoted paper parameter without any extra conversion."
            ),
        )
    )
    a: float = 0.8
    r: float = 0.3
    theta12_deg: float = 115.0
    theta23_deg: float = 65.0
    theta13_deg: float = 70.0
    delta: float = 0.6

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_EQ3_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        if not isinstance(self.source, Paper07101869SourceRef):
            raise ValueError("source must be a Paper07101869SourceRef")
        object.__setattr__(self, "a", require_positive_finite("a", self.a))
        object.__setattr__(self, "r", require_positive_finite("r", self.r))
        object.__setattr__(
            self,
            "theta12_deg",
            _require_finite_float("theta12_deg", self.theta12_deg),
        )
        object.__setattr__(
            self,
            "theta23_deg",
            _require_finite_float("theta23_deg", self.theta23_deg),
        )
        object.__setattr__(
            self,
            "theta13_deg",
            _require_finite_float("theta13_deg", self.theta13_deg),
        )
        object.__setattr__(self, "delta", _require_finite_float("delta", self.delta))

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "source": self.source.as_dict(),
            "a": self.a,
            "r": self.r,
            "theta12_deg": self.theta12_deg,
            "theta23_deg": self.theta23_deg,
            "theta13_deg": self.theta13_deg,
            "delta": self.delta,
        }


def default_paper_0710_1869_table_i_inputs() -> Paper07101869TableIInputs:
    """Return the frozen Table I eigenvalue inputs for PR1."""
    q_source = Paper07101869SourceRef(
        source_kind="table",
        locator_label="Table I",
        detail="c_Q/f_Q eigenvalues",
        citation="arXiv:0710.1869, Table I",
    )
    u_source = Paper07101869SourceRef(
        source_kind="table",
        locator_label="Table I",
        detail="c_u/f_u eigenvalues",
        citation="arXiv:0710.1869, Table I",
    )
    d_source = Paper07101869SourceRef(
        source_kind="table",
        locator_label="Table I",
        detail="c_d/f_d eigenvalues",
        citation="arXiv:0710.1869, Table I",
    )
    return Paper07101869TableIInputs(
        sectors=(
            Paper07101869SectorTableInputs.from_pairs(
                sector_id="Q",
                source=q_source,
                pairs=((0.64, 0.002), (0.59, 0.01), (0.46, 0.2)),
            ),
            Paper07101869SectorTableInputs.from_pairs(
                sector_id="u",
                source=u_source,
                pairs=((0.68, 0.0007), (0.53, 0.06), (-0.06, 0.8)),
            ),
            Paper07101869SectorTableInputs.from_pairs(
                sector_id="d",
                source=d_source,
                pairs=((0.65, 0.002), (0.60, 0.008), (0.58, 0.02)),
            ),
        )
    )


def default_paper_0710_1869_eq3_example() -> Paper07101869Eq3Example:
    """Return the frozen Eq. (3) example parameters for PR1."""
    return Paper07101869Eq3Example()
