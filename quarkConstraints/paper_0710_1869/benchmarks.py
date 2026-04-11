"""Paper-owned benchmark definitions for the dedicated ``paper_0710_1869`` mode."""

from __future__ import annotations

from dataclasses import dataclass, field

from .conventions import PAPER_0710_1869_MODE_ID, PAPER_0710_1869_PAPER_ID
from .inputs import (
    Paper07101869Eq3Example,
    Paper07101869TableIInputs,
    default_paper_0710_1869_eq3_example,
    default_paper_0710_1869_table_i_inputs,
)
from .validation import require_known_schema_id, require_member, require_nonempty_identifier

PAPER_0710_1869_BENCHMARK_SCHEMA_ID = "quarkConstraints.paper_0710_1869.benchmarks.v1"


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


def default_paper_0710_1869_pr1_benchmark() -> Paper07101869Benchmark:
    """Return the PR1 structural benchmark assembled from frozen paper inputs."""
    return Paper07101869Benchmark()


def paper_0710_1869_pr1_benchmarks() -> tuple[Paper07101869Benchmark, ...]:
    """Return the closed set of sourced PR1 benchmark definitions."""
    return (default_paper_0710_1869_pr1_benchmark(),)
