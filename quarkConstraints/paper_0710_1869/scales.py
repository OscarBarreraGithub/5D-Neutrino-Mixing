"""Explicit scale objects for the dedicated paper-facing 0710.1869 mode."""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

from .conventions import PAPER_0710_1869_MODE_ID
from .validation import (
    normalize_optional_positive_finite,
    require_known_schema_id,
    require_member,
    require_nonempty_identifier,
    require_positive_finite,
)

PAPER_0710_1869_SCALES_SCHEMA_ID = "quarkConstraints.paper_0710_1869.scales.v1"


@dataclass(frozen=True)
class Paper07101869ScalePoint:
    """Explicit, non-overloaded scale object for one paper-mode point."""

    schema_id: str = PAPER_0710_1869_SCALES_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    label: str = "central"
    Lambda_IR_GeV: float = 3000.0
    m_g1_GeV: float = 3000.0
    xi_g: float | None = None
    mu_match_GeV: float = 3000.0
    mu_gs_GeV: float = 3000.0
    m_KK_eff_GeV: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_SCALES_SCHEMA_ID,
            ),
        )
        object.__setattr__(self, "label", require_nonempty_identifier("label", self.label))
        object.__setattr__(
            self,
            "Lambda_IR_GeV",
            require_positive_finite("Lambda_IR_GeV", self.Lambda_IR_GeV),
        )
        object.__setattr__(self, "m_g1_GeV", require_positive_finite("m_g1_GeV", self.m_g1_GeV))
        derived_xi_g = float(self.m_g1_GeV / self.Lambda_IR_GeV)
        resolved_xi_g = (
            derived_xi_g if self.xi_g is None else require_positive_finite("xi_g", self.xi_g)
        )
        if not math.isclose(resolved_xi_g, derived_xi_g, rel_tol=1e-12, abs_tol=0.0):
            raise ValueError("xi_g must equal m_g1_GeV / Lambda_IR_GeV")
        object.__setattr__(self, "xi_g", resolved_xi_g)
        object.__setattr__(
            self,
            "mu_match_GeV",
            require_positive_finite("mu_match_GeV", self.mu_match_GeV),
        )
        object.__setattr__(self, "mu_gs_GeV", require_positive_finite("mu_gs_GeV", self.mu_gs_GeV))
        object.__setattr__(
            self,
            "m_KK_eff_GeV",
            normalize_optional_positive_finite("m_KK_eff_GeV", self.m_KK_eff_GeV),
        )
        object.__setattr__(self, "mode_id", require_nonempty_identifier("mode_id", self.mode_id))
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))

    @property
    def propagator_mass_GeV(self) -> float:
        """Mass scale used in the heavy propagator."""
        return self.m_g1_GeV if self.m_KK_eff_GeV is None else self.m_KK_eff_GeV

    @property
    def has_explicit_effective_kk_scale(self) -> bool:
        """Whether a separate effective KK scale is carried explicitly."""
        return self.m_KK_eff_GeV is not None

    def as_dict(self) -> dict[str, float | str | None]:
        """Return a stable mapping representation."""
        return asdict(self)


def default_paper_0710_1869_scales() -> Paper07101869ScalePoint:
    """Return the default explicit scale point."""
    return Paper07101869ScalePoint()
