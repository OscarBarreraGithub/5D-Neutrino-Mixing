"""Base types for the flavor-catalog constraint plugin system.

This module defines:

- :class:`ConstraintLevel` — PRIMARY / SECONDARY / DEFERRED tier.
- :class:`Severity` — HARD / SOFT / INFO failure semantics.
- :class:`ParameterPointExtras` — frozen TypedDict registry of valid
  pre-computed helper keys exposed to constraints. **Adding a new key
  requires editing this TypedDict.** It is the only legitimate
  cross-constraint coordination point in the scaffold.
- :class:`ParameterPoint` — the per-draw input handed to
  :meth:`ConstraintBase.evaluate`.
- :class:`ConstraintResult` — the per-constraint return type.
- :class:`ConstraintBase` — the structural Protocol every constraint
  class must satisfy.

The design rationale is documented in
``.orchestration/PHASE3_SCAFFOLDING_PLAN.md`` §B.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping, Protocol, TypedDict, runtime_checkable


class ConstraintLevel(str, Enum):
    """Tier classification for a constraint."""

    PRIMARY = "PRIMARY"
    SECONDARY = "SECONDARY"
    DEFERRED = "DEFERRED"


class Severity(str, Enum):
    """How a constraint's ``passes=False`` affects a parameter point."""

    HARD = "HARD"   # observed bound, NP must fit inside experimental error
    SOFT = "SOFT"   # SM-vs-exp tension / projection / theory-dominated bound
    INFO = "INFO"   # informational only, never vetoes a point


# --------------------------------------------------------------------- #
# Frozen extras-key registry.
#
# Adding a new key REQUIRES editing this TypedDict. This is the only
# legitimate cross-constraint coordination point in the scaffold.
#
# ``total=False`` because not every point pre-computes every helper —
# constraints that need a key check for its presence via
# ``point.extras.get(...)``.
# --------------------------------------------------------------------- #
class ParameterPointExtras(TypedDict, total=False):
    """Frozen registry of valid ``ParameterPoint.extras`` keys.

    The TypedDict serves three purposes:

    1. Static documentation of which extras are produced by
       :func:`flavor_catalog_constraints.point_builder.build_parameter_point`.
    2. mypy/IDE auto-completion for constraint authors.
    3. A grep-able registry used by ``test_extras_keys_are_declared``
       in the global contract test to catch typos.
    """

    quark_mass_basis_couplings: Any   # QuarkMassBasisCouplings
    lepton_mass_basis_couplings: Any  # LeptonMassBasisCouplings (future)
    kk_gluon_mass_gev: float
    kk_w_mass_gev: float
    kk_z_mass_gev: float
    deltaf2_wilsons: Any              # tuple[DeltaF2WilsonCoefficients, ...]
    # New keys go here. Constraints reference them by literal string.


@dataclass(frozen=True)
class ParameterPoint:
    """A single RS parameter point handed to every constraint.

    Parameters
    ----------
    raw
        Whatever the scan driver produces (warpConfig.Setup + scanParams
        row, opaque to constraints). Constraints should NOT rely on
        ``raw`` having any particular shape; they read from ``extras``.
    extras
        Pre-computed intermediates (see :class:`ParameterPointExtras`).
        Produced exclusively by
        :func:`flavor_catalog_constraints.point_builder.build_parameter_point`.
    """

    raw: Any = None
    extras: ParameterPointExtras = field(default_factory=lambda: ParameterPointExtras())


@dataclass(frozen=True)
class ConstraintResult:
    """Standardised return value of ``ConstraintBase.evaluate``.

    All numeric fields are ``float`` (real). Complex amplitudes belong
    in ``diagnostics`` (see R1 N-5 in the plan).
    """

    process_id: str
    passes: bool
    predicted: float | None              # NP prediction, observable units (real)
    sm_prediction: float | None          # SM central value where known
    experimental: float | None           # central experimental value (or bound)
    ratio: float | None                  # |predicted| / budget (or similar)
    budget: float | None                 # bound used in the ratio
    severity: Severity
    notes: str = ""
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


@runtime_checkable
class ConstraintBase(Protocol):
    """Structural protocol every constraint class satisfies.

    A constraint module declares a *plain class* (no inheritance) with
    the required attributes and an ``evaluate`` method, and decorates
    it with
    :func:`flavor_catalog_constraints.registry.register_constraint`.

    The decorator pins ``level`` and ``family`` as instance attributes
    derived from the module path; constraints must NOT declare them as
    class attributes.
    """

    process_id: str            # "K001"  (REQUIRED class attr)
    severity: Severity         # HARD / SOFT / INFO  (REQUIRED class attr)
    observable: str            # human-readable observable label

    # The following are pinned by the decorator from the module path:
    #   level: ConstraintLevel
    #   family: str

    def evaluate(self, point: ParameterPoint) -> ConstraintResult: ...
