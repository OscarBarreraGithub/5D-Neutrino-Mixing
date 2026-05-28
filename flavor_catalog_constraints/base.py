"""Base types for the ``flavor_catalog_constraints`` plugin system.

This module is the *contract*. Everything else in the package
(registry, anchors, adapters, the constraint files themselves) is
written against the types declared here, and the contract test pins
every registered constraint to :class:`ConstraintProtocol`.

Layered design (each layer depends only on the ones above it):

1. ``base``                 — types only, no I/O, no physics imports.
2. ``anchors``              — schema-flex YAML loader -> typed views.
3. ``physics_adapters.*``   — the only modules that import physics code.
4. ``registry``             — auto-discovery + exception-isolated bulk eval.
5. ``primary/`` & ``secondary/`` — one self-contained constraint per file.

Numeric discipline
-------------------
Every numeric field on :class:`ConstraintResult` is a real ``float``
(or ``None``). Complex amplitudes, matrices, intermediate Wilson
coefficients, etc. go in ``diagnostics`` (a free-form mapping) — never
in ``predicted`` / ``sm_prediction`` / ``experimental`` / ``ratio`` /
``budget``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Mapping, Protocol, runtime_checkable

__all__ = [
    "ConstraintLevel",
    "Severity",
    "ParameterPoint",
    "ConstraintResult",
    "ConstraintProtocol",
]


class ConstraintLevel(str, Enum):
    """Catalog tier a constraint belongs to.

    Derived from the constraint module's directory (``primary`` vs
    ``secondary``); never declared by the constraint author. ``str``
    mix-in so values serialize cleanly and compare to plain strings.
    """

    PRIMARY = "PRIMARY"
    SECONDARY = "SECONDARY"


class Severity(str, Enum):
    """How a ``passes=False`` result is meant to be interpreted.

    HARD
        Observed experimental bound. New-physics contribution must fit
        inside the experimental error budget. A failure is a veto: the
        parameter point is excluded.
    SOFT
        SM-vs-experiment tension, projection, or theory-dominated
        bound. A failure flags tension but a scan may choose to keep
        the point (it is the caller's policy, not the constraint's).
    INFO
        Informational only. ``passes`` is reported for completeness but
        an INFO constraint must NEVER veto a point.
    """

    HARD = "HARD"
    SOFT = "SOFT"
    INFO = "INFO"


@dataclass(frozen=True)
class ParameterPoint:
    """One RS model draw, as handed to every ``evaluate(point)`` call.

    The point is a *frozen* container so a constraint cannot mutate
    state shared with sibling constraints (isolation). It carries two
    things:

    ``raw``
        The opaque object the scan driver produced for this draw
        (e.g. a warp/geometry setup plus the sampled bulk-mass /
        Yukawa row). Constraints SHOULD NOT poke at ``raw`` directly;
        it is kept only so the point-builder and diagnostics can trace
        provenance.

    ``extras``
        Pre-computed physics intermediates, keyed by string. This is
        the supported way for a constraint to obtain inputs it needs
        (e.g. ``"quark_mass_basis_couplings"``). The set of valid keys
        is documented in :mod:`flavor_catalog_constraints.point_builder`
        (``KNOWN_EXTRA_KEYS``) and pinned by the contract test, so a
        typo'd key is caught rather than silently returning ``None``.

    A constraint reads an extra with ``point.get_extra("key")`` and is
    expected to degrade gracefully (return an INFO/None-ratio result
    with a clear note) when a needed extra is absent, rather than
    raising. Raising is also safe — the registry isolates it — but a
    clean note is friendlier for diagnostics.
    """

    raw: Any = None
    extras: Mapping[str, Any] = field(default_factory=dict)

    def get_extra(self, key: str, default: Any = None) -> Any:
        """Return ``extras[key]`` or ``default`` if absent.

        Convenience accessor so constraints do not assume ``extras`` is
        a plain ``dict`` (it is typed as a ``Mapping``).
        """
        return self.extras.get(key, default)


@dataclass(frozen=True)
class ConstraintResult:
    """Standardized return value of every ``evaluate(point)`` call.

    Parameters
    ----------
    process_id
        Catalog ID, e.g. ``"K001"``. Must match the owning constraint.
    severity
        Copied from the constraint so a caller can decide vetoing
        policy without a second lookup.
    passes
        ``True`` if the point satisfies the constraint. For an INFO
        constraint this is advisory only and must not veto.
    predicted, sm_prediction, experimental
        Real-valued observable quantities (NP prediction, SM central
        value, experimental central value / bound). ``None`` where not
        applicable.
    ratio, budget
        ``ratio`` is the dimensionless saturation of the allowed
        budget (typically ``|predicted| / budget``); ``budget`` is the
        bound used. A point with ``ratio <= 1`` passes a HARD
        constraint. ``None`` where not applicable.
    notes
        Short human-readable explanation (formula used, why it
        passed/failed, why an input was missing).
    diagnostics
        Free-form mapping for anything that is not a clean real scalar:
        complex amplitudes, per-operator Wilson coefficients, the names
        of missing extras, exception metadata, etc.
    """

    process_id: str
    severity: Severity
    passes: bool
    predicted: float | None = None
    sm_prediction: float | None = None
    experimental: float | None = None
    ratio: float | None = None
    budget: float | None = None
    notes: str = ""
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Enforce the real-number contract at construction. bool/int/float
        # and None are fine; a complex amplitude is a contract violation and
        # belongs in ``diagnostics``. Raise (not coerce) so the author fixes
        # the call site rather than silently dropping the imaginary part.
        for name in ("predicted", "sm_prediction", "experimental", "ratio", "budget"):
            value = getattr(self, name)
            if isinstance(value, complex):
                raise TypeError(
                    f"ConstraintResult.{name} must be a real number, got complex "
                    f"{value!r}; put complex amplitudes in `diagnostics` instead."
                )


@runtime_checkable
class ConstraintProtocol(Protocol):
    """Structural type every registered constraint instance satisfies.

    A constraint author writes a *plain class* named ``Constraint``
    with these three class attributes and one method, then decorates it
    with :func:`flavor_catalog_constraints.registry.register`. The
    decorator pins ``level`` and ``family`` from the module path, so
    the author never sets them (and must not).

    Required class attributes
    -------------------------
    ``process_id`` : str
        Catalog ID matching ``^[A-Z]+[0-9]+$`` and the yaml sidecar.
    ``severity`` : Severity
        HARD / SOFT / INFO.
    ``observable`` : str
        Human-readable label for the observable.

    Pinned by the registry decorator (do NOT declare these yourself)
    ----------------------------------------------------------------
    ``level`` : ConstraintLevel
    ``family`` : str
    """

    process_id: str
    severity: Severity
    observable: str
    level: ConstraintLevel
    family: str

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        ...
