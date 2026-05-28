"""flavor_catalog_constraints — plugin registry for catalogued constraints.

This package implements ~102 flavor-physics constraints documented under
``flavor_catalog/processes/``. Each constraint lives in its own module
(one file per constraint) and is registered via the
:func:`flavor_catalog_constraints.registry.register_constraint`
decorator.

The package intentionally does NOT eagerly walk submodules at import
time (R1 REC-1). Discovery is lazy: callers either invoke
``ConstraintRegistry.discover()`` once at startup (production drivers,
the session-autouse pytest fixture in ``tests/constraints/conftest.py``)
or rely on the lookup helpers (``get``, ``all``, ``filter``,
``evaluate_all``) which trigger discovery on first use.

See ``.orchestration/PHASE3_SCAFFOLDING_PLAN.md`` for the full design.
"""

from __future__ import annotations

from .base import (
    ConstraintBase,
    ConstraintLevel,
    ConstraintResult,
    ParameterPoint,
    ParameterPointExtras,
    Severity,
)
from .registry import ConstraintRegistry, register_constraint

__all__ = [
    "ConstraintBase",
    "ConstraintLevel",
    "ConstraintResult",
    "ParameterPoint",
    "ParameterPointExtras",
    "Severity",
    "ConstraintRegistry",
    "register_constraint",
    "discover",
]


def discover() -> None:
    """Convenience alias for :meth:`ConstraintRegistry.discover`."""
    ConstraintRegistry.discover()
