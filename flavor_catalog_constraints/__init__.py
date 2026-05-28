"""Flavor-catalog constraint plugin system.

Public surface::

    from flavor_catalog_constraints import (
        ConstraintLevel, Severity, ParameterPoint, ConstraintResult,
        ConstraintProtocol, register, evaluate_all, all_constraints,
        make_point, empty_point,
    )

Each constraint is one self-contained file under ``primary/<family>/``
or ``secondary/<family>/`` that decorates a ``Constraint`` class with
:func:`register`. See ``README.md`` and ``TEMPLATE.py``.
"""

from __future__ import annotations

from .anchors import Anchor, AnchorError, load_anchor, load_pdg_block
from .base import (
    ConstraintLevel,
    ConstraintProtocol,
    ConstraintResult,
    ParameterPoint,
    Severity,
)
from .point_builder import KNOWN_EXTRA_KEYS, build_from_quark_couplings, empty_point, make_point
from .registry import (
    all_constraints,
    by_family,
    by_level,
    by_severity,
    discover,
    evaluate_all,
    get,
    import_failures,
    register,
)

__all__ = [
    "ConstraintLevel",
    "Severity",
    "ParameterPoint",
    "ConstraintResult",
    "ConstraintProtocol",
    "Anchor",
    "AnchorError",
    "load_anchor",
    "load_pdg_block",
    "register",
    "discover",
    "evaluate_all",
    "all_constraints",
    "get",
    "by_family",
    "by_level",
    "by_severity",
    "import_failures",
    "make_point",
    "empty_point",
    "build_from_quark_couplings",
    "KNOWN_EXTRA_KEYS",
]
