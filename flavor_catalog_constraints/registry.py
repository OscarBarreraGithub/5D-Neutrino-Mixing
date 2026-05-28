"""Plugin registry for ``flavor_catalog_constraints``.

Responsibilities:

- The ``@register_constraint`` decorator (one line in each constraint
  module) instantiates the class, validates ``process_id``, derives
  ``level`` and ``family`` from the module path, checks the yaml
  sidecar exists, and installs the instance into the global registry.
- :class:`ConstraintRegistry` exposes lazy discovery, lookup, filter,
  and bulk-evaluation helpers. Discovery is idempotent and runs at
  most once per Python process. Import errors in any individual
  constraint module are captured into ``_IMPORT_FAILURES`` rather than
  propagating, so one broken module does not poison the whole
  registry (R1 M-5).

See ``.orchestration/PHASE3_SCAFFOLDING_PLAN.md`` §C for the full
design.
"""

from __future__ import annotations

import importlib
import pkgutil
import re
from typing import Any, Dict

from .base import (
    ConstraintBase,
    ConstraintLevel,
    ConstraintResult,
    ParameterPoint,
    Severity,
)

_REGISTRY: Dict[str, ConstraintBase] = {}
_IMPORT_FAILURES: Dict[str, BaseException] = {}
_DISCOVERED: bool = False

_ID_RE = re.compile(r"^[A-Z]+[0-9]+$")
_VALID_FAMILIES = {
    "beauty",
    "charged_lepton",
    "charm",
    "collider_rs",
    "edm_neutrino",
    "kaon",
    "top_higgs_ew",
}


def _derive_tier_and_family_from_module(module_name: str) -> tuple[ConstraintLevel, str]:
    """Derive ``(level, family)`` from a constraint module's dotted name.

    Examples
    --------
    ``flavor_catalog_constraints.primary.kaon.K001``     -> (PRIMARY, "kaon")
    ``flavor_catalog_constraints.secondary.beauty.B007`` -> (SECONDARY, "beauty")
    """
    parts = module_name.split(".")
    # Expect: ["flavor_catalog_constraints", "<tier>", "<family>", "<ID>"]
    if len(parts) < 4 or parts[0] != "flavor_catalog_constraints":
        raise RuntimeError(
            f"Cannot derive tier/family from module name {module_name!r}: "
            "constraint modules must live at "
            "flavor_catalog_constraints.<tier>.<family>.<ID>"
        )
    tier_str, family = parts[1], parts[2]
    tier_map = {
        "primary": ConstraintLevel.PRIMARY,
        "secondary": ConstraintLevel.SECONDARY,
    }
    if tier_str not in tier_map:
        raise RuntimeError(
            f"Unknown tier '{tier_str}' in module {module_name!r}; "
            f"expected one of {sorted(tier_map)}"
        )
    if family not in _VALID_FAMILIES:
        raise RuntimeError(
            f"Unknown family '{family}' in module {module_name!r}; "
            f"expected one of {sorted(_VALID_FAMILIES)}"
        )
    return tier_map[tier_str], family


def register_constraint(cls: type) -> type:
    """Decorator: validate, instantiate, pin metadata, and register.

    The decorator runs at module import time. It performs:

    1. ``process_id`` format check (``^[A-Z]+[0-9]+$``).
    2. Tier + family derivation from the module path, pinned as
       instance attributes via :func:`object.__setattr__` (which works
       on frozen dataclasses and plain classes alike).
    3. Yaml sidecar existence check.
    4. Duplicate-ID guard.
    5. Installation into the global ``_REGISTRY``.

    Returns the original class unchanged so the decorator is
    composable with type-checkers and IDEs.
    """
    instance = cls()

    # 1. Validate process_id
    pid = getattr(instance, "process_id", None)
    if not isinstance(pid, str) or not _ID_RE.match(pid):
        raise RuntimeError(
            f"{cls.__module__}: process_id {pid!r} does not match {_ID_RE.pattern}"
        )

    # 2. Derive tier + family from module path; pin as instance attributes
    level, family = _derive_tier_and_family_from_module(cls.__module__)
    object.__setattr__(instance, "level", level)
    object.__setattr__(instance, "family", family)

    # 3. Check yaml sidecar exists (REC-3)
    from .anchors import yaml_path_for

    yaml_path = yaml_path_for(pid, family=family, tier=level)
    if not yaml_path.is_file():
        raise RuntimeError(
            f"{cls.__module__}: yaml sidecar missing at {yaml_path}"
        )

    # 4. Duplicate guard
    if pid in _REGISTRY:
        prev = _REGISTRY[pid].__class__.__module__
        raise RuntimeError(
            f"Duplicate constraint registration for {pid}: "
            f"already provided by {prev}, now also by {cls.__module__}"
        )

    # 5. Install
    _REGISTRY[pid] = instance
    return cls


class ConstraintRegistry:
    """Static-method facade over the module-global ``_REGISTRY`` dict.

    All public methods trigger lazy ``discover()`` on first call.
    """

    @staticmethod
    def discover() -> None:
        """Walk ``primary`` and ``secondary`` subpackages and import every leaf module.

        Idempotent: subsequent calls are no-ops. Per-module exceptions
        are captured into ``_IMPORT_FAILURES`` so one broken
        constraint module does not poison the rest of the registry.
        """
        global _DISCOVERED
        if _DISCOVERED:
            return
        import flavor_catalog_constraints.primary as primary
        import flavor_catalog_constraints.secondary as secondary

        for pkg in (primary, secondary):
            for _, modname, ispkg in pkgutil.walk_packages(
                pkg.__path__, pkg.__name__ + "."
            ):
                if ispkg:
                    continue  # only leaf modules carry constraints
                try:
                    importlib.import_module(modname)
                except Exception as exc:  # noqa: BLE001
                    _IMPORT_FAILURES[modname] = exc
        _DISCOVERED = True

    @staticmethod
    def import_failures() -> Dict[str, BaseException]:
        """Return a shallow copy of the per-module import-failure map."""
        return dict(_IMPORT_FAILURES)

    @staticmethod
    def reset_for_tests() -> None:
        """Test-only escape hatch to reset the discovery flag.

        Used by ``tests/constraints/conftest.py`` if discovery state
        needs to be re-bootstrapped between test sessions. Not part of
        the production API.
        """
        global _DISCOVERED
        _DISCOVERED = False
        _IMPORT_FAILURES.clear()

    @staticmethod
    def all() -> Dict[str, ConstraintBase]:
        """Return a shallow copy of the full ``{process_id: instance}`` map."""
        ConstraintRegistry.discover()
        return dict(_REGISTRY)

    @staticmethod
    def get(process_id: str) -> ConstraintBase:
        """Return the registered instance for ``process_id``."""
        ConstraintRegistry.discover()
        return _REGISTRY[process_id]

    @staticmethod
    def filter(
        *,
        family: str | None = None,
        level: ConstraintLevel | None = None,
        severity: Severity | None = None,
    ) -> Dict[str, ConstraintBase]:
        """Return the registry filtered by family / level / severity (intersection)."""
        ConstraintRegistry.discover()
        out: Dict[str, ConstraintBase] = {}
        for pid, c in _REGISTRY.items():
            if family is not None and getattr(c, "family", None) != family:
                continue
            if level is not None and getattr(c, "level", None) != level:
                continue
            if severity is not None and getattr(c, "severity", None) != severity:
                continue
            out[pid] = c
        return out

    @staticmethod
    def evaluate_all(
        point: ParameterPoint,
        *,
        include_deferred: bool = False,
    ) -> Dict[str, ConstraintResult]:
        """Evaluate every registered constraint at ``point``.

        Each ``constraint.evaluate(point)`` call is wrapped in a
        try/except (R1 M-4). A constraint that raises is reported as
        ``passes=False`` with a diagnostic note; the other 101
        constraints still run.

        DEFERRED constraints are skipped by default (R1 N-4);
        ``include_deferred=True`` re-enables them for diagnostics.
        """
        ConstraintRegistry.discover()
        out: Dict[str, ConstraintResult] = {}
        for pid, c in _REGISTRY.items():
            level = getattr(c, "level", None)
            if level == ConstraintLevel.DEFERRED and not include_deferred:
                continue
            try:
                out[pid] = c.evaluate(point)
            except Exception as exc:  # noqa: BLE001
                out[pid] = ConstraintResult(
                    process_id=pid,
                    passes=False,
                    predicted=None,
                    sm_prediction=None,
                    experimental=None,
                    ratio=None,
                    budget=None,
                    severity=getattr(c, "severity", Severity.HARD),
                    notes=f"evaluate() raised {type(exc).__name__}: {exc}",
                    diagnostics={"exception_type": type(exc).__name__},
                )
        return out


__all__ = [
    "ConstraintRegistry",
    "register_constraint",
]
