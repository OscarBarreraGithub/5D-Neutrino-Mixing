"""Auto-discovery plugin registry with per-constraint exception isolation.

A constraint registers itself by decorating its ``Constraint`` class
with :func:`register`. Discovery walks the ``primary`` and ``secondary``
subpackages and imports every leaf module; the decorator fires at import
time and installs an instance into the global registry.

Isolation guarantees (the whole point of the design):

- **One broken module cannot poison discovery.** Each leaf module is
  imported in its own ``try/except``; an import-time failure (syntax
  error, bad anchor, decorator raise) is recorded in
  :func:`import_failures` instead of propagating.
- **One raising ``evaluate`` cannot abort the sweep.** Every call in
  :func:`evaluate_all` is wrapped; a constraint that raises yields a
  failing :class:`ConstraintResult` with the exception in
  ``diagnostics`` and the other constraints still run.

The registry holds no hand-maintained list of IDs: the set of
constraints *is* the set of files under ``primary/`` and ``secondary/``.
"""

from __future__ import annotations

import importlib
import pkgutil
import re
from typing import Dict, Optional

from .anchors import yaml_path_for
from .base import (
    ConstraintLevel,
    ConstraintProtocol,
    ConstraintResult,
    ParameterPoint,
    Severity,
)

__all__ = [
    "register",
    "discover",
    "all_constraints",
    "get",
    "by_family",
    "by_level",
    "by_severity",
    "evaluate_all",
    "import_failures",
    "reset_for_tests",
    "VALID_FAMILIES",
]

VALID_FAMILIES = frozenset(
    {
        "beauty",
        "charged_lepton",
        "charm",
        "collider_rs",
        "edm_neutrino",
        "kaon",
        "top_higgs_ew",
    }
)

_ID_RE = re.compile(r"^[A-Z]+[0-9]+$")

# Module-global state. Discovery is idempotent; populated once per process.
_REGISTRY: Dict[str, ConstraintProtocol] = {}
_IMPORT_FAILURES: Dict[str, BaseException] = {}
_DISCOVERED = False


def _derive_level_and_family(module_name: str) -> tuple[ConstraintLevel, str]:
    """Derive ``(level, family)`` from a constraint module's dotted path.

    Expected layout::

        flavor_catalog_constraints.<tier>.<family>.<ID>

    e.g. ``flavor_catalog_constraints.primary.kaon.K001`` ->
    ``(PRIMARY, "kaon")``.
    """
    parts = module_name.split(".")
    if len(parts) != 4 or parts[0] != "flavor_catalog_constraints":
        raise RuntimeError(
            f"constraint module {module_name!r} must live at "
            "flavor_catalog_constraints.<tier>.<family>.<ID>"
        )
    tier_str, family = parts[1], parts[2]
    tier_map = {"primary": ConstraintLevel.PRIMARY, "secondary": ConstraintLevel.SECONDARY}
    if tier_str not in tier_map:
        raise RuntimeError(f"unknown tier {tier_str!r} in {module_name!r}")
    if family not in VALID_FAMILIES:
        raise RuntimeError(
            f"unknown family {family!r} in {module_name!r}; "
            f"expected one of {sorted(VALID_FAMILIES)}"
        )
    return tier_map[tier_str], family


def register(cls: type) -> type:
    """Class decorator: validate, instantiate, pin metadata, register.

    Runs at import time and performs, in order:

    1. ``process_id`` format check (``^[A-Z]+[0-9]+$``).
    2. ``severity`` is a :class:`Severity`.
    3. Tier + family derivation from the module path, pinned onto the
       instance as ``level`` / ``family`` (so the author never sets
       them — keeping the file portable across tiers).
    4. Yaml sidecar existence check (a constraint without its catalog
       entry is rejected loudly).
    5. Duplicate-ID guard.

    Returns the class unchanged so type-checkers and IDEs still see the
    original definition.
    """
    instance = cls()

    pid = getattr(instance, "process_id", None)
    if not isinstance(pid, str) or not _ID_RE.match(pid):
        raise RuntimeError(f"{cls.__module__}: process_id {pid!r} must match {_ID_RE.pattern}")

    if not isinstance(getattr(instance, "severity", None), Severity):
        raise RuntimeError(f"{cls.__module__}: severity must be a Severity enum")

    level, family = _derive_level_and_family(cls.__module__)
    # object.__setattr__ works on both plain and frozen-dataclass instances.
    object.__setattr__(instance, "level", level)
    object.__setattr__(instance, "family", family)

    sidecar = yaml_path_for(pid, family=family, tier=level)
    if not sidecar.is_file():
        raise RuntimeError(f"{cls.__module__}: catalog sidecar missing at {sidecar}")

    # Duplicate process_id guard: the FIRST module to import wins (it is
    # already in _REGISTRY); the second raises here and discover() isolates
    # it as an import failure. Which one survives is therefore import-order
    # dependent (pkgutil walk order), not deterministic by content.
    if pid in _REGISTRY:
        prev = type(_REGISTRY[pid]).__module__
        raise RuntimeError(
            f"duplicate constraint id {pid!r}: registered by {prev}, "
            f"now also by {cls.__module__}"
        )

    _REGISTRY[pid] = instance  # type: ignore[assignment]
    return cls


def discover() -> None:
    """Import every leaf constraint module under primary/ and secondary/.

    Idempotent. Per-module import failures are captured in
    :func:`import_failures`, never raised, so a single broken file does
    not break discovery of the rest.
    """
    global _DISCOVERED
    if _DISCOVERED:
        return
    # Import the tier subpackages. These must exist (empty __init__.py).
    from . import primary, secondary

    for pkg in (primary, secondary):
        for _, modname, ispkg in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
            if ispkg:
                continue  # family subpackages carry no constraints themselves
            try:
                importlib.import_module(modname)
            except Exception as exc:  # noqa: BLE001 — isolation is the point
                _IMPORT_FAILURES[modname] = exc
    _DISCOVERED = True


def import_failures() -> Dict[str, BaseException]:
    """Return a copy of the ``{module_name: exception}`` import-failure map."""
    discover()
    return dict(_IMPORT_FAILURES)


def reset_for_tests() -> None:
    """Clear discovery state. Test-only; not part of the production API."""
    global _DISCOVERED
    _REGISTRY.clear()
    _IMPORT_FAILURES.clear()
    _DISCOVERED = False


def all_constraints() -> Dict[str, ConstraintProtocol]:
    """Return a shallow copy of the ``{process_id: instance}`` registry."""
    discover()
    return dict(_REGISTRY)


def get(process_id: str) -> ConstraintProtocol:
    """Return the registered instance for ``process_id`` (KeyError if absent)."""
    discover()
    return _REGISTRY[process_id]


def by_family(family: str) -> Dict[str, ConstraintProtocol]:
    discover()
    return {pid: c for pid, c in _REGISTRY.items() if getattr(c, "family", None) == family}


def by_level(level: ConstraintLevel) -> Dict[str, ConstraintProtocol]:
    discover()
    return {pid: c for pid, c in _REGISTRY.items() if getattr(c, "level", None) == level}


def by_severity(severity: Severity) -> Dict[str, ConstraintProtocol]:
    discover()
    return {pid: c for pid, c in _REGISTRY.items() if getattr(c, "severity", None) == severity}


def evaluate_all(
    point: ParameterPoint,
    *,
    only_family: Optional[str] = None,
) -> Dict[str, ConstraintResult]:
    """Evaluate every registered constraint at ``point``, isolating failures.

    Each ``evaluate`` is wrapped: a constraint that raises is reported
    as a failing :class:`ConstraintResult` carrying the exception in
    ``diagnostics``, and the remaining constraints still run.

    ``only_family`` optionally restricts evaluation to one family.
    """
    discover()
    out: Dict[str, ConstraintResult] = {}
    for pid, c in _REGISTRY.items():
        if only_family is not None and getattr(c, "family", None) != only_family:
            continue
        try:
            out[pid] = c.evaluate(point)
        except Exception as exc:  # noqa: BLE001 — isolation is the point
            out[pid] = ConstraintResult(
                process_id=pid,
                severity=getattr(c, "severity", Severity.HARD),
                passes=False,
                notes=f"evaluate() raised {type(exc).__name__}: {exc}",
                diagnostics={"exception_type": type(exc).__name__, "exception": str(exc)},
            )
    return out
