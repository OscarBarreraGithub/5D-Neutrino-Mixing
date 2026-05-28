"""Single source of truth for :class:`ParameterPoint` construction.

The constraint zoo consumes ``ParameterPoint.extras`` keys declared in
:class:`flavor_catalog_constraints.base.ParameterPointExtras`. This
module is the **only** place those extras are populated.

Status (Phase 3 scaffold)
-------------------------
The full builder requires a stable scan-row shape and a
``build_quark_mass_basis_couplings(scan_row)`` helper on the scan-driver
side. Both are TBD per ``.orchestration/PHASE3_SCAFFOLDING_PLAN.md``
§K (last open item). Until those are pinned, the builder is a
skeleton:

- :func:`build_parameter_point` raises :class:`NotImplementedError`
  with a pointer to the plan.
- :func:`make_empty_point` and :func:`make_point` are usable for
  tests and unit experimentation. They do **not** populate the
  ``quark_mass_basis_couplings`` extra — tests that need it must
  inject one explicitly (see ``tests/constraints/conftest.py``).
"""

from __future__ import annotations

from typing import Any

from .base import ParameterPoint, ParameterPointExtras

__all__ = [
    "build_parameter_point",
    "make_empty_point",
    "make_point",
]


def make_empty_point() -> ParameterPoint:
    """Return a :class:`ParameterPoint` with no extras.

    Useful for constraint smoke tests that don't depend on any
    pre-computed helper. Most realistic constraints will need
    :func:`make_point` (or the future
    :func:`build_parameter_point`) instead.
    """
    return ParameterPoint(raw=None, extras=ParameterPointExtras())


def make_point(
    raw: Any = None,
    *,
    extras: ParameterPointExtras | None = None,
) -> ParameterPoint:
    """Hand-build a :class:`ParameterPoint` from caller-supplied pieces.

    The ``extras`` argument must be a mapping whose keys are declared
    in :class:`ParameterPointExtras`. The contract test
    ``test_extras_keys_are_declared`` statically catches stray keys
    referenced inside constraint modules; ``make_point`` itself does
    not enforce key membership at runtime so that test fixtures can
    inject mocked helpers freely.
    """
    if extras is None:
        extras = ParameterPointExtras()
    return ParameterPoint(raw=raw, extras=extras)


def build_parameter_point(scan_row: Any) -> ParameterPoint:
    """Construct a :class:`ParameterPoint` from a scan-driver row.

    **TBD.** The exact shape of ``scan_row`` is owned by
    ``scripts/run_rs_anarchy.py`` and friends. Once that integration
    happens, this function must:

    1. Call ``build_quark_mass_basis_couplings(scan_row)`` from
       ``quarkConstraints.couplings`` (helper does not yet exist —
       see plan §K).
    2. Cache the Delta F=2 Wilson coefficients via
       ``compute_delta_f2_wilsons(couplings, system="epsilon_k")``.
    3. Materialise KK gauge-boson masses from ``scan_row``.
    4. Populate :class:`ParameterPointExtras` and return.

    Until then this raises :class:`NotImplementedError` with a pointer
    so a caller does not silently get a half-built point.

    See ``.orchestration/PHASE3_SCAFFOLDING_PLAN.md`` §K for the open
    items and the scan-driver integration plan.
    """
    raise NotImplementedError(
        "build_parameter_point() is a scaffold stub. See "
        ".orchestration/PHASE3_SCAFFOLDING_PLAN.md §K — the "
        "scan-row shape and the build_quark_mass_basis_couplings() "
        "helper are TBD."
    )
