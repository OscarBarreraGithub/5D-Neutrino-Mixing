"""Schema-flex loader for catalog YAML sidecars -> typed anchor views.

Why this module exists
----------------------
Each constraint needs the experimental (and often SM) anchor value plus
its provenance, which live in the catalog sidecar
``flavor_catalog/processes/<family>/<ID>.yaml`` under the
``pdg_or_equivalent`` block. That block's *inner* shape varies a lot
across the ~102 entries: some are flat (``value``/``uncertainty`` at the
top), some nest under ``canonical_experimental_value``, some under
``pdg_fit_assuming_cpt`` / ``canonical_limit`` / ``primary_current_limit``
/ etc.

A prior scaffold tried to paper over this by grabbing
``next(iter(block))`` and silently falling back to hardcoded default
strings. That is the failure mode we explicitly reject: a missing or
mismatched anchor must fail **loudly**, not default.

Design
------
- :func:`load_pdg_block` returns the parsed ``pdg_or_equivalent`` mapping
  verbatim and raises :class:`AnchorError` if the file or that key is
  missing.
- :func:`find_block` resolves a value sub-block by trying an explicit,
  ordered list of candidate keys (the constraint author states exactly
  which keys *its* yaml uses). Unknown/empty -> :class:`AnchorError`.
- :func:`build_anchor` extracts a TYPED :class:`Anchor` (real
  ``value`` float, plus optional uncertainty and provenance) from a
  named sub-block, coercing numeric fields and raising loudly on a
  missing/non-numeric value.

A constraint that needs a bespoke multi-field view (e.g. exp + SM
together) composes these helpers; the common single-anchor case is one
call to :func:`load_anchor`.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import yaml

from .base import ConstraintLevel

__all__ = [
    "AnchorError",
    "Anchor",
    "yaml_path_for",
    "load_full_yaml",
    "load_pdg_block",
    "find_block",
    "build_anchor",
    "load_anchor",
    "CATALOG_ROOT",
]

# This file lives at <repo>/flavor_catalog_constraints/anchors.py, so the
# repo root is two parents up.
REPO_ROOT = Path(__file__).resolve().parents[1]
CATALOG_ROOT = REPO_ROOT / "flavor_catalog" / "processes"


class AnchorError(Exception):
    """Raised when an anchor cannot be located or coerced.

    Deliberately a *loud* failure: a constraint whose anchor is missing
    or has the wrong shape is broken and must not silently fall back to
    a default value.
    """


def yaml_path_for(
    process_id: str,
    *,
    family: str,
    tier: ConstraintLevel | str = ConstraintLevel.PRIMARY,
) -> Path:
    """Return the absolute sidecar path for ``process_id``.

    PRIMARY  -> ``flavor_catalog/processes/<family>/<ID>.yaml``
    SECONDARY -> ``flavor_catalog/processes/secondary/<family>/<ID>.yaml``
    """
    tier_str = tier.value.lower() if isinstance(tier, ConstraintLevel) else str(tier).lower()
    if tier_str == "primary":
        return CATALOG_ROOT / family / f"{process_id}.yaml"
    if tier_str == "secondary":
        return CATALOG_ROOT / "secondary" / family / f"{process_id}.yaml"
    raise AnchorError(f"Unknown tier {tier!r} for {process_id}")


@functools.lru_cache(maxsize=None)
def _load_cached(path_str: str) -> Mapping[str, Any]:
    path = Path(path_str)
    if not path.is_file():
        raise AnchorError(f"catalog sidecar not found: {path}")
    with open(path) as f:
        data = yaml.safe_load(f)
    if data is None:
        raise AnchorError(f"catalog sidecar is empty: {path}")
    if not isinstance(data, Mapping):
        raise AnchorError(f"catalog sidecar is not a mapping: {path}")
    return data


def load_full_yaml(
    process_id: str,
    *,
    family: str,
    tier: ConstraintLevel | str = ConstraintLevel.PRIMARY,
) -> Mapping[str, Any]:
    """Return the full parsed yaml document (cached). Loud on missing/empty."""
    return _load_cached(str(yaml_path_for(process_id, family=family, tier=tier)))


def load_pdg_block(
    process_id: str,
    *,
    family: str,
    tier: ConstraintLevel | str = ConstraintLevel.PRIMARY,
) -> Mapping[str, Any]:
    """Return the ``pdg_or_equivalent`` mapping verbatim. Loud if absent."""
    data = load_full_yaml(process_id, family=family, tier=tier)
    block = data.get("pdg_or_equivalent")
    if block is None:
        raise AnchorError(
            f"{process_id}: missing top-level 'pdg_or_equivalent' key in sidecar"
        )
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: 'pdg_or_equivalent' is not a mapping (got {type(block).__name__})"
        )
    return block


def find_block(
    block: Mapping[str, Any],
    candidates: Sequence[str],
    *,
    process_id: str = "?",
) -> Mapping[str, Any]:
    """Return the first sub-mapping in ``block`` whose key is in ``candidates``.

    ``candidates`` is an explicit, ordered list supplied by the
    constraint author — the keys *its* yaml is known to use. This is the
    schema-flex hook: we never guess via ``next(iter(...))``; the author
    declares the allowed shapes, and an unmatched block raises loudly so
    a yaml refactor that drops the expected key is caught immediately.

    Raises
    ------
    AnchorError
        If none of ``candidates`` is present, or the matched value is
        not itself a mapping.
    """
    for key in candidates:
        if key in block:
            sub = block[key]
            if not isinstance(sub, Mapping):
                raise AnchorError(
                    f"{process_id}: anchor sub-block '{key}' is not a mapping "
                    f"(got {type(sub).__name__})"
                )
            return sub
    raise AnchorError(
        f"{process_id}: none of the expected anchor keys {list(candidates)} "
        f"found in pdg_or_equivalent (present keys: {sorted(block)})"
    )


@dataclass(frozen=True)
class Anchor:
    """Typed view over one experimental/SM value block.

    All numeric fields are real ``float`` (or ``None``). Provenance
    fields are strings carried verbatim from the yaml for traceability.
    """

    process_id: str
    block_key: str                    # which sub-block this came from
    value: float                      # the central value / bound (REQUIRED)
    uncertainty: float | None = None  # symmetric uncertainty if present
    value_id: str | None = None
    observable: str | None = None
    units: str | None = None
    confidence_level: str | None = None
    source: str | None = None
    source_url: str | None = None
    year: int | None = None
    snapshot_path: str | None = None

    def with_block(self, sub: Mapping[str, Any]) -> "Anchor":  # pragma: no cover
        """Reserved for future extension (kept for API symmetry)."""
        raise NotImplementedError


def _coerce_float(x: Any, *, process_id: str, field_name: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: anchor field '{field_name}'={x!r} is not numeric"
        ) from exc


def _opt_float(x: Any, *, process_id: str, field_name: str) -> float | None:
    if x is None:
        return None
    return _coerce_float(x, process_id=process_id, field_name=field_name)


def build_anchor(
    sub: Mapping[str, Any],
    *,
    process_id: str,
    block_key: str,
    value_key: str = "value",
    uncertainty_key: str = "uncertainty",
) -> Anchor:
    """Build a typed :class:`Anchor` from a resolved value sub-block.

    Coerces ``value`` to a real float (REQUIRED — loud on missing or
    non-numeric) and pulls common provenance fields when present.
    ``uncertainty`` is optional (some bound-style entries have none).
    """
    if value_key not in sub:
        raise AnchorError(
            f"{process_id}: anchor block '{block_key}' has no '{value_key}' field "
            f"(keys: {sorted(sub)})"
        )
    return Anchor(
        process_id=process_id,
        block_key=block_key,
        value=_coerce_float(sub[value_key], process_id=process_id, field_name=value_key),
        uncertainty=_opt_float(
            sub.get(uncertainty_key), process_id=process_id, field_name=uncertainty_key
        ),
        value_id=_as_str(sub.get("value_id")),
        observable=_as_str(sub.get("observable")),
        units=_as_str(sub.get("units")),
        confidence_level=_as_str(_first_present(sub, "confidence_level", "cl")),
        source=_as_str(sub.get("source")),
        source_url=_as_str(sub.get("source_url")),
        year=_as_int(sub.get("year")),
        snapshot_path=_as_str(sub.get("snapshot_path")),
    )


def _as_str(x: Any) -> str | None:
    return None if x is None else str(x)


def _first_present(block: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in block:
            return block[key]
    return None


def _as_int(x: Any) -> int | None:
    if x is None:
        return None
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def load_anchor(
    process_id: str,
    *,
    family: str,
    candidates: Sequence[str],
    tier: ConstraintLevel | str = ConstraintLevel.PRIMARY,
    value_key: str = "value",
    uncertainty_key: str = "uncertainty",
    expected_value_id: str | None = None,
    expected_block_key: str | None = None,
    expected_units: str | None = None,
    expected_confidence_level: str | None = None,
) -> Anchor:
    """One-call convenience for the common single-anchor case.

    Loads the sidecar's ``pdg_or_equivalent`` block, resolves the first
    matching sub-block from ``candidates``, and returns a typed
    :class:`Anchor`. Any failure (missing file, missing block, missing
    value) raises :class:`AnchorError`.

    Optional ``expected_*`` arguments are additive validation guards. They
    preserve existing call behavior when omitted, and raise
    :class:`AnchorError` when a loaded sidecar silently drifts away from the
    caller's intended value id, block key, units, or confidence level.
    """
    block = load_pdg_block(process_id, family=family, tier=tier)
    matched_key = next((k for k in candidates if k in block), None)
    sub = find_block(block, candidates, process_id=process_id)
    anchor = build_anchor(
        sub,
        process_id=process_id,
        block_key=matched_key or "?",
        value_key=value_key,
        uncertainty_key=uncertainty_key,
    )
    _validate_expected(
        anchor,
        expected_value_id=expected_value_id,
        expected_block_key=expected_block_key,
        expected_units=expected_units,
        expected_confidence_level=expected_confidence_level,
    )
    return anchor


def _validate_expected(
    anchor: Anchor,
    *,
    expected_value_id: str | None,
    expected_block_key: str | None,
    expected_units: str | None,
    expected_confidence_level: str | None,
) -> None:
    checks = (
        ("value_id", expected_value_id, anchor.value_id),
        ("block_key", expected_block_key, anchor.block_key),
        ("units", expected_units, anchor.units),
        ("confidence_level", expected_confidence_level, anchor.confidence_level),
    )
    for field_name, expected, actual in checks:
        if expected is not None and actual != expected:
            raise AnchorError(
                f"{anchor.process_id}: anchor {field_name} mismatch: "
                f"expected {expected!r}, got {actual!r}"
            )
