"""Schema-flex YAML sidecar loader for catalogued constraints.

Round-1 of the scaffolding plan hard-coded
``data["pdg_or_equivalent"]["canonical_experimental_value"]``, a path
that exists in **exactly 1 of 102 yamls**. The aggregate count of
top-level keys inside ``pdg_or_equivalent`` is ~30 distinct shapes
(``canonical_average``, ``canonical_limit``, ``primary_current_limit``,
``values:`` list, family-prefixed keys, flat ``source:``/``year:``,
…).

The R2 design (this module) therefore returns the parsed
``pdg_or_equivalent`` mapping **verbatim**. Each constraint defines
its own typed view inline via a private ``_Anchor`` dataclass and a
``_build_anchor(raw)`` function that knows the keys appropriate to
*its own* yaml. ``test_anchor_matches_yaml`` in each
``tests/constraints/.../test_<ID>.py`` pins the per-constraint
extraction back against the yaml.

See ``.orchestration/PHASE3_SCAFFOLDING_PLAN.md`` §E for the worked
patterns A through D.
"""

from __future__ import annotations

import functools
from pathlib import Path
from typing import Any, Mapping

import yaml

from .base import ConstraintLevel

# Repository root: this file lives at
#   <root>/flavor_catalog_constraints/anchors.py
# so ``parents[1]`` is the repo root.
REPO_ROOT = Path(__file__).resolve().parents[1]
CATALOG_ROOT = REPO_ROOT / "flavor_catalog" / "processes"


def yaml_path_for(
    process_id: str,
    *,
    family: str,
    tier: ConstraintLevel | str,
) -> Path:
    """Return the absolute path to a constraint's yaml sidecar.

    PRIMARY constraints live at
    ``flavor_catalog/processes/<family>/<process_id>.yaml``.
    SECONDARY constraints live at
    ``flavor_catalog/processes/secondary/<family>/<process_id>.yaml``.
    """
    tier_str = tier.value.lower() if isinstance(tier, ConstraintLevel) else str(tier).lower()
    if tier_str == "primary":
        return CATALOG_ROOT / family / f"{process_id}.yaml"
    if tier_str == "secondary":
        return CATALOG_ROOT / "secondary" / family / f"{process_id}.yaml"
    raise ValueError(f"Unknown tier {tier!r}")


@functools.lru_cache(maxsize=None)
def _load_full_yaml(path_str: str) -> Mapping[str, Any]:
    """Cache the parsed yaml document keyed by absolute path string."""
    with open(path_str) as f:
        return yaml.safe_load(f)


def load_raw(
    process_id: str,
    *,
    family: str,
    tier: ConstraintLevel | str = ConstraintLevel.PRIMARY,
) -> Mapping[str, Any]:
    """Return the parsed ``pdg_or_equivalent`` dict verbatim.

    The schema of this dict varies across the catalog. Each constraint
    is responsible for picking the field(s) it needs. This loader makes
    no assumptions beyond the presence of the top-level
    ``pdg_or_equivalent`` key (verified to be present in all 102+
    catalog yamls).

    Raises
    ------
    KeyError
        If the yaml does not contain a ``pdg_or_equivalent`` top-level
        key (loud failure — constraint cannot be evaluated).
    """
    path = yaml_path_for(process_id, family=family, tier=tier)
    data = _load_full_yaml(str(path))
    if data is None:
        raise KeyError(f"{path}: yaml document is empty")
    raw = data.get("pdg_or_equivalent")
    if raw is None:
        raise KeyError(
            f"{path}: missing top-level 'pdg_or_equivalent' key"
        )
    return raw


def load_full(
    process_id: str,
    *,
    family: str,
    tier: ConstraintLevel | str = ConstraintLevel.PRIMARY,
) -> Mapping[str, Any]:
    """Return the full parsed yaml document.

    Mostly used by the global contract test to cross-check fields the
    constraint extracts against the canonical yaml.
    """
    path = yaml_path_for(process_id, family=family, tier=tier)
    return _load_full_yaml(str(path))
