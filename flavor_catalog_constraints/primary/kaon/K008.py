"""K008 - Long-lived neutral-kaon rare electron mode ``K_L -> pi0 e+ e-``.

Severity rationale
------------------
INFO placeholder. This observable is a Delta S = 1 rare kaon decay, not
a Delta F = 2 kaon-mixing constraint like K001/K002, not the
nonleptonic eps'/eps observable in K003, and not a K -> pi nu nubar
mode like K004/K005. The measured branching-fraction limit constrains
a mode with short-distance + indirect-CP-violation + CP-conserving
photon-mediated contributions. A faithful prediction
requires ``s -> d e+ e-`` Wilson matching, CKM inputs, SM--NP
interference conventions, the K_S -> pi0 e+ e- indirect-CP input, and
photon-mediated CP-conserving amplitudes. The current physics core has
no such machinery, so this constraint records the experimental anchor
and returns a non-vetoing deferred result.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K008.yaml`` (schema-flex anchor:
flat ``pdg_or_equivalent`` block with value ``"<2.8e-10"``).

Physics core
------------
No Delta S = 1 rare-kaon electron-mode implementation exists yet. The
placeholder reaches
:mod:`flavor_catalog_constraints.physics_adapters.rare_kaon_decays`,
whose result carries ``deferred=True``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from flavor_catalog_constraints.anchors import load_raw
from flavor_catalog_constraints.base import (
    ConstraintResult,
    ParameterPoint,
    Severity,
)
from flavor_catalog_constraints.physics_adapters.rare_kaon_decays import (
    evaluate_klong_pi0ee_placeholder,
)
from flavor_catalog_constraints.registry import register_constraint

_SM_BRANCHING_RATIO_ORDER = 3.0e-11


@dataclass(frozen=True)
class _Anchor:
    """Typed view over ``K008.yaml``'s flat ``pdg_or_equivalent`` block."""

    year: int
    source: str
    observable: str
    raw_value: str
    upper_limit: float
    uncertainty: Any
    units: str
    confidence_level: str
    source_url: str
    access_date: str
    snapshot_path: str
    sha256_of_text_snapshot: str


def _parse_upper_limit(raw_value: str) -> float:
    """Parse catalog upper-limit strings like ``<2.8e-10``."""
    text = str(raw_value).strip()
    if not text.startswith("<"):
        raise ValueError(f"Expected upper-limit string starting with '<': {text!r}")
    return float(text[1:].strip())


def _build_anchor(raw: dict[str, Any]) -> _Anchor:
    """K008 uses a flat ``pdg_or_equivalent`` upper-limit anchor."""
    raw_value = str(raw["value"])
    return _Anchor(
        year=int(raw["year"]),
        source=str(raw["source"]),
        observable=str(raw["observable"]),
        raw_value=raw_value,
        upper_limit=_parse_upper_limit(raw_value),
        uncertainty=raw["uncertainty"],
        units=str(raw["units"]),
        confidence_level=str(raw["confidence_level"]),
        source_url=str(raw["source_url"]),
        access_date=str(raw["access_date"]),
        snapshot_path=str(raw["snapshot_path"]),
        sha256_of_text_snapshot=str(raw["sha256_of_text_snapshot"]),
    )


@register_constraint
class Constraint:
    """Catalogued K_L -> pi0 e+ e- constraint (process_id ``K008``)."""

    process_id = "K008"
    severity = Severity.INFO
    observable = "BR(K_L -> pi0 e+ e-)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.snapshot_path,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_klong_pi0ee_placeholder()
        return ConstraintResult(
            process_id=self.process_id,
            passes=bool(result.passes),
            predicted=result.predicted,
            sm_prediction=None,
            experimental=self.anchor.upper_limit,
            ratio=None,
            budget=None,
            severity=self.severity,
            notes=(
                "Evaluation deferred: K008 is a Delta S=1 rare kaon "
                "electron mode with short-distance + indirect-CP-violation "
                "+ CP-conserving photon-mediated contributions. It must "
                "not reuse the Delta F=2 K001/K002 mixing machinery, the "
                "K003 eps'/eps placeholder, or the K004/K005 K -> pi nu "
                "nubar branching-ratio machinery."
            ),
            diagnostics={
                "deferred": bool(result.deferred),
                "evaluation_deferred": result.reason,
                "confidence_level": self.anchor.confidence_level,
                "upper_limit": self.anchor.upper_limit,
                "sm_branching_ratio_order": _SM_BRANCHING_RATIO_ORDER,
                "sm_branching_ratio_order_source": (
                    "KTeV 2004 arXiv hep-ex/0309072 introduction"
                ),
                "rate_components": (
                    "short-distance",
                    "indirect-CP-violation",
                    "CP-conserving photon-mediated",
                ),
                "component_summary": (
                    "short-distance + indirect-CP-violation + "
                    "CP-conserving photon-mediated"
                ),
                "required_machinery": (
                    "s -> d e+ e- Wilson matching, CKM elements, SM-NP "
                    "interference conventions, K_S -> pi0 e+ e- indirect-CP "
                    "input, and CP-conserving photon-mediated amplitudes"
                ),
            },
        )
