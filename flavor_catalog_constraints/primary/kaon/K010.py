"""K010 - Short-lived neutral-kaon rare electron mode ``K_S -> pi0 e+ e-``.

Severity rationale
------------------
INFO placeholder. This observable is a Delta S = 1 rare kaon decay, not
a Delta F = 2 kaon-mixing constraint like K001/K002, not the
nonleptonic eps'/eps observable in K003, and not a K -> pi nu nubar
mode like K004/K005. The measured branching fraction is a short-lived
neutral-kaon electron-mode observable. A faithful prediction requires
``s -> d e+ e-`` Wilson
matching, CKM inputs, SM--NP interference conventions, form-factor and
phase-space treatment, and branching-ratio machinery. The current
physics core has no such machinery, so this constraint records the
experimental anchor and returns a non-vetoing deferred result.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K010.yaml`` (schema-flex anchor:
``pdg_or_equivalent.extrapolated_total_rate`` with value ``5.8e-9``).

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
    evaluate_kshort_pi0ee_placeholder,
)
from flavor_catalog_constraints.registry import register_constraint

_SM_BRANCHING_RATIO_ORDER = 5.8e-9


@dataclass(frozen=True)
class _Anchor:
    """Typed view over ``K010.yaml``'s nested ``pdg_or_equivalent`` block."""

    year: int
    source: str
    observable: str
    branching_ratio: float
    uncertainty: Any
    units: str
    assumptions: str
    confidence_level: str | None
    source_url: str
    access_date: str
    snapshot_path: str
    sha256_of_text_snapshot: str
    canonical_observable: str
    canonical_branching_ratio: float
    canonical_uncertainty: Any


def _build_anchor(raw: dict[str, Any]) -> _Anchor:
    """K010 uses nested canonical and extrapolated branching-ratio anchors."""
    canonical = raw["canonical_measurement"]
    extrapolated = raw["extrapolated_total_rate"]
    return _Anchor(
        year=int(extrapolated["year"]),
        source=str(raw["source"]),
        observable=str(extrapolated["observable"]),
        branching_ratio=float(extrapolated["value"]),
        uncertainty=extrapolated["uncertainty"],
        units=str(extrapolated["units"]),
        assumptions=str(extrapolated["assumptions"]),
        confidence_level=canonical["confidence_level"],
        source_url=str(extrapolated["source_url"]),
        access_date=str(extrapolated["access_date"]),
        snapshot_path=str(extrapolated["snapshot_path"]),
        sha256_of_text_snapshot=str(extrapolated["sha256_of_text_snapshot"]),
        canonical_observable=str(canonical["observable"]),
        canonical_branching_ratio=float(canonical["value"]),
        canonical_uncertainty=canonical["uncertainty"],
    )


@register_constraint
class Constraint:
    """Catalogued K_S -> pi0 e+ e- constraint (process_id ``K010``)."""

    process_id = "K010"
    severity = Severity.INFO
    observable = "BR(K_S -> pi0 e+ e-)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.snapshot_path,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_kshort_pi0ee_placeholder()
        return ConstraintResult(
            process_id=self.process_id,
            passes=bool(result.passes),
            predicted=result.predicted,
            sm_prediction=None,
            experimental=self.anchor.branching_ratio,
            ratio=None,
            budget=None,
            severity=self.severity,
            notes=(
                "Evaluation deferred: K010 is a Delta S=1 rare kaon "
                "electron mode for K_S -> pi0 e+ e-. It must not reuse "
                "the Delta F=2 K001/K002 mixing machinery, the K003 "
                "eps'/eps placeholder, or the K004/K005 K -> pi nu nubar "
                "branching-ratio machinery."
            ),
            diagnostics={
                "deferred": bool(result.deferred),
                "evaluation_deferred": result.reason,
                "confidence_level": self.anchor.confidence_level,
                "branching_ratio": self.anchor.branching_ratio,
                "canonical_partial_region_branching_ratio": (
                    self.anchor.canonical_branching_ratio
                ),
                "canonical_partial_region_observable": (
                    self.anchor.canonical_observable
                ),
                "extrapolation_assumptions": self.anchor.assumptions,
                "sm_branching_ratio_order": _SM_BRANCHING_RATIO_ORDER,
                "sm_branching_ratio_order_source": (
                    "PDG 2026 S012.10 full-region extrapolation"
                ),
                "rate_components": (
                    "short-distance",
                    "long-distance form-factor",
                    "phase-space extrapolation",
                ),
                "component_summary": (
                    "short-distance + long-distance form-factor + "
                    "phase-space extrapolation"
                ),
                "required_machinery": (
                    "s -> d e+ e- Wilson matching, CKM elements, SM-NP "
                    "interference conventions, K_S -> pi0 e+ e- form-factor "
                    "treatment, and branching-ratio machinery"
                ),
            },
        )
