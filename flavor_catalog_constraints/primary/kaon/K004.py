"""K004 - Charged rare kaon decay ``K+ -> pi+ nu nubar``.

Severity rationale
------------------
INFO placeholder. This observable is not a Delta F = 2 kaon-mixing
constraint like K001/K002, and it is not the nonleptonic eps'/eps
observable in K003. A faithful prediction requires Delta S = 1
``s -> d nu nubar`` Wilson matching, CKM inputs, SM--NP interference
conventions, and the charged-mode branching-ratio formula. The current
physics core has no such machinery, so this constraint records the
experimental anchor and returns a non-vetoing deferred result.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K004.yaml`` (schema-flex anchor:
nested ``pdg_or_equivalent.latest_experimental_value`` block).

Physics core
------------
No Delta S = 1 rare-kaon implementation exists yet. The placeholder
reaches :mod:`flavor_catalog_constraints.physics_adapters.rare_kaon_decays`,
whose result carries ``deferred=True``.
"""

from __future__ import annotations

from dataclasses import dataclass

from flavor_catalog_constraints.anchors import load_raw
from flavor_catalog_constraints.base import (
    ConstraintResult,
    ParameterPoint,
    Severity,
)
from flavor_catalog_constraints.physics_adapters.rare_kaon_decays import (
    evaluate_kplus_piplus_nunu_placeholder,
)
from flavor_catalog_constraints.registry import register_constraint


@dataclass(frozen=True)
class _Anchor:
    """Typed view over ``K004.yaml``'s latest experimental-value block."""

    year: int
    source: str
    raw_value: str
    raw_uncertainty: str
    scale: float
    value: float
    uncertainty_up: float
    uncertainty_down: float
    units: str
    dataset: str
    value_summary: str
    source_url: str
    access_date: str
    snapshot_path: str
    sha256_of_local_snapshot: str


def _parse_asymmetric_uncertainty(raw: str, scale: float) -> tuple[float, float]:
    """Parse strings like ``+1.9/-1.8`` into positive scaled errors."""
    up_raw, down_raw = str(raw).split("/")
    return abs(float(up_raw)) * scale, abs(float(down_raw)) * scale


def _build_anchor(raw) -> _Anchor:
    """K004 uses ``latest_experimental_value`` nested under ``pdg_or_equivalent``."""
    latest = raw["latest_experimental_value"]
    scale = float(latest["scale"])
    value = float(latest["value"]) * scale
    uncertainty_up, uncertainty_down = _parse_asymmetric_uncertainty(
        latest["uncertainty"], scale
    )
    return _Anchor(
        year=int(latest["year"]),
        source=str(latest["source"]),
        raw_value=str(latest["value"]),
        raw_uncertainty=str(latest["uncertainty"]),
        scale=scale,
        value=value,
        uncertainty_up=uncertainty_up,
        uncertainty_down=uncertainty_down,
        units=str(latest["units"]),
        dataset=str(latest["dataset"]),
        value_summary=str(latest["value_summary"]),
        source_url=str(latest["source_url"]),
        access_date=str(latest["access_date"]),
        snapshot_path=str(latest["snapshot_path"]),
        sha256_of_local_snapshot=str(latest["sha256_of_local_snapshot"]),
    )


@register_constraint
class Constraint:
    """Catalogued K+ -> pi+ nu nubar constraint (process_id ``K004``)."""

    process_id = "K004"
    severity = Severity.INFO
    observable = "BR(K+ -> pi+ nu nubar)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.snapshot_path,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_kplus_piplus_nunu_placeholder()
        return ConstraintResult(
            process_id=self.process_id,
            passes=bool(result.passes),
            predicted=result.predicted,
            sm_prediction=None,
            experimental=self.anchor.value,
            ratio=None,
            budget=None,
            severity=self.severity,
            notes=(
                "Evaluation deferred: K004 is a Delta S=1 rare kaon "
                "decay and must not reuse the Delta F=2 K001/K002 "
                "mixing machinery or the K003 eps'/eps placeholder."
            ),
            diagnostics={
                "deferred": bool(result.deferred),
                "evaluation_deferred": result.reason,
                "required_machinery": (
                    "K -> pi nu nubar Wilson matching, BMU-operator "
                    "mapping, CKM elements, SM-NP interference "
                    "conventions, and charged-mode branching-ratio "
                    "evaluation"
                ),
            },
        )
