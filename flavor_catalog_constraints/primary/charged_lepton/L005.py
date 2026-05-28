"""L005 - Muon-to-electron conversion in titanium (placeholder)."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Any

from flavor_catalog_constraints.anchors import load_raw
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.lepton import evaluate_mue_conversion_placeholder
from flavor_catalog_constraints.registry import register_constraint


@dataclass(frozen=True)
class _Anchor:
    raw: dict[str, Any]
    observable: str
    source: str


def _build_anchor(raw) -> _Anchor:
    # schema-flex: find the first nested dict with observable/source fields
    inner = raw
    if isinstance(raw, dict):
        for v in raw.values():
            if isinstance(v, dict) and ("observable" in v or "source" in v):
                inner = v
                break
    obs = str(inner.get("observable", "R_mue^Ti")) if isinstance(inner, dict) else "R_mue^Ti"
    src = str(inner.get("source", "PDG / lab publication")) if isinstance(inner, dict) else "PDG / lab publication"
    return _Anchor(raw=raw, observable=obs, source=src)


@register_constraint
class Constraint:
    process_id = "L005"
    severity = Severity.INFO
    observable = "R_mue^Ti"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="charged_lepton", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.source,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_mue_conversion_placeholder("Ti")
        return ConstraintResult(
            process_id=self.process_id, passes=bool(result.passes),
            predicted=result.predicted, sm_prediction=None, experimental=None,
            ratio=None, budget=None, severity=self.severity,
            notes="Evaluation deferred: L005 requires LFV/trident machinery not in core.",
            diagnostics={"deferred": bool(result.deferred),
                          "evaluation_deferred": result.reason},
        )
