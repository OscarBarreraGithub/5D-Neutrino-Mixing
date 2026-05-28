"""T005 - BR(t->cg) (placeholder)."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Any

from flavor_catalog_constraints.anchors import load_raw
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.top_higgs_ew import evaluate_top_higgs_ew_placeholder
from flavor_catalog_constraints.registry import register_constraint


@dataclass(frozen=True)
class _Anchor:
    raw: dict[str, Any]
    observable: str
    source: str


def _build_anchor(raw) -> _Anchor:
    inner = raw
    if isinstance(raw, dict):
        for v in raw.values():
            if isinstance(v, dict) and ("observable" in v or "source" in v):
                inner = v; break
    obs = str(inner.get("observable", "BR(t->cg)")) if isinstance(inner, dict) else "BR(t->cg)"
    src = str(inner.get("source", "PDG / ATLAS / CMS")) if isinstance(inner, dict) else "PDG / ATLAS / CMS"
    return _Anchor(raw=raw, observable=obs, source=src)


@register_constraint
class Constraint:
    process_id = "T005"
    severity = Severity.INFO
    observable = "BR(t->cg)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="top_higgs_ew", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.source,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_top_higgs_ew_placeholder("top FCNC gluonic")
        return ConstraintResult(
            process_id=self.process_id, passes=bool(result.passes),
            predicted=result.predicted, sm_prediction=None, experimental=None,
            ratio=None, budget=None, severity=self.severity,
            notes="Evaluation deferred: T005 requires top/EW/LFV machinery not in core.",
            diagnostics={"deferred": bool(result.deferred),
                          "evaluation_deferred": result.reason},
        )
