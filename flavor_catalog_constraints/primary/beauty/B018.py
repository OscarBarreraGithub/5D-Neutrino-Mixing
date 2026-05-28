"""B018 - LFU R_K = BR(B+->K+μμ)/BR(B+->K+ee) (placeholder)."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Any

from flavor_catalog_constraints.anchors import load_raw
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.beauty import evaluate_beauty_placeholder
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
                inner = v
                break
    obs = str(inner.get("observable", "LFU R_K = BR(B+->K+μμ)/BR(B+->K+ee)")) if isinstance(inner, dict) else "LFU R_K = BR(B+->K+μμ)/BR(B+->K+ee)"
    src = str(inner.get("source", "PDG / HFLAV")) if isinstance(inner, dict) else "PDG / HFLAV"
    return _Anchor(raw=raw, observable=obs, source=src)


@register_constraint
class Constraint:
    process_id = "B018"
    severity = Severity.INFO
    observable = "LFU R_K = BR(B+->K+μμ)/BR(B+->K+ee)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="beauty", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.source,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_beauty_placeholder("LFU in b->s ll")
        return ConstraintResult(
            process_id=self.process_id, passes=bool(result.passes),
            predicted=result.predicted, sm_prediction=None, experimental=None,
            ratio=None, budget=None, severity=self.severity,
            notes="Evaluation deferred: B018 requires B-sector machinery not in core.",
            diagnostics={"deferred": bool(result.deferred),
                          "evaluation_deferred": result.reason},
        )
