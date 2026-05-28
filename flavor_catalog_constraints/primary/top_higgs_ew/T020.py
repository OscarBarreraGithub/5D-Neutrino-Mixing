"""T020 - BR(h->eμ) (placeholder)."""

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
    obs = str(inner.get("observable", "BR(h->eμ)")) if isinstance(inner, dict) else "BR(h->eμ)"
    src = str(inner.get("source", "PDG / ATLAS / CMS")) if isinstance(inner, dict) else "PDG / ATLAS / CMS"
    return _Anchor(raw=raw, observable=obs, source=src)


@register_constraint
class Constraint:
    process_id = "T020"
    severity = Severity.INFO
    observable = "BR(h->eμ)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="top_higgs_ew", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.source,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_top_higgs_ew_placeholder("LFV Higgs decay")
        return ConstraintResult(
            process_id=self.process_id, passes=bool(result.passes),
            predicted=result.predicted, sm_prediction=None, experimental=None,
            ratio=None, budget=None, severity=self.severity,
            notes="Evaluation deferred: T020 requires top/EW/LFV machinery not in core.",
            diagnostics={"deferred": bool(result.deferred),
                          "evaluation_deferred": result.reason},
        )
