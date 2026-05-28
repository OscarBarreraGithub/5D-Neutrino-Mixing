"""B008 - Rare tauonic Bd/Bs->ττ (SECONDARY placeholder)."""

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
                inner = v; break
            if isinstance(v, list) and v and isinstance(v[0], dict):
                inner = v[0]; break
    obs = str(inner.get("observable", "Rare tauonic Bd/Bs->ττ")) if isinstance(inner, dict) else "Rare tauonic Bd/Bs->ττ"
    src = str(inner.get("source", "PDG")) if isinstance(inner, dict) else "PDG"
    return _Anchor(raw=raw, observable=obs, source=src)


@register_constraint
class Constraint:
    process_id = "B008"
    severity = Severity.INFO
    observable = "Rare tauonic Bd/Bs->ττ"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="beauty", tier="secondary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.source,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_beauty_placeholder("ΔB=1 rare leptonic")
        return ConstraintResult(
            process_id=self.process_id, passes=bool(result.passes),
            predicted=result.predicted, sm_prediction=None, experimental=None,
            ratio=None, budget=None, severity=self.severity,
            notes="Evaluation deferred: B008 SECONDARY placeholder.",
            diagnostics={"deferred": bool(result.deferred),
                          "evaluation_deferred": result.reason},
        )
