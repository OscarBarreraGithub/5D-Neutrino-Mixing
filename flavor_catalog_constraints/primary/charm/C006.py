"""C006 - LFV neutral-charm decay D0 -> e mu (placeholder)."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Any

from flavor_catalog_constraints.anchors import load_raw
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.charm import evaluate_charm_placeholder
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
    obs = str(inner.get("observable", "BR(D0->emu)")) if isinstance(inner, dict) else "BR(D0->emu)"
    src = str(inner.get("source", "PDG / HFLAV")) if isinstance(inner, dict) else "PDG / HFLAV"
    return _Anchor(raw=raw, observable=obs, source=src)


@register_constraint
class Constraint:
    process_id = "C006"
    severity = Severity.INFO
    observable = "BR(D0->emu)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="charm", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.source,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_charm_placeholder("ΔC=1 LFV")
        return ConstraintResult(
            process_id=self.process_id, passes=bool(result.passes),
            predicted=result.predicted, sm_prediction=None, experimental=None,
            ratio=None, budget=None, severity=self.severity,
            notes="Evaluation deferred: C006 requires charm-sector NP machinery not in core.",
            diagnostics={"deferred": bool(result.deferred),
                          "evaluation_deferred": result.reason},
        )
