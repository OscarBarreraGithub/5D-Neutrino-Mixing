"""E002 — Muon electric dipole moment.

Severity rationale
------------------
INFO placeholder. EDM observables require CP-violating dipole Wilson
matching, 5D-RS one-loop contributions, and (for hadronic/atomic
species) hadronic matrix elements that are not yet in the physics core.

Catalog sidecar
---------------
``flavor_catalog/processes/edm_neutrino/E002.yaml`` (schema-flex anchor:
``pdg_or_equivalent.canonical_direct_limit`` block).
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
from flavor_catalog_constraints.physics_adapters.edm import evaluate_lepton_edm_placeholder
from flavor_catalog_constraints.registry import register_constraint


@dataclass(frozen=True)
class _Anchor:
    """Typed view over the ``E002.yaml`` anchor block (schema-flex)."""

    raw: dict[str, Any]
    observable: str
    source: str
    year: int | None


def _build_anchor(raw) -> _Anchor:
    """Build a minimal typed view; full schema varies across EDM yamls."""
    block = raw.get("canonical_direct_limit", raw)
    if isinstance(block, dict):
        # Find the first nested dict with an "observable" or "source" field
        if "observable" in block or "source" in block:
            inner = block
        else:
            inner = next(iter(block.values())) if block else {}
    else:
        inner = raw
    return _Anchor(
        raw=raw,
        observable=str(inner.get("observable", "Muon EDM |d_mu|")),
        source=str(inner.get("source", "PDG / lab publication")),
        year=int(inner["year"]) if isinstance(inner.get("year"), (int, str)) and str(inner.get("year")).isdigit() else None,
    )


@register_constraint
class Constraint:
    """Catalogued Muon electric dipole moment constraint (process_id ``E002``)."""

    process_id = "E002"
    severity = Severity.INFO
    observable = "Muon EDM |d_mu|"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="edm_neutrino", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.source,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_lepton_edm_placeholder("muon")
        return ConstraintResult(
            process_id=self.process_id,
            passes=bool(result.passes),
            predicted=result.predicted,
            sm_prediction=None,
            experimental=None,
            ratio=None,
            budget=None,
            severity=self.severity,
            notes=(
                "Evaluation deferred: E002 is an EDM observable requiring "
                "CP-violating dipole Wilson matching and dedicated 5D RS "
                "loop machinery not yet in the physics core."
            ),
            diagnostics={
                "deferred": bool(result.deferred),
                "evaluation_deferred": result.reason,
            },
        )
