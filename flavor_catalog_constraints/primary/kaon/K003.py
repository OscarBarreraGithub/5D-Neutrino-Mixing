"""K003 - Direct CP violation in K -> pi pi (``Re(epsilon'/epsilon)``).

Severity rationale
------------------
INFO placeholder. This observable is not a Delta F = 2 kaon-mixing
constraint like K001/K002. A faithful NP prediction requires Delta S = 1
Wilson coefficients, QCD/electroweak-penguin running, and K -> pi pi
hadronic matrix elements. The current physics core has no such
machinery, so this constraint records the experimental anchor and
returns a non-vetoing deferred result.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K003.yaml`` (schema-flex anchor:
flat ``pdg_or_equivalent`` block).

Physics core
------------
No Delta S = 1 implementation exists yet. The placeholder reaches
:mod:`flavor_catalog_constraints.physics_adapters.eps_prime`, whose
result carries ``deferred=True``.
"""

from __future__ import annotations

from dataclasses import dataclass

from flavor_catalog_constraints.anchors import load_raw
from flavor_catalog_constraints.base import (
    ConstraintResult,
    ParameterPoint,
    Severity,
)
from flavor_catalog_constraints.physics_adapters.eps_prime import (
    evaluate_eps_prime_placeholder,
)
from flavor_catalog_constraints.registry import register_constraint


@dataclass(frozen=True)
class _Anchor:
    """Typed view over ``K003.yaml``'s flat ``pdg_or_equivalent`` block."""

    year: int
    source: str
    observable: str
    value: float
    uncertainty: float
    units: str
    display_value: str
    source_url: str
    access_date: str
    snapshot_path: str
    sha256_of_text_snapshot: str


def _build_anchor(raw) -> _Anchor:
    """K003 uses a flat ``pdg_or_equivalent`` anchor."""
    return _Anchor(
        year=int(raw["year"]),
        source=str(raw["source"]),
        observable=str(raw["observable"]),
        value=float(raw["value"]),
        uncertainty=float(raw["uncertainty"]),
        units=str(raw["units"]),
        display_value=str(raw["display_value"]),
        source_url=str(raw["source_url"]),
        access_date=str(raw["access_date"]),
        snapshot_path=str(raw["snapshot_path"]),
        sha256_of_text_snapshot=str(raw["sha256_of_text_snapshot"]),
    )


@register_constraint
class Constraint:
    """Catalogued Re(epsilon'/epsilon) constraint (process_id ``K003``)."""

    process_id = "K003"
    severity = Severity.INFO
    observable = "Re(epsilon'/epsilon)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.snapshot_path,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_eps_prime_placeholder()
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
                "Evaluation deferred: K003 is a Delta S=1 K -> pi pi "
                "observable and must not reuse the Delta F=2 K001/K002 "
                "mixing machinery."
            ),
            diagnostics={
                "deferred": bool(result.deferred),
                "evaluation_deferred": result.reason,
                "required_machinery": (
                    "Delta S=1 RG, QCD/electroweak-penguin Wilson "
                    "coefficients, K -> pi pi matrix elements, and RS NP "
                    "matching"
                ),
            },
        )
