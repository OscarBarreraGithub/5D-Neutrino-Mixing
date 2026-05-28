"""K006 - Long-lived neutral-kaon dimuon decay ``K_L -> mu+ mu-``.

Severity rationale
------------------
INFO placeholder. This observable is a Delta S = 1 rare kaon decay, not
a Delta F = 2 kaon-mixing constraint like K001/K002, not the
nonleptonic eps'/eps observable in K003, and not a K -> pi nu nubar
mode like K004/K005. The measured branching fraction is long-distance
dominated by the two-photon intermediate state. A faithful
short-distance prediction requires ``s -> d mu+ mu-`` Wilson matching,
CKM inputs, SM--NP interference conventions, and either long-distance
two-photon treatment or a documented short-distance subtraction. The
current physics core has no such machinery, so this constraint records
the experimental anchor and returns a non-vetoing deferred result.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K006.yaml`` (schema-flex anchor:
flat ``pdg_or_equivalent`` block).

Physics core
------------
No Delta S = 1 rare-kaon dimuon implementation exists yet. The
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
    evaluate_klong_mumu_placeholder,
)
from flavor_catalog_constraints.registry import register_constraint

_SHORT_DISTANCE_SM_REFERENCE = 9.0e-10


@dataclass(frozen=True)
class _Anchor:
    """Typed view over ``K006.yaml``'s flat ``pdg_or_equivalent`` block."""

    year: int
    source: str
    observable: str
    value: float
    uncertainty: float
    units: str
    confidence_level: Any
    source_url: str
    access_date: str
    snapshot_path: str
    sha256_of_text_snapshot: str


def _build_anchor(raw) -> _Anchor:
    """K006 uses a flat ``pdg_or_equivalent`` anchor."""
    return _Anchor(
        year=int(raw["year"]),
        source=str(raw["source"]),
        observable=str(raw["observable"]),
        value=float(raw["value"]),
        uncertainty=float(raw["uncertainty"]),
        units=str(raw["units"]),
        confidence_level=raw["confidence_level"],
        source_url=str(raw["source_url"]),
        access_date=str(raw["access_date"]),
        snapshot_path=str(raw["snapshot_path"]),
        sha256_of_text_snapshot=str(raw["sha256_of_text_snapshot"]),
    )


@register_constraint
class Constraint:
    """Catalogued K_L -> mu+ mu- constraint (process_id ``K006``)."""

    process_id = "K006"
    severity = Severity.INFO
    observable = "BR(K_L -> mu+ mu-)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.snapshot_path,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_klong_mumu_placeholder()
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
                "Evaluation deferred: K006 is a Delta S=1 rare kaon "
                "dimuon decay whose total rate is long-distance dominated "
                "by the two-photon intermediate state. It must not reuse "
                "the Delta F=2 K001/K002 mixing machinery, the K003 "
                "eps'/eps placeholder, or the K004/K005 K -> pi nu nubar "
                "branching-ratio machinery."
            ),
            diagnostics={
                "deferred": bool(result.deferred),
                "evaluation_deferred": result.reason,
                "short_distance_sm_reference": _SHORT_DISTANCE_SM_REFERENCE,
                "short_distance_sm_reference_source": "Buras et al.",
                "long_distance_dominated": True,
                "dominant_long_distance_state": "two-photon",
                "required_machinery": (
                    "s -> d mu+ mu- Wilson matching, CKM elements, "
                    "SM-NP interference conventions, and long-distance "
                    "two-photon treatment or documented short-distance "
                    "subtraction"
                ),
            },
        )
