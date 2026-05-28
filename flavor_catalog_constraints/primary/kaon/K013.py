"""K013 - Radiative neutral-kaon decay ``K_L -> pi0 gamma gamma``.

Severity rationale
------------------
INFO placeholder. This observable is a radiative Delta S = 1 kaon
decay, not a Delta F = 2 kaon-mixing constraint like K001/K002, not
the nonleptonic eps'/eps observable in K003, and not a K -> pi nu
nubar mode like K004/K005. The Standard Model rate is dominated by
chiral perturbation theory O(p^4) plus O(p^6) contributions and photon
matrix elements, with BR_SM ~ 1.4e-6. The PDG experimental branching
fraction is about 1.27e-6. The current physics core has no radiative
ChPT amplitude machinery, so this constraint records the experimental
anchor and returns a non-vetoing deferred result.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K013.yaml`` (schema-flex anchor:
flat ``pdg_or_equivalent`` block).

Physics core
------------
No radiative ChPT rare-kaon implementation exists yet. The placeholder
reaches
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
    evaluate_klong_pi0gg_placeholder,
)
from flavor_catalog_constraints.registry import register_constraint

_SM_BRANCHING_RATIO_REFERENCE = 1.4e-6


@dataclass(frozen=True)
class _Anchor:
    """Typed view over ``K013.yaml``'s flat ``pdg_or_equivalent`` block."""

    year: int
    source: str
    observable: str
    value: float
    uncertainty: float
    units: str
    display: str
    source_url: str
    access_date: str
    snapshot_path: str
    sha256_of_text_snapshot: str


def _build_anchor(raw: dict[str, Any]) -> _Anchor:
    """K013 uses a flat ``pdg_or_equivalent`` anchor."""
    return _Anchor(
        year=int(raw["year"]),
        source=str(raw["source"]),
        observable=str(raw["observable"]),
        value=float(raw["value"]),
        uncertainty=float(raw["uncertainty"]),
        units=str(raw["units"]),
        display=str(raw["display"]),
        source_url=str(raw["source_url"]),
        access_date=str(raw["access_date"]),
        snapshot_path=str(raw["snapshot_path"]),
        sha256_of_text_snapshot=str(raw["sha256"]),
    )


@register_constraint
class Constraint:
    """Catalogued K_L -> pi0 gamma gamma constraint (process_id ``K013``)."""

    process_id = "K013"
    severity = Severity.INFO
    observable = "BR(K_L -> pi0 gamma gamma)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.snapshot_path,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_klong_pi0gg_placeholder()
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
                "Evaluation deferred: K013 is a radiative kaon decay "
                "whose SM rate is dominated by chiral perturbation theory "
                "O(p^4)+O(p^6) contributions and photon matrix elements. "
                "It must not reuse the Delta F=2 K001/K002 mixing "
                "machinery, the K003 eps'/eps placeholder, or a Delta "
                "S=1 Wilson-matching-only treatment."
            ),
            diagnostics={
                "deferred": bool(result.deferred),
                "evaluation_deferred": result.reason,
                "standard_notation": "K_L -> pi0 gamma gamma",
                "sm_branching_ratio_reference": _SM_BRANCHING_RATIO_REFERENCE,
                "sm_branching_ratio_reference_source": (
                    "K013 placeholder task; radiative-kaon ChPT benchmark"
                ),
                "radiative_decay": True,
                "dominant_sm_theory": (
                    "chiral perturbation theory O(p^4)+O(p^6)"
                ),
                "requires_photon_matrix_element": True,
                "required_machinery": (
                    "chiral perturbation theory O(p^4)+O(p^6) amplitudes, "
                    "photon matrix elements, and branching-ratio machinery; "
                    "Delta S=1 Wilson matching alone is insufficient"
                ),
            },
        )
