"""K018 - Kaon semileptonic ``K_l3`` extraction of ``|V_us| f_+(0)``.

Severity rationale
------------------
INFO placeholder. ``K_l3`` decays are tree-level Standard Model
charged-current processes that provide the source-level
``|V_us| f_+(0)`` input. They are not Delta S = 1 rare kaon decays and
must not reuse the neutral-kaon Delta F = 2 machinery. New physics would
enter through modified W-quark-quark couplings in the 5D RS model or new
tree-level four-fermion contact operators. The current physics core has
generic CKM-fit targets but no source-level ``K_l3`` extraction,
radiative/isospin correction, or BSM charged-current matching machinery,
so this constraint records the catalog anchor and returns a non-vetoing
deferred result.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K018.yaml`` (schema-flex anchor:
list-valued ``pdg_or_equivalent`` block, using
``PDG2025:K018:average_fplus_vus``).

Physics core
------------
No source-level semileptonic-kaon implementation exists yet. The
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
    evaluate_kl3_vus_placeholder,
)
from flavor_catalog_constraints.registry import register_constraint

_ANCHOR_VALUE_ID = "PDG2025:K018:average_fplus_vus"


@dataclass(frozen=True)
class _Anchor:
    """Typed view over the selected ``K018.yaml`` ``pdg_or_equivalent`` record."""

    value_id: str
    observable: str
    year: int
    value: float
    uncertainty: float
    units: str
    source: str
    source_url: str
    access_date: str
    snapshot_path: str
    sha256_of_local_snapshot: str
    value_summary: str


def _find_record(raw: list[dict[str, Any]], value_id: str) -> dict[str, Any]:
    for record in raw:
        if record.get("value_id") == value_id:
            return record
    raise KeyError(f"K018 anchor record {value_id!r} not found")


def _build_anchor(raw: list[dict[str, Any]]) -> _Anchor:
    """K018 uses the PDG ``|V_us| f_+(0)`` average record."""
    record = _find_record(raw, _ANCHOR_VALUE_ID)
    return _Anchor(
        value_id=str(record["value_id"]),
        observable=str(record["observable"]),
        year=int(record["year"]),
        value=float(record["value"]),
        uncertainty=float(record["uncertainty"]),
        units=str(record["units"]),
        source=str(record["source"]),
        source_url=str(record["source_url"]),
        access_date=str(record["access_date"]),
        snapshot_path=str(record["snapshot_path"]),
        sha256_of_local_snapshot=str(record["sha256"]),
        value_summary=str(record["value_summary"]),
    )


@register_constraint
class Constraint:
    """Catalogued K_l3 ``|V_us| f_+(0)`` constraint (process_id ``K018``)."""

    process_id = "K018"
    severity = Severity.INFO
    observable = "K_l3 |V_us| f_+(0)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.snapshot_path,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_kl3_vus_placeholder()
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
                "Evaluation deferred: K018 is the tree-level charged-current "
                "K_l3 extraction of |V_us| f_+(0), not a Delta S=1 rare "
                "kaon decay. It requires source-level |V_us| extraction, "
                "5D RS effects on the W-quark-quark coupling, radiative "
                "corrections, and possible tree-level four-fermion contact "
                "operators."
            ),
            diagnostics={
                "deferred": bool(result.deferred),
                "evaluation_deferred": result.reason,
                "standard_notation": (
                    "K_l3: K -> pi e nu and K -> pi mu nu; |V_us| f_+(0)"
                ),
                "tree_level_sm_observable": True,
                "charged_current_semileptonic": True,
                "not_delta_s1_rare_decay": True,
                "value_id": self.anchor.value_id,
                "experimental_uncertainty": self.anchor.uncertainty,
                "vus_fplus_product": self.anchor.value,
                "radiative_corrections_required": True,
                "np_entry": (
                    "modified W-quark-quark coupling or new tree-level "
                    "four-fermion contact operators"
                ),
                "required_machinery": (
                    "|V_us| extraction, 5D RS effects on W-quark vertex, "
                    "radiative/isospin corrections, lattice f_+(0) "
                    "provenance, and four-fermion contact-operator matching"
                ),
            },
        )
