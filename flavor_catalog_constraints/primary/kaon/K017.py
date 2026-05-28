"""K017 - Charged-kaon leptonic LFU ratio ``R_K``.

Severity rationale
------------------
INFO placeholder. ``R_K = Gamma(K+ -> e+ nu) / Gamma(K+ -> mu+ nu)``
is a tree-level Standard Model charged-current lepton-universality
ratio, not a Delta S = 1 rare kaon decay. Its SM prediction is very
precise once radiative QED corrections are included. New physics would
enter through modified W-lepton couplings in the 5D RS model or through
new tree-level four-fermion contact operators. The current physics core
has no such charged-current LFU machinery, so this constraint records
the experimental and SM anchors and returns a non-vetoing deferred
result.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K017.yaml`` (schema-flex anchor:
nested ``pdg_or_equivalent.canonical_ratio`` and ``sm_prediction``
blocks).

Physics core
------------
No charged-kaon LFU ratio implementation exists yet. The placeholder
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
    evaluate_rk_lfu_placeholder,
)
from flavor_catalog_constraints.registry import register_constraint


@dataclass(frozen=True)
class _Anchor:
    """Typed view over selected ``K017.yaml`` ``pdg_or_equivalent`` records."""

    year: int
    source: str
    edition: str
    citation: str
    raw_value: str
    raw_uncertainty: str
    scale: float
    value: float
    uncertainty: float
    units: str
    value_summary: str
    pdg_identifier: str
    source_url: str
    access_date: str
    snapshot_path: str
    sha256_of_local_snapshot: str
    sm_source: str
    sm_arxiv_id: str
    sm_year: int
    sm_raw_value: str
    sm_raw_uncertainty: str
    sm_scale: float
    sm_prediction: float
    sm_uncertainty: float
    sm_units: str
    sm_value_summary: str
    sm_source_url: str
    sm_access_date: str
    sm_snapshot_path: str
    sm_sha256_of_local_snapshot: str


def _scaled(record: dict[str, Any], key: str) -> float:
    """Return a string-valued numeric field scaled by ``record['scale']``."""
    return float(record[key]) * float(record["scale"])


def _build_anchor(raw: dict[str, Any]) -> _Anchor:
    """K017 uses nested experimental-ratio and SM-prediction records."""
    canonical = raw["canonical_ratio"]
    sm_prediction = raw["sm_prediction"]
    return _Anchor(
        year=int(canonical["year"]),
        source=str(canonical["source"]),
        edition=str(canonical["edition"]),
        citation=str(canonical["citation"]),
        raw_value=str(canonical["value"]),
        raw_uncertainty=str(canonical["uncertainty"]),
        scale=float(canonical["scale"]),
        value=_scaled(canonical, "value"),
        uncertainty=_scaled(canonical, "uncertainty"),
        units=str(canonical["units"]),
        value_summary=str(canonical["value_summary"]),
        pdg_identifier=str(canonical["pdg_identifier"]),
        source_url=str(canonical["source_url"]),
        access_date=str(canonical["access_date"]),
        snapshot_path=str(canonical["snapshot_path"]),
        sha256_of_local_snapshot=str(canonical["sha256_of_local_snapshot"]),
        sm_source=str(sm_prediction["source"]),
        sm_arxiv_id=str(sm_prediction["arxiv_id"]),
        sm_year=int(sm_prediction["year"]),
        sm_raw_value=str(sm_prediction["value"]),
        sm_raw_uncertainty=str(sm_prediction["uncertainty"]),
        sm_scale=float(sm_prediction["scale"]),
        sm_prediction=_scaled(sm_prediction, "value"),
        sm_uncertainty=_scaled(sm_prediction, "uncertainty"),
        sm_units=str(sm_prediction["units"]),
        sm_value_summary=str(sm_prediction["value_summary"]),
        sm_source_url=str(sm_prediction["source_url"]),
        sm_access_date=str(sm_prediction["access_date"]),
        sm_snapshot_path=str(sm_prediction["snapshot_path"]),
        sm_sha256_of_local_snapshot=str(
            sm_prediction["sha256_of_local_snapshot"]
        ),
    )


@register_constraint
class Constraint:
    """Catalogued charged-kaon LFU ratio constraint (process_id ``K017``)."""

    process_id = "K017"
    severity = Severity.INFO
    observable = "R_K = Gamma(K+ -> e+ nu) / Gamma(K+ -> mu+ nu)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (
            self.anchor.snapshot_path,
            self.anchor.sm_snapshot_path,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_rk_lfu_placeholder()
        return ConstraintResult(
            process_id=self.process_id,
            passes=bool(result.passes),
            predicted=result.predicted,
            sm_prediction=self.anchor.sm_prediction,
            experimental=self.anchor.value,
            ratio=None,
            budget=None,
            severity=self.severity,
            notes=(
                "Evaluation deferred: K017 is the charged-kaon "
                "lepton-universality ratio R_K, not a Delta S=1 rare "
                "kaon decay. It requires a charged-current treatment "
                "of modified W-lepton couplings or new tree-level "
                "four-fermion contact operators, with radiative QED "
                "corrections matched to the inclusive experimental "
                "definition."
            ),
            diagnostics={
                "deferred": bool(result.deferred),
                "evaluation_deferred": result.reason,
                "standard_notation": (
                    "R_K = Gamma(K+ -> e+ nu_e(gamma)) / "
                    "Gamma(K+ -> mu+ nu_mu(gamma))"
                ),
                "tree_level_sm_ratio": True,
                "not_delta_s1_rare_decay": True,
                "pdg_identifier": self.anchor.pdg_identifier,
                "experimental_uncertainty": self.anchor.uncertainty,
                "sm_uncertainty": self.anchor.sm_uncertainty,
                "radiative_qed_corrections_required": True,
                "np_entry": (
                    "modified W-lepton couplings or new tree-level "
                    "four-fermion contact operators"
                ),
                "required_machinery": (
                    "modified W-lepton coupling treatment in 5D RS, "
                    "new tree-level four-fermion contact-operator "
                    "matching, and radiative QED corrections"
                ),
            },
        )
