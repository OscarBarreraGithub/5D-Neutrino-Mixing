"""K012 - Short-lived neutral-kaon dimuon decay ``K_S -> mu+ mu-``.

Severity rationale
------------------
INFO placeholder. This observable is a Delta S = 1 rare kaon decay, not
a Delta F = 2 kaon-mixing constraint like K001/K002, not the
nonleptonic eps'/eps observable in K003, and not a K -> pi nu nubar
mode like K004/K005. The experimental anchor is an upper limit on the
short-lived dimuon branching fraction. A faithful prediction requires
``s -> d mu+ mu-`` Wilson matching, CKM inputs, SM--NP interference
conventions, and K_S-specific long-distance and time-dependent
interference treatment. The current physics core has no such machinery,
so this constraint records the experimental anchor and returns a
non-vetoing deferred result.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K012.yaml`` (schema-flex anchor:
``pdg_or_equivalent`` list containing the headline upper limit and SM
estimate records).

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
    evaluate_kshort_mumu_placeholder,
)
from flavor_catalog_constraints.registry import register_constraint

_SM_BRANCHING_RATIO_REFERENCE = 5.1e-12


@dataclass(frozen=True)
class _Anchor:
    """Typed view over selected ``K012.yaml`` ``pdg_or_equivalent`` records."""

    year: int
    source: str
    observable: str
    raw_value: str
    upper_limit: float
    uncertainty: Any
    units: str
    confidence_level: Any
    source_url: str
    access_date: str
    snapshot_path: str
    sha256_of_text_snapshot: str
    sm_observable: str
    sm_raw_value: str
    catalog_sm_estimate: float
    sm_source: str
    sm_source_url: str


def _by_value_id(raw: list[dict[str, Any]], value_id: str) -> dict[str, Any]:
    """Return the K012 record with the requested ``value_id``."""
    for record in raw:
        if record.get("value_id") == value_id:
            return record
    raise KeyError(f"K012.yaml pdg_or_equivalent missing value_id {value_id!r}")


def _parse_upper_limit(raw_value: str) -> float:
    """Parse catalog upper-limit strings like ``<2.1e-10``."""
    text = str(raw_value).strip()
    if not text.startswith("<"):
        raise ValueError(f"Expected upper-limit string starting with '<': {text!r}")
    return float(text[1:].strip())


def _parse_approximate_value(raw_value: str) -> float:
    """Parse catalog strings like ``approximately 5e-12``."""
    text = str(raw_value).strip().lower()
    prefix = "approximately "
    if text.startswith(prefix):
        text = text[len(prefix) :]
    return float(text)


def _build_anchor(raw: list[dict[str, Any]]) -> _Anchor:
    """K012 uses list records for the limit and auxiliary SM estimate."""
    headline = _by_value_id(raw, "PDG2025:K012:headline_limit")
    sm_estimate = _by_value_id(raw, "Chobanova2018:K012:sm_estimate")
    raw_value = str(headline["value"])
    return _Anchor(
        year=int(headline["year"]),
        source=str(headline["source"]),
        observable=str(headline["observable"]),
        raw_value=raw_value,
        upper_limit=_parse_upper_limit(raw_value),
        uncertainty=headline["uncertainty"],
        units=str(headline["units"]),
        confidence_level=headline["confidence_level"],
        source_url=str(headline["source_url"]),
        access_date=str(headline["access_date"]),
        snapshot_path=str(headline["snapshot_path"]),
        sha256_of_text_snapshot=str(headline["sha256"]),
        sm_observable=str(sm_estimate["observable"]),
        sm_raw_value=str(sm_estimate["value"]),
        catalog_sm_estimate=_parse_approximate_value(str(sm_estimate["value"])),
        sm_source=str(sm_estimate["source"]),
        sm_source_url=str(sm_estimate["source_url"]),
    )


@register_constraint
class Constraint:
    """Catalogued K_S -> mu+ mu- constraint (process_id ``K012``)."""

    process_id = "K012"
    severity = Severity.INFO
    observable = "BR(K_S -> mu+ mu-)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.snapshot_path,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_kshort_mumu_placeholder()
        return ConstraintResult(
            process_id=self.process_id,
            passes=bool(result.passes),
            predicted=result.predicted,
            sm_prediction=None,
            experimental=self.anchor.upper_limit,
            ratio=None,
            budget=None,
            severity=self.severity,
            notes=(
                "Evaluation deferred: K012 is a Delta S=1 rare kaon "
                "dimuon decay for K_S -> mu+ mu-. It must not reuse the "
                "Delta F=2 K001/K002 mixing machinery, the K003 eps'/eps "
                "placeholder, or the K004/K005 K -> pi nu nubar "
                "branching-ratio machinery."
            ),
            diagnostics={
                "deferred": bool(result.deferred),
                "evaluation_deferred": result.reason,
                "confidence_level": self.anchor.confidence_level,
                "upper_limit": self.anchor.upper_limit,
                "sm_branching_ratio_reference": _SM_BRANCHING_RATIO_REFERENCE,
                "sm_branching_ratio_reference_source": (
                    "K_S -> mu+ mu- SM order quoted in K012 implementation task"
                ),
                "catalog_sm_estimate": self.anchor.catalog_sm_estimate,
                "catalog_sm_estimate_source": self.anchor.sm_source,
                "long_distance_dominated": True,
                "dominant_long_distance_state": "two-photon",
                "required_machinery": (
                    "s -> d mu+ mu- Wilson matching, CKM elements, "
                    "SM-NP interference conventions, K_S-specific "
                    "long-distance treatment, time-dependent neutral-kaon "
                    "interference conventions, and branching-ratio machinery"
                ),
            },
        )
