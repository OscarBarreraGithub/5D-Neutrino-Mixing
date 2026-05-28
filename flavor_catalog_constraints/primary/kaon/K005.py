"""K005 - Neutral rare kaon decay ``K_L -> pi0 nu nubar``.

Severity rationale
------------------
INFO placeholder. This observable is a Delta S = 1 rare kaon decay, not
a Delta F = 2 kaon-mixing constraint like K001/K002 and not the
nonleptonic eps'/eps observable in K003. A faithful prediction requires
``s -> d nu nubar`` Wilson matching, CKM inputs, SM--NP interference
conventions, the CP-odd K_L projection, and the neutral-mode
branching-ratio formula. The current physics core has no such
machinery, so this constraint records the experimental upper-limit
anchor and returns a non-vetoing deferred result.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K005.yaml`` (schema-flex anchor:
one-entry list under ``pdg_or_equivalent`` with value ``"<2.2e-9"``).

Physics core
------------
No Delta S = 1 rare-kaon implementation exists yet. The placeholder
reaches :mod:`flavor_catalog_constraints.physics_adapters.rare_kaon_decays`,
whose result carries ``deferred=True``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

from flavor_catalog_constraints.anchors import load_raw
from flavor_catalog_constraints.base import (
    ConstraintResult,
    ParameterPoint,
    Severity,
)
from flavor_catalog_constraints.physics_adapters.rare_kaon_decays import (
    evaluate_klong_pi0_nunu_placeholder,
)
from flavor_catalog_constraints.registry import register_constraint


@dataclass(frozen=True)
class _Anchor:
    """Typed view over ``K005.yaml``'s one-entry ``pdg_or_equivalent`` list."""

    observable: str
    year: int
    raw_value: str
    upper_limit: float
    uncertainty: Any
    units: str
    confidence_level: str
    source: str
    source_url: str
    access_date: str
    snapshot_path: str
    sha256_of_text_snapshot: str
    notes: str


def _parse_upper_limit(raw_value: str) -> float:
    """Parse catalog upper-limit strings like ``<2.2e-9``."""
    text = str(raw_value).strip()
    if not text.startswith("<"):
        raise ValueError(f"Expected upper-limit string starting with '<': {text!r}")
    return float(text[1:].strip())


def _build_anchor(raw: Sequence[Any]) -> _Anchor:
    """K005 uses a one-entry list under ``pdg_or_equivalent``."""
    if not isinstance(raw, list) or len(raw) != 1:
        raise ValueError("K005 expects pdg_or_equivalent to be a one-entry list")

    entry = raw[0]
    raw_value = str(entry["value"])
    return _Anchor(
        observable=str(entry["observable"]),
        year=int(entry["year"]),
        raw_value=raw_value,
        upper_limit=_parse_upper_limit(raw_value),
        uncertainty=entry["uncertainty"],
        units=str(entry["units"]),
        confidence_level=str(entry["confidence_level"]),
        source=str(entry["source"]),
        source_url=str(entry["source_url"]),
        access_date=str(entry["access_date"]),
        snapshot_path=str(entry["snapshot_path"]),
        sha256_of_text_snapshot=str(entry["sha256_of_text_snapshot"]),
        notes=str(entry["notes"]),
    )


@register_constraint
class Constraint:
    """Catalogued K_L -> pi0 nu nubar constraint (process_id ``K005``)."""

    process_id = "K005"
    severity = Severity.INFO
    observable = "BR(K_L -> pi0 nu nubar)"

    def __init__(self) -> None:
        raw = load_raw(self.process_id, family="kaon", tier="primary")
        self.anchor = _build_anchor(raw)
        self.references = (self.anchor.snapshot_path,)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = evaluate_klong_pi0_nunu_placeholder()
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
                "Evaluation deferred: K005 is a Delta S=1 rare kaon "
                "decay with a CP-violating K_L projection and must not "
                "reuse the Delta F=2 K001/K002 mixing machinery, the "
                "K003 eps'/eps placeholder, or the K004 charged-mode "
                "branching-ratio formula without neutral-mode support."
            ),
            diagnostics={
                "deferred": bool(result.deferred),
                "evaluation_deferred": result.reason,
                "confidence_level": self.anchor.confidence_level,
                "upper_limit": self.anchor.upper_limit,
                "required_machinery": (
                    "K -> pi nu nubar Wilson matching, BMU-operator "
                    "mapping, CKM elements, SM-NP interference "
                    "conventions, CP-odd K_L projection, and "
                    "neutral-mode branching-ratio evaluation"
                ),
            },
        )
