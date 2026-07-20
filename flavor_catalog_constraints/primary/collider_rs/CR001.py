"""CR001 - KK-gluon resonance in ``t tbar``.

Physics
-------
The observable is the first KK-gluon mass compared with the catalogued
ATLAS/CMS ``pp -> g_KK -> t tbar`` benchmark exclusion.  The reusable
mass-limit and future ``sigma*BR`` comparison machinery lives in
``quarkConstraints.collider_resonance`` and is reached only through the
``flavor_catalog_constraints.physics_adapters.collider_resonance`` boundary.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete collider recast needs
``sigma(pp -> g_KK) * BR(g_KK -> t tbar)``, total width, interference,
acceptance, and the experiment's mass-dependent limit curve.  The current
``ParameterPoint`` carries the KK scale but not that collider model, so CR001
uses the documented benchmark mass-exclusion proxy.

Severity
--------
HARD.  The active budget is the current CMS-B2G-25-009 95% CL KK-gluon
benchmark exclusion edge loaded from ``CR001.yaml``.  The ratio is
``m_limit / m_KK`` so ``ratio <= 1`` passes.

Catalog sidecar
---------------
``flavor_catalog/processes/collider_rs/CR001.yaml`` is the source of truth for
the CMS/ATLAS mass-exclusion values.  CR001 stores its numbers in a
``pdg_or_equivalent.values`` list, so this module adapts selected list entries
into the scaffold ``load_anchor`` path and fails loudly if the expected value
IDs or units are missing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.collider_resonance import (
    COLLIDER_RESONANCE_MASS_PROXY_ASSUMPTION_V1,
    MASS_LOWER_BOUND,
    ColliderResonanceLimit,
    evaluate_collider_resonance_limit,
    kk_gluon_prediction_from_m_kk_gev,
    resolve_kk_gluon_mass_gev,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "collider_rs"
_MASS_EXTRA = "kk_gluon_mass_gev"
_COUPLINGS_EXTRA = "quark_mass_basis_couplings"
_PARENT_KEY = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__cr001_uncertainty_is_not_used__"
_EXPECTED_UNITS = "TeV"
_ACTIVE_CMS_2026 = "CMS2026:CR001:gkk_ttbar_mass_exclusion"
_CMS_2019 = "CMS2019:CR001:gkk_ttbar_mass_exclusion"
_PDG_2026 = "PDGLive2026:CR001:s071kkg_aaboud2018bi"
_ATLAS_2013 = "ATLAS2013:CR001:gkk_ttbar_mass_exclusion"
_ATLAS_2012 = "ATLAS2012:CR001:gkk_ttbar_mass_exclusion"
_BUDGET_POLICY = (
    "active HARD budget is the CMS2026 current canonical RS KK-gluon "
    "ttbar benchmark mass-exclusion upper edge in CR001.yaml, interpreted as "
    "a lower bound above the excluded interval for that benchmark."
)


@dataclass(frozen=True)
class CR001ValueAnchor:
    """Typed CR001 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    experiment: str | None
    source_key: str | None
    mass_interval_low: float | None
    benchmark_model: str | None

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def budget(self) -> float:
        return self.anchor.value

    @property
    def units(self) -> str | None:
        return self.anchor.units

    @property
    def source(self) -> str | None:
        return self.anchor.source

    @property
    def source_url(self) -> str | None:
        return self.anchor.source_url

    @property
    def snapshot_path(self) -> str | None:
        return self.anchor.snapshot_path


@dataclass(frozen=True)
class CR001Anchor:
    """YAML-loaded CR001 mass-limit bundle."""

    active_limit: CR001ValueAnchor
    cms_2019: CR001ValueAnchor
    pdg_2026: CR001ValueAnchor
    atlas_2013: CR001ValueAnchor
    atlas_2012: CR001ValueAnchor
    source_summary: str | None
    budget_policy: str

    @property
    def value(self) -> float:
        return self.active_limit.value

    @property
    def budget(self) -> float:
        return self.active_limit.budget

    @property
    def source_url(self) -> str | None:
        return self.active_limit.source_url


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: CR001 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR001 field {field_name!r}={value!r} is not finite"
        )
    return number


def _optional_positive_float(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> float | None:
    if value is None:
        return None
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: {field_name} must be positive")
    return number


def _parent(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    parent = data.get(_PARENT_KEY)
    if not isinstance(parent, Mapping):
        raise AnchorError(f"{process_id}: expected mapping-shaped {_PARENT_KEY}")
    return parent


def _value_entries(process_id: str) -> Sequence[Mapping[str, Any]]:
    parent = _parent(process_id)
    values = parent.get("values")
    if not isinstance(values, list) or not values:
        raise AnchorError(f"{process_id}: expected non-empty {_PARENT_KEY}.values")
    for index, entry in enumerate(values):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: {_PARENT_KEY}.values[{index}] is not a mapping"
            )
    return values


def _entry_by_value_id(
    process_id: str,
    value_id: str,
) -> tuple[int, Mapping[str, Any]]:
    matches: list[tuple[int, Mapping[str, Any]]] = []
    for index, entry in enumerate(_value_entries(process_id)):
        if entry.get("value_id") == value_id:
            matches.append((index, entry))
    if not matches:
        present = [str(entry.get("value_id")) for entry in _value_entries(process_id)]
        raise AnchorError(
            f"{process_id}: value_id {value_id!r} not found in "
            f"{_PARENT_KEY}.values (present: {present})"
        )
    if len(matches) > 1:
        raise AnchorError(f"{process_id}: duplicate value_id {value_id!r}")
    return matches[0]


def _benchmark_model(entry: Mapping[str, Any]) -> str | None:
    assumptions = entry.get("benchmark_assumptions")
    if not isinstance(assumptions, Mapping):
        return None
    return _optional_str(assumptions.get("model"))


def _load_mass_limit_anchor(value_id: str, *, process_id: str) -> CR001ValueAnchor:
    index, entry = _entry_by_value_id(process_id, value_id)
    value = _required_float(
        entry.get("value"),
        process_id=process_id,
        field_name=f"{value_id}.value",
    )
    if value <= 0.0:
        raise AnchorError(f"{process_id}: {value_id}.value must be positive")
    if entry.get("units") != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_UNITS!r} for "
            f"{value_id}, got {entry.get('units')!r}"
        )

    block_key = f"{_PARENT_KEY}.values[{index}]"
    virtual_entry = dict(entry)
    virtual_entry["value"] = value
    virtual_block = {block_key: virtual_entry}
    original_load_pdg_block = anchor_scaffold.load_pdg_block

    def _load_virtual_pdg_block(
        request_process_id: str,
        **kwargs: Any,
    ) -> Mapping[str, Any]:
        if request_process_id == process_id and kwargs.get("family") == _FAMILY:
            return virtual_block
        return original_load_pdg_block(request_process_id, **kwargs)

    anchor_scaffold.load_pdg_block = _load_virtual_pdg_block
    try:
        scaffold_anchor = load_anchor(
            process_id,
            family=_FAMILY,
            candidates=(block_key,),
            uncertainty_key=_SCAFFOLD_UNCERTAINTY_KEY,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for CR001 value_id {value_id!r}"
        )

    mass_interval_low = _optional_positive_float(
        entry.get("benchmark_assumptions", {}).get("mass_interval_low")
        if isinstance(entry.get("benchmark_assumptions"), Mapping)
        else None,
        process_id=process_id,
        field_name=f"{value_id}.benchmark_assumptions.mass_interval_low",
    )
    return CR001ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        source_key=_optional_str(entry.get("source_key")),
        mass_interval_low=mass_interval_low,
        benchmark_model=_benchmark_model(entry),
    )


def _load_cr001_anchor(process_id: str) -> CR001Anchor:
    active = _load_mass_limit_anchor(_ACTIVE_CMS_2026, process_id=process_id)
    anchor = CR001Anchor(
        active_limit=active,
        cms_2019=_load_mass_limit_anchor(_CMS_2019, process_id=process_id),
        pdg_2026=_load_mass_limit_anchor(_PDG_2026, process_id=process_id),
        atlas_2013=_load_mass_limit_anchor(_ATLAS_2013, process_id=process_id),
        atlas_2012=_load_mass_limit_anchor(_ATLAS_2012, process_id=process_id),
        source_summary=_optional_str(_parent(process_id).get("source_summary")),
        budget_policy=_BUDGET_POLICY,
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: CR001 active budget must be positive")
    return anchor


def _limit_from_anchor(anchor: CR001ValueAnchor, *, process_id: str) -> ColliderResonanceLimit:
    return ColliderResonanceLimit(
        process_id=process_id,
        resonance="g_KK^(1)",
        final_state="t tbar",
        limit_kind=MASS_LOWER_BOUND,
        value=float(anchor.value),
        units=str(anchor.units),
        cl=anchor.cl,
        source=anchor.source,
        source_url=anchor.source_url,
        limit_type=anchor.limit_type,
        benchmark_model=anchor.benchmark_model,
        mass_interval_low=anchor.mass_interval_low,
        diagnostics={
            "value_id": anchor.value_id,
            "entry_index": anchor.entry_index,
            "display": anchor.display,
            "experiment": anchor.experiment,
            "source_key": anchor.source_key,
            "snapshot_path": anchor.snapshot_path,
        },
    )


@register
class Constraint:
    """Catalogued direct KK-gluon ``ttbar`` resonance mass constraint."""

    process_id = "CR001"
    severity = Severity.HARD
    observable = "m(g_KK -> t tbar)"

    def __init__(self) -> None:
        self.anchor = _load_cr001_anchor(self.process_id)
        self.limit = _limit_from_anchor(
            self.anchor.active_limit,
            process_id=self.process_id,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        try:
            m_kk_gev, mass_source = resolve_kk_gluon_mass_gev(
                mass_extra=point.get_extra(_MASS_EXTRA),
                couplings=point.get_extra(_COUPLINGS_EXTRA),
            )
        except (TypeError, ValueError, AttributeError) as exc:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=False,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=f"invalid KK-gluon mass input; CR001 not evaluated: {exc}",
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "required_mass_sources": (_MASS_EXTRA, f"{_COUPLINGS_EXTRA}.M_KK"),
                },
            )

        if m_kk_gev is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extras {_MASS_EXTRA!r} and {_COUPLINGS_EXTRA!r} absent; "
                    "KK-gluon ttbar resonance constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extras": (_MASS_EXTRA, _COUPLINGS_EXTRA),
                    "budget_policy": self.anchor.budget_policy,
                    "active_value_id": self.anchor.active_limit.value_id,
                },
            )

        prediction = kk_gluon_prediction_from_m_kk_gev(float(m_kk_gev))
        comparison = evaluate_collider_resonance_limit(prediction, self.limit)
        diagnostics = dict(comparison.diagnostics)
        diagnostics.update(
            {
                "m_kk_gev": float(m_kk_gev),
                "mass_source": mass_source,
                "active_value_id": self.anchor.active_limit.value_id,
                "active_experiment": self.anchor.active_limit.experiment,
                "active_display": self.anchor.active_limit.display,
                "budget_policy": self.anchor.budget_policy,
                "source_summary": self.anchor.source_summary,
                "historical_mass_limits_tev": {
                    "CMS2019": float(self.anchor.cms_2019.value),
                    "PDGLive2026": float(self.anchor.pdg_2026.value),
                    "ATLAS2013": float(self.anchor.atlas_2013.value),
                    "ATLAS2012": float(self.anchor.atlas_2012.value),
                },
                "needs_human_physics": COLLIDER_RESONANCE_MASS_PROXY_ASSUMPTION_V1,
                "sigma_times_br_proxy_pb": comparison.predicted_sigma_times_br,
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(comparison.passes),
            predicted=float(comparison.predicted_mass_tev),
            experimental=float(comparison.experimental_limit),
            ratio=float(comparison.ratio_to_budget),
            budget=float(comparison.budget),
            notes=(
                "KK-gluon ttbar resonance uses the documented benchmark "
                "mass-exclusion proxy. HARD ratio is m_limit/m_KK; "
                "full sigma*BR recast is marked NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
