"""CR002 - vector-like ``T_5/3`` pair-production mass limit.

Physics
-------
The observable is the vector-like charge-5/3 top-partner mass compared with
the catalogued ATLAS/CMS pair-production exclusions in
``flavor_catalog/processes/collider_rs/CR002.yaml``.  The reusable
mass-limit comparison machinery is the CR001 path in
``quarkConstraints.collider_resonance``, reached only through the
``flavor_catalog_constraints.physics_adapters.collider_resonance`` adapter.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete collider recast needs
``sigma(pp -> T_5/3 Tbar_5/3) * BR(T_5/3 -> t W)^2``, the branching-ratio
pattern, width assumptions, acceptance, and the experiment's mass-dependent
limit curve.  The current ``ParameterPoint`` carries the common quark-sector
``M_KK`` scale, not a dedicated custodial VLQ spectrum, so CR002 uses the
documented proxy ``m_VLQ = M_KK`` and compares that mass to the active
benchmark mass lower bound.

Severity
--------
HARD.  The active budget is the PDG/ATLAS 2023 exclusive ``B/X -> Wt``
pair-production lower bound loaded from ``CR002.yaml``.  The ratio is
``m_limit / m_VLQ`` so ``ratio <= 1`` passes.

Catalog sidecar
---------------
``CR002.yaml`` stores its numerical limits in a ``pdg_or_equivalent.values``
list.  This module adapts selected list entries into the scaffold
``load_anchor`` path and fails loudly if the expected value IDs or units are
missing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.collider_resonance import (
    MASS_LOWER_BOUND,
    VLQ_PAIR_MASS_PROXY_ASSUMPTION_V1,
    ColliderResonanceLimit,
    evaluate_collider_resonance_limit,
    resolve_vlq_mkk_gev,
    vlq_pair_prediction_from_m_kk_gev,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "collider_rs"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_PARENT_KEY = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__cr002_uncertainty_is_not_used__"
_EXPECTED_UNITS = "GeV"
_LIMIT_UNITS = "TeV"
_GEV_PER_TEV = 1000.0
_RESONANCE = "T_5/3 pair"
_FINAL_STATE = "tW tW"

_ACTIVE_ATLAS_PDG_2026 = "PDG2026:CR002:ATLAS2023:X53_pair_Wt"
_ATLAS_2023_DEGENERATE = "ATLAS2023:CR002:mass_degenerate_doublet"
_CMS_2019_RH = "CMS2019:CR002:X53_RH"
_CMS_2019_LH = "CMS2019:CR002:X53_LH"
_ATLAS_2018_PAIR_ONLY = "ATLAS2018:CR002:T53_pair_only"
_ATLAS_2018_PAIR_PLUS_SINGLE = "ATLAS2018:CR002:T53_pair_plus_single_unit_coupling"
_CMS_2017_RH = "CMS2017:CR002:X53_RH"
_CMS_2017_LH = "CMS2017:CR002:X53_LH"
_CMS_2014_SAME_SIGN = "CMS2014:CR002:X53_same_sign"
_KNOWN_VALUE_IDS = (
    _ACTIVE_ATLAS_PDG_2026,
    _ATLAS_2023_DEGENERATE,
    _CMS_2019_RH,
    _CMS_2019_LH,
    _ATLAS_2018_PAIR_ONLY,
    _ATLAS_2018_PAIR_PLUS_SINGLE,
    _CMS_2017_RH,
    _CMS_2017_LH,
    _CMS_2014_SAME_SIGN,
)
_BUDGET_POLICY = (
    "active HARD budget is the PDG2026/ATLAS2023 exclusive B/X -> Wt "
    "pair-production mass lower limit in CR002.yaml. The stronger ATLAS2018 "
    "pair-plus-single interpretation is kept as diagnostics only because it "
    "assumes a unit t' t W single-production coupling."
)


@dataclass(frozen=True)
class CR002ValueAnchor:
    """Typed CR002 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    confidence_level: str | None
    cl: str | None
    experiment: str | None
    source_key: str | None
    supporting_source_key: str | None
    supporting_source_url: str | None
    model_assumptions: tuple[str, ...]

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def value_gev(self) -> float:
        return self.anchor.value

    @property
    def value_tev(self) -> float:
        return float(self.anchor.value / _GEV_PER_TEV)

    @property
    def budget(self) -> float:
        return self.value_tev

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
class CR002Anchor:
    """YAML-loaded CR002 VLQ mass-limit bundle."""

    active_limit: CR002ValueAnchor
    atlas_2023_degenerate: CR002ValueAnchor
    cms_2019_rh: CR002ValueAnchor
    cms_2019_lh: CR002ValueAnchor
    atlas_2018_pair_only: CR002ValueAnchor
    atlas_2018_pair_plus_single: CR002ValueAnchor
    all_limits: tuple[CR002ValueAnchor, ...]
    parent_source: str | None
    parent_source_key: str | None
    parent_source_url: str | None
    budget_policy: str

    @property
    def value(self) -> float:
        return self.active_limit.value_tev

    @property
    def value_gev(self) -> float:
        return self.active_limit.value_gev

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
            f"{process_id}: CR002 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR002 field {field_name!r}={value!r} is not finite"
        )
    return number


def _string_tuple(value: Any, *, process_id: str, field_name: str) -> tuple[str, ...]:
    if value is None:
        return ()
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise AnchorError(f"{process_id}: {field_name} must be a list of strings")
    return tuple(str(item) for item in value)


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


def _load_mass_limit_anchor(value_id: str, *, process_id: str) -> CR002ValueAnchor:
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
            f"expected {block_key!r} for CR002 value_id {value_id!r}"
        )

    return CR002ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        confidence_level=_optional_str(entry.get("confidence_level")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        source_key=_optional_str(entry.get("source_key")),
        supporting_source_key=_optional_str(entry.get("supporting_source_key")),
        supporting_source_url=_optional_str(entry.get("supporting_source_url")),
        model_assumptions=_string_tuple(
            entry.get("model_assumptions"),
            process_id=process_id,
            field_name=f"{value_id}.model_assumptions",
        ),
    )


def _load_cr002_anchor(process_id: str) -> CR002Anchor:
    limits = tuple(
        _load_mass_limit_anchor(value_id, process_id=process_id)
        for value_id in _KNOWN_VALUE_IDS
    )
    by_id = {limit.value_id: limit for limit in limits}
    anchor = CR002Anchor(
        active_limit=by_id[_ACTIVE_ATLAS_PDG_2026],
        atlas_2023_degenerate=by_id[_ATLAS_2023_DEGENERATE],
        cms_2019_rh=by_id[_CMS_2019_RH],
        cms_2019_lh=by_id[_CMS_2019_LH],
        atlas_2018_pair_only=by_id[_ATLAS_2018_PAIR_ONLY],
        atlas_2018_pair_plus_single=by_id[_ATLAS_2018_PAIR_PLUS_SINGLE],
        all_limits=limits,
        parent_source=_optional_str(_parent(process_id).get("source")),
        parent_source_key=_optional_str(_parent(process_id).get("source_key")),
        parent_source_url=_optional_str(_parent(process_id).get("source_url")),
        budget_policy=_BUDGET_POLICY,
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: CR002 active budget must be positive")
    return anchor


def _limit_from_anchor(anchor: CR002ValueAnchor, *, process_id: str) -> ColliderResonanceLimit:
    return ColliderResonanceLimit(
        process_id=process_id,
        resonance=_RESONANCE,
        final_state=_FINAL_STATE,
        limit_kind=MASS_LOWER_BOUND,
        value=float(anchor.value_tev),
        units=_LIMIT_UNITS,
        cl=anchor.cl,
        source=anchor.source,
        source_url=anchor.source_url,
        limit_type=anchor.limit_type,
        benchmark_model="; ".join(anchor.model_assumptions) or None,
        diagnostics={
            "value_id": anchor.value_id,
            "entry_index": anchor.entry_index,
            "display": anchor.display,
            "experiment": anchor.experiment,
            "source_key": anchor.source_key,
            "supporting_source_key": anchor.supporting_source_key,
            "supporting_source_url": anchor.supporting_source_url,
            "snapshot_path": anchor.snapshot_path,
            "yaml_value_gev": float(anchor.value_gev),
            "yaml_units": anchor.units,
        },
    )


@register
class Constraint:
    """Catalogued vector-like ``T_5/3`` pair-production mass constraint."""

    process_id = "CR002"
    severity = Severity.HARD
    observable = "m(T_5/3 pair -> tW tW)"

    def __init__(self) -> None:
        self.anchor = _load_cr002_anchor(self.process_id)
        self.limit = _limit_from_anchor(
            self.anchor.active_limit,
            process_id=self.process_id,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        try:
            m_kk_gev, mass_source = resolve_vlq_mkk_gev(
                couplings=point.get_extra(_REQUIRED_EXTRA),
            )
        except (TypeError, ValueError, AttributeError) as exc:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=False,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=f"invalid VLQ proxy mass input; CR002 not evaluated: {exc}",
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "required_mass_source": f"{_REQUIRED_EXTRA}.M_KK",
                    "needs_human_physics": VLQ_PAIR_MASS_PROXY_ASSUMPTION_V1,
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
                    f"extra {_REQUIRED_EXTRA!r} absent; vector-like T_5/3 "
                    "pair-production constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "required_mass_source": f"{_REQUIRED_EXTRA}.M_KK",
                    "budget_policy": self.anchor.budget_policy,
                    "active_value_id": self.anchor.active_limit.value_id,
                    "needs_human_physics": VLQ_PAIR_MASS_PROXY_ASSUMPTION_V1,
                },
            )

        prediction = vlq_pair_prediction_from_m_kk_gev(float(m_kk_gev))
        comparison = evaluate_collider_resonance_limit(prediction, self.limit)
        diagnostics = dict(comparison.diagnostics)
        diagnostics.update(
            {
                "m_kk_gev": float(m_kk_gev),
                "m_vlq_proxy_gev": float(m_kk_gev),
                "m_vlq_proxy_tev": float(comparison.predicted_mass_tev),
                "mass_proxy": "m_VLQ = M_KK",
                "mass_source": mass_source,
                "active_value_id": self.anchor.active_limit.value_id,
                "active_experiment": self.anchor.active_limit.experiment,
                "active_display": self.anchor.active_limit.display,
                "active_limit_gev": float(self.anchor.active_limit.value_gev),
                "budget_policy": self.anchor.budget_policy,
                "parent_source": self.anchor.parent_source,
                "parent_source_key": self.anchor.parent_source_key,
                "parent_source_url": self.anchor.parent_source_url,
                "all_mass_limits_gev": {
                    limit.value_id: float(limit.value_gev)
                    for limit in self.anchor.all_limits
                },
                "atlas2018_pair_plus_single_not_used_gev": float(
                    self.anchor.atlas_2018_pair_plus_single.value_gev
                ),
                "needs_human_physics": VLQ_PAIR_MASS_PROXY_ASSUMPTION_V1,
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
                "Vector-like T_5/3 pair production uses the documented mass "
                "proxy m_VLQ = M_KK and the active ATLAS/PDG pair-production "
                "mass lower bound. HARD ratio is m_limit/m_VLQ; full "
                "sigma*BR recast is marked NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
