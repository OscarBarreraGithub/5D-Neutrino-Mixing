"""CR003 - charge-2/3 vector-like top-partner pair-production mass limit.

Physics
-------
The observable is the custodial charge-2/3 top-partner mass compared with the
catalogued ATLAS/CMS pair-production exclusions in
``flavor_catalog/processes/collider_rs/CR003.yaml``.  The reusable
mass-limit comparison machinery is the CR001 path in
``quarkConstraints.collider_resonance``, reached only through the
``flavor_catalog_constraints.physics_adapters.collider_resonance`` adapter.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete collider recast needs
``sigma(pp -> T Tbar) * BR(T -> Wb/Zt/Ht)^2``, the branching-fraction point,
width assumptions, acceptance, and the experiment's mass-dependent limit
curve.  The current ``ParameterPoint`` carries the common quark-sector
``M_KK`` scale, not a dedicated custodial charge-2/3 VLQ spectrum, so CR003
uses the documented proxy ``m_T = M_KK``.

Severity
--------
HARD.  Because ``ParameterPoint`` does not carry VLQ branching fractions, the
active budget is the CMS 2023 all-third-generation decay-mixture lower edge
loaded from ``CR003.yaml``.  Pure ``Wb``, ``Zt``, ``Ht``, and singlet
benchmark limits remain in diagnostics.  The ratio is ``m_limit / m_T`` so
``ratio <= 1`` passes.

Catalog sidecar
---------------
``CR003.yaml`` stores its numerical limits in a ``pdg_or_equivalent.values``
list.  This module adapts selected list entries into the scaffold
``load_anchor`` path and fails loudly if the expected value IDs or units are
missing.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.collider_resonance import (
    MASS_LOWER_BOUND,
    VLQ_T_PAIR_MASS_PROXY_ASSUMPTION_V1,
    ColliderResonanceLimit,
    evaluate_collider_resonance_limit,
    resolve_charge_two_thirds_vlq_mkk_gev,
    vlq_t_pair_prediction_from_m_kk_gev,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "collider_rs"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_PARENT_KEY = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__cr003_uncertainty_is_not_used__"
_EXPECTED_UNITS = "TeV"
_RESONANCE = "T pair"
_FINAL_STATE = "Wb/Zt/Ht"

_PURE_WB_ATLAS_2024 = "PDGLive2026:CR003:T_Wb_pair_ATLAS2024"
_SINGLET_ATLAS_2024 = "PDGLive2026:CR003:T_singlet_pair_ATLAS2024"
_PURE_ZT_ATLAS_2023 = "PDGLive2026:CR003:T_Zt_pair_ATLAS2023"
_PURE_HT_CMS_2023 = "PDGLive2026:CR003:T_Ht_pair_CMS2023"
_ACTIVE_CMS_ALL_MIXTURES = "CMS2023:CR003:T_all_third_generation_decay_mixtures"
_KNOWN_VALUE_IDS = (
    _PURE_WB_ATLAS_2024,
    _SINGLET_ATLAS_2024,
    _PURE_ZT_ATLAS_2023,
    _PURE_HT_CMS_2023,
    _ACTIVE_CMS_ALL_MIXTURES,
)
_BUDGET_POLICY = (
    "active HARD budget is the CMS2023 all-third-generation Wb/Zt/Ht "
    "decay-mixture mass lower edge in CR003.yaml. Stronger pure-channel "
    "limits are retained as diagnostics because ParameterPoint does not "
    "provide VLQ branching fractions."
)


@dataclass(frozen=True)
class CR003ValueAnchor:
    """Typed CR003 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    confidence_level: str | None
    cl: str | None
    experiment: str | None
    source_key: str | None
    production_mode: str | None
    decay_mode: str | None
    branching_assumption: str | None
    normalized_value_gev: float | None

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def value_tev(self) -> float:
        return self.anchor.value

    @property
    def value_gev(self) -> float:
        return float(self.anchor.value * 1000.0)

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
class CR003Anchor:
    """YAML-loaded CR003 charge-2/3 top-partner mass-limit bundle."""

    active_limit: CR003ValueAnchor
    pure_wb_atlas_2024: CR003ValueAnchor
    singlet_atlas_2024: CR003ValueAnchor
    pure_zt_atlas_2023: CR003ValueAnchor
    pure_ht_cms_2023: CR003ValueAnchor
    all_limits: tuple[CR003ValueAnchor, ...]
    parent_canonical_source: str | None
    parent_source_key: str | None
    parent_source_url: str | None
    budget_policy: str

    @property
    def value(self) -> float:
        return self.active_limit.value_tev

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
            f"{process_id}: CR003 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR003 field {field_name!r}={value!r} is not finite"
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


def _load_mass_limit_anchor(value_id: str, *, process_id: str) -> CR003ValueAnchor:
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
            f"expected {block_key!r} for CR003 value_id {value_id!r}"
        )

    normalized_value_gev = _optional_positive_float(
        entry.get("normalized_value_GeV"),
        process_id=process_id,
        field_name=f"{value_id}.normalized_value_GeV",
    )
    if normalized_value_gev is not None and not math.isclose(
        normalized_value_gev,
        value * 1000.0,
        rel_tol=0.0,
        abs_tol=1.0e-9,
    ):
        raise AnchorError(
            f"{process_id}: {value_id}.normalized_value_GeV is inconsistent "
            "with the TeV value"
        )

    return CR003ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        confidence_level=_optional_str(entry.get("confidence_level")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        source_key=_optional_str(entry.get("source_key")),
        production_mode=_optional_str(entry.get("production_mode")),
        decay_mode=_optional_str(entry.get("decay_mode")),
        branching_assumption=_optional_str(entry.get("branching_assumption")),
        normalized_value_gev=normalized_value_gev,
    )


def _load_cr003_anchor(process_id: str) -> CR003Anchor:
    limits = tuple(
        _load_mass_limit_anchor(value_id, process_id=process_id)
        for value_id in _KNOWN_VALUE_IDS
    )
    by_id = {limit.value_id: limit for limit in limits}
    anchor = CR003Anchor(
        active_limit=by_id[_ACTIVE_CMS_ALL_MIXTURES],
        pure_wb_atlas_2024=by_id[_PURE_WB_ATLAS_2024],
        singlet_atlas_2024=by_id[_SINGLET_ATLAS_2024],
        pure_zt_atlas_2023=by_id[_PURE_ZT_ATLAS_2023],
        pure_ht_cms_2023=by_id[_PURE_HT_CMS_2023],
        all_limits=limits,
        parent_canonical_source=_optional_str(_parent(process_id).get("canonical_source")),
        parent_source_key=_optional_str(_parent(process_id).get("source_key")),
        parent_source_url=_optional_str(_parent(process_id).get("source_url")),
        budget_policy=_BUDGET_POLICY,
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: CR003 active budget must be positive")
    return anchor


def _limit_from_anchor(anchor: CR003ValueAnchor, *, process_id: str) -> ColliderResonanceLimit:
    return ColliderResonanceLimit(
        process_id=process_id,
        resonance=_RESONANCE,
        final_state=_FINAL_STATE,
        limit_kind=MASS_LOWER_BOUND,
        value=float(anchor.value_tev),
        units=_EXPECTED_UNITS,
        cl=anchor.cl,
        source=anchor.source,
        source_url=anchor.source_url,
        limit_type=anchor.limit_type,
        benchmark_model=anchor.branching_assumption,
        diagnostics={
            "value_id": anchor.value_id,
            "entry_index": anchor.entry_index,
            "display": anchor.display,
            "experiment": anchor.experiment,
            "source_key": anchor.source_key,
            "production_mode": anchor.production_mode,
            "decay_mode": anchor.decay_mode,
            "branching_assumption": anchor.branching_assumption,
            "snapshot_path": anchor.snapshot_path,
            "normalized_value_gev": anchor.normalized_value_gev,
            "yaml_units": anchor.units,
        },
    )


@register
class Constraint:
    """Catalogued charge-2/3 vector-like top-partner mass constraint."""

    process_id = "CR003"
    severity = Severity.HARD
    observable = "m(T pair -> Wb/Zt/Ht)"

    def __init__(self) -> None:
        self.anchor = _load_cr003_anchor(self.process_id)
        self.limit = _limit_from_anchor(
            self.anchor.active_limit,
            process_id=self.process_id,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        try:
            m_kk_gev, mass_source = resolve_charge_two_thirds_vlq_mkk_gev(
                couplings=point.get_extra(_REQUIRED_EXTRA),
            )
        except (TypeError, ValueError, AttributeError) as exc:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=False,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=f"invalid charge-2/3 VLQ proxy mass input; CR003 not evaluated: {exc}",
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "required_mass_source": f"{_REQUIRED_EXTRA}.M_KK",
                    "needs_human_physics": VLQ_T_PAIR_MASS_PROXY_ASSUMPTION_V1,
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
                    f"extra {_REQUIRED_EXTRA!r} absent; charge-2/3 T "
                    "pair-production constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "required_mass_source": f"{_REQUIRED_EXTRA}.M_KK",
                    "budget_policy": self.anchor.budget_policy,
                    "active_value_id": self.anchor.active_limit.value_id,
                    "needs_human_physics": VLQ_T_PAIR_MASS_PROXY_ASSUMPTION_V1,
                },
            )

        prediction = vlq_t_pair_prediction_from_m_kk_gev(
            float(m_kk_gev),
            resonance=_RESONANCE,
            final_state=_FINAL_STATE,
        )
        comparison = evaluate_collider_resonance_limit(prediction, self.limit)
        diagnostics = dict(comparison.diagnostics)
        diagnostics.update(
            {
                "m_kk_gev": float(m_kk_gev),
                "m_t_partner_proxy_gev": float(m_kk_gev),
                "m_t_partner_proxy_tev": float(comparison.predicted_mass_tev),
                "mass_proxy": "m_T = M_KK",
                "mass_source": mass_source,
                "active_value_id": self.anchor.active_limit.value_id,
                "active_experiment": self.anchor.active_limit.experiment,
                "active_display": self.anchor.active_limit.display,
                "active_limit_gev": float(self.anchor.active_limit.value_gev),
                "active_branching_assumption": (
                    self.anchor.active_limit.branching_assumption
                ),
                "budget_policy": self.anchor.budget_policy,
                "parent_canonical_source": self.anchor.parent_canonical_source,
                "parent_source_key": self.anchor.parent_source_key,
                "parent_source_url": self.anchor.parent_source_url,
                "all_mass_limits_tev": {
                    limit.value_id: float(limit.value_tev)
                    for limit in self.anchor.all_limits
                },
                "pure_wb_limit_tev": float(self.anchor.pure_wb_atlas_2024.value_tev),
                "singlet_limit_tev": float(self.anchor.singlet_atlas_2024.value_tev),
                "pure_zt_limit_tev": float(self.anchor.pure_zt_atlas_2023.value_tev),
                "pure_ht_limit_tev": float(self.anchor.pure_ht_cms_2023.value_tev),
                "needs_human_physics": VLQ_T_PAIR_MASS_PROXY_ASSUMPTION_V1,
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
                "Charge-2/3 vector-like T pair production uses the documented "
                "mass proxy m_T = M_KK and the active CMS all-decay-mixture "
                "mass lower bound. HARD ratio is m_limit/m_T; full sigma*BR "
                "and branching-simplex recast is marked NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
