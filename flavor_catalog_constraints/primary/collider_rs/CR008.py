"""CR008 - non-custodial singlet vector-like ``T`` pair-production limit.

Physics
-------
The observable is the singlet charge-2/3 vector-like ``T`` mass compared with
the catalogued ATLAS pair-production exclusion for

    pp -> T Tbar,  B(T -> Wb:Ht:Zt) = 1/2:1/4:1/4.

The reusable mass-vs-limit recast is the ``collider_resonance`` machinery in
``quarkConstraints.collider_resonance``, reached only through
``flavor_catalog_constraints.physics_adapters.collider_resonance``.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete recast needs
``sigma(pp -> T Tbar) * BR(T -> Wb/Ht/Zt)^2``, the actual singlet top-partner
spectrum, width assumptions, acceptance, and the experiment's mass-dependent
limit curve.  The current ``ParameterPoint`` carries the common quark-sector
``M_KK`` scale, not a dedicated non-custodial singlet VLQ spectrum, so CR008
uses the documented proxy ``m_T = M_KK``.

Severity
--------
HARD.  The active budget is the ATLAS 2024 singlet benchmark lower limit
loaded from ``CR008.yaml``.  The CMS 2023 all-third-generation decay envelope
is retained as diagnostics because it assumes a broader branching-fraction
scan than the CR008 singlet channel.  The ratio is ``m_limit / m_T`` so
``ratio <= 1`` passes.

Catalog sidecar
---------------
``CR008.yaml`` stores numerical limits in a ``pdg_or_equivalent.values`` list.
This module adapts selected list entries into the scaffold ``load_anchor`` path
and fails loudly if the expected value IDs or units are missing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
)
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
_SCAFFOLD_UNCERTAINTY_KEY = "__cr008_uncertainty_is_not_used__"
_EXPECTED_UNITS = "TeV"
_RESONANCE = "T pair"
_FINAL_STATE = "Wb/Ht/Zt singlet"

_ACTIVE_ATLAS_SINGLET = "ATLAS2024:CR008:T_singlet_pair_mass_limit"
_CMS_ALL_THIRD_GENERATION = "CMS2023:CR008:T_pair_all_third_generation_decays_envelope"
_KNOWN_VALUE_IDS = (_ACTIVE_ATLAS_SINGLET, _CMS_ALL_THIRD_GENERATION)
_BUDGET_POLICY = (
    "active HARD budget is the ATLAS2024 isospin-singlet T pair-production "
    "mass lower limit in CR008.yaml for B(T -> Wb:Ht:Zt)=1/2:1/4:1/4. "
    "The CMS2023 all-third-generation decay envelope is kept as diagnostics "
    "because it is a broader branching-fraction scan."
)


@dataclass(frozen=True)
class CR008ValueAnchor:
    """Typed CR008 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    experiment: str | None
    sqrt_s: str | None
    source_key: str | None
    model_assumptions: Mapping[str, str]
    pdg_crosscheck_source_key: str | None
    pdg_crosscheck_source_url: str | None

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
class CR008Anchor:
    """YAML-loaded CR008 singlet-T mass-limit bundle."""

    active_limit: CR008ValueAnchor
    cms_all_third_generation: CR008ValueAnchor
    all_limits: tuple[CR008ValueAnchor, ...]
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
            f"{process_id}: CR008 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR008 field {field_name!r}={value!r} is not finite"
        )
    return number


def _string_mapping(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> Mapping[str, str]:
    if not isinstance(value, Mapping):
        raise AnchorError(f"{process_id}: {field_name} must be a mapping")
    return {str(key): str(item) for key, item in value.items()}


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


def _load_mass_limit_anchor(value_id: str, *, process_id: str) -> CR008ValueAnchor:
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
            f"expected {block_key!r} for CR008 value_id {value_id!r}"
        )

    pdg_crosscheck = entry.get("pdg_crosscheck")
    if pdg_crosscheck is not None and not isinstance(pdg_crosscheck, Mapping):
        raise AnchorError(f"{process_id}: {value_id}.pdg_crosscheck must be a mapping")

    return CR008ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        sqrt_s=_optional_str(entry.get("sqrt_s")),
        source_key=_optional_str(entry.get("source_key")),
        model_assumptions=_string_mapping(
            entry.get("model_assumptions"),
            process_id=process_id,
            field_name=f"{value_id}.model_assumptions",
        ),
        pdg_crosscheck_source_key=(
            _optional_str(pdg_crosscheck.get("source_key"))
            if isinstance(pdg_crosscheck, Mapping)
            else None
        ),
        pdg_crosscheck_source_url=(
            _optional_str(pdg_crosscheck.get("source_url"))
            if isinstance(pdg_crosscheck, Mapping)
            else None
        ),
    )


def _load_cr008_anchor(process_id: str) -> CR008Anchor:
    limits = tuple(
        _load_mass_limit_anchor(value_id, process_id=process_id)
        for value_id in _KNOWN_VALUE_IDS
    )
    by_id = {limit.value_id: limit for limit in limits}
    anchor = CR008Anchor(
        active_limit=by_id[_ACTIVE_ATLAS_SINGLET],
        cms_all_third_generation=by_id[_CMS_ALL_THIRD_GENERATION],
        all_limits=limits,
        budget_policy=_BUDGET_POLICY,
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: CR008 active budget must be positive")
    return anchor


def _limit_from_anchor(anchor: CR008ValueAnchor, *, process_id: str) -> ColliderResonanceLimit:
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
        benchmark_model=anchor.model_assumptions.get("branching_fractions"),
        diagnostics={
            "value_id": anchor.value_id,
            "entry_index": anchor.entry_index,
            "display": anchor.display,
            "experiment": anchor.experiment,
            "sqrt_s": anchor.sqrt_s,
            "source_key": anchor.source_key,
            "model_assumptions": dict(anchor.model_assumptions),
            "pdg_crosscheck_source_key": anchor.pdg_crosscheck_source_key,
            "pdg_crosscheck_source_url": anchor.pdg_crosscheck_source_url,
            "snapshot_path": anchor.snapshot_path,
            "yaml_units": anchor.units,
        },
    )


@register
class Constraint:
    """Catalogued singlet vector-like ``T`` pair-production mass constraint."""

    process_id = "CR008"
    severity = Severity.HARD
    observable = "m(T singlet pair -> Wb/Ht/Zt)"

    def __init__(self) -> None:
        self.anchor = _load_cr008_anchor(self.process_id)
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
                notes=f"invalid singlet-T proxy mass input; CR008 not evaluated: {exc}",
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "required_mass_source": f"{_REQUIRED_EXTRA}.M_KK",
                    "needs_human_physics": VLQ_T_PAIR_MASS_PROXY_ASSUMPTION_V1,
                    "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
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
                    f"extra {_REQUIRED_EXTRA!r} absent; singlet T "
                    "pair-production constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "required_mass_source": f"{_REQUIRED_EXTRA}.M_KK",
                    "budget_policy": self.anchor.budget_policy,
                    "active_value_id": self.anchor.active_limit.value_id,
                    "needs_human_physics": VLQ_T_PAIR_MASS_PROXY_ASSUMPTION_V1,
                    "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
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
                "m_t_singlet_proxy_gev": float(m_kk_gev),
                "m_t_singlet_proxy_tev": float(comparison.predicted_mass_tev),
                "mass_proxy": "m_T = M_KK",
                "mass_source": mass_source,
                "active_value_id": self.anchor.active_limit.value_id,
                "active_experiment": self.anchor.active_limit.experiment,
                "active_display": self.anchor.active_limit.display,
                "active_limit_gev": float(self.anchor.active_limit.value_gev),
                "active_branching_fractions": (
                    self.anchor.active_limit.model_assumptions.get(
                        "branching_fractions"
                    )
                ),
                "active_topology": self.anchor.active_limit.model_assumptions.get(
                    "topology"
                ),
                "budget_policy": self.anchor.budget_policy,
                "all_mass_limits_tev": {
                    limit.value_id: float(limit.value_tev)
                    for limit in self.anchor.all_limits
                },
                "cms_all_third_generation_limit_tev": float(
                    self.anchor.cms_all_third_generation.value_tev
                ),
                "needs_human_physics": VLQ_T_PAIR_MASS_PROXY_ASSUMPTION_V1,
                "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
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
                "Singlet vector-like T pair production uses the documented "
                "mass proxy m_T = M_KK and the active ATLAS singlet benchmark "
                "mass lower bound. HARD ratio is m_limit/m_T; full sigma*BR "
                "and detector recast is marked NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
