"""CR006 - charged-current electroweak KK/W' resonance mass limit.

Physics
-------
The observable is the first charged electroweak KK gauge-boson mass in

    pp -> W_KK^(1) -> ell nu, tb

compared with the catalogued high-mass W' benchmark exclusions in
``flavor_catalog/processes/collider_rs/CR006.yaml``.  The reusable
collider-resonance mass-limit comparison machinery lives in
``quarkConstraints.collider_resonance`` and is reached only through the
``flavor_catalog_constraints.physics_adapters.collider_resonance`` boundary.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete recast needs
``sigma(pp -> W_KK) * BR(W_KK -> ell nu, tb)``, light-quark, lepton, and
third-generation quark couplings, total width, interference, acceptance, and
the experiment's mass-dependent limit curve.  The current ``ParameterPoint``
may carry the EW KK mass but not that collider signal model.  CR006 therefore
keeps the documented benchmark mass-exclusion proxy as an advisory diagnostic
only.  The ``sigma*BR`` recast status is reported only as diagnostics.

Severity
--------
INFO.  The active budget is the PDG 2025 / ATLAS 2019 95% CL SSM-W'
``e nu`` benchmark lower mass limit loaded from ``CR006.yaml``.  The advisory
ratio is ``m_limit / m_KK_EW`` so ``ratio <= 1`` passes the SSM proxy, but
failures must not veto RS points.

Catalog sidecar
---------------
``CR006.yaml`` stores its numerical limits in a ``pdg_or_equivalent.values``
list.  This module adapts selected list entries into the scaffold
``load_anchor`` path and fails loudly if the expected value IDs or units are
missing.
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
    KK_CHARGED_CURRENT_MASS_PROXY_ASSUMPTION_V1,
    MASS_LOWER_BOUND,
    ColliderResonanceLimit,
    evaluate_collider_resonance_limit,
    kk_charged_current_prediction_from_m_kk_gev,
    resolve_kk_ew_mass_gev,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "collider_rs"
_MASS_EXTRA = "kk_ew_mass_gev"
_COUPLINGS_EXTRA = "quark_mass_basis_couplings"
_PARENT_KEY = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__cr006_uncertainty_is_not_used__"
_MASS_UNITS = "TeV"
_RESONANCE = "W_KK^(1)"
_FINAL_STATE = "ell nu"

_ACTIVE_PDG_ATLAS_ENU = "PDG2025:CR006:Wprime_SSM_enu_mass_lower_bound"
_PDG_CMS_MUNU = "PDG2025:CR006:Wprime_SSM_munu_mass_lower_bound"
_CMS_COMBINED_LNU = "CMS2022:CR006:Wprime_SSM_enu_munu_combined_mass_lower_bound"
_CMS_TB_R = "CMS2024:CR006:Wprime_R_tb_mass_lower_bound"
_CMS_TB_L = "CMS2024:CR006:Wprime_L_tb_mass_lower_bound"
_KNOWN_VALUE_IDS = (
    _ACTIVE_PDG_ATLAS_ENU,
    _PDG_CMS_MUNU,
    _CMS_COMBINED_LNU,
    _CMS_TB_R,
    _CMS_TB_L,
)
_BUDGET_POLICY = (
    "active INFO budget is the PDG2025/ATLAS2019 SSM W' -> e nu mass "
    "lower limit in CR006.yaml. Per M-26 this raw SSM mass edge is "
    "advisory only for RS electroweak KK bosons because light-quark "
    "production couplings are sqrt(L)-suppressed. CMS combined e/mu and "
    "CMS tb limits are retained as diagnostics because ParameterPoint does "
    "not provide the charged-vector production, branching-fraction, width, "
    "and acceptance model needed to choose a channel-specific RS limit."
)
_M26_SEVERITY_POLICY = "INFO/non-vetoing until an RS sigma*BR recast is supplied"
_RS_VOLUME_LOG_FOR_SUPPRESSION_NOTE = 35.0
_LIGHT_QUARK_COUPLING_RATIO_TO_SSM = 1.0 / math.sqrt(_RS_VOLUME_LOG_FOR_SUPPRESSION_NOTE)
_PRODUCTION_RATE_SUPPRESSION_ESTIMATE = _LIGHT_QUARK_COUPLING_RATIO_TO_SSM**2


@dataclass(frozen=True)
class CR006ValueAnchor:
    """Typed CR006 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    experiment: str | None
    result_year: int | None
    source_key: str | None
    primary_source_key: str | None
    primary_source_url: str | None
    model_assumption: str | None

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def value_tev(self) -> float:
        return self.anchor.value

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
class CR006Anchor:
    """YAML-loaded CR006 charged-current mass-limit bundle."""

    active_limit: CR006ValueAnchor
    pdg_cms_munu: CR006ValueAnchor
    cms_combined_lnu: CR006ValueAnchor
    cms_tb_r: CR006ValueAnchor
    cms_tb_l: CR006ValueAnchor
    all_limits: tuple[CR006ValueAnchor, ...]
    summary: str | None
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


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: CR006 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR006 field {field_name!r}={value!r} is not finite"
        )
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


def _load_value_anchor(value_id: str, *, process_id: str) -> CR006ValueAnchor:
    index, entry = _entry_by_value_id(process_id, value_id)
    value = _required_float(
        entry.get("value"),
        process_id=process_id,
        field_name=f"{value_id}.value",
    )
    if value <= 0.0:
        raise AnchorError(f"{process_id}: {value_id}.value must be positive")
    if entry.get("units") != _MASS_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_MASS_UNITS!r} for "
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
            f"expected {block_key!r} for CR006 value_id {value_id!r}"
        )

    return CR006ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        result_year=_optional_int(entry.get("result_year")),
        source_key=_optional_str(entry.get("source_key")),
        primary_source_key=_optional_str(entry.get("primary_source_key")),
        primary_source_url=_optional_str(entry.get("primary_source_url")),
        model_assumption=_optional_str(entry.get("model_assumption")),
    )


def _load_cr006_anchor(process_id: str) -> CR006Anchor:
    limits = tuple(
        _load_value_anchor(value_id, process_id=process_id)
        for value_id in _KNOWN_VALUE_IDS
    )
    by_id = {limit.value_id: limit for limit in limits}
    anchor = CR006Anchor(
        active_limit=by_id[_ACTIVE_PDG_ATLAS_ENU],
        pdg_cms_munu=by_id[_PDG_CMS_MUNU],
        cms_combined_lnu=by_id[_CMS_COMBINED_LNU],
        cms_tb_r=by_id[_CMS_TB_R],
        cms_tb_l=by_id[_CMS_TB_L],
        all_limits=limits,
        summary=_optional_str(_parent(process_id).get("summary")),
        budget_policy=_BUDGET_POLICY,
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: CR006 active budget must be positive")
    return anchor


def _limit_from_anchor(
    anchor: CR006ValueAnchor,
    *,
    process_id: str,
) -> ColliderResonanceLimit:
    return ColliderResonanceLimit(
        process_id=process_id,
        resonance=_RESONANCE,
        final_state=_FINAL_STATE,
        limit_kind=MASS_LOWER_BOUND,
        value=float(anchor.value_tev),
        units=_MASS_UNITS,
        cl=anchor.cl,
        source=anchor.source,
        source_url=anchor.source_url,
        limit_type=anchor.limit_type,
        benchmark_model=anchor.model_assumption or anchor.anchor.observable,
        diagnostics={
            "value_id": anchor.value_id,
            "entry_index": anchor.entry_index,
            "display": anchor.display,
            "experiment": anchor.experiment,
            "result_year": anchor.result_year,
            "source_key": anchor.source_key,
            "primary_source_key": anchor.primary_source_key,
            "primary_source_url": anchor.primary_source_url,
            "snapshot_path": anchor.snapshot_path,
            "yaml_units": anchor.units,
            "yaml_observable": anchor.anchor.observable,
        },
    )


@register
class Constraint:
    """Catalogued charged-current electroweak KK/W' advisory record."""

    process_id = "CR006"
    severity = Severity.INFO
    observable = "m(W_KK^(1) -> ell nu, tb)"

    def __init__(self) -> None:
        self.anchor = _load_cr006_anchor(self.process_id)
        self.limit = _limit_from_anchor(
            self.anchor.active_limit,
            process_id=self.process_id,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        try:
            m_kk_gev, mass_source = resolve_kk_ew_mass_gev(
                mass_extra=point.get_extra(_MASS_EXTRA),
                couplings=point.get_extra(_COUPLINGS_EXTRA),
            )
        except (TypeError, ValueError, AttributeError) as exc:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "invalid charged-EW KK mass input; non-vetoing CR006 INFO "
                    f"record kept without RS recast: {exc}"
                ),
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "required_mass_sources": (
                        _MASS_EXTRA,
                        f"{_COUPLINGS_EXTRA}.M_KK",
                    ),
                    "needs_human_physics": (
                        KK_CHARGED_CURRENT_MASS_PROXY_ASSUMPTION_V1
                    ),
                    "severity_policy": _M26_SEVERITY_POLICY,
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
                    "charged-current EW KK resonance constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extras": (_MASS_EXTRA, _COUPLINGS_EXTRA),
                    "required_mass_sources": (
                        _MASS_EXTRA,
                        f"{_COUPLINGS_EXTRA}.M_KK",
                    ),
                    "budget_policy": self.anchor.budget_policy,
                    "active_value_id": self.anchor.active_limit.value_id,
                    "needs_human_physics": (
                        KK_CHARGED_CURRENT_MASS_PROXY_ASSUMPTION_V1
                    ),
                    "severity_policy": _M26_SEVERITY_POLICY,
                },
            )

        prediction = kk_charged_current_prediction_from_m_kk_gev(
            float(m_kk_gev),
            resonance=_RESONANCE,
            final_state=_FINAL_STATE,
        )
        comparison = evaluate_collider_resonance_limit(prediction, self.limit)
        diagnostics = dict(comparison.diagnostics)
        diagnostics.update(
            {
                "m_kk_gev": float(m_kk_gev),
                "m_wkk_proxy_gev": float(m_kk_gev),
                "m_wkk_proxy_tev": float(comparison.predicted_mass_tev),
                "mass_proxy": "m_WKK = kk_ew_mass_gev or M_KK",
                "mass_source": mass_source,
                "active_value_id": self.anchor.active_limit.value_id,
                "active_experiment": self.anchor.active_limit.experiment,
                "active_display": self.anchor.active_limit.display,
                "active_model_assumption": (
                    self.anchor.active_limit.model_assumption
                ),
                "budget_policy": self.anchor.budget_policy,
                "all_mass_limits_tev": {
                    limit.value_id: float(limit.value_tev)
                    for limit in self.anchor.all_limits
                },
                "leptonic_mass_limits_tev": {
                    self.anchor.active_limit.value_id: float(
                        self.anchor.active_limit.value_tev
                    ),
                    self.anchor.pdg_cms_munu.value_id: float(
                        self.anchor.pdg_cms_munu.value_tev
                    ),
                    self.anchor.cms_combined_lnu.value_id: float(
                        self.anchor.cms_combined_lnu.value_tev
                    ),
                },
                "tb_mass_limits_tev": {
                    self.anchor.cms_tb_r.value_id: float(
                        self.anchor.cms_tb_r.value_tev
                    ),
                    self.anchor.cms_tb_l.value_id: float(
                        self.anchor.cms_tb_l.value_tev
                    ),
                },
                "needs_human_physics": (
                    KK_CHARGED_CURRENT_MASS_PROXY_ASSUMPTION_V1
                ),
                "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
                "sigma_times_br_proxy_fb": comparison.predicted_sigma_times_br,
                "severity_policy": _M26_SEVERITY_POLICY,
                "m26_ssm_benchmark_status": "advisory_not_hard_veto",
                "rs_light_quark_coupling_ratio_to_ssm": (
                    _LIGHT_QUARK_COUPLING_RATIO_TO_SSM
                ),
                "rs_light_quark_production_rate_suppression_estimate": (
                    _PRODUCTION_RATE_SUPPRESSION_ESTIMATE
                ),
                "rs_sigma_times_br_suppression_note": (
                    "M-26: production is already suppressed by roughly 1/L from "
                    "g_q^RS/g_q^SSM ~ 1/sqrt(L); branching and width effects can "
                    "reduce sigma*BR further, so the raw SSM mass edge is not an "
                    "RS recast."
                ),
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
                "Charged-current electroweak KK/W' resonance uses the documented "
                "benchmark mass-exclusion proxy and the active PDG/ATLAS SSM "
                "e nu mass lower bound as an INFO advisory. Per M-26 this is "
                "not a HARD veto on sqrt(L)-suppressed RS KK electroweak bosons; "
                "full sigma*BR recast is marked NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
