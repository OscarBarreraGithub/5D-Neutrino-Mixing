"""CR010 - weak-isospin ``(T,B)`` vector-like-quark doublet limit.

Physics
-------
The observable is the common weak-isospin vector-like-quark doublet mass
compared with the catalogued simultaneous ATLAS ``(T,B)`` pair-production
exclusion,

    m_T, m_B > 1.37 TeV  for a weak-isospin (T,B) doublet.

The reusable mass-vs-limit comparison is the ``collider_resonance`` machinery
in ``quarkConstraints.collider_resonance``, reached only through
``flavor_catalog_constraints.physics_adapters.collider_resonance``.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete recast needs
``sigma(pp -> T Tbar/B Bbar) * BR`` over the doublet branching-fraction
surface, the actual ``T`` and ``B`` spectrum, widths, mass splittings,
relative T/B rates, acceptance, nonstandard RS cascade decays, and the
experiment's mass-dependent limit surface.  The current ``ParameterPoint``
carries the common quark-sector ``M_KK`` scale, not a dedicated weak-isospin
doublet spectrum, so CR010 uses the documented proxy ``m_T = m_B = M_KK``.

Severity
--------
HARD.  The active budget is the simultaneous ATLAS/PDG 2025 ``(T,B)`` doublet
lower limit loaded from ``CR010.yaml``.  Individual one-particle endpoints
from the YAML are retained as diagnostics because their own
``model_assumption`` fields say they are not simultaneous ``(T,B)`` doublet
exclusions.  The ratio is ``m_limit / M_KK`` so ``ratio <= 1`` passes.

Catalog sidecar
---------------
``CR010.yaml`` stores numerical limits in a ``pdg_or_equivalent.values`` list.
This module adapts selected list entries into the scaffold ``load_anchor`` path
and fails loudly if the expected value IDs or units are missing.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
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
    VLQ_TB_DOUBLET_MASS_PROXY_ASSUMPTION_V1,
    ColliderResonanceLimit,
    evaluate_collider_resonance_limit,
    resolve_vlq_tb_doublet_mkk_gev,
    vlq_tb_doublet_prediction_from_m_kk_gev,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "collider_rs"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_PARENT_KEY = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__cr010_uncertainty_is_not_used__"
_EXPECTED_UNITS = "TeV"
_RESONANCE = "(T,B) doublet pair"
_FINAL_STATE = "mixed W/Z/H third-generation final states"

_ACTIVE_T_DOUBLET = "PDG2025:CR010:TB_doublet_T_mass_lower_limit"
_ACTIVE_B_DOUBLET = "PDG2025:CR010:TB_doublet_B_mass_lower_limit"
_T_WB_100PCT = "PDG2025:CR010:T_to_Wb_100pct_mass_lower_limit"
_B_HB_100PCT = "PDG2025:CR010:B_to_Hb_100pct_mass_lower_limit"
_B_WT_100PCT = "PDG2025:CR010:B_to_Wt_100pct_mass_lower_limit"
_B_ZB_100PCT = "PDG2025:CR010:B_to_Zb_100pct_mass_lower_limit"
_CMS_T_UNIFORM = "CMS2023:CR010:T_all_third_generation_decays_uniform_lower_limit"
_KNOWN_VALUE_IDS = (
    _ACTIVE_T_DOUBLET,
    _ACTIVE_B_DOUBLET,
    _T_WB_100PCT,
    _B_HB_100PCT,
    _B_WT_100PCT,
    _B_ZB_100PCT,
    _CMS_T_UNIFORM,
)
_BUDGET_POLICY = (
    "active HARD budget is the simultaneous ATLAS/PDG2025 weak-isospin "
    "(T,B)-doublet mass lower limit in CR010.yaml. The individual T and B "
    "pure-channel endpoints are retained as diagnostics because the YAML "
    "model_assumption fields mark them as non-simultaneous doublet exclusions."
)


@dataclass(frozen=True)
class CR010ValueAnchor:
    """Typed CR010 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    experiment: str | None
    collider: str | None
    source_key: str | None
    original_result_url: str | None
    model_assumption: str | None
    sha256: str | None

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
    def year(self) -> int | None:
        return self.anchor.year

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
class CR010Anchor:
    """YAML-loaded CR010 weak-isospin ``(T,B)`` mass-limit bundle."""

    t_doublet_limit: CR010ValueAnchor
    b_doublet_limit: CR010ValueAnchor
    t_wb_100pct: CR010ValueAnchor
    b_hb_100pct: CR010ValueAnchor
    b_wt_100pct: CR010ValueAnchor
    b_zb_100pct: CR010ValueAnchor
    cms_t_uniform: CR010ValueAnchor
    all_limits: tuple[CR010ValueAnchor, ...]
    budget_policy: str

    @property
    def active_limits(self) -> tuple[CR010ValueAnchor, CR010ValueAnchor]:
        return (self.t_doublet_limit, self.b_doublet_limit)

    @property
    def active_value_ids(self) -> tuple[str, str]:
        return tuple(limit.value_id for limit in self.active_limits)

    @property
    def value(self) -> float:
        """Common degenerate-mass lower bound: both active limits must pass."""

        return max(limit.value_tev for limit in self.active_limits)

    @property
    def value_gev(self) -> float:
        return float(self.value * 1000.0)

    @property
    def budget(self) -> float:
        return self.value

    @property
    def source_url(self) -> str | None:
        return self.t_doublet_limit.source_url


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: CR010 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR010 field {field_name!r}={value!r} is not finite"
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


def _load_mass_limit_anchor(value_id: str, *, process_id: str) -> CR010ValueAnchor:
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
            f"expected {block_key!r} for CR010 value_id {value_id!r}"
        )

    return CR010ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        collider=_optional_str(entry.get("collider")),
        source_key=_optional_str(entry.get("source_key")),
        original_result_url=_optional_str(entry.get("original_result_url")),
        model_assumption=_optional_str(entry.get("model_assumption")),
        sha256=_optional_str(entry.get("sha256")),
    )


def _load_cr010_anchor(process_id: str) -> CR010Anchor:
    limits = tuple(
        _load_mass_limit_anchor(value_id, process_id=process_id)
        for value_id in _KNOWN_VALUE_IDS
    )
    by_id = {limit.value_id: limit for limit in limits}
    anchor = CR010Anchor(
        t_doublet_limit=by_id[_ACTIVE_T_DOUBLET],
        b_doublet_limit=by_id[_ACTIVE_B_DOUBLET],
        t_wb_100pct=by_id[_T_WB_100PCT],
        b_hb_100pct=by_id[_B_HB_100PCT],
        b_wt_100pct=by_id[_B_WT_100PCT],
        b_zb_100pct=by_id[_B_ZB_100PCT],
        cms_t_uniform=by_id[_CMS_T_UNIFORM],
        all_limits=limits,
        budget_policy=_BUDGET_POLICY,
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: CR010 active budget must be positive")
    return anchor


def _combined_text(values: Sequence[str | None]) -> str | None:
    present = tuple(value for value in values if value)
    if not present:
        return None
    if len(set(present)) == 1:
        return present[0]
    return "; ".join(present)


def _limit_from_anchor(anchor: CR010Anchor, *, process_id: str) -> ColliderResonanceLimit:
    return ColliderResonanceLimit(
        process_id=process_id,
        resonance=_RESONANCE,
        final_state=_FINAL_STATE,
        limit_kind=MASS_LOWER_BOUND,
        value=float(anchor.value),
        units=_EXPECTED_UNITS,
        cl=_combined_text([limit.cl for limit in anchor.active_limits]),
        source=_combined_text([limit.source for limit in anchor.active_limits]),
        source_url=_combined_text([limit.source_url for limit in anchor.active_limits]),
        limit_type=_combined_text([limit.limit_type for limit in anchor.active_limits]),
        benchmark_model=_combined_text(
            [limit.model_assumption for limit in anchor.active_limits]
        ),
        diagnostics={
            "active_value_ids": anchor.active_value_ids,
            "active_displays": {
                limit.value_id: limit.display for limit in anchor.active_limits
            },
            "active_experiments": {
                limit.value_id: limit.experiment for limit in anchor.active_limits
            },
            "active_source_keys": {
                limit.value_id: limit.source_key for limit in anchor.active_limits
            },
            "active_model_assumptions": {
                limit.value_id: limit.model_assumption for limit in anchor.active_limits
            },
            "active_limits_tev": {
                limit.value_id: float(limit.value_tev)
                for limit in anchor.active_limits
            },
            "active_entry_indices": {
                limit.value_id: limit.entry_index for limit in anchor.active_limits
            },
            "active_snapshot_paths": {
                limit.value_id: limit.snapshot_path for limit in anchor.active_limits
            },
            "active_yaml_units": {
                limit.value_id: limit.units for limit in anchor.active_limits
            },
        },
    )


@register
class Constraint:
    """Catalogued weak-isospin ``(T,B)`` doublet pair-production constraint."""

    process_id = "CR010"
    severity = Severity.HARD
    observable = "m((T,B) doublet pair -> mixed third-generation final states)"

    def __init__(self) -> None:
        self.anchor = _load_cr010_anchor(self.process_id)
        self.limit = _limit_from_anchor(self.anchor, process_id=self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        try:
            m_kk_gev, mass_source = resolve_vlq_tb_doublet_mkk_gev(
                couplings=point.get_extra(_REQUIRED_EXTRA),
            )
        except (TypeError, ValueError, AttributeError) as exc:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=False,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=f"invalid (T,B)-doublet proxy mass input; CR010 not evaluated: {exc}",
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "required_mass_source": f"{_REQUIRED_EXTRA}.M_KK",
                    "needs_human_physics": VLQ_TB_DOUBLET_MASS_PROXY_ASSUMPTION_V1,
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
                    f"extra {_REQUIRED_EXTRA!r} absent; (T,B)-doublet "
                    "pair-production constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "required_mass_source": f"{_REQUIRED_EXTRA}.M_KK",
                    "budget_policy": self.anchor.budget_policy,
                    "active_value_ids": self.anchor.active_value_ids,
                    "needs_human_physics": VLQ_TB_DOUBLET_MASS_PROXY_ASSUMPTION_V1,
                    "sigma_times_br_recast_status": "NEEDS-HUMAN-PHYSICS",
                },
            )

        prediction = vlq_tb_doublet_prediction_from_m_kk_gev(
            float(m_kk_gev),
            resonance=_RESONANCE,
            final_state=_FINAL_STATE,
        )
        comparison = evaluate_collider_resonance_limit(prediction, self.limit)
        diagnostics = dict(comparison.diagnostics)
        diagnostics.update(
            {
                "m_kk_gev": float(m_kk_gev),
                "m_t_doublet_proxy_gev": float(m_kk_gev),
                "m_b_doublet_proxy_gev": float(m_kk_gev),
                "m_tb_doublet_proxy_tev": float(comparison.predicted_mass_tev),
                "mass_proxy": "m_T = m_B = M_KK",
                "mass_source": mass_source,
                "active_value_ids": self.anchor.active_value_ids,
                "active_limit_gev": float(self.anchor.value_gev),
                "active_t_limit_tev": float(self.anchor.t_doublet_limit.value_tev),
                "active_b_limit_tev": float(self.anchor.b_doublet_limit.value_tev),
                "active_t_model_assumption": self.anchor.t_doublet_limit.model_assumption,
                "active_b_model_assumption": self.anchor.b_doublet_limit.model_assumption,
                "budget_policy": self.anchor.budget_policy,
                "all_mass_limits_tev": {
                    limit.value_id: float(limit.value_tev)
                    for limit in self.anchor.all_limits
                },
                "individual_endpoint_limits_tev": {
                    self.anchor.t_wb_100pct.value_id: float(
                        self.anchor.t_wb_100pct.value_tev
                    ),
                    self.anchor.b_hb_100pct.value_id: float(
                        self.anchor.b_hb_100pct.value_tev
                    ),
                    self.anchor.b_wt_100pct.value_id: float(
                        self.anchor.b_wt_100pct.value_tev
                    ),
                    self.anchor.b_zb_100pct.value_id: float(
                        self.anchor.b_zb_100pct.value_tev
                    ),
                    self.anchor.cms_t_uniform.value_id: float(
                        self.anchor.cms_t_uniform.value_tev
                    ),
                },
                "needs_human_physics": VLQ_TB_DOUBLET_MASS_PROXY_ASSUMPTION_V1,
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
                "Weak-isospin (T,B) doublet pair production uses the documented "
                "mass proxy m_T = m_B = M_KK and the active simultaneous "
                "ATLAS/PDG doublet mass lower bound. HARD ratio is "
                "m_limit/M_KK; full sigma*BR and branching-surface recast is "
                "marked NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
