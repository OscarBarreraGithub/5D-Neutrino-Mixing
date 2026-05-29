"""CR004 - charge -1/3 vector-like bottom-partner pair-production limit.

Physics
-------
The observable is the custodial charge -1/3 bottom-partner mass compared with
the PDG 2025 ``b'(-1/3)`` pair-production limits in
``flavor_catalog/processes/collider_rs/CR004.yaml``.  The reusable
mass-limit comparison machinery is the CR001 path in
``quarkConstraints.collider_resonance``, reached only through the
``flavor_catalog_constraints.physics_adapters.collider_resonance`` adapter.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete collider recast needs
``sigma(pp -> B Bbar) * BR(B -> tW/bZ/bH)^2``, the branching-fraction point,
width assumptions, acceptance, and the experiment's mass-dependent limit
curve.  The current ``ParameterPoint`` carries the common quark-sector
``M_KK`` scale, not a dedicated custodial charge -1/3 ``B`` spectrum, so CR004
uses the documented proxy ``m_B = M_KK``.

Severity
--------
HARD.  Because ``ParameterPoint`` does not carry bottom-partner branching
fractions, the active budget is the strongest PDG 2025 pure-channel pair
production lower bound loaded from ``CR004.yaml``: ``B(B -> H b)=1``.  The
``B -> Z b`` and ``B -> W t`` limits remain in diagnostics.  The ratio is
``m_limit / m_B`` so ``ratio <= 1`` passes.

Catalog sidecar
---------------
``CR004.yaml`` stores its numerical limits in a ``pdg_or_equivalent.values``
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
    VLQ_B_PAIR_MASS_PROXY_ASSUMPTION_V1,
    ColliderResonanceLimit,
    evaluate_collider_resonance_limit,
    resolve_charge_minus_one_third_vlq_mkk_gev,
    vlq_b_pair_prediction_from_m_kk_gev,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "collider_rs"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_PARENT_KEY = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__cr004_uncertainty_is_not_used__"
_EXPECTED_UNITS = "GeV"
_LIMIT_UNITS = "TeV"
_GEV_PER_TEV = 1000.0
_RESONANCE = "B pair"
_FINAL_STATE = "tW/bZ/bH"

_ACTIVE_BH_PDG_2025 = "PDG2025:CR004:B_bH_pair_mass_limit"
_BZ_PDG_2025 = "PDG2025:CR004:B_bZ_pair_mass_limit"
_TW_PDG_2025 = "PDG2025:CR004:B_tW_pair_mass_limit"
_KNOWN_VALUE_IDS = (
    _ACTIVE_BH_PDG_2025,
    _BZ_PDG_2025,
    _TW_PDG_2025,
)
_BUDGET_POLICY = (
    "active HARD budget is the strongest PDG2025 pure-channel bottom-partner "
    "pair-production mass lower limit in CR004.yaml, B(B -> H b)=1. The "
    "B -> Z b and B -> W t pure-channel limits are retained as diagnostics "
    "because ParameterPoint does not provide B branching fractions."
)


@dataclass(frozen=True)
class CR004ValueAnchor:
    """Typed CR004 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    experiment: str | None
    source_key: str | None
    underlying_source_key: str | None
    production_assumption: str
    branching_assumption: str

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
class CR004Anchor:
    """YAML-loaded CR004 charge -1/3 bottom-partner mass-limit bundle."""

    active_limit: CR004ValueAnchor
    b_h: CR004ValueAnchor
    b_z: CR004ValueAnchor
    t_w: CR004ValueAnchor
    all_limits: tuple[CR004ValueAnchor, ...]
    parent_canonical_source: str | None
    parent_note: str | None
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
            f"{process_id}: CR004 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR004 field {field_name!r}={value!r} is not finite"
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


def _required_assumption(
    entry: Mapping[str, Any],
    key: str,
    *,
    process_id: str,
    value_id: str,
) -> str:
    assumptions = entry.get("assumptions")
    if not isinstance(assumptions, Mapping):
        raise AnchorError(f"{process_id}: {value_id}.assumptions must be a mapping")
    value = assumptions.get(key)
    if value is None:
        raise AnchorError(f"{process_id}: {value_id}.assumptions.{key} is missing")
    return str(value)


def _load_mass_limit_anchor(value_id: str, *, process_id: str) -> CR004ValueAnchor:
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
            f"expected {block_key!r} for CR004 value_id {value_id!r}"
        )

    return CR004ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        source_key=_optional_str(entry.get("source_key")),
        underlying_source_key=_optional_str(entry.get("underlying_source_key")),
        production_assumption=_required_assumption(
            entry,
            "production",
            process_id=process_id,
            value_id=value_id,
        ),
        branching_assumption=_required_assumption(
            entry,
            "branching_fraction",
            process_id=process_id,
            value_id=value_id,
        ),
    )


def _load_cr004_anchor(process_id: str) -> CR004Anchor:
    limits = tuple(
        _load_mass_limit_anchor(value_id, process_id=process_id)
        for value_id in _KNOWN_VALUE_IDS
    )
    by_id = {limit.value_id: limit for limit in limits}
    anchor = CR004Anchor(
        active_limit=by_id[_ACTIVE_BH_PDG_2025],
        b_h=by_id[_ACTIVE_BH_PDG_2025],
        b_z=by_id[_BZ_PDG_2025],
        t_w=by_id[_TW_PDG_2025],
        all_limits=limits,
        parent_canonical_source=_optional_str(_parent(process_id).get("canonical_source")),
        parent_note=_optional_str(_parent(process_id).get("note")),
        budget_policy=_BUDGET_POLICY,
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: CR004 active budget must be positive")
    return anchor


def _limit_from_anchor(anchor: CR004ValueAnchor, *, process_id: str) -> ColliderResonanceLimit:
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
        benchmark_model=anchor.branching_assumption,
        diagnostics={
            "value_id": anchor.value_id,
            "entry_index": anchor.entry_index,
            "display": anchor.display,
            "experiment": anchor.experiment,
            "source_key": anchor.source_key,
            "underlying_source_key": anchor.underlying_source_key,
            "production_assumption": anchor.production_assumption,
            "branching_assumption": anchor.branching_assumption,
            "snapshot_path": anchor.snapshot_path,
            "yaml_value_gev": float(anchor.value_gev),
            "yaml_units": anchor.units,
        },
    )


@register
class Constraint:
    """Catalogued charge -1/3 vector-like bottom-partner mass constraint."""

    process_id = "CR004"
    severity = Severity.HARD
    observable = "m(B pair -> tW/bZ/bH)"

    def __init__(self) -> None:
        self.anchor = _load_cr004_anchor(self.process_id)
        self.limit = _limit_from_anchor(
            self.anchor.active_limit,
            process_id=self.process_id,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        try:
            m_kk_gev, mass_source = resolve_charge_minus_one_third_vlq_mkk_gev(
                couplings=point.get_extra(_REQUIRED_EXTRA),
            )
        except (TypeError, ValueError, AttributeError) as exc:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=False,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=f"invalid charge -1/3 VLQ proxy mass input; CR004 not evaluated: {exc}",
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "required_mass_source": f"{_REQUIRED_EXTRA}.M_KK",
                    "needs_human_physics": VLQ_B_PAIR_MASS_PROXY_ASSUMPTION_V1,
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
                    f"extra {_REQUIRED_EXTRA!r} absent; charge -1/3 B "
                    "pair-production constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "required_mass_source": f"{_REQUIRED_EXTRA}.M_KK",
                    "budget_policy": self.anchor.budget_policy,
                    "active_value_id": self.anchor.active_limit.value_id,
                    "needs_human_physics": VLQ_B_PAIR_MASS_PROXY_ASSUMPTION_V1,
                },
            )

        prediction = vlq_b_pair_prediction_from_m_kk_gev(
            float(m_kk_gev),
            resonance=_RESONANCE,
            final_state=_FINAL_STATE,
        )
        comparison = evaluate_collider_resonance_limit(prediction, self.limit)
        diagnostics = dict(comparison.diagnostics)
        diagnostics.update(
            {
                "m_kk_gev": float(m_kk_gev),
                "m_b_partner_proxy_gev": float(m_kk_gev),
                "m_b_partner_proxy_tev": float(comparison.predicted_mass_tev),
                "mass_proxy": "m_B = M_KK",
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
                "parent_note": self.anchor.parent_note,
                "all_mass_limits_gev": {
                    limit.value_id: float(limit.value_gev)
                    for limit in self.anchor.all_limits
                },
                "b_h_limit_gev": float(self.anchor.b_h.value_gev),
                "b_z_limit_gev": float(self.anchor.b_z.value_gev),
                "t_w_limit_gev": float(self.anchor.t_w.value_gev),
                "needs_human_physics": VLQ_B_PAIR_MASS_PROXY_ASSUMPTION_V1,
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
                "Charge -1/3 vector-like B pair production uses the documented "
                "mass proxy m_B = M_KK and the active PDG2025 B -> H b "
                "pure-channel mass lower bound. HARD ratio is m_limit/m_B; "
                "full sigma*BR and branching-simplex recast is marked "
                "NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
