"""T002 - top FCNC decay ``t -> u Z``.

Physics
-------
The observable is the branching fraction ``BR(t -> u Z)``.  The rigorous
decay-width part is evaluated by the shared top-FCNC module,
``quarkConstraints.top_fcnc``, reached only through the
``flavor_catalog_constraints.physics_adapters.top_fcnc`` boundary.  The
effective vector-coupling convention is

    L = (g / 2 c_W) Z_mu ubar gamma^mu (X_L P_L + X_R P_R) t + h.c.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete RS prediction needs electroweak KK/Z/Z'
mixing, top-sector rotations, correlated SMEFT operators, and the collider
interpretation of ATLAS tensor-benchmark limits.  The current ParameterPoint
only carries quark mass-basis KK-gluon-style couplings, so v1 uses the
documented top-Z proxy from ``quarkConstraints.top_fcnc``:
``X_ut = (g_ut / g_s) * m_Z^2 / M_KK^2``.

Severity
--------
HARD.  The SM branching fraction is negligible at the current collider limit,
so the HARD ratio is the RS-proxy NP branching fraction divided by the active
95% CL upper limit loaded from T002.yaml.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/T002.yaml`` is the source of truth for
the PDG/ATLAS limits and the SM reference value.  T002 stores its values in a
``pdg_or_equivalent.values`` list and writes several values as strings with
``<``; this module adapts selected list entries into the scaffold
``load_anchor`` path after parsing the numeric bound, and fails loudly if an
expected value ID is missing or malformed.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.top_fcnc import (
    TOP_FCNC_RS_Z_MATCHING_ASSUMPTION_V1,
    TOP_FCNC_Z_VECTOR_CONVENTION,
    t_to_q_z_from_couplings,
    top_fcnc_default_sm_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_LIGHT_QUARK = "u"
_LIGHT_UP_INDEX = 0

_PDG_TZU_LEFT = "PDG2025:T002:tZu_left"
_PDG_TZU_RIGHT = "PDG2025:T002:tZu_right"
_ATLAS_TZU_LEFT = "ATLAS2023:T002:tZu_left"
_ATLAS_TZU_RIGHT = "ATLAS2023:T002:tZu_right"
_CMS_TZU = "CMS2017:T002:tZu"
_SM_AGUILAR_SAAVEDRA = "AguilarSaavedra2004:T002:SM"
_EXPECTED_LIMIT_UNITS = "branching fraction, 95% CL upper limit"
_BUDGET_POLICY = (
    "active budget is the numerically strongest current PDG2025 t->Zu "
    "benchmark limit in T002.yaml; left/right benchmark values are retained "
    "separately in diagnostics."
)
_NUMBER_RE = re.compile(
    r"(?P<number>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(?P<pct>%)?"
)


@dataclass(frozen=True)
class T002ValueAnchor:
    """Typed T002 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    raw_value: str
    is_upper_limit: bool

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def block_key(self) -> str:
        return self.anchor.block_key

    @property
    def source_url(self) -> str | None:
        return self.anchor.source_url

    @property
    def year(self) -> int | None:
        return self.anchor.year


@dataclass(frozen=True)
class T002Anchor:
    """All YAML-loaded T002 anchors used by the constraint."""

    pdg_tzu_left: T002ValueAnchor
    pdg_tzu_right: T002ValueAnchor
    atlas_tzu_left: T002ValueAnchor
    atlas_tzu_right: T002ValueAnchor
    cms_tzu: T002ValueAnchor
    sm_aguilar_saavedra: T002ValueAnchor
    active_limit: T002ValueAnchor
    budget_policy: str

    @property
    def value(self) -> float:
        return self.active_limit.value

    @property
    def budget(self) -> float:
        return self.active_limit.value

    @property
    def sm_value(self) -> float:
        return self.sm_aguilar_saavedra.value

    @property
    def source_url(self) -> str | None:
        return self.active_limit.source_url


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: T002 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: T002 field {field_name!r}={value!r} is not finite"
        )
    return number


def _parse_numeric_value(value: Any, *, process_id: str, field_name: str) -> float:
    if isinstance(value, (int, float)):
        return _required_float(value, process_id=process_id, field_name=field_name)
    text = str(value).strip()
    match = _NUMBER_RE.search(text)
    if match is None:
        raise AnchorError(
            f"{process_id}: cannot parse numeric value from {field_name}={value!r}"
        )
    number = _required_float(
        match.group("number"),
        process_id=process_id,
        field_name=field_name,
    )
    if match.group("pct"):
        number *= 1.0e-2
    if number < 0.0:
        raise AnchorError(f"{process_id}: {field_name} must be non-negative")
    return float(number)


def _pdg_values(process_id: str) -> Sequence[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg = data.get("pdg_or_equivalent")
    if not isinstance(pdg, Mapping):
        raise AnchorError(f"{process_id}: expected mapping-shaped pdg_or_equivalent")
    values = pdg.get("values")
    if not isinstance(values, list) or not values:
        raise AnchorError(f"{process_id}: expected non-empty pdg_or_equivalent.values")
    for index, entry in enumerate(values):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent.values[{index}] is not a mapping"
            )
    return values


def _pdg_parent_metadata(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg = data.get("pdg_or_equivalent")
    if not isinstance(pdg, Mapping):
        raise AnchorError(f"{process_id}: expected mapping-shaped pdg_or_equivalent")
    return pdg


def _entry_by_value_id(
    process_id: str,
    value_id: str,
) -> tuple[int, Mapping[str, Any]]:
    matches: list[tuple[int, Mapping[str, Any]]] = []
    for index, entry in enumerate(_pdg_values(process_id)):
        if entry.get("value_id") == value_id:
            matches.append((index, entry))
    if not matches:
        present = [str(entry.get("value_id")) for entry in _pdg_values(process_id)]
        raise AnchorError(
            f"{process_id}: value_id {value_id!r} not found in "
            f"pdg_or_equivalent.values (present: {present})"
        )
    if len(matches) > 1:
        raise AnchorError(f"{process_id}: duplicate value_id {value_id!r}")
    return matches[0]


def _load_scaffold_value_anchor(
    value_id: str,
    *,
    process_id: str,
    require_upper_limit: bool,
) -> T002ValueAnchor:
    index, entry = _entry_by_value_id(process_id, value_id)
    parent = _pdg_parent_metadata(process_id)
    raw_value = str(entry.get("value"))
    value_source = entry.get("normalized_value", entry.get("value"))
    numeric_value = _parse_numeric_value(
        value_source,
        process_id=process_id,
        field_name=f"{value_id}.value",
    )
    is_upper_limit = raw_value.strip().startswith("<") or str(
        entry.get("normalized_value", "")
    ).strip().startswith("<")
    if require_upper_limit and not is_upper_limit:
        raise AnchorError(f"{process_id}: {value_id} must be an upper-limit entry")
    if require_upper_limit and entry.get("units") != _EXPECTED_LIMIT_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_LIMIT_UNITS!r} for "
            f"{value_id}, got {entry.get('units')!r}"
        )

    block_key = f"pdg_or_equivalent.values[{index}]"
    virtual_entry = dict(entry)
    virtual_entry["value"] = numeric_value
    virtual_entry.setdefault("source", parent.get("source"))
    virtual_entry.setdefault("snapshot_path", parent.get("snapshot_path"))
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
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for T002 value_id {value_id!r}"
        )
    return T002ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        raw_value=raw_value,
        is_upper_limit=is_upper_limit,
    )


def _load_limit_anchor(value_id: str, *, process_id: str) -> T002ValueAnchor:
    return _load_scaffold_value_anchor(
        value_id,
        process_id=process_id,
        require_upper_limit=True,
    )


def _load_reference_anchor(value_id: str, *, process_id: str) -> T002ValueAnchor:
    return _load_scaffold_value_anchor(
        value_id,
        process_id=process_id,
        require_upper_limit=False,
    )


def _load_t002_anchor(process_id: str) -> T002Anchor:
    pdg_left = _load_limit_anchor(_PDG_TZU_LEFT, process_id=process_id)
    pdg_right = _load_limit_anchor(_PDG_TZU_RIGHT, process_id=process_id)
    atlas_left = _load_limit_anchor(_ATLAS_TZU_LEFT, process_id=process_id)
    atlas_right = _load_limit_anchor(_ATLAS_TZU_RIGHT, process_id=process_id)
    cms = _load_limit_anchor(_CMS_TZU, process_id=process_id)
    sm_aguilar = _load_reference_anchor(_SM_AGUILAR_SAAVEDRA, process_id=process_id)
    active = min((pdg_left, pdg_right), key=lambda item: item.value)
    if active.value <= 0.0 or sm_aguilar.value <= 0.0:
        raise AnchorError(f"{process_id}: T002 budget and SM estimate must be positive")
    return T002Anchor(
        pdg_tzu_left=pdg_left,
        pdg_tzu_right=pdg_right,
        atlas_tzu_left=atlas_left,
        atlas_tzu_right=atlas_right,
        cms_tzu=cms,
        sm_aguilar_saavedra=sm_aguilar,
        active_limit=active,
        budget_policy=_BUDGET_POLICY,
    )


@register
class Constraint:
    """Catalogued ``BR(t -> u Z)`` top-FCNC constraint (process_id T002)."""

    process_id = "T002"
    severity = Severity.HARD
    observable = "BR(t -> u Z)"

    def __init__(self) -> None:
        self.anchor = _load_t002_anchor(self.process_id)
        self.sm_inputs = top_fcnc_default_sm_inputs()
        self.sm_value = self.anchor.sm_value

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_value),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; t -> u Z constraint "
                    "was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "needs_human_physics": TOP_FCNC_RS_Z_MATCHING_ASSUMPTION_V1,
                    "budget_policy": self.anchor.budget_policy,
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        result = t_to_q_z_from_couplings(
            couplings,
            light_quark=_LIGHT_QUARK,
            light_up_index=_LIGHT_UP_INDEX,
            m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
            inputs=self.sm_inputs,
        )
        predicted_np = float(result.branching_fraction)
        budget = float(self.anchor.budget)
        ratio = float(predicted_np / budget) if budget > 0.0 else float("inf")
        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "sm_anchor_branching_fraction": float(self.sm_value),
                "total_branching_fraction_including_sm": float(
                    predicted_np + self.sm_value
                ),
                "np_branching_fraction": predicted_np,
                "light_quark": _LIGHT_QUARK,
                "light_up_index": _LIGHT_UP_INDEX,
                "active_limit_value_id": self.anchor.active_limit.value_id,
                "active_limit_block": self.anchor.active_limit.block_key,
                "pdg_tzu_left_limit": float(self.anchor.pdg_tzu_left.value),
                "pdg_tzu_right_limit": float(self.anchor.pdg_tzu_right.value),
                "atlas_tzu_left_limit": float(self.anchor.atlas_tzu_left.value),
                "atlas_tzu_right_limit": float(self.anchor.atlas_tzu_right.value),
                "cms_tzu_limit": float(self.anchor.cms_tzu.value),
                "budget_policy": self.anchor.budget_policy,
                "operator_convention": TOP_FCNC_Z_VECTOR_CONVENTION,
                "rs_matching_assumption": TOP_FCNC_RS_Z_MATCHING_ASSUMPTION_V1,
                "needs_human_physics": TOP_FCNC_RS_Z_MATCHING_ASSUMPTION_V1,
                # Structured tag-class hint consumed by run_full_catalog_scan
                # tag_result (M5): the RS t->qZ contribution is a documented
                # overlap proxy.  This replaces brittle prose substring matching
                # (the legacy "proxy" in needs_text missed the plural "proxies").
                "tag_class": "proxy",
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "sm_is_negligible_for_limit": True,
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(ratio <= 1.0),
            predicted=predicted_np,
            sm_prediction=float(self.sm_value),
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "BR(t -> u Z) uses the shared top-FCNC vector-width formula. "
                "The RS contribution is a documented top-Z overlap proxy and "
                "is flagged NEEDS-HUMAN-PHYSICS; the HARD budget is the active "
                "95% CL upper limit from T002.yaml."
            ),
            diagnostics=diagnostics,
        )
