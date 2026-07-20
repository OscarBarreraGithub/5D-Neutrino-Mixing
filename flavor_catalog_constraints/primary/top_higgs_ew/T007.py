"""T007 - top FCNC decay ``t -> H c``.

Physics
-------
The observable is the branching fraction ``BR(t -> H c)``.  The rigorous
two-body scalar-Yukawa width is reused from the shared top-FCNC module,
``quarkConstraints.top_fcnc``, reached only through the
``flavor_catalog_constraints.physics_adapters.top_fcnc`` boundary.  The
effective Yukawa convention is

    L = -H cbar (y_L P_L + y_R P_R) t + h.c.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete RS prediction needs the Higgs-sector Yukawa
misalignment, top-sector rotations, SMEFT normalization, electroweak KK/Z/Z'
mixing, and a collider prescription combining associated ``tH`` production
with ``ttbar`` FCNC decays.  The current ``ParameterPoint`` only carries quark
mass-basis KK-gluon-style couplings, so v1 uses the documented scalar-Yukawa
proxy ``y_ct = (g_ct / g_s) * m_t^2 / M_KK^2``.  This is not a human-approved
RS top-Higgs FCNC calculation.

Severity
--------
HARD.  The SM contribution is negligible at the current collider limit, so the
HARD ratio is the RS-proxy pure-NP branching fraction divided by the PDG 2026
headline 95% CL upper limit loaded from T007.yaml.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/T007.yaml`` is the source of truth for
the PDG, ATLAS, and CMS limits.  T007 stores its numeric limits in
``pdg_or_equivalent.values``; this module adapts selected list entries into the
scaffold ``load_anchor`` path after parsing numeric values, and fails loudly if
an expected value ID is missing or malformed.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.top_fcnc import (
    TOP_FCNC_HIGGS_SCALAR_CONVENTION,
    TOP_FCNC_RS_HIGGS_YUKAWA_PROXY_ASSUMPTION_V1,
    t_to_q_higgs_from_couplings,
    top_fcnc_default_sm_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_HIGGS_MASS_EXTRA = "kk_ew_mass_gev"
_LIGHT_QUARK = "c"
_LIGHT_UP_INDEX = 1
_PARENT = "pdg_or_equivalent"

_PDG_THC_COMBINED = "PDG2026:T007:tHc_combined"
_ATLAS_THC_MULTILEPTON = "ATLAS2024:T007:tHc_multilepton"
_ATLAS_THC_COMBINED = "ATLAS2024:T007:tHc_combined"
_CMS_THC_SAMESIGN = "CMS2025:T007:tHc_samesign"
_CMS_THC_COMBINED = "CMS2025:T007:tHc_combined"

_EXPECTED_LIMIT_UNITS = "branching fraction, 95% CL upper limit"
_SCAFFOLD_UNCERTAINTY_KEY = "__t007_uncertainty_is_not_used__"
_BUDGET_POLICY = (
    "active budget is the PDG2026 headline t->Hc combined 95% CL upper "
    "limit in T007.yaml; ATLAS and CMS sub/comparison limits are retained "
    "as diagnostics."
)
_SM_POLICY = (
    "T007 is evaluated as a pure-NP collider-limit bound. The zero-effective-"
    "Yukawa rate from the top-FCNC scalar formula is used as the negligible-SM "
    "reference because T007.yaml does not provide a numeric SM central value."
)
_NUMBER_RE = re.compile(
    r"(?P<number>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(?P<pct>%)?"
)


@dataclass(frozen=True)
class T007ValueAnchor:
    """Typed T007 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    parent_key: str
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

    @property
    def units(self) -> str | None:
        return self.anchor.units

    @property
    def snapshot_path(self) -> str | None:
        return self.anchor.snapshot_path


@dataclass(frozen=True)
class T007Anchor:
    """All YAML-loaded T007 anchors used by the constraint."""

    pdg_thc_combined: T007ValueAnchor
    atlas_thc_multilepton: T007ValueAnchor
    atlas_thc_combined: T007ValueAnchor
    cms_thc_samesign: T007ValueAnchor
    cms_thc_combined: T007ValueAnchor
    active_limit: T007ValueAnchor
    budget_policy: str

    @property
    def value(self) -> float:
        return self.active_limit.value

    @property
    def budget(self) -> float:
        return self.active_limit.value

    @property
    def source_url(self) -> str | None:
        return self.active_limit.source_url


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: T007 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: T007 field {field_name!r}={value!r} is not finite"
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


def _value_entries(process_id: str) -> Sequence[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY)
    parent = data.get(_PARENT)
    if not isinstance(parent, Mapping):
        raise AnchorError(f"{process_id}: expected mapping-shaped {_PARENT}")
    values = parent.get("values")
    if not isinstance(values, list) or not values:
        raise AnchorError(f"{process_id}: expected non-empty {_PARENT}.values")
    for index, entry in enumerate(values):
        if not isinstance(entry, Mapping):
            raise AnchorError(f"{process_id}: {_PARENT}.values[{index}] is not a mapping")
    return values


def _parent_metadata(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    parent = data.get(_PARENT)
    if not isinstance(parent, Mapping):
        raise AnchorError(f"{process_id}: expected mapping-shaped {_PARENT}")
    return parent


def _entry_by_value_id(process_id: str, value_id: str) -> tuple[int, Mapping[str, Any]]:
    matches: list[tuple[int, Mapping[str, Any]]] = []
    for index, entry in enumerate(_value_entries(process_id)):
        if entry.get("value_id") == value_id:
            matches.append((index, entry))
    if not matches:
        present = [str(entry.get("value_id")) for entry in _value_entries(process_id)]
        raise AnchorError(
            f"{process_id}: value_id {value_id!r} not found in "
            f"{_PARENT}.values (present: {present})"
        )
    if len(matches) > 1:
        raise AnchorError(f"{process_id}: duplicate value_id {value_id!r}")
    return matches[0]


def _load_scaffold_value_anchor(
    value_id: str,
    *,
    process_id: str,
    require_upper_limit: bool,
    expected_units: str,
) -> T007ValueAnchor:
    index, entry = _entry_by_value_id(process_id, value_id)
    parent = _parent_metadata(process_id)
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
    if entry.get("units") != expected_units:
        raise AnchorError(
            f"{process_id}: expected units {expected_units!r} for "
            f"{value_id}, got {entry.get('units')!r}"
        )

    block_key = f"{_PARENT}.values[{index}]"
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
            uncertainty_key=_SCAFFOLD_UNCERTAINTY_KEY,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for T007 value_id {value_id!r}"
        )
    return T007ValueAnchor(
        anchor=scaffold_anchor,
        parent_key=_PARENT,
        value_id=value_id,
        entry_index=index,
        raw_value=raw_value,
        is_upper_limit=is_upper_limit,
    )


def _load_limit_anchor(
    value_id: str,
    *,
    process_id: str,
) -> T007ValueAnchor:
    return _load_scaffold_value_anchor(
        value_id,
        process_id=process_id,
        require_upper_limit=True,
        expected_units=_EXPECTED_LIMIT_UNITS,
    )


def _load_t007_anchor(process_id: str) -> T007Anchor:
    pdg_thc = _load_limit_anchor(_PDG_THC_COMBINED, process_id=process_id)
    atlas_multilepton = _load_limit_anchor(
        _ATLAS_THC_MULTILEPTON,
        process_id=process_id,
    )
    atlas_combined = _load_limit_anchor(
        _ATLAS_THC_COMBINED,
        process_id=process_id,
    )
    cms_samesign = _load_limit_anchor(_CMS_THC_SAMESIGN, process_id=process_id)
    cms_combined = _load_limit_anchor(_CMS_THC_COMBINED, process_id=process_id)
    if any(
        item.value <= 0.0
        for item in (pdg_thc, atlas_multilepton, atlas_combined, cms_samesign, cms_combined)
    ):
        raise AnchorError(f"{process_id}: T007 collider limits must be positive")
    return T007Anchor(
        pdg_thc_combined=pdg_thc,
        atlas_thc_multilepton=atlas_multilepton,
        atlas_thc_combined=atlas_combined,
        cms_thc_samesign=cms_samesign,
        cms_thc_combined=cms_combined,
        active_limit=pdg_thc,
        budget_policy=_BUDGET_POLICY,
    )


@register
class Constraint:
    """Catalogued ``BR(t -> H c)`` top-Higgs-FCNC constraint (process_id T007)."""

    process_id = "T007"
    severity = Severity.HARD
    observable = "BR(t -> H c)"

    def __init__(self) -> None:
        self.anchor = _load_t007_anchor(self.process_id)
        self.sm_inputs = top_fcnc_default_sm_inputs()
        self.sm_result = t_to_q_higgs_from_couplings(
            None,
            light_quark=_LIGHT_QUARK,
            light_up_index=_LIGHT_UP_INDEX,
            inputs=self.sm_inputs,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_result.branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; t -> H c constraint "
                    "was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "active_limit_value_id": self.anchor.active_limit.value_id,
                    "active_limit_block": self.anchor.active_limit.block_key,
                    "needs_human_physics": (
                        TOP_FCNC_RS_HIGGS_YUKAWA_PROXY_ASSUMPTION_V1
                    ),
                    "budget_policy": self.anchor.budget_policy,
                    "sm_prediction_policy": _SM_POLICY,
                },
            )

        kk_higgs_mass = point.get_extra(_OPTIONAL_HIGGS_MASS_EXTRA)
        result = t_to_q_higgs_from_couplings(
            couplings,
            light_quark=_LIGHT_QUARK,
            light_up_index=_LIGHT_UP_INDEX,
            m_kk_gev=None if kk_higgs_mass is None else float(kk_higgs_mass),
            inputs=self.sm_inputs,
        )
        predicted_np = float(result.branching_fraction)
        budget = float(self.anchor.budget)
        ratio = float(predicted_np / budget) if budget > 0.0 else float("inf")
        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "zero_yukawa_branching_fraction": float(
                    self.sm_result.branching_fraction
                ),
                "np_branching_fraction": predicted_np,
                "total_branching_fraction_including_sm_proxy": float(
                    predicted_np + self.sm_result.branching_fraction
                ),
                "light_quark": _LIGHT_QUARK,
                "light_up_index": _LIGHT_UP_INDEX,
                "active_limit_value_id": self.anchor.active_limit.value_id,
                "active_limit_block": self.anchor.active_limit.block_key,
                "active_limit_source_url": self.anchor.active_limit.source_url,
                "pdg_thc_combined_limit": float(self.anchor.pdg_thc_combined.value),
                "atlas_thc_multilepton_limit": float(
                    self.anchor.atlas_thc_multilepton.value
                ),
                "atlas_thc_combined_limit": float(self.anchor.atlas_thc_combined.value),
                "cms_thc_samesign_limit": float(self.anchor.cms_thc_samesign.value),
                "cms_thc_combined_limit": float(self.anchor.cms_thc_combined.value),
                "budget_policy": self.anchor.budget_policy,
                "operator_convention": TOP_FCNC_HIGGS_SCALAR_CONVENTION,
                "rs_matching_assumption": (
                    TOP_FCNC_RS_HIGGS_YUKAWA_PROXY_ASSUMPTION_V1
                ),
                "needs_human_physics": (
                    TOP_FCNC_RS_HIGGS_YUKAWA_PROXY_ASSUMPTION_V1
                ),
                "kk_ew_mass_extra_used": kk_higgs_mass is not None,
                "sm_is_negligible_for_limit": True,
                "sm_prediction_policy": _SM_POLICY,
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(ratio <= 1.0),
            predicted=predicted_np,
            sm_prediction=float(self.sm_result.branching_fraction),
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "BR(t -> H c) uses the shared top-FCNC scalar-Yukawa width. "
                "The RS contribution is a documented Higgs-Yukawa proxy and "
                "is flagged NEEDS-HUMAN-PHYSICS; the HARD budget is the PDG "
                "2026 headline 95% CL upper limit from T007.yaml."
            ),
            diagnostics=diagnostics,
        )
