"""T005 - top FCNC decay ``t -> c g``.

Physics
-------
The observable is the branching fraction ``BR(t -> c g)``.  The rigorous
two-body chromomagnetic-dipole width is evaluated by the shared top-FCNC
module, ``quarkConstraints.top_fcnc``, reached only through the
``flavor_catalog_constraints.physics_adapters.top_fcnc`` boundary.  The
effective dipole convention is

    L = g_s G^a_mu cbar T^a i sigma^{mu nu} k_nu / m_t
        (zeta_L P_L + zeta_R P_R) t + h.c.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete RS prediction needs the chromomagnetic
operator basis, loop/Yukawa matching, QCD running, top-sector rotations,
SMEFT ``C_uG`` normalization, and the production-vs-decay collider recast.
The current ``ParameterPoint`` only carries quark mass-basis KK-gluon-style
couplings, so v1 uses the documented diagnostic proxy
``zeta_ct = (g_ct / g_s) * m_t^2 / M_KK^2``.  This is not a human-approved RS
chromomagnetic dipole calculation.

Severity
--------
HARD.  The SM branching fraction is negligible at the current collider limit,
so the HARD ratio is the RS-proxy NP branching fraction divided by the active
95% CL upper limit loaded from T005.yaml.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/T005.yaml`` is the source of truth for
the PDG/ATLAS/CMS limits and the Aguilar-Saavedra SM reference value.  T005
stores its values in a ``pdg_or_equivalent.values`` list; this module adapts
selected list entries into the scaffold ``load_anchor`` path after parsing
numeric values, and fails loudly if an expected value ID is missing or
malformed.
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
    TOP_FCNC_GLUON_DIPOLE_CONVENTION,
    TOP_FCNC_RS_GLUON_DIPOLE_PROXY_ASSUMPTION_V1,
    t_to_q_gluon_from_couplings,
    top_fcnc_default_sm_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_DIPOLE_MASS_EXTRA = "kk_gluon_mass_gev"
_LIGHT_QUARK = "c"
_LIGHT_UP_INDEX = 1
_PARENT = "pdg_or_equivalent"

_ATLAS_TCG = "PDG2025:T005:ATLAS2022:t_cg"
_ATLAS_CUGCT = "PDG2025:T005:ATLAS2022:CuGct"
_ATLAS_CG_TO_T = "PDG2025:T005:ATLAS2022:cg_to_t_cross_section"
_CMS_TCG = "CMS2017:T005:t_cg"
_CMS_KAPPA_TCG = "CMS2017:T005:kappa_tcg"
_SM_AGUILAR_SAAVEDRA = "AguilarSaavedra2004:T005:SM"

_EXPECTED_LIMIT_UNITS = "branching fraction, 95% CL upper limit"
_EXPECTED_CUG_UNITS = "TeV^-2, 95% CL upper limit"
_EXPECTED_XSEC_UNITS = "pb, 95% CL upper limit"
_EXPECTED_KAPPA_UNITS = "TeV^-1, 95% CL upper limit"
_EXPECTED_SM_UNITS = "branching fraction"
_SCAFFOLD_UNCERTAINTY_KEY = "__t005_uncertainty_is_not_used__"
_BUDGET_POLICY = (
    "active budget is the numerically strongest t->c g branching-fraction "
    "95% CL upper limit in T005.yaml; ATLAS C_uG, cg->t, and older CMS "
    "normalizations are retained as context only."
)
_NUMBER_RE = re.compile(
    r"(?P<number>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(?P<pct>%)?"
)


@dataclass(frozen=True)
class T005ValueAnchor:
    """Typed T005 value entry routed through the scaffold ``load_anchor``."""

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
class T005Anchor:
    """All YAML-loaded T005 anchors used by the constraint."""

    atlas_tcg: T005ValueAnchor
    atlas_cugct: T005ValueAnchor
    atlas_cg_to_t_cross_section: T005ValueAnchor
    cms_tcg: T005ValueAnchor
    cms_kappa_tcg: T005ValueAnchor
    sm_aguilar_saavedra: T005ValueAnchor
    active_limit: T005ValueAnchor
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
            f"{process_id}: T005 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: T005 field {field_name!r}={value!r} is not finite"
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
) -> T005ValueAnchor:
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
            f"expected {block_key!r} for T005 value_id {value_id!r}"
        )
    return T005ValueAnchor(
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
    expected_units: str = _EXPECTED_LIMIT_UNITS,
) -> T005ValueAnchor:
    return _load_scaffold_value_anchor(
        value_id,
        process_id=process_id,
        require_upper_limit=True,
        expected_units=expected_units,
    )


def _load_reference_anchor(
    value_id: str,
    *,
    process_id: str,
    expected_units: str,
) -> T005ValueAnchor:
    return _load_scaffold_value_anchor(
        value_id,
        process_id=process_id,
        require_upper_limit=False,
        expected_units=expected_units,
    )


def _load_t005_anchor(process_id: str) -> T005Anchor:
    atlas_tcg = _load_limit_anchor(_ATLAS_TCG, process_id=process_id)
    atlas_cugct = _load_limit_anchor(
        _ATLAS_CUGCT,
        process_id=process_id,
        expected_units=_EXPECTED_CUG_UNITS,
    )
    atlas_cg_to_t = _load_limit_anchor(
        _ATLAS_CG_TO_T,
        process_id=process_id,
        expected_units=_EXPECTED_XSEC_UNITS,
    )
    cms_tcg = _load_limit_anchor(_CMS_TCG, process_id=process_id)
    cms_kappa = _load_limit_anchor(
        _CMS_KAPPA_TCG,
        process_id=process_id,
        expected_units=_EXPECTED_KAPPA_UNITS,
    )
    sm = _load_reference_anchor(
        _SM_AGUILAR_SAAVEDRA,
        process_id=process_id,
        expected_units=_EXPECTED_SM_UNITS,
    )
    active = min((atlas_tcg, cms_tcg), key=lambda item: item.value)
    if active.value <= 0.0 or sm.value <= 0.0:
        raise AnchorError(f"{process_id}: T005 budget and SM estimate must be positive")
    return T005Anchor(
        atlas_tcg=atlas_tcg,
        atlas_cugct=atlas_cugct,
        atlas_cg_to_t_cross_section=atlas_cg_to_t,
        cms_tcg=cms_tcg,
        cms_kappa_tcg=cms_kappa,
        sm_aguilar_saavedra=sm,
        active_limit=active,
        budget_policy=_BUDGET_POLICY,
    )


@register
class Constraint:
    """Catalogued ``BR(t -> c g)`` top-FCNC constraint (process_id T005)."""

    process_id = "T005"
    severity = Severity.HARD
    observable = "BR(t -> c g)"

    def __init__(self) -> None:
        self.anchor = _load_t005_anchor(self.process_id)
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
                    f"extra {_REQUIRED_EXTRA!r} absent; t -> c g constraint "
                    "was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "needs_human_physics": (
                        TOP_FCNC_RS_GLUON_DIPOLE_PROXY_ASSUMPTION_V1
                    ),
                    "budget_policy": self.anchor.budget_policy,
                },
            )

        kk_dipole_mass = point.get_extra(_OPTIONAL_DIPOLE_MASS_EXTRA)
        result = t_to_q_gluon_from_couplings(
            couplings,
            light_quark=_LIGHT_QUARK,
            light_up_index=_LIGHT_UP_INDEX,
            m_kk_gev=None if kk_dipole_mass is None else float(kk_dipole_mass),
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
                "active_limit_source_url": self.anchor.active_limit.source_url,
                "atlas_tcg_limit": float(self.anchor.atlas_tcg.value),
                "atlas_cugct_limit_tev_minus2": float(self.anchor.atlas_cugct.value),
                "atlas_cg_to_t_cross_section_limit_pb": float(
                    self.anchor.atlas_cg_to_t_cross_section.value
                ),
                "cms_tcg_limit": float(self.anchor.cms_tcg.value),
                "cms_kappa_tcg_limit_tev_minus1": float(
                    self.anchor.cms_kappa_tcg.value
                ),
                "budget_policy": self.anchor.budget_policy,
                "operator_convention": TOP_FCNC_GLUON_DIPOLE_CONVENTION,
                "rs_matching_assumption": (
                    TOP_FCNC_RS_GLUON_DIPOLE_PROXY_ASSUMPTION_V1
                ),
                "needs_human_physics": (
                    TOP_FCNC_RS_GLUON_DIPOLE_PROXY_ASSUMPTION_V1
                ),
                "kk_gluon_mass_extra_used": kk_dipole_mass is not None,
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
                "BR(t -> c g) uses the shared top-FCNC chromomagnetic-dipole "
                "width. The RS contribution is a documented dipole proxy and "
                "is flagged NEEDS-HUMAN-PHYSICS; the HARD budget is the active "
                "95% CL upper limit from T005.yaml."
            ),
            diagnostics=diagnostics,
        )
