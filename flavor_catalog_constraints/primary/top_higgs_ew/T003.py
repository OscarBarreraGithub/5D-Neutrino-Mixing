"""T003 - top FCNC decay ``t -> c gamma``.

Physics
-------
The observable is the branching fraction ``BR(t -> c gamma)``.  The rigorous
two-body photon-dipole width is evaluated by the shared top-FCNC module,
``quarkConstraints.top_fcnc``, reached only through the
``flavor_catalog_constraints.physics_adapters.top_fcnc`` boundary.  The
effective dipole convention is

    L = e A_mu cbar i sigma^{mu nu} k_nu / m_t
        (lambda_L P_L + lambda_R P_R) t + h.c.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete RS prediction needs the electromagnetic
dipole operator basis, loop/Yukawa matching, electroweak KK photon/Z mixing,
top-sector rotations, and a collider interpretation of the production-plus-
decay searches.  The current ParameterPoint only carries quark mass-basis
KK-gluon-style couplings, so v1 uses the documented diagnostic proxy
``lambda_ct = (g_ct / g_s) * m_t^2 / M_KK^2``.  This is not a human-approved
RS electromagnetic dipole calculation.

Severity
--------
HARD.  The SM branching fraction is negligible at the current collider limit,
so the HARD ratio is the RS-proxy NP branching fraction divided by the
CMS2024 observed 95% CL upper limit loaded from T003.yaml.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/T003.yaml`` is the source of truth for
the CMS/ATLAS/PDG limits, the Aguilar-Saavedra SM reference value, and the
paper-era dipole normalization.  T003 stores its values in value lists under
``pdg_or_equivalent`` and ``paper_era_reference``; this module adapts selected
list entries into the scaffold ``load_anchor`` path after parsing numeric
values, and fails loudly if an expected value ID is missing or malformed.
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
    TOP_FCNC_PHOTON_DIPOLE_CONVENTION,
    TOP_FCNC_RS_PHOTON_DIPOLE_PROXY_ASSUMPTION_V1,
    t_to_q_gamma_from_couplings,
    top_fcnc_default_sm_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_DIPOLE_MASS_EXTRA = "kk_ew_mass_gev"
_LIGHT_QUARK = "c"
_LIGHT_UP_INDEX = 1
_PDG_PARENT = "pdg_or_equivalent"
_PAPER_PARENT = "paper_era_reference"

_CMS_TCGAMMA = "CMS2024:T003:tcgamma"
_ATLAS_TCGAMMA_LEFT = "ATLAS2023:T003:tcgamma_left"
_ATLAS_TCGAMMA_RIGHT = "ATLAS2023:T003:tcgamma_right"
_PDG_GAMMAQ_SUMMARY = "PDG2026:T003:gammaq_summary"
_SM_AGUILAR_SAAVEDRA = "AguilarSaavedra2004:T003:SM"
_DIPOLE_NORMALIZATION = "AguilarSaavedra2004:T003:dipole_normalization"

_EXPECTED_LIMIT_UNITS = "branching fraction, 95% CL upper limit"
_EXPECTED_GENERIC_RATIO_UNITS = "ratio, 95% CL upper limit"
_EXPECTED_SM_UNITS = "branching fraction"
_EXPECTED_DIPOLE_NORMALIZATION_UNITS = (
    "branching-fraction normalization in the source convention"
)
_SCAFFOLD_UNCERTAINTY_KEY = "__t003_uncertainty_is_not_used__"
_BUDGET_POLICY = (
    "active budget is the CMS2024 observed charm-separated t->c gamma "
    "single-nonzero-coupling limit in T003.yaml; ATLAS chirality benchmarks "
    "and the PDG generic gamma-q ratio are retained as context only."
)
_NUMBER_RE = re.compile(
    r"(?P<number>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(?P<pct>%)?"
)


@dataclass(frozen=True)
class T003ValueAnchor:
    """Typed T003 value entry routed through the scaffold ``load_anchor``."""

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
class T003Anchor:
    """All YAML-loaded T003 anchors used by the constraint."""

    cms_tcgamma: T003ValueAnchor
    atlas_tcgamma_left: T003ValueAnchor
    atlas_tcgamma_right: T003ValueAnchor
    pdg_gammaq_summary: T003ValueAnchor
    sm_aguilar_saavedra: T003ValueAnchor
    dipole_normalization: T003ValueAnchor
    active_limit: T003ValueAnchor
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
            f"{process_id}: T003 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: T003 field {field_name!r}={value!r} is not finite"
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


def _value_entries(process_id: str, parent_key: str) -> Sequence[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY)
    parent = data.get(parent_key)
    if not isinstance(parent, Mapping):
        raise AnchorError(f"{process_id}: expected mapping-shaped {parent_key}")
    values = parent.get("values")
    if not isinstance(values, list) or not values:
        raise AnchorError(f"{process_id}: expected non-empty {parent_key}.values")
    for index, entry in enumerate(values):
        if not isinstance(entry, Mapping):
            raise AnchorError(f"{process_id}: {parent_key}.values[{index}] is not a mapping")
    return values


def _parent_metadata(process_id: str, parent_key: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    parent = data.get(parent_key)
    if not isinstance(parent, Mapping):
        raise AnchorError(f"{process_id}: expected mapping-shaped {parent_key}")
    return parent


def _entry_by_value_id(
    process_id: str,
    value_id: str,
    *,
    parent_key: str,
) -> tuple[int, Mapping[str, Any]]:
    matches: list[tuple[int, Mapping[str, Any]]] = []
    for index, entry in enumerate(_value_entries(process_id, parent_key)):
        if entry.get("value_id") == value_id:
            matches.append((index, entry))
    if not matches:
        present = [
            str(entry.get("value_id"))
            for entry in _value_entries(process_id, parent_key)
        ]
        raise AnchorError(
            f"{process_id}: value_id {value_id!r} not found in "
            f"{parent_key}.values (present: {present})"
        )
    if len(matches) > 1:
        raise AnchorError(f"{process_id}: duplicate value_id {value_id!r}")
    return matches[0]


def _load_scaffold_value_anchor(
    value_id: str,
    *,
    process_id: str,
    parent_key: str,
    require_upper_limit: bool,
    expected_units: str,
) -> T003ValueAnchor:
    index, entry = _entry_by_value_id(process_id, value_id, parent_key=parent_key)
    parent = _parent_metadata(process_id, parent_key)
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

    block_key = f"{parent_key}.values[{index}]"
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
            f"expected {block_key!r} for T003 value_id {value_id!r}"
        )
    return T003ValueAnchor(
        anchor=scaffold_anchor,
        parent_key=parent_key,
        value_id=value_id,
        entry_index=index,
        raw_value=raw_value,
        is_upper_limit=is_upper_limit,
    )


def _load_limit_anchor(
    value_id: str,
    *,
    process_id: str,
    parent_key: str = _PDG_PARENT,
    expected_units: str = _EXPECTED_LIMIT_UNITS,
) -> T003ValueAnchor:
    return _load_scaffold_value_anchor(
        value_id,
        process_id=process_id,
        parent_key=parent_key,
        require_upper_limit=True,
        expected_units=expected_units,
    )


def _load_reference_anchor(
    value_id: str,
    *,
    process_id: str,
    parent_key: str,
    expected_units: str,
) -> T003ValueAnchor:
    return _load_scaffold_value_anchor(
        value_id,
        process_id=process_id,
        parent_key=parent_key,
        require_upper_limit=False,
        expected_units=expected_units,
    )


def _load_t003_anchor(process_id: str) -> T003Anchor:
    cms_tcgamma = _load_limit_anchor(_CMS_TCGAMMA, process_id=process_id)
    atlas_left = _load_limit_anchor(_ATLAS_TCGAMMA_LEFT, process_id=process_id)
    atlas_right = _load_limit_anchor(_ATLAS_TCGAMMA_RIGHT, process_id=process_id)
    pdg_gammaq = _load_limit_anchor(
        _PDG_GAMMAQ_SUMMARY,
        process_id=process_id,
        expected_units=_EXPECTED_GENERIC_RATIO_UNITS,
    )
    sm = _load_reference_anchor(
        _SM_AGUILAR_SAAVEDRA,
        process_id=process_id,
        parent_key=_PAPER_PARENT,
        expected_units=_EXPECTED_SM_UNITS,
    )
    dipole_norm = _load_reference_anchor(
        _DIPOLE_NORMALIZATION,
        process_id=process_id,
        parent_key=_PAPER_PARENT,
        expected_units=_EXPECTED_DIPOLE_NORMALIZATION_UNITS,
    )
    if cms_tcgamma.value <= 0.0 or sm.value <= 0.0 or dipole_norm.value <= 0.0:
        raise AnchorError(f"{process_id}: T003 budget, SM, and dipole norm must be positive")
    return T003Anchor(
        cms_tcgamma=cms_tcgamma,
        atlas_tcgamma_left=atlas_left,
        atlas_tcgamma_right=atlas_right,
        pdg_gammaq_summary=pdg_gammaq,
        sm_aguilar_saavedra=sm,
        dipole_normalization=dipole_norm,
        active_limit=cms_tcgamma,
        budget_policy=_BUDGET_POLICY,
    )


@register
class Constraint:
    """Catalogued ``BR(t -> c gamma)`` top-FCNC constraint (process_id T003)."""

    process_id = "T003"
    severity = Severity.HARD
    observable = "BR(t -> c gamma)"

    def __init__(self) -> None:
        self.anchor = _load_t003_anchor(self.process_id)
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
                    f"extra {_REQUIRED_EXTRA!r} absent; t -> c gamma constraint "
                    "was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "needs_human_physics": (
                        TOP_FCNC_RS_PHOTON_DIPOLE_PROXY_ASSUMPTION_V1
                    ),
                    "budget_policy": self.anchor.budget_policy,
                },
            )

        kk_dipole_mass = point.get_extra(_OPTIONAL_DIPOLE_MASS_EXTRA)
        result = t_to_q_gamma_from_couplings(
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
                "paper_dipole_normalization_coefficient": float(
                    self.anchor.dipole_normalization.value
                ),
                "light_quark": _LIGHT_QUARK,
                "light_up_index": _LIGHT_UP_INDEX,
                "active_limit_value_id": self.anchor.active_limit.value_id,
                "active_limit_block": self.anchor.active_limit.block_key,
                "cms_tcgamma_limit": float(self.anchor.cms_tcgamma.value),
                "atlas_tcgamma_left_limit": float(
                    self.anchor.atlas_tcgamma_left.value
                ),
                "atlas_tcgamma_right_limit": float(
                    self.anchor.atlas_tcgamma_right.value
                ),
                "pdg_gammaq_generic_ratio_limit_context": float(
                    self.anchor.pdg_gammaq_summary.value
                ),
                "pdg_gammaq_generic_units": self.anchor.pdg_gammaq_summary.units,
                "budget_policy": self.anchor.budget_policy,
                "operator_convention": TOP_FCNC_PHOTON_DIPOLE_CONVENTION,
                "rs_matching_assumption": (
                    TOP_FCNC_RS_PHOTON_DIPOLE_PROXY_ASSUMPTION_V1
                ),
                "needs_human_physics": (
                    TOP_FCNC_RS_PHOTON_DIPOLE_PROXY_ASSUMPTION_V1
                ),
                "kk_ew_mass_extra_used": kk_dipole_mass is not None,
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
                "BR(t -> c gamma) uses the shared top-FCNC photon-dipole "
                "width. The RS contribution is a documented dipole proxy and "
                "is flagged NEEDS-HUMAN-PHYSICS; the HARD budget is the "
                "CMS2024 observed 95% CL upper limit from T003.yaml."
            ),
            diagnostics=diagnostics,
        )
