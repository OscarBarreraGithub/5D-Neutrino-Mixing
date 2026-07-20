"""T004 - top FCNC decay ``t -> u gamma``.

Physics
-------
The observable is the branching fraction ``BR(t -> u gamma)``.  The rigorous
two-body photon-dipole width is reused from the shared top-FCNC module,
``quarkConstraints.top_fcnc``, reached only through the
``flavor_catalog_constraints.physics_adapters.top_fcnc`` boundary.  The
effective dipole convention is

    L = e A_mu ubar i sigma^{mu nu} k_nu / m_t
        (lambda_L P_L + lambda_R P_R) t + h.c.

This is the up-quark analogue of T003 and deliberately reuses the same
photon-dipole machinery.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete RS prediction needs the electromagnetic
dipole operator basis, loop/Yukawa matching, electroweak KK photon/Z mixing,
top-sector rotations, and a collider interpretation of the production-plus-
decay searches.  The current ParameterPoint only carries quark mass-basis
KK-gluon-style couplings, so v1 uses the documented diagnostic proxy
``lambda_ut = (g_ut / g_s) * m_t^2 / M_KK^2``.  This is not a human-approved
RS electromagnetic dipole calculation.

Severity
--------
HARD.  The SM contribution is negligible at the current collider limit and
T004.yaml provides only textual SM context, so the HARD ratio is the RS-proxy
pure-NP branching fraction divided by the active 95% CL collider upper limit
loaded from T004.yaml.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/T004.yaml`` is the source of truth for
the PDG, ATLAS, and CMS limits plus the paper-era SM/projection context.
T004 stores its numeric limits in a value list under ``pdg_or_equivalent``;
this module adapts selected list entries into the scaffold ``load_anchor``
path after parsing numeric values, and fails loudly if an expected value ID is
missing or malformed.
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
    TOP_FCNC_PHOTON_DIPOLE_CONVENTION,
    TOP_FCNC_RS_PHOTON_DIPOLE_PROXY_ASSUMPTION_V1,
    t_to_q_gamma_from_couplings,
    top_fcnc_default_sm_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_DIPOLE_MASS_EXTRA = "kk_ew_mass_gev"
_LIGHT_QUARK = "u"
_LIGHT_UP_INDEX = 0
_PDG_PARENT = "pdg_or_equivalent"
_PAPER_REFERENCES_PARENT = "paper_era_reference"

_PDG_GAMMAQ_COMBINED = "PDG2026:T004:tgammaq_combined"
_ATLAS_TUGAMMA_LEFT = "ATLAS2023:T004:tugamma_left"
_ATLAS_TUGAMMA_RIGHT = "ATLAS2023:T004:tugamma_right"
_CMS_TUGAMMA_OBSERVED = "CMS2024:T004:tugamma_observed"
_SM_PROJECTION_CONTEXT = "AguilarSaavedra2017:T004:SM_and_projection_context"

_EXPECTED_LIMIT_UNITS = "branching fraction, 95% CL upper limit"
_SCAFFOLD_UNCERTAINTY_KEY = "__t004_uncertainty_is_not_used__"
_BUDGET_POLICY = (
    "conservative chirality-tight policy: active budget intentionally uses the "
    "ATLAS2023 left-handed t-u-gamma benchmark, the tightest flavor-separated "
    "T004 collider limit in T004.yaml, for all proxy points including pure-RH "
    "points; ATLAS right-handed, CMS observed, and PDG generic gamma-q entries "
    "are retained separately in diagnostics."
)
_NUMBER_RE = re.compile(
    r"(?P<number>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(?P<pct>%)?"
)


@dataclass(frozen=True)
class T004ValueAnchor:
    """Typed T004 value entry routed through the scaffold ``load_anchor``."""

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
class T004TextReference:
    """Typed view over a non-numeric paper-era context entry."""

    value_id: str
    year: int | None
    source: str | None
    role: str | None
    values: tuple[str, ...]
    snapshot_path: str | None
    source_url: str | None
    sha256: str | None


@dataclass(frozen=True)
class T004Anchor:
    """All YAML-loaded T004 anchors used by the constraint."""

    pdg_gammaq_combined: T004ValueAnchor
    atlas_tugamma_left: T004ValueAnchor
    atlas_tugamma_right: T004ValueAnchor
    cms_tugamma_observed: T004ValueAnchor
    sm_projection_context: T004TextReference
    active_limit: T004ValueAnchor
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
            f"{process_id}: T004 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: T004 field {field_name!r}={value!r} is not finite"
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
) -> T004ValueAnchor:
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
            f"expected {block_key!r} for T004 value_id {value_id!r}"
        )
    return T004ValueAnchor(
        anchor=scaffold_anchor,
        parent_key=parent_key,
        value_id=value_id,
        entry_index=index,
        raw_value=raw_value,
        is_upper_limit=is_upper_limit,
    )


def _load_limit_anchor(value_id: str, *, process_id: str) -> T004ValueAnchor:
    return _load_scaffold_value_anchor(
        value_id,
        process_id=process_id,
        parent_key=_PDG_PARENT,
        require_upper_limit=True,
        expected_units=_EXPECTED_LIMIT_UNITS,
    )


def _paper_reference_entry(
    process_id: str,
    value_id: str,
) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    entries = data.get(_PAPER_REFERENCES_PARENT)
    if not isinstance(entries, list) or not entries:
        raise AnchorError(
            f"{process_id}: expected non-empty list-shaped {_PAPER_REFERENCES_PARENT}"
        )
    matches: list[Mapping[str, Any]] = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: {_PAPER_REFERENCES_PARENT}[{index}] is not a mapping"
            )
        if entry.get("value_id") == value_id:
            matches.append(entry)
    if not matches:
        present = [str(entry.get("value_id")) for entry in entries if isinstance(entry, Mapping)]
        raise AnchorError(
            f"{process_id}: value_id {value_id!r} not found in "
            f"{_PAPER_REFERENCES_PARENT} (present: {present})"
        )
    if len(matches) > 1:
        raise AnchorError(f"{process_id}: duplicate value_id {value_id!r}")
    return matches[0]


def _load_text_reference(value_id: str, *, process_id: str) -> T004TextReference:
    entry = _paper_reference_entry(process_id, value_id)
    values = entry.get("values")
    if not isinstance(values, list) or not values:
        raise AnchorError(f"{process_id}: {value_id} must provide non-empty values")
    return T004TextReference(
        value_id=str(entry.get("value_id")),
        year=_optional_int(entry.get("year")),
        source=_optional_str(entry.get("source")),
        role=_optional_str(entry.get("role")),
        values=tuple(str(value) for value in values),
        snapshot_path=_optional_str(entry.get("snapshot_path")),
        source_url=_optional_str(entry.get("source_url")),
        sha256=_optional_str(entry.get("sha256")),
    )


def _load_t004_anchor(process_id: str) -> T004Anchor:
    pdg_gammaq = _load_limit_anchor(_PDG_GAMMAQ_COMBINED, process_id=process_id)
    atlas_left = _load_limit_anchor(_ATLAS_TUGAMMA_LEFT, process_id=process_id)
    atlas_right = _load_limit_anchor(_ATLAS_TUGAMMA_RIGHT, process_id=process_id)
    cms_observed = _load_limit_anchor(_CMS_TUGAMMA_OBSERVED, process_id=process_id)
    sm_context = _load_text_reference(_SM_PROJECTION_CONTEXT, process_id=process_id)
    if (
        pdg_gammaq.value <= 0.0
        or atlas_left.value <= 0.0
        or atlas_right.value <= 0.0
        or cms_observed.value <= 0.0
    ):
        raise AnchorError(f"{process_id}: T004 collider limits must be positive")
    return T004Anchor(
        pdg_gammaq_combined=pdg_gammaq,
        atlas_tugamma_left=atlas_left,
        atlas_tugamma_right=atlas_right,
        cms_tugamma_observed=cms_observed,
        sm_projection_context=sm_context,
        active_limit=atlas_left,
        budget_policy=_BUDGET_POLICY,
    )


@register
class Constraint:
    """Catalogued ``BR(t -> u gamma)`` top-FCNC constraint (process_id T004)."""

    process_id = "T004"
    severity = Severity.HARD
    observable = "BR(t -> u gamma)"

    def __init__(self) -> None:
        self.anchor = _load_t004_anchor(self.process_id)
        self.sm_inputs = top_fcnc_default_sm_inputs()

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=None,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; t -> u gamma constraint "
                    "was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "active_limit_value_id": self.anchor.active_limit.value_id,
                    "active_limit_block": self.anchor.active_limit.block_key,
                    "atlas_tugamma_left_limit": float(
                        self.anchor.atlas_tugamma_left.value
                    ),
                    "atlas_tugamma_right_limit": float(
                        self.anchor.atlas_tugamma_right.value
                    ),
                    "cms_tugamma_observed_limit": float(
                        self.anchor.cms_tugamma_observed.value
                    ),
                    "pdg_gammaq_combined_limit_context": float(
                        self.anchor.pdg_gammaq_combined.value
                    ),
                    "needs_human_physics": (
                        TOP_FCNC_RS_PHOTON_DIPOLE_PROXY_ASSUMPTION_V1
                    ),
                    "budget_policy": self.anchor.budget_policy,
                    "sm_prediction_policy": (
                        "T004 is evaluated as a pure-NP collider-limit bound; "
                        "T004.yaml gives textual SM-negligibility context but "
                        "no numeric SM central value."
                    ),
                    "sm_context_values": self.anchor.sm_projection_context.values,
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
                "np_branching_fraction": predicted_np,
                "light_quark": _LIGHT_QUARK,
                "light_up_index": _LIGHT_UP_INDEX,
                "active_limit_value_id": self.anchor.active_limit.value_id,
                "active_limit_block": self.anchor.active_limit.block_key,
                "atlas_tugamma_left_limit": float(
                    self.anchor.atlas_tugamma_left.value
                ),
                "atlas_tugamma_right_limit": float(
                    self.anchor.atlas_tugamma_right.value
                ),
                "cms_tugamma_observed_limit": float(
                    self.anchor.cms_tugamma_observed.value
                ),
                "pdg_gammaq_combined_limit_context": float(
                    self.anchor.pdg_gammaq_combined.value
                ),
                "pdg_gammaq_combined_units": self.anchor.pdg_gammaq_combined.units,
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
                "sm_prediction_policy": (
                    "T004 is evaluated as a pure-NP collider-limit bound; "
                    "T004.yaml gives textual SM-negligibility context but "
                    "no numeric SM central value."
                ),
                "sm_context_source": self.anchor.sm_projection_context.source,
                "sm_context_values": self.anchor.sm_projection_context.values,
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(ratio <= 1.0),
            predicted=predicted_np,
            sm_prediction=None,
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "BR(t -> u gamma) uses the shared top-FCNC photon-dipole "
                "width. The RS contribution is a documented dipole proxy and "
                "is flagged NEEDS-HUMAN-PHYSICS; the HARD budget is the "
                "active 95% CL collider limit from T004.yaml."
            ),
            diagnostics=diagnostics,
        )
