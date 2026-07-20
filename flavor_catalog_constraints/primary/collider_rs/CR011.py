"""CR011 - longitudinal vector-boson scattering (VBS).

Physics
-------
CR011 records the fiducial same-sign longitudinal-WW VBS upper limit

    sigma_fid(pp -> jj W_L^+/- W_L^+/-) < 0.45 fb at 95% CL

from the catalog sidecar.  The implementation is an honest hard stub:
there is no first-principles SM/EFT likelihood, no EFT coefficient recast, no
unitarity/form-factor convention, and no RS strong-EWSB matching calculation.
The only optional numerical comparison is against an explicitly
human-supplied fiducial cross-section prediction; the scan does not infer one
from ``M_KK`` or from KK electroweak masses.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS on both sides.  The SM/EFT side needs a
polarization-aware VBS likelihood/recast with basis, profiling, unitarity, and
form-factor choices.  The RS strong-EWSB side needs the custodial
gauge/scalar resonance spectrum, widths, VV branching fractions, SM
interference, and acceptance.  These gaps are reported separately in
``ConstraintResult.diagnostics``.

Severity
--------
INFO.  The measured VBS bound is recorded as a non-vetoing collider-sector
diagnostic.  If a human-supplied fiducial prediction is present, ``passes`` is
the advisory comparison ``sigma_fid / 0.45 fb <= 1``; callers must still treat
INFO as non-vetoing.

Catalog sidecar
---------------
``flavor_catalog/processes/collider_rs/CR011.yaml`` is the source of truth.
It stores numerical limits in ``pdg_or_equivalent.values``; this module routes
selected entries through the scaffold ``load_anchor`` path and fails loudly if
the expected value IDs or units are missing.
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
from flavor_catalog_constraints.physics_adapters.vbs_longitudinal import (
    HUMAN_SUPPLIED_SIGMA_RAW_KEY,
    VBS_LONGITUDINAL_RS_MATCHING_GAP_V1,
    VBS_LONGITUDINAL_SM_EFT_GAP_V1,
    VBSFiducialLimit,
    compare_vbs_fiducial_limit,
    human_prediction_from_raw,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "collider_rs"
_PARENT_KEY = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__cr011_uncertainty_is_not_used__"
_EXPECTED_UNITS = "fb"
_ACTIVE_ATLAS_2025 = "ATLAS2025:CR011:fiducial_sigma_WLWL_same_sign_upper_limit"
_CMS_2020 = "CMS2020:CR011:fiducial_sigma_WLWL_same_sign_upper_limit"
_KNOWN_VALUE_IDS = (_ACTIVE_ATLAS_2025, _CMS_2020)
_BUDGET_POLICY = (
    "active INFO budget is the ATLAS 2025 observed 95% CL fiducial "
    "longitudinal same-sign WW upper limit loaded from CR011.yaml. It is "
    "recorded as a non-vetoing diagnostic unless a human-supplied fiducial "
    "prediction is explicitly attached to ParameterPoint.raw."
)


@dataclass(frozen=True)
class CR011ValueAnchor:
    """Typed CR011 value entry routed through the scaffold ``load_anchor``."""

    anchor: Anchor
    value_id: str
    entry_index: int
    display: str | None
    limit_type: str | None
    cl: str | None
    experiment: str | None
    expected_limit_value: float | None
    expected_limit_units: str | None
    expected_limit_cl: str | None
    arxiv_id: str | None
    doi: str | None
    note: str | None

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def value_fb(self) -> float:
        return self.anchor.value

    @property
    def budget(self) -> float:
        return self.value_fb

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
class CR011Anchor:
    """YAML-loaded CR011 longitudinal-VBS limit bundle."""

    active_limit: CR011ValueAnchor
    cms_2020: CR011ValueAnchor
    all_limits: tuple[CR011ValueAnchor, ...]
    summary: str | None
    budget_policy: str

    @property
    def value(self) -> float:
        return self.active_limit.value_fb

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
            f"{process_id}: CR011 field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: CR011 field {field_name!r}={value!r} is not finite"
        )
    return number


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, field_name=field_name)


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


def _expected_limit(entry: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = entry.get("expected_limit")
    if expected is None:
        return {}
    if not isinstance(expected, Mapping):
        raise AnchorError("CR011: expected_limit must be a mapping when present")
    return expected


def _load_value_anchor(value_id: str, *, process_id: str) -> CR011ValueAnchor:
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
            f"expected {block_key!r} for CR011 value_id {value_id!r}"
        )

    expected = _expected_limit(entry)
    return CR011ValueAnchor(
        anchor=scaffold_anchor,
        value_id=value_id,
        entry_index=index,
        display=_optional_str(entry.get("display")),
        limit_type=_optional_str(entry.get("limit_type")),
        cl=_optional_str(entry.get("cl")),
        experiment=_optional_str(entry.get("experiment")),
        expected_limit_value=_optional_float(
            expected.get("value"),
            process_id=process_id,
            field_name=f"{value_id}.expected_limit.value",
        ),
        expected_limit_units=_optional_str(expected.get("units")),
        expected_limit_cl=_optional_str(expected.get("cl")),
        arxiv_id=_optional_str(entry.get("arxiv_id")),
        doi=_optional_str(entry.get("doi")),
        note=_optional_str(entry.get("note")),
    )


def _load_cr011_anchor(process_id: str) -> CR011Anchor:
    limits = tuple(
        _load_value_anchor(value_id, process_id=process_id)
        for value_id in _KNOWN_VALUE_IDS
    )
    by_id = {limit.value_id: limit for limit in limits}
    anchor = CR011Anchor(
        active_limit=by_id[_ACTIVE_ATLAS_2025],
        cms_2020=by_id[_CMS_2020],
        all_limits=limits,
        summary=_optional_str(_parent(process_id).get("edition")),
        budget_policy=_BUDGET_POLICY,
    )
    if anchor.budget <= 0.0:
        raise AnchorError(f"{process_id}: CR011 active budget must be positive")
    return anchor


def _limit_from_anchor(
    anchor: CR011ValueAnchor,
    *,
    process_id: str,
) -> VBSFiducialLimit:
    return VBSFiducialLimit(
        process_id=process_id,
        value_fb=float(anchor.value_fb),
        cl=anchor.cl,
        experiment=anchor.experiment,
        source=anchor.source,
        source_url=anchor.source_url,
        diagnostics={
            "value_id": anchor.value_id,
            "entry_index": anchor.entry_index,
            "display": anchor.display,
            "limit_type": anchor.limit_type,
            "snapshot_path": anchor.snapshot_path,
            "arxiv_id": anchor.arxiv_id,
            "doi": anchor.doi,
            "yaml_units": anchor.units,
            "yaml_observable": anchor.anchor.observable,
        },
    )


@register
class Constraint:
    """Catalogued non-vetoing longitudinal-VBS fiducial-limit record."""

    process_id = "CR011"
    severity = Severity.INFO
    observable = "sigma_fid(pp -> jj W_L^+/- W_L^+/-)"

    def __init__(self) -> None:
        self.anchor = _load_cr011_anchor(self.process_id)
        self.limit = _limit_from_anchor(
            self.anchor.active_limit,
            process_id=self.process_id,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        try:
            prediction = human_prediction_from_raw(point.raw)
            comparison = compare_vbs_fiducial_limit(self.limit, prediction)
        except (TypeError, ValueError, AttributeError) as exc:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "invalid human-supplied CR011 fiducial prediction; "
                    f"non-vetoing INFO record kept without recast: {exc}"
                ),
                diagnostics={
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                    "human_prediction_raw_key": HUMAN_SUPPLIED_SIGMA_RAW_KEY,
                    "budget_policy": self.anchor.budget_policy,
                    "needs_human_physics": (
                        VBS_LONGITUDINAL_SM_EFT_GAP_V1,
                        VBS_LONGITUDINAL_RS_MATCHING_GAP_V1,
                    ),
                    "needs_human_physics_sm_eft": VBS_LONGITUDINAL_SM_EFT_GAP_V1,
                    "needs_human_physics_rs_matching": (
                        VBS_LONGITUDINAL_RS_MATCHING_GAP_V1
                    ),
                },
            )

        diagnostics = dict(comparison.diagnostics)
        diagnostics.update(
            {
                "active_value_id": self.anchor.active_limit.value_id,
                "active_experiment": self.anchor.active_limit.experiment,
                "active_display": self.anchor.active_limit.display,
                "active_expected_limit_fb": (
                    self.anchor.active_limit.expected_limit_value
                ),
                "cms_2020_limit_fb": float(self.anchor.cms_2020.value_fb),
                "all_limits_fb": {
                    limit.value_id: float(limit.value_fb)
                    for limit in self.anchor.all_limits
                },
                "budget_policy": self.anchor.budget_policy,
                "severity_policy": "INFO/non-vetoing even when advisory passes=False",
                "human_prediction_raw_key": HUMAN_SUPPLIED_SIGMA_RAW_KEY,
                "needs_human_physics": (
                    VBS_LONGITUDINAL_SM_EFT_GAP_V1,
                    VBS_LONGITUDINAL_RS_MATCHING_GAP_V1,
                ),
                "needs_human_physics_sm_eft": VBS_LONGITUDINAL_SM_EFT_GAP_V1,
                "needs_human_physics_rs_matching": (
                    VBS_LONGITUDINAL_RS_MATCHING_GAP_V1
                ),
            }
        )

        if comparison.predicted_fiducial_sigma_fb is None:
            notes = (
                "CR011 records the ATLAS longitudinal same-sign WW fiducial "
                "upper limit only. No SM/EFT likelihood or RS strong-EWSB "
                "matching is available on ParameterPoint, so no recast is "
                "performed; result is INFO/non-vetoing."
            )
        else:
            notes = (
                "CR011 compares an explicitly human-supplied fiducial "
                "longitudinal-WW prediction with the ATLAS 2025 upper limit. "
                "This is still INFO/non-vetoing and is not an RS-derived "
                "mass/EFT recast."
            )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(comparison.passes),
            predicted=(
                None
                if comparison.predicted_fiducial_sigma_fb is None
                else float(comparison.predicted_fiducial_sigma_fb)
            ),
            experimental=float(comparison.experimental_limit_fb),
            ratio=(
                None
                if comparison.ratio_to_budget is None
                else float(comparison.ratio_to_budget)
            ),
            budget=float(comparison.budget_fb),
            notes=notes,
            diagnostics=diagnostics,
        )
