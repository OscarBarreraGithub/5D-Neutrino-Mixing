"""T019 - lepton-flavor-violating ``h -> e tau`` branching fraction.

Physics
-------
The Standard Model contribution is zero for catalog purposes.  This constraint
therefore applies a pure-NP HARD upper bound to the charge-summed branching
fraction ``BR(h -> e+- tau-+)``.  The rate uses the shared Higgs LFV
effective-Yukawa convention

    Gamma(h -> e tau) = m_h (|Y_e_tau|^2 + |Y_tau_e|^2) / (8 pi),
    BR = Gamma / Gamma_h^total.

RS matching status
------------------
The production path reads the Phase-6b ``rs_higgs_yukawas`` extra and routes
its mass-basis Higgs-Yukawa matrix through the shared Higgs LFV adapter.  With
the repo v1 diagonal charged-lepton fit this gives exactly zero off-diagonal
tree Higgs-LFV entries; nonzero generic RS Higgs-LFV requires a non-diagonal
charged-lepton Yukawa/rotation structure supplied by a future lepton fit.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/T019.yaml`` is the source of truth for
the 95% CL upper limit and provenance.  The numeric branching-fraction limit is
parsed from the YAML entry and routed through the scaffold ``load_anchor`` path;
no experimental number is hardcoded here.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.higgs_lfv import (
    higgs_lfv_branching_fraction_from_yukawas,
    higgs_lfv_branching_fraction_with_proxy,
    higgs_lfv_default_inputs,
    higgs_lfv_effective_yukawa_limit,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_REQUIRED_EXTRA = "rs_higgs_yukawas"
_PDG_ATLAS_LIMIT_VALUE_ID = "PDG2025:T019:atlas_run2"
_PDG_CMS_LIMIT_VALUE_ID = "PDG2025:T019:cms_run2"
_ATLAS_LIMIT_VALUE_ID = "ATLAS2023:T019:etau_limit"
_ATLAS_DATASET_VALUE_ID = "ATLAS2023:T019:dataset"
_CMS_LIMIT_VALUE_ID = "CMS2021:T019:etau_limit"
_CMS_DATASET_VALUE_ID = "CMS2021:T019:dataset"
_ATLAS_2019_LIMIT_VALUE_ID = "ATLAS2019:T019:etau_limit"
_ATLAS_2019_DATASET_VALUE_ID = "ATLAS2019:T019:dataset"
_HKZ_EFT_VALUE_ID = "HarnikKoppZupan2012:T019:LFV_Higgs_EFT_allowance"
_NUMBER_RE = re.compile(
    r"^\s*<?\s*(?P<number>[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"
    r"(?:[eE][+-]?[0-9]+)?)\s*(?P<percent>%?)\s*$"
)
_UNEVALUATED_REASON = (
    "no off-diagonal Higgs-Yukawa e-tau prediction available "
    "(rs_higgs_yukawas not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"
_NEEDS_HUMAN_PHYSICS = (
    "Phase-6b v1 evaluates the minimal-RS tree/ZMA Higgs-LFV matrix supplied "
    "on rs_higgs_yukawas. Generic nonzero RS Higgs-LFV still requires a "
    "future non-diagonal charged-lepton fit; the current diagonal builder is "
    "rigorous tree-level zero."
)


@dataclass(frozen=True)
class HiggsLFVLimitEntry:
    """Typed view over one value entry in the T019 YAML ``values`` list."""

    value_id: str
    observable: str | None
    year: int | None
    limit: float
    expected_limit: float | None
    confidence_level: str | None
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    value_summary: str | None


@dataclass(frozen=True)
class T019Anchor:
    """Typed T019 anchor bundle: ATLAS budget plus CMS/PDG context."""

    experimental: Anchor
    atlas_limit: HiggsLFVLimitEntry
    pdg_atlas_limit: HiggsLFVLimitEntry
    pdg_cms_limit: HiggsLFVLimitEntry
    cms_limit: HiggsLFVLimitEntry
    atlas_2019_limit: HiggsLFVLimitEntry
    atlas_dataset: Mapping[str, Any]
    cms_dataset: Mapping[str, Any]
    atlas_2019_dataset: Mapping[str, Any]
    harnik_kopp_zupan_eft_allowance: Mapping[str, Any]
    effective_yukawa_limit: float

    @property
    def value(self) -> float:
        return self.experimental.value

    @property
    def budget(self) -> float:
        return self.experimental.value

    @property
    def source_url(self) -> str | None:
        return self.experimental.source_url


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _parse_limit_value(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> float:
    if isinstance(value, (int, float)):
        number = float(value)
    elif isinstance(value, str):
        match = _NUMBER_RE.match(value)
        if match is None:
            raise AnchorError(
                f"{process_id}: T019 field {field_name!r}={value!r} is not "
                "a numeric upper-limit string"
            )
        number = float(match.group("number"))
        if match.group("percent"):
            number /= 100.0
    else:
        raise AnchorError(
            f"{process_id}: T019 field {field_name!r}={value!r} is not numeric"
        )
    if not math.isfinite(number) or number <= 0.0:
        raise AnchorError(f"{process_id}: T019 field {field_name!r} must be positive")
    return number


def _pdg_values(process_id: str) -> Sequence[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg = data.get("pdg_or_equivalent")
    if not isinstance(pdg, Mapping):
        raise AnchorError(f"{process_id}: expected pdg_or_equivalent mapping")
    values = pdg.get("values")
    if not isinstance(values, list) or not values:
        raise AnchorError(f"{process_id}: expected non-empty pdg_or_equivalent.values")
    for index, entry in enumerate(values):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent.values[{index}] is not a mapping"
            )
    return values


def _find_value_entry(
    process_id: str,
    value_id: str,
) -> tuple[int, Mapping[str, Any]]:
    for index, entry in enumerate(_pdg_values(process_id)):
        if entry.get("value_id") == value_id:
            return index, entry
    present = [str(entry.get("value_id")) for entry in _pdg_values(process_id)]
    raise AnchorError(
        f"{process_id}: value_id {value_id!r} not found in T019 values "
        f"(present: {present})"
    )


def _load_scaffold_value_anchor(
    value_id: str,
    *,
    process_id: str,
) -> Anchor:
    index, entry = _find_value_entry(process_id, value_id)
    block_key = f"pdg_or_equivalent.values[{index}]"
    numeric_entry = dict(entry)
    numeric_entry["value"] = _parse_limit_value(
        entry.get("normalized_value", entry.get("value")),
        process_id=process_id,
        field_name=f"{value_id}.normalized_value",
    )
    if "snapshot_path" not in numeric_entry:
        pdg = load_full_yaml(process_id, family=_FAMILY)["pdg_or_equivalent"]
        if isinstance(pdg, Mapping) and "snapshot_path" in pdg:
            numeric_entry["snapshot_path"] = pdg["snapshot_path"]
    virtual_block = {block_key: numeric_entry}
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
            f"expected {block_key!r} for value_id {value_id!r}"
        )
    return scaffold_anchor


def _entry_view(
    process_id: str,
    value_id: str,
    *,
    value_summary: str | None,
) -> HiggsLFVLimitEntry:
    _, entry = _find_value_entry(process_id, value_id)
    expected_raw = entry.get("expected_normalized_value", entry.get("expected_value"))
    return HiggsLFVLimitEntry(
        value_id=value_id,
        observable=_optional_str(entry.get("observable")),
        year=_optional_int(entry.get("year")),
        limit=_parse_limit_value(
            entry.get("normalized_value", entry.get("value")),
            process_id=process_id,
            field_name=f"{value_id}.normalized_value",
        ),
        expected_limit=(
            None
            if expected_raw is None
            else _parse_limit_value(
                expected_raw,
                process_id=process_id,
                field_name=f"{value_id}.expected_value",
            )
        ),
        confidence_level=_optional_str(entry.get("confidence_level")),
        units=_optional_str(entry.get("units")),
        source_url=_optional_str(entry.get("source_url")),
        snapshot_path=_optional_str(entry.get("snapshot_path")),
        value_summary=value_summary,
    )


def _dataset_entry(process_id: str, value_id: str) -> Mapping[str, Any]:
    _, entry = _find_value_entry(process_id, value_id)
    return dict(entry)


def _load_t019_anchor(process_id: str) -> T019Anchor:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg = data.get("pdg_or_equivalent")
    if not isinstance(pdg, Mapping):
        raise AnchorError(f"{process_id}: expected pdg_or_equivalent mapping")
    value_summary = _optional_str(pdg.get("value_summary"))
    experimental = _load_scaffold_value_anchor(
        _ATLAS_LIMIT_VALUE_ID,
        process_id=process_id,
    )
    if experimental.units != "branching fraction":
        raise AnchorError(
            f"{process_id}: ATLAS limit units must be 'branching fraction', "
            f"got {experimental.units!r}"
        )
    atlas = _entry_view(process_id, _ATLAS_LIMIT_VALUE_ID, value_summary=value_summary)
    pdg_atlas = _entry_view(
        process_id,
        _PDG_ATLAS_LIMIT_VALUE_ID,
        value_summary=value_summary,
    )
    pdg_cms = _entry_view(
        process_id,
        _PDG_CMS_LIMIT_VALUE_ID,
        value_summary=value_summary,
    )
    cms = _entry_view(process_id, _CMS_LIMIT_VALUE_ID, value_summary=value_summary)
    atlas_2019 = _entry_view(
        process_id,
        _ATLAS_2019_LIMIT_VALUE_ID,
        value_summary=value_summary,
    )
    if not math.isclose(experimental.value, atlas.limit, rel_tol=0.0, abs_tol=0.0):
        raise AnchorError(f"{process_id}: scaffold anchor and ATLAS entry disagree")
    if not math.isclose(experimental.value, pdg_atlas.limit, rel_tol=0.0, abs_tol=0.0):
        raise AnchorError(f"{process_id}: ATLAS entry and PDG ATLAS entry disagree")
    if experimental.value > cms.limit:
        raise AnchorError(
            f"{process_id}: selected ATLAS budget is weaker than CMS limit"
        )
    return T019Anchor(
        experimental=experimental,
        atlas_limit=atlas,
        pdg_atlas_limit=pdg_atlas,
        pdg_cms_limit=pdg_cms,
        cms_limit=cms,
        atlas_2019_limit=atlas_2019,
        atlas_dataset=_dataset_entry(process_id, _ATLAS_DATASET_VALUE_ID),
        cms_dataset=_dataset_entry(process_id, _CMS_DATASET_VALUE_ID),
        atlas_2019_dataset=_dataset_entry(process_id, _ATLAS_2019_DATASET_VALUE_ID),
        harnik_kopp_zupan_eft_allowance=_dataset_entry(process_id, _HKZ_EFT_VALUE_ID),
        effective_yukawa_limit=higgs_lfv_effective_yukawa_limit(
            experimental.value,
            inputs=higgs_lfv_default_inputs(),
        ),
    )


@register
class Constraint:
    """Catalogued ``h -> e tau`` LFV pure-NP branching-fraction constraint."""

    process_id = "T019"
    severity = Severity.HARD
    observable = "BR(h -> e tau)"

    def __init__(self) -> None:
        self.anchor = _load_t019_anchor(self.process_id)
        self.sm_inputs = higgs_lfv_default_inputs()
        self.sm_result = higgs_lfv_branching_fraction_from_yukawas(
            initial_flavor="e",
            final_flavor="tau",
            yukawa_ij=0.0j,
            yukawa_ji=0.0j,
            br_limit=self.anchor.budget,
            inputs=self.sm_inputs,
        )

    def _unevaluated_result(self, diagnostics: Mapping[str, Any]) -> ConstraintResult:
        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=True,
            predicted=None,
            sm_prediction=0.0,
            experimental=float(self.anchor.value),
            ratio=None,
            budget=float(self.anchor.budget),
            notes=_UNEVALUATED_NOTES,
            diagnostics={
                "evaluated": False,
                "unevaluated_reason": _UNEVALUATED_REASON,
                "passes_semantics": (
                    "non-vetoing only; no BR(h -> e tau) NP prediction was evaluated"
                ),
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "budget_source": self.anchor.source_url,
                "experimental_block": self.anchor.experimental.block_key,
                "sm_branching_fraction": 0.0,
                **dict(diagnostics),
            },
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        higgs_input = point.get_extra(_REQUIRED_EXTRA)
        if higgs_input is None:
            return self._unevaluated_result({"missing_extra": _REQUIRED_EXTRA})

        try:
            result, proxy = higgs_lfv_branching_fraction_with_proxy(
                higgs_input,
                initial_flavor="e",
                final_flavor="tau",
                br_limit=self.anchor.budget,
                inputs=self.sm_inputs,
            )
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            return self._unevaluated_result(
                {
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception": str(exc),
                },
            )

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                "sm_prediction_policy": (
                    "Charged-LFV Higgs decay has zero SM rate for catalog "
                    "purposes; the HARD budget is applied to the pure-NP "
                    "branching fraction."
                ),
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "rs_matching_assumption": str(proxy.matching_assumption),
                "rs_higgs_yukawas_units": getattr(higgs_input, "units", None),
                "includes_fermion_kk_mixing": bool(
                    getattr(higgs_input, "includes_fermion_kk_mixing", False)
                ),
                "rs_higgs_yukawas_diagnostics": dict(
                    getattr(higgs_input, "diagnostics", {})
                ),
                "higgs_lfv_proxy": dict(proxy.diagnostics),
                "y_e_tau": complex(result.yukawa_ij),
                "y_tau_e": complex(result.yukawa_ji),
                "effective_yukawa_norm": float(result.yukawa_norm),
                "effective_yukawa_norm_squared": float(result.yukawa_norm_squared),
                "effective_yukawa_limit": float(self.anchor.effective_yukawa_limit),
                "atlas_limit_value_id": self.anchor.atlas_limit.value_id,
                "atlas_expected_limit": self.anchor.atlas_limit.expected_limit,
                "pdg_atlas_limit": float(self.anchor.pdg_atlas_limit.limit),
                "pdg_cms_limit": float(self.anchor.pdg_cms_limit.limit),
                "cms_limit": float(self.anchor.cms_limit.limit),
                "atlas_2019_limit": float(self.anchor.atlas_2019_limit.limit),
                "atlas_dataset": dict(self.anchor.atlas_dataset),
                "cms_dataset": dict(self.anchor.cms_dataset),
                "atlas_2019_dataset": dict(self.anchor.atlas_2019_dataset),
                "harnik_kopp_zupan_eft_allowance": dict(
                    self.anchor.harnik_kopp_zupan_eft_allowance
                ),
                "experimental_block": self.anchor.experimental.block_key,
                "budget_source": self.anchor.source_url,
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(result.passes),
            predicted=float(result.branching_fraction),
            sm_prediction=0.0,
            experimental=float(self.anchor.value),
            ratio=float(result.ratio_to_limit),
            budget=float(result.br_limit),
            notes=(
                "Pure-NP BR(h -> e tau) bound using the shared Higgs LFV "
                "effective-Yukawa width formula and the Phase-6b "
                "rs_higgs_yukawas tree matrix. The diagonal charged-lepton "
                "v1 builder gives rigorous zero off-diagonal entries; "
                "non-diagonal charged-lepton fits are deferred."
            ),
            diagnostics=diagnostics,
        )
