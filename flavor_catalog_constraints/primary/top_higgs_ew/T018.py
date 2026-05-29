"""T018 - lepton-flavor-violating ``h -> mu tau`` branching fraction.

Physics
-------
The Standard Model contribution is zero for catalog purposes.  This constraint
therefore applies a pure-NP HARD upper bound to the charge-summed branching
fraction ``BR(h -> mu+- tau-+)``.  The rate uses the shared Higgs LFV
effective-Yukawa convention

    Gamma(h -> mu tau) = m_h (|Y_mu_tau|^2 + |Y_tau_mu|^2) / (8 pi),
    BR = Gamma / Gamma_h^total.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A rigorous prediction needs the off-diagonal physical
charged-lepton Higgs-Yukawa matrix after RS Higgs localization, KK-fermion
mixing, and charged-lepton mass-basis rotation.  That matrix is not present on
``ParameterPoint``.  If a caller supplies ``lepton_mass_basis_couplings``, this
v1 implementation uses a documented effective-Yukawa proxy through the Higgs
LFV adapter: ``Y_mu_tau`` and ``Y_tau_mu`` are taken directly from the supplied
proxy or Higgs-Yukawa matrix.

Catalog sidecar
---------------
``flavor_catalog/processes/top_higgs_ew/T018.yaml`` is the source of truth for
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
    HIGGS_LFV_RS_PROXY_V1,
    higgs_lfv_branching_fraction_from_yukawas,
    higgs_lfv_branching_fraction_with_proxy,
    higgs_lfv_default_inputs,
    higgs_lfv_effective_yukawa_limit,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "top_higgs_ew"
_REQUIRED_EXTRA = "lepton_mass_basis_couplings"
_CMS_LIMIT_VALUE_ID = "CMS2021:T018:mutau_limit"
_PDG_CMS_LIMIT_VALUE_ID = "PDG2025:T018:cms_run2"
_PDG_ATLAS_LIMIT_VALUE_ID = "PDG2025:T018:atlas_run2"
_ATLAS_LIMIT_VALUE_ID = "ATLAS2023:T018:mutau_limit"
_CMS_DATASET_VALUE_ID = "CMS2021:T018:dataset"
_ATLAS_DATASET_VALUE_ID = "ATLAS2023:T018:dataset"
_CMS_2015_YUKAWA_VALUE_ID = "CMS2015:T018:yukawa_limit"
_NUMBER_RE = re.compile(
    r"^\s*<?\s*(?P<number>[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"
    r"(?:[eE][+-]?[0-9]+)?)\s*(?P<percent>%?)\s*$"
)
_UNEVALUATED_REASON = (
    "no off-diagonal Higgs-Yukawa mu-tau prediction available "
    "(charged-lepton Higgs-Yukawa FCNC not on ParameterPoint)"
)
_UNEVALUATED_NOTES = f"NOT EVALUATED - {_UNEVALUATED_REASON}"
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: off-diagonal charged-lepton Higgs Yukawas are not "
    "on ParameterPoint; this v1 uses a documented effective-Yukawa proxy "
    "through the shared Higgs LFV module."
)


@dataclass(frozen=True)
class HiggsLFVLimitEntry:
    """Typed view over one value entry in the T018 YAML ``values`` list."""

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
class T018Anchor:
    """Typed T018 anchor bundle: CMS budget plus ATLAS/PDG context."""

    experimental: Anchor
    cms_limit: HiggsLFVLimitEntry
    pdg_cms_limit: HiggsLFVLimitEntry
    pdg_atlas_limit: HiggsLFVLimitEntry
    atlas_limit: HiggsLFVLimitEntry
    cms_dataset: Mapping[str, Any]
    atlas_dataset: Mapping[str, Any]
    cms_2015_yukawa_limit: HiggsLFVLimitEntry
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
                f"{process_id}: T018 field {field_name!r}={value!r} is not "
                "a numeric upper-limit string"
            )
        number = float(match.group("number"))
        if match.group("percent"):
            number /= 100.0
    else:
        raise AnchorError(
            f"{process_id}: T018 field {field_name!r}={value!r} is not numeric"
        )
    if not math.isfinite(number) or number <= 0.0:
        raise AnchorError(f"{process_id}: T018 field {field_name!r} must be positive")
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
        f"{process_id}: value_id {value_id!r} not found in T018 values "
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


def _load_t018_anchor(process_id: str) -> T018Anchor:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg = data.get("pdg_or_equivalent")
    if not isinstance(pdg, Mapping):
        raise AnchorError(f"{process_id}: expected pdg_or_equivalent mapping")
    value_summary = _optional_str(pdg.get("value_summary"))
    experimental = _load_scaffold_value_anchor(
        _CMS_LIMIT_VALUE_ID,
        process_id=process_id,
    )
    if experimental.units != "branching fraction":
        raise AnchorError(
            f"{process_id}: CMS limit units must be 'branching fraction', "
            f"got {experimental.units!r}"
        )
    cms = _entry_view(process_id, _CMS_LIMIT_VALUE_ID, value_summary=value_summary)
    pdg_cms = _entry_view(
        process_id,
        _PDG_CMS_LIMIT_VALUE_ID,
        value_summary=value_summary,
    )
    pdg_atlas = _entry_view(
        process_id,
        _PDG_ATLAS_LIMIT_VALUE_ID,
        value_summary=value_summary,
    )
    atlas = _entry_view(process_id, _ATLAS_LIMIT_VALUE_ID, value_summary=value_summary)
    cms_2015_yukawa = _entry_view(
        process_id,
        _CMS_2015_YUKAWA_VALUE_ID,
        value_summary=value_summary,
    )
    if not math.isclose(experimental.value, cms.limit, rel_tol=0.0, abs_tol=0.0):
        raise AnchorError(f"{process_id}: scaffold anchor and CMS entry disagree")
    if experimental.value > atlas.limit:
        raise AnchorError(
            f"{process_id}: selected CMS budget is weaker than ATLAS limit"
        )
    return T018Anchor(
        experimental=experimental,
        cms_limit=cms,
        pdg_cms_limit=pdg_cms,
        pdg_atlas_limit=pdg_atlas,
        atlas_limit=atlas,
        cms_dataset=_dataset_entry(process_id, _CMS_DATASET_VALUE_ID),
        atlas_dataset=_dataset_entry(process_id, _ATLAS_DATASET_VALUE_ID),
        cms_2015_yukawa_limit=cms_2015_yukawa,
        effective_yukawa_limit=higgs_lfv_effective_yukawa_limit(
            experimental.value,
            inputs=higgs_lfv_default_inputs(),
        ),
    )


@register
class Constraint:
    """Catalogued ``h -> mu tau`` LFV pure-NP branching-fraction constraint."""

    process_id = "T018"
    severity = Severity.HARD
    observable = "BR(h -> mu tau)"

    def __init__(self) -> None:
        self.anchor = _load_t018_anchor(self.process_id)
        self.sm_inputs = higgs_lfv_default_inputs()
        self.sm_result = higgs_lfv_branching_fraction_from_yukawas(
            initial_flavor="mu",
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
                    "non-vetoing only; no BR(h -> mu tau) NP prediction was evaluated"
                ),
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "budget_source": self.anchor.source_url,
                "experimental_block": self.anchor.experimental.block_key,
                "sm_branching_fraction": 0.0,
                **dict(diagnostics),
            },
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        lepton_input = point.get_extra(_REQUIRED_EXTRA)
        if lepton_input is None:
            return self._unevaluated_result({"missing_extra": _REQUIRED_EXTRA})

        try:
            result, proxy = higgs_lfv_branching_fraction_with_proxy(
                lepton_input,
                initial_flavor="mu",
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
                "rs_matching_assumption": HIGGS_LFV_RS_PROXY_V1,
                "higgs_lfv_proxy": dict(proxy.diagnostics),
                "y_mu_tau": complex(result.yukawa_ij),
                "y_tau_mu": complex(result.yukawa_ji),
                "effective_yukawa_norm": float(result.yukawa_norm),
                "effective_yukawa_limit": float(self.anchor.effective_yukawa_limit),
                "cms_limit_value_id": self.anchor.cms_limit.value_id,
                "cms_expected_limit": self.anchor.cms_limit.expected_limit,
                "pdg_cms_limit": float(self.anchor.pdg_cms_limit.limit),
                "pdg_atlas_limit": float(self.anchor.pdg_atlas_limit.limit),
                "atlas_limit": float(self.anchor.atlas_limit.limit),
                "cms_2015_yukawa_limit": float(
                    self.anchor.cms_2015_yukawa_limit.limit
                ),
                "cms_dataset": dict(self.anchor.cms_dataset),
                "atlas_dataset": dict(self.anchor.atlas_dataset),
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
                "Pure-NP BR(h -> mu tau) bound using the shared Higgs LFV "
                "effective-Yukawa width formula. The off-diagonal Higgs "
                "Yukawas are documented proxy inputs and are flagged "
                "NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
