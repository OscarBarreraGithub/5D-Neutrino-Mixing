"""K008 - ``K_L -> pi0 e+ e-``.

Physics
-------
The total rate is conventionally decomposed as

    BR = [C_mix |a_S|^2 +/- C_int |a_S| A
          + C_dir A^2 + C_CPC] x N,

with ``A = Im(lambda_t) / 1e-4`` in the SM shorthand.  The live K008 path uses
the ISU/BMS Wilson-level form instead: direct CP depends separately on
``y7V`` and ``y7A``, while the direct/indirect interference depends only on the
vector ``y7V`` amplitude.

The indirect-CP, interference, and CPC pieces are parametrized from the
catalogued Isidori-Smith-Unterdorfer coefficients, but their sign/CPC treatment
is diagnostic only and marked ``NEEDS-HUMAN-PHYSICS``.

Severity
--------
HARD.  The veto compares the direct-CP short-distance branching fraction to
the observed PDG/KTeV 90% CL upper limit from ``K008.yaml``.  This is the
constrained NP part requested for the RS scan; total-rate envelopes are
reported in diagnostics because the ChPT interference sign is not fixed here.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K008.yaml`` is the source of truth for the
PDG limit, KTeV/NA48 supporting numbers, and the ISU coefficient bundle.  This
module adapts the flat and list-shaped YAML entries into the scaffold
``load_anchor`` path rather than hardcoding value-bearing anchors.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_kaon_dilepton import (
    KLongPi0EEChPTInputs,
    RARE_KAON_PI0EE_INTERFERENCE_LIMITATION_V1,
    RARE_KAON_PI0EE_PARAMETRIZATION_CITATION,
    RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1,
    RARE_KAON_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1,
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    klong_pi0ee_y7_direct_cp_from_rs_semileptonic_wilsons,
    klong_pi0ee_y7_direct_cp_sm,
    rare_kaon_dilepton_default_sm_inputs,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_PDG_LIMIT_BLOCK = "pdg_or_equivalent"
_VALUES_SECTION = "pdg_or_equivalent.values"
_SCAFFOLD_UNCERTAINTY_KEY = "__k008_uncertainty_is_parsed_below__"
_EXPECTED_BRANCHING_UNITS = "branching fraction"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/kaon/K008.yaml pdg_or_equivalent "
    "PDG/KTeV 90% CL upper limit"
)

_KTEV_STANDALONE_LIMIT_OBS = "KTeV 1999-2000 standalone BR(K_L -> pi0 e+ e-) limit"
_KTEV_SM_ORDER_OBS = "KTeV theory-context Standard-Model branching-ratio order"
_KTEV_TOTAL_CPV_RANGE_OBS = "KTeV theory-context total CP-violating branching-ratio range"
_KTEV_DIRECT_CPV_RANGE_OBS = "KTeV theory-context direct-CP branching-ratio estimate"
_CPC_COMPONENT_OBS = "PDG-listed CP-conserving part inferred from K_L -> pi0 gamma gamma"
_KS_PARTIAL_OBS = "PDG / NA48 BR(K_S -> pi0 e+ e-) for m_ee > 0.165 GeV"
_KS_FULL_OBS = "NA48 full-region extrapolated BR(K_S -> pi0 e+ e-)"
_FORMULA_NORMALIZATION_OBS = "Isidori-Smith-Unterdorfer formula normalization"
_C_MIX_OBS = "Isidori-Smith-Unterdorfer electron C_mix coefficient"
_C_INT_OBS = "Isidori-Smith-Unterdorfer electron C_int coefficient"
_C_DIR_OBS = "Isidori-Smith-Unterdorfer electron C_dir coefficient"
_C_CPC_OBS = "Isidori-Smith-Unterdorfer electron C_CPC coefficient"
_A_S_OBS = "Isidori-Smith-Unterdorfer |a_S| input"
_SM_IM_LAMBDA_OBS = "Isidori-Smith-Unterdorfer SM Im(lambda_t)/1e-4 input"

_UPPER_LIMIT_RE = re.compile(r"^\s*<\s*(?P<value>[0-9.eE+-]+)\s*$")
_RANGE_RE = re.compile(
    r"^\s*(?P<low>[0-9.eE+-]+)\s+to\s+(?P<high>[0-9.eE+-]+)\s*$"
)


@dataclass(frozen=True)
class CatalogValueAnchor:
    """Value-bearing K008 YAML entry routed through ``load_anchor``."""

    section_key: str
    entry_index: int
    block_key: str
    observable: str | None
    source: str | None
    source_key: str | None
    year: int | None
    value: float
    uncertainty: float | None
    is_upper_limit: bool
    range_low: float | None
    range_high: float | None
    confidence_level: str | None
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    value_raw: str | None


@dataclass(frozen=True)
class K008BudgetBand:
    """K008 upper-limit budget for the HARD direct-CP veto."""

    source: str
    construction: str
    hard_veto_budget: float
    confidence_level: str | None
    sm_direct_cp_branching_fraction: float
    sm_constructive_total_branching_fraction: float
    sm_destructive_total_branching_fraction: float
    ktev_direct_cpv_low: float
    ktev_direct_cpv_high: float
    ktev_total_cpv_low: float
    ktev_total_cpv_high: float
    sm_subtracted: bool


@dataclass(frozen=True)
class K008Anchor:
    """Typed K008 anchor: limit, supporting values, and ChPT coefficients."""

    experimental_limit: CatalogValueAnchor
    ktev_standalone_limit: CatalogValueAnchor
    ktev_sm_order: CatalogValueAnchor
    ktev_total_cpv_range: CatalogValueAnchor
    ktev_direct_cpv_range: CatalogValueAnchor
    cpc_component: CatalogValueAnchor
    ks_partial_rate: CatalogValueAnchor
    ks_full_rate: CatalogValueAnchor
    formula_normalization: CatalogValueAnchor
    c_mix: CatalogValueAnchor
    c_int: CatalogValueAnchor
    c_dir: CatalogValueAnchor
    c_cpc: CatalogValueAnchor
    a_s_abs: CatalogValueAnchor
    sm_im_lambda_t_over_1e_minus4: CatalogValueAnchor
    budget_band: K008BudgetBand | None = None

    @property
    def value(self) -> float:
        return self.experimental_limit.value

    @property
    def uncertainty(self) -> float | None:
        return None

    @property
    def budget(self) -> float:
        if self.budget_band is None:
            return self.experimental_limit.value
        return self.budget_band.hard_veto_budget

    def chpt_inputs(self) -> KLongPi0EEChPTInputs:
        return KLongPi0EEChPTInputs(
            normalization=float(self.formula_normalization.value),
            c_mix=float(self.c_mix.value),
            c_int=float(self.c_int.value),
            c_dir=float(self.c_dir.value),
            c_cpc=float(self.c_cpc.value),
            a_s_abs=float(self.a_s_abs.value),
            sm_im_lambda_t_over_1e_minus4=float(
                self.sm_im_lambda_t_over_1e_minus4.value
            ),
            citation=RARE_KAON_PI0EE_PARAMETRIZATION_CITATION,
        )


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: K008 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: K008 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    try:
        return _required_float(value, process_id=process_id, field_name=field_name)
    except AnchorError:
        return None


def _parse_value(
    value: Any,
    *,
    process_id: str,
    field_name: str,
    expect_upper_limit: bool = False,
    expect_range: bool = False,
) -> tuple[float, bool, str | None, float | None, float | None]:
    raw = _optional_str(value)
    if isinstance(value, str):
        limit_match = _UPPER_LIMIT_RE.match(value)
        if limit_match is not None:
            number = _required_float(
                limit_match.group("value"),
                process_id=process_id,
                field_name=field_name,
            )
            return number, True, raw, None, None
        range_match = _RANGE_RE.match(value)
        if range_match is not None:
            low = _required_float(
                range_match.group("low"),
                process_id=process_id,
                field_name=f"{field_name}.low",
            )
            high = _required_float(
                range_match.group("high"),
                process_id=process_id,
                field_name=f"{field_name}.high",
            )
            if high <= low:
                raise AnchorError(f"{process_id}: K008 range {field_name} is invalid")
            return 0.5 * (low + high), False, raw, low, high
    if expect_upper_limit:
        raise AnchorError(
            f"{process_id}: expected an upper-limit value in {field_name}, got {value!r}"
        )
    if expect_range:
        raise AnchorError(
            f"{process_id}: expected a range value in {field_name}, got {value!r}"
        )
    return (
        _required_float(value, process_id=process_id, field_name=field_name),
        False,
        raw,
        None,
        None,
    )


def _load_virtual_anchor(
    process_id: str,
    *,
    block_key: str,
    entry: Mapping[str, Any],
    expect_upper_limit: bool = False,
    expect_range: bool = False,
) -> tuple[Anchor, bool, str | None, float | None, float | None]:
    value, is_upper_limit, raw, range_low, range_high = _parse_value(
        entry.get("value"),
        process_id=process_id,
        field_name=f"{block_key}.value",
        expect_upper_limit=expect_upper_limit,
        expect_range=expect_range,
    )
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

    return scaffold_anchor, is_upper_limit, raw, range_low, range_high


def _catalog_anchor_from_entry(
    *,
    scaffold_anchor: Anchor,
    entry: Mapping[str, Any],
    section_key: str,
    entry_index: int,
    is_upper_limit: bool,
    value_raw: str | None,
    range_low: float | None,
    range_high: float | None,
    process_id: str,
) -> CatalogValueAnchor:
    return CatalogValueAnchor(
        section_key=section_key,
        entry_index=entry_index,
        block_key=scaffold_anchor.block_key,
        observable=_optional_str(scaffold_anchor.observable),
        source=_optional_str(scaffold_anchor.source),
        source_key=_optional_str(entry.get("source_key")),
        year=scaffold_anchor.year,
        value=float(scaffold_anchor.value),
        uncertainty=_optional_float(
            entry.get("uncertainty"),
            process_id=process_id,
            field_name=f"{scaffold_anchor.block_key}.uncertainty",
        ),
        is_upper_limit=is_upper_limit,
        range_low=range_low,
        range_high=range_high,
        confidence_level=_optional_str(entry.get("confidence_level")),
        units=scaffold_anchor.units,
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        value_raw=value_raw,
    )


def _load_flat_limit_anchor(process_id: str) -> CatalogValueAnchor:
    data = load_full_yaml(process_id, family=_FAMILY)
    entry = data.get("pdg_or_equivalent")
    if not isinstance(entry, Mapping):
        raise AnchorError(f"{process_id}: expected flat pdg_or_equivalent mapping")
    scaffold_anchor, is_upper_limit, raw, range_low, range_high = _load_virtual_anchor(
        process_id,
        block_key=_PDG_LIMIT_BLOCK,
        entry=entry,
        expect_upper_limit=True,
    )
    if scaffold_anchor.units != _EXPECTED_BRANCHING_UNITS:
        raise AnchorError(
            f"{process_id}: expected branching-fraction units for PDG limit, "
            f"got {scaffold_anchor.units!r}"
        )
    return _catalog_anchor_from_entry(
        scaffold_anchor=scaffold_anchor,
        entry=entry,
        section_key="pdg_or_equivalent",
        entry_index=-1,
        is_upper_limit=is_upper_limit,
        value_raw=raw,
        range_low=range_low,
        range_high=range_high,
        process_id=process_id,
    )


def _values_entries(data: Mapping[str, Any], *, process_id: str) -> Sequence[Mapping[str, Any]]:
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


def _value_entry_for_observable(
    process_id: str,
    observable: str,
) -> tuple[int, Mapping[str, Any]]:
    entries = _values_entries(load_full_yaml(process_id, family=_FAMILY), process_id=process_id)
    for index, entry in enumerate(entries):
        if entry.get("observable") == observable:
            return index, entry
    raise AnchorError(f"{process_id}: no pdg_or_equivalent.values entry {observable!r}")


def _load_value_anchor(
    process_id: str,
    observable: str,
    *,
    expect_upper_limit: bool = False,
    expect_range: bool = False,
    expected_units: str | None = None,
) -> CatalogValueAnchor:
    index, entry = _value_entry_for_observable(process_id, observable)
    block_key = f"pdg_or_equivalent.values[{index}]"
    scaffold_anchor, is_upper_limit, raw, range_low, range_high = _load_virtual_anchor(
        process_id,
        block_key=block_key,
        entry=entry,
        expect_upper_limit=expect_upper_limit,
        expect_range=expect_range,
    )
    if scaffold_anchor.observable != observable:
        raise AnchorError(
            f"{process_id}: selected {block_key} observable "
            f"{scaffold_anchor.observable!r}, expected {observable!r}"
        )
    if expected_units is not None and scaffold_anchor.units != expected_units:
        raise AnchorError(
            f"{process_id}: expected units {expected_units!r} for {observable!r}, "
            f"got {scaffold_anchor.units!r}"
        )
    return _catalog_anchor_from_entry(
        scaffold_anchor=scaffold_anchor,
        entry=entry,
        section_key=_VALUES_SECTION,
        entry_index=index,
        is_upper_limit=is_upper_limit,
        value_raw=raw,
        range_low=range_low,
        range_high=range_high,
        process_id=process_id,
    )


def _build_budget_band(
    anchor: K008Anchor,
    *,
    sm_direct_cp: float,
    sm_constructive_total: float,
    sm_destructive_total: float,
) -> K008BudgetBand:
    limit = float(anchor.experimental_limit.value)
    if not anchor.experimental_limit.is_upper_limit or limit <= 0.0:
        raise AnchorError("K008: PDG limit must be a positive upper limit")
    direct_low = anchor.ktev_direct_cpv_range.range_low
    direct_high = anchor.ktev_direct_cpv_range.range_high
    total_low = anchor.ktev_total_cpv_range.range_low
    total_high = anchor.ktev_total_cpv_range.range_high
    if None in (direct_low, direct_high, total_low, total_high):
        raise AnchorError("K008: KTeV direct/total CPV ranges are required")
    if not (float(direct_low) <= sm_direct_cp <= float(direct_high)):
        raise AnchorError(
            f"K008: SM direct-CP value {sm_direct_cp:.6e} is outside "
            f"the KTeV direct-CP range [{direct_low:.6e}, {direct_high:.6e}]"
        )
    if not (float(total_low) <= sm_constructive_total <= float(total_high)):
        raise AnchorError(
            f"K008: SM constructive total {sm_constructive_total:.6e} is outside "
            f"the KTeV CPV range [{total_low:.6e}, {total_high:.6e}]"
        )
    return K008BudgetBand(
        source=_BUDGET_SOURCE,
        construction=(
            "Observed 90% CL upper limit on BR(K_L -> pi0 e+ e-); HARD ratio "
            "uses the direct-CP short-distance prediction divided by this limit. "
            "No SM subtraction is applied."
        ),
        hard_veto_budget=limit,
        confidence_level=anchor.experimental_limit.confidence_level,
        sm_direct_cp_branching_fraction=float(sm_direct_cp),
        sm_constructive_total_branching_fraction=float(sm_constructive_total),
        sm_destructive_total_branching_fraction=float(sm_destructive_total),
        ktev_direct_cpv_low=float(direct_low),
        ktev_direct_cpv_high=float(direct_high),
        ktev_total_cpv_low=float(total_low),
        ktev_total_cpv_high=float(total_high),
        sm_subtracted=False,
    )


def _load_k008_anchor(
    process_id: str,
    *,
    sm_direct_cp: float | None = None,
    sm_constructive_total: float | None = None,
    sm_destructive_total: float | None = None,
) -> K008Anchor:
    anchor = K008Anchor(
        experimental_limit=_load_flat_limit_anchor(process_id),
        ktev_standalone_limit=_load_value_anchor(
            process_id,
            _KTEV_STANDALONE_LIMIT_OBS,
            expect_upper_limit=True,
            expected_units=_EXPECTED_BRANCHING_UNITS,
        ),
        ktev_sm_order=_load_value_anchor(
            process_id,
            _KTEV_SM_ORDER_OBS,
            expected_units=_EXPECTED_BRANCHING_UNITS,
        ),
        ktev_total_cpv_range=_load_value_anchor(
            process_id,
            _KTEV_TOTAL_CPV_RANGE_OBS,
            expect_range=True,
            expected_units=_EXPECTED_BRANCHING_UNITS,
        ),
        ktev_direct_cpv_range=_load_value_anchor(
            process_id,
            _KTEV_DIRECT_CPV_RANGE_OBS,
            expect_range=True,
            expected_units=_EXPECTED_BRANCHING_UNITS,
        ),
        cpc_component=_load_value_anchor(
            process_id,
            _CPC_COMPONENT_OBS,
            expected_units=_EXPECTED_BRANCHING_UNITS,
        ),
        ks_partial_rate=_load_value_anchor(
            process_id,
            _KS_PARTIAL_OBS,
            expected_units=_EXPECTED_BRANCHING_UNITS,
        ),
        ks_full_rate=_load_value_anchor(
            process_id,
            _KS_FULL_OBS,
            expected_units=_EXPECTED_BRANCHING_UNITS,
        ),
        formula_normalization=_load_value_anchor(
            process_id,
            _FORMULA_NORMALIZATION_OBS,
            expected_units="branching fraction normalization",
        ),
        c_mix=_load_value_anchor(process_id, _C_MIX_OBS),
        c_int=_load_value_anchor(process_id, _C_INT_OBS),
        c_dir=_load_value_anchor(process_id, _C_DIR_OBS),
        c_cpc=_load_value_anchor(process_id, _C_CPC_OBS),
        a_s_abs=_load_value_anchor(
            process_id,
            _A_S_OBS,
            expected_units="dimensionless",
        ),
        sm_im_lambda_t_over_1e_minus4=_load_value_anchor(
            process_id,
            _SM_IM_LAMBDA_OBS,
            expected_units="dimensionless ratio",
        ),
    )
    if sm_direct_cp is None or sm_constructive_total is None or sm_destructive_total is None:
        return anchor
    budget = _build_budget_band(
        anchor,
        sm_direct_cp=sm_direct_cp,
        sm_constructive_total=sm_constructive_total,
        sm_destructive_total=sm_destructive_total,
    )
    return K008Anchor(**{**anchor.__dict__, "budget_band": budget})


@register
class Constraint:
    """Catalogued ``K_L -> pi0 e+ e-`` direct-CP constraint (K008)."""

    process_id = "K008"
    severity = Severity.HARD
    observable = "BR(K_L -> pi0 e+ e-)_direct CP"

    def __init__(self) -> None:
        preliminary_anchor = _load_k008_anchor(self.process_id)
        self.chpt_inputs = preliminary_anchor.chpt_inputs()
        self.sm_inputs = rare_kaon_dilepton_default_sm_inputs()
        self.sm_result = klong_pi0ee_y7_direct_cp_sm(
            self.chpt_inputs,
            inputs=self.sm_inputs,
        )
        self.anchor = _load_k008_anchor(
            self.process_id,
            sm_direct_cp=float(self.sm_result.direct_cp_branching_fraction),
            sm_constructive_total=float(
                self.sm_result.constructive_total_branching_fraction
            ),
            sm_destructive_total=float(
                self.sm_result.destructive_total_branching_fraction
            ),
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        rs_wilsons = point.get_extra(_REQUIRED_EXTRA)
        if rs_wilsons is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_result.direct_cp_branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; K_L -> pi0 e+ e- "
                    "direct-CP constraint was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_quark_mass_basis_couplings_present": (
                        point.get_extra("quark_mass_basis_couplings") is not None
                    ),
                    "budget_source": self.anchor.budget_band.source,
                    "semileptonic_qcd_running_applied": False,
                    "semileptonic_qcd_running_multiplicative_factor": 1.0,
                    "semileptonic_qcd_running_effect_fraction": 0.0,
                    "semileptonic_qcd_running_diagnostic": (
                        RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1
                    ),
                    "sm_direct_cp_branching_fraction": float(
                        self.sm_result.direct_cp_branching_fraction
                    ),
                    "sm_constructive_total_branching_fraction": float(
                        self.sm_result.constructive_total_branching_fraction
                        + self.anchor.cpc_component.value
                        - self.sm_result.cpc_branching_fraction
                    ),
                    "sm_destructive_total_branching_fraction": float(
                        self.sm_result.destructive_total_branching_fraction
                        + self.anchor.cpc_component.value
                        - self.sm_result.cpc_branching_fraction
                    ),
                    "cpc_branching_fraction": float(self.anchor.cpc_component.value),
                    "cpc_branching_fraction_source": (
                        "K008.yaml PDG-listed CP-conserving part inferred from "
                        "K_L -> pi0 gamma gamma"
                    ),
                    "needs_human_physics": (
                        RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1
                        + " "
                        + RARE_KAON_PI0EE_INTERFERENCE_LIMITATION_V1
                    ),
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = klong_pi0ee_y7_direct_cp_from_rs_semileptonic_wilsons(
                rs_wilsons,
                chpt_inputs=self.chpt_inputs,
                lepton="e",
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_result.direct_cp_branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_semileptonic_wilsons for "
                    "K_L -> pi0 e+ e- direct CP"
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "budget_source": self.anchor.budget_band.source,
                },
            )
        predicted = float(result.direct_cp_branching_fraction)
        budget = float(self.anchor.budget)
        ratio = predicted / budget if budget > 0.0 else float("inf")
        passes = bool(ratio <= 1.0)

        diagnostics = dict(result.diagnostics)
        cpc_yaml = float(self.anchor.cpc_component.value)
        cpc_core = float(result.cpc_branching_fraction)
        constructive_with_yaml_cpc = (
            float(result.constructive_total_branching_fraction) - cpc_core + cpc_yaml
        )
        destructive_with_yaml_cpc = (
            float(result.destructive_total_branching_fraction) - cpc_core + cpc_yaml
        )
        sm_constructive_with_yaml_cpc = (
            float(result.sm_constructive_total_branching_fraction)
            - cpc_core
            + cpc_yaml
        )
        sm_destructive_with_yaml_cpc = (
            float(result.sm_destructive_total_branching_fraction)
            - cpc_core
            + cpc_yaml
        )
        diagnostics.update(
            {
                "experimental_limit": float(self.anchor.experimental_limit.value),
                "experimental_limit_block": self.anchor.experimental_limit.block_key,
                "experimental_confidence_level": (
                    self.anchor.experimental_limit.confidence_level
                ),
                "ktev_standalone_limit": float(
                    self.anchor.ktev_standalone_limit.value
                ),
                "budget": float(budget),
                "budget_source": self.anchor.budget_band.source,
                "budget_construction": self.anchor.budget_band.construction,
                "budget_confidence_level": self.anchor.budget_band.confidence_level,
                "budget_sm_subtracted": self.anchor.budget_band.sm_subtracted,
                "parametrization_citation": RARE_KAON_PI0EE_PARAMETRIZATION_CITATION,
                "rs_matching_assumption": RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
                "needs_human_physics": (
                    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1
                    + " "
                    + RARE_KAON_PI0EE_INTERFERENCE_LIMITATION_V1
                ),
                "rs_semileptonic_vector_matching_status": (
                    RARE_KAON_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1
                ),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "formula_normalization": float(
                    self.anchor.formula_normalization.value
                ),
                "c_mix": float(self.anchor.c_mix.value),
                "c_int": float(self.anchor.c_int.value),
                "c_dir": float(self.anchor.c_dir.value),
                "c_cpc": float(self.anchor.c_cpc.value),
                "a_s_abs": float(self.anchor.a_s_abs.value),
                "yaml_sm_im_lambda_t_over_1e_minus4": float(
                    self.anchor.sm_im_lambda_t_over_1e_minus4.value
                ),
                "direct_cp_amplitude_ratio": float(
                    result.direct_cp_amplitude_ratio
                ),
                "sm_direct_cp_amplitude_ratio": float(
                    result.sm_direct_cp_amplitude_ratio
                ),
                "direct_cp_branching_fraction": float(
                    result.direct_cp_branching_fraction
                ),
                "indirect_cp_branching_fraction": float(
                    result.indirect_cp_branching_fraction
                ),
                "interference_abs_branching_fraction": float(
                    result.interference_abs_branching_fraction
                ),
                "cpc_branching_fraction": cpc_yaml,
                "cpc_branching_fraction_source": (
                    "K008.yaml PDG-listed CP-conserving part inferred from "
                    "K_L -> pi0 gamma gamma"
                ),
                "cpc_isu_coefficient_branching_fraction": cpc_core,
                "constructive_total_branching_fraction": float(
                    constructive_with_yaml_cpc
                ),
                "destructive_total_branching_fraction": float(
                    destructive_with_yaml_cpc
                ),
                "sm_direct_cp_branching_fraction": float(
                    result.sm_direct_cp_branching_fraction
                ),
                "sm_constructive_total_branching_fraction": float(
                    sm_constructive_with_yaml_cpc
                ),
                "sm_destructive_total_branching_fraction": float(
                    sm_destructive_with_yaml_cpc
                ),
                "ktev_sm_order_branching_fraction": float(
                    self.anchor.ktev_sm_order.value
                ),
                "ktev_direct_cpv_range_low": float(
                    self.anchor.budget_band.ktev_direct_cpv_low
                ),
                "ktev_direct_cpv_range_high": float(
                    self.anchor.budget_band.ktev_direct_cpv_high
                ),
                "ktev_total_cpv_range_low": float(
                    self.anchor.budget_band.ktev_total_cpv_low
                ),
                "ktev_total_cpv_range_high": float(
                    self.anchor.budget_band.ktev_total_cpv_high
                ),
                "cpc_component_yaml_branching_fraction": float(
                    self.anchor.cpc_component.value
                ),
                "ks_partial_branching_fraction": float(
                    self.anchor.ks_partial_rate.value
                ),
                "ks_full_branching_fraction": float(self.anchor.ks_full_rate.value),
                "lambda_wolfenstein": float(result.lambda_wolfenstein),
                "lambda_t": complex(result.lambda_t),
                "y_t": float(result.y_t),
                "y_eff": complex(result.y_eff),
                "y_np_total": complex(result.y_np_total),
                "np_shift_branching_fraction": float(
                    result.direct_cp_np_shift_branching_fraction
                ),
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=predicted,
            sm_prediction=float(result.sm_direct_cp_branching_fraction),
            experimental=float(self.anchor.value),
            ratio=float(ratio),
            budget=budget,
            notes=(
                "BR(K_L -> pi0 e+ e-)_direct CP uses the "
                "Isidori-Smith-Unterdorfer/Buras-Mescia-Smith y7V/y7A "
                "Wilson structure: direct CP is vector plus axial, while "
                "the interference diagnostic is vector-only. The HARD "
                "budget is the K008.yaml PDG/KTeV 90% CL total-rate limit; "
                "Phase-3a RS C9/C10 Wilsons enter additively through y7V/y7A; "
                "indirect CP, interference, and CPC totals are diagnostics."
            ),
            diagnostics=diagnostics,
        )
