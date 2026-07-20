"""K009 - ``K_L -> pi0 mu+ mu-``.

Physics
-------
The total rate is conventionally decomposed as

    BR = [C_mix |a_S|^2 +/- C_int |a_S| A
          + C_dir A^2 + C_CPC] x N,

with ``A = Im(lambda_t) / 1e-4`` in the SM shorthand.  K009 reuses the K008
``y7V/y7A`` direct-CP machinery for the muon mode: direct CP depends on vector
and axial amplitudes, while the direct/indirect interference diagnostic uses
the vector amplitude only.  The muon phase-space coefficients are loaded from
``K009.yaml`` and converted to the y7 basis in the adapter.

Severity
--------
HARD.  The veto compares the direct-CP short-distance branching fraction to
the observed PDG/KTeV 90% CL upper limit from ``K009.yaml``.  Total-rate
envelopes are reported in diagnostics because the ChPT interference sign and
CPC treatment need human physics review.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K009.yaml`` is the source of truth for the
PDG limit, KTeV/NA48 supporting numbers, and the Isidori-Smith-Unterdorfer
muon coefficient bundle.  Value-bearing entries are routed through the scaffold
``load_anchor`` path, including the list-shaped YAML entries.

NEEDS-HUMAN-PHYSICS
-------------------
The short-distance RS y7V/y7A slots are filled from Phase-3a
``rs_semileptonic_wilsons.s_to_d_ll`` C9/C10/C9p/C10p.  A complete K009
prediction still needs the long-distance interference/CPC treatment and
non-light-Z neutral-current pieces that are not part of Phase 3a.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_kaon_dilepton import (
    RARE_KAON_PI0EE_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1,
    RARE_KAON_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1,
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    KLongPi0EEChPTInputs,
    rare_kaon_dilepton_default_sm_inputs,
)
from flavor_catalog_constraints.physics_adapters.rare_kaon_dilepton_muon import (
    RARE_KAON_PI0MUMU_PARAMETRIZATION_CITATION,
    RARE_KAON_PI0MUMU_RS_MATCHING_ASSUMPTION_V1,
    klong_pi0mumu_y7_direct_cp_from_rs_semileptonic_wilsons,
    klong_pi0mumu_y7_direct_cp_sm,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_PDG_LIMIT_BLOCK = "pdg_or_equivalent"
_VALUES_SECTION = "pdg_or_equivalent.values"
_SCAFFOLD_UNCERTAINTY_KEY = "__k009_uncertainty_is_parsed_below__"
_EXPECTED_BRANCHING_UNITS = "branching fraction"
_SM_VALIDATION_SIGMA_TOLERANCE = 2.0
_BUDGET_SOURCE = (
    "flavor_catalog/processes/kaon/K009.yaml pdg_or_equivalent "
    "PDG/KTeV 90% CL upper limit"
)

_KTEV_LIMIT_OBS = "KTeV 1997-data BR(K_L -> pi0 mu+ mu-) limit"
_KTEV_OBSERVED_OBS = "KTeV 1997-data observed candidates for K_L -> pi0 mu+ mu-"
_KTEV_BACKGROUND_OBS = "KTeV 1997-data expected background for K_L -> pi0 mu+ mu-"
_KS_RATE_OBS = "BR(K_S -> pi0 mu+ mu-) supporting indirect-CP input"
_NA48_OBSERVED_OBS = "NA48 observed K_S -> pi0 mu+ mu- candidates"
_NA48_BACKGROUND_OBS = "NA48 expected background for K_S -> pi0 mu+ mu-"
_FORMULA_NORMALIZATION_OBS = (
    "Isidori-Smith-Unterdorfer rate decomposition normalization for "
    "K_L -> pi0 mu+ mu-"
)
_C_MIX_OBS = "C_mix^mu coefficient in K_L -> pi0 mu+ mu- decomposition"
_C_INT_OBS = "C_int^mu coefficient in K_L -> pi0 mu+ mu- decomposition"
_C_DIR_OBS = "C_dir^mu coefficient in K_L -> pi0 mu+ mu- decomposition"
_C_CPC_OBS = "C_CPC^mu coefficient in K_L -> pi0 mu+ mu- decomposition"
_A_S_OBS = "|a_S| input for K_L -> pi0 mu+ mu- decomposition"
_SM_CONSTRUCTIVE_OBS = (
    "Constructive-interference SM expectation for BR(K_L -> pi0 mu+ mu-)"
)
_SM_IM_LAMBDA_THEORY_KEY = "Im_lambda_t_over_1e_minus4_SM"
_SM_DESTRUCTIVE_THEORY_KEY = "SM_destructive_interference"

_UPPER_LIMIT_RE = re.compile(r"^\s*<\s*(?P<value>[0-9.eE+-]+)\s*$")
_RANGE_RE = re.compile(
    r"^\s*(?P<low>[0-9.eE+-]+)\s+to\s+(?P<high>[0-9.eE+-]+)\s*$"
)
_VALUE_WITH_UNCERT_TIMES_RE = re.compile(
    r"^\s*\(?\s*(?P<value>[0-9.eE+-]+)\s*\+/-\s*"
    r"(?P<uncertainty>[0-9.eE+-]+)\s*\)?\s*x\s*10\^(?P<exponent>[+-]?[0-9]+)\s*$"
)
_VALUE_WITH_UNCERT_RE = re.compile(
    r"^\s*\(?\s*(?P<value>[0-9.eE+-]+)\s*\+/-\s*"
    r"(?P<uncertainty>[0-9.eE+-]+)\s*\)?\s*$"
)
_NUMERIC_PREFIX_RE = re.compile(
    r"^\s*(?P<value>[+-]?[0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)(?:\s+.*)?$"
)


@dataclass(frozen=True)
class CatalogValueAnchor:
    """Value-bearing K009 YAML entry routed through ``load_anchor``."""

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
class K009BudgetBand:
    """K009 upper-limit budget and SM validation diagnostics."""

    source: str
    construction: str
    hard_veto_budget: float
    confidence_level: str | None
    sm_direct_cp_branching_fraction: float
    sm_constructive_total_branching_fraction: float
    sm_destructive_total_branching_fraction: float
    sm_constructive_anchor_value: float
    sm_constructive_anchor_uncertainty: float
    sm_destructive_anchor_value: float
    sm_destructive_anchor_uncertainty: float
    sm_constructive_anchor_pull: float
    sm_destructive_anchor_pull: float
    sm_subtracted: bool


@dataclass(frozen=True)
class K009Anchor:
    """Typed K009 anchor: limit, supporting values, and ChPT coefficients."""

    experimental_limit: CatalogValueAnchor
    ktev_limit: CatalogValueAnchor
    ktev_observed_candidates: CatalogValueAnchor
    ktev_expected_background: CatalogValueAnchor
    ks_supporting_rate: CatalogValueAnchor
    na48_observed_candidates: CatalogValueAnchor
    na48_expected_background: CatalogValueAnchor
    formula_normalization: CatalogValueAnchor
    c_mix: CatalogValueAnchor
    c_int: CatalogValueAnchor
    c_dir: CatalogValueAnchor
    c_cpc: CatalogValueAnchor
    a_s_abs: CatalogValueAnchor
    sm_im_lambda_t_over_1e_minus4: CatalogValueAnchor
    sm_constructive_total: CatalogValueAnchor
    sm_destructive_total: CatalogValueAnchor
    budget_band: K009BudgetBand | None = None

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
            citation=RARE_KAON_PI0MUMU_PARAMETRIZATION_CITATION,
        )


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: K009 anchor field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: K009 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _value_with_uncertainty_times_10(value: str) -> tuple[float, float] | None:
    match = _VALUE_WITH_UNCERT_TIMES_RE.match(value)
    if match is None:
        plain_match = _VALUE_WITH_UNCERT_RE.match(value)
        if plain_match is None:
            return None
        return float(plain_match.group("value")), float(
            plain_match.group("uncertainty")
        )
    scale = 10.0 ** int(match.group("exponent"))
    return float(match.group("value")) * scale, float(match.group("uncertainty")) * scale


def _numeric_prefix(value: str) -> float | None:
    stripped = value.strip()
    if stripped.startswith("+") and ("-" in stripped[1:] or "/" in stripped):
        return None
    match = _NUMERIC_PREFIX_RE.match(value)
    if match is None:
        return None
    return float(match.group("value"))


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        parsed = _value_with_uncertainty_times_10(value)
        if parsed is not None:
            return float(parsed[1])
        prefix = _numeric_prefix(value)
        if prefix is not None:
            return float(prefix)
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
) -> tuple[float, bool, str | None, float | None, float | None, float | None]:
    raw = _optional_str(value)
    parsed_uncertainty = None
    if isinstance(value, str):
        limit_match = _UPPER_LIMIT_RE.match(value)
        if limit_match is not None:
            number = _required_float(
                limit_match.group("value"),
                process_id=process_id,
                field_name=field_name,
            )
            return number, True, raw, None, None, parsed_uncertainty
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
                raise AnchorError(f"{process_id}: K009 range {field_name} is invalid")
            return 0.5 * (low + high), False, raw, low, high, parsed_uncertainty
        value_with_uncertainty = _value_with_uncertainty_times_10(value)
        if value_with_uncertainty is not None:
            number, parsed_uncertainty = value_with_uncertainty
            return float(number), False, raw, None, None, float(parsed_uncertainty)
        prefix = _numeric_prefix(value)
        if prefix is not None:
            return float(prefix), False, raw, None, None, parsed_uncertainty
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
        parsed_uncertainty,
    )


def _load_virtual_anchor(
    process_id: str,
    *,
    block_key: str,
    entry: Mapping[str, Any],
    expect_upper_limit: bool = False,
    expect_range: bool = False,
) -> tuple[Anchor, bool, str | None, float | None, float | None, float | None]:
    value, is_upper_limit, raw, range_low, range_high, parsed_uncertainty = _parse_value(
        entry.get("value"),
        process_id=process_id,
        field_name=f"{block_key}.value",
        expect_upper_limit=expect_upper_limit,
        expect_range=expect_range,
    )
    virtual_entry = dict(entry)
    virtual_entry["value"] = value
    if parsed_uncertainty is not None and virtual_entry.get("uncertainty") is None:
        virtual_entry["uncertainty"] = parsed_uncertainty
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

    return scaffold_anchor, is_upper_limit, raw, range_low, range_high, parsed_uncertainty


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
    parsed_uncertainty: float | None,
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
        uncertainty=(
            parsed_uncertainty
            if parsed_uncertainty is not None
            else _optional_float(
                entry.get("uncertainty"),
                process_id=process_id,
                field_name=f"{scaffold_anchor.block_key}.uncertainty",
            )
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
    (
        scaffold_anchor,
        is_upper_limit,
        raw,
        range_low,
        range_high,
        parsed_uncertainty,
    ) = _load_virtual_anchor(
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
        parsed_uncertainty=parsed_uncertainty,
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
    (
        scaffold_anchor,
        is_upper_limit,
        raw,
        range_low,
        range_high,
        parsed_uncertainty,
    ) = _load_virtual_anchor(
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
        parsed_uncertainty=parsed_uncertainty,
        process_id=process_id,
    )


def _load_theory_muon_anchor(process_id: str, key: str) -> CatalogValueAnchor:
    data = load_full_yaml(process_id, family=_FAMILY)
    theory = data.get("theory_decomposition")
    if not isinstance(theory, Mapping):
        raise AnchorError(f"{process_id}: expected theory_decomposition mapping")
    coefficients = theory.get("muon_coefficients")
    if not isinstance(coefficients, Mapping):
        raise AnchorError(f"{process_id}: expected theory_decomposition.muon_coefficients")
    if key not in coefficients:
        raise AnchorError(
            f"{process_id}: no theory_decomposition.muon_coefficients[{key!r}]"
        )
    units = "dimensionless ratio" if key == _SM_IM_LAMBDA_THEORY_KEY else "branching fraction"
    entry = {
        "observable": key,
        "year": theory.get("year"),
        "value": coefficients[key],
        "uncertainty": None,
        "units": units,
        "source": theory.get("source"),
        "source_key": "IsidoriSmithUnterdorfer2004:KLPi0MuMu",
        "source_url": theory.get("source_url"),
        "snapshot_path": theory.get("snapshot_path"),
    }
    block_key = f"theory_decomposition.muon_coefficients.{key}"
    (
        scaffold_anchor,
        is_upper_limit,
        raw,
        range_low,
        range_high,
        parsed_uncertainty,
    ) = _load_virtual_anchor(
        process_id,
        block_key=block_key,
        entry=entry,
    )
    return _catalog_anchor_from_entry(
        scaffold_anchor=scaffold_anchor,
        entry=entry,
        section_key="theory_decomposition.muon_coefficients",
        entry_index=-1,
        is_upper_limit=is_upper_limit,
        value_raw=raw,
        range_low=range_low,
        range_high=range_high,
        parsed_uncertainty=parsed_uncertainty,
        process_id=process_id,
    )


def _require_uncertainty(anchor: CatalogValueAnchor, *, process_id: str) -> float:
    if anchor.uncertainty is None or anchor.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: {anchor.block_key} requires a positive uncertainty"
        )
    return float(anchor.uncertainty)


def _validate_sm_total(
    *,
    process_id: str,
    label: str,
    predicted: float,
    anchor: CatalogValueAnchor,
) -> float:
    sigma = _require_uncertainty(anchor, process_id=process_id)
    pull = (float(predicted) - float(anchor.value)) / sigma
    if abs(pull) > _SM_VALIDATION_SIGMA_TOLERANCE:
        raise AnchorError(
            f"{process_id}: SM {label} total {predicted:.6e} differs from YAML "
            f"anchor {anchor.value:.6e} by {pull:.3g} sigma"
        )
    return float(pull)


def _build_budget_band(
    anchor: K009Anchor,
    *,
    sm_direct_cp: float,
    sm_constructive_total: float,
    sm_destructive_total: float,
) -> K009BudgetBand:
    limit = float(anchor.experimental_limit.value)
    if not anchor.experimental_limit.is_upper_limit or limit <= 0.0:
        raise AnchorError("K009: PDG limit must be a positive upper limit")
    constructive_pull = _validate_sm_total(
        process_id="K009",
        label="constructive-interference",
        predicted=sm_constructive_total,
        anchor=anchor.sm_constructive_total,
    )
    destructive_pull = _validate_sm_total(
        process_id="K009",
        label="destructive-interference",
        predicted=sm_destructive_total,
        anchor=anchor.sm_destructive_total,
    )
    return K009BudgetBand(
        source=_BUDGET_SOURCE,
        construction=(
            "Observed 90% CL upper limit on BR(K_L -> pi0 mu+ mu-); HARD ratio "
            "uses the direct-CP short-distance prediction divided by this limit. "
            "No SM subtraction is applied."
        ),
        hard_veto_budget=limit,
        confidence_level=anchor.experimental_limit.confidence_level,
        sm_direct_cp_branching_fraction=float(sm_direct_cp),
        sm_constructive_total_branching_fraction=float(sm_constructive_total),
        sm_destructive_total_branching_fraction=float(sm_destructive_total),
        sm_constructive_anchor_value=float(anchor.sm_constructive_total.value),
        sm_constructive_anchor_uncertainty=float(anchor.sm_constructive_total.uncertainty),
        sm_destructive_anchor_value=float(anchor.sm_destructive_total.value),
        sm_destructive_anchor_uncertainty=float(anchor.sm_destructive_total.uncertainty),
        sm_constructive_anchor_pull=constructive_pull,
        sm_destructive_anchor_pull=destructive_pull,
        sm_subtracted=False,
    )


def _load_k009_anchor(
    process_id: str,
    *,
    sm_direct_cp: float | None = None,
    sm_constructive_total: float | None = None,
    sm_destructive_total: float | None = None,
) -> K009Anchor:
    anchor = K009Anchor(
        experimental_limit=_load_flat_limit_anchor(process_id),
        ktev_limit=_load_value_anchor(
            process_id,
            _KTEV_LIMIT_OBS,
            expect_upper_limit=True,
            expected_units=_EXPECTED_BRANCHING_UNITS,
        ),
        ktev_observed_candidates=_load_value_anchor(
            process_id,
            _KTEV_OBSERVED_OBS,
            expected_units="events",
        ),
        ktev_expected_background=_load_value_anchor(
            process_id,
            _KTEV_BACKGROUND_OBS,
            expected_units="events",
        ),
        ks_supporting_rate=_load_value_anchor(
            process_id,
            _KS_RATE_OBS,
            expected_units=_EXPECTED_BRANCHING_UNITS,
        ),
        na48_observed_candidates=_load_value_anchor(
            process_id,
            _NA48_OBSERVED_OBS,
            expected_units="events",
        ),
        na48_expected_background=_load_value_anchor(
            process_id,
            _NA48_BACKGROUND_OBS,
            expected_units="events",
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
        sm_im_lambda_t_over_1e_minus4=_load_theory_muon_anchor(
            process_id,
            _SM_IM_LAMBDA_THEORY_KEY,
        ),
        sm_constructive_total=_load_value_anchor(
            process_id,
            _SM_CONSTRUCTIVE_OBS,
            expected_units=_EXPECTED_BRANCHING_UNITS,
        ),
        sm_destructive_total=_load_theory_muon_anchor(
            process_id,
            _SM_DESTRUCTIVE_THEORY_KEY,
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
    return K009Anchor(**{**anchor.__dict__, "budget_band": budget})


@register
class Constraint:
    """Catalogued ``K_L -> pi0 mu+ mu-`` direct-CP constraint (K009)."""

    process_id = "K009"
    severity = Severity.HARD
    observable = "BR(K_L -> pi0 mu+ mu-)_direct CP"

    def __init__(self) -> None:
        preliminary_anchor = _load_k009_anchor(self.process_id)
        self.chpt_inputs = preliminary_anchor.chpt_inputs()
        self.sm_inputs = rare_kaon_dilepton_default_sm_inputs()
        self.sm_result = klong_pi0mumu_y7_direct_cp_sm(
            self.chpt_inputs,
            inputs=self.sm_inputs,
        )
        self.anchor = _load_k009_anchor(
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
                    f"extra {_REQUIRED_EXTRA!r} absent; K_L -> pi0 mu+ mu- "
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
                    ),
                    "sm_destructive_total_branching_fraction": float(
                        self.sm_result.destructive_total_branching_fraction
                    ),
                    "cpc_branching_fraction": float(
                        self.sm_result.cpc_branching_fraction
                    ),
                    "needs_human_physics": (
                        RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1
                        + " "
                        + RARE_KAON_PI0MUMU_RS_MATCHING_ASSUMPTION_V1
                    ),
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = klong_pi0mumu_y7_direct_cp_from_rs_semileptonic_wilsons(
                rs_wilsons,
                chpt_inputs=self.chpt_inputs,
                lepton="mu",
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
                    "K_L -> pi0 mu+ mu- direct CP"
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
        diagnostics.update(
            {
                "experimental_limit": float(self.anchor.experimental_limit.value),
                "experimental_limit_block": self.anchor.experimental_limit.block_key,
                "experimental_confidence_level": (
                    self.anchor.experimental_limit.confidence_level
                ),
                "ktev_limit": float(self.anchor.ktev_limit.value),
                "ktev_observed_candidates": float(
                    self.anchor.ktev_observed_candidates.value
                ),
                "ktev_expected_background": float(
                    self.anchor.ktev_expected_background.value
                ),
                "ks_supporting_branching_fraction": float(
                    self.anchor.ks_supporting_rate.value
                ),
                "na48_observed_candidates": float(
                    self.anchor.na48_observed_candidates.value
                ),
                "na48_expected_background": float(
                    self.anchor.na48_expected_background.value
                ),
                "budget": float(budget),
                "budget_source": self.anchor.budget_band.source,
                "budget_construction": self.anchor.budget_band.construction,
                "budget_confidence_level": self.anchor.budget_band.confidence_level,
                "budget_sm_subtracted": self.anchor.budget_band.sm_subtracted,
                "parametrization_citation": RARE_KAON_PI0MUMU_PARAMETRIZATION_CITATION,
                "rs_matching_assumption": RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
                "needs_human_physics": (
                    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1
                    + " "
                    + RARE_KAON_PI0MUMU_RS_MATCHING_ASSUMPTION_V1
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
                "direct_cp_vector_amplitude_ratio": float(
                    result.direct_cp_vector_amplitude_ratio
                ),
                "direct_cp_axial_amplitude_ratio": float(
                    result.direct_cp_axial_amplitude_ratio
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
                "cpc_branching_fraction": float(result.cpc_branching_fraction),
                "constructive_total_branching_fraction": float(
                    result.constructive_total_branching_fraction
                ),
                "destructive_total_branching_fraction": float(
                    result.destructive_total_branching_fraction
                ),
                "sm_direct_cp_branching_fraction": float(
                    result.sm_direct_cp_branching_fraction
                ),
                "sm_constructive_total_branching_fraction": float(
                    result.sm_constructive_total_branching_fraction
                ),
                "sm_destructive_total_branching_fraction": float(
                    result.sm_destructive_total_branching_fraction
                ),
                "yaml_sm_constructive_total_branching_fraction": float(
                    self.anchor.sm_constructive_total.value
                ),
                "yaml_sm_constructive_total_uncertainty": float(
                    self.anchor.sm_constructive_total.uncertainty
                ),
                "yaml_sm_destructive_total_branching_fraction": float(
                    self.anchor.sm_destructive_total.value
                ),
                "yaml_sm_destructive_total_uncertainty": float(
                    self.anchor.sm_destructive_total.uncertainty
                ),
                "sm_constructive_anchor_pull": float(
                    self.anchor.budget_band.sm_constructive_anchor_pull
                ),
                "sm_destructive_anchor_pull": float(
                    self.anchor.budget_band.sm_destructive_anchor_pull
                ),
                "lambda_wolfenstein": float(result.lambda_wolfenstein),
                "lambda_t": complex(result.lambda_t),
                "lambda_y7v_eff": complex(result.lambda_y7v_eff),
                "lambda_y7a_eff": complex(result.lambda_y7a_eff),
                "lambda_y7v_np_proxy": complex(result.lambda_y7v_np_proxy),
                "lambda_y7a_np_proxy": complex(result.lambda_y7a_np_proxy),
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
                "BR(K_L -> pi0 mu+ mu-)_direct CP reuses the K008 "
                "Isidori-Smith-Unterdorfer/Buras-Mescia-Smith y7V/y7A "
                "machinery with K009 muon phase-space coefficients. The HARD "
                "budget is the K009.yaml PDG/KTeV 90% CL total-rate limit; "
                "Phase-3a RS C9/C10 Wilsons enter additively through y7V/y7A; "
                "indirect CP, interference, and CPC totals are diagnostics."
            ),
            diagnostics=diagnostics,
        )
