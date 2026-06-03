"""K012 - short-distance component of ``K_S -> mu+ mu-``.

Physics
-------
The observed ``K_S -> mu+ mu-`` total rate is long-distance dominated.  For the
RS scan, K012 constrains only the short-distance dimuon ell=0 piece by reusing
the K006 ``Y`` effective-input core with Phase-3a
``rs_semileptonic_wilsons.s_to_d_ll`` C10/C10p and applying the K_S imaginary
CP projection supplied by the K012 adapter,

    BR(K_S -> mu+mu-)_SD,ell=0
        = (tau_KS / tau_KL) kappa_mu [
            Im(lambda_t Y_t + Y_NP) / lambda^5
            - Im(lambda_c) P_c(Y) / lambda
        ]^2.

This is intentionally only the constrained short-distance piece, not a
complete total-rate prediction.

Severity
--------
HARD.  The short-distance component is compared with the K012.yaml PDG/LHCb 90% CL
upper limit.  The sidecar's SM entry is an approximate total-rate context
number, not a short-distance anchor; it is carried diagnostically and is not
used for a formula validation or subtraction.  The default SM short-distance
formula gives about 1.86e-13.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K012.yaml`` is the source of truth for the
PDG/LHCb upper limits, approximate SM total-rate context, and hadronic
uncertainty context.  Its ``pdg_or_equivalent`` block is list-shaped, so K012
adapts selected list entries into the scaffold ``load_anchor`` path rather than
hardcoding value-bearing anchors.

NEEDS-HUMAN-PHYSICS
-------------------
The semileptonic short-distance input is the Phase-3a light-Z Wilson path.
Complete K012 rigor still needs the K_S/K_L CP and time-dependent extraction
treatment, the long-distance two-photon amplitude, and non-light-Z or
muon-sector effects outside Phase 3a.
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
    RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION,
    RARE_KAON_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1,
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    rare_kaon_dilepton_default_sm_inputs,
)
from flavor_catalog_constraints.physics_adapters.rare_kaon_dilepton_kshort_mumu import (
    RARE_KAON_KSHORT_MUMU_PARAMETRIZATION_CITATION,
    RARE_KAON_KSHORT_MUMU_RS_MATCHING_ASSUMPTION_V1,
    kshort_mumu_lifetime_inputs_default,
    kshort_mumu_short_distance_from_rs_semileptonic_wilsons,
    kshort_mumu_short_distance_sm,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_VALUES_SECTION = "pdg_or_equivalent"
_SCAFFOLD_UNCERTAINTY_KEY = "__k012_uncertainty_is_parsed_below__"
_EXPECTED_BRANCHING_UNITS = "branching fraction"
_EXPECTED_RELATIVE_UNITS = "relative hadronic uncertainty"
_HEADLINE_LIMIT_VALUE_ID = "PDG2025:K012:headline_limit"
_LHCB_STANDALONE_VALUE_ID = "LHCb2020:K012:standalone_2016_2018_limit"
_LHCB_COMBINED_VALUE_ID = "LHCb2020:K012:combined_limit"
_RUN1_VALUE_ID = "LHCb2017:K012:run1_limit"
_SM_ESTIMATE_VALUE_ID = "Chobanova2018:K012:sm_estimate"
_HADRONIC_UNCERTAINTY_VALUE_ID = (
    "DeryGhoshGrossmanSchacht2021:K012:ell0_hadronic_uncertainty"
)
_BUDGET_SOURCE = (
    "flavor_catalog/processes/kaon/K012.yaml "
    "PDG2025:K012:headline_limit 90% CL upper limit"
)
_NUMBER_PATTERN = (
    r"(?P<number>[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"
    r"(?:[eE][+-]?[0-9]+)?)"
)
_UPPER_LIMIT_RE = re.compile(rf"^\s*<\s*{_NUMBER_PATTERN}(?:\s+.*)?$")
_APPROXIMATE_RE = re.compile(rf"^\s*approximately\s+{_NUMBER_PATTERN}\s*$")
_PERCENT_LIMIT_RE = re.compile(rf"^\s*<\s*{_NUMBER_PATTERN}\s*%\s*$")


@dataclass(frozen=True)
class CatalogValueAnchor:
    """Value-bearing K012 YAML list entry routed through ``load_anchor``."""

    section_key: str
    entry_index: int
    block_key: str
    value_id: str
    observable: str | None
    source: str | None
    year: int | None
    value: float
    uncertainty: float | None
    is_upper_limit: bool
    is_approximate: bool
    confidence_level: str | None
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    value_raw: str | None


@dataclass(frozen=True)
class K012BudgetBand:
    """K012 upper-limit budget for the HARD short-distance component."""

    source: str
    construction: str
    hard_veto_budget: float
    confidence_level: str | None
    sm_sd_formula_value: float
    sm_total_estimate_value: float
    sm_total_estimate_is_short_distance: bool
    limit_minus_sm_sd_formula: float
    sm_subtracted: bool
    long_distance_dominated: bool


@dataclass(frozen=True)
class K012Anchor:
    """Typed K012 anchor bundle."""

    headline_limit: CatalogValueAnchor
    lhcb_standalone_limit: CatalogValueAnchor
    lhcb_combined_limit: CatalogValueAnchor
    run1_limit: CatalogValueAnchor
    standard_model_total_estimate: CatalogValueAnchor
    hadronic_uncertainty_target: CatalogValueAnchor
    budget_band: K012BudgetBand

    @property
    def value(self) -> float:
        return self.headline_limit.value

    @property
    def uncertainty(self) -> float | None:
        return None

    @property
    def budget(self) -> float:
        return self.budget_band.hard_veto_budget

    @property
    def source_url(self) -> str | None:
        return self.headline_limit.source_url


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        raise AnchorError(f"{process_id}: K012 field {field_name!r} is not finite")
    return number


def _required_positive(value: float, *, process_id: str, field_name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise AnchorError(f"{process_id}: K012 field {field_name!r} must be positive")
    return number


def _entries(process_id: str) -> Sequence[Mapping[str, Any]]:
    data = load_full_yaml(process_id, family=_FAMILY)
    entries = data.get("pdg_or_equivalent")
    if not isinstance(entries, list) or not entries:
        raise AnchorError(f"{process_id}: expected non-empty pdg_or_equivalent list")
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent[{index}] is not a mapping"
            )
    return entries


def _find_entry(process_id: str, value_id: str) -> tuple[int, Mapping[str, Any]]:
    entries = _entries(process_id)
    for index, entry in enumerate(entries):
        if entry.get("value_id") == value_id:
            return index, entry
    present = [str(entry.get("value_id")) for entry in entries]
    raise AnchorError(
        f"{process_id}: value_id {value_id!r} not found in K012 "
        f"pdg_or_equivalent entries (present: {present})"
    )


def _parse_entry_value(
    value: Any,
    *,
    process_id: str,
    field_name: str,
    expect_upper_limit: bool = False,
    expect_approximate: bool = False,
    expect_percent_limit: bool = False,
) -> tuple[float, bool, bool, str | None]:
    raw = _optional_str(value)
    if isinstance(value, (int, float)):
        if expect_upper_limit or expect_approximate or expect_percent_limit:
            raise AnchorError(
                f"{process_id}: K012 field {field_name!r}={value!r} does not "
                "carry the expected limit/approximate marker"
            )
        number = _required_positive(
            float(value),
            process_id=process_id,
            field_name=field_name,
        )
        return number, False, False, raw
    if isinstance(value, str):
        if expect_percent_limit:
            match = _PERCENT_LIMIT_RE.match(value)
            if match is None:
                raise AnchorError(
                    f"{process_id}: expected a percent upper limit in "
                    f"{field_name}, got {value!r}"
                )
            number = _required_positive(
                float(match.group("number")) / 100.0,
                process_id=process_id,
                field_name=field_name,
            )
            return number, True, False, raw

        limit_match = _UPPER_LIMIT_RE.match(value)
        if limit_match is not None:
            if expect_approximate:
                raise AnchorError(
                    f"{process_id}: expected approximate value in {field_name}, "
                    f"got upper limit {value!r}"
                )
            number = _required_positive(
                float(limit_match.group("number")),
                process_id=process_id,
                field_name=field_name,
            )
            return number, True, False, raw

        approx_match = _APPROXIMATE_RE.match(value)
        if approx_match is not None:
            if expect_upper_limit:
                raise AnchorError(
                    f"{process_id}: expected upper-limit value in {field_name}, "
                    f"got approximate value {value!r}"
                )
            number = _required_positive(
                float(approx_match.group("number")),
                process_id=process_id,
                field_name=field_name,
            )
            return number, False, True, raw

    if expect_upper_limit:
        raise AnchorError(
            f"{process_id}: expected upper-limit value in {field_name}, "
            f"got {value!r}"
        )
    if expect_approximate:
        raise AnchorError(
            f"{process_id}: expected approximate value in {field_name}, "
            f"got {value!r}"
        )
    raise AnchorError(f"{process_id}: K012 field {field_name!r} is not numeric")


def _load_value_anchor(
    process_id: str,
    value_id: str,
    *,
    expect_upper_limit: bool = False,
    expect_approximate: bool = False,
    expect_percent_limit: bool = False,
    expected_units: str | None = None,
) -> CatalogValueAnchor:
    index, entry = _find_entry(process_id, value_id)
    block_key = f"pdg_or_equivalent[{index}]"
    value, is_upper_limit, is_approximate, raw = _parse_entry_value(
        entry.get("value"),
        process_id=process_id,
        field_name=f"{value_id}.value",
        expect_upper_limit=expect_upper_limit,
        expect_approximate=expect_approximate,
        expect_percent_limit=expect_percent_limit,
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

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for value_id {value_id!r}"
        )
    if scaffold_anchor.observable != entry.get("observable"):
        raise AnchorError(
            f"{process_id}: selected {block_key} observable "
            f"{scaffold_anchor.observable!r}, expected {entry.get('observable')!r}"
        )
    if expected_units is not None and scaffold_anchor.units != expected_units:
        raise AnchorError(
            f"{process_id}: expected units {expected_units!r} for {value_id}, "
            f"got {scaffold_anchor.units!r}"
        )

    return CatalogValueAnchor(
        section_key=_VALUES_SECTION,
        entry_index=index,
        block_key=scaffold_anchor.block_key,
        value_id=value_id,
        observable=_optional_str(scaffold_anchor.observable),
        source=_optional_str(scaffold_anchor.source),
        year=scaffold_anchor.year,
        value=float(scaffold_anchor.value),
        uncertainty=_optional_float(
            entry.get("uncertainty"),
            process_id=process_id,
            field_name=f"{block_key}.uncertainty",
        ),
        is_upper_limit=is_upper_limit,
        is_approximate=is_approximate,
        confidence_level=_optional_str(entry.get("confidence_level")),
        units=scaffold_anchor.units,
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        value_raw=raw,
    )


def _build_budget_band(
    *,
    headline_limit: CatalogValueAnchor,
    standard_model_total_estimate: CatalogValueAnchor,
    sm_formula_value: float,
) -> K012BudgetBand:
    budget = float(headline_limit.value)
    if not headline_limit.is_upper_limit or budget <= 0.0:
        raise AnchorError("K012: headline budget must be a positive upper limit")
    if not standard_model_total_estimate.is_approximate:
        raise AnchorError("K012: SM context entry must be marked approximate")
    return K012BudgetBand(
        source=_BUDGET_SOURCE,
        construction=(
            "Compare the K006-Wilson-derived K_S ell=0 Im-projected "
            "short-distance component directly with the K012.yaml PDG/LHCb "
            "total-rate upper limit. No SM or long-distance subtraction is "
            "applied."
        ),
        hard_veto_budget=budget,
        confidence_level=headline_limit.confidence_level,
        sm_sd_formula_value=float(sm_formula_value),
        sm_total_estimate_value=float(standard_model_total_estimate.value),
        sm_total_estimate_is_short_distance=False,
        limit_minus_sm_sd_formula=float(budget - sm_formula_value),
        sm_subtracted=False,
        long_distance_dominated=True,
    )


def _load_k012_anchor(process_id: str, *, sm_formula_value: float) -> K012Anchor:
    headline = _load_value_anchor(
        process_id,
        _HEADLINE_LIMIT_VALUE_ID,
        expect_upper_limit=True,
        expected_units=_EXPECTED_BRANCHING_UNITS,
    )
    standalone = _load_value_anchor(
        process_id,
        _LHCB_STANDALONE_VALUE_ID,
        expect_upper_limit=True,
        expected_units=_EXPECTED_BRANCHING_UNITS,
    )
    combined = _load_value_anchor(
        process_id,
        _LHCB_COMBINED_VALUE_ID,
        expect_upper_limit=True,
        expected_units=_EXPECTED_BRANCHING_UNITS,
    )
    run1 = _load_value_anchor(
        process_id,
        _RUN1_VALUE_ID,
        expect_upper_limit=True,
        expected_units=_EXPECTED_BRANCHING_UNITS,
    )
    sm_estimate = _load_value_anchor(
        process_id,
        _SM_ESTIMATE_VALUE_ID,
        expect_approximate=True,
        expected_units=_EXPECTED_BRANCHING_UNITS,
    )
    hadronic = _load_value_anchor(
        process_id,
        _HADRONIC_UNCERTAINTY_VALUE_ID,
        expect_percent_limit=True,
        expected_units=_EXPECTED_RELATIVE_UNITS,
    )
    return K012Anchor(
        headline_limit=headline,
        lhcb_standalone_limit=standalone,
        lhcb_combined_limit=combined,
        run1_limit=run1,
        standard_model_total_estimate=sm_estimate,
        hadronic_uncertainty_target=hadronic,
        budget_band=_build_budget_band(
            headline_limit=headline,
            standard_model_total_estimate=sm_estimate,
            sm_formula_value=sm_formula_value,
        ),
    )


@register
class Constraint:
    """Catalogued short-distance ``K_S -> mu+ mu-`` constraint (K012)."""

    process_id = "K012"
    severity = Severity.HARD
    observable = "BR(K_S -> mu+ mu-)_SD"

    def __init__(self) -> None:
        self.sm_inputs = rare_kaon_dilepton_default_sm_inputs()
        self.lifetime_inputs = kshort_mumu_lifetime_inputs_default()
        self.sm_result = kshort_mumu_short_distance_sm(
            self.sm_inputs,
            lifetime_inputs=self.lifetime_inputs,
        )
        self.anchor = _load_k012_anchor(
            self.process_id,
            sm_formula_value=float(self.sm_result.branching_fraction),
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        rs_wilsons = point.get_extra(_REQUIRED_EXTRA)
        if rs_wilsons is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_result.branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; K_S -> mu+ mu- "
                    "short-distance constraint was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_quark_mass_basis_couplings_present": (
                        point.get_extra("quark_mass_basis_couplings") is not None
                    ),
                    "pdg_headline_limit": float(self.anchor.headline_limit.value),
                    "lhcb_combined_limit": float(
                        self.anchor.lhcb_combined_limit.value
                    ),
                    "sm_context_total_branching_fraction": float(
                        self.anchor.standard_model_total_estimate.value
                    ),
                    "sm_context_is_short_distance": False,
                    "sm_formula_short_distance_branching_fraction": float(
                        self.sm_result.branching_fraction
                    ),
                    "kshort_ell0_projection_source": (
                        "Im[-lambda_c Y_c + lambda_t C10]^2"
                    ),
                    "uses_imaginary_ks_ell0_projection": True,
                    "sm_context_total_not_sd_anchor": True,
                    "budget_source": self.anchor.budget_band.source,
                    "long_distance_dominated_total_rate": True,
                    "needs_human_physics": (
                        RARE_KAON_KSHORT_MUMU_RS_MATCHING_ASSUMPTION_V1
                    ),
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = kshort_mumu_short_distance_from_rs_semileptonic_wilsons(
                rs_wilsons,
                lepton="mu",
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
                lifetime_inputs=self.lifetime_inputs,
            )
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_result.branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_semileptonic_wilsons for "
                    "K_S -> mu+ mu- short-distance"
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "budget_source": self.anchor.budget_band.source,
                    "uses_imaginary_ks_ell0_projection": True,
                },
            )
        predicted = float(result.branching_fraction)
        budget = float(self.anchor.budget)
        ratio = predicted / budget if budget > 0.0 else float("inf")
        passes = bool(ratio <= 1.0)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "pdg_headline_limit": float(self.anchor.headline_limit.value),
                "pdg_headline_limit_block": self.anchor.headline_limit.block_key,
                "lhcb_standalone_limit": float(
                    self.anchor.lhcb_standalone_limit.value
                ),
                "lhcb_combined_limit": float(self.anchor.lhcb_combined_limit.value),
                "run1_limit_90cl": float(self.anchor.run1_limit.value),
                "sm_context_total_branching_fraction": float(
                    self.anchor.standard_model_total_estimate.value
                ),
                "sm_context_value_raw": self.anchor.standard_model_total_estimate.value_raw,
                "sm_context_is_short_distance": False,
                "hadronic_uncertainty_target": float(
                    self.anchor.hadronic_uncertainty_target.value
                ),
                "sm_formula_short_distance_branching_fraction": float(
                    self.sm_result.branching_fraction
                ),
                "sm_formula_minus_sm_context_total": float(
                    self.sm_result.branching_fraction
                    - self.anchor.standard_model_total_estimate.value
                ),
                "budget_bound_minus_sm_formula": float(
                    self.anchor.budget_band.limit_minus_sm_sd_formula
                ),
                "budget_sm_subtracted": self.anchor.budget_band.sm_subtracted,
                "budget_long_distance_dominated": (
                    self.anchor.budget_band.long_distance_dominated
                ),
                "budget_source": self.anchor.budget_band.source,
                "budget_construction": self.anchor.budget_band.construction,
                "parametrization_citation": (
                    RARE_KAON_KSHORT_MUMU_PARAMETRIZATION_CITATION
                ),
                "k006_parametrization_citation": (
                    RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION
                ),
                "rs_matching_assumption": RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
                "k006_matching_assumption_reused": RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
                "needs_human_physics": (
                    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1
                    + " "
                    + RARE_KAON_KSHORT_MUMU_RS_MATCHING_ASSUMPTION_V1
                ),
                "rs_semileptonic_vector_matching_status": (
                    RARE_KAON_RS_SEMILEPTONIC_VECTOR_MATCHING_STATUS_V1
                ),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "kappa_mu": float(result.kappa_mu),
                "p_c_y": float(result.p_c_y),
                "y_t": float(result.y_t),
                "lambda_wolfenstein": float(result.lambda_wolfenstein),
                "lambda_c": complex(result.lambda_c),
                "lambda_t": complex(result.lambda_t),
                "y_eff": complex(result.y_eff),
                "y_np_total": complex(result.y_np_total),
                "short_distance_amplitude": float(result.short_distance_amplitude),
                "np_shift_branching_fraction": float(
                    result.np_shift_branching_fraction
                ),
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=predicted,
            sm_prediction=float(result.sm_branching_fraction),
            experimental=float(self.anchor.value),
            ratio=float(ratio),
            budget=budget,
            notes=(
                "BR(K_S -> mu+ mu-)_SD,ell=0 reuses the K006 Buras/Isidori "
                "s -> d mu+mu- C10/C10p Wilson path with the K_S imaginary "
                "CP projection Im[-lambda_c Y_c + lambda_t C10]^2. The HARD "
                "budget is the K012.yaml PDG/LHCb 90% CL total-rate limit, "
                "applied to the constrained short-distance component. The total "
                "rate is long-distance dominated and the RS matching is "
                "flagged NEEDS-HUMAN-PHYSICS."
            ),
            diagnostics=diagnostics,
        )
