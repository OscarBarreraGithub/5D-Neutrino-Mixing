"""B022 - charged rare beauty decay ``B+ -> K+ nu nubar``.

Physics
-------
The short-distance SM piece is evaluated with the shared ``b -> s nu nubar``
core in ``quarkConstraints.rare_b_nunu`` through the adapter boundary.  The
core uses the Inami-Lim top function ``X_t`` and the standard
``C_L^SM = -X_t / sin^2(theta_W)`` convention.  The exclusive
``B+ -> K+`` SM branching fraction is not hardcoded here: it is loaded from
``B022.yaml`` and supplied to the core as the total SM normalization.  The
charged-mode long-distance term is kept fixed while new physics rescales only
the short-distance remainder.  New physics is the Phase-4d
``rs_semileptonic_wilsons.b_to_s_nunu`` block mapped as ``X_NP=C/g_SM^2``.

Severity
--------
HARD.  The total predicted branching fraction is compared with the Belle II
evidence measurement recorded in ``B022.yaml``.  The direction-aware one-sigma
budget combines the Belle II statistical/systematic uncertainty and the HPQCD
SM-normalization uncertainty in quadrature.  Because Belle II is above the SM
central prediction, the SM-limit point is a measured tension rather than a
guaranteed passing point.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B022.yaml`` is the source of truth for the
Belle II measurement and HPQCD SM reference.  Numeric values below are loaded
from that sidecar, not hardcoded in this constraint.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_full_yaml,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_b_nunu import (
    RARE_B_NUNU_RS_SEMILEPTONIC_MATCHING_STATUS_V1,
    RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
    bplus_kplus_nunu_from_rs_semileptonic_wilsons,
    rare_b_nunu_default_sm_inputs,
    rare_b_nunu_sm_branching_fraction,
    rare_b_nunu_sm_inputs_with_bplus_kplus_normalization,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPECTED_UNITS = "branching fraction"
_BELLE_II_OBSERVABLE_NAME = "Belle II evidence measurement"
_HPQCD_SM_VALUE_ID = "HPQCD2023:B022:sm_prediction"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B022.yaml "
    "pdg_or_equivalent.observables[Belle II evidence measurement] "
    "+ pdg_or_equivalent.values[HPQCD2023:B022:sm_prediction]"
)
_PARAMETRIZATION_CITATION = (
    "Inami-Lim X_t short-distance b -> s nu nubar response; "
    "HPQCD 2023 arXiv:2207.13371 SM B+ -> K+ nu nubar normalization"
)
_TREE_LEVEL_STATUS = (
    "Phase-4d rigorous light-Z active-neutrino Wilson from "
    "rs_semileptonic_wilsons.b_to_s_nunu; old one-Z-like proxy resolved."
)
_UNEVALUATED_REASON = "rs_semileptonic_wilsons.b_to_s_nunu not provided"
_SCAFFOLD_NO_SCALAR_UNCERTAINTY_KEY = "__b022_no_scalar_uncertainty__"


@dataclass(frozen=True)
class BranchingFractionAnchor:
    """Branching-fraction value with symmetric or asymmetric uncertainties."""

    block_key: str
    entry_id: str
    source: str | None
    year: int | None
    value: float
    uncertainty_upper: float
    uncertainty_lower: float
    units: str | None
    source_key: str | None
    source_url: str | None
    snapshot_path: str | None
    note: str | None = None

    @property
    def uncertainty(self) -> float:
        return float(0.5 * (self.uncertainty_upper + self.uncertainty_lower))


@dataclass(frozen=True)
class RareBBudgetBand:
    """Direction-aware B022 total-branching-fraction budget."""

    source: str
    central_residual: float
    experimental_sigma_upper: float
    experimental_sigma_lower: float
    sm_theory_sigma: float
    combined_sigma_upper: float
    combined_sigma_lower: float
    hard_veto_budget: float
    lower_edge: float
    upper_edge: float


@dataclass(frozen=True)
class B022Anchor:
    """Typed B022 anchor: Belle II experiment, HPQCD SM reference, and budget."""

    experimental: BranchingFractionAnchor
    standard_model: BranchingFractionAnchor
    budget_band: RareBBudgetBand

    @property
    def value(self) -> float:
        return self.experimental.value

    @property
    def uncertainty(self) -> float:
        return self.experimental.uncertainty

    @property
    def sm_value(self) -> float:
        return self.standard_model.value

    @property
    def budget(self) -> float:
        return self.budget_band.hard_veto_budget


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
            f"{process_id}: B022 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B022 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: B022 anchor field {field_name!r} <= 0")
    return number


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    block = data.get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: expected mapping-shaped 'pdg_or_equivalent' for B022, "
            f"got {type(block).__name__}"
        )
    return block


def _entry_list(
    pdg_block: Mapping[str, Any],
    key: str,
    *,
    process_id: str,
) -> tuple[Mapping[str, Any], ...]:
    entries = pdg_block.get(key)
    if not isinstance(entries, list | tuple):
        raise AnchorError(f"{process_id}: pdg_or_equivalent.{key} must be a list")
    out: list[Mapping[str, Any]] = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: pdg_or_equivalent.{key}[{index}] is not a mapping"
            )
        out.append(entry)
    return tuple(out)


def _find_entry_index(
    entries: tuple[Mapping[str, Any], ...],
    *,
    key: str,
    expected: str,
    process_id: str,
    list_name: str,
) -> tuple[int, Mapping[str, Any]]:
    for index, entry in enumerate(entries):
        if str(entry.get(key)) == expected:
            return index, entry
    raise AnchorError(
        f"{process_id}: no pdg_or_equivalent.{list_name} entry with "
        f"{key}={expected!r}"
    )


def _validate_units(
    entry: Mapping[str, Any],
    *,
    process_id: str,
    field_name: str,
) -> str:
    units = _optional_str(entry.get("units"))
    if units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: {field_name} must use units {_EXPECTED_UNITS!r}, "
            f"got {units!r}"
        )
    return units


def _load_scaffold_list_anchor(
    section_key: str,
    *,
    match_key: str,
    expected: str,
    process_id: str,
    uncertainty_key: str = "uncertainty",
) -> tuple[Anchor, Mapping[str, Any], Mapping[str, Any]]:
    pdg_block = _pdg_block(process_id)
    index, entry = _find_entry_index(
        _entry_list(pdg_block, section_key, process_id=process_id),
        key=match_key,
        expected=expected,
        process_id=process_id,
        list_name=section_key,
    )
    block_key = f"pdg_or_equivalent.{section_key}[{index}]"
    virtual_block = {block_key: dict(entry)}
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
        scaffold_anchor = anchor_scaffold.load_anchor(
            process_id,
            family=_FAMILY,
            candidates=(block_key,),
            uncertainty_key=uncertainty_key,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    if scaffold_anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {scaffold_anchor.block_key!r}, "
            f"expected {block_key!r} for B022 {section_key} entry {expected!r}"
        )
    return scaffold_anchor, entry, pdg_block


def _load_belle_ii_measurement(
    *,
    process_id: str,
) -> BranchingFractionAnchor:
    scaffold_anchor, entry, pdg_block = _load_scaffold_list_anchor(
        "observables",
        match_key="name",
        expected=_BELLE_II_OBSERVABLE_NAME,
        process_id=process_id,
        uncertainty_key=_SCAFFOLD_NO_SCALAR_UNCERTAINTY_KEY,
    )
    value = _positive_float(
        scaffold_anchor.value,
        process_id=process_id,
        field_name="Belle II evidence measurement.value",
    )
    stat = _positive_float(
        entry.get("statistical_uncertainty"),
        process_id=process_id,
        field_name="Belle II evidence measurement.statistical_uncertainty",
    )
    syst_up = _positive_float(
        entry.get("systematic_uncertainty_positive"),
        process_id=process_id,
        field_name="Belle II evidence measurement.systematic_uncertainty_positive",
    )
    syst_down = _positive_float(
        entry.get("systematic_uncertainty_negative"),
        process_id=process_id,
        field_name="Belle II evidence measurement.systematic_uncertainty_negative",
    )
    return BranchingFractionAnchor(
        block_key=scaffold_anchor.block_key,
        entry_id=_BELLE_II_OBSERVABLE_NAME,
        source=_optional_str(pdg_block.get("source")),
        year=_optional_int(pdg_block.get("year")),
        value=float(value),
        uncertainty_upper=float(math.hypot(stat, syst_up)),
        uncertainty_lower=float(math.hypot(stat, syst_down)),
        units=_validate_units(
            entry,
            process_id=process_id,
            field_name="Belle II evidence measurement",
        ),
        source_key=_optional_str(entry.get("source_key")),
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        note=_optional_str(entry.get("sm_excess")),
    )


def _load_hpqcd_sm_reference(
    *,
    process_id: str,
) -> BranchingFractionAnchor:
    scaffold_anchor, entry, _ = _load_scaffold_list_anchor(
        "values",
        match_key="value_id",
        expected=_HPQCD_SM_VALUE_ID,
        process_id=process_id,
    )
    value = _positive_float(
        scaffold_anchor.value,
        process_id=process_id,
        field_name=f"{_HPQCD_SM_VALUE_ID}.value",
    )
    uncertainty = _positive_float(
        scaffold_anchor.uncertainty,
        process_id=process_id,
        field_name=f"{_HPQCD_SM_VALUE_ID}.uncertainty",
    )
    return BranchingFractionAnchor(
        block_key=scaffold_anchor.block_key,
        entry_id=_HPQCD_SM_VALUE_ID,
        source=_optional_str(entry.get("source_key")),
        year=scaffold_anchor.year,
        value=float(value),
        uncertainty_upper=float(uncertainty),
        uncertainty_lower=float(uncertainty),
        units=_validate_units(entry, process_id=process_id, field_name=_HPQCD_SM_VALUE_ID),
        source_key=_optional_str(entry.get("source_key")),
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        note=_optional_str(entry.get("note")),
    )


def _build_budget_band(
    *,
    experimental: BranchingFractionAnchor,
    standard_model: BranchingFractionAnchor,
) -> RareBBudgetBand:
    combined_upper = math.sqrt(
        experimental.uncertainty_upper**2 + standard_model.uncertainty_upper**2
    )
    combined_lower = math.sqrt(
        experimental.uncertainty_lower**2 + standard_model.uncertainty_lower**2
    )
    lower_edge = experimental.value - combined_lower
    upper_edge = experimental.value + combined_upper
    hard_budget = max(combined_upper, combined_lower)
    if hard_budget <= 0.0 or lower_edge <= 0.0:
        raise AnchorError("B022: constructed branching-ratio budget is invalid")
    return RareBBudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(experimental.value - standard_model.value),
        experimental_sigma_upper=float(experimental.uncertainty_upper),
        experimental_sigma_lower=float(experimental.uncertainty_lower),
        sm_theory_sigma=float(standard_model.uncertainty),
        combined_sigma_upper=float(combined_upper),
        combined_sigma_lower=float(combined_lower),
        hard_veto_budget=float(hard_budget),
        lower_edge=float(lower_edge),
        upper_edge=float(upper_edge),
    )


def _load_b022_anchor(process_id: str) -> B022Anchor:
    experimental = _load_belle_ii_measurement(process_id=process_id)
    standard_model = _load_hpqcd_sm_reference(process_id=process_id)
    return B022Anchor(
        experimental=experimental,
        standard_model=standard_model,
        budget_band=_build_budget_band(
            experimental=experimental,
            standard_model=standard_model,
        ),
    )


def _selected_budget(
    predicted: float,
    anchor: B022Anchor,
) -> tuple[float, float, bool]:
    pull = float(predicted - anchor.value)
    budget = (
        anchor.budget_band.combined_sigma_upper
        if pull >= 0.0
        else anchor.budget_band.combined_sigma_lower
    )
    ratio = abs(pull) / budget if budget > 0.0 else float("inf")
    return float(budget), float(ratio), bool(ratio <= 1.0)


@register
class Constraint:
    """Catalogued charged rare-B invisible branching-ratio constraint (B022)."""

    process_id = "B022"
    severity = Severity.HARD
    observable = "BR(B+ -> K+ nu nubar)"

    def __init__(self) -> None:
        self.anchor = _load_b022_anchor(self.process_id)
        base_inputs = rare_b_nunu_default_sm_inputs()
        self.sm_inputs = rare_b_nunu_sm_inputs_with_bplus_kplus_normalization(
            self.anchor.sm_value,
            inputs=base_inputs,
        )
        self.sm_result = rare_b_nunu_sm_branching_fraction(self.sm_inputs)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                predicted=None,
                sm_prediction=float(self.sm_result.branching_fraction),
                experimental=float(self.anchor.value),
                ratio=None,
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - rs_semileptonic_wilsons.b_to_s_nunu "
                    "absent; B+ -> K+ nu nubar constraint is non-vetoing."
                ),
                diagnostics={
                    "evaluated": False,
                    "unevaluated_reason": _UNEVALUATED_REASON,
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_quark_mass_basis_couplings_present": (
                        point.get_extra("quark_mass_basis_couplings") is not None
                    ),
                    "sm_anchor_branching_fraction": float(self.anchor.sm_value),
                    "sm_formula_branching_fraction": float(
                        self.sm_result.branching_fraction
                    ),
                    "budget_source": self.anchor.budget_band.source,
                    "tree_level_status": _TREE_LEVEL_STATUS,
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = bplus_kplus_nunu_from_rs_semileptonic_wilsons(
                couplings,
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
                inputs=self.sm_inputs,
            )
        except (AttributeError, KeyError, TypeError, ValueError) as exc:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                predicted=None,
                sm_prediction=float(self.sm_result.branching_fraction),
                experimental=float(self.anchor.value),
                ratio=None,
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_semileptonic_wilsons for "
                    "B+ -> K+ nu nubar"
                ),
                diagnostics={
                    "evaluated": False,
                    "unevaluated_reason": _UNEVALUATED_REASON,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "tree_level_status": _TREE_LEVEL_STATUS,
                },
            )
        predicted = float(result.branching_fraction)
        budget, ratio, passes = _selected_budget(predicted, self.anchor)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                "sm_anchor_branching_fraction": float(self.anchor.sm_value),
                "sm_formula_branching_fraction": float(
                    result.sm_branching_fraction
                ),
                "sm_formula_minus_anchor": float(
                    result.sm_branching_fraction - self.anchor.sm_value
                ),
                "experimental_block": self.anchor.experimental.block_key,
                "experimental_entry_id": self.anchor.experimental.entry_id,
                "sm_block": self.anchor.standard_model.block_key,
                "sm_entry_id": self.anchor.standard_model.entry_id,
                "experimental_sigma_upper": float(
                    self.anchor.budget_band.experimental_sigma_upper
                ),
                "experimental_sigma_lower": float(
                    self.anchor.budget_band.experimental_sigma_lower
                ),
                "sm_theory_sigma": float(self.anchor.budget_band.sm_theory_sigma),
                "budget_combined_sigma_upper": float(
                    self.anchor.budget_band.combined_sigma_upper
                ),
                "budget_combined_sigma_lower": float(
                    self.anchor.budget_band.combined_sigma_lower
                ),
                "budget_lower_edge": float(self.anchor.budget_band.lower_edge),
                "budget_upper_edge": float(self.anchor.budget_band.upper_edge),
                "budget_source": self.anchor.budget_band.source,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
                "tree_level_status": _TREE_LEVEL_STATUS,
                "rs_semileptonic_nunu_matching_status": (
                    RARE_B_NUNU_RS_SEMILEPTONIC_MATCHING_STATUS_V1
                ),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "lambda_wolfenstein": float(result.lambda_wolfenstein),
                "lambda_t_bs": complex(result.lambda_t_bs),
                "c_l_sm": complex(result.c_l_sm),
                "c_l_total": complex(result.c_l_total),
                "c_r_total": complex(result.c_r_total),
                "x_eff_left": complex(result.x_eff_left),
                "x_eff_right": complex(result.x_eff_right),
                "x_t": float(result.x_t),
                "epsilon": float(result.epsilon),
                "eta": float(result.eta),
                "r_k": float(result.r_k),
                "r_kstar": float(result.r_kstar),
                "np_shift_branching_fraction": float(
                    result.np_shift_branching_fraction
                ),
                "down_sector_indices": (1, 2),
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=predicted,
            sm_prediction=float(result.sm_branching_fraction),
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "BR(B+ -> K+ nu nubar) uses the shared b -> s nu nubar "
                "short-distance X_t response normalized to the HPQCD SM "
                "branching fraction from B022.yaml. RS contribution is the "
                "Phase-4d b_to_s_nunu Wilson block mapped as X_NP=C/g_SM^2; "
                "the old one-Z-like proxy is not used. Budget combines Belle "
                "II and SM uncertainties in quadrature."
            ),
            diagnostics=diagnostics,
        )
