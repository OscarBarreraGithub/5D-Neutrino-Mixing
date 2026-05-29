"""B006 - rare leptonic decay ``B0 -> mu+ mu-``.

Physics
-------
``BR(B0 -> mu+ mu-)`` is evaluated with the shared ``b -> q l l`` core
introduced for B005, reached only through
``flavor_catalog_constraints.physics_adapters.rare_b_meson``.  This module
selects the ``b_d`` entry point, i.e. ``lambda_t^d = V_tb V_td^*`` and the
``B_d`` hadronic/lifetime inputs.  The time-integrated branching fraction uses
the same amplitude-dependent ``A_DeltaGamma`` treatment as B005; for ``B_d``
the width-difference input is ``y_d ~= 0``, so the SM-like time factor is one.

Severity
--------
HARD.  The experimental object is an upper limit, not a measured central
value.  The veto is therefore direction-aware: only upward excursions above
the YAML SM anchor are compared with the YAML CMS/PDG room
``BR_limit - BR_SM(anchor)``.  Large downward interference is not excluded by
an upper bound.  The LHCb direct limit is loaded from the same sidecar and
kept as a documented comparison; the active budget is the strongest CMS/PDG
90%-CL limit recorded in ``B006.yaml``.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B006.yaml`` is the source of truth for the
PDG/CMS upper limit, the LHCb comparison limit, the HFLAV ratio provenance,
and the Bobeth et al. SM prediction.  Numeric values below are loaded from
that sidecar, not hardcoded here.

NEEDS-HUMAN-PHYSICS
-------------------
The RS contribution uses the shared documented Z/KK-penguin proxy because the
``ParameterPoint`` does not provide the full electroweak KK/Z/Z', muon,
Higgs/radion, scalar, or pseudoscalar matching inputs needed for a rigorous
RS ``b -> d mu mu`` calculation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_b_meson import (
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    bd_mumu_from_couplings,
    rare_b_dilepton_default_sm_inputs,
    rare_b_dilepton_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPECTED_BRANCHING_UNITS = "branching fraction"
_EXPECTED_RATIO_UNITS = "dimensionless ratio"

_EXPERIMENTAL_LIMIT_CANDIDATES = ("canonical_experimental_limit",)
_HFLAV_RATIO_CANDIDATES = ("hflav_equivalent_ratio",)
_SM_ANCHOR_CANDIDATES = ("standard_model_prediction",)
_INPUT_MEASUREMENTS_KEY = "input_measurements"
_CMS_SOURCE_KEY = "CMS2023:BdMuMu"
_LHCB_SOURCE_KEY = "LHCb2022:BdMuMu"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B006.yaml "
    "canonical_experimental_limit/CMS2023:BdMuMu + LHCb2022:BdMuMu "
    "compared with standard_model_prediction"
)
_PARAMETRIZATION_CITATION = (
    "Buras b->d l l effective Hamiltonian; "
    "Bobeth et al. arXiv:1311.0903 SM B_d -> mu mu validation"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: full RS electroweak KK/Z/Z', muon, "
    "Higgs/radion, scalar and pseudoscalar b->d mu mu matching is not "
    "available on ParameterPoint; v1 uses the documented Z/KK-penguin "
    "C9/C10 proxy."
)


@dataclass(frozen=True)
class BranchingLimitAnchor:
    """Upper-limit anchor loaded through the scaffold ``load_anchor`` path."""

    anchor: Anchor
    confidence_level: float | None
    source_key: str | None
    raw_block: Mapping[str, Any]

    @property
    def block_key(self) -> str:
        return self.anchor.block_key

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def units(self) -> str | None:
        return self.anchor.units

    @property
    def source_url(self) -> str | None:
        return self.anchor.source_url

    @property
    def snapshot_path(self) -> str | None:
        return self.anchor.snapshot_path


@dataclass(frozen=True)
class BDMuMuBudgetBand:
    """Direction-aware B006 upward-NP budget in branching-fraction units."""

    source: str
    active_limit: float
    active_limit_block: str
    active_limit_source_key: str | None
    active_confidence_level: float | None
    cms_limit: float
    cms_limit_block: str
    lhcb_limit: float
    lhcb_limit_block: str
    sm_anchor_value: float
    sm_theory_sigma: float
    hard_veto_budget: float
    limit_minus_formula_sm: float
    construction: str


@dataclass(frozen=True)
class BDMuMuAnchor:
    """Typed B006 anchor: limits, HFLAV ratio provenance, SM, and budget."""

    experimental_limit: BranchingLimitAnchor
    cms_limit: BranchingLimitAnchor
    lhcb_limit: BranchingLimitAnchor
    hflav_ratio_limit: BranchingLimitAnchor
    standard_model: Anchor
    budget_band: BDMuMuBudgetBand

    @property
    def value(self) -> float:
        return self.experimental_limit.value

    @property
    def uncertainty(self) -> None:
        return None

    @property
    def sm_value(self) -> float:
        return self.standard_model.value

    @property
    def budget(self) -> float:
        return self.budget_band.hard_veto_budget


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: B006 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B006 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    block = data.get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: expected mapping-shaped 'pdg_or_equivalent' for B006, "
            f"got {type(block).__name__}"
        )
    return block


def _validate_positive_anchor(
    anchor: Anchor,
    *,
    process_id: str,
    label: str,
    expected_units: str | None,
    require_uncertainty: bool = False,
) -> None:
    if expected_units is not None and anchor.units != expected_units:
        raise AnchorError(
            f"{process_id}: {label} anchor {anchor.block_key!r} must use units "
            f"{expected_units!r}, got {anchor.units!r}"
        )
    if anchor.value <= 0.0 or not math.isfinite(anchor.value):
        raise AnchorError(
            f"{process_id}: {label} anchor {anchor.block_key!r} must have a "
            "positive finite value"
        )
    if require_uncertainty and (
        anchor.uncertainty is None
        or anchor.uncertainty <= 0.0
        or not math.isfinite(anchor.uncertainty)
    ):
        raise AnchorError(
            f"{process_id}: {label} uncertainty is required for the B006 budget"
        )


def _limit_anchor_from_top_level(
    process_id: str,
    *,
    candidates: tuple[str, ...],
    value_key: str,
    expected_units: str,
    label: str,
) -> BranchingLimitAnchor:
    anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=candidates,
        value_key=value_key,
    )
    _validate_positive_anchor(
        anchor,
        process_id=process_id,
        label=label,
        expected_units=expected_units,
    )
    block = _pdg_block(process_id).get(anchor.block_key)
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: selected B006 anchor block {anchor.block_key!r} "
            "is not a mapping"
        )
    return BranchingLimitAnchor(
        anchor=anchor,
        confidence_level=_optional_float(
            block.get("confidence_level"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.confidence_level",
        ),
        source_key=_optional_str(block.get("source_key")),
        raw_block=block,
    )


def _measurement_entry(
    process_id: str,
    *,
    source_key: str,
) -> tuple[int, Mapping[str, Any]]:
    entries = _pdg_block(process_id).get(_INPUT_MEASUREMENTS_KEY)
    if not isinstance(entries, list) or not entries:
        raise AnchorError(
            f"{process_id}: pdg_or_equivalent.{_INPUT_MEASUREMENTS_KEY} "
            "must be a non-empty list"
        )
    for index, entry in enumerate(entries):
        if not isinstance(entry, Mapping):
            raise AnchorError(
                f"{process_id}: {_INPUT_MEASUREMENTS_KEY}[{index}] is not a mapping"
            )
        if str(entry.get("source")) == source_key:
            return index, entry
    raise AnchorError(
        f"{process_id}: no B006 input_measurements entry with source={source_key!r}"
    )


def _limit_anchor_from_measurement(
    process_id: str,
    *,
    source_key: str,
    value_key: str,
) -> BranchingLimitAnchor:
    index, entry = _measurement_entry(process_id, source_key=source_key)
    block_key = f"{_INPUT_MEASUREMENTS_KEY}[{index}]"
    virtual_entry = dict(entry)
    virtual_entry["source_key"] = source_key
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
        anchor = load_anchor(
            process_id,
            family=_FAMILY,
            candidates=(block_key,),
            value_key=value_key,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block

    _validate_positive_anchor(
        anchor,
        process_id=process_id,
        label=source_key,
        expected_units=None,
    )
    return BranchingLimitAnchor(
        anchor=anchor,
        confidence_level=0.90 if value_key.endswith("_90cl") else None,
        source_key=source_key,
        raw_block=entry,
    )


def _build_budget_band(
    *,
    process_id: str,
    experimental_limit: BranchingLimitAnchor,
    cms_limit: BranchingLimitAnchor,
    lhcb_limit: BranchingLimitAnchor,
    standard_model: Anchor,
    formula_sm: float,
) -> BDMuMuBudgetBand:
    _validate_positive_anchor(
        standard_model,
        process_id=process_id,
        label="standard_model",
        expected_units=_EXPECTED_BRANCHING_UNITS,
        require_uncertainty=True,
    )
    limits = (experimental_limit, cms_limit, lhcb_limit)
    active = min(limits, key=lambda item: item.value)
    hard_budget = float(active.value - standard_model.value)
    if hard_budget <= 0.0 or not math.isfinite(hard_budget):
        raise AnchorError(
            f"{process_id}: B006 limit-SM budget must be positive and finite"
        )
    return BDMuMuBudgetBand(
        source=_BUDGET_SOURCE,
        active_limit=float(active.value),
        active_limit_block=active.block_key,
        active_limit_source_key=active.source_key,
        active_confidence_level=active.confidence_level,
        cms_limit=float(cms_limit.value),
        cms_limit_block=cms_limit.block_key,
        lhcb_limit=float(lhcb_limit.value),
        lhcb_limit_block=lhcb_limit.block_key,
        sm_anchor_value=float(standard_model.value),
        sm_theory_sigma=float(standard_model.uncertainty),
        hard_veto_budget=hard_budget,
        limit_minus_formula_sm=float(active.value - formula_sm),
        construction=(
            "Direction-aware upper-limit budget: upward excess above the "
            "B006.yaml SM anchor is compared with the strongest CMS/PDG "
            "90%-CL upper-limit room, BR_limit - BR_SM(anchor)."
        ),
    )


def _load_bd_mumu_anchor(process_id: str, *, formula_sm: float) -> BDMuMuAnchor:
    experimental_limit = _limit_anchor_from_top_level(
        process_id,
        candidates=_EXPERIMENTAL_LIMIT_CANDIDATES,
        value_key="upper_limit",
        expected_units=_EXPECTED_BRANCHING_UNITS,
        label="experimental_limit",
    )
    hflav_ratio_limit = _limit_anchor_from_top_level(
        process_id,
        candidates=_HFLAV_RATIO_CANDIDATES,
        value_key="upper_limit",
        expected_units=_EXPECTED_RATIO_UNITS,
        label="hflav_ratio_limit",
    )
    standard_model = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_SM_ANCHOR_CANDIDATES,
    )
    cms_limit = _limit_anchor_from_measurement(
        process_id,
        source_key=_CMS_SOURCE_KEY,
        value_key="upper_limit_90cl",
    )
    lhcb_limit = _limit_anchor_from_measurement(
        process_id,
        source_key=_LHCB_SOURCE_KEY,
        value_key="upper_limit_90cl",
    )
    return BDMuMuAnchor(
        experimental_limit=experimental_limit,
        cms_limit=cms_limit,
        lhcb_limit=lhcb_limit,
        hflav_ratio_limit=hflav_ratio_limit,
        standard_model=standard_model,
        budget_band=_build_budget_band(
            process_id=process_id,
            experimental_limit=experimental_limit,
            cms_limit=cms_limit,
            lhcb_limit=lhcb_limit,
            standard_model=standard_model,
            formula_sm=formula_sm,
        ),
    )


def _upper_limit_ratio(predicted: float, anchor: BDMuMuAnchor) -> tuple[float, float, bool]:
    upward_excess = max(0.0, float(predicted) - float(anchor.sm_value))
    budget = float(anchor.budget)
    ratio = upward_excess / budget if budget > 0.0 else float("inf")
    total_limit_passes = float(predicted) <= float(anchor.budget_band.active_limit)
    return float(upward_excess), float(ratio), bool(ratio <= 1.0 and total_limit_passes)


@register
class Constraint:
    """Catalogued ``B0 -> mu+ mu-`` branching-ratio constraint (B006)."""

    process_id = "B006"
    severity = Severity.HARD
    observable = "BR(B0 -> mu+ mu-)"

    def __init__(self) -> None:
        self.sm_inputs = rare_b_dilepton_default_sm_inputs()
        self.sm_result = rare_b_dilepton_sm_branching_fraction("b_d", self.sm_inputs)
        self.anchor = _load_bd_mumu_anchor(
            self.process_id,
            formula_sm=float(self.sm_result.branching_fraction),
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.sm_result.branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; B0 -> mu+ mu- "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "sm_anchor_branching_fraction": float(self.anchor.sm_value),
                    "sm_formula_branching_fraction": float(
                        self.sm_result.branching_fraction
                    ),
                    "sm_a_delta_gamma": float(self.sm_result.a_delta_gamma),
                    "sm_time_integrated_factor": float(
                        self.sm_result.time_integrated_factor
                    ),
                    "budget_source": self.anchor.budget_band.source,
                    "budget_construction": self.anchor.budget_band.construction,
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        result = bd_mumu_from_couplings(
            couplings,
            m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
            inputs=self.sm_inputs,
        )
        predicted = float(result.branching_fraction)
        upward_excess, ratio, passes = _upper_limit_ratio(predicted, self.anchor)
        total_minus_limit = predicted - float(self.anchor.budget_band.active_limit)
        total_limit_ratio = predicted / float(self.anchor.budget_band.active_limit)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "sm_anchor_branching_fraction": float(self.anchor.sm_value),
                "sm_formula_branching_fraction": float(
                    result.sm_branching_fraction
                ),
                "sm_formula_minus_anchor": float(
                    result.sm_branching_fraction - self.anchor.sm_value
                ),
                "experimental_block": self.anchor.experimental_limit.block_key,
                "experimental_source_key": self.anchor.experimental_limit.source_key,
                "experimental_confidence_level": (
                    self.anchor.experimental_limit.confidence_level
                ),
                "active_limit": float(self.anchor.budget_band.active_limit),
                "active_limit_block": self.anchor.budget_band.active_limit_block,
                "active_limit_source_key": (
                    self.anchor.budget_band.active_limit_source_key
                ),
                "cms_90cl_limit": float(self.anchor.budget_band.cms_limit),
                "cms_limit_block": self.anchor.budget_band.cms_limit_block,
                "lhcb_90cl_limit": float(self.anchor.budget_band.lhcb_limit),
                "lhcb_limit_block": self.anchor.budget_band.lhcb_limit_block,
                "hflav_bd_over_bs_ratio_90cl_limit": float(
                    self.anchor.hflav_ratio_limit.value
                ),
                "hflav_ratio_block": self.anchor.hflav_ratio_limit.block_key,
                "sm_block": self.anchor.standard_model.block_key,
                "sm_theory_sigma": float(self.anchor.budget_band.sm_theory_sigma),
                "upward_excess_over_sm_anchor": float(upward_excess),
                "hard_veto_np_shift_budget": float(
                    self.anchor.budget_band.hard_veto_budget
                ),
                "limit_minus_formula_sm": float(
                    self.anchor.budget_band.limit_minus_formula_sm
                ),
                "total_minus_active_limit": float(total_minus_limit),
                "total_limit_ratio": float(total_limit_ratio),
                "budget_source": self.anchor.budget_band.source,
                "budget_construction": self.anchor.budget_band.construction,
                "budget_policy": (
                    "max(0, BR_total - BR_SM(anchor)) / "
                    "(BR_limit - BR_SM(anchor)); additionally require "
                    "BR_total <= BR_limit"
                ),
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "wilson_coefficients": (
                    {}
                    if result.wilsons is None
                    else {
                        key: complex(value)
                        for key, value in result.wilsons.wilsons.items()
                    }
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
            budget=float(self.anchor.budget),
            notes=(
                "BR(B0 -> mu+ mu-) reuses the shared Buras b->qll "
                "C10-dominant formula with the b_d CKM factor and "
                "amplitude-dependent A_DeltaGamma time integration. The RS "
                "contribution is the documented Z/KK-penguin C9/C10 proxy "
                "from mass-basis b-d couplings; full EW/lepton/scalar "
                "matching is marked NEEDS-HUMAN-PHYSICS. The HARD ratio is "
                "the upward excess above the YAML SM anchor over the "
                "CMS/PDG upper-limit room from B006.yaml."
            ),
            diagnostics=diagnostics,
        )
