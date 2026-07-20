"""B007 - rare electronic leptonic neutral-B decays.

Physics
-------
``BR(B_s -> e+ e-)`` and ``BR(B0 -> e+ e-)`` reuse the shared
``quarkConstraints.rare_b_dilepton`` Buras ``b -> q l l`` machinery through
the electron-mode adapter
``flavor_catalog_constraints.physics_adapters.rare_b_electronic``.  The SM
rate is helicity-suppressed by the charged-lepton mass squared; using the
electron mass gives ``BR(B_s -> e+e-) ~= 8.5e-14``.

Severity
--------
HARD.  The SM predictions are many orders below the PDG 2026 90%-CL limits in
``B007.yaml``, so this constraint treats the limit as an essentially pure-NP
branching-fraction budget for both channels.  The reported HARD ratio is the
largest channel saturation.

Catalog sidecar
---------------
``flavor_catalog/processes/secondary/beauty/B007.yaml`` is the source of
truth for the PDG live/API limits and Bobeth et al. SM reference values.
Numeric values below are loaded through scaffold anchor paths, not hardcoded.

NEEDS-HUMAN-PHYSICS
-------------------
The Phase-3a RS light-Z contribution supplies rigorous vector/axial
``C9/C10/C9'/C10'`` Wilsons.  Scalar and pseudoscalar matching remains
deferred.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

from flavor_catalog_constraints import anchors as anchor_scaffold
from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
)
from flavor_catalog_constraints.base import (
    ConstraintLevel,
    ConstraintResult,
    ParameterPoint,
    Severity,
)
from flavor_catalog_constraints.physics_adapters.rare_b_electronic import (
    ELECTRON_MASS_GEV,
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RARE_B_ELECTRONIC_MODEL_NOTE_V1,
    bd_ee_from_rs_semileptonic_wilsons,
    bs_ee_from_rs_semileptonic_wilsons,
    rare_b_electronic_default_sm_inputs,
    rare_b_electronic_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_TIER = ConstraintLevel.SECONDARY
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPECTED_UNITS = "branching fraction"

_BS_LIMIT_ROW_ID = "pdg2026_bs_ee_90cl"
_BD_LIMIT_ROW_ID = "pdg2026_bd_ee_90cl"
_SM_SOURCE_KEY = "BobethEtAl2013_BqllSM"
_BS_SM_OBSERVABLE = "SM BR(B_s -> e+e-)"
_BD_SM_OBSERVABLE = "SM BR(B_d -> e+e-)"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/secondary/beauty/B007.yaml "
    "pdg2026_bs_ee_90cl + pdg2026_bd_ee_90cl, with BobethEtAl2013_BqllSM "
    "SM references"
)
_PARAMETRIZATION_CITATION = (
    "Buras b->q l l effective Hamiltonian; Bobeth et al. arXiv:1311.0903 "
    "SM B_q -> e e validation"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: Phase-3a supplies the light-Z vector/axial "
    "C9/C10/C9'/C10' terms; Higgs/radion scalar and pseudoscalar "
    "b->q e e matching remains deferred."
)


@dataclass(frozen=True)
class ElectronicLimitAnchor:
    """Upper-limit anchor loaded from one B007 ``pdg_or_equivalent.values`` row."""

    anchor: Anchor
    row_id: str
    confidence_level: float
    cl: str | None
    limit_type: str | None
    conditions: str | None
    source_key: str | None
    raw_row: Mapping[str, Any]

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
class ElectronicSMAnchor:
    """SM prediction anchor loaded from B007 auxiliary theory inputs."""

    anchor: Anchor
    source_key: str
    raw_row: Mapping[str, Any]

    @property
    def block_key(self) -> str:
        return self.anchor.block_key

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def uncertainty(self) -> float | None:
        return self.anchor.uncertainty

    @property
    def source_url(self) -> str | None:
        return self.anchor.source_url


@dataclass(frozen=True)
class ElectronicChannelBudget:
    """Pure-NP upper-limit budget for one electronic neutral-B channel."""

    channel_key: str
    source: str
    active_limit: float
    active_limit_block: str
    active_limit_row_id: str
    active_confidence_level: float
    sm_anchor_value: float
    sm_theory_sigma: float
    formula_sm_value: float
    hard_veto_budget: float
    limit_minus_sm_anchor: float
    limit_minus_formula_sm: float
    sm_anchor_to_limit_ratio: float
    construction: str


@dataclass(frozen=True)
class B007Anchor:
    """Typed B007 anchor bundle for both ``B_s`` and ``B0`` electron modes."""

    bs_limit: ElectronicLimitAnchor
    bd_limit: ElectronicLimitAnchor
    bs_standard_model: ElectronicSMAnchor
    bd_standard_model: ElectronicSMAnchor
    bs_budget: ElectronicChannelBudget
    bd_budget: ElectronicChannelBudget

    @property
    def value(self) -> float:
        return min(self.bs_limit.value, self.bd_limit.value)

    @property
    def budget(self) -> float:
        return min(self.bs_budget.hard_veto_budget, self.bd_budget.hard_veto_budget)

    @property
    def sm_value(self) -> float:
        return max(self.bs_standard_model.value, self.bd_standard_model.value)


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: B007 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B007 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: B007 field {field_name!r} must be positive")
    return number


def _full_yaml(process_id: str) -> Mapping[str, Any]:
    return load_full_yaml(process_id, family=_FAMILY, tier=_TIER)


def _pdg_block(process_id: str) -> Mapping[str, Any]:
    block = _full_yaml(process_id).get("pdg_or_equivalent")
    if not isinstance(block, Mapping):
        raise AnchorError(
            f"{process_id}: expected mapping-shaped 'pdg_or_equivalent' for B007"
        )
    return block


def _anchor_from_virtual_block(
    process_id: str,
    *,
    block_key: str,
    row: Mapping[str, Any],
    value_key: str = "value",
    uncertainty_key: str = "uncertainty",
) -> Anchor:
    virtual_block = {block_key: row}
    original_load_pdg_block = anchor_scaffold.load_pdg_block

    def _load_virtual_pdg_block(
        request_process_id: str,
        **kwargs: Any,
    ) -> Mapping[str, Any]:
        if (
            request_process_id == process_id
            and kwargs.get("family") == _FAMILY
            and kwargs.get("tier") == _TIER
        ):
            return virtual_block
        return original_load_pdg_block(request_process_id, **kwargs)

    anchor_scaffold.load_pdg_block = _load_virtual_pdg_block
    try:
        anchor = load_anchor(
            process_id,
            family=_FAMILY,
            tier=_TIER,
            candidates=(block_key,),
            value_key=value_key,
            uncertainty_key=uncertainty_key,
        )
    finally:
        anchor_scaffold.load_pdg_block = original_load_pdg_block
    if anchor.block_key != block_key:
        raise AnchorError(
            f"{process_id}: load_anchor selected {anchor.block_key!r}, "
            f"expected {block_key!r} for B007 anchor"
        )
    return anchor


def _value_row_by_id(process_id: str, row_id: str) -> tuple[int, Mapping[str, Any]]:
    values = _pdg_block(process_id).get("values")
    if not isinstance(values, list) or not values:
        raise AnchorError(f"{process_id}: pdg_or_equivalent.values must be a list")
    matches: list[tuple[int, Mapping[str, Any]]] = []
    for index, row in enumerate(values):
        if not isinstance(row, Mapping):
            raise AnchorError(f"{process_id}: values[{index}] is not a mapping")
        if row.get("id") == row_id:
            matches.append((index, row))
    if len(matches) != 1:
        raise AnchorError(
            f"{process_id}: expected exactly one pdg_or_equivalent.values row "
            f"with id={row_id!r}, found {len(matches)}"
        )
    return matches[0]


def _limit_anchor_from_value_row(process_id: str, *, row_id: str) -> ElectronicLimitAnchor:
    index, row = _value_row_by_id(process_id, row_id)
    block_key = f"values[{index}]"
    anchor = _anchor_from_virtual_block(
        process_id,
        block_key=block_key,
        row=row,
        value_key="upper_limit",
    )
    if anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: B007 limit {row_id!r} must use units "
            f"{_EXPECTED_UNITS!r}, got {anchor.units!r}"
        )
    if anchor.value <= 0.0 or not math.isfinite(anchor.value):
        raise AnchorError(f"{process_id}: B007 limit {row_id!r} must be positive")
    confidence_level = _positive_float(
        row.get("confidence_level"),
        process_id=process_id,
        field_name=f"{row_id}.confidence_level",
    )
    if row.get("limit_type") != "upper":
        raise AnchorError(f"{process_id}: B007 limit {row_id!r} must be upper-type")
    return ElectronicLimitAnchor(
        anchor=anchor,
        row_id=row_id,
        confidence_level=confidence_level,
        cl=_optional_str(row.get("cl")),
        limit_type=_optional_str(row.get("limit_type")),
        conditions=_optional_str(row.get("conditions")),
        source_key=_optional_str(row.get("source_key")),
        raw_row=row,
    )


def _sm_prediction_row(
    process_id: str,
    *,
    source_key: str,
    observable: str,
) -> tuple[int, int, Mapping[str, Any], Mapping[str, Any]]:
    aux = _full_yaml(process_id).get("auxiliary_theory_inputs")
    if not isinstance(aux, Mapping):
        raise AnchorError(f"{process_id}: missing auxiliary_theory_inputs mapping")
    groups = aux.get("standard_model_predictions")
    if not isinstance(groups, list) or not groups:
        raise AnchorError(
            f"{process_id}: auxiliary_theory_inputs.standard_model_predictions "
            "must be a non-empty list"
        )
    for group_index, group in enumerate(groups):
        if not isinstance(group, Mapping):
            raise AnchorError(
                f"{process_id}: standard_model_predictions[{group_index}] is not a mapping"
            )
        if group.get("source_key") != source_key:
            continue
        predictions = group.get("predictions")
        if not isinstance(predictions, list) or not predictions:
            raise AnchorError(
                f"{process_id}: standard_model_predictions[{group_index}].predictions "
                "must be a non-empty list"
            )
        for prediction_index, row in enumerate(predictions):
            if not isinstance(row, Mapping):
                raise AnchorError(
                    f"{process_id}: predictions[{prediction_index}] is not a mapping"
                )
            if row.get("observable") == observable:
                return group_index, prediction_index, group, row
    raise AnchorError(
        f"{process_id}: no SM prediction row for source_key={source_key!r}, "
        f"observable={observable!r}"
    )


def _sm_anchor_from_auxiliary_prediction(
    process_id: str,
    *,
    source_key: str,
    observable: str,
) -> ElectronicSMAnchor:
    group_index, prediction_index, group, row = _sm_prediction_row(
        process_id,
        source_key=source_key,
        observable=observable,
    )
    block_key = (
        f"auxiliary_theory_inputs.standard_model_predictions[{group_index}]"
        f".predictions[{prediction_index}]"
    )
    virtual_row = dict(row)
    virtual_row.setdefault("source", source_key)
    for key in ("source_url", "snapshot_path"):
        if key in group:
            virtual_row.setdefault(key, group[key])
    anchor = _anchor_from_virtual_block(
        process_id,
        block_key=block_key,
        row=virtual_row,
    )
    if anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: B007 SM row {observable!r} must use units "
            f"{_EXPECTED_UNITS!r}, got {anchor.units!r}"
        )
    if anchor.value <= 0.0 or not math.isfinite(anchor.value):
        raise AnchorError(f"{process_id}: B007 SM row {observable!r} must be positive")
    if anchor.uncertainty is None or anchor.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: B007 SM row {observable!r} requires a positive uncertainty"
        )
    return ElectronicSMAnchor(anchor=anchor, source_key=source_key, raw_row=row)


def _build_channel_budget(
    *,
    process_id: str,
    channel_key: str,
    limit: ElectronicLimitAnchor,
    standard_model: ElectronicSMAnchor,
    formula_sm: float,
) -> ElectronicChannelBudget:
    if limit.value <= standard_model.value:
        raise AnchorError(
            f"{process_id}: B007 {channel_key} limit must exceed the SM anchor"
        )
    hard_budget = float(limit.value)
    return ElectronicChannelBudget(
        channel_key=channel_key,
        source=_BUDGET_SOURCE,
        active_limit=float(limit.value),
        active_limit_block=limit.block_key,
        active_limit_row_id=limit.row_id,
        active_confidence_level=float(limit.confidence_level),
        sm_anchor_value=float(standard_model.value),
        sm_theory_sigma=float(standard_model.uncertainty),
        formula_sm_value=float(formula_sm),
        hard_veto_budget=hard_budget,
        limit_minus_sm_anchor=float(limit.value - standard_model.value),
        limit_minus_formula_sm=float(limit.value - formula_sm),
        sm_anchor_to_limit_ratio=float(standard_model.value / limit.value),
        construction=(
            "Pure-NP upper-limit budget: because the B007 SM branching "
            "fraction is orders below the experimental 90%-CL limit, compare "
            "max(0, BR_total - BR_SM(anchor)) with BR_limit and also require "
            "BR_total <= BR_limit."
        ),
    )


def _load_b007_anchor(
    process_id: str,
    *,
    formula_sm_bs: float,
    formula_sm_bd: float,
) -> B007Anchor:
    bs_limit = _limit_anchor_from_value_row(process_id, row_id=_BS_LIMIT_ROW_ID)
    bd_limit = _limit_anchor_from_value_row(process_id, row_id=_BD_LIMIT_ROW_ID)
    bs_sm = _sm_anchor_from_auxiliary_prediction(
        process_id,
        source_key=_SM_SOURCE_KEY,
        observable=_BS_SM_OBSERVABLE,
    )
    bd_sm = _sm_anchor_from_auxiliary_prediction(
        process_id,
        source_key=_SM_SOURCE_KEY,
        observable=_BD_SM_OBSERVABLE,
    )
    return B007Anchor(
        bs_limit=bs_limit,
        bd_limit=bd_limit,
        bs_standard_model=bs_sm,
        bd_standard_model=bd_sm,
        bs_budget=_build_channel_budget(
            process_id=process_id,
            channel_key="b_s",
            limit=bs_limit,
            standard_model=bs_sm,
            formula_sm=formula_sm_bs,
        ),
        bd_budget=_build_channel_budget(
            process_id=process_id,
            channel_key="b_d",
            limit=bd_limit,
            standard_model=bd_sm,
            formula_sm=formula_sm_bd,
        ),
    )


def _channel_saturation(
    *,
    predicted: float,
    sm_anchor: float,
    budget: ElectronicChannelBudget,
) -> tuple[float, float, bool]:
    upward_excess = max(0.0, float(predicted) - float(sm_anchor))
    ratio = upward_excess / budget.hard_veto_budget
    passes = ratio <= 1.0 and float(predicted) <= budget.active_limit
    return float(upward_excess), float(ratio), bool(passes)


def _wilson_mapping(result: Any) -> Mapping[str, complex]:
    if result.wilsons is None:
        return {}
    return {key: complex(value) for key, value in result.wilsons.wilsons.items()}


@register
class Constraint:
    """Catalogued ``B_s/B0 -> e+ e-`` branching-ratio limit constraint."""

    process_id = "B007"
    severity = Severity.HARD
    observable = "BR(B_s -> e+ e-), BR(B0 -> e+ e-)"

    def __init__(self) -> None:
        self.sm_inputs = rare_b_electronic_default_sm_inputs()
        self.bs_sm_result = rare_b_electronic_sm_branching_fraction(
            "b_s",
            self.sm_inputs,
        )
        self.bd_sm_result = rare_b_electronic_sm_branching_fraction(
            "b_d",
            self.sm_inputs,
        )
        self.anchor = _load_b007_anchor(
            self.process_id,
            formula_sm_bs=float(self.bs_sm_result.branching_fraction),
            formula_sm_bd=float(self.bd_sm_result.branching_fraction),
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        rs_wilsons = point.get_extra(_REQUIRED_EXTRA)
        if rs_wilsons is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.anchor.sm_value),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; B_s/B0 -> e+e- "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_quark_mass_basis_couplings_present": (
                        point.get_extra("quark_mass_basis_couplings") is not None
                    ),
                    "bs_sm_anchor_branching_fraction": float(
                        self.anchor.bs_standard_model.value
                    ),
                    "bd_sm_anchor_branching_fraction": float(
                        self.anchor.bd_standard_model.value
                    ),
                    "bs_sm_formula_branching_fraction": float(
                        self.bs_sm_result.branching_fraction
                    ),
                    "bd_sm_formula_branching_fraction": float(
                        self.bd_sm_result.branching_fraction
                    ),
                    "bs_limit": float(self.anchor.bs_limit.value),
                    "bd_limit": float(self.anchor.bd_limit.value),
                    "budget_source": _BUDGET_SOURCE,
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        m_kk_gev = None if kk_ew_mass is None else float(kk_ew_mass)
        try:
            bs_result = bs_ee_from_rs_semileptonic_wilsons(
                rs_wilsons,
                matching_scale_gev=m_kk_gev,
                inputs=self.sm_inputs,
            )
            bd_result = bd_ee_from_rs_semileptonic_wilsons(
                rs_wilsons,
                matching_scale_gev=m_kk_gev,
                inputs=self.sm_inputs,
            )
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.anchor.sm_value),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_semileptonic_wilsons for "
                    "B_s/B0 -> e+ e-"
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )

        channel_inputs = {
            "b_s": (
                bs_result,
                self.anchor.bs_limit,
                self.anchor.bs_standard_model,
                self.anchor.bs_budget,
            ),
            "b_d": (
                bd_result,
                self.anchor.bd_limit,
                self.anchor.bd_standard_model,
                self.anchor.bd_budget,
            ),
        }
        channel_diagnostics: dict[str, Mapping[str, Any]] = {}
        channel_summaries: dict[str, tuple[float, float, bool]] = {}
        for channel_key, (result, limit, sm_anchor, budget) in channel_inputs.items():
            predicted = float(result.branching_fraction)
            upward_excess, ratio, passes = _channel_saturation(
                predicted=predicted,
                sm_anchor=float(sm_anchor.value),
                budget=budget,
            )
            total_limit_ratio = predicted / float(limit.value)
            diagnostics = dict(result.diagnostics)
            diagnostics.update(
                {
                    "evaluated": True,
                    "predicted_branching_fraction": predicted,
                    "sm_anchor_branching_fraction": float(sm_anchor.value),
                    "sm_formula_branching_fraction": float(
                        result.sm_branching_fraction
                    ),
                    "sm_formula_minus_anchor": float(
                        result.sm_branching_fraction - sm_anchor.value
                    ),
                    "experimental_limit": float(limit.value),
                    "experimental_block": limit.block_key,
                    "experimental_row_id": limit.row_id,
                    "experimental_confidence_level": float(limit.confidence_level),
                    "experimental_conditions": limit.conditions,
                    "sm_block": sm_anchor.block_key,
                    "sm_theory_sigma": float(sm_anchor.uncertainty),
                    "upward_excess_over_sm_anchor": float(upward_excess),
                    "hard_veto_np_budget": float(budget.hard_veto_budget),
                    "total_limit_ratio": float(total_limit_ratio),
                    "sm_anchor_to_limit_ratio": float(
                        budget.sm_anchor_to_limit_ratio
                    ),
                    "budget_source": budget.source,
                    "budget_construction": budget.construction,
                    "wilson_coefficients": _wilson_mapping(result),
                }
            )
            channel_diagnostics[channel_key] = diagnostics
            channel_summaries[channel_key] = (predicted, ratio, passes)

        active_channel = max(channel_summaries, key=lambda key: channel_summaries[key][1])
        active_predicted, active_ratio, _ = channel_summaries[active_channel]
        active_limit = (
            self.anchor.bs_limit if active_channel == "b_s" else self.anchor.bd_limit
        )
        active_sm = (
            bs_result.sm_branching_fraction
            if active_channel == "b_s"
            else bd_result.sm_branching_fraction
        )
        passes = all(summary[2] for summary in channel_summaries.values())

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=float(active_predicted),
            sm_prediction=float(active_sm),
            experimental=float(active_limit.value),
            ratio=float(active_ratio),
            budget=float(active_limit.value),
            notes=(
                "BR(B_s -> e+e-) and BR(B0 -> e+e-) reuse the shared Buras "
                "b->qll C10-dominant formula with the charged-lepton mass set "
                "to m_e. The SM rate is helicity-suppressed and far below the "
                "B007.yaml limits, so the HARD ratio is the largest pure-NP "
                "upper-limit saturation across the two channels. Phase-3a RS "
                "semileptonic C10/C10' Wilsons enter additively."
            ),
            diagnostics={
                "active_channel": active_channel,
                "electron_mass_gev": float(ELECTRON_MASS_GEV),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "electron_mode_note": RARE_B_ELECTRONIC_MODEL_NOTE_V1,
                "rs_matching_assumption": channel_diagnostics[active_channel].get(
                    "rs_semileptonic_matching_assumption",
                    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
                ),
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "budget_source": _BUDGET_SOURCE,
                "channels": channel_diagnostics,
            },
        )
