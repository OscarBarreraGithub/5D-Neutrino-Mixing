"""K001 - indirect CP violation in neutral-kaon mixing.

Physics
-------
``epsilon_K`` measures indirect CP violation in ``K0-K0bar`` mixing.  The
new-physics contribution is evaluated with the audited Delta F = 2 core,

    epsilon_K^NP = kappa_epsilon Im(M12^NP) / (sqrt(2) Delta m_K),

where ``M12^NP`` is built from the kaon Wilson coefficients produced by
``quarkConstraints.deltaf2``.  This module reaches that core only through the
``flavor_catalog_constraints.physics_adapters.deltaf2`` boundary.

Severity
--------
HARD.  ``epsilon_K`` is observed, and any RS new-physics contribution must fit
inside the uncertainty-aware NP budget band documented in
``docs/audits/epsilon_k_sm_decision.md:34-37,71-100``.  The central residual
``epsilon_exp - epsilon_SM`` is kept as a diagnostic; the HARD veto uses the
shared Delta-F=2 direction-aware one-sigma budget.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K001.yaml`` is the source of truth for the
experimental ``|epsilon|`` value, the BGS2020 SM reference, and FLAG bag-input
provenance.  The numbers below are loaded through the scaffold anchor loader,
not hardcoded in this constraint.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    find_block,
    load_anchor,
    load_pdg_block,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.deltaf2 import (
    EpsilonKBudgetPolicy,
    epsilon_k_budget_policy,
    epsilon_k_from_wilsons_with_running,
    epsilon_k_wilsons_from_couplings,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"

_EXPERIMENTAL_ANCHOR_CANDIDATES = ("canonical_experimental_value",)
_SM_ANCHOR_CANDIDATES = ("standard_model_reference",)
_FLAG_BAG_CANDIDATES = ("flag_bag_parameters",)
_BUDGET_DOC_CITATION = "docs/audits/epsilon_k_sm_decision.md:34-37,71-100"
_MU_HAD_GEV = 3.0


@dataclass(frozen=True)
class EpsilonKBudgetBand:
    """Shared epsilon_K NP budget policy plus K001 YAML provenance.

    The numbers are owned by ``quarkConstraints.deltaf2``.  K001 only validates
    that its YAML anchor agrees with the shared policy and carries the grouped
    BGS uncertainty components as provenance diagnostics.
    """

    policy_id: str
    confidence_level: str
    doc_citation: str
    central_budget: float
    tight_budget: float
    loose_budget: float
    hard_veto_budget: float
    budget_lowers_epsilon_k: float
    budget_raises_epsilon_k: float
    signed_lower_edge: float
    signed_upper_edge: float
    sm_theory_sigma: float
    experimental_sigma: float
    combined_sigma: float
    sm_choice_sensitivity: float
    bgs_grouped_uncertainties: Mapping[str, float]


@dataclass(frozen=True)
class FlagBagInputs:
    """Typed provenance view over the FLAG bag-parameter block."""

    block_key: str
    source: str | None
    year: int | None
    preferred_theory: str | None
    hat_b_k_nf21: str | None
    b_k_msbar_2gev_nf21: str | None
    b_k_msbar_3gev_nf21: str | None
    hat_b_k_nf211: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class EpsilonKAnchor:
    """Typed K001 anchor: experiment, SM reference, FLAG provenance."""

    experimental: Anchor
    standard_model: Anchor
    flag_bag_parameters: FlagBagInputs
    budget_band: EpsilonKBudgetBand

    @property
    def value(self) -> float:
        """Experimental central value, kept for common Anchor-like access."""
        return self.experimental.value

    @property
    def uncertainty(self) -> float | None:
        """Experimental uncertainty, kept for common Anchor-like access."""
        return self.experimental.uncertainty

    @property
    def source_url(self) -> str | None:
        """Experimental source URL, kept for common Anchor-like access."""
        return self.experimental.source_url

    @property
    def sm_value(self) -> float:
        """BGS2020 SM central value used to construct the NP budget band."""
        return self.standard_model.value

    @property
    def central_budget(self) -> float:
        """Central residual ``|epsilon_K^exp - epsilon_K^SM|``."""
        return self.budget_band.central_budget

    @property
    def budget(self) -> float:
        """Largest scalar edge for legacy display; actual veto is sign-selected."""
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
        return float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: epsilon_K budget field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc


def _load_bgs_grouped_uncertainties(
    sub: Mapping[str, Any],
    *,
    process_id: str,
) -> dict[str, float]:
    grouped = sub.get("grouped_uncertainties")
    if not isinstance(grouped, Mapping):
        raise AnchorError(
            f"{process_id}: standard_model_reference must provide "
            "'grouped_uncertainties' for the documented epsilon_K budget band"
        )
    required = ("parametric", "non_perturbative", "perturbative")
    return {
        key: _required_float(
            grouped.get(key),
            process_id=process_id,
            field_name=f"grouped_uncertainties.{key}",
        )
        for key in required
    }


def _build_budget_band(
    *,
    process_id: str,
    experimental: Anchor,
    standard_model: Anchor,
    standard_model_sub: Mapping[str, Any],
) -> EpsilonKBudgetBand:
    policy: EpsilonKBudgetPolicy = epsilon_k_budget_policy()
    if experimental.uncertainty is None or experimental.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: epsilon_K experimental uncertainty is required "
            "for the documented budget band"
        )
    if not math.isclose(
        experimental.value,
        policy.experimental_value,
        rel_tol=0.0,
        abs_tol=1e-15,
    ):
        raise AnchorError(
            f"{process_id}: experimental epsilon_K value {experimental.value} "
            f"does not match shared policy {policy.experimental_value}"
        )
    if not math.isclose(
        experimental.uncertainty,
        policy.sigma_exp,
        rel_tol=0.0,
        abs_tol=1e-15,
    ):
        raise AnchorError(
            f"{process_id}: experimental epsilon_K uncertainty "
            f"{experimental.uncertainty} does not match shared policy {policy.sigma_exp}"
        )
    if not math.isclose(
        standard_model.value,
        policy.sm_value,
        rel_tol=0.0,
        abs_tol=1e-15,
    ):
        raise AnchorError(
            f"{process_id}: SM epsilon_K value {standard_model.value} "
            f"does not match shared policy {policy.sm_value}"
        )

    grouped = _load_bgs_grouped_uncertainties(
        standard_model_sub,
        process_id=process_id,
    )
    if (
        policy.central_budget <= 0.0
        or policy.tight_budget <= 0.0
        or policy.loose_budget <= 0.0
    ):
        raise AnchorError(
            f"{process_id}: epsilon_K NP budgets must be positive "
            f"(central={policy.central_budget}, loose={policy.loose_budget}, "
            f"tight={policy.tight_budget})"
        )

    return EpsilonKBudgetBand(
        policy_id=policy.policy_id,
        confidence_level=policy.confidence_level,
        doc_citation=policy.doc_citation,
        central_budget=policy.central_budget,
        tight_budget=policy.tight_budget,
        loose_budget=policy.loose_budget,
        hard_veto_budget=policy.loose_budget,
        budget_lowers_epsilon_k=policy.budget_lowers_epsilon_k,
        budget_raises_epsilon_k=policy.budget_raises_epsilon_k,
        signed_lower_edge=policy.lower_signed_edge,
        signed_upper_edge=policy.upper_signed_edge,
        sm_theory_sigma=policy.sigma_bgs,
        experimental_sigma=policy.sigma_exp,
        combined_sigma=policy.primary_combined_sigma,
        sm_choice_sensitivity=policy.sm_choice_sensitivity,
        bgs_grouped_uncertainties=grouped,
    )


def _load_flag_bag_inputs(process_id: str) -> FlagBagInputs:
    pdg_block = load_pdg_block(process_id, family=_FAMILY)
    sub = find_block(pdg_block, _FLAG_BAG_CANDIDATES, process_id=process_id)
    block_key = next(key for key in _FLAG_BAG_CANDIDATES if key in pdg_block)
    return FlagBagInputs(
        block_key=block_key,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        preferred_theory=_optional_str(sub.get("preferred_theory")),
        hat_b_k_nf21=_optional_str(sub.get("hat_B_K_Nf21")),
        b_k_msbar_2gev_nf21=_optional_str(sub.get("B_K_MSbar_2GeV_Nf21")),
        b_k_msbar_3gev_nf21=_optional_str(sub.get("B_K_MSbar_3GeV_Nf21")),
        hat_b_k_nf211=_optional_str(sub.get("hat_B_K_Nf211")),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _load_epsilon_k_anchor(process_id: str) -> EpsilonKAnchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_EXPERIMENTAL_ANCHOR_CANDIDATES,
    )
    standard_model = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_SM_ANCHOR_CANDIDATES,
    )
    pdg_block = load_pdg_block(process_id, family=_FAMILY)
    standard_model_sub = find_block(
        pdg_block,
        _SM_ANCHOR_CANDIDATES,
        process_id=process_id,
    )
    anchor = EpsilonKAnchor(
        experimental=experimental,
        standard_model=standard_model,
        flag_bag_parameters=_load_flag_bag_inputs(process_id),
        budget_band=_build_budget_band(
            process_id=process_id,
            experimental=experimental,
            standard_model=standard_model,
            standard_model_sub=standard_model_sub,
        ),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(
            f"{process_id}: epsilon_K NP budget must be positive, got {anchor.budget}"
        )
    return anchor


def _complex_mapping(values: Mapping[str, complex]) -> dict[str, complex]:
    return {str(key): complex(value) for key, value in values.items()}


@register
class Constraint:
    """Catalogued epsilon_K Delta F=2 constraint (process_id K001)."""

    process_id = "K001"
    severity = Severity.HARD
    observable = "epsilon_K"

    def __init__(self) -> None:
        self.anchor = _load_epsilon_k_anchor(self.process_id)
        self.sm_value = self.anchor.sm_value

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.anchor.sm_value),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; epsilon_K constraint "
                    "was not evaluated. NOTE: catalog-wide single-CL migration "
                    "is a separate policy decision."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "budget_policy_id": self.anchor.budget_band.policy_id,
                    "confidence_level": self.anchor.budget_band.confidence_level,
                    "catalog_confidence_level_note": (
                        "Delta-F=2 epsilon_K reports its own one-sigma "
                        "sensitivity label; no global catalog CL is implied."
                    ),
                },
            )

        wilsons = epsilon_k_wilsons_from_couplings(couplings)
        result = epsilon_k_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
        )

        predicted = float(result.epsilon_k_np_signed)
        ratio = float(result.ratio_to_budget)
        budget = float(result.epsilon_k_np_budget)

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(result.passes),
            predicted=predicted,
            sm_prediction=float(self.anchor.sm_value),
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "Signed epsilon_K^NP uses kappa_epsilon Im(M12^NP)/(sqrt(2) Delta m_K); "
                "Wilsons are QCD-evolved to 3 GeV; HARD budget is selected by "
                "the sign of epsilon_K^NP from the shared one-sigma policy. "
                "NOTE: catalog-wide single-CL migration is a separate policy decision."
            ),
            diagnostics={
                "im_m12_np_gev": float(result.im_m12_np),
                "epsilon_k_np_signed": float(result.epsilon_k_np_signed),
                "epsilon_k_np_abs": float(result.epsilon_k_np),
                "epsilon_k_np_is_absolute": False,
                "epsilon_k_selected_budget_direction": result.selected_budget_direction,
                "epsilon_k_selected_signed_budget": float(result.epsilon_k_np_budget),
                "qcd_running_applied": True,
                "hadronic_scale_gev": _MU_HAD_GEV,
                "running_order": result.running_order,
                "running_bias_note": result.running_bias_note,
                "matching_scale_gev": float(wilsons.matching_scale),
                "m_kk_gev": float(wilsons.M_KK),
                "m_kk_physical_gev": float(
                    wilsons.m_kk_physical_gev
                    if wilsons.m_kk_physical_gev is not None
                    else wilsons.M_KK
                ),
                "lambda_ir_gev": (
                    None if wilsons.lambda_ir_gev is None else float(wilsons.lambda_ir_gev)
                ),
                "mass_convention_id": wilsons.mass_convention_id,
                "coupling_policy_id": wilsons.coupling_policy_id,
                "operator_convention_id": wilsons.operator_convention_id,
                "g_s_4d": None if wilsons.g_s_4d is None else float(wilsons.g_s_4d),
                "g_eff": None if wilsons.g_eff is None else float(wilsons.g_eff),
                "g_s_multiplier": (
                    None if wilsons.g_s_multiplier is None else float(wilsons.g_s_multiplier)
                ),
                "left_sd_coupling": complex(wilsons.left_coupling),
                "right_sd_coupling": complex(wilsons.right_coupling),
                "wilson_coefficients": _complex_mapping(wilsons.wilsons),
                "budget_policy_id": self.anchor.budget_band.policy_id,
                "confidence_level": self.anchor.budget_band.confidence_level,
                "budget_doc_citation": self.anchor.budget_band.doc_citation,
                "central_np_budget": self.anchor.budget_band.central_budget,
                "tight_band_np_budget": self.anchor.budget_band.tight_budget,
                "loose_band_np_budget": self.anchor.budget_band.loose_budget,
                "hard_veto_np_budget": float(result.epsilon_k_np_budget),
                "budget_lowers_epsilon_k": (
                    self.anchor.budget_band.budget_lowers_epsilon_k
                ),
                "budget_raises_epsilon_k": (
                    self.anchor.budget_band.budget_raises_epsilon_k
                ),
                "signed_lower_edge": self.anchor.budget_band.signed_lower_edge,
                "signed_upper_edge": self.anchor.budget_band.signed_upper_edge,
                "budget_combined_sigma": self.anchor.budget_band.combined_sigma,
                "budget_primary_combined_sigma": (
                    self.anchor.budget_band.combined_sigma
                ),
                "budget_sm_theory_sigma": self.anchor.budget_band.sm_theory_sigma,
                "budget_experimental_sigma": (
                    self.anchor.budget_band.experimental_sigma
                ),
                "sm_choice_sensitivity": (
                    self.anchor.budget_band.sm_choice_sensitivity
                ),
                "budget_sm_choice_sensitivity": (
                    self.anchor.budget_band.sm_choice_sensitivity
                ),
                "sm_choice_sensitivity_in_hard_gate": False,
                "catalog_confidence_level_note": (
                    "Delta-F=2 epsilon_K reports its own one-sigma "
                    "sensitivity label; no global catalog CL is implied."
                ),
                "bgs_grouped_uncertainties": dict(
                    self.anchor.budget_band.bgs_grouped_uncertainties
                ),
                "kappa_epsilon_uncertainty_policy": (
                    "BGS kappa_epsilon=0.94(2) is part of the published "
                    "sigma_BGS=0.18e-3 used by the shared hard budget policy."
                ),
                "experimental_block": self.anchor.experimental.block_key,
                "sm_block": self.anchor.standard_model.block_key,
                "flag_bag_block": self.anchor.flag_bag_parameters.block_key,
                "flag_b_k_msbar_2gev_nf21": (
                    self.anchor.flag_bag_parameters.b_k_msbar_2gev_nf21
                ),
                "flag_b_k_msbar_3gev_nf21": (
                    self.anchor.flag_bag_parameters.b_k_msbar_3gev_nf21
                ),
                "flag_hat_b_k_nf21": self.anchor.flag_bag_parameters.hat_b_k_nf21,
            },
        )
