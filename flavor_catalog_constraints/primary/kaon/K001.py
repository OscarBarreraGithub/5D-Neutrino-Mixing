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
``|epsilon_exp - epsilon_SM|`` is kept as a diagnostic; the HARD veto uses the
loose edge of the documented 1-sigma band.

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
_SM_CHOICE_SENSITIVITY = 0.15e-3
_MU_HAD_GEV = 2.0


@dataclass(frozen=True)
class EpsilonKBudgetBand:
    """Documented epsilon_K NP budget band for the K001 HARD veto.

    The construction follows ``docs/audits/epsilon_k_sm_decision.md``:
    combine the BGS grouped SM uncertainty, PDG experimental uncertainty, and
    SM-choice sensitivity in quadrature; keep the central residual as a
    diagnostic; use the loose edge as the veto budget for catalog exclusions.
    The BGS grouped uncertainty includes the published kappa_epsilon=0.94(2)
    contribution through the non-perturbative group.
    """

    doc_citation: str
    central_budget: float
    tight_budget: float
    loose_budget: float
    hard_veto_budget: float
    sm_theory_sigma: float
    experimental_sigma: float
    sm_choice_sigma: float
    combined_sigma: float
    sm_at_loose_edge: float
    sm_at_tight_edge: float
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
        """Uncertainty-aware HARD veto budget for ``|epsilon_K^NP|``."""
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
    if experimental.uncertainty is None or experimental.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: epsilon_K experimental uncertainty is required "
            "for the documented budget band"
        )

    grouped = _load_bgs_grouped_uncertainties(
        standard_model_sub,
        process_id=process_id,
    )
    sm_theory_sigma = math.sqrt(sum(value * value for value in grouped.values()))
    combined_sigma = math.sqrt(
        sm_theory_sigma * sm_theory_sigma
        + experimental.uncertainty * experimental.uncertainty
        + _SM_CHOICE_SENSITIVITY * _SM_CHOICE_SENSITIVITY
    )
    central_budget = abs(experimental.value - standard_model.value)
    sm_at_loose_edge = standard_model.value - combined_sigma
    sm_at_tight_edge = standard_model.value + combined_sigma
    loose_budget = abs(experimental.value - sm_at_loose_edge)

    # The audit prescribes the lower/tighter band edge as the experimental
    # error floor once the SM-shifted central value crosses the experiment.
    tight_budget = float(experimental.uncertainty)

    if central_budget <= 0.0 or loose_budget <= 0.0 or tight_budget <= 0.0:
        raise AnchorError(
            f"{process_id}: epsilon_K NP budgets must be positive "
            f"(central={central_budget}, loose={loose_budget}, tight={tight_budget})"
        )

    return EpsilonKBudgetBand(
        doc_citation=_BUDGET_DOC_CITATION,
        central_budget=central_budget,
        tight_budget=tight_budget,
        loose_budget=loose_budget,
        hard_veto_budget=loose_budget,
        sm_theory_sigma=sm_theory_sigma,
        experimental_sigma=float(experimental.uncertainty),
        sm_choice_sigma=_SM_CHOICE_SENSITIVITY,
        combined_sigma=combined_sigma,
        sm_at_loose_edge=sm_at_loose_edge,
        sm_at_tight_edge=sm_at_tight_edge,
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
                    "was not evaluated."
                ),
                diagnostics={"missing_extra": _REQUIRED_EXTRA},
            )

        wilsons = epsilon_k_wilsons_from_couplings(couplings)
        result = epsilon_k_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
            epsilon_k_np_budget=self.anchor.budget,
        )

        predicted = float(result.epsilon_k_np)
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
                "|epsilon_K^NP| uses kappa_epsilon Im(M12^NP)/(sqrt(2) Delta m_K); "
                "Wilsons are QCD-evolved to 2 GeV; HARD budget is the "
                "uncertainty-aware loose band edge from "
                f"{_BUDGET_DOC_CITATION}."
            ),
            diagnostics={
                "im_m12_np_gev": float(result.im_m12_np),
                "epsilon_k_np_is_absolute": True,
                "qcd_running_applied": True,
                "hadronic_scale_gev": _MU_HAD_GEV,
                "matching_scale_gev": float(wilsons.matching_scale),
                "m_kk_gev": float(wilsons.M_KK),
                "left_sd_coupling": complex(wilsons.left_coupling),
                "right_sd_coupling": complex(wilsons.right_coupling),
                "wilson_coefficients": _complex_mapping(wilsons.wilsons),
                "budget_doc_citation": self.anchor.budget_band.doc_citation,
                "central_np_budget": self.anchor.budget_band.central_budget,
                "tight_band_np_budget": self.anchor.budget_band.tight_budget,
                "loose_band_np_budget": self.anchor.budget_band.loose_budget,
                "hard_veto_np_budget": self.anchor.budget_band.hard_veto_budget,
                "budget_combined_sigma": self.anchor.budget_band.combined_sigma,
                "budget_sm_theory_sigma": self.anchor.budget_band.sm_theory_sigma,
                "budget_experimental_sigma": (
                    self.anchor.budget_band.experimental_sigma
                ),
                "budget_sm_choice_sigma": self.anchor.budget_band.sm_choice_sigma,
                "budget_sm_at_loose_edge": self.anchor.budget_band.sm_at_loose_edge,
                "budget_sm_at_tight_edge": self.anchor.budget_band.sm_at_tight_edge,
                "bgs_grouped_uncertainties": dict(
                    self.anchor.budget_band.bgs_grouped_uncertainties
                ),
                "kappa_epsilon_uncertainty_policy": (
                    "BGS kappa_epsilon=0.94(2) is folded into the grouped "
                    "SM uncertainty used by the documented budget band."
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
