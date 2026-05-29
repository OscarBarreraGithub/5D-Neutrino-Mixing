"""B011 - inclusive radiative decay ``Bbar -> X_s gamma``.

Physics
-------
The SM inclusive rate is normalized to the B011 YAML NNLO anchor at
``E_gamma > 1.6 GeV`` and evaluated with the shared C7 dipole machinery,

    BR = BR_SM * (|C7_SM + C7_NP|^2 + |C7p_SM + C7p_NP|^2)
              / (|C7_SM|^2 + |C7p_SM|^2).

The low-level C7 formula lives in ``quarkConstraints.bsgamma`` and is reached
only through the ``flavor_catalog_constraints.physics_adapters.bsgamma``
boundary so B012 can reuse the same dipole coefficient machinery.

Severity
--------
HARD.  The inclusive branching fraction is observed and SM-predicted.  The
HARD veto uses the YAML-measured-vs-SM room for the NP shift,

    |BR_total - BR_SM| <= |BR_exp - BR_SM| + sqrt(sigma_exp^2 + sigma_SM^2).

RS matching is a documented C7 coefficient proxy, not a full RS loop match, and
is explicitly flagged NEEDS-HUMAN-PHYSICS in the result diagnostics.  The
proxy dipoles are evolved from ``M_KK`` to ``mu_b`` with the shared
``b -> s gamma`` leading-log ``C7``-``C8`` running before they are combined
with the low-scale SM ``C7_SM(mu_b)``.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B011.yaml`` is the source of truth for the
HFLAV inclusive average and Misiak-Rehman-Steinhauser SM prediction.  Numeric
branching-fraction values below are loaded from that sidecar, not hardcoded in
this constraint.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.bsgamma import (
    BSGAMMA_RS_MATCHING_ASSUMPTION_V1,
    bsgamma_default_sm_inputs,
    inclusive_bsgamma_from_couplings,
    inclusive_bsgamma_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_DIPOLE_MASS_EXTRA = "kk_ew_mass_gev"
_EXPERIMENTAL_ANCHOR_CANDIDATES = ("canonical_experimental_average",)
_SM_ANCHOR_CANDIDATES = ("sm_prediction_current",)
_SM_VALIDATION_CANDIDATES = ("sm_prediction_pdg_review_quote",)
_EXPECTED_UNITS = "branching fraction"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B011.yaml "
    "canonical_experimental_average + sm_prediction_current"
)
_PARAMETRIZATION_CITATION = (
    "Misiak-Rehman-Steinhauser arXiv:2002.01548 for the SM BR anchor; "
    "C7-normalized inclusive-rate proxy in quarkConstraints.bsgamma"
)


@dataclass(frozen=True)
class BsgammaBudgetBand:
    """Uncertainty-aware B011 NP-shift budget for the HARD veto."""

    source: str
    central_residual: float
    experimental_sigma: float
    sm_theory_sigma: float
    combined_sigma: float
    hard_veto_budget: float
    lower_total_edge: float
    upper_total_edge: float


@dataclass(frozen=True)
class B011Anchor:
    """Typed B011 anchor: HFLAV experiment, SM references, and budget."""

    experimental: Anchor
    standard_model: Anchor
    validation_standard_model: Anchor
    budget_band: BsgammaBudgetBand

    @property
    def value(self) -> float:
        return self.experimental.value

    @property
    def uncertainty(self) -> float | None:
        return self.experimental.uncertainty

    @property
    def source_url(self) -> str | None:
        return self.experimental.source_url

    @property
    def sm_value(self) -> float:
        return self.standard_model.value

    @property
    def budget(self) -> float:
        return self.budget_band.hard_veto_budget


def _validate_branching_anchor(anchor: Anchor, *, process_id: str, label: str) -> None:
    if anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: {label} anchor {anchor.block_key!r} must use units "
            f"{_EXPECTED_UNITS!r}, got {anchor.units!r}"
        )
    if anchor.value <= 0.0 or not math.isfinite(anchor.value):
        raise AnchorError(
            f"{process_id}: {label} anchor {anchor.block_key!r} must have a "
            "positive finite value"
        )
    if anchor.uncertainty is None or anchor.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: {label} anchor {anchor.block_key!r} must provide "
            "a positive uncertainty"
        )


def _build_budget_band(
    *,
    experimental: Anchor,
    standard_model: Anchor,
) -> BsgammaBudgetBand:
    _validate_branching_anchor(experimental, process_id="B011", label="experimental")
    _validate_branching_anchor(standard_model, process_id="B011", label="SM")
    exp_sigma = float(experimental.uncertainty)
    sm_sigma = float(standard_model.uncertainty)
    combined_sigma = math.sqrt(exp_sigma * exp_sigma + sm_sigma * sm_sigma)
    central_residual = abs(experimental.value - standard_model.value)
    hard_budget = central_residual + combined_sigma
    if hard_budget <= 0.0 or not math.isfinite(hard_budget):
        raise AnchorError("B011: constructed b -> s gamma budget is invalid")
    return BsgammaBudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(central_residual),
        experimental_sigma=exp_sigma,
        sm_theory_sigma=sm_sigma,
        combined_sigma=float(combined_sigma),
        hard_veto_budget=float(hard_budget),
        lower_total_edge=float(standard_model.value - hard_budget),
        upper_total_edge=float(standard_model.value + hard_budget),
    )


def _load_b011_anchor(process_id: str) -> B011Anchor:
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
    validation_standard_model = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_SM_VALIDATION_CANDIDATES,
    )
    _validate_branching_anchor(
        validation_standard_model,
        process_id=process_id,
        label="validation SM",
    )
    return B011Anchor(
        experimental=experimental,
        standard_model=standard_model,
        validation_standard_model=validation_standard_model,
        budget_band=_build_budget_band(
            experimental=experimental,
            standard_model=standard_model,
        ),
    )


def _budget_result(predicted: float, anchor: B011Anchor) -> tuple[float, float, float, bool]:
    np_shift = float(predicted - anchor.sm_value)
    budget = float(anchor.budget)
    ratio = abs(np_shift) / budget if budget > 0.0 else float("inf")
    return np_shift, budget, float(ratio), bool(ratio <= 1.0)


def _complex_mapping(mapping: Any) -> dict[str, complex]:
    return {str(key): complex(value) for key, value in dict(mapping).items()}


@register
class Constraint:
    """Catalogued inclusive ``Bbar -> X_s gamma`` branching-fraction constraint."""

    process_id = "B011"
    severity = Severity.HARD
    observable = "BR(Bbar -> X_s gamma)"

    def __init__(self) -> None:
        self.anchor = _load_b011_anchor(self.process_id)
        self.sm_inputs = bsgamma_default_sm_inputs()
        self.sm_result = inclusive_bsgamma_sm_branching_fraction(
            sm_branching_fraction=self.anchor.sm_value,
            inputs=self.sm_inputs,
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
                    f"extra {_REQUIRED_EXTRA!r} absent; Bbar -> X_s gamma "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "sm_anchor_branching_fraction": float(self.anchor.sm_value),
                    "sm_formula_branching_fraction": float(
                        self.sm_result.branching_fraction
                    ),
                    "budget_source": self.anchor.budget_band.source,
                },
            )

        kk_mass = point.get_extra(_OPTIONAL_DIPOLE_MASS_EXTRA)
        try:
            result = inclusive_bsgamma_from_couplings(
                couplings,
                sm_branching_fraction=self.anchor.sm_value,
                m_kk_gev=None if kk_mass is None else float(kk_mass),
                inputs=self.sm_inputs,
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
                    "NOT EVALUATED - invalid quark_mass_basis_couplings for "
                    "the Bbar -> X_s gamma C7 proxy"
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "budget_source": self.anchor.budget_band.source,
                },
            )

        predicted = float(result.branching_fraction)
        np_shift, budget, ratio, passes = _budget_result(predicted, self.anchor)
        diagnostics = dict(result.diagnostics)
        wilsons = result.wilsons.wilsons if result.wilsons is not None else {}
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
                "validation_sm_pdg_review_branching_fraction": float(
                    self.anchor.validation_standard_model.value
                ),
                "experimental_block": self.anchor.experimental.block_key,
                "sm_block": self.anchor.standard_model.block_key,
                "validation_sm_block": self.anchor.validation_standard_model.block_key,
                "np_shift_branching_fraction": float(np_shift),
                "ratio_to_sm_c7_power": float(result.ratio_to_sm),
                "c7_sm_eff": complex(result.c7_sm_eff),
                "c7p_sm_eff": complex(result.c7p_sm_eff),
                "c7_total": complex(result.c7_total),
                "c7p_total": complex(result.c7p_total),
                "c7_np": complex(result.c7_np),
                "c7p_np": complex(result.c7p_np),
                "wilson_coefficients": _complex_mapping(wilsons),
                "experimental_sigma": float(self.anchor.budget_band.experimental_sigma),
                "sm_theory_sigma": float(self.anchor.budget_band.sm_theory_sigma),
                "budget_combined_sigma": float(self.anchor.budget_band.combined_sigma),
                "budget_central_residual": float(
                    self.anchor.budget_band.central_residual
                ),
                "hard_veto_np_shift_budget": float(
                    self.anchor.budget_band.hard_veto_budget
                ),
                "budget_lower_total_edge": float(
                    self.anchor.budget_band.lower_total_edge
                ),
                "budget_upper_total_edge": float(
                    self.anchor.budget_band.upper_total_edge
                ),
                "budget_source": self.anchor.budget_band.source,
                "budget_policy": (
                    "|BR_total - BR_SM| compared with "
                    "|BR_exp - BR_SM| + sqrt(sigma_exp^2 + sigma_SM^2)"
                ),
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": BSGAMMA_RS_MATCHING_ASSUMPTION_V1,
                "needs_human_physics": (
                    "NEEDS-HUMAN-PHYSICS: full RS loop matching for C7/C8 "
                    "requires KK fermion, Higgs/Goldstone, charged-current, "
                    "and finite matching inputs not available on ParameterPoint; "
                    "v1 uses documented b-s overlap C7/C8 proxies with LL "
                    "QCD running to mu_b."
                ),
                "kk_ew_mass_extra_used": kk_mass is not None,
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
                "BR(Bbar -> X_s gamma) uses a C7-normalized inclusive-rate "
                "formula with the B011 YAML SM anchor. The RS contribution is "
                "a documented b-s overlap proxy for C7/C7p and is marked "
                "NEEDS-HUMAN-PHYSICS; the HARD budget is the measured-vs-SM "
                "NP-shift room from B011.yaml."
            ),
            diagnostics=diagnostics,
        )
