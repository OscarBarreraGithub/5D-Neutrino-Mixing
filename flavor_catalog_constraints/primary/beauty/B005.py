"""B005 - rare leptonic decay ``B_s -> mu+ mu-``.

Physics
-------
``BR(B_s -> mu+ mu-)`` is evaluated with the shared ``b -> q l l`` core built
for this process.  The SM part uses the Buras effective-Hamiltonian
normalization with dominant ``C10``:

    BR(B_s -> mu mu) proportional to f_Bs^2 tau_Bs m_mu^2
    |V_tb V_ts^*|^2 |C10 - C10'|^2.

The time-integrated branching fraction uses ``A_DeltaGamma`` computed from the
complex scalar and pseudoscalar amplitudes, reducing to ``+1`` in the SM-like
pure-C10 limit.

The underlying core also computes a documented ``C9``/``C10`` RS proxy for
reuse by later ``b -> s l l`` siblings.  ``C9`` is diagnostic-only for this
pure leptonic mode.

Severity
--------
HARD.  The observed branching fraction is compared to the SM anchor by giving
new physics the uncertainty-aware room
``|BR_exp - BR_SM(anchor)| + sqrt(sigma_exp^2 + sigma_SM^2)``.  This mirrors
the existing Delta-F=2 loose-edge budget policy: the SM point is not vetoed by
the current one-sigma SM-vs-data offset, but a large RS shift is excluded.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B005.yaml`` is the source of truth for the
PDG live/API experimental average and the Czaja-Misiak 2024 SM prediction.
Numeric budget inputs below are loaded from that sidecar, not hardcoded here.

NEEDS-HUMAN-PHYSICS
-------------------
The Phase-3a RS light-Z contribution supplies rigorous vector/axial
``C9/C10/C9'/C10'`` Wilsons.  Higgs/radion scalar and pseudoscalar matching
remains deferred.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_b_meson import (
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    bs_mumu_from_rs_semileptonic_wilsons,
    rare_b_dilepton_default_sm_inputs,
    rare_b_dilepton_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "rs_semileptonic_wilsons"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPECTED_UNITS = "branching fraction"

_EXPERIMENTAL_ANCHOR_CANDIDATES = ("canonical_experimental_average",)
_HISTORICAL_ANCHOR_CANDIDATES = ("hflav_historical_average",)
_SM_ANCHOR_CANDIDATES = ("standard_model_prediction",)
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B005.yaml "
    "canonical_experimental_average + standard_model_prediction"
)
_PARAMETRIZATION_CITATION = (
    "Buras b->s l l effective Hamiltonian; "
    "Czaja-Misiak arXiv:2407.03810 SM validation"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: Phase-3a supplies the light-Z vector/axial "
    "C9/C10/C9'/C10' terms; Higgs/radion scalar and pseudoscalar "
    "b->s mu mu matching remains deferred."
)


@dataclass(frozen=True)
class BSMuMuBudgetBand:
    """Uncertainty-aware B005 NP-shift budget in branching-fraction units."""

    source: str
    central_residual: float
    experimental_sigma: float
    sm_theory_sigma: float
    combined_sigma: float
    hard_veto_budget: float
    experimental_lower_edge: float
    experimental_upper_edge: float


@dataclass(frozen=True)
class BSMuMuAnchor:
    """Typed B005 anchor: experiment, historical average, SM, and budget."""

    experimental: Anchor
    historical_average: Anchor
    standard_model: Anchor
    budget_band: BSMuMuBudgetBand

    @property
    def value(self) -> float:
        return self.experimental.value

    @property
    def uncertainty(self) -> float | None:
        return self.experimental.uncertainty

    @property
    def sm_value(self) -> float:
        return self.standard_model.value

    @property
    def budget(self) -> float:
        return self.budget_band.hard_veto_budget


def _validate_branching_anchor(anchor: Anchor, *, process_id: str, name: str) -> None:
    if anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: {name} anchor {anchor.block_key!r} must use units "
            f"{_EXPECTED_UNITS!r}, got {anchor.units!r}"
        )
    if anchor.value <= 0.0 or not math.isfinite(anchor.value):
        raise AnchorError(
            f"{process_id}: {name} branching fraction must be positive and finite"
        )
    if anchor.uncertainty is None or anchor.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: {name} uncertainty is required for the B005 budget"
        )


def _build_budget_band(
    *,
    process_id: str,
    experimental: Anchor,
    standard_model: Anchor,
) -> BSMuMuBudgetBand:
    _validate_branching_anchor(experimental, process_id=process_id, name="experimental")
    _validate_branching_anchor(standard_model, process_id=process_id, name="standard_model")
    exp_sigma = float(experimental.uncertainty)
    sm_sigma = float(standard_model.uncertainty)
    combined = math.sqrt(exp_sigma * exp_sigma + sm_sigma * sm_sigma)
    central = abs(float(experimental.value) - float(standard_model.value))
    budget = central + combined
    if budget <= 0.0 or not math.isfinite(budget):
        raise AnchorError(f"{process_id}: B005 NP budget must be positive and finite")
    return BSMuMuBudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(central),
        experimental_sigma=exp_sigma,
        sm_theory_sigma=sm_sigma,
        combined_sigma=float(combined),
        hard_veto_budget=float(budget),
        experimental_lower_edge=float(experimental.value - combined),
        experimental_upper_edge=float(experimental.value + combined),
    )


def _load_bs_mumu_anchor(process_id: str) -> BSMuMuAnchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_EXPERIMENTAL_ANCHOR_CANDIDATES,
    )
    historical = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_HISTORICAL_ANCHOR_CANDIDATES,
    )
    standard_model = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_SM_ANCHOR_CANDIDATES,
    )
    _validate_branching_anchor(historical, process_id=process_id, name="historical")
    return BSMuMuAnchor(
        experimental=experimental,
        historical_average=historical,
        standard_model=standard_model,
        budget_band=_build_budget_band(
            process_id=process_id,
            experimental=experimental,
            standard_model=standard_model,
        ),
    )


@register
class Constraint:
    """Catalogued ``B_s -> mu+ mu-`` branching-ratio constraint (B005)."""

    process_id = "B005"
    severity = Severity.HARD
    observable = "BR(B_s -> mu+ mu-)"

    def __init__(self) -> None:
        self.anchor = _load_bs_mumu_anchor(self.process_id)
        self.sm_inputs = rare_b_dilepton_default_sm_inputs()
        self.sm_result = rare_b_dilepton_sm_branching_fraction("b_s", self.sm_inputs)

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
                    f"extra {_REQUIRED_EXTRA!r} absent; B_s -> mu+ mu- "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "legacy_quark_mass_basis_couplings_present": (
                        point.get_extra("quark_mass_basis_couplings") is not None
                    ),
                    "sm_anchor_branching_fraction": float(self.anchor.sm_value),
                    "sm_formula_branching_fraction": float(
                        self.sm_result.branching_fraction
                    ),
                    "sm_a_delta_gamma": float(self.sm_result.a_delta_gamma),
                    "sm_time_integrated_factor": float(
                        self.sm_result.time_integrated_factor
                    ),
                    "budget_source": self.anchor.budget_band.source,
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        try:
            result = bs_mumu_from_rs_semileptonic_wilsons(
                rs_wilsons,
                lepton="mu",
                matching_scale_gev=None if kk_ew_mass is None else float(kk_ew_mass),
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
                    "NOT EVALUATED - invalid rs_semileptonic_wilsons for "
                    "B_s -> mu+ mu-"
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )
        predicted = float(result.branching_fraction)
        np_shift = float(result.np_shift_branching_fraction)
        budget = float(self.anchor.budget)
        ratio = abs(np_shift) / budget if budget > 0.0 else float("inf")
        total_exp_pull = predicted - float(self.anchor.value)
        total_exp_pull_ratio = (
            abs(total_exp_pull) / self.anchor.budget_band.combined_sigma
        )

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
                "historical_hflav_branching_fraction": float(
                    self.anchor.historical_average.value
                ),
                "experimental_block": self.anchor.experimental.block_key,
                "sm_block": self.anchor.standard_model.block_key,
                "experimental_sigma": float(
                    self.anchor.budget_band.experimental_sigma
                ),
                "sm_theory_sigma": float(self.anchor.budget_band.sm_theory_sigma),
                "budget_combined_sigma": float(
                    self.anchor.budget_band.combined_sigma
                ),
                "budget_central_residual": float(
                    self.anchor.budget_band.central_residual
                ),
                "budget_source": self.anchor.budget_band.source,
                "total_minus_experiment": float(total_exp_pull),
                "total_minus_experiment_ratio_to_combined_sigma": float(
                    total_exp_pull_ratio
                ),
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": (
                    result.wilsons.matching_assumption
                    if result.wilsons is not None
                    else RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1
                ),
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
            passes=bool(ratio <= 1.0),
            predicted=predicted,
            sm_prediction=float(result.sm_branching_fraction),
            experimental=float(self.anchor.value),
            ratio=float(ratio),
            budget=budget,
            notes=(
                "BR(B_s -> mu+ mu-) uses the Buras b->sll C10-dominant "
                "short-distance formula with amplitude-dependent A_DeltaGamma "
                "time integration. Phase-3a RS semileptonic C10/C10' Wilsons "
                "enter additively; scalar matching remains deferred. "
                "The HARD ratio is the absolute NP branching-fraction shift "
                "over the YAML exp-vs-SM loose-edge budget."
            ),
            diagnostics=diagnostics,
        )
