"""B002 - CP phase in B_d -> J/psi K_S.

Physics
-------
``S_psiK_S`` is the mixing-induced CP asymmetry in
``B_d -> J/psi K_S``.  For this RS quark-flavor plugin the new physics is
restricted to the B_d Delta F = 2 mixing amplitude,

    phi_d^NP = arg(1 + M12^NP / M12^SM),
    S_psiK_S = sin(2 beta + phi_d^NP).

The complex ``M12^NP`` is evaluated through the Delta F = 2 adapter after
QCD-running the B_d Wilson coefficients to ``mu_had = 2 GeV``.  The SM
amplitude convention is ``M12^SM = Delta m_Bd^SM / 2`` from the same B_d core
inputs used by B001.  The ``2 beta`` phase is computed in core from the
repo-owned CKM target via the rephasing-invariant unitarity-triangle angle.

Severity
--------
HARD.  ``S_psiK_S`` is measured and the total predicted value must stay within
an uncertainty-aware dimensionless budget,

    sqrt(sigma_exp^2 + sigma_sin2beta^2 + sigma_penguin^2),

where the first term comes from the HFLAV average, the second uses the
documented HFLAV physical-beta uncertainty in ``B002.yaml`` around the in-core
CKM central value, and the penguin term converts the documented
``|Delta phi_d| <= 0.68 deg`` bound into an S shift using ``|cos(2 beta)|``.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B002.yaml`` is the source of truth for the
experimental HFLAV average, physical beta solution, and penguin phase bound.
"""

from __future__ import annotations

from dataclasses import dataclass
import cmath
import math
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.deltaf2 import (
    bd_mixing_core_inputs,
    bd_mixing_from_wilsons_with_running,
    bd_mixing_m12_np_from_wilsons_with_running,
    bd_mixing_wilsons_from_couplings,
)
from flavor_catalog_constraints.physics_adapters.ckm_extraction import (
    repo_default_ckm_phases,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_EXPECTED_OBSERVABLE = "S_psiK_S"
_MU_HAD_GEV = 2.0

_EXPERIMENTAL_ANCHOR_CANDIDATES = ("canonical_experimental_average",)
_MODE_SPECIFIC_ANCHOR_CANDIDATES = ("mode_specific_jpsi_ks",)
_BETA_SOLUTION_CANDIDATES = ("beta_physical_solution",)
_PENGUIN_PHASE_CANDIDATES = ("penguin_phase_bound",)
_BUDGET_DOC_CITATION = (
    "flavor_catalog/processes/beauty/B002.yaml:"
    "pdg_or_equivalent.canonical_experimental_average,"
    "beta_physical_solution,penguin_phase_bound; "
    "flavor_catalog/processes/beauty/B002.tex:"
    "Constraint validity and model dependence"
)
_SM_PHASE_SOURCE_POLICY = (
    "SM 2 beta computed in core from the repo-owned ModernDefaultCKMTarget; "
    "B002.yaml beta_physical_solution supplies only the uncertainty anchor."
)


@dataclass(frozen=True)
class BetaSolutionInputs:
    """Typed HFLAV physical beta solution from the B002 sidecar."""

    block_key: str
    source: str | None
    year: int | None
    beta_degrees: float
    upper_uncertainty_degrees: float
    lower_uncertainty_degrees: float
    sigma_beta_degrees: float
    two_beta_radians: float
    two_beta_degrees: float
    sigma_two_beta_radians: float
    alternate_solution_degrees: float | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class PenguinPhaseBound:
    """Typed long-distance penguin phase bound from the B002 sidecar."""

    block_key: str
    source: str | None
    year: int | None
    bound_abs_delta_phi_d_degrees: float
    bound_abs_delta_phi_d_radians: float
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class SpsiKsBudgetBand:
    """Uncertainty-aware S_psiK_S residual budget."""

    doc_citation: str
    central_residual: float
    hard_veto_budget: float
    experimental_sigma: float
    sm_sin2beta_sigma: float
    penguin_phase_sigma: float
    sm_sin2beta: float
    ckm_phase_source: str
    beta_radians: float
    beta_degrees: float
    two_beta_radians: float
    two_beta_degrees: float
    sigma_two_beta_radians: float
    penguin_bound_radians: float
    m12_sm_gev: float
    delta_m_bd_sm_gev: float
    core_delta_m_bd_exp_gev: float
    core_legacy_m12_budget: float


@dataclass(frozen=True)
class SpsiKsAnchor:
    """Typed B002 anchor: HFLAV average, beta solution, and budget."""

    experimental: Anchor
    mode_specific_jpsi_ks: Anchor
    beta_solution: BetaSolutionInputs
    penguin_phase_bound: PenguinPhaseBound
    budget_band: SpsiKsBudgetBand

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
        return self.budget_band.sm_sin2beta

    @property
    def budget(self) -> float:
        return self.budget_band.hard_veto_budget


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: S_psiK_S field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc


def _load_beta_solution(process_id: str) -> BetaSolutionInputs:
    beta_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_BETA_SOLUTION_CANDIDATES,
        value_key="value_degrees",
        uncertainty_key="upper_uncertainty_degrees",
    )
    block_key = beta_anchor.block_key
    beta_degrees = beta_anchor.value
    upper = abs(
        _required_float(
            beta_anchor.uncertainty,
            process_id=process_id,
            field_name=f"{block_key}.upper_uncertainty_degrees",
        )
    )
    lower_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_BETA_SOLUTION_CANDIDATES,
        value_key="lower_uncertainty_degrees",
    )
    alternate_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_BETA_SOLUTION_CANDIDATES,
        value_key="alternate_solution_degrees",
    )
    lower = abs(lower_anchor.value)
    sigma_beta_degrees = max(upper, lower)
    if beta_degrees <= 0.0 or sigma_beta_degrees <= 0.0:
        raise AnchorError(
            f"{process_id}: beta physical solution and uncertainty must be positive"
        )

    two_beta_degrees = 2.0 * beta_degrees
    return BetaSolutionInputs(
        block_key=block_key,
        source=beta_anchor.source,
        year=beta_anchor.year,
        beta_degrees=beta_degrees,
        upper_uncertainty_degrees=upper,
        lower_uncertainty_degrees=lower,
        sigma_beta_degrees=sigma_beta_degrees,
        two_beta_radians=math.radians(two_beta_degrees),
        two_beta_degrees=two_beta_degrees,
        sigma_two_beta_radians=math.radians(2.0 * sigma_beta_degrees),
        alternate_solution_degrees=alternate_anchor.value,
        source_url=beta_anchor.source_url,
        snapshot_path=beta_anchor.snapshot_path,
    )


def _load_penguin_phase_bound(process_id: str) -> PenguinPhaseBound:
    bound_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_PENGUIN_PHASE_CANDIDATES,
        value_key="bound_abs_delta_phi_d_degrees",
    )
    block_key = bound_anchor.block_key
    bound_degrees = abs(bound_anchor.value)
    if bound_degrees <= 0.0:
        raise AnchorError(
            f"{process_id}: penguin phase bound must be positive, got {bound_degrees}"
        )
    return PenguinPhaseBound(
        block_key=block_key,
        source=bound_anchor.source,
        year=bound_anchor.year,
        bound_abs_delta_phi_d_degrees=bound_degrees,
        bound_abs_delta_phi_d_radians=math.radians(bound_degrees),
        source_url=bound_anchor.source_url,
        snapshot_path=bound_anchor.snapshot_path,
    )


def _build_budget_band(
    *,
    process_id: str,
    experimental: Anchor,
    beta_solution: BetaSolutionInputs,
    penguin_phase_bound: PenguinPhaseBound,
) -> SpsiKsBudgetBand:
    if experimental.uncertainty is None or experimental.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: S_psiK_S experimental uncertainty is required "
            "for the uncertainty-aware budget band"
        )
    if experimental.units != "dimensionless":
        raise AnchorError(
            f"{process_id}: S_psiK_S experimental anchor must be dimensionless, "
            f"got units={experimental.units!r}"
        )

    core_inputs = bd_mixing_core_inputs()
    delta_m_bd_sm = core_inputs["delta_m_bd_sm_gev"]
    m12_sm = delta_m_bd_sm / 2.0
    if m12_sm <= 0.0:
        raise AnchorError(
            f"{process_id}: core B_d M12^SM must be positive, got {m12_sm}"
        )

    ckm_phase = repo_default_ckm_phases()
    sm_sin2beta = ckm_phase.sin_2beta
    cos_two_beta_abs = abs(math.cos(ckm_phase.two_beta))
    sm_sigma = cos_two_beta_abs * beta_solution.sigma_two_beta_radians
    penguin_sigma = (
        cos_two_beta_abs * penguin_phase_bound.bound_abs_delta_phi_d_radians
    )
    budget = math.sqrt(
        float(experimental.uncertainty) * float(experimental.uncertainty)
        + sm_sigma * sm_sigma
        + penguin_sigma * penguin_sigma
    )
    central_residual = abs(float(experimental.value) - sm_sin2beta)
    if budget <= 0.0 or not math.isfinite(budget):
        raise AnchorError(
            f"{process_id}: S_psiK_S budget must be positive and finite, got {budget}"
        )

    return SpsiKsBudgetBand(
        doc_citation=_BUDGET_DOC_CITATION,
        central_residual=central_residual,
        hard_veto_budget=budget,
        experimental_sigma=float(experimental.uncertainty),
        sm_sin2beta_sigma=sm_sigma,
        penguin_phase_sigma=penguin_sigma,
        sm_sin2beta=sm_sin2beta,
        ckm_phase_source=ckm_phase.source,
        beta_radians=ckm_phase.beta,
        beta_degrees=ckm_phase.beta_degrees,
        two_beta_radians=ckm_phase.two_beta,
        two_beta_degrees=ckm_phase.two_beta_degrees,
        sigma_two_beta_radians=beta_solution.sigma_two_beta_radians,
        penguin_bound_radians=penguin_phase_bound.bound_abs_delta_phi_d_radians,
        m12_sm_gev=m12_sm,
        delta_m_bd_sm_gev=delta_m_bd_sm,
        core_delta_m_bd_exp_gev=core_inputs["delta_m_bd_exp_gev"],
        core_legacy_m12_budget=core_inputs["core_m12_budget_gev"],
    )


def _load_spsi_ks_anchor(process_id: str) -> SpsiKsAnchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_EXPERIMENTAL_ANCHOR_CANDIDATES,
    )
    mode_specific = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_MODE_SPECIFIC_ANCHOR_CANDIDATES,
    )
    beta_solution = _load_beta_solution(process_id)
    penguin_phase_bound = _load_penguin_phase_bound(process_id)
    anchor = SpsiKsAnchor(
        experimental=experimental,
        mode_specific_jpsi_ks=mode_specific,
        beta_solution=beta_solution,
        penguin_phase_bound=penguin_phase_bound,
        budget_band=_build_budget_band(
            process_id=process_id,
            experimental=experimental,
            beta_solution=beta_solution,
            penguin_phase_bound=penguin_phase_bound,
        ),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(
            f"{process_id}: S_psiK_S NP budget must be positive, got {anchor.budget}"
        )
    return anchor


def _complex_mapping(values: Mapping[str, complex]) -> dict[str, complex]:
    return {str(key): complex(value) for key, value in values.items()}


@register
class Constraint:
    """Catalogued S_psiK_S / sin(2 beta) constraint (process_id B002)."""

    process_id = "B002"
    severity = Severity.HARD
    observable = _EXPECTED_OBSERVABLE

    def __init__(self) -> None:
        self.anchor = _load_spsi_ks_anchor(self.process_id)
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
                    f"extra {_REQUIRED_EXTRA!r} absent; S_psiK_S constraint "
                    "was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "sm_phase_source_policy": _SM_PHASE_SOURCE_POLICY,
                    "ckm_phase_source": self.anchor.budget_band.ckm_phase_source,
                },
            )

        wilsons = bd_mixing_wilsons_from_couplings(couplings)
        m12_np = bd_mixing_m12_np_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
        )
        core_magnitude = bd_mixing_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
        )

        m12_sm = self.anchor.budget_band.m12_sm_gev
        m12_ratio = complex(m12_np / m12_sm)
        total_mixing_ratio = 1.0 + m12_ratio
        phi_d_np = float(cmath.phase(total_mixing_ratio))
        predicted = float(
            math.sin(self.anchor.budget_band.two_beta_radians + phi_d_np)
        )
        residual = abs(predicted - float(self.anchor.value))
        budget = float(self.anchor.budget)
        ratio = residual / budget

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(ratio <= 1.0),
            predicted=predicted,
            sm_prediction=float(self.anchor.sm_value),
            experimental=float(self.anchor.value),
            ratio=float(ratio),
            budget=budget,
            notes=(
                "S_psiK_S = sin(2 beta + arg(1 + M12_Bd^NP/M12_Bd^SM)); "
                "M12_Bd^NP is the complex QCD-running Delta F=2 amplitude at "
                "2 GeV. 2 beta is computed in core from the repo CKM target."
            ),
            diagnostics={
                "m12_np_gev": complex(m12_np),
                "abs_m12_np_gev": float(abs(m12_np)),
                "core_abs_m12_np_gev": float(core_magnitude.abs_m12_np),
                "m12_sm_gev": float(m12_sm),
                "m12_np_over_m12_sm": m12_ratio,
                "re_m12_np_over_m12_sm": float(m12_ratio.real),
                "im_m12_np_over_m12_sm": float(m12_ratio.imag),
                "total_mixing_ratio": total_mixing_ratio,
                "phi_d_np_rad": phi_d_np,
                "phi_d_np_deg": math.degrees(phi_d_np),
                "s_psi_ks_total": predicted,
                "s_psi_ks_residual": residual,
                "qcd_running_applied": True,
                "hadronic_scale_gev": _MU_HAD_GEV,
                "matching_scale_gev": float(wilsons.matching_scale),
                "m_kk_gev": float(wilsons.M_KK),
                "core_input_key": wilsons.input.key,
                "down_sector_indices": tuple(wilsons.input.generations),
                "left_db_coupling": complex(wilsons.left_coupling),
                "right_db_coupling": complex(wilsons.right_coupling),
                "wilson_coefficients": _complex_mapping(wilsons.wilsons),
                "phase_formula": "phi_d_np = arg(1 + M12_NP / M12_SM)",
                "phase_uses_complex_m12_not_abs": True,
                "budget_doc_citation": self.anchor.budget_band.doc_citation,
                "budget_construction": (
                    "sqrt(sigma_exp^2 + sigma_sin2beta_SM^2 + "
                    "sigma_penguin^2)"
                ),
                "central_sm_exp_residual": (
                    self.anchor.budget_band.central_residual
                ),
                "hard_veto_budget": self.anchor.budget_band.hard_veto_budget,
                "budget_experimental_sigma": (
                    self.anchor.budget_band.experimental_sigma
                ),
                "budget_sm_sin2beta_sigma": (
                    self.anchor.budget_band.sm_sin2beta_sigma
                ),
                "budget_penguin_phase_sigma": (
                    self.anchor.budget_band.penguin_phase_sigma
                ),
                "two_beta_rad": self.anchor.budget_band.two_beta_radians,
                "two_beta_deg": self.anchor.budget_band.two_beta_degrees,
                "beta_rad": self.anchor.budget_band.beta_radians,
                "beta_deg": self.anchor.budget_band.beta_degrees,
                "ckm_phase_source": self.anchor.budget_band.ckm_phase_source,
                "sigma_two_beta_rad": (
                    self.anchor.budget_band.sigma_two_beta_radians
                ),
                "penguin_bound_rad": self.anchor.budget_band.penguin_bound_radians,
                "delta_m_bd_sm_gev": self.anchor.budget_band.delta_m_bd_sm_gev,
                "core_delta_m_bd_exp_gev": (
                    self.anchor.budget_band.core_delta_m_bd_exp_gev
                ),
                "core_legacy_m12_budget": (
                    self.anchor.budget_band.core_legacy_m12_budget
                ),
                "core_legacy_ratio_to_budget": float(
                    core_magnitude.ratio_to_budget
                ),
                "sm_phase_source_policy": _SM_PHASE_SOURCE_POLICY,
                "experimental_block": self.anchor.experimental.block_key,
                "mode_specific_block": (
                    self.anchor.mode_specific_jpsi_ks.block_key
                ),
                "mode_specific_jpsi_ks_value": (
                    self.anchor.mode_specific_jpsi_ks.value
                ),
                "mode_specific_jpsi_ks_uncertainty": (
                    self.anchor.mode_specific_jpsi_ks.uncertainty
                ),
                "beta_solution_block": self.anchor.beta_solution.block_key,
                "beta_solution_degrees": self.anchor.beta_solution.beta_degrees,
                "beta_solution_sigma_degrees": (
                    self.anchor.beta_solution.sigma_beta_degrees
                ),
                "beta_solution_source_url": self.anchor.beta_solution.source_url,
                "penguin_phase_block": (
                    self.anchor.penguin_phase_bound.block_key
                ),
                "penguin_phase_bound_degrees": (
                    self.anchor.penguin_phase_bound.bound_abs_delta_phi_d_degrees
                ),
                "penguin_phase_source_url": (
                    self.anchor.penguin_phase_bound.source_url
                ),
            },
        )
