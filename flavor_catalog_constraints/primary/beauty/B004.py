"""B004 - CP phase in B_s -> J/psi phi.

Physics
-------
``phi_s`` is the mixing-induced CP phase in tree-dominated
``B_s -> J/psi phi`` / ``B_s -> J/psi K+K-`` analyses.  For this RS
quark-flavor plugin the new physics is restricted to the complex
``B_s`` Delta F = 2 mixing amplitude,

    phi_s^NP = arg(1 + M12^NP / M12^SM),
    phi_s = phi_s^SM + phi_s^NP.

The complex ``M12^NP`` is evaluated through the Delta F = 2 adapter after
QCD-running the B_s Wilson coefficients to ``mu_had = m_b``.  The SM
amplitude has magnitude ``|M12^SM| = Delta m_Bs^SM / 2`` from the same B_s
core inputs used by B003 and phase ``arg((V_ts* V_tb)^2)`` from the fitted
point CKM matrix.  The SM ``phi_s`` reference phase is computed in core from
the repo-owned CKM target as ``phi_s = -2 beta_s``.

Severity
--------
HARD.  ``phi_s`` is measured and the total predicted phase must stay within
an uncertainty-aware radian budget that combines the HFLAV experimental
uncertainty and the asymmetric SM ``-2 beta_s`` uncertainty from ``B004.yaml``
in quadrature.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B004.yaml`` is the source of truth for the
HFLAV phase average, the mode-specific ``J/psi K+K-`` average, the latest LHCb
input, and the SM ``-2 beta_s`` uncertainty.  The central SM phase is computed
from the in-core CKM matrix target.
"""

from __future__ import annotations

from dataclasses import dataclass
import cmath
import math
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.ckm_extraction import (
    neutral_b_mixing_sm_amplitude,
    repo_default_ckm_phases,
)
from flavor_catalog_constraints.physics_adapters.deltaf2 import (
    bs_mixing_core_inputs,
    bs_mixing_from_wilsons_with_running,
    bs_mixing_m12_np_from_wilsons_with_running,
    bs_mixing_wilsons_from_couplings,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_EXPECTED_UNITS = "rad"
_EXPECTED_OBSERVABLE = "phi_s"
_PHASE_OBSERVABLE_TOKENS = ("phi_s", "beta_s")
_MU_HAD_GEV = 4.18

_EXPERIMENTAL_ANCHOR_CANDIDATES = ("canonical_hflav_average",)
_MODE_SPECIFIC_ANCHOR_CANDIDATES = ("mode_specific_jpsi_kk_average",)
_LATEST_LHCB_ANCHOR_CANDIDATES = ("latest_lhcb_jpsi_kk_input",)
_SM_REFERENCE_CANDIDATES = ("standard_model_reference",)
_BUDGET_DOC_CITATION = (
    "flavor_catalog/processes/beauty/B004.yaml:"
    "pdg_or_equivalent.canonical_hflav_average,standard_model_reference; "
    "flavor_catalog/processes/beauty/B004.tex:"
    "Constraint validity / model dependence"
)
_SM_PHASE_SOURCE_POLICY = (
    "SM phi_s = -2 beta_s computed in core from the repo-owned "
    "ModernDefaultCKMTarget; B004.yaml standard_model_reference supplies only "
    "the uncertainty anchor."
)
_SM_BOX_CKM_SOURCE_POLICY = (
    "M12_SM uses the fitted QuarkMassBasisCouplings.ckm_matrix when present; "
    "legacy/manual coupling objects fall back to the repo-owned "
    "ModernDefaultCKMTarget. The phase convention is (conj(V_ts) V_tb)^2."
)


@dataclass(frozen=True)
class PhiSStandardModelReference:
    """Typed SM ``-2 beta_s`` phase reference and uncertainty provenance."""

    block_key: str
    source: str | None
    year: int | None
    value: float
    yaml_value: float
    upper_uncertainty: float
    lower_uncertainty: float
    sigma: float
    ckm_phase_source: str
    beta_s: float
    beta_s_degrees: float
    units: str | None
    observable: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class PhiSBudgetBand:
    """Uncertainty-aware B004 phase budget in radians."""

    doc_citation: str
    central_residual: float
    hard_veto_budget: float
    experimental_sigma: float
    sm_phase_sigma_upper: float
    sm_phase_sigma_lower: float
    combined_sigma_upper: float
    combined_sigma_lower: float
    sm_phi_s: float
    m12_sm_gev: float
    delta_m_bs_sm_gev: float
    core_delta_m_bs_exp_gev: float
    core_legacy_m12_budget: float


@dataclass(frozen=True)
class PhiSAnchor:
    """Typed B004 anchor: HFLAV experiment, SM phase, and budget."""

    experimental: Anchor
    mode_specific_jpsi_kk: Anchor
    latest_lhcb_jpsi_kk: Anchor
    standard_model: PhiSStandardModelReference
    budget_band: PhiSBudgetBand

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


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: phi_s field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: phi_s field {field_name!r}={value!r} is not finite"
        )
    return number


def _wrap_phase_radians(value: float) -> float:
    """Wrap a phase to the principal interval ``[-pi, pi]``."""
    return float(math.atan2(math.sin(value), math.cos(value)))


def _validate_phase_anchor(anchor: Anchor, *, process_id: str, name: str) -> None:
    if anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: {name} anchor {anchor.block_key!r} must use units "
            f"{_EXPECTED_UNITS!r}, got {anchor.units!r}"
        )
    if anchor.observable is None or not any(
        token in anchor.observable for token in _PHASE_OBSERVABLE_TOKENS
    ):
        raise AnchorError(
            f"{process_id}: {name} anchor {anchor.block_key!r} must describe "
            f"a phi_s/beta_s phase, got observable={anchor.observable!r}"
        )
    _required_float(anchor.value, process_id=process_id, field_name=f"{name}.value")


def _load_sm_reference(process_id: str) -> PhiSStandardModelReference:
    sm_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_SM_REFERENCE_CANDIDATES,
        uncertainty_key="upper_uncertainty",
    )
    _validate_phase_anchor(sm_anchor, process_id=process_id, name="standard_model")
    if sm_anchor.uncertainty is None:
        raise AnchorError(
            f"{process_id}: standard_model_reference.upper_uncertainty is required"
        )
    lower_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_SM_REFERENCE_CANDIDATES,
        value_key="lower_uncertainty",
    )
    upper = abs(
        _required_float(
            sm_anchor.uncertainty,
            process_id=process_id,
            field_name="standard_model_reference.upper_uncertainty",
        )
    )
    lower = abs(
        _required_float(
            lower_anchor.value,
            process_id=process_id,
            field_name="standard_model_reference.lower_uncertainty",
        )
    )
    if upper <= 0.0 or lower <= 0.0:
        raise AnchorError(
            f"{process_id}: SM phi_s asymmetric uncertainties must be positive"
        )
    ckm_phase = repo_default_ckm_phases()
    return PhiSStandardModelReference(
        block_key=sm_anchor.block_key,
        source=sm_anchor.source,
        year=sm_anchor.year,
        value=ckm_phase.phi_s,
        yaml_value=sm_anchor.value,
        upper_uncertainty=upper,
        lower_uncertainty=lower,
        sigma=max(upper, lower),
        ckm_phase_source=ckm_phase.source,
        beta_s=ckm_phase.beta_s,
        beta_s_degrees=ckm_phase.beta_s_degrees,
        units=sm_anchor.units,
        observable=sm_anchor.observable,
        source_url=sm_anchor.source_url,
        snapshot_path=sm_anchor.snapshot_path,
    )


def _build_budget_band(
    *,
    process_id: str,
    experimental: Anchor,
    standard_model: PhiSStandardModelReference,
) -> PhiSBudgetBand:
    if experimental.uncertainty is None or experimental.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: phi_s experimental uncertainty is required for "
            "the uncertainty-aware budget band"
        )

    core_inputs = bs_mixing_core_inputs()
    delta_m_bs_sm = core_inputs["delta_m_bs_sm_gev"]
    m12_sm = delta_m_bs_sm / 2.0
    if m12_sm <= 0.0 or not math.isfinite(m12_sm):
        raise AnchorError(
            f"{process_id}: core B_s M12^SM must be positive and finite, "
            f"got {m12_sm}"
        )

    exp_sigma = float(experimental.uncertainty)
    combined_upper = math.sqrt(
        exp_sigma * exp_sigma
        + standard_model.upper_uncertainty * standard_model.upper_uncertainty
    )
    combined_lower = math.sqrt(
        exp_sigma * exp_sigma
        + standard_model.lower_uncertainty * standard_model.lower_uncertainty
    )
    central_residual = abs(
        _wrap_phase_radians(standard_model.value - experimental.value)
    )
    hard_budget = max(combined_upper, combined_lower)
    if hard_budget <= 0.0 or not math.isfinite(hard_budget):
        raise AnchorError(
            f"{process_id}: phi_s phase budget must be positive and finite, "
            f"got {hard_budget}"
        )

    return PhiSBudgetBand(
        doc_citation=_BUDGET_DOC_CITATION,
        central_residual=central_residual,
        hard_veto_budget=hard_budget,
        experimental_sigma=exp_sigma,
        sm_phase_sigma_upper=standard_model.upper_uncertainty,
        sm_phase_sigma_lower=standard_model.lower_uncertainty,
        combined_sigma_upper=combined_upper,
        combined_sigma_lower=combined_lower,
        sm_phi_s=standard_model.value,
        m12_sm_gev=m12_sm,
        delta_m_bs_sm_gev=delta_m_bs_sm,
        core_delta_m_bs_exp_gev=core_inputs["delta_m_bs_exp_gev"],
        core_legacy_m12_budget=core_inputs["core_m12_budget_gev"],
    )


def _load_phi_s_anchor(process_id: str) -> PhiSAnchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_EXPERIMENTAL_ANCHOR_CANDIDATES,
    )
    _validate_phase_anchor(experimental, process_id=process_id, name="experimental")
    mode_specific = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_MODE_SPECIFIC_ANCHOR_CANDIDATES,
    )
    _validate_phase_anchor(
        mode_specific,
        process_id=process_id,
        name="mode_specific_jpsi_kk",
    )
    latest_lhcb = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_LATEST_LHCB_ANCHOR_CANDIDATES,
    )
    _validate_phase_anchor(latest_lhcb, process_id=process_id, name="latest_lhcb")
    standard_model = _load_sm_reference(process_id)
    anchor = PhiSAnchor(
        experimental=experimental,
        mode_specific_jpsi_kk=mode_specific,
        latest_lhcb_jpsi_kk=latest_lhcb,
        standard_model=standard_model,
        budget_band=_build_budget_band(
            process_id=process_id,
            experimental=experimental,
            standard_model=standard_model,
        ),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(
            f"{process_id}: phi_s budget must be positive, got {anchor.budget}"
        )
    return anchor


def _complex_mapping(values: Mapping[str, complex]) -> dict[str, complex]:
    return {str(key): complex(value) for key, value in values.items()}


def _selected_budget(
    predicted: float,
    anchor: PhiSAnchor,
) -> tuple[float, float, float, bool]:
    residual_signed = _wrap_phase_radians(predicted - anchor.value)
    budget = (
        anchor.budget_band.combined_sigma_upper
        if residual_signed >= 0.0
        else anchor.budget_band.combined_sigma_lower
    )
    residual = abs(residual_signed)
    ratio = residual / budget if budget > 0.0 else float("inf")
    return float(residual), float(budget), float(ratio), bool(ratio <= 1.0)


@register
class Constraint:
    """Catalogued phi_s / B_s Delta F=2 phase constraint (process_id B004)."""

    process_id = "B004"
    severity = Severity.HARD
    observable = "phi_s(B_s -> J/psi phi)"

    def __init__(self) -> None:
        self.anchor = _load_phi_s_anchor(self.process_id)
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
                    f"extra {_REQUIRED_EXTRA!r} absent; phi_s constraint was "
                    "not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "sm_phase_source_policy": _SM_PHASE_SOURCE_POLICY,
                    "ckm_phase_source": self.anchor.standard_model.ckm_phase_source,
                },
            )

        wilsons = bs_mixing_wilsons_from_couplings(couplings)
        m12_np = bs_mixing_m12_np_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
        )
        core_magnitude = bs_mixing_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
        )

        sm_box = neutral_b_mixing_sm_amplitude(
            delta_m_sm_gev=self.anchor.budget_band.delta_m_bs_sm_gev,
            light_down_index=1,
            ckm=getattr(couplings, "ckm_matrix", None),
            ckm_source=getattr(couplings, "ckm_source", None),
        )
        m12_sm = sm_box.amplitude_gev
        m12_ratio = complex(m12_np / m12_sm)
        total_mixing_ratio = 1.0 + m12_ratio
        phi_s_np = float(cmath.phase(total_mixing_ratio))
        predicted = _wrap_phase_radians(self.anchor.sm_value + phi_s_np)
        residual, budget, ratio, passes = _selected_budget(predicted, self.anchor)

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=predicted,
            sm_prediction=float(self.anchor.sm_value),
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "phi_s = phi_s^SM + arg(1 + M12_Bs^NP/M12_Bs^SM); "
                "M12_Bs^NP is the complex QCD-running Delta F=2 amplitude at "
                "m_b=4.18 GeV. phi_s^SM is computed in core from the repo CKM target."
            ),
            diagnostics={
                "m12_np_gev": complex(m12_np),
                "abs_m12_np_gev": float(abs(m12_np)),
                "im_m12_np_gev": float(m12_np.imag),
                "core_abs_m12_np_gev": float(core_magnitude.abs_m12_np),
                "m12_sm_gev": float(sm_box.magnitude_gev),
                "m12_sm_box_gev": complex(m12_sm),
                "m12_sm_box_phase_rad": float(cmath.phase(m12_sm)),
                "m12_sm_box_ckm_factor": complex(sm_box.ckm_factor),
                "m12_sm_box_ckm_source": sm_box.ckm_source,
                "m12_np_over_m12_sm": m12_ratio,
                "re_m12_np_over_m12_sm": float(m12_ratio.real),
                "im_m12_np_over_m12_sm": float(m12_ratio.imag),
                "total_mixing_ratio": total_mixing_ratio,
                "phi_s_np_rad": phi_s_np,
                "phi_s_np_deg": math.degrees(phi_s_np),
                "phi_s_total_rad": predicted,
                "phi_s_residual_rad": residual,
                "qcd_running_applied": True,
                "hadronic_scale_gev": _MU_HAD_GEV,
                "matching_scale_gev": float(wilsons.matching_scale),
                "m_kk_gev": float(wilsons.M_KK),
                "core_input_key": wilsons.input.key,
                "down_sector_indices": tuple(wilsons.input.generations),
                "left_sb_coupling": complex(wilsons.left_coupling),
                "right_sb_coupling": complex(wilsons.right_coupling),
                "wilson_coefficients": _complex_mapping(wilsons.wilsons),
                "phase_formula": "phi_s_np = arg(1 + M12_NP / M12_SM)",
                "phase_uses_complex_m12_not_abs": True,
                "m12_sm_convention": (
                    "|M12_SM| = Delta m_Bs^SM / 2; "
                    "phase = arg((conj(V_ts) V_tb)^2)"
                ),
                "budget_doc_citation": self.anchor.budget_band.doc_citation,
                "budget_construction": (
                    "sqrt(sigma_exp^2 + sigma_phi_s_SM^2), direction-aware "
                    "for asymmetric SM uncertainty"
                ),
                "central_sm_exp_residual": (
                    self.anchor.budget_band.central_residual
                ),
                "hard_veto_budget": self.anchor.budget_band.hard_veto_budget,
                "budget_combined_sigma_upper": (
                    self.anchor.budget_band.combined_sigma_upper
                ),
                "budget_combined_sigma_lower": (
                    self.anchor.budget_band.combined_sigma_lower
                ),
                "budget_experimental_sigma": (
                    self.anchor.budget_band.experimental_sigma
                ),
                "budget_sm_phase_sigma_upper": (
                    self.anchor.budget_band.sm_phase_sigma_upper
                ),
                "budget_sm_phase_sigma_lower": (
                    self.anchor.budget_band.sm_phase_sigma_lower
                ),
                "delta_m_bs_sm_gev": self.anchor.budget_band.delta_m_bs_sm_gev,
                "core_delta_m_bs_exp_gev": (
                    self.anchor.budget_band.core_delta_m_bs_exp_gev
                ),
                "core_legacy_m12_budget": (
                    self.anchor.budget_band.core_legacy_m12_budget
                ),
                "core_legacy_ratio_to_budget": float(
                    core_magnitude.ratio_to_budget
                ),
                "sm_phase_source_policy": _SM_PHASE_SOURCE_POLICY,
                "sm_box_ckm_source_policy": _SM_BOX_CKM_SOURCE_POLICY,
                "ckm_phase_source": self.anchor.standard_model.ckm_phase_source,
                "beta_s_rad": self.anchor.standard_model.beta_s,
                "beta_s_deg": self.anchor.standard_model.beta_s_degrees,
                "experimental_block": self.anchor.experimental.block_key,
                "mode_specific_block": (
                    self.anchor.mode_specific_jpsi_kk.block_key
                ),
                "mode_specific_jpsi_kk_value": (
                    self.anchor.mode_specific_jpsi_kk.value
                ),
                "mode_specific_jpsi_kk_uncertainty": (
                    self.anchor.mode_specific_jpsi_kk.uncertainty
                ),
                "latest_lhcb_block": self.anchor.latest_lhcb_jpsi_kk.block_key,
                "latest_lhcb_jpsi_kk_value": (
                    self.anchor.latest_lhcb_jpsi_kk.value
                ),
                "sm_reference_block": self.anchor.standard_model.block_key,
                "sm_reference_yaml_value": self.anchor.standard_model.yaml_value,
                "sm_reference_source_url": self.anchor.standard_model.source_url,
            },
        )
