"""B025 - charged-current LFU ratio ``R_D`` in ``B -> D tau nu``.

Physics
-------
The rigorous SM normalization is the HFLAV CKM-2025 SM reference in
``B025.yaml``.  This constraint does not hardcode the SM ratio; it loads the
YAML ``sm_reference_used_by_hflav.rd_value`` anchor and passes it into the
shared semileptonic-LFU core.  The v1 point-dependent prediction is

    R_D^proxy = R_D^SM |1 + C_tau^proxy|^2,

where ``C_tau^proxy`` is the documented charged-current proxy in
``quarkConstraints.semileptonic_lfu`` and scales like
``m_b m_tau / M_KK^2``.

Severity
--------
HARD.  ``R_D`` is an observed LFU ratio.  The HARD budget is built entirely
from ``B025.yaml`` as
``|R_D^exp - R_D^SM| + sqrt(sigma_exp^2 + sigma_SM^2)``.  The verdict compares
the total predicted ratio with the HFLAV average using that envelope, matching
the LFU-ratio treatment used by B018/B019.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B025.yaml`` is the source of truth for the
HFLAV average, SM reference, joint ``R_D``/``R_D*`` context, and provenance.

NEEDS-HUMAN-PHYSICS
-------------------
Full RS ``R_D`` matching requires the charged electroweak KK/W' tower,
charged-Higgs or leptoquark sector, lepton and neutrino couplings, a full
``b -> c tau nu`` WET basis, and ``B -> D`` form-factor integration.  Those
inputs are not available on ``ParameterPoint``; the implemented NP term is a
documented stress proxy, not a complete model prediction.
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
from flavor_catalog_constraints.physics_adapters.semileptonic_lfu import (
    SEMILEPTONIC_LFU_RS_MATCHING_ASSUMPTION_V1,
    rd_lfu_ratio_from_couplings,
    semileptonic_lfu_inputs_with_sm_ratio,
    semileptonic_lfu_sm_ratio,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPERIMENTAL_ANCHOR_CANDIDATES = ("canonical_average",)
_SM_ANCHOR_CANDIDATES = ("sm_reference_used_by_hflav",)
_JOINT_CONTEXT_CANDIDATES = ("joint_fit_context",)
_EXPECTED_UNITS = "dimensionless ratio"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B025.yaml canonical_average "
    "+ sm_reference_used_by_hflav"
)
_PARAMETRIZATION_CITATION = (
    "HFLAV CKM 2025 R(D) average and SM arithmetic-reference row from "
    "B025.yaml; charged-current scalar amplitude proxy in "
    "quarkConstraints.semileptonic_lfu"
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: full RS R_D matching needs charged EW KK/W', "
    "charged-Higgs or leptoquark sectors, lepton/neutrino couplings, the full "
    "b -> c tau nu WET operator basis, and B -> D form-factor integration; "
    "v1 uses the documented m_b m_tau/M_KK^2 charged-current proxy."
)


@dataclass(frozen=True)
class RDJointFitContext:
    """HFLAV joint ``R_D``/``R_D*`` context retained for diagnostics."""

    block_key: str
    rd_value: float | None
    rd_uncertainty: float | None
    rdstar_value: float
    rdstar_uncertainty: float
    correlation_rd_rdstar: float
    chi2: float | None
    dof: float | None
    confidence_level: float | None
    source: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class RDBudgetBand:
    """B025 measured-vs-SM envelope for the HARD verdict."""

    source: str
    central_residual: float
    experimental_sigma: float
    sm_theory_sigma: float
    combined_sigma: float
    hard_veto_budget: float
    lower_edge: float
    upper_edge: float
    construction: str


@dataclass(frozen=True)
class B025Anchor:
    """Typed B025 anchor: HFLAV average, SM reference, and budget."""

    experimental: Anchor
    standard_model: Anchor
    joint_fit_context: RDJointFitContext
    budget_band: RDBudgetBand

    @property
    def value(self) -> float:
        return self.experimental.value

    @property
    def uncertainty(self) -> float:
        if self.experimental.uncertainty is None:
            raise AnchorError("B025: experimental anchor lacks uncertainty")
        return float(self.experimental.uncertainty)

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
    return _required_float(value, process_id=process_id, field_name=field_name)


def _optional_positive_float(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> float | None:
    if value is None:
        return None
    return _positive_float(value, process_id=process_id, field_name=field_name)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: B025 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: B025 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(f"{process_id}: B025 anchor field {field_name!r} <= 0")
    return number


def _load_joint_context(process_id: str) -> RDJointFitContext:
    pdg_block = load_pdg_block(process_id, family=_FAMILY)
    sub = find_block(pdg_block, _JOINT_CONTEXT_CANDIDATES, process_id=process_id)
    block_key = next(key for key in _JOINT_CONTEXT_CANDIDATES if key in pdg_block)
    return RDJointFitContext(
        block_key=block_key,
        rd_value=_optional_positive_float(
            sub.get("rd_value"),
            process_id=process_id,
            field_name=f"{block_key}.rd_value",
        ),
        rd_uncertainty=_optional_positive_float(
            sub.get("rd_uncertainty"),
            process_id=process_id,
            field_name=f"{block_key}.rd_uncertainty",
        ),
        rdstar_value=_positive_float(
            sub.get("rdstar_value"),
            process_id=process_id,
            field_name=f"{block_key}.rdstar_value",
        ),
        rdstar_uncertainty=_positive_float(
            sub.get("rdstar_uncertainty"),
            process_id=process_id,
            field_name=f"{block_key}.rdstar_uncertainty",
        ),
        correlation_rd_rdstar=_required_float(
            sub.get("correlation_rd_rdstar"),
            process_id=process_id,
            field_name=f"{block_key}.correlation_rd_rdstar",
        ),
        chi2=_optional_float(
            sub.get("chi2"),
            process_id=process_id,
            field_name=f"{block_key}.chi2",
        ),
        dof=_optional_float(
            sub.get("dof"),
            process_id=process_id,
            field_name=f"{block_key}.dof",
        ),
        confidence_level=_optional_float(
            sub.get("confidence_level"),
            process_id=process_id,
            field_name=f"{block_key}.confidence_level",
        ),
        source=_optional_str(sub.get("source")),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _build_budget_band(*, experimental: Anchor, standard_model: Anchor) -> RDBudgetBand:
    if experimental.uncertainty is None or experimental.uncertainty <= 0.0:
        raise AnchorError("B025: experimental uncertainty is required")
    if standard_model.uncertainty is None or standard_model.uncertainty <= 0.0:
        raise AnchorError("B025: SM uncertainty is required")
    central = abs(experimental.value - standard_model.value)
    combined_sigma = math.sqrt(
        experimental.uncertainty**2 + standard_model.uncertainty**2
    )
    budget = central + combined_sigma
    if budget <= 0.0 or not math.isfinite(budget):
        raise AnchorError("B025: constructed R_D budget is invalid")
    return RDBudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(central),
        experimental_sigma=float(experimental.uncertainty),
        sm_theory_sigma=float(standard_model.uncertainty),
        combined_sigma=float(combined_sigma),
        hard_veto_budget=float(budget),
        lower_edge=float(experimental.value - budget),
        upper_edge=float(experimental.value + budget),
        construction=(
            "|R_D^exp - R_D^SM| + sqrt(sigma_exp^2 + sigma_SM^2) from "
            "B025.yaml HFLAV average and SM reference"
        ),
    )


def _load_b025_anchor(process_id: str) -> B025Anchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_EXPERIMENTAL_ANCHOR_CANDIDATES,
    )
    if experimental.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: expected experimental units {_EXPECTED_UNITS!r}, "
            f"got {experimental.units!r}"
        )
    standard_model = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_SM_ANCHOR_CANDIDATES,
        value_key="rd_value",
        uncertainty_key="rd_uncertainty",
    )
    if standard_model.value <= 0.0:
        raise AnchorError(f"{process_id}: SM R_D value must be positive")
    return B025Anchor(
        experimental=experimental,
        standard_model=standard_model,
        joint_fit_context=_load_joint_context(process_id),
        budget_band=_build_budget_band(
            experimental=experimental,
            standard_model=standard_model,
        ),
    )


def _budget_result(predicted: float, anchor: B025Anchor) -> tuple[float, float, bool]:
    budget = float(anchor.budget)
    pull = float(predicted - anchor.value)
    ratio = abs(pull) / budget if budget > 0.0 else float("inf")
    return float(budget), float(ratio), bool(ratio <= 1.0)


def _anchor_diagnostics(anchor: B025Anchor) -> dict[str, float | str | None]:
    return {
        "experimental_block": anchor.experimental.block_key,
        "sm_block": anchor.standard_model.block_key,
        "joint_context_block": anchor.joint_fit_context.block_key,
        "budget_source": anchor.budget_band.source,
        "budget_construction": anchor.budget_band.construction,
        "budget_central_residual": float(anchor.budget_band.central_residual),
        "budget_experimental_sigma": float(anchor.budget_band.experimental_sigma),
        "budget_sm_theory_sigma": float(anchor.budget_band.sm_theory_sigma),
        "budget_combined_sigma": float(anchor.budget_band.combined_sigma),
        "budget_lower_edge": float(anchor.budget_band.lower_edge),
        "budget_upper_edge": float(anchor.budget_band.upper_edge),
        "joint_fit_rdstar_value": float(anchor.joint_fit_context.rdstar_value),
        "joint_fit_rdstar_uncertainty": float(
            anchor.joint_fit_context.rdstar_uncertainty
        ),
        "joint_fit_correlation_rd_rdstar": float(
            anchor.joint_fit_context.correlation_rd_rdstar
        ),
        "joint_fit_chi2": anchor.joint_fit_context.chi2,
        "joint_fit_dof": anchor.joint_fit_context.dof,
        "joint_fit_confidence_level": anchor.joint_fit_context.confidence_level,
    }


@register
class Constraint:
    """Catalogued ``R_D`` charged-current LFU ratio constraint (B025)."""

    process_id = "B025"
    severity = Severity.HARD
    observable = "R_D = Gamma(B -> D tau nu) / Gamma(B -> D l nu)"

    def __init__(self) -> None:
        self.anchor = _load_b025_anchor(self.process_id)
        self.sm_inputs = semileptonic_lfu_inputs_with_sm_ratio(
            self.anchor.sm_value,
            mode="B->D",
        )
        self.sm_result = semileptonic_lfu_sm_ratio(self.sm_inputs)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                predicted=None,
                sm_prediction=float(self.sm_result.ratio),
                experimental=float(self.anchor.value),
                ratio=None,
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; R_D constraint was not "
                    "evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                    "rs_matching_assumption": SEMILEPTONIC_LFU_RS_MATCHING_ASSUMPTION_V1,
                    **_anchor_diagnostics(self.anchor),
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        result = rd_lfu_ratio_from_couplings(
            couplings,
            m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
            inputs=self.sm_inputs,
        )
        predicted = float(result.ratio)
        budget, ratio, passes = _budget_result(predicted, self.anchor)
        diagnostics: dict[str, Any] = dict(result.diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                "sm_anchor_rd": float(self.anchor.sm_value),
                "sm_formula_rd": float(result.sm_ratio),
                "sm_formula_minus_anchor": float(result.sm_ratio - self.anchor.sm_value),
                "np_shift_rd": float(result.np_shift_ratio),
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "rs_matching_assumption": SEMILEPTONIC_LFU_RS_MATCHING_ASSUMPTION_V1,
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                **_anchor_diagnostics(self.anchor),
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=passes,
            predicted=predicted,
            sm_prediction=float(result.sm_ratio),
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "R_D uses the YAML HFLAV SM reference and the semileptonic-LFU "
                "amplitude proxy R_D = R_D^SM |1 + C_tau^proxy|^2. The proxy "
                "scales as m_b m_tau/M_KK^2 and is marked NEEDS-HUMAN-PHYSICS. "
                "The HARD budget is the B025.yaml measured-vs-SM envelope."
            ),
            diagnostics=diagnostics,
        )
