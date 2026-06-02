"""B014 - exclusive radiative ``b -> d gamma`` decays.

Physics
-------
The exclusive ``B+ -> rho+ gamma``, ``B0 -> rho0 gamma``, and
``B0 -> omega gamma`` branching fractions are evaluated with the shared C7/C7'
dipole response used by B011/B012/B013,

    BR(B -> rho/omega gamma) = BR_norm^(b->d) *
        (|C7_SM + C7_NP^(d)|^2 + |C7p_SM + C7p_NP^(d)|^2)
        / (|C7_SM|^2 + |C7p_SM|^2).

The C7/C8 matching proxy and leading-log QCD running to ``mu_b`` live behind
``flavor_catalog_constraints.physics_adapters.bsgamma``.  B014 uses the
append-only b->d adapter wrapper there: the down-sector d-b entry ``(0, 2)``
is selected and passed through the unchanged bsgamma C7/C8 machinery.  The
repo-default CKM target is recorded in diagnostics, including the expected
``|V_td/V_ts|^2`` suppression relative to b->s gamma.

Severity
--------
HARD.  ``B014.yaml`` has no theory-only exclusive SM prediction/form-factor
block.  The PDG rows are therefore used only as measurement-normalization
proxies, with per-mode measurement-consistency bands,

    |BR_total(mode) - BR_ref_measurement(mode)| <= sigma_measurement(mode).

The HARD decision uses the largest per-mode saturation.  This is deliberately
not presented as SM-theory-vs-data room.  NEEDS-HUMAN-PHYSICS: a rigorous RS
exclusive constraint needs full b->d dipole matching plus B->rho/omega tensor
form factors, weak-annihilation/spectator amplitudes, and a likelihood for the
branching fractions, isospin, and CP/helicity observables.

Catalog sidecar
---------------
``flavor_catalog/processes/secondary/beauty/B014.yaml`` is the source of truth
for the PDG/HFLAV branching-fraction rows and supporting b->d gamma
observables.  Numeric branching-fraction values below are loaded from that
sidecar, not hardcoded in this constraint.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_pdg_block,
)
from flavor_catalog_constraints.base import (
    ConstraintLevel,
    ConstraintResult,
    ParameterPoint,
    Severity,
)
from flavor_catalog_constraints.physics_adapters.bsgamma import (
    BSGAMMA_RS_MATCHING_ASSUMPTION_V1,
    bsgamma_default_sm_inputs,
    bdgamma_ckm_factors,
    exclusive_bdrhogamma_from_couplings,
    exclusive_bdrhogamma_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_TIER = ConstraintLevel.SECONDARY
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_DIPOLE_MASS_EXTRA = "kk_ew_mass_gev"
_CHARGED_RHO_CANDIDATES = ("pdg2025_bp_to_rhop_gamma",)
_NEUTRAL_RHO_CANDIDATES = ("pdg2025_b0_to_rho0_gamma",)
_OMEGA_CANDIDATES = ("pdg2025_b0_to_omega_gamma",)
_HFLAV_RHO_CANDIDATES = ("hflav2024_b_to_rho_gamma",)
_HFLAV_RHO_OMEGA_CANDIDATES = ("hflav2024_b_to_rho_omega_gamma",)
_ISOSPIN_CANDIDATES = ("hflav2024_isospin_rho_gamma",)
_XD_GAMMA_CANDIDATES = ("hflav2024_b_to_xd_gamma",)
_LHCB_RATIO_CANDIDATES = ("lhcb2025_ratio_b0rho0_to_b0kstar0_gamma",)
_EXPECTED_BRANCHING_UNITS = "branching fraction"
_EXPECTED_DIMENSIONLESS_UNITS = "dimensionless"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/secondary/beauty/B014.yaml "
    "measurement-consistency bands from PDG 2025 B+ -> rho+ gamma, "
    "B0 -> rho0 gamma, and B0 -> omega gamma rows"
)
_NORMALIZATION_LABEL = "B014 PDG exclusive b -> d gamma measurement proxy"
_PARAMETRIZATION_CITATION = (
    "Exclusive C7-normalized B -> rho/omega gamma proxy: B014.yaml supplies "
    "PDG measurement normalization proxies, while C7/C8 matching and LL QCD "
    "running are reused from quarkConstraints.bsgamma through the bsgamma "
    "adapter's b->d wrapper."
)
_BUDGET_POLICY = (
    "measurement-consistency band: each exclusive total BR is compared with "
    "its own PDG measurement proxy and measurement uncertainty; this is not a "
    "SM-theory-vs-data budget because B014.yaml has no theory-only exclusive "
    "SM prediction."
)
_SM_THEORY_GAP_NOTE = (
    "NEEDS-HUMAN-PHYSICS: B014.yaml contains no theory-only QCDF/SCET/LCSR "
    "exclusive SM prediction for BR(B+ -> rho+ gamma), BR(B0 -> rho0 gamma), "
    "or BR(B0 -> omega gamma). A production constraint needs B->rho/omega "
    "form factors, weak-annihilation/spectator treatment, and correlated "
    "exclusive theory uncertainties."
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: full RS exclusive b -> d gamma matching requires "
    "KK fermion, Higgs/Goldstone, charged-current, chromomagnetic, B->rho/"
    "omega form-factor, weak-annihilation, spectator-amplitude, isospin, and "
    "right-handed-photon inputs not available on ParameterPoint; v1 reuses "
    "the documented C7/C8 dipole proxy with the d-b down-sector coupling and "
    "LL QCD running to mu_b. "
    + _SM_THEORY_GAP_NOTE
)


@dataclass(frozen=True)
class B014ModeBudgetBand:
    """Per-mode B014 measurement-consistency NP-shift budget."""

    source: str
    mode_key: str
    experimental_sigma: float
    hard_veto_budget: float
    lower_total_edge: float
    upper_total_edge: float


@dataclass(frozen=True)
class B014ModeAnchor:
    """One exclusive B014 branching-fraction normalization row."""

    mode_key: str
    mode_label: str
    anchor: Anchor
    raw: Mapping[str, Any]
    budget_band: B014ModeBudgetBand

    @property
    def value(self) -> float:
        return self.anchor.value

    @property
    def uncertainty(self) -> float:
        return self.budget_band.experimental_sigma

    @property
    def reference_value(self) -> float:
        return self.anchor.value

    @property
    def budget(self) -> float:
        return self.budget_band.hard_veto_budget

    @property
    def source_url(self) -> str | None:
        return self.anchor.source_url


@dataclass(frozen=True)
class B014Anchor:
    """Typed B014 anchors: exclusive modes plus supporting diagnostics."""

    charged_rho: B014ModeAnchor
    neutral_rho: B014ModeAnchor
    omega: B014ModeAnchor
    hflav_rho: Anchor
    hflav_rho_omega: Anchor
    isospin_rho: Anchor
    xd_gamma: Anchor
    lhcb_ratio_rho0_kstar0: Anchor
    budget_source: str

    @property
    def modes(self) -> tuple[B014ModeAnchor, ...]:
        return (self.neutral_rho, self.charged_rho, self.omega)

    @property
    def value(self) -> float:
        return self.neutral_rho.value

    @property
    def reference_value(self) -> float:
        return self.neutral_rho.reference_value

    @property
    def budget(self) -> float:
        return self.neutral_rho.budget


@dataclass(frozen=True)
class B014ModeEvaluation:
    """Per-mode evaluated branching fraction and budget saturation."""

    mode: B014ModeAnchor
    predicted: float
    sm_prediction: float
    np_shift: float
    budget: float
    ratio: float
    passes: bool
    adapter_result: Any


def _validate_anchor_block_key(
    anchor: Anchor,
    *,
    process_id: str,
    label: str,
    candidates: tuple[str, ...],
) -> None:
    if anchor.block_key not in candidates:
        raise AnchorError(
            f"{process_id}: load_anchor selected {anchor.block_key!r} for {label}, "
            f"expected one of {candidates!r}"
        )


def _validate_branching_anchor(anchor: Anchor, *, process_id: str, label: str) -> None:
    if anchor.units != _EXPECTED_BRANCHING_UNITS:
        raise AnchorError(
            f"{process_id}: {label} anchor {anchor.block_key!r} must use units "
            f"{_EXPECTED_BRANCHING_UNITS!r}, got {anchor.units!r}"
        )
    if anchor.value <= 0.0 or not math.isfinite(anchor.value):
        raise AnchorError(
            f"{process_id}: {label} anchor {anchor.block_key!r} must have a "
            "positive finite value"
        )


def _validate_dimensionless_anchor(
    anchor: Anchor,
    *,
    process_id: str,
    label: str,
) -> None:
    if anchor.units != _EXPECTED_DIMENSIONLESS_UNITS:
        raise AnchorError(
            f"{process_id}: {label} anchor {anchor.block_key!r} must use units "
            f"{_EXPECTED_DIMENSIONLESS_UNITS!r}, got {anchor.units!r}"
        )
    if not math.isfinite(anchor.value):
        raise AnchorError(
            f"{process_id}: {label} anchor {anchor.block_key!r} must have a "
            "finite value"
        )


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number) or number <= 0.0:
        raise AnchorError(
            f"{process_id}: field {field_name!r} must be positive and finite"
        )
    return number


def _mode_uncertainty(
    anchor: Anchor,
    raw: Mapping[str, Any],
    *,
    process_id: str,
) -> float:
    if "uncertainty_plus" in raw or "uncertainty_minus" in raw:
        plus = _positive_float(
            raw.get("uncertainty_plus"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.uncertainty_plus",
        )
        minus = _positive_float(
            raw.get("uncertainty_minus"),
            process_id=process_id,
            field_name=f"{anchor.block_key}.uncertainty_minus",
        )
        return float(max(plus, minus))
    if anchor.uncertainty is None:
        raise AnchorError(
            f"{process_id}: {anchor.block_key!r} must provide a measurement "
            "uncertainty"
        )
    return _positive_float(
        anchor.uncertainty,
        process_id=process_id,
        field_name=f"{anchor.block_key}.uncertainty",
    )


def _build_budget_band(
    *,
    process_id: str,
    mode_key: str,
    anchor: Anchor,
    raw: Mapping[str, Any],
) -> B014ModeBudgetBand:
    sigma = _mode_uncertainty(anchor, raw, process_id=process_id)
    return B014ModeBudgetBand(
        source=f"{_BUDGET_SOURCE}; mode={mode_key}",
        mode_key=mode_key,
        experimental_sigma=sigma,
        hard_veto_budget=sigma,
        lower_total_edge=float(anchor.value - sigma),
        upper_total_edge=float(anchor.value + sigma),
    )


def _load_mode_anchor(
    process_id: str,
    *,
    pdg_block: Mapping[str, Any],
    candidates: tuple[str, ...],
    mode_key: str,
    mode_label: str,
    uncertainty_key: str = "uncertainty",
) -> B014ModeAnchor:
    anchor = load_anchor(
        process_id,
        family=_FAMILY,
        tier=_TIER,
        candidates=candidates,
        uncertainty_key=uncertainty_key,
    )
    _validate_anchor_block_key(
        anchor,
        process_id=process_id,
        label=mode_label,
        candidates=candidates,
    )
    _validate_branching_anchor(anchor, process_id=process_id, label=mode_label)
    raw = dict(pdg_block[anchor.block_key])
    return B014ModeAnchor(
        mode_key=mode_key,
        mode_label=mode_label,
        anchor=anchor,
        raw=raw,
        budget_band=_build_budget_band(
            process_id=process_id,
            mode_key=mode_key,
            anchor=anchor,
            raw=raw,
        ),
    )


def _load_b014_anchor(process_id: str) -> B014Anchor:
    pdg_block = load_pdg_block(process_id, family=_FAMILY, tier=_TIER)
    charged_rho = _load_mode_anchor(
        process_id,
        pdg_block=pdg_block,
        candidates=_CHARGED_RHO_CANDIDATES,
        mode_key="bp_to_rhop_gamma",
        mode_label="B+ -> rho+ gamma",
    )
    neutral_rho = _load_mode_anchor(
        process_id,
        pdg_block=pdg_block,
        candidates=_NEUTRAL_RHO_CANDIDATES,
        mode_key="b0_to_rho0_gamma",
        mode_label="B0 -> rho0 gamma",
    )
    omega = _load_mode_anchor(
        process_id,
        pdg_block=pdg_block,
        candidates=_OMEGA_CANDIDATES,
        mode_key="b0_to_omega_gamma",
        mode_label="B0 -> omega gamma",
        uncertainty_key="uncertainty_plus",
    )
    hflav_rho = load_anchor(
        process_id,
        family=_FAMILY,
        tier=_TIER,
        candidates=_HFLAV_RHO_CANDIDATES,
    )
    hflav_rho_omega = load_anchor(
        process_id,
        family=_FAMILY,
        tier=_TIER,
        candidates=_HFLAV_RHO_OMEGA_CANDIDATES,
    )
    isospin_rho = load_anchor(
        process_id,
        family=_FAMILY,
        tier=_TIER,
        candidates=_ISOSPIN_CANDIDATES,
    )
    xd_gamma = load_anchor(
        process_id,
        family=_FAMILY,
        tier=_TIER,
        candidates=_XD_GAMMA_CANDIDATES,
    )
    lhcb_ratio = load_anchor(
        process_id,
        family=_FAMILY,
        tier=_TIER,
        candidates=_LHCB_RATIO_CANDIDATES,
    )
    for label, anchor in (
        ("HFLAV B -> rho gamma", hflav_rho),
        ("HFLAV B -> rho/omega gamma", hflav_rho_omega),
        ("HFLAV B -> Xd gamma", xd_gamma),
    ):
        _validate_branching_anchor(anchor, process_id=process_id, label=label)
    for label, anchor in (
        ("HFLAV rho isospin", isospin_rho),
        ("LHCb rho0/Kstar0 ratio", lhcb_ratio),
    ):
        _validate_dimensionless_anchor(anchor, process_id=process_id, label=label)
    return B014Anchor(
        charged_rho=charged_rho,
        neutral_rho=neutral_rho,
        omega=omega,
        hflav_rho=hflav_rho,
        hflav_rho_omega=hflav_rho_omega,
        isospin_rho=isospin_rho,
        xd_gamma=xd_gamma,
        lhcb_ratio_rho0_kstar0=lhcb_ratio,
        budget_source=_BUDGET_SOURCE,
    )


def _budget_result(
    predicted: float,
    mode: B014ModeAnchor,
) -> tuple[float, float, float, bool]:
    np_shift = float(predicted - mode.reference_value)
    budget = float(mode.budget)
    ratio = abs(np_shift) / budget if budget > 0.0 else float("inf")
    return np_shift, budget, float(ratio), bool(ratio <= 1.0)


def _complex_mapping(mapping: Any) -> dict[str, complex]:
    return {str(key): complex(value) for key, value in dict(mapping).items()}


def _mode_budget_diagnostics(anchor: B014Anchor) -> dict[str, Any]:
    out: dict[str, Any] = {
        "budget_source": anchor.budget_source,
        "budget_policy": _BUDGET_POLICY,
        "budget_sm_theory_sigma": None,
        "sm_theory_prediction_available": False,
        "sm_theory_prediction_gap": _SM_THEORY_GAP_NOTE,
    }
    for mode in anchor.modes:
        prefix = mode.mode_key
        out.update(
            {
                f"{prefix}_experimental": float(mode.value),
                f"{prefix}_experimental_sigma": float(mode.uncertainty),
                f"{prefix}_reference_measurement": float(mode.reference_value),
                f"{prefix}_hard_veto_np_shift_budget": float(mode.budget),
                f"{prefix}_budget_lower_total_edge": float(
                    mode.budget_band.lower_total_edge
                ),
                f"{prefix}_budget_upper_total_edge": float(
                    mode.budget_band.upper_total_edge
                ),
                f"{prefix}_anchor_block": mode.anchor.block_key,
            }
        )
    return out


def _supporting_diagnostics(anchor: B014Anchor) -> dict[str, Any]:
    return {
        "hflav_rho_branching_fraction": float(anchor.hflav_rho.value),
        "hflav_rho_branching_fraction_uncertainty": float(
            anchor.hflav_rho.uncertainty
        ),
        "hflav_rho_omega_branching_fraction": float(anchor.hflav_rho_omega.value),
        "hflav_rho_omega_branching_fraction_uncertainty": float(
            anchor.hflav_rho_omega.uncertainty
        ),
        "hflav_isospin_rho_gamma": float(anchor.isospin_rho.value),
        "hflav_isospin_rho_gamma_uncertainty": float(anchor.isospin_rho.uncertainty),
        "hflav_isospin_rho_gamma_evaluated": False,
        "hflav_xd_gamma_branching_fraction": float(anchor.xd_gamma.value),
        "hflav_xd_gamma_branching_fraction_uncertainty": float(
            anchor.xd_gamma.uncertainty
        ),
        "lhcb_ratio_b0rho0_to_b0kstar0_gamma": float(
            anchor.lhcb_ratio_rho0_kstar0.value
        ),
        "supporting_observable_note": (
            "HFLAV aggregate, isospin, semi-inclusive X_d gamma, and LHCb "
            "rho0/K*0 ratio rows are loaded for provenance but not evaluated "
            "by this branching-fraction proxy."
        ),
    }


def _ckm_diagnostics() -> dict[str, Any]:
    ckm = bdgamma_ckm_factors()
    return {
        "ckm_input_bundle": ckm.input_bundle,
        "lambda_t_bd": complex(ckm.lambda_t_bd),
        "lambda_t_bs": complex(ckm.lambda_t_bs),
        "abs_vtd_over_vts": float(ckm.abs_vtd_over_vts),
        "abs_lambda_t_bd_over_lambda_t_bs": float(
            ckm.abs_lambda_t_bd_over_lambda_t_bs
        ),
        "ckm_power_suppression_vtd_over_vts_squared": float(
            ckm.ckm_power_suppression
        ),
    }


def _mode_evaluation_diagnostics(evaluations: tuple[B014ModeEvaluation, ...]) -> dict[str, Any]:
    return {
        "mode_predictions": {
            evaluation.mode.mode_key: {
                "mode_label": evaluation.mode.mode_label,
                "predicted_branching_fraction": float(evaluation.predicted),
                "normalization_proxy_branching_fraction": float(
                    evaluation.sm_prediction
                ),
                "experimental_branching_fraction": float(evaluation.mode.value),
                "np_shift_branching_fraction": float(evaluation.np_shift),
                "budget": float(evaluation.budget),
                "ratio": float(evaluation.ratio),
                "passes": bool(evaluation.passes),
            }
            for evaluation in evaluations
        },
        "max_mode_ratio": float(max(evaluation.ratio for evaluation in evaluations)),
        "all_mode_ratios_pass": all(evaluation.passes for evaluation in evaluations),
    }


@register
class Constraint:
    """Catalogued exclusive ``b -> d gamma`` branching-fraction constraint."""

    process_id = "B014"
    severity = Severity.HARD
    observable = "BR(B+ -> rho+ gamma), BR(B0 -> rho0/omega gamma)"

    def __init__(self) -> None:
        self.anchor = _load_b014_anchor(self.process_id)
        self.sm_inputs = bsgamma_default_sm_inputs()
        self.normalization_results = {
            mode.mode_key: exclusive_bdrhogamma_sm_branching_fraction(
                exclusive_sm_branching_fraction=mode.reference_value,
                inputs=self.sm_inputs,
                normalization_label=_NORMALIZATION_LABEL,
                normalization_source=mode.source_url,
                mode_label=mode.mode_label,
            )
            for mode in self.anchor.modes
        }

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(
                    self.normalization_results[
                        self.anchor.neutral_rho.mode_key
                    ].branching_fraction
                ),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; exclusive b -> d gamma "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    **_mode_budget_diagnostics(self.anchor),
                    **_supporting_diagnostics(self.anchor),
                    **_ckm_diagnostics(),
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "exclusive_normalization_branching_fraction": float(
                        self.anchor.reference_value
                    ),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )

        kk_mass = point.get_extra(_OPTIONAL_DIPOLE_MASS_EXTRA)
        try:
            evaluations = tuple(
                self._evaluate_mode(mode, couplings, kk_mass)
                for mode in self.anchor.modes
            )
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(
                    self.normalization_results[
                        self.anchor.neutral_rho.mode_key
                    ].branching_fraction
                ),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid quark_mass_basis_couplings for "
                    "the exclusive b -> d gamma C7 proxy"
                ),
                diagnostics={
                    **_mode_budget_diagnostics(self.anchor),
                    **_supporting_diagnostics(self.anchor),
                    **_ckm_diagnostics(),
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )

        worst = max(evaluations, key=lambda item: item.ratio)
        adapter_result = worst.adapter_result
        diagnostics = dict(adapter_result.diagnostics)
        wilsons = adapter_result.wilsons.wilsons if adapter_result.wilsons is not None else {}
        diagnostics.update(
            {
                **_mode_budget_diagnostics(self.anchor),
                **_supporting_diagnostics(self.anchor),
                **_mode_evaluation_diagnostics(evaluations),
                "evaluated": True,
                "dominant_mode": worst.mode.mode_key,
                "dominant_mode_label": worst.mode.mode_label,
                "normalization_anchor_policy": (
                    "B014.yaml has no separate theory-only branching-fraction "
                    "block; the PDG rows are used as C7 no-NP/exclusive "
                    "normalization proxies, not as theory predictions."
                ),
                "normalization_formula_branching_fraction": float(
                    adapter_result.sm_branching_fraction
                ),
                "normalization_formula_minus_reference_measurement": float(
                    adapter_result.sm_branching_fraction - worst.mode.reference_value
                ),
                "np_shift_branching_fraction": float(worst.np_shift),
                "ratio_to_sm_c7_power": float(adapter_result.ratio_to_sm),
                "c7_sm_eff": complex(adapter_result.c7_sm_eff),
                "c7p_sm_eff": complex(adapter_result.c7p_sm_eff),
                "c7_total": complex(adapter_result.c7_total),
                "c7p_total": complex(adapter_result.c7p_total),
                "c7_np": complex(adapter_result.c7_np),
                "c7p_np": complex(adapter_result.c7p_np),
                "wilson_coefficients": _complex_mapping(wilsons),
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": BSGAMMA_RS_MATCHING_ASSUMPTION_V1,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "kk_ew_mass_extra_used": kk_mass is not None,
                "down_sector_indices": (0, 2),
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=all(evaluation.passes for evaluation in evaluations),
            predicted=float(worst.predicted),
            sm_prediction=float(worst.sm_prediction),
            experimental=float(worst.mode.value),
            ratio=float(worst.ratio),
            budget=float(worst.budget),
            notes=(
                "BR(B+ -> rho+ gamma), BR(B0 -> rho0 gamma), and "
                "BR(B0 -> omega gamma) use the shared C7-normalized dipole "
                "formula with B014 YAML measurement normalization proxies. "
                "The RS contribution is a documented d-b C7/C7p proxy and is "
                "marked NEEDS-HUMAN-PHYSICS; the HARD budget is the largest "
                "per-mode measurement-consistency saturation, not "
                "SM-theory-vs-data room."
            ),
            diagnostics=diagnostics,
        )

    def _evaluate_mode(
        self,
        mode: B014ModeAnchor,
        couplings: Any,
        kk_mass: Any,
    ) -> B014ModeEvaluation:
        result = exclusive_bdrhogamma_from_couplings(
            couplings,
            exclusive_sm_branching_fraction=mode.reference_value,
            m_kk_gev=None if kk_mass is None else float(kk_mass),
            inputs=self.sm_inputs,
            normalization_label=_NORMALIZATION_LABEL,
            normalization_source=mode.source_url,
            mode_label=mode.mode_label,
        )
        predicted = float(result.branching_fraction)
        np_shift, budget, ratio, passes = _budget_result(predicted, mode)
        return B014ModeEvaluation(
            mode=mode,
            predicted=predicted,
            sm_prediction=float(result.sm_branching_fraction),
            np_shift=np_shift,
            budget=budget,
            ratio=ratio,
            passes=passes,
            adapter_result=result,
        )
