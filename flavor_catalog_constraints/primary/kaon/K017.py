"""K017 - charged-kaon leptonic LFU ratio ``R_K``.

Physics
-------
``R_K = Gamma(K+ -> e+ nu_e(gamma)) / Gamma(K+ -> mu+ nu_mu(gamma))`` is
helicity-suppressed and clean because the common ``f_K``, ``V_us``, and kaon
normalization cancel.  The SM ratio is evaluated with the reusable
charged-leptonic tree kernel in
``flavor_catalog_constraints.physics_adapters.leptonic_tree``; the
process-specific radiative multiplier is fixed by the K017 YAML
Cirigliano-Rosell SM anchor.

RS matching
-----------
For points carrying ``rs_charged_current``, K017 consumes the stored
minimal-left-handed epsilons and evaluates

    R_K = R_K^SM |1 + epsilon_us^e|^2 / |1 + epsilon_us^mu|^2.

Charged-Higgs, heavy-neutrino, scalar/RH-current, and radiative nonminimal
effects remain outside this minimal vector matching.  No mass proxy is used.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K017.yaml`` is the source of truth for the
PDG average, Cirigliano-Rosell SM value, uncertainties, scale factors, and
provenance.  Numeric anchor values below are loaded and scale-converted from
that sidecar, not hardcoded in this constraint.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.charged_current import (
    CHARGED_CURRENT_NONMINIMAL_NEEDS_HUMAN,
    charged_current_epsilon,
    charged_current_source_diagnostics,
    shifted_lfu_ratio,
)
from flavor_catalog_constraints.physics_adapters.leptonic_tree import (
    leptonic_tree_kplus_enu_over_munu_inputs_from_sm_ratio_anchor,
    leptonic_tree_sm_lfu_ratio,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_REQUIRED_EXTRA = "rs_charged_current"
_EXPERIMENTAL_ANCHOR_CANDIDATES = ("canonical_ratio",)
_SM_ANCHOR_CANDIDATES = ("sm_prediction",)
_DOMINANT_INPUT_CANDIDATES = ("dominant_experimental_input",)
_SUPPORTING_INPUT_CANDIDATES = ("supporting_experimental_input_kloe",)
_EXPECTED_UNITS = "dimensionless decay-rate ratio"
_BUDGET_SOURCE = "flavor_catalog/processes/kaon/K017.yaml canonical_ratio + sm_prediction"
_PARAMETRIZATION_CITATION = (
    "charged-pseudoscalar leptonic tree kernel with K017.yaml "
    "Cirigliano-Rosell radiative SM ratio anchor; minimal-LH epsilon_us "
    "charged-current ratio from rs_charged_current"
)
_NONMINIMAL_STATUS = (
    "NEEDS-HUMAN-PHYSICS: charged-Higgs, heavy-neutrino, scalar/RH-current, "
    "and radiative nonminimal effects are not included in the minimal "
    "left-handed epsilon_us LFU ratio."
)


@dataclass(frozen=True)
class ScaledRatioAnchor:
    """Dimensionless ratio anchor converted from YAML display units."""

    anchor: Anchor
    value: float
    uncertainty: float
    scale: float

    @property
    def block_key(self) -> str:
        return self.anchor.block_key

    @property
    def source(self) -> str | None:
        return self.anchor.source

    @property
    def source_url(self) -> str | None:
        return self.anchor.source_url

    @property
    def snapshot_path(self) -> str | None:
        return self.anchor.snapshot_path


@dataclass(frozen=True)
class ExperimentalInputSummary:
    """Compact NA62/KLOE provenance row."""

    block_key: str
    source: str | None
    year: int | None
    value: float
    stat_uncertainty: float | None
    syst_uncertainty: float | None
    combined_uncertainty: float | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class K017BudgetBand:
    """Uncertainty-aware K017 NP-shift budget in ratio units."""

    source: str
    central_residual: float
    experimental_sigma: float
    sm_theory_sigma: float
    combined_sigma: float
    hard_veto_budget: float
    lower_total_edge: float
    upper_total_edge: float


@dataclass(frozen=True)
class K017Anchor:
    """Typed K017 anchor: PDG ratio, SM prediction, inputs, and budget."""

    experimental: ScaledRatioAnchor
    standard_model: ScaledRatioAnchor
    dominant_experimental_input: ExperimentalInputSummary
    supporting_experimental_input: ExperimentalInputSummary
    budget_band: K017BudgetBand

    @property
    def value(self) -> float:
        return self.experimental.value

    @property
    def uncertainty(self) -> float:
        return self.experimental.uncertainty

    @property
    def sm_value(self) -> float:
        return self.standard_model.value

    @property
    def budget(self) -> float:
        return self.budget_band.hard_veto_budget


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: K017 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: K017 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(
            f"{process_id}: K017 anchor field {field_name!r} must be positive"
        )
    return number


def _pdg_subblock_for_anchor(
    anchor: Anchor,
    *,
    process_id: str,
) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg_block = data.get("pdg_or_equivalent")
    if not isinstance(pdg_block, Mapping):
        raise AnchorError(f"{process_id}: 'pdg_or_equivalent' is not a mapping")
    sub = pdg_block.get(anchor.block_key)
    if not isinstance(sub, Mapping):
        present = sorted(str(key) for key in pdg_block)
        raise AnchorError(
            f"{process_id}: load_anchor selected {anchor.block_key!r}, but that "
            f"block is not available as a mapping (present keys: {present})"
        )
    return sub


def _load_scaled_ratio_anchor(
    process_id: str,
    candidates: tuple[str, ...],
) -> ScaledRatioAnchor:
    anchor = load_anchor(process_id, family=_FAMILY, candidates=candidates)
    if anchor.uncertainty is None:
        raise AnchorError(
            f"{process_id}: anchor {anchor.block_key!r} must provide uncertainty"
        )
    if anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: anchor {anchor.block_key!r} has unsupported units "
            f"{anchor.units!r}"
        )
    if anchor.value <= 0.0 or anchor.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: anchor {anchor.block_key!r} value and uncertainty "
            "must be positive"
        )
    sub = _pdg_subblock_for_anchor(anchor, process_id=process_id)
    scale = _positive_float(
        sub.get("scale", 1.0),
        process_id=process_id,
        field_name=f"{anchor.block_key}.scale",
    )
    return ScaledRatioAnchor(
        anchor=anchor,
        value=float(anchor.value * scale),
        uncertainty=float(anchor.uncertainty * scale),
        scale=float(scale),
    )


def _optional_scaled_float(
    value: Any,
    *,
    scale: float,
    process_id: str,
    field_name: str,
) -> float | None:
    if value is None:
        return None
    return float(
        _positive_float(value, process_id=process_id, field_name=field_name) * scale
    )


def _load_experimental_input_summary(
    process_id: str,
    candidates: tuple[str, ...],
) -> ExperimentalInputSummary:
    anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=candidates,
        uncertainty_key="__k017_input_uncertainty_parsed_below__",
    )
    if anchor.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: input {anchor.block_key!r} has unsupported units "
            f"{anchor.units!r}"
        )
    sub = _pdg_subblock_for_anchor(anchor, process_id=process_id)
    scale = _positive_float(
        sub.get("scale", 1.0),
        process_id=process_id,
        field_name=f"{anchor.block_key}.scale",
    )
    stat = _optional_scaled_float(
        sub.get("stat_uncertainty"),
        scale=scale,
        process_id=process_id,
        field_name=f"{anchor.block_key}.stat_uncertainty",
    )
    syst = _optional_scaled_float(
        sub.get("syst_uncertainty"),
        scale=scale,
        process_id=process_id,
        field_name=f"{anchor.block_key}.syst_uncertainty",
    )
    combined = None
    if stat is not None and syst is not None:
        combined = float(math.sqrt(stat * stat + syst * syst))
    return ExperimentalInputSummary(
        block_key=anchor.block_key,
        source=_optional_str(anchor.source),
        year=anchor.year,
        value=float(anchor.value * scale),
        stat_uncertainty=stat,
        syst_uncertainty=syst,
        combined_uncertainty=combined,
        source_url=_optional_str(anchor.source_url),
        snapshot_path=_optional_str(anchor.snapshot_path),
    )


def _build_budget_band(
    *,
    experimental: ScaledRatioAnchor,
    standard_model: ScaledRatioAnchor,
) -> K017BudgetBand:
    combined_sigma = math.sqrt(
        experimental.uncertainty**2 + standard_model.uncertainty**2
    )
    central_residual = abs(experimental.value - standard_model.value)
    budget = central_residual + combined_sigma
    if budget <= 0.0 or not math.isfinite(budget):
        raise AnchorError("K017: constructed LFU-ratio budget is invalid")
    return K017BudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(central_residual),
        experimental_sigma=float(experimental.uncertainty),
        sm_theory_sigma=float(standard_model.uncertainty),
        combined_sigma=float(combined_sigma),
        hard_veto_budget=float(budget),
        lower_total_edge=float(standard_model.value - budget),
        upper_total_edge=float(standard_model.value + budget),
    )


def _load_k017_anchor(process_id: str) -> K017Anchor:
    experimental = _load_scaled_ratio_anchor(
        process_id,
        _EXPERIMENTAL_ANCHOR_CANDIDATES,
    )
    standard_model = _load_scaled_ratio_anchor(
        process_id,
        _SM_ANCHOR_CANDIDATES,
    )
    return K017Anchor(
        experimental=experimental,
        standard_model=standard_model,
        dominant_experimental_input=_load_experimental_input_summary(
            process_id,
            _DOMINANT_INPUT_CANDIDATES,
        ),
        supporting_experimental_input=_load_experimental_input_summary(
            process_id,
            _SUPPORTING_INPUT_CANDIDATES,
        ),
        budget_band=_build_budget_band(
            experimental=experimental,
            standard_model=standard_model,
        ),
    )


def _theory_input_citation(anchor: K017Anchor) -> str:
    return (
        "K017.yaml SM anchor: "
        f"{anchor.standard_model.block_key} R_K^SM={anchor.sm_value:.6g} "
        f"({anchor.standard_model.source}); radiative multiplier is normalized "
        "so the reusable helicity-suppressed tree ratio reproduces this anchor."
    )


@register
class Constraint:
    """Catalogued charged-kaon leptonic LFU-ratio constraint (K017)."""

    process_id = "K017"
    severity = Severity.HARD
    observable = "R_K = Gamma(K -> e nu) / Gamma(K -> mu nu)"

    def __init__(self) -> None:
        self.anchor = _load_k017_anchor(self.process_id)
        self.ratio_inputs = leptonic_tree_kplus_enu_over_munu_inputs_from_sm_ratio_anchor(
            sm_ratio=self.anchor.sm_value,
            constants_citation=_theory_input_citation(self.anchor),
        )
        self.sm_result = leptonic_tree_sm_lfu_ratio(self.ratio_inputs)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = self.sm_result
        common_diagnostics: dict[str, Any] = {
            "sm_anchor_ratio": float(self.anchor.sm_value),
            "sm_formula_ratio": float(result.sm_ratio),
            "sm_formula_minus_anchor": float(result.sm_ratio - self.anchor.sm_value),
            "tree_ratio_without_radiation": float(
                result.diagnostics["tree_ratio_without_radiation"]
            ),
            "radiative_correction_multiplier": float(
                result.diagnostics["radiative_correction_multiplier"]
            ),
            "experimental_block": self.anchor.experimental.block_key,
            "sm_block": self.anchor.standard_model.block_key,
            "dominant_experimental_input_block": (
                self.anchor.dominant_experimental_input.block_key
            ),
            "supporting_experimental_input_block": (
                self.anchor.supporting_experimental_input.block_key
            ),
            "dominant_experimental_input_ratio": float(
                self.anchor.dominant_experimental_input.value
            ),
            "supporting_experimental_input_ratio": float(
                self.anchor.supporting_experimental_input.value
            ),
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
                "|R_total - R_SM(formula)| compared with "
                "|R_exp - R_SM(anchor)| + sqrt(sigma_exp^2 + sigma_SM^2)"
            ),
            "parametrization_citation": _PARAMETRIZATION_CITATION,
            "nonminimal_charged_current_status": _NONMINIMAL_STATUS,
            "needs_human_physics": CHARGED_CURRENT_NONMINIMAL_NEEDS_HUMAN,
        }

        charged = point.get_extra(_REQUIRED_EXTRA)
        if charged is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                predicted=None,
                sm_prediction=float(result.sm_ratio),
                experimental=float(self.anchor.value),
                ratio=None,
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; R_K charged-current "
                    "epsilon shift was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "np_shift_ratio": 0.0,
                    **common_diagnostics,
                },
            )

        try:
            eps_e = charged_current_epsilon(
                charged,
                up="u",
                down="s",
                lepton="e",
            )
            eps_mu = charged_current_epsilon(
                charged,
                up="u",
                down="s",
                lepton="mu",
            )
            predicted = shifted_lfu_ratio(float(result.sm_ratio), eps_e, eps_mu)
            source_diag = charged_current_source_diagnostics(charged)
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                predicted=float(self.sm_result.ratio),
                sm_prediction=float(self.sm_result.sm_ratio),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_charged_current extra for K017 "
                    "epsilon_us LFU ratio."
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "np_shift_ratio": 0.0,
                    **common_diagnostics,
                },
            )

        np_shift = float(predicted - result.sm_ratio)
        budget = float(self.anchor.budget)
        ratio = abs(np_shift) / budget if budget > 0.0 else float("inf")
        diagnostics = dict(common_diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                "epsilon_us_e": eps_e.epsilon,
                "epsilon_us_mu": eps_mu.epsilon,
                "numerator_rate_multiplier": float(eps_e.rate_multiplier),
                "denominator_rate_multiplier": float(eps_mu.rate_multiplier),
                "np_shift_ratio": float(np_shift),
                "rs_matching_formula": (
                    "R_K = R_K^SM |1+epsilon_us^e|^2/"
                    "|1+epsilon_us^mu|^2"
                ),
                **source_diag,
            }
        )

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(ratio <= 1.0),
            predicted=predicted,
            sm_prediction=float(result.sm_ratio),
            experimental=float(self.anchor.value),
            ratio=float(ratio),
            budget=budget,
            notes=(
                "R_K uses the charged-leptonic tree helicity ratio with the "
                "K017 YAML SM radiative normalization and minimal-LH "
                "rs_charged_current epsilons. The HARD ratio is the vector "
                "LFU NP shift over the PDG-vs-SM YAML budget."
            ),
            diagnostics=diagnostics,
        )
