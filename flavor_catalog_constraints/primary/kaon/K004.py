"""K004 - charged rare kaon decay ``K+ -> pi+ nu nubar``.

Physics
-------
The SM short-distance prediction uses the Buras/Brod-Gorbahn-Stamou
parametrization

    BR = kappa_+ (1 + Delta_EM) [
           (Im(lambda_t X_t) / lambda^5)^2
           + (Re(lambda_c) P_c / lambda
              + Re(lambda_t X_t) / lambda^5)^2 ],

with the corresponding RS proxy entering as a complex shift
``lambda_t X_t -> lambda_t X_t + X_NP``.  The low-level formula and the
documented RS matching assumption live in ``quarkConstraints.rare_kaon_snd``
and are reached only through the
``flavor_catalog_constraints.physics_adapters.rare_kaon`` boundary.

Severity
--------
HARD.  The predicted total branching fraction is compared with the latest
NA62 K004 anchor using a direction-aware one-sigma budget that combines the
asymmetric experimental uncertainty and the catalogued Buras-Venturini SM
theory uncertainty in quadrature.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K004.yaml`` is the source of truth for the
NA62 branching-ratio anchor and SM reference values.  Numeric values below are
loaded and scale-converted from that sidecar, not hardcoded here.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.rare_kaon import (
    RARE_KAON_RS_MATCHING_ASSUMPTION_V1,
    kplus_piplus_nunu_from_couplings,
    rare_kaon_default_sm_inputs,
    rare_kaon_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_EW_MASS_EXTRA = "kk_ew_mass_gev"
_EXPERIMENTAL_ANCHOR_CANDIDATES = ("latest_experimental_value",)
_SM_ANCHOR_CANDIDATES = ("sm_prediction_buras_venturini_2022",)
_SM_VALIDATION_CANDIDATES = ("sm_prediction_brod_gorbahn_stamou_2021",)
_SCAFFOLD_UNCERTAINTY_KEY = "__k004_uncertainty_is_parsed_below__"
_EXPECTED_UNITS = "dimensionless branching fraction"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/kaon/K004.yaml latest_experimental_value "
    "+ sm_prediction_buras_venturini_2022"
)
_PARAMETRIZATION_CITATION = (
    "Buras et al. JHEP 11 (2015) 033, arXiv:1503.02693; "
    "Brod-Gorbahn-Stamou arXiv:2105.02868"
)


@dataclass(frozen=True)
class BranchingFractionAnchor:
    """Scaled branching-fraction value with symmetric or asymmetric errors."""

    block_key: str
    source: str | None
    year: int | None
    value: float
    uncertainty_upper: float
    uncertainty_lower: float
    scale: float
    units: str | None
    source_url: str | None
    snapshot_path: str | None
    value_summary: str | None

    @property
    def uncertainty(self) -> float:
        return float(0.5 * (self.uncertainty_upper + self.uncertainty_lower))


@dataclass(frozen=True)
class RareKaonBudgetBand:
    """Direction-aware K004 total-BR budget."""

    source: str
    central_residual: float
    experimental_sigma_upper: float
    experimental_sigma_lower: float
    sm_theory_sigma: float
    combined_sigma_upper: float
    combined_sigma_lower: float
    hard_veto_budget: float
    lower_edge: float
    upper_edge: float


@dataclass(frozen=True)
class K004Anchor:
    """Typed K004 anchor: NA62 experiment, SM references, and BR budget."""

    experimental: BranchingFractionAnchor
    standard_model: BranchingFractionAnchor
    validation_standard_model: BranchingFractionAnchor
    budget_band: RareKaonBudgetBand

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
            f"{process_id}: K004 anchor field {field_name!r}={value!r} is not numeric"
        ) from exc
    if not math.isfinite(number):
        raise AnchorError(
            f"{process_id}: K004 anchor field {field_name!r}={value!r} is not finite"
        )
    return number


def _positive_float(value: Any, *, process_id: str, field_name: str) -> float:
    number = _required_float(value, process_id=process_id, field_name=field_name)
    if number <= 0.0:
        raise AnchorError(
            f"{process_id}: K004 anchor field {field_name!r} must be positive"
        )
    return number


_ASYMMETRIC_UNCERTAINTY_RE = re.compile(
    r"^\s*\+?(?P<upper>[0-9.eE+-]+)\s*/\s*-(?P<lower>[0-9.eE+-]+)\s*$"
)


def _parse_uncertainty_pair(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> tuple[float, float]:
    if isinstance(value, str):
        match = _ASYMMETRIC_UNCERTAINTY_RE.match(value)
        if match is not None:
            upper = _positive_float(
                match.group("upper"),
                process_id=process_id,
                field_name=f"{field_name}.upper",
            )
            lower = _positive_float(
                match.group("lower"),
                process_id=process_id,
                field_name=f"{field_name}.lower",
            )
            return upper, lower
    sigma = _positive_float(value, process_id=process_id, field_name=field_name)
    return sigma, sigma


def _pdg_subblock_for_anchor(
    anchor: Anchor,
    *,
    process_id: str,
) -> Mapping[str, Any]:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg_block = data.get("pdg_or_equivalent")
    if not isinstance(pdg_block, Mapping):
        raise AnchorError(
            f"{process_id}: 'pdg_or_equivalent' is not a mapping while loading "
            f"{anchor.block_key}"
        )
    sub = pdg_block.get(anchor.block_key)
    if not isinstance(sub, Mapping):
        present = sorted(str(key) for key in pdg_block)
        raise AnchorError(
            f"{process_id}: load_anchor selected {anchor.block_key!r}, but that "
            f"block is not available as a mapping (present keys: {present})"
        )
    return sub


def _load_branching_anchor(
    candidates: tuple[str, ...],
    *,
    process_id: str,
) -> BranchingFractionAnchor:
    scaffold_anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=candidates,
        # K004 stores asymmetric experimental errors as "+up/-down" strings.
        # The scaffold selects and types the block; K004 parses that richer
        # uncertainty shape after selection.
        uncertainty_key=_SCAFFOLD_UNCERTAINTY_KEY,
    )
    sub = _pdg_subblock_for_anchor(scaffold_anchor, process_id=process_id)
    block_key = scaffold_anchor.block_key
    units = scaffold_anchor.units
    if units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_UNITS!r} for {block_key}, "
            f"got {units!r}"
        )
    scale = _positive_float(
        sub.get("scale", 1.0),
        process_id=process_id,
        field_name=f"{block_key}.scale",
    )
    value = _required_float(
        scaffold_anchor.value,
        process_id=process_id,
        field_name=f"{block_key}.value",
    )
    uncertainty_upper, uncertainty_lower = _parse_uncertainty_pair(
        sub.get("uncertainty"),
        process_id=process_id,
        field_name=f"{block_key}.uncertainty",
    )
    return BranchingFractionAnchor(
        block_key=block_key,
        source=_optional_str(scaffold_anchor.source),
        year=scaffold_anchor.year,
        value=float(value * scale),
        uncertainty_upper=float(uncertainty_upper * scale),
        uncertainty_lower=float(uncertainty_lower * scale),
        scale=float(scale),
        units=units,
        source_url=_optional_str(scaffold_anchor.source_url),
        snapshot_path=_optional_str(scaffold_anchor.snapshot_path),
        value_summary=_optional_str(sub.get("value_summary")),
    )


def _build_budget_band(
    *,
    experimental: BranchingFractionAnchor,
    standard_model: BranchingFractionAnchor,
) -> RareKaonBudgetBand:
    combined_upper = math.sqrt(
        experimental.uncertainty_upper**2 + standard_model.uncertainty_upper**2
    )
    combined_lower = math.sqrt(
        experimental.uncertainty_lower**2 + standard_model.uncertainty_lower**2
    )
    lower_edge = experimental.value - combined_lower
    upper_edge = experimental.value + combined_upper
    hard_budget = max(combined_upper, combined_lower)
    if hard_budget <= 0.0 or lower_edge <= 0.0:
        raise AnchorError("K004: constructed branching-ratio budget is invalid")
    return RareKaonBudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(experimental.value - standard_model.value),
        experimental_sigma_upper=float(experimental.uncertainty_upper),
        experimental_sigma_lower=float(experimental.uncertainty_lower),
        sm_theory_sigma=float(standard_model.uncertainty),
        combined_sigma_upper=float(combined_upper),
        combined_sigma_lower=float(combined_lower),
        hard_veto_budget=float(hard_budget),
        lower_edge=float(lower_edge),
        upper_edge=float(upper_edge),
    )


def _load_k004_anchor(process_id: str) -> K004Anchor:
    experimental = _load_branching_anchor(
        _EXPERIMENTAL_ANCHOR_CANDIDATES,
        process_id=process_id,
    )
    standard_model = _load_branching_anchor(
        _SM_ANCHOR_CANDIDATES,
        process_id=process_id,
    )
    validation_standard_model = _load_branching_anchor(
        _SM_VALIDATION_CANDIDATES,
        process_id=process_id,
    )
    return K004Anchor(
        experimental=experimental,
        standard_model=standard_model,
        validation_standard_model=validation_standard_model,
        budget_band=_build_budget_band(
            experimental=experimental,
            standard_model=standard_model,
        ),
    )


def _selected_budget(
    predicted: float,
    anchor: K004Anchor,
) -> tuple[float, float, bool]:
    pull = float(predicted - anchor.value)
    budget = (
        anchor.budget_band.combined_sigma_upper
        if pull >= 0.0
        else anchor.budget_band.combined_sigma_lower
    )
    ratio = abs(pull) / budget if budget > 0.0 else float("inf")
    return float(budget), float(ratio), bool(ratio <= 1.0)


@register
class Constraint:
    """Catalogued charged rare-kaon branching-ratio constraint (K004)."""

    process_id = "K004"
    severity = Severity.HARD
    observable = "BR(K+ -> pi+ nu nubar)"

    def __init__(self) -> None:
        self.anchor = _load_k004_anchor(self.process_id)
        self.sm_inputs = rare_kaon_default_sm_inputs()
        self.sm_result = rare_kaon_sm_branching_fraction(self.sm_inputs)

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
                    f"extra {_REQUIRED_EXTRA!r} absent; K+ -> pi+ nu nubar "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    "missing_extra": _REQUIRED_EXTRA,
                    "sm_anchor_branching_fraction": float(self.anchor.sm_value),
                    "sm_formula_branching_fraction": float(
                        self.sm_result.branching_fraction
                    ),
                    "budget_source": self.anchor.budget_band.source,
                },
            )

        kk_ew_mass = point.get_extra(_OPTIONAL_EW_MASS_EXTRA)
        result = kplus_piplus_nunu_from_couplings(
            couplings,
            m_kk_gev=None if kk_ew_mass is None else float(kk_ew_mass),
            inputs=self.sm_inputs,
        )
        predicted = float(result.branching_fraction)
        budget, ratio, passes = _selected_budget(predicted, self.anchor)

        diagnostics = dict(result.diagnostics)
        diagnostics.update(
            {
                "sm_anchor_branching_fraction": float(self.anchor.sm_value),
                "sm_formula_branching_fraction": float(
                    result.sm_branching_fraction
                ),
                "sm_formula_minus_anchor": float(
                    result.sm_branching_fraction - self.anchor.sm_value
                ),
                "validation_sm_bgs2021_branching_fraction": float(
                    self.anchor.validation_standard_model.value
                ),
                "experimental_block": self.anchor.experimental.block_key,
                "sm_block": self.anchor.standard_model.block_key,
                "experimental_sigma_upper": float(
                    self.anchor.budget_band.experimental_sigma_upper
                ),
                "experimental_sigma_lower": float(
                    self.anchor.budget_band.experimental_sigma_lower
                ),
                "sm_theory_sigma": float(self.anchor.budget_band.sm_theory_sigma),
                "budget_combined_sigma_upper": float(
                    self.anchor.budget_band.combined_sigma_upper
                ),
                "budget_combined_sigma_lower": float(
                    self.anchor.budget_band.combined_sigma_lower
                ),
                "budget_lower_edge": float(self.anchor.budget_band.lower_edge),
                "budget_upper_edge": float(self.anchor.budget_band.upper_edge),
                "budget_source": self.anchor.budget_band.source,
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": RARE_KAON_RS_MATCHING_ASSUMPTION_V1,
                "needs_human_physics": (
                    "NEEDS-HUMAN-PHYSICS: full RS electroweak KK/Z/Z' tower "
                    "and neutrino-coupling matching are not available on "
                    "ParameterPoint; v1 uses the documented Z-like proxy."
                ),
                "kk_ew_mass_extra_used": kk_ew_mass is not None,
                "kappa_plus": float(result.kappa_plus),
                "p_c": float(result.p_c),
                "x_t": float(result.x_t),
                "lambda_wolfenstein": float(result.lambda_wolfenstein),
                "lambda_c": complex(result.lambda_c),
                "lambda_t": complex(result.lambda_t),
                "x_eff_top": complex(result.x_eff_top),
                "x_np_total": complex(result.x_np_total),
                "np_shift_branching_fraction": float(
                    result.np_shift_branching_fraction
                ),
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
                "BR(K+ -> pi+ nu nubar) uses the Buras/Brod-Gorbahn-Stamou "
                "short-distance parametrization. RS contribution is a "
                "documented Z-like X-function shift from mass-basis s-d "
                "couplings; full EW/lepton matching is marked "
                "NEEDS-HUMAN-PHYSICS. Budget combines NA62 and SM theory "
                "uncertainties in quadrature."
            ),
            diagnostics=diagnostics,
        )
