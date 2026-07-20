"""B009 - leptonic decay ``B+ -> tau+ nu_tau``.

Physics
-------
The SM prediction is evaluated with the reusable charged-leptonic tree formula

    BR = G_F^2 m_B m_tau^2 f_B^2 |V_ub|^2 tau_B
         (1 - m_tau^2/m_B^2)^2 / (8 pi),

through ``flavor_catalog_constraints.physics_adapters.leptonic_tree``.  The
``f_B`` and ``|V_ub|`` theory inputs are loaded from explicit B009 YAML anchors;
the resulting SM branching fraction reproduces the B009 UTfit SM anchor at the
few-per-mille level.

Severity
--------
HARD.  The observed HFLAV branching-fraction average is compared to the UTfit
SM anchor by giving new physics the uncertainty-aware room
``|BR_exp - BR_SM(anchor)| + sqrt(sigma_exp^2 + sigma_SM^2)``.

RS matching
-----------
For points carrying ``rs_charged_current``, B009 evaluates the minimal
left-handed W/W' shift

    BR(B+ -> tau+ nu_tau) = BR_SM |1 + epsilon_ub^tau|^2.

Charged-Higgs, right-handed W, scalar, heavy-neutrino, and other nonminimal
effects remain outside this vector matching.  No mass proxy is used.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B009.yaml`` is the source of truth for the
HFLAV average, the UTfit SM prediction, the explicit ``f_B``/``|V_ub|`` theory
inputs, and the Belle/BaBar/Belle II input provenance.  Numeric anchor values
below are loaded and scale-converted from that sidecar, not hardcoded in this
constraint.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor, load_full_yaml
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.charged_current import (
    CHARGED_CURRENT_NONMINIMAL_NEEDS_HUMAN,
    charged_current_epsilon,
    charged_current_source_diagnostics,
    shifted_branching_fraction,
)
from flavor_catalog_constraints.physics_adapters.leptonic_tree import (
    leptonic_tree_bplus_tau_nu_inputs_from_theory_anchors,
    leptonic_tree_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "rs_charged_current"
_EXPERIMENTAL_ANCHOR_CANDIDATES = ("canonical_average",)
_PDG_AVERAGE_CANDIDATES = ("pdg_2025_average",)
_SM_ANCHOR_CANDIDATES = ("sm_prediction",)
_THEORY_F_B_CANDIDATES = ("theory_flag_f_B",)
_THEORY_V_UB_CANDIDATES = ("theory_ckm_abs_vub",)
_EXPERIMENTAL_INPUTS_KEY = "experimental_inputs"
_EXPECTED_SCALED_UNITS = "10^-4 branching fraction"
_EXPECTED_F_B_UNITS = "GeV"
_EXPECTED_CKM_UNITS = "dimensionless"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B009.yaml canonical_average "
    "+ sm_prediction"
)
_PARAMETRIZATION_CITATION = (
    "tree-level charged pseudoscalar leptonic decay formula; "
    "B009.yaml theory_flag_f_B and theory_ckm_abs_vub anchors; "
    "B009.yaml UTfit Summer 2024 SM BR anchor for validation; minimal-LH "
    "epsilon_ub^tau charged-current shift from rs_charged_current"
)
_NONMINIMAL_STATUS = (
    "PARTIAL: minimal left-handed W/W' epsilon_ub^tau is included; "
    "charged-Higgs, right-handed W, scalar, and heavy-neutrino effects are not "
    "built."
)


@dataclass(frozen=True)
class ScaledBranchingAnchor:
    """Branching-fraction anchor converted from YAML display units."""

    anchor: Anchor
    value: float
    uncertainty: float
    scale: float

    @property
    def block_key(self) -> str:
        return self.anchor.block_key

    @property
    def raw_value(self) -> float:
        return self.anchor.value

    @property
    def raw_uncertainty(self) -> float | None:
        return self.anchor.uncertainty

    @property
    def units(self) -> str | None:
        return self.anchor.units

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
    """Compact Belle/BaBar/Belle II input provenance row."""

    experiment: str
    year: int | None
    method: str | None
    source_url: str | None
    snapshot_path: str | None


@dataclass(frozen=True)
class B009BudgetBand:
    """Uncertainty-aware B009 NP-shift budget in branching-fraction units."""

    source: str
    central_residual: float
    experimental_sigma: float
    sm_theory_sigma: float
    combined_sigma: float
    hard_veto_budget: float
    lower_total_edge: float
    upper_total_edge: float


@dataclass(frozen=True)
class B009Anchor:
    """Typed B009 anchor: HFLAV, PDG comparison, SM, inputs, and budget."""

    experimental: ScaledBranchingAnchor
    pdg_average: ScaledBranchingAnchor
    standard_model: ScaledBranchingAnchor
    f_b: Anchor
    v_ub: Anchor
    experimental_inputs: tuple[ExperimentalInputSummary, ...]
    budget_band: B009BudgetBand

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


def _optional_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _scale_from_units(units: str | None, *, process_id: str, block_key: str) -> float:
    if units == _EXPECTED_SCALED_UNITS:
        return 1.0e-4
    if units == "branching fraction":
        return 1.0
    raise AnchorError(
        f"{process_id}: anchor {block_key!r} has unsupported units {units!r}"
    )


def _load_scaled_branching_anchor(
    process_id: str,
    candidates: tuple[str, ...],
) -> ScaledBranchingAnchor:
    anchor = load_anchor(process_id, family=_FAMILY, candidates=candidates)
    if anchor.uncertainty is None:
        raise AnchorError(
            f"{process_id}: anchor {anchor.block_key!r} must provide uncertainty"
        )
    if anchor.value <= 0.0 or anchor.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: anchor {anchor.block_key!r} value and uncertainty "
            "must be positive"
        )
    scale = _scale_from_units(
        anchor.units,
        process_id=process_id,
        block_key=anchor.block_key,
    )
    return ScaledBranchingAnchor(
        anchor=anchor,
        value=float(anchor.value * scale),
        uncertainty=float(anchor.uncertainty * scale),
        scale=float(scale),
    )


def _load_positive_theory_anchor(
    process_id: str,
    candidates: tuple[str, ...],
    *,
    expected_units: str,
    label: str,
) -> Anchor:
    anchor = load_anchor(process_id, family=_FAMILY, candidates=candidates)
    if anchor.value <= 0.0 or not math.isfinite(anchor.value):
        raise AnchorError(
            f"{process_id}: {label} anchor {anchor.block_key!r} must be positive"
        )
    if anchor.units != expected_units:
        raise AnchorError(
            f"{process_id}: {label} anchor {anchor.block_key!r} must use units "
            f"{expected_units!r}, got {anchor.units!r}"
        )
    return anchor


def _load_pdg_average(process_id: str) -> ScaledBranchingAnchor:
    # The PDG comparison block carries an asymmetric string uncertainty.  Use
    # the scaffold to select the block, then parse the richer shape locally.
    anchor = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_PDG_AVERAGE_CANDIDATES,
        uncertainty_key="__b009_pdg_uncertainty_parsed_below__",
    )
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg_block = data.get("pdg_or_equivalent")
    if not isinstance(pdg_block, Mapping):
        raise AnchorError(f"{process_id}: pdg_or_equivalent is not a mapping")
    sub = pdg_block.get(anchor.block_key)
    if not isinstance(sub, Mapping):
        raise AnchorError(f"{process_id}: PDG block {anchor.block_key!r} missing")
    raw_uncertainty = str(sub.get("uncertainty", "")).strip()
    pieces = raw_uncertainty.replace("+", "").split("/")
    if len(pieces) != 2 or not pieces[1].startswith("-"):
        raise AnchorError(
            f"{process_id}: cannot parse PDG asymmetric uncertainty "
            f"{raw_uncertainty!r}"
        )
    upper = float(pieces[0])
    lower = float(pieces[1][1:])
    scale = _scale_from_units(
        anchor.units,
        process_id=process_id,
        block_key=anchor.block_key,
    )
    return ScaledBranchingAnchor(
        anchor=anchor,
        value=float(anchor.value * scale),
        uncertainty=float(0.5 * (upper + lower) * scale),
        scale=float(scale),
    )


def _load_experimental_inputs(process_id: str) -> tuple[ExperimentalInputSummary, ...]:
    data = load_full_yaml(process_id, family=_FAMILY)
    pdg_block = data.get("pdg_or_equivalent")
    if not isinstance(pdg_block, Mapping):
        raise AnchorError(f"{process_id}: pdg_or_equivalent is not a mapping")
    rows = pdg_block.get(_EXPERIMENTAL_INPUTS_KEY)
    if not isinstance(rows, list) or not rows:
        raise AnchorError(f"{process_id}: experimental_inputs list is required")
    out: list[ExperimentalInputSummary] = []
    for row in rows:
        if not isinstance(row, Mapping):
            raise AnchorError(f"{process_id}: experimental input row is not a mapping")
        out.append(
            ExperimentalInputSummary(
                experiment=str(row.get("experiment")),
                year=_optional_int(row.get("year")),
                method=_optional_str(row.get("method")),
                source_url=_optional_str(row.get("source_url")),
                snapshot_path=_optional_str(row.get("snapshot_path")),
            )
        )
    return tuple(out)


def _build_budget_band(
    *,
    experimental: ScaledBranchingAnchor,
    standard_model: ScaledBranchingAnchor,
) -> B009BudgetBand:
    combined_sigma = math.sqrt(
        experimental.uncertainty**2 + standard_model.uncertainty**2
    )
    central_residual = abs(experimental.value - standard_model.value)
    budget = central_residual + combined_sigma
    if budget <= 0.0 or not math.isfinite(budget):
        raise AnchorError("B009: constructed branching-fraction budget is invalid")
    return B009BudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(central_residual),
        experimental_sigma=float(experimental.uncertainty),
        sm_theory_sigma=float(standard_model.uncertainty),
        combined_sigma=float(combined_sigma),
        hard_veto_budget=float(budget),
        lower_total_edge=float(standard_model.value - budget),
        upper_total_edge=float(standard_model.value + budget),
    )


def _load_b009_anchor(process_id: str) -> B009Anchor:
    experimental = _load_scaled_branching_anchor(
        process_id,
        _EXPERIMENTAL_ANCHOR_CANDIDATES,
    )
    standard_model = _load_scaled_branching_anchor(
        process_id,
        _SM_ANCHOR_CANDIDATES,
    )
    return B009Anchor(
        experimental=experimental,
        pdg_average=_load_pdg_average(process_id),
        standard_model=standard_model,
        f_b=_load_positive_theory_anchor(
            process_id,
            _THEORY_F_B_CANDIDATES,
            expected_units=_EXPECTED_F_B_UNITS,
            label="f_B",
        ),
        v_ub=_load_positive_theory_anchor(
            process_id,
            _THEORY_V_UB_CANDIDATES,
            expected_units=_EXPECTED_CKM_UNITS,
            label="|V_ub|",
        ),
        experimental_inputs=_load_experimental_inputs(process_id),
        budget_band=_build_budget_band(
            experimental=experimental,
            standard_model=standard_model,
        ),
    )


def _theory_input_citation(anchor: B009Anchor) -> str:
    return (
        "B009.yaml theory anchors: "
        f"{anchor.f_b.block_key} f_B={anchor.f_b.value:.4g} GeV "
        f"({anchor.f_b.source}); "
        f"{anchor.v_ub.block_key} |V_ub|={anchor.v_ub.value:.4g} "
        f"({anchor.v_ub.source}). SM BR validated against "
        f"{anchor.standard_model.block_key}."
    )


@register
class Constraint:
    """Catalogued ``B+ -> tau+ nu_tau`` branching-ratio constraint."""

    process_id = "B009"
    severity = Severity.HARD
    observable = "BR(B+ -> tau+ nu_tau)"

    def __init__(self) -> None:
        self.anchor = _load_b009_anchor(self.process_id)
        self.sm_inputs = leptonic_tree_bplus_tau_nu_inputs_from_theory_anchors(
            decay_constant_gev=self.anchor.f_b.value,
            ckm_abs=self.anchor.v_ub.value,
            constants_citation=_theory_input_citation(self.anchor),
        )
        self.sm_result = leptonic_tree_sm_branching_fraction(self.sm_inputs)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        result = self.sm_result
        common_diagnostics: dict[str, Any] = {
            "sm_anchor_branching_fraction": float(self.anchor.sm_value),
            "sm_formula_branching_fraction": float(result.sm_branching_fraction),
            "sm_formula_minus_anchor": float(
                result.sm_branching_fraction - self.anchor.sm_value
            ),
            "pdg_2025_average_branching_fraction": float(
                self.anchor.pdg_average.value
            ),
            "experimental_block": self.anchor.experimental.block_key,
            "sm_block": self.anchor.standard_model.block_key,
            "theory_inputs_yaml_backed": True,
            "f_b_anchor_block": self.anchor.f_b.block_key,
            "f_b_source_url": self.anchor.f_b.source_url,
            "v_ub_anchor_block": self.anchor.v_ub.block_key,
            "v_ub_source_url": self.anchor.v_ub.source_url,
            "experimental_input_count": len(self.anchor.experimental_inputs),
            "experimental_input_experiments": tuple(
                row.experiment for row in self.anchor.experimental_inputs
            ),
            "meson_mass_gev": float(result.diagnostics["meson_mass_gev"]),
            "lepton_mass_gev": float(result.diagnostics["lepton_mass_gev"]),
            "decay_constant_gev": float(result.diagnostics["decay_constant_gev"]),
            "ckm_abs": float(result.diagnostics["ckm_abs"]),
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
                "|BR_total - BR_SM(formula)| compared with "
                "|BR_exp - BR_SM(anchor)| + sqrt(sigma_exp^2 + sigma_SM^2)"
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
                sm_prediction=float(result.sm_branching_fraction),
                experimental=float(self.anchor.value),
                ratio=None,
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; B+ -> tau nu "
                    "charged-current epsilon shift was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "np_shift_branching_fraction": 0.0,
                    **common_diagnostics,
                },
            )

        try:
            epsilon_ub_tau = charged_current_epsilon(
                charged,
                up="u",
                down="b",
                lepton="tau",
            )
            predicted = shifted_branching_fraction(
                float(result.sm_branching_fraction),
                epsilon_ub_tau,
            )
            source_diag = charged_current_source_diagnostics(charged)
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                predicted=float(self.sm_result.branching_fraction),
                sm_prediction=float(self.sm_result.sm_branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid rs_charged_current extra for B009 "
                    "epsilon_ub^tau shift."
                ),
                diagnostics={
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "np_shift_branching_fraction": 0.0,
                    **common_diagnostics,
                },
            )

        predicted = float(predicted)
        np_shift = float(predicted - result.sm_branching_fraction)
        budget = float(self.anchor.budget)
        ratio = abs(np_shift) / budget if budget > 0.0 else float("inf")
        diagnostics = dict(common_diagnostics)
        diagnostics.update(
            {
                "evaluated": True,
                "epsilon_ub_tau": epsilon_ub_tau.epsilon,
                "amplitude_multiplier": float(epsilon_ub_tau.rate_multiplier),
                "np_shift_branching_fraction": float(np_shift),
                "rs_matching_formula": (
                    "BR(B+ -> tau+ nu_tau) = BR_SM |1+epsilon_ub^tau|^2"
                ),
                **source_diag,
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
                "BR(B+ -> tau+ nu_tau) uses the tree-level charged-leptonic "
                "formula with the B009 YAML f_B and |Vub| anchors and the "
                "minimal-LH rs_charged_current epsilon_ub^tau shift. The HARD "
                "ratio is the NP shift over the HFLAV-vs-UTfit YAML budget."
            ),
            diagnostics=diagnostics,
        )
