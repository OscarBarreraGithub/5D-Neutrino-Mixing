"""K002 - neutral-kaon mass splitting.

Physics
-------
``Delta m_K`` is the ``K_L-K_S`` mass difference.  The new-physics
contribution is evaluated with the audited Delta F = 2 core through
``flavor_catalog_constraints.physics_adapters.deltaf2``:

    Delta m_K^NP = 2 |M12^NP|.

The HARD veto is expressed in the equivalent core convention
``|M12^NP| <= Delta m_K^exp / 2``.  The SM short-distance prediction is not
subtracted because the observable is long-distance dominated.

Catalog sidecar
---------------
``flavor_catalog/processes/kaon/K002.yaml`` stores the PDG CPT-assuming fit in
``10^10 hbar s^-1``.  This module loads that value through the scaffold anchor
loader and converts it to GeV before building the ``Delta m_K^exp / 2`` budget.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

from scipy import constants

from flavor_catalog_constraints.anchors import Anchor, AnchorError, load_anchor
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.deltaf2 import (
    delta_mk_core_inputs,
    delta_mk_from_wilsons_with_running,
    delta_mk_wilsons_from_couplings,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "kaon"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_EXPERIMENTAL_ANCHOR_CANDIDATES = ("pdg_fit_assuming_cpt",)
_EXPECTED_UNITS = "10^10 hbar s^-1"
_UNIT_SCALE_HBAR_PER_SECOND = 1.0e10
_HBAR_GEV_SECONDS = float(constants.hbar / constants.electron_volt * 1.0e-9)
_MU_HAD_GEV = 2.0
_BUDGET_DOC_CITATION = (
    "docs/quark_scan_assumptions_compact.tex:464-466; "
    "quarkConstraints/deltaf2.py:809-825"
)
_INPUT_AUDIT_CITATION = "docs/audits/bag_param_inventory.md:20"


@dataclass(frozen=True)
class DeltaMKAnchor:
    """Typed K002 anchor: PDG fit converted to GeV plus NP budget."""

    experimental: Anchor
    value: float
    uncertainty: float | None
    budget: float
    hbar_gev_seconds: float
    unit_scale_hbar_per_second: float
    budget_doc_citation: str
    input_audit_citation: str

    @property
    def source_url(self) -> str | None:
        """Experimental source URL, kept for common Anchor-like access."""
        return self.experimental.source_url


def _convert_pdg_hbar_rate_to_gev(
    value: float,
    *,
    process_id: str,
    field_name: str,
) -> float:
    converted = float(value) * _UNIT_SCALE_HBAR_PER_SECOND * _HBAR_GEV_SECONDS
    if not math.isfinite(converted) or converted <= 0.0:
        raise AnchorError(
            f"{process_id}: converted {field_name} must be positive and finite, "
            f"got {converted!r}"
        )
    return converted


def _load_delta_mk_anchor(process_id: str) -> DeltaMKAnchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_EXPERIMENTAL_ANCHOR_CANDIDATES,
    )
    if experimental.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: expected units {_EXPECTED_UNITS!r} for Delta m_K, "
            f"got {experimental.units!r}"
        )

    value_gev = _convert_pdg_hbar_rate_to_gev(
        experimental.value,
        process_id=process_id,
        field_name="value",
    )
    uncertainty_gev = (
        None
        if experimental.uncertainty is None
        else _convert_pdg_hbar_rate_to_gev(
            experimental.uncertainty,
            process_id=process_id,
            field_name="uncertainty",
        )
    )
    budget = value_gev / 2.0
    if not math.isfinite(budget) or budget <= 0.0:
        raise AnchorError(f"{process_id}: Delta m_K NP budget must be positive")

    return DeltaMKAnchor(
        experimental=experimental,
        value=value_gev,
        uncertainty=uncertainty_gev,
        budget=budget,
        hbar_gev_seconds=_HBAR_GEV_SECONDS,
        unit_scale_hbar_per_second=_UNIT_SCALE_HBAR_PER_SECOND,
        budget_doc_citation=_BUDGET_DOC_CITATION,
        input_audit_citation=_INPUT_AUDIT_CITATION,
    )


def _complex_mapping(values: Mapping[str, complex]) -> dict[str, complex]:
    return {str(key): complex(value) for key, value in values.items()}


def _budget_pull_sigma(
    *,
    catalog_value: float,
    catalog_uncertainty: float | None,
    core_value: float,
) -> float | None:
    if catalog_uncertainty is None or catalog_uncertainty <= 0.0:
        return None
    return float((core_value - catalog_value) / catalog_uncertainty)


@register
class Constraint:
    """Catalogued Delta m_K Delta F=2 constraint (process_id K002)."""

    process_id = "K002"
    severity = Severity.HARD
    observable = "Delta m_K"

    def __init__(self) -> None:
        self.anchor = _load_delta_mk_anchor(self.process_id)

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; Delta m_K constraint "
                    "was not evaluated."
                ),
                diagnostics={"missing_extra": _REQUIRED_EXTRA},
            )

        wilsons = delta_mk_wilsons_from_couplings(couplings)
        result = delta_mk_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
            m12_np_budget=self.anchor.budget,
        )
        core_inputs = delta_mk_core_inputs()
        core_budget = float(core_inputs["core_m12_budget_gev"])
        predicted = float(result.abs_m12_np)
        ratio = float(result.ratio_to_exp)
        budget = float(self.anchor.budget)
        core_ratio = predicted / core_budget if core_budget > 0.0 else float("inf")

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(result.passes),
            predicted=predicted,
            experimental=float(self.anchor.value),
            ratio=ratio,
            budget=budget,
            notes=(
                "|M12^NP| is compared to Delta m_K^exp / 2 with no SM "
                "short-distance subtraction; Wilsons are QCD-evolved to 2 GeV."
            ),
            diagnostics={
                "abs_m12_np_gev": predicted,
                "qcd_running_applied": True,
                "hadronic_scale_gev": _MU_HAD_GEV,
                "matching_scale_gev": float(wilsons.matching_scale),
                "m_kk_gev": float(wilsons.M_KK),
                "left_sd_coupling": complex(wilsons.left_coupling),
                "right_sd_coupling": complex(wilsons.right_coupling),
                "wilson_coefficients": _complex_mapping(wilsons.wilsons),
                "sm_subtracted": False,
                "long_distance_dominated": True,
                "budget_construction": "Delta m_K^exp / 2",
                "budget_doc_citation": self.anchor.budget_doc_citation,
                "input_audit_citation": self.anchor.input_audit_citation,
                "experimental_block": self.anchor.experimental.block_key,
                "experimental_units": self.anchor.experimental.units,
                "experimental_value_raw": float(self.anchor.experimental.value),
                "experimental_uncertainty_raw": (
                    None
                    if self.anchor.experimental.uncertainty is None
                    else float(self.anchor.experimental.uncertainty)
                ),
                "delta_m_k_exp_gev": float(self.anchor.value),
                "delta_m_k_uncertainty_gev": (
                    None if self.anchor.uncertainty is None else float(self.anchor.uncertainty)
                ),
                "hbar_gev_seconds": self.anchor.hbar_gev_seconds,
                "unit_scale_hbar_per_second": self.anchor.unit_scale_hbar_per_second,
                "core_delta_m_k_exp_gev": float(core_inputs["delta_m_k_exp_gev"]),
                "core_m12_budget_gev": core_budget,
                "core_ratio_to_exp_without_catalog_override": float(core_ratio),
                "core_catalog_budget_pull_sigma": _budget_pull_sigma(
                    catalog_value=self.anchor.value,
                    catalog_uncertainty=self.anchor.uncertainty,
                    core_value=float(core_inputs["delta_m_k_exp_gev"]),
                ),
                "catalog_core_budget_relative_difference": float(
                    (budget - core_budget) / core_budget
                ),
            },
        )
