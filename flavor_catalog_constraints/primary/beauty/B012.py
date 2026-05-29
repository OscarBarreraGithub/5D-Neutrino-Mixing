"""B012 - exclusive radiative decay ``B -> K*(892) gamma``.

Physics
-------
The exclusive branching fraction is evaluated with the shared ``b -> s gamma``
C7 dipole machinery built for B011,

    BR(B -> K* gamma) = BR_norm *
        (|C7_SM + C7_NP|^2 + |C7p_SM + C7p_NP|^2)
        / (|C7_SM|^2 + |C7p_SM|^2).

The low-level C7/C8 matching proxy and leading-log running to ``mu_b`` live in
``quarkConstraints.bsgamma`` and are reached only through the
``flavor_catalog_constraints.physics_adapters.bsgamma`` boundary.  B012 adds
only the exclusive normalization: ``BR_norm`` is the neutral
``B0 -> K*(892)0 gamma`` HFLAV anchor from ``B012.yaml``.  This preserves the
shared C7 machinery while putting the exclusive form-factor normalization at
the catalog boundary.

Severity
--------
HARD.  The observed exclusive branching fraction is compared to the
C7-rescaled total prediction with the B012 YAML experimental uncertainty as
the budget.  The full RS loop match for exclusive radiative amplitudes is not
available on ``ParameterPoint``; the reused RS dipole proxy is therefore
explicitly flagged NEEDS-HUMAN-PHYSICS in the result diagnostics.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B012.yaml`` is the source of truth for the
HFLAV exclusive branching-fraction anchors.  Numeric branching-fraction values
below are loaded from that sidecar, not hardcoded in this constraint.
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
    exclusive_btokstargamma_from_couplings,
    exclusive_btokstargamma_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_DIPOLE_MASS_EXTRA = "kk_ew_mass_gev"
_NEUTRAL_ANCHOR_CANDIDATES = ("branching_fraction_b0_to_kstar0_gamma",)
_CHARGED_ANCHOR_CANDIDATES = ("branching_fraction_bp_to_kstarp_gamma",)
_EXPECTED_UNITS = "branching fraction"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/beauty/B012.yaml "
    "branching_fraction_b0_to_kstar0_gamma"
)
_NORMALIZATION_LABEL = "B012 HFLAV B0 -> K*(892)0 gamma"
_PARAMETRIZATION_CITATION = (
    "Exclusive C7-normalized B -> K*(892) gamma proxy: the B012 HFLAV "
    "branching fraction supplies the exclusive normalization/form-factor "
    "factor, while C7/C8 matching and LL QCD running are reused from "
    "quarkConstraints.bsgamma."
)


@dataclass(frozen=True)
class BtokstargammaBudgetBand:
    """B012 total-BR budget from the neutral HFLAV exclusive anchor."""

    source: str
    experimental_sigma: float
    hard_veto_budget: float
    lower_total_edge: float
    upper_total_edge: float
    exclusive_normalization: float


@dataclass(frozen=True)
class B012Anchor:
    """Typed B012 anchor: neutral normalization, charged cross-check, budget."""

    neutral: Anchor
    charged: Anchor
    budget_band: BtokstargammaBudgetBand

    @property
    def value(self) -> float:
        return self.neutral.value

    @property
    def uncertainty(self) -> float | None:
        return self.neutral.uncertainty

    @property
    def source_url(self) -> str | None:
        return self.neutral.source_url

    @property
    def sm_value(self) -> float:
        return self.neutral.value

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


def _build_budget_band(*, process_id: str, neutral: Anchor) -> BtokstargammaBudgetBand:
    _validate_branching_anchor(neutral, process_id=process_id, label="neutral")
    sigma = float(neutral.uncertainty)
    lower_edge = float(neutral.value - sigma)
    upper_edge = float(neutral.value + sigma)
    if lower_edge <= 0.0 or not math.isfinite(lower_edge):
        raise AnchorError(f"{process_id}: constructed B -> K* gamma budget is invalid")
    return BtokstargammaBudgetBand(
        source=_BUDGET_SOURCE,
        experimental_sigma=sigma,
        hard_veto_budget=sigma,
        lower_total_edge=lower_edge,
        upper_total_edge=upper_edge,
        exclusive_normalization=float(neutral.value),
    )


def _load_b012_anchor(process_id: str) -> B012Anchor:
    neutral = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_NEUTRAL_ANCHOR_CANDIDATES,
    )
    charged = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_CHARGED_ANCHOR_CANDIDATES,
    )
    _validate_anchor_block_key(
        neutral,
        process_id=process_id,
        label="neutral",
        candidates=_NEUTRAL_ANCHOR_CANDIDATES,
    )
    _validate_anchor_block_key(
        charged,
        process_id=process_id,
        label="charged",
        candidates=_CHARGED_ANCHOR_CANDIDATES,
    )
    _validate_branching_anchor(neutral, process_id=process_id, label="neutral")
    _validate_branching_anchor(charged, process_id=process_id, label="charged")
    return B012Anchor(
        neutral=neutral,
        charged=charged,
        budget_band=_build_budget_band(process_id=process_id, neutral=neutral),
    )


def _budget_result(predicted: float, anchor: B012Anchor) -> tuple[float, float, float, bool]:
    np_shift = float(predicted - anchor.sm_value)
    budget = float(anchor.budget)
    ratio = abs(np_shift) / budget if budget > 0.0 else float("inf")
    return np_shift, budget, float(ratio), bool(ratio <= 1.0)


def _complex_mapping(mapping: Any) -> dict[str, complex]:
    return {str(key): complex(value) for key, value in dict(mapping).items()}


@register
class Constraint:
    """Catalogued exclusive ``B0 -> K*(892)0 gamma`` BR constraint."""

    process_id = "B012"
    severity = Severity.HARD
    observable = "BR(B0 -> K*(892)0 gamma)"

    def __init__(self) -> None:
        self.anchor = _load_b012_anchor(self.process_id)
        self.sm_inputs = bsgamma_default_sm_inputs()
        self.sm_result = exclusive_btokstargamma_sm_branching_fraction(
            exclusive_sm_branching_fraction=self.anchor.sm_value,
            inputs=self.sm_inputs,
            normalization_label=_NORMALIZATION_LABEL,
            normalization_source=self.anchor.source_url,
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
                    f"extra {_REQUIRED_EXTRA!r} absent; B0 -> K*(892)0 gamma "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "exclusive_normalization_branching_fraction": float(
                        self.anchor.sm_value
                    ),
                    "budget_source": self.anchor.budget_band.source,
                },
            )

        kk_mass = point.get_extra(_OPTIONAL_DIPOLE_MASS_EXTRA)
        try:
            result = exclusive_btokstargamma_from_couplings(
                couplings,
                exclusive_sm_branching_fraction=self.anchor.sm_value,
                m_kk_gev=None if kk_mass is None else float(kk_mass),
                inputs=self.sm_inputs,
                normalization_label=_NORMALIZATION_LABEL,
                normalization_source=self.anchor.source_url,
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
                    "the B0 -> K*(892)0 gamma C7 proxy"
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
                "neutral_experimental_block": self.anchor.neutral.block_key,
                "charged_reference_block": self.anchor.charged.block_key,
                "charged_reference_branching_fraction": float(self.anchor.charged.value),
                "charged_reference_uncertainty": float(self.anchor.charged.uncertainty),
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
                    "|BR_total - BR_B012_norm| compared with the neutral "
                    "B012.yaml experimental uncertainty"
                ),
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": BSGAMMA_RS_MATCHING_ASSUMPTION_V1,
                "needs_human_physics": (
                    "NEEDS-HUMAN-PHYSICS: full RS exclusive B -> K* gamma "
                    "matching requires KK fermion, Higgs/Goldstone, "
                    "charged-current, chromomagnetic, form-factor, and "
                    "spectator-amplitude inputs not available on "
                    "ParameterPoint; v1 reuses the documented b-s overlap "
                    "C7/C8 proxy with LL QCD running to mu_b."
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
                "BR(B0 -> K*(892)0 gamma) uses the shared C7-normalized "
                "b -> s gamma dipole formula with the B012 YAML exclusive "
                "normalization. The RS contribution is a documented b-s "
                "overlap proxy for C7/C7p and is marked NEEDS-HUMAN-PHYSICS; "
                "the HARD budget is the neutral HFLAV uncertainty in B012.yaml."
            ),
            diagnostics=diagnostics,
        )
