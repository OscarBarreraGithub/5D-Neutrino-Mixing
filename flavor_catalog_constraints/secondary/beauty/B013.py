"""B013 - exclusive radiative decay ``B_s0 -> phi(1020) gamma``.

Physics
-------
The exclusive branching fraction is evaluated with the shared ``b -> s gamma``
C7/C7' dipole machinery built for B011 and reused by B012,

    BR(B_s -> phi gamma) = BR_norm *
        (|C7_SM + C7_NP|^2 + |C7p_SM + C7p_NP|^2)
        / (|C7_SM|^2 + |C7p_SM|^2).

The low-level C7/C8 matching proxy and leading-log running to ``mu_b`` live in
``quarkConstraints.bsgamma`` and are reached only through the
``flavor_catalog_constraints.physics_adapters.bsgamma`` boundary.  B013 adds
only the exclusive ``B_s -> phi`` normalization: ``BR_norm`` is the catalog
PDG measurement normalization proxy loaded from ``B013.yaml``.

Severity
--------
HARD.  ``B013.yaml`` has no theory-only exclusive SM branching-fraction
prediction.  The observed exclusive branching fraction is therefore compared to
the C7-rescaled total prediction with a conservative measurement-consistency
NP-shift band,

    |BR_total - BR_ref| <=
        |BR_exp - BR_ref| + sqrt(sigma_exp^2 + sigma_ref^2).

This is deliberately not presented as SM-theory-vs-data room: ``BR_ref`` and
``sigma_ref`` are the PDG measurement row used as the normalization proxy.  The
full RS loop match for exclusive radiative amplitudes is not available on
``ParameterPoint``; the reused RS dipole proxy is explicitly flagged
NEEDS-HUMAN-PHYSICS in the result diagnostics.  The photon-helicity/
time-dependent CP observables ``A_Delta`` and ``S_phi gamma`` are loaded for
provenance but are not evaluated because the present adapter does not compute
the full right-handed-photon likelihood.

Catalog sidecar
---------------
``flavor_catalog/processes/secondary/beauty/B013.yaml`` is the source of
truth for the B_s -> phi gamma branching-fraction rows and the LHCb helicity
observables.  Numeric values below are loaded from that sidecar, not hardcoded
in this constraint.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
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
    exclusive_bsphigamma_from_couplings,
    exclusive_bsphigamma_sm_branching_fraction,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_TIER = ConstraintLevel.SECONDARY
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_OPTIONAL_DIPOLE_MASS_EXTRA = "kk_ew_mass_gev"
_EXPERIMENTAL_ANCHOR_CANDIDATES = ("branching_fraction_hflav_2024",)
_REFERENCE_MEASUREMENT_ANCHOR_CANDIDATES = ("branching_fraction_pdg_2025",)
_S_PHIGAMMA_CANDIDATES = ("s_phigamma_lhcb_2019",)
_A_DELTA_CANDIDATES = ("a_delta_phigamma_lhcb_2019",)
_EXPECTED_BRANCHING_UNITS = "branching fraction"
_EXPECTED_DIMENSIONLESS_UNITS = "dimensionless"
_BUDGET_SOURCE = (
    "flavor_catalog/processes/secondary/beauty/B013.yaml "
    "measurement-consistency band from branching_fraction_hflav_2024 + "
    "branching_fraction_pdg_2025"
)
_NORMALIZATION_LABEL = (
    "B013 B_s0 -> phi(1020) gamma PDG measurement normalization proxy"
)
_PARAMETRIZATION_CITATION = (
    "Exclusive C7-normalized B_s -> phi(1020) gamma proxy: B013.yaml "
    "supplies a PDG measurement normalization proxy, while C7/C8 matching "
    "and LL QCD running are reused from quarkConstraints.bsgamma."
)
_BUDGET_POLICY = (
    "measurement-consistency band: |BR_total - BR_ref_measurement| compared "
    "with |BR_exp - BR_ref_measurement| + "
    "sqrt(sigma_exp^2 + sigma_ref_measurement^2); this is not a "
    "SM-theory-vs-data budget because B013.yaml has no theory-only "
    "exclusive SM prediction."
)
_SM_THEORY_GAP_NOTE = (
    "NEEDS-HUMAN-PHYSICS: B013.yaml contains no theory-only SCET/QCDF/"
    "form-factor SM prediction for BR(B_s0 -> phi(1020) gamma). A rigorous "
    "exclusive theory prediction and uncertainty are required to set a true "
    "SM-theory-vs-data NP budget."
)
_NEEDS_HUMAN_PHYSICS = (
    "NEEDS-HUMAN-PHYSICS: full RS exclusive B_s -> phi gamma matching "
    "requires KK fermion, Higgs/Goldstone, charged-current, chromomagnetic, "
    "B_s -> phi form-factor, spectator-amplitude, and time-dependent "
    "right-handed-photon inputs not available on ParameterPoint; v1 reuses "
    "the documented b-s overlap C7/C8 proxy with LL QCD running to mu_b. "
    + _SM_THEORY_GAP_NOTE
    + " The A_Delta and S_phi gamma photon-helicity observables are loaded but "
    "not evaluated by this branching-fraction proxy."
)


@dataclass(frozen=True)
class BsphigammaBudgetBand:
    """Measurement-consistency B013 NP-shift budget for the HARD veto."""

    source: str
    central_residual: float
    experimental_sigma: float
    reference_measurement_sigma: float
    combined_sigma: float
    hard_veto_budget: float
    lower_total_edge: float
    upper_total_edge: float


@dataclass(frozen=True)
class B013Anchor:
    """Typed B013 anchors: branching fraction, helicity rows, and budget."""

    experimental: Anchor
    reference_measurement: Anchor
    s_phigamma: Anchor
    a_delta_phigamma: Anchor
    raw_s_phigamma: Mapping[str, Any]
    raw_a_delta_phigamma: Mapping[str, Any]
    budget_band: BsphigammaBudgetBand

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
    def reference_value(self) -> float:
        return self.reference_measurement.value

    @property
    def sm_value(self) -> float:
        """Compatibility alias for the C7 no-NP normalization proxy."""

        return self.reference_value

    @property
    def budget(self) -> float:
        return self.budget_band.hard_veto_budget


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
    if anchor.uncertainty is None or anchor.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: {label} anchor {anchor.block_key!r} must provide "
            "a positive uncertainty"
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


def _positive_optional_float(row: Mapping[str, Any], key: str) -> float | None:
    value = row.get(key)
    if value is None:
        return None
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        return None
    return number


def _combined_uncertainty(row: Mapping[str, Any], *keys: str) -> float | None:
    pieces = [_positive_optional_float(row, key) for key in keys]
    if any(piece is None for piece in pieces):
        return None
    return float(math.sqrt(sum(float(piece) ** 2 for piece in pieces)))


def _build_budget_band(
    *,
    process_id: str,
    experimental: Anchor,
    reference_measurement: Anchor,
) -> BsphigammaBudgetBand:
    _validate_branching_anchor(
        experimental,
        process_id=process_id,
        label="experimental",
    )
    _validate_branching_anchor(
        reference_measurement,
        process_id=process_id,
        label="reference measurement",
    )
    exp_sigma = float(experimental.uncertainty)
    ref_sigma = float(reference_measurement.uncertainty)
    combined_sigma = math.sqrt(exp_sigma * exp_sigma + ref_sigma * ref_sigma)
    central_residual = abs(experimental.value - reference_measurement.value)
    hard_budget = central_residual + combined_sigma
    if hard_budget <= 0.0 or not math.isfinite(hard_budget):
        raise AnchorError(f"{process_id}: constructed B_s -> phi gamma budget is invalid")
    return BsphigammaBudgetBand(
        source=_BUDGET_SOURCE,
        central_residual=float(central_residual),
        experimental_sigma=exp_sigma,
        reference_measurement_sigma=ref_sigma,
        combined_sigma=float(combined_sigma),
        hard_veto_budget=float(hard_budget),
        lower_total_edge=float(reference_measurement.value - hard_budget),
        upper_total_edge=float(reference_measurement.value + hard_budget),
    )


def _load_b013_anchor(process_id: str) -> B013Anchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        tier=_TIER,
        candidates=_EXPERIMENTAL_ANCHOR_CANDIDATES,
    )
    reference_measurement = load_anchor(
        process_id,
        family=_FAMILY,
        tier=_TIER,
        candidates=_REFERENCE_MEASUREMENT_ANCHOR_CANDIDATES,
    )
    s_phigamma = load_anchor(
        process_id,
        family=_FAMILY,
        tier=_TIER,
        candidates=_S_PHIGAMMA_CANDIDATES,
        uncertainty_key="uncertainty_stat",
    )
    a_delta_phigamma = load_anchor(
        process_id,
        family=_FAMILY,
        tier=_TIER,
        candidates=_A_DELTA_CANDIDATES,
        uncertainty_key="uncertainty_stat_plus",
    )
    _validate_anchor_block_key(
        experimental,
        process_id=process_id,
        label="experimental",
        candidates=_EXPERIMENTAL_ANCHOR_CANDIDATES,
    )
    _validate_anchor_block_key(
        reference_measurement,
        process_id=process_id,
        label="reference measurement",
        candidates=_REFERENCE_MEASUREMENT_ANCHOR_CANDIDATES,
    )
    _validate_anchor_block_key(
        s_phigamma,
        process_id=process_id,
        label="S_phi gamma",
        candidates=_S_PHIGAMMA_CANDIDATES,
    )
    _validate_anchor_block_key(
        a_delta_phigamma,
        process_id=process_id,
        label="A_Delta phi gamma",
        candidates=_A_DELTA_CANDIDATES,
    )
    _validate_branching_anchor(
        experimental,
        process_id=process_id,
        label="experimental",
    )
    _validate_branching_anchor(
        reference_measurement,
        process_id=process_id,
        label="reference measurement",
    )
    _validate_dimensionless_anchor(
        s_phigamma,
        process_id=process_id,
        label="S_phi gamma",
    )
    _validate_dimensionless_anchor(
        a_delta_phigamma,
        process_id=process_id,
        label="A_Delta phi gamma",
    )
    pdg_block = load_pdg_block(process_id, family=_FAMILY, tier=_TIER)
    return B013Anchor(
        experimental=experimental,
        reference_measurement=reference_measurement,
        s_phigamma=s_phigamma,
        a_delta_phigamma=a_delta_phigamma,
        raw_s_phigamma=dict(pdg_block[s_phigamma.block_key]),
        raw_a_delta_phigamma=dict(pdg_block[a_delta_phigamma.block_key]),
        budget_band=_build_budget_band(
            process_id=process_id,
            experimental=experimental,
            reference_measurement=reference_measurement,
        ),
    )


def _budget_result(predicted: float, anchor: B013Anchor) -> tuple[float, float, float, bool]:
    np_shift = float(predicted - anchor.reference_value)
    budget = float(anchor.budget)
    ratio = abs(np_shift) / budget if budget > 0.0 else float("inf")
    return np_shift, budget, float(ratio), bool(ratio <= 1.0)


def _complex_mapping(mapping: Any) -> dict[str, complex]:
    return {str(key): complex(value) for key, value in dict(mapping).items()}


def _budget_diagnostics(anchor: B013Anchor) -> dict[str, Any]:
    return {
        "experimental_sigma": float(anchor.budget_band.experimental_sigma),
        "budget_reference_measurement_sigma": float(
            anchor.budget_band.reference_measurement_sigma
        ),
        "budget_sm_theory_sigma": None,
        "budget_combined_sigma": float(anchor.budget_band.combined_sigma),
        "budget_central_residual": float(anchor.budget_band.central_residual),
        "hard_veto_np_shift_budget": float(anchor.budget_band.hard_veto_budget),
        "budget_lower_total_edge": float(anchor.budget_band.lower_total_edge),
        "budget_upper_total_edge": float(anchor.budget_band.upper_total_edge),
        "budget_source": anchor.budget_band.source,
        "budget_policy": _BUDGET_POLICY,
        "sm_theory_prediction_available": False,
        "sm_theory_prediction_gap": _SM_THEORY_GAP_NOTE,
    }


def _helicity_diagnostics(anchor: B013Anchor) -> dict[str, Any]:
    s_combined = _combined_uncertainty(
        anchor.raw_s_phigamma,
        "uncertainty_stat",
        "uncertainty_syst",
    )
    a_delta_plus = _combined_uncertainty(
        anchor.raw_a_delta_phigamma,
        "uncertainty_stat_plus",
        "uncertainty_syst",
    )
    a_delta_minus = _combined_uncertainty(
        anchor.raw_a_delta_phigamma,
        "uncertainty_stat_minus",
        "uncertainty_syst",
    )
    return {
        "helicity_observables_evaluated": False,
        "s_phigamma_block": anchor.s_phigamma.block_key,
        "s_phigamma_experimental": float(anchor.s_phigamma.value),
        "s_phigamma_uncertainty_stat": float(anchor.s_phigamma.uncertainty),
        "s_phigamma_uncertainty_syst": float(
            anchor.raw_s_phigamma["uncertainty_syst"]
        ),
        "s_phigamma_uncertainty_combined": s_combined,
        "a_delta_phigamma_block": anchor.a_delta_phigamma.block_key,
        "a_delta_phigamma_experimental": float(anchor.a_delta_phigamma.value),
        "a_delta_phigamma_uncertainty_stat_plus": float(
            anchor.raw_a_delta_phigamma["uncertainty_stat_plus"]
        ),
        "a_delta_phigamma_uncertainty_stat_minus": float(
            anchor.raw_a_delta_phigamma["uncertainty_stat_minus"]
        ),
        "a_delta_phigamma_uncertainty_syst": float(
            anchor.raw_a_delta_phigamma["uncertainty_syst"]
        ),
        "a_delta_phigamma_uncertainty_combined_plus": a_delta_plus,
        "a_delta_phigamma_uncertainty_combined_minus": a_delta_minus,
        "helicity_observable_note": (
            "A_Delta and S_phi gamma constrain right-handed photon "
            "polarization/time-dependent CP, but this B013 implementation "
            "only evaluates the branching fraction."
        ),
    }


@register
class Constraint:
    """Catalogued exclusive ``B_s0 -> phi(1020) gamma`` BR constraint."""

    process_id = "B013"
    severity = Severity.HARD
    observable = "BR(B_s0 -> phi(1020) gamma)"

    def __init__(self) -> None:
        self.anchor = _load_b013_anchor(self.process_id)
        self.sm_inputs = bsgamma_default_sm_inputs()
        self.normalization_result = exclusive_bsphigamma_sm_branching_fraction(
            exclusive_sm_branching_fraction=self.anchor.reference_value,
            inputs=self.sm_inputs,
            normalization_label=_NORMALIZATION_LABEL,
            normalization_source=self.anchor.reference_measurement.source_url,
        )

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.normalization_result.branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; B_s0 -> phi(1020) gamma "
                    "constraint was not evaluated."
                ),
                diagnostics={
                    **_budget_diagnostics(self.anchor),
                    **_helicity_diagnostics(self.anchor),
                    "evaluated": False,
                    "missing_extra": _REQUIRED_EXTRA,
                    "reference_measurement_branching_fraction": float(
                        self.anchor.reference_value
                    ),
                    "exclusive_normalization_branching_fraction": float(
                        self.anchor.reference_value
                    ),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )

        kk_mass = point.get_extra(_OPTIONAL_DIPOLE_MASS_EXTRA)
        try:
            result = exclusive_bsphigamma_from_couplings(
                couplings,
                exclusive_sm_branching_fraction=self.anchor.reference_value,
                m_kk_gev=None if kk_mass is None else float(kk_mass),
                inputs=self.sm_inputs,
                normalization_label=_NORMALIZATION_LABEL,
                normalization_source=self.anchor.reference_measurement.source_url,
            )
        except Exception as exc:  # noqa: BLE001 - constraints degrade cleanly
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.normalization_result.branching_fraction),
                experimental=float(self.anchor.value),
                budget=float(self.anchor.budget),
                notes=(
                    "NOT EVALUATED - invalid quark_mass_basis_couplings for "
                    "the B_s0 -> phi(1020) gamma C7 proxy"
                ),
                diagnostics={
                    **_budget_diagnostics(self.anchor),
                    **_helicity_diagnostics(self.anchor),
                    "evaluated": False,
                    "invalid_extra": _REQUIRED_EXTRA,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                },
            )

        predicted = float(result.branching_fraction)
        np_shift, budget, ratio, passes = _budget_result(predicted, self.anchor)
        diagnostics = dict(result.diagnostics)
        wilsons = result.wilsons.wilsons if result.wilsons is not None else {}
        diagnostics.update(
            {
                "evaluated": True,
                "reference_measurement_branching_fraction": float(
                    self.anchor.reference_value
                ),
                "exclusive_normalization_branching_fraction": float(
                    self.anchor.reference_value
                ),
                "normalization_formula_branching_fraction": float(
                    result.sm_branching_fraction
                ),
                "normalization_formula_minus_reference_measurement": float(
                    result.sm_branching_fraction - self.anchor.reference_value
                ),
                "experimental_block": self.anchor.experimental.block_key,
                "reference_measurement_block": (
                    self.anchor.reference_measurement.block_key
                ),
                "normalization_anchor_policy": (
                    "B013.yaml has no separate theory-only branching-fraction "
                    "block; the PDG 2025 3.4e-5 row is used as the catalog "
                    "C7 no-NP/exclusive normalization proxy, not as a "
                    "theory prediction."
                ),
                "np_shift_branching_fraction": float(np_shift),
                "ratio_to_sm_c7_power": float(result.ratio_to_sm),
                "c7_sm_eff": complex(result.c7_sm_eff),
                "c7p_sm_eff": complex(result.c7p_sm_eff),
                "c7_total": complex(result.c7_total),
                "c7p_total": complex(result.c7p_total),
                "c7_np": complex(result.c7_np),
                "c7p_np": complex(result.c7p_np),
                "wilson_coefficients": _complex_mapping(wilsons),
                "parametrization_citation": _PARAMETRIZATION_CITATION,
                "rs_matching_assumption": BSGAMMA_RS_MATCHING_ASSUMPTION_V1,
                "needs_human_physics": _NEEDS_HUMAN_PHYSICS,
                "kk_ew_mass_extra_used": kk_mass is not None,
                "down_sector_indices": (1, 2),
            }
        )
        diagnostics.update(_budget_diagnostics(self.anchor))
        diagnostics.update(_helicity_diagnostics(self.anchor))

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
                "BR(B_s0 -> phi(1020) gamma) uses the shared C7-normalized "
                "b -> s gamma dipole formula with the B013 YAML exclusive "
                "PDG measurement normalization proxy. The RS contribution is "
                "a documented b-s overlap proxy for C7/C7p and is marked "
                "NEEDS-HUMAN-PHYSICS; the HARD budget is a conservative "
                "measurement-consistency band, not SM-theory-vs-data room. "
                "A_Delta and S_phi gamma are not evaluated."
            ),
            diagnostics=diagnostics,
        )
