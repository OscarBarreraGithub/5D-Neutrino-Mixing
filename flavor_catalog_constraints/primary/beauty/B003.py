"""B003 - B_s mixing mass splitting.

Physics
-------
``Delta m_s`` is the oscillation frequency / mass splitting in
``B_s^0-B_sbar^0`` mixing.  The new-physics contribution is evaluated through
the audited Delta F = 2 core,

    Delta m_s^NP = 2 |M12^NP|,

with ``M12^NP`` built from the B_s Wilson coefficients produced by
``quarkConstraints.deltaf2``.  This module reaches that core only through the
``flavor_catalog_constraints.physics_adapters.deltaf2`` boundary, and uses the
QCD-running path that evolves Wilsons to ``mu_had = m_b`` before applying the
core B_s matrix elements.

Severity
--------
HARD.  ``Delta m_s`` is observed and SM-predicted.  The HARD veto uses an
uncertainty-aware NP room in the M12 amplitude,

    (|Delta m_exp - Delta m_SM| + sigma_combined) / 2,

where ``sigma_combined`` combines the HFLAV experimental uncertainty and the
FLAG ``f_Bs sqrt(Bhat_Bs)`` normalization uncertainty in quadrature.  The core
and catalog now use this same promoted B003 budget; the old full measured
splitting envelope is retained only as a diagnostic comparison.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B003.yaml`` is the source of truth for the
HFLAV experimental value and FLAG B_s mixing normalization provenance.  The
SM central value is the Delta F = 2 core convention documented in
``docs/audits/bag_param_inventory.md:38-45``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from flavor_catalog_constraints.anchors import (
    Anchor,
    AnchorError,
    load_anchor,
    load_full_yaml,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.deltaf2 import (
    bs_mixing_core_inputs,
    bs_mixing_from_wilsons_with_running,
    bs_mixing_wilsons_from_couplings,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"
_EXPECTED_OBSERVABLE = "Delta m_s"
_EXPECTED_UNITS = "ps^-1"

_EXPERIMENTAL_ANCHOR_CANDIDATES = ("canonical_hflav_recommended",)
_AUXILIARY_THEORY_KEY = "auxiliary_theory_inputs"
_FLAG_BMIXING_CANDIDATES = ("flag_2024_bmixing",)
_BUDGET_DOC_CITATION = (
    "flavor_catalog/processes/beauty/B003.yaml:"
    "pdg_or_equivalent.canonical_hflav_recommended,"
    "auxiliary_theory_inputs.flag_2024_bmixing; "
    "docs/audits/bag_param_inventory.md:38-45"
)
_MU_HAD_GEV = 4.18
_HBAR_GEV_PER_PS = 6.582119569e-13


@dataclass(frozen=True)
class FlagBsMixingInputs:
    """Typed view over the FLAG B_s-mixing auxiliary theory block."""

    block_key: str
    source: str | None
    year: int | None
    n_f: str | None
    f_bs_sqrt_bhat_bs_value_mev: float
    f_bs_sqrt_bhat_bs_uncertainty_mev: float
    f_bs_sqrt_bhat_bs_units: str | None
    xi_value: float | None
    xi_uncertainty: float | None
    source_url: str | None
    snapshot_path: str | None

    @property
    def f_bs_sqrt_bhat_relative_sigma(self) -> float:
        return self.f_bs_sqrt_bhat_bs_uncertainty_mev / self.f_bs_sqrt_bhat_bs_value_mev


@dataclass(frozen=True)
class BsMixingBudgetBand:
    """Uncertainty-aware B_s M12 NP budget for the B003 HARD veto."""

    doc_citation: str
    central_budget: float
    loose_budget: float
    hard_veto_budget: float
    experimental_delta_m_gev: float
    experimental_delta_m_ps_inv: float
    experimental_m12_gev: float
    sm_delta_m_gev: float
    sm_delta_m_ps_inv: float
    sm_m12_gev: float
    experimental_sigma_delta_m_gev: float
    sm_theory_sigma_delta_m_gev: float
    combined_sigma_delta_m_gev: float
    flag_f_bs_sqrt_bhat_relative_sigma: float
    core_delta_m_exp_gev: float
    core_delta_m_sm_gev: float
    core_default_m12_budget: float
    core_budget_policy_id: str
    core_budget_confidence_level: str
    core_legacy_m12_budget: float


@dataclass(frozen=True)
class BsMixingAnchor:
    """Typed B003 anchor: HFLAV experiment, FLAG inputs, NP budget band."""

    experimental: Anchor
    flag_bmixing: FlagBsMixingInputs
    budget_band: BsMixingBudgetBand

    @property
    def value(self) -> float:
        """Experimental Delta m_s central value in ps^-1."""
        return self.experimental.value

    @property
    def uncertainty(self) -> float | None:
        """Experimental Delta m_s uncertainty in ps^-1."""
        return self.experimental.uncertainty

    @property
    def source_url(self) -> str | None:
        return self.experimental.source_url

    @property
    def sm_value(self) -> float:
        """Core SM M12 amplitude convention in GeV."""
        return self.budget_band.sm_m12_gev

    @property
    def central_budget(self) -> float:
        """Central residual ``|Delta m_exp - Delta m_SM| / 2`` in GeV."""
        return self.budget_band.central_budget

    @property
    def budget(self) -> float:
        """Uncertainty-aware HARD veto budget for ``|M12^NP|`` in GeV."""
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


def _required_mapping(
    value: Any,
    *,
    process_id: str,
    field_name: str,
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise AnchorError(
            f"{process_id}: B_s mixing field {field_name!r} is not a mapping "
            f"(got {type(value).__name__})"
        )
    return value


def _find_auxiliary_theory_block(
    auxiliary: Mapping[str, Any],
    candidates: tuple[str, ...],
    *,
    process_id: str,
) -> tuple[str, Mapping[str, Any]]:
    for key in candidates:
        if key in auxiliary:
            sub = auxiliary[key]
            if not isinstance(sub, Mapping):
                raise AnchorError(
                    f"{process_id}: auxiliary_theory_inputs sub-block {key!r} "
                    f"is not a mapping (got {type(sub).__name__})"
                )
            return key, sub
    present = sorted(str(key) for key in auxiliary)
    raise AnchorError(
        f"{process_id}: none of the expected auxiliary theory keys "
        f"{list(candidates)} found in auxiliary_theory_inputs "
        f"(present keys: {present})"
    )


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: B_s mixing field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, field_name=field_name)


def _validate_delta_ms_anchor(experimental: Anchor, *, process_id: str) -> None:
    if experimental.units != _EXPECTED_UNITS:
        raise AnchorError(
            f"{process_id}: Delta m_s anchor {experimental.block_key!r} must "
            f"use units {_EXPECTED_UNITS!r}, got {experimental.units!r}"
        )
    if experimental.observable != _EXPECTED_OBSERVABLE:
        raise AnchorError(
            f"{process_id}: Delta m_s anchor {experimental.block_key!r} must "
            f"have observable {_EXPECTED_OBSERVABLE!r}, got "
            f"{experimental.observable!r}"
        )


def _load_flag_bmixing_inputs(process_id: str) -> FlagBsMixingInputs:
    data = load_full_yaml(process_id, family=_FAMILY)
    auxiliary = _required_mapping(
        data.get(_AUXILIARY_THEORY_KEY),
        process_id=process_id,
        field_name=_AUXILIARY_THEORY_KEY,
    )
    block_key, sub = _find_auxiliary_theory_block(
        auxiliary,
        _FLAG_BMIXING_CANDIDATES,
        process_id=process_id,
    )
    f_block = _required_mapping(
        sub.get("f_Bs_sqrt_Bhat_Bs"),
        process_id=process_id,
        field_name=f"{block_key}.f_Bs_sqrt_Bhat_Bs",
    )
    xi_block = _required_mapping(
        sub.get("xi"),
        process_id=process_id,
        field_name=f"{block_key}.xi",
    )

    f_value = _required_float(
        f_block.get("value"),
        process_id=process_id,
        field_name=f"{block_key}.f_Bs_sqrt_Bhat_Bs.value",
    )
    f_uncertainty = _required_float(
        f_block.get("uncertainty"),
        process_id=process_id,
        field_name=f"{block_key}.f_Bs_sqrt_Bhat_Bs.uncertainty",
    )
    if f_value <= 0.0 or f_uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: FLAG f_Bs sqrt(Bhat_Bs) value and uncertainty "
            f"must be positive, got value={f_value}, uncertainty={f_uncertainty}"
        )

    return FlagBsMixingInputs(
        block_key=block_key,
        source=_optional_str(sub.get("source")),
        year=_optional_int(sub.get("year")),
        n_f=_optional_str(sub.get("n_f")),
        f_bs_sqrt_bhat_bs_value_mev=f_value,
        f_bs_sqrt_bhat_bs_uncertainty_mev=f_uncertainty,
        f_bs_sqrt_bhat_bs_units=_optional_str(f_block.get("units")),
        xi_value=_optional_float(
            xi_block.get("value"),
            process_id=process_id,
            field_name=f"{block_key}.xi.value",
        ),
        xi_uncertainty=_optional_float(
            xi_block.get("uncertainty"),
            process_id=process_id,
            field_name=f"{block_key}.xi.uncertainty",
        ),
        source_url=_optional_str(sub.get("source_url")),
        snapshot_path=_optional_str(sub.get("snapshot_path")),
    )


def _build_budget_band(
    *,
    process_id: str,
    experimental: Anchor,
    flag_bmixing: FlagBsMixingInputs,
) -> BsMixingBudgetBand:
    if experimental.uncertainty is None or experimental.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: Delta m_s experimental uncertainty is required "
            "for the B_s SM-vs-exp budget band"
        )

    core_inputs = bs_mixing_core_inputs()
    experimental_delta_m_gev = experimental.value * _HBAR_GEV_PER_PS
    experimental_sigma_delta_m_gev = experimental.uncertainty * _HBAR_GEV_PER_PS
    sm_delta_m_gev = core_inputs["delta_m_bs_sm_gev"]
    sm_theory_sigma_delta_m_gev = (
        abs(sm_delta_m_gev) * 2.0 * flag_bmixing.f_bs_sqrt_bhat_relative_sigma
    )
    combined_sigma_delta_m_gev = math.sqrt(
        experimental_sigma_delta_m_gev * experimental_sigma_delta_m_gev
        + sm_theory_sigma_delta_m_gev * sm_theory_sigma_delta_m_gev
    )
    residual_delta_m_gev = abs(experimental_delta_m_gev - sm_delta_m_gev)
    central_budget = residual_delta_m_gev / 2.0
    loose_budget = (residual_delta_m_gev + combined_sigma_delta_m_gev) / 2.0

    if central_budget <= 0.0 or loose_budget <= 0.0:
        raise AnchorError(
            f"{process_id}: B_s NP budgets must be positive "
            f"(central={central_budget}, loose={loose_budget})"
        )

    return BsMixingBudgetBand(
        doc_citation=_BUDGET_DOC_CITATION,
        central_budget=central_budget,
        loose_budget=loose_budget,
        hard_veto_budget=loose_budget,
        experimental_delta_m_gev=experimental_delta_m_gev,
        experimental_delta_m_ps_inv=float(experimental.value),
        experimental_m12_gev=experimental_delta_m_gev / 2.0,
        sm_delta_m_gev=sm_delta_m_gev,
        sm_delta_m_ps_inv=sm_delta_m_gev / _HBAR_GEV_PER_PS,
        sm_m12_gev=sm_delta_m_gev / 2.0,
        experimental_sigma_delta_m_gev=experimental_sigma_delta_m_gev,
        sm_theory_sigma_delta_m_gev=sm_theory_sigma_delta_m_gev,
        combined_sigma_delta_m_gev=combined_sigma_delta_m_gev,
        flag_f_bs_sqrt_bhat_relative_sigma=(
            flag_bmixing.f_bs_sqrt_bhat_relative_sigma
        ),
        core_delta_m_exp_gev=core_inputs["delta_m_bs_exp_gev"],
        core_delta_m_sm_gev=core_inputs["delta_m_bs_sm_gev"],
        core_default_m12_budget=core_inputs["core_m12_budget_gev"],
        core_budget_policy_id=str(core_inputs["core_budget_policy_id"]),
        core_budget_confidence_level=str(core_inputs["core_budget_confidence_level"]),
        core_legacy_m12_budget=core_inputs["legacy_full_delta_m_m12_budget_gev"],
    )


def _load_bs_mixing_anchor(process_id: str) -> BsMixingAnchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_EXPERIMENTAL_ANCHOR_CANDIDATES,
    )
    _validate_delta_ms_anchor(experimental, process_id=process_id)
    flag_bmixing = _load_flag_bmixing_inputs(process_id)
    anchor = BsMixingAnchor(
        experimental=experimental,
        flag_bmixing=flag_bmixing,
        budget_band=_build_budget_band(
            process_id=process_id,
            experimental=experimental,
            flag_bmixing=flag_bmixing,
        ),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(
            f"{process_id}: B_s mixing NP budget must be positive, got {anchor.budget}"
        )
    return anchor


def _complex_mapping(values: Mapping[str, complex]) -> dict[str, complex]:
    return {str(key): complex(value) for key, value in values.items()}


@register
class Constraint:
    """Catalogued Delta m_s / B_s Delta F=2 constraint (process_id B003)."""

    process_id = "B003"
    severity = Severity.HARD
    observable = _EXPECTED_OBSERVABLE

    def __init__(self) -> None:
        self.anchor = _load_bs_mixing_anchor(self.process_id)
        self.sm_value = self.anchor.sm_value

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.anchor.budget_band.sm_delta_m_gev),
                experimental=float(self.anchor.budget_band.experimental_delta_m_gev),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; Delta m_s constraint "
                    "was not evaluated. sm_prediction/experimental are "
                    "Delta m_s in GeV; predicted/budget use |M12^NP| in GeV."
                ),
                diagnostics={"missing_extra": _REQUIRED_EXTRA},
            )

        wilsons = bs_mixing_wilsons_from_couplings(couplings)
        result = bs_mixing_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
            m12_np_budget=self.anchor.budget,
        )

        predicted = float(result.abs_m12_np)
        ratio = float(result.ratio_to_budget)
        budget = float(result.budget)
        core_budget = float(self.anchor.budget_band.core_legacy_m12_budget)
        core_default_budget = float(self.anchor.budget_band.core_default_m12_budget)

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(result.passes),
            predicted=predicted,
            sm_prediction=float(self.anchor.budget_band.sm_delta_m_gev),
            experimental=float(self.anchor.budget_band.experimental_delta_m_gev),
            ratio=ratio,
            budget=budget,
            notes=(
                "|M12^NP(B_s)| is evaluated with Wilsons QCD-evolved to m_b=4.18 GeV; "
                "predicted and budget are |M12^NP| in GeV, while "
                "sm_prediction/experimental are Delta m_s in GeV; "
                "the HARD budget is the one-sigma SM-vs-HFLAV M12 room using "
                f"FLAG f_Bs sqrt(Bhat_Bs), cited in {_BUDGET_DOC_CITATION}."
            ),
            diagnostics={
                "abs_m12_np_gev": predicted,
                "delta_m_np_gev": 2.0 * predicted,
                "qcd_running_applied": True,
                "hadronic_scale_gev": _MU_HAD_GEV,
                "matching_scale_gev": float(wilsons.matching_scale),
                "m_kk_gev": float(wilsons.M_KK),
                "left_sb_coupling": complex(wilsons.left_coupling),
                "right_sb_coupling": complex(wilsons.right_coupling),
                "wilson_coefficients": _complex_mapping(wilsons.wilsons),
                "system": result.system,
                "budget_doc_citation": self.anchor.budget_band.doc_citation,
                "budget_policy_id": self.anchor.budget_band.core_budget_policy_id,
                "confidence_level": self.anchor.budget_band.core_budget_confidence_level,
                "central_np_budget": self.anchor.budget_band.central_budget,
                "loose_band_np_budget": self.anchor.budget_band.loose_budget,
                "hard_veto_np_budget": self.anchor.budget_band.hard_veto_budget,
                "experimental_delta_m_ps_inv": (
                    self.anchor.budget_band.experimental_delta_m_ps_inv
                ),
                "experimental_delta_m_gev": (
                    self.anchor.budget_band.experimental_delta_m_gev
                ),
                "experimental_m12_gev": (
                    self.anchor.budget_band.experimental_m12_gev
                ),
                "sm_delta_m_ps_inv": self.anchor.budget_band.sm_delta_m_ps_inv,
                "sm_delta_m_gev": self.anchor.budget_band.sm_delta_m_gev,
                "sm_m12_gev": self.anchor.budget_band.sm_m12_gev,
                "budget_combined_sigma_delta_m_gev": (
                    self.anchor.budget_band.combined_sigma_delta_m_gev
                ),
                "budget_sm_theory_sigma_delta_m_gev": (
                    self.anchor.budget_band.sm_theory_sigma_delta_m_gev
                ),
                "budget_experimental_sigma_delta_m_gev": (
                    self.anchor.budget_band.experimental_sigma_delta_m_gev
                ),
                "flag_f_bs_sqrt_bhat_relative_sigma": (
                    self.anchor.budget_band.flag_f_bs_sqrt_bhat_relative_sigma
                ),
                "flag_f_bs_sqrt_bhat_bs_mev": (
                    self.anchor.flag_bmixing.f_bs_sqrt_bhat_bs_value_mev
                ),
                "flag_f_bs_sqrt_bhat_bs_uncertainty_mev": (
                    self.anchor.flag_bmixing.f_bs_sqrt_bhat_bs_uncertainty_mev
                ),
                "core_delta_m_exp_gev": self.anchor.budget_band.core_delta_m_exp_gev,
                "core_delta_m_sm_gev": self.anchor.budget_band.core_delta_m_sm_gev,
                "core_default_m12_budget": core_default_budget,
                "core_legacy_m12_budget": core_budget,
                "core_legacy_ratio_to_budget": (
                    predicted / core_budget if core_budget > 0.0 else float("inf")
                ),
                "core_default_ratio_to_budget": (
                    predicted / core_default_budget
                    if core_default_budget > 0.0
                    else float("inf")
                ),
                "core_budget_policy": (
                    "Delta F=2 core default is the B003 one-sigma SM-vs-HFLAV "
                    "M12 room; the legacy diagnostic is Delta m_exp/2 only."
                ),
                "experimental_block": self.anchor.experimental.block_key,
                "auxiliary_theory_block": self.anchor.flag_bmixing.block_key,
                "flag_n_f": self.anchor.flag_bmixing.n_f,
                "flag_source_url": self.anchor.flag_bmixing.source_url,
            },
        )
