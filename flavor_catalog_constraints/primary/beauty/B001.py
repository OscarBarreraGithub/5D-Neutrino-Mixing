"""B001 - neutral B_d mixing mass difference.

Physics
-------
``Delta m_d`` is the mass splitting in ``B_d0-B_dbar0`` mixing.  The
new-physics contribution is evaluated with the audited Delta F = 2 core as
``|M12^NP|`` after QCD-running the Wilson coefficients to ``mu_had = 2 GeV``.

Severity
--------
HARD.  The oscillation frequency is observed, and the RS new-physics
contribution must fit inside an uncertainty-aware SM-vs-experiment room.  The
core's legacy B_d convention uses
``max(Delta m_exp / 2, |Delta m_exp - Delta m_SM| / 2)``; this plugin keeps
that core budget as a diagnostic and applies the requested catalog budget,

    (|Delta m_exp - Delta m_SM| + sigma_combined) / 2,

where ``sigma_combined`` combines the HFLAV/PDG experimental uncertainty with
the larger HPQCD 2019 asymmetric SM uncertainty in quadrature.

Catalog sidecar
---------------
``flavor_catalog/processes/beauty/B001.yaml`` is the source of truth for the
experimental ``Delta m_d`` average, the HPQCD SM-prediction uncertainty, and
the code-side GeV conversion anchor.  The SM central value is read through the
Delta F = 2 adapter from the same core evaluator convention used for the
Wilson/matrix-element calculation.
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
    load_full_yaml,
    load_pdg_block,
)
from flavor_catalog_constraints.base import ConstraintResult, ParameterPoint, Severity
from flavor_catalog_constraints.physics_adapters.deltaf2 import (
    bd_mixing_core_inputs,
    bd_mixing_from_wilsons_with_running,
    bd_mixing_wilsons_from_couplings,
)
from flavor_catalog_constraints.registry import register

_FAMILY = "beauty"
_REQUIRED_EXTRA = "quark_mass_basis_couplings"

_EXPERIMENTAL_ANCHOR_CANDIDATES = ("canonical_experimental_average",)
_SM_PREDICTION_ANCHOR_CANDIDATES = ("standard_model_prediction",)
_AUXILIARY_CODE_CANDIDATES = ("deltaf2_bd_constants",)
_BUDGET_DOC_CITATION = (
    "flavor_catalog/processes/beauty/B001.yaml:"
    "pdg_or_equivalent.standard_model_prediction; "
    "docs/quark_scan_assumptions_compact.tex:445-472; "
    "docs/audits/bag_param_inventory.md:36-37; "
    "HPQCD 2019 arXiv:1907.01025"
)
_MU_HAD_GEV = 2.0


@dataclass(frozen=True)
class BdCodeInputs:
    """Typed view over B001's auxiliary code-input provenance block."""

    block_key: str
    source: str | None
    lines: tuple[str, ...]
    delta_m_bd_exp_gev: float
    f_bd_gev: float | None
    m_bd_gev: float | None
    b_1_bd: float | None
    b_4_bd: float | None
    b_5_bd: float | None


@dataclass(frozen=True)
class BdMixingBudgetBand:
    """Uncertainty-aware B_d mixing NP budget for the B001 HARD veto."""

    doc_citation: str
    central_budget: float
    loose_budget: float
    tight_budget: float
    hard_veto_budget: float
    combined_delta_m_sigma_gev: float
    experimental_delta_m_sigma_gev: float
    sm_delta_m_sigma_gev: float
    experimental_delta_m_gev: float
    sm_delta_m_gev: float
    experimental_delta_m_ps_inv: float
    sm_delta_m_ps_inv: float
    experimental_delta_m_sigma_ps_inv: float
    sm_delta_m_sigma_ps_inv: float
    sm_delta_m_sigma_policy: str | None
    gev_per_ps_inverse: float
    core_default_budget: float


@dataclass(frozen=True)
class BdMixingAnchor:
    """Typed B001 anchor: experiment, SM reference, code inputs, and budget."""

    experimental: Anchor
    standard_model: Anchor
    code_inputs: BdCodeInputs
    budget_band: BdMixingBudgetBand

    @property
    def value(self) -> float:
        """Experimental central value in the YAML units, ps^-1."""
        return self.experimental.value

    @property
    def uncertainty(self) -> float | None:
        """Experimental uncertainty in the YAML units, ps^-1."""
        return self.experimental.uncertainty

    @property
    def source_url(self) -> str | None:
        return self.experimental.source_url

    @property
    def delta_m_experimental_gev(self) -> float:
        return self.budget_band.experimental_delta_m_gev

    @property
    def sm_value(self) -> float:
        """SM central prediction for Delta m_d in GeV."""
        return self.budget_band.sm_delta_m_gev

    @property
    def central_budget(self) -> float:
        """Central residual room for ``|M12^NP|`` in GeV."""
        return self.budget_band.central_budget

    @property
    def budget(self) -> float:
        """Uncertainty-aware HARD veto budget for ``|M12^NP|`` in GeV."""
        return self.budget_band.hard_veto_budget


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _required_float(value: Any, *, process_id: str, field_name: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorError(
            f"{process_id}: B_d mixing field {field_name!r}={value!r} "
            "is not numeric"
        ) from exc


def _optional_float(value: Any, *, process_id: str, field_name: str) -> float | None:
    if value is None:
        return None
    return _required_float(value, process_id=process_id, field_name=field_name)


def _load_code_inputs(process_id: str) -> BdCodeInputs:
    data = load_full_yaml(process_id, family=_FAMILY)
    aux = data.get("auxiliary_code_inputs")
    if not isinstance(aux, Mapping):
        raise AnchorError(
            f"{process_id}: missing mapping top-level 'auxiliary_code_inputs'"
        )
    sub = find_block(aux, _AUXILIARY_CODE_CANDIDATES, process_id=process_id)
    block_key = next(key for key in _AUXILIARY_CODE_CANDIDATES if key in aux)
    values = sub.get("values")
    if not isinstance(values, Mapping):
        raise AnchorError(
            f"{process_id}: auxiliary_code_inputs.{block_key} must provide values"
        )
    lines = sub.get("lines", ())
    if not isinstance(lines, list | tuple):
        lines = ()
    return BdCodeInputs(
        block_key=block_key,
        source=_optional_str(sub.get("source")),
        lines=tuple(str(line) for line in lines),
        delta_m_bd_exp_gev=_required_float(
            values.get("DELTA_M_BD_EXP_GeV"),
            process_id=process_id,
            field_name="DELTA_M_BD_EXP_GeV",
        ),
        f_bd_gev=_optional_float(
            values.get("F_BD_GeV"),
            process_id=process_id,
            field_name="F_BD_GeV",
        ),
        m_bd_gev=_optional_float(
            values.get("M_BD_GeV"),
            process_id=process_id,
            field_name="M_BD_GeV",
        ),
        b_1_bd=_optional_float(
            values.get("B_1_BD"),
            process_id=process_id,
            field_name="B_1_BD",
        ),
        b_4_bd=_optional_float(
            values.get("B_4_BD"),
            process_id=process_id,
            field_name="B_4_BD",
        ),
        b_5_bd=_optional_float(
            values.get("B_5_BD"),
            process_id=process_id,
            field_name="B_5_BD",
        ),
    )


def _load_sm_prediction_sub(process_id: str) -> Mapping[str, Any]:
    pdg_block = load_pdg_block(process_id, family=_FAMILY)
    return find_block(
        pdg_block,
        _SM_PREDICTION_ANCHOR_CANDIDATES,
        process_id=process_id,
    )


def _validated_sm_theory_sigma_ps_inv(
    *,
    process_id: str,
    standard_model: Anchor,
    standard_model_sub: Mapping[str, Any],
) -> float:
    block_key = standard_model.block_key
    if standard_model.uncertainty is None or standard_model.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: SM Delta m_d uncertainty anchor is required "
            f"in {block_key}.uncertainty"
        )
    if standard_model.value <= 0.0:
        raise AnchorError(
            f"{process_id}: SM Delta m_d central anchor must be positive, "
            f"got {standard_model.value}"
        )
    if standard_model.units != "ps^-1":
        raise AnchorError(
            f"{process_id}: SM Delta m_d uncertainty anchor must use ps^-1, "
            f"got units={standard_model.units!r}"
        )

    asymmetric_plus = abs(
        _required_float(
            standard_model_sub.get("asymmetric_uncertainty_plus"),
            process_id=process_id,
            field_name=f"{block_key}.asymmetric_uncertainty_plus",
        )
    )
    asymmetric_minus = abs(
        _required_float(
            standard_model_sub.get("asymmetric_uncertainty_minus"),
            process_id=process_id,
            field_name=f"{block_key}.asymmetric_uncertainty_minus",
        )
    )
    second_uncertainty = abs(
        _required_float(
            standard_model_sub.get("second_uncertainty"),
            process_id=process_id,
            field_name=f"{block_key}.second_uncertainty",
        )
    )
    expected_sigma = math.hypot(
        max(asymmetric_plus, asymmetric_minus),
        second_uncertainty,
    )
    sigma = float(standard_model.uncertainty)
    if not math.isclose(sigma, expected_sigma, rel_tol=5.0e-3, abs_tol=5.0e-4):
        raise AnchorError(
            f"{process_id}: SM Delta m_d uncertainty anchor {sigma} does not "
            f"match component quadrature {expected_sigma}"
        )
    return sigma


def _build_budget_band(
    *,
    process_id: str,
    experimental: Anchor,
    standard_model: Anchor,
    standard_model_sub: Mapping[str, Any],
    code_inputs: BdCodeInputs,
) -> BdMixingBudgetBand:
    if experimental.uncertainty is None or experimental.uncertainty <= 0.0:
        raise AnchorError(
            f"{process_id}: Delta m_d experimental uncertainty is required "
            "for the uncertainty-aware budget band"
        )
    if experimental.value <= 0.0 or code_inputs.delta_m_bd_exp_gev <= 0.0:
        raise AnchorError(f"{process_id}: Delta m_d anchors must be positive")

    core_inputs = bd_mixing_core_inputs()
    core_exp_gev = core_inputs["delta_m_bd_exp_gev"]
    if not math.isclose(
        code_inputs.delta_m_bd_exp_gev,
        core_exp_gev,
        rel_tol=1.0e-12,
        abs_tol=0.0,
    ):
        raise AnchorError(
            f"{process_id}: YAML DELTA_M_BD_EXP_GeV="
            f"{code_inputs.delta_m_bd_exp_gev} does not match core "
            f"DELTA_M_BD_EXP={core_exp_gev}"
        )

    gev_per_ps_inverse = code_inputs.delta_m_bd_exp_gev / experimental.value
    experimental_sigma_gev = float(experimental.uncertainty) * gev_per_ps_inverse
    sm_sigma_ps_inv = _validated_sm_theory_sigma_ps_inv(
        process_id=process_id,
        standard_model=standard_model,
        standard_model_sub=standard_model_sub,
    )
    sm_sigma_gev = sm_sigma_ps_inv * gev_per_ps_inverse
    combined_sigma = math.sqrt(
        experimental_sigma_gev * experimental_sigma_gev
        + sm_sigma_gev * sm_sigma_gev
    )
    sm_delta_m_gev = core_inputs["delta_m_bd_sm_gev"]
    central_gap = abs(code_inputs.delta_m_bd_exp_gev - sm_delta_m_gev)
    central_budget = central_gap / 2.0
    loose_budget = (central_gap + combined_sigma) / 2.0
    tight_budget = max(central_gap - combined_sigma, experimental_sigma_gev) / 2.0

    if central_budget <= 0.0 or loose_budget <= 0.0 or tight_budget <= 0.0:
        raise AnchorError(
            f"{process_id}: B_d NP budgets must be positive "
            f"(central={central_budget}, loose={loose_budget}, tight={tight_budget})"
        )

    return BdMixingBudgetBand(
        doc_citation=_BUDGET_DOC_CITATION,
        central_budget=central_budget,
        loose_budget=loose_budget,
        tight_budget=tight_budget,
        hard_veto_budget=loose_budget,
        combined_delta_m_sigma_gev=combined_sigma,
        experimental_delta_m_sigma_gev=experimental_sigma_gev,
        sm_delta_m_sigma_gev=sm_sigma_gev,
        experimental_delta_m_gev=code_inputs.delta_m_bd_exp_gev,
        sm_delta_m_gev=sm_delta_m_gev,
        experimental_delta_m_ps_inv=experimental.value,
        sm_delta_m_ps_inv=sm_delta_m_gev / gev_per_ps_inverse,
        experimental_delta_m_sigma_ps_inv=float(experimental.uncertainty),
        sm_delta_m_sigma_ps_inv=sm_sigma_ps_inv,
        sm_delta_m_sigma_policy=_optional_str(
            standard_model_sub.get("uncertainty_policy")
        ),
        gev_per_ps_inverse=gev_per_ps_inverse,
        core_default_budget=core_inputs["core_m12_budget_gev"],
    )


def _load_bd_mixing_anchor(process_id: str) -> BdMixingAnchor:
    experimental = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_EXPERIMENTAL_ANCHOR_CANDIDATES,
    )
    standard_model = load_anchor(
        process_id,
        family=_FAMILY,
        candidates=_SM_PREDICTION_ANCHOR_CANDIDATES,
    )
    standard_model_sub = _load_sm_prediction_sub(process_id)
    code_inputs = _load_code_inputs(process_id)
    anchor = BdMixingAnchor(
        experimental=experimental,
        standard_model=standard_model,
        code_inputs=code_inputs,
        budget_band=_build_budget_band(
            process_id=process_id,
            experimental=experimental,
            standard_model=standard_model,
            standard_model_sub=standard_model_sub,
            code_inputs=code_inputs,
        ),
    )
    if anchor.budget <= 0.0:
        raise AnchorError(
            f"{process_id}: B_d mixing NP budget must be positive, got {anchor.budget}"
        )
    return anchor


def _complex_mapping(values: Mapping[str, complex]) -> dict[str, complex]:
    return {str(key): complex(value) for key, value in values.items()}


@register
class Constraint:
    """Catalogued Delta m_d Delta F=2 constraint (process_id B001)."""

    process_id = "B001"
    severity = Severity.HARD
    observable = "Delta m_d"

    def __init__(self) -> None:
        self.anchor = _load_bd_mixing_anchor(self.process_id)
        self.sm_value = self.anchor.sm_value

    def evaluate(self, point: ParameterPoint) -> ConstraintResult:
        couplings = point.get_extra(_REQUIRED_EXTRA)
        if couplings is None:
            return ConstraintResult(
                process_id=self.process_id,
                severity=self.severity,
                passes=True,
                sm_prediction=float(self.anchor.sm_value),
                experimental=float(self.anchor.delta_m_experimental_gev),
                budget=float(self.anchor.budget),
                notes=(
                    f"extra {_REQUIRED_EXTRA!r} absent; Delta m_d constraint "
                    "was not evaluated."
                ),
                diagnostics={"missing_extra": _REQUIRED_EXTRA},
            )

        wilsons = bd_mixing_wilsons_from_couplings(couplings)
        core_result = bd_mixing_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
        )
        result = bd_mixing_from_wilsons_with_running(
            wilsons,
            mu_had=_MU_HAD_GEV,
            m12_np_budget=self.anchor.budget,
        )

        predicted = float(result.abs_m12_np)
        ratio = float(result.ratio_to_budget)
        budget = float(result.budget)

        return ConstraintResult(
            process_id=self.process_id,
            severity=self.severity,
            passes=bool(result.passes),
            predicted=predicted,
            sm_prediction=float(self.anchor.sm_value),
            experimental=float(self.anchor.delta_m_experimental_gev),
            ratio=ratio,
            budget=budget,
            notes=(
                "|M12_Bd^NP| is evaluated by the QCD-running Delta F=2 core "
                "at 2 GeV; the HARD budget is the uncertainty-aware "
                "SM-vs-experiment room documented in diagnostics."
            ),
            diagnostics={
                "abs_m12_np_gev": predicted,
                "qcd_running_applied": True,
                "hadronic_scale_gev": _MU_HAD_GEV,
                "matching_scale_gev": float(wilsons.matching_scale),
                "m_kk_gev": float(wilsons.M_KK),
                "core_input_key": wilsons.input.key,
                "down_sector_indices": tuple(wilsons.input.generations),
                "left_db_coupling": complex(wilsons.left_coupling),
                "right_db_coupling": complex(wilsons.right_coupling),
                "wilson_coefficients": _complex_mapping(wilsons.wilsons),
                "budget_doc_citation": self.anchor.budget_band.doc_citation,
                "budget_construction": (
                    "(|Delta m_exp - Delta m_SM| + "
                    "sqrt(sigma_exp^2 + sigma_SM^2)) / 2"
                ),
                "sm_theory_sigma_policy": (
                    self.anchor.budget_band.sm_delta_m_sigma_policy
                ),
                "central_np_budget": self.anchor.budget_band.central_budget,
                "tight_band_np_budget": self.anchor.budget_band.tight_budget,
                "loose_band_np_budget": self.anchor.budget_band.loose_budget,
                "hard_veto_np_budget": self.anchor.budget_band.hard_veto_budget,
                "budget_combined_delta_m_sigma_gev": (
                    self.anchor.budget_band.combined_delta_m_sigma_gev
                ),
                "budget_experimental_delta_m_sigma_gev": (
                    self.anchor.budget_band.experimental_delta_m_sigma_gev
                ),
                "budget_sm_delta_m_sigma_gev": (
                    self.anchor.budget_band.sm_delta_m_sigma_gev
                ),
                "experimental_delta_m_ps_inv": (
                    self.anchor.budget_band.experimental_delta_m_ps_inv
                ),
                "sm_delta_m_ps_inv": self.anchor.budget_band.sm_delta_m_ps_inv,
                "experimental_delta_m_sigma_ps_inv": (
                    self.anchor.budget_band.experimental_delta_m_sigma_ps_inv
                ),
                "sm_delta_m_sigma_ps_inv": (
                    self.anchor.budget_band.sm_delta_m_sigma_ps_inv
                ),
                "gev_per_ps_inverse": self.anchor.budget_band.gev_per_ps_inverse,
                "core_default_budget": self.anchor.budget_band.core_default_budget,
                "core_default_ratio": float(core_result.ratio_to_budget),
                "core_default_passes": bool(core_result.passes),
                "experimental_block": self.anchor.experimental.block_key,
                "sm_prediction_block": self.anchor.standard_model.block_key,
                "sm_prediction_source": self.anchor.standard_model.source,
                "sm_prediction_source_url": self.anchor.standard_model.source_url,
                "sm_prediction_snapshot_path": (
                    self.anchor.standard_model.snapshot_path
                ),
                "auxiliary_code_block": self.anchor.code_inputs.block_key,
                "auxiliary_code_source": self.anchor.code_inputs.source,
                "auxiliary_code_lines": self.anchor.code_inputs.lines,
                "f_bd_gev": self.anchor.code_inputs.f_bd_gev,
                "m_bd_gev": self.anchor.code_inputs.m_bd_gev,
                "b_1_bd": self.anchor.code_inputs.b_1_bd,
                "b_4_bd": self.anchor.code_inputs.b_4_bd,
                "b_5_bd": self.anchor.code_inputs.b_5_bd,
            },
        )
