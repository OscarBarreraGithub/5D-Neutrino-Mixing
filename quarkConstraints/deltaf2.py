"""Tree-level KK-gluon-inspired ``Delta F = 2`` benchmark observables.

This module intentionally implements a compact, repo-owned v1 matching layer:

- exact fitted quark mass-basis couplings come from ``quarkConstraints.couplings``
- KK-gluon exchange is represented by a small fixed operator basis
- hadronic/bound inputs live in a deterministic in-repo bundle

The goal is to turn the current MFV fit machinery into a first exclusion slice,
not to provide a general-purpose EFT/RG package.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

from .couplings import (
    COUPLING_POLICY_PERTURBATIVE_4D_LEGACY,
    OPERATOR_CONVENTION_PERTURBATIVE_4D_LEGACY,
    QuarkMassBasisCouplings,
    compute_quark_kk_gluon_couplings,
)
from .fit import QuarkFitResult
from .qcd_running import evolve_deltaf2_wilsons
from .scales import DEFAULT_QUARK_TARGET_SCALE_GEV, DEFAULT_QUARK_XI_KK

DELTA_F2_MODEL_V1 = OPERATOR_CONVENTION_PERTURBATIVE_4D_LEGACY
DELTA_F2_INPUT_BUNDLE_V1 = "deltaf2_inputs_mu3tev_v1"
DELTA_F2_OPERATOR_CONVENTION = DELTA_F2_MODEL_V1
DELTA_F2_INPUT_BUNDLE = DELTA_F2_INPUT_BUNDLE_V1
DELTA_F2_CATALOG_CL_NOTE = (
    "NOTE: Delta-F=2 epsilon_K uses a 68.27% one-sigma sensitivity budget; "
    "catalog-wide single-CL harmonization is a separate policy decision."
)
DELTA_F2_RUNNING_ORDER = "LO"
DELTA_F2_LO_RUNNING_BIAS_NOTE = (
    "Delta-F=2 running is LO; the LR (C4) enhancement is likely low by "
    "~10-15% in M_KK vs NLO -- eps_K floors are therefore mildly "
    "CONSERVATIVE (understated) in this direction."
)
DELTA_F2_MU_HAD_KAON_GEV = 3.0
DELTA_F2_MU_HAD_B_GEV = 4.18
DELTA_F2_MU_HAD_D_GEV = 3.0

# HFLAV/B001-B003 mass-splitting budget inputs. These are intentionally kept
# near DEFAULT_DELTA_F2_INPUTS_V1 because the input bundle's ``bound`` field is
# a public diagnostic, while the hadronic constants live further below.
HBAR_GEV_PER_PS = 6.582119569e-13

DELTA_M_BD_EXP = 3.334e-13
DELTA_M_BD_SM = 3.6e-13
DELTA_M_BD_EXP_PS_INV = 0.5069
DELTA_M_BD_EXP_SIGMA_PS_INV = 0.0019
DELTA_M_BD_SM_SIGMA_PS_INV = 0.062
DELTA_M_BD_GEV_PER_PS_INV = DELTA_M_BD_EXP / DELTA_M_BD_EXP_PS_INV
DELTA_M_BD_EXP_SIGMA = DELTA_M_BD_EXP_SIGMA_PS_INV * DELTA_M_BD_GEV_PER_PS_INV
DELTA_M_BD_SM_SIGMA = DELTA_M_BD_SM_SIGMA_PS_INV * DELTA_M_BD_GEV_PER_PS_INV
DELTA_M_BD_COMBINED_SIGMA = math.hypot(
    DELTA_M_BD_EXP_SIGMA,
    DELTA_M_BD_SM_SIGMA,
)
DELTA_M_BD_CENTRAL_BUDGET = abs(DELTA_M_BD_EXP - DELTA_M_BD_SM) / 2.0
DELTA_M_BD_BUDGET = (
    abs(DELTA_M_BD_EXP - DELTA_M_BD_SM) + DELTA_M_BD_COMBINED_SIGMA
) / 2.0

DELTA_M_BS_EXP_PS_INV = 17.766
DELTA_M_BS_EXP_SIGMA_PS_INV = 0.006
DELTA_M_BS_EXP = DELTA_M_BS_EXP_PS_INV * HBAR_GEV_PER_PS
DELTA_M_BS_EXP_SIGMA = DELTA_M_BS_EXP_SIGMA_PS_INV * HBAR_GEV_PER_PS
DELTA_M_BS_SM = 1.17e-11
DELTA_M_BS_FLAG_F_BS_SQRT_BHAT_BS_MEV = 256.1
DELTA_M_BS_FLAG_F_BS_SQRT_BHAT_BS_SIGMA_MEV = 5.7
DELTA_M_BS_SM_SIGMA = (
    abs(DELTA_M_BS_SM)
    * 2.0
    * DELTA_M_BS_FLAG_F_BS_SQRT_BHAT_BS_SIGMA_MEV
    / DELTA_M_BS_FLAG_F_BS_SQRT_BHAT_BS_MEV
)
DELTA_M_BS_COMBINED_SIGMA = math.hypot(
    DELTA_M_BS_EXP_SIGMA,
    DELTA_M_BS_SM_SIGMA,
)
DELTA_M_BS_CENTRAL_BUDGET = abs(DELTA_M_BS_EXP - DELTA_M_BS_SM) / 2.0
DELTA_M_BS_BUDGET = (
    abs(DELTA_M_BS_EXP - DELTA_M_BS_SM) + DELTA_M_BS_COMBINED_SIGMA
) / 2.0

DELTA_M_D_EXP = 6.25e-15
DELTA_M_D_BUDGET = DELTA_M_D_EXP / 2.0


def delta_f2_default_mu_had_for_system(key: str) -> float:
    """Return the default hadronic RG endpoint for one Delta-F=2 system."""

    if key == "epsilon_k":
        return DELTA_F2_MU_HAD_KAON_GEV
    if key in {"b_d", "b_s"}:
        return DELTA_F2_MU_HAD_B_GEV
    if key == "d":
        return DELTA_F2_MU_HAD_D_GEV
    raise KeyError(key)


def _resolve_mu_had_for_system(
    key: str,
    mu_had: float | Mapping[str, float] | None,
) -> float:
    if mu_had is None:
        return delta_f2_default_mu_had_for_system(key)
    if isinstance(mu_had, Mapping):
        return float(mu_had.get(key, delta_f2_default_mu_had_for_system(key)))
    return float(mu_had)


@dataclass(frozen=True)
class DeltaF2Input:
    """Repo-owned v1 input for one neutral-meson mixing system."""

    key: str
    display_name: str
    column_name: str
    reject_reason: str
    sector: str
    generations: tuple[int, int]
    bound: float
    ll_weight: float = 1.0
    rr_weight: float = 1.0
    lr1_weight: float = 7.0
    lr2_weight: float = 2.0
    reference_scale: float = DEFAULT_QUARK_TARGET_SCALE_GEV
    note: str = ""

    def __post_init__(self) -> None:
        if self.sector not in {"down", "up"}:
            raise ValueError("sector must be 'down' or 'up'")
        if len(self.generations) != 2:
            raise ValueError("generations must be a pair of flavor indices")
        i, j = self.generations
        if i == j or not (0 <= i < 3) or not (0 <= j < 3):
            raise ValueError("generations must be distinct indices in {0,1,2}")
        if self.bound <= 0.0:
            raise ValueError("bound must be positive")
        if self.reference_scale <= 0.0:
            raise ValueError("reference_scale must be positive")

    @property
    def system(self) -> str:
        aliases = {
            "epsilon_k": "K",
            "b_d": "B_d",
            "b_s": "B_s",
            "d": "D",
        }
        return aliases.get(self.key, self.display_name)

    @property
    def flavor_indices(self) -> tuple[int, int]:
        return self.generations


@dataclass(frozen=True)
class DeltaF2WilsonCoefficients:
    """Wilson coefficients for one neutral-meson system."""

    input: DeltaF2Input
    M_KK: float
    matching_scale: float
    left_coupling: complex
    right_coupling: complex
    c1_vll: complex
    c1_vrr: complex
    c4_lr: complex
    c5_lr: complex
    lambda_ir_gev: float | None = None
    m_kk_physical_gev: float | None = None
    mass_convention_id: str | None = None
    g_s_4d: float | None = None
    g_eff: float | None = None
    g_s_multiplier: float | None = None
    coupling_policy_id: str = COUPLING_POLICY_PERTURBATIVE_4D_LEGACY
    operator_convention_id: str = DELTA_F2_MODEL_V1

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {
            "C1_VLL": self.c1_vll,
            "C1_VRR": self.c1_vrr,
            "C4_LR": self.c4_lr,
            "C5_LR": self.c5_lr,
        }


@dataclass(frozen=True)
class EpsilonKBudgetPolicy:
    """Shared epsilon_K NP-budget policy for core and catalog consumers."""

    policy_id: str
    confidence_level: str
    doc_citation: str
    experimental_value: float
    sm_value: float
    central_residual_signed: float
    central_budget: float
    sigma_bgs: float
    sigma_exp: float
    primary_combined_sigma: float
    sm_choice_sensitivity: float
    lower_signed_edge: float
    upper_signed_edge: float
    budget_lowers_epsilon_k: float
    budget_raises_epsilon_k: float
    tight_budget: float
    loose_budget: float

    def selected_signed_budget(self, epsilon_k_np_signed: float) -> tuple[float, str]:
        """Return the direction-aware HARD budget for a signed NP shift."""

        if epsilon_k_np_signed >= 0.0:
            return self.budget_raises_epsilon_k, "raises_epsilon_k"
        return self.budget_lowers_epsilon_k, "lowers_epsilon_k"

    def as_diagnostics(self) -> dict[str, Any]:
        """Diagnostic payload shared by core and catalog reports."""

        return {
            "budget_policy_id": self.policy_id,
            "confidence_level": self.confidence_level,
            "budget_doc_citation": self.doc_citation,
            "central_np_budget": float(self.central_budget),
            "central_residual_signed": float(self.central_residual_signed),
            "tight_band_np_budget": float(self.tight_budget),
            "loose_band_np_budget": float(self.loose_budget),
            "signed_lower_edge": float(self.lower_signed_edge),
            "signed_upper_edge": float(self.upper_signed_edge),
            "budget_lowers_epsilon_k": float(self.budget_lowers_epsilon_k),
            "budget_raises_epsilon_k": float(self.budget_raises_epsilon_k),
            "hard_veto_np_budget": "direction_selected",
            "budget_combined_sigma": float(self.primary_combined_sigma),
            "budget_primary_combined_sigma": float(self.primary_combined_sigma),
            "budget_sm_theory_sigma": float(self.sigma_bgs),
            "budget_experimental_sigma": float(self.sigma_exp),
            "sm_choice_sensitivity": float(self.sm_choice_sensitivity),
            "budget_sm_choice_sensitivity": float(self.sm_choice_sensitivity),
            "sm_choice_sensitivity_in_hard_gate": False,
            "catalog_confidence_level_note": DELTA_F2_CATALOG_CL_NOTE,
            "running_order": DELTA_F2_RUNNING_ORDER,
            "running_bias_note": DELTA_F2_LO_RUNNING_BIAS_NOTE,
        }


@dataclass(frozen=True)
class DeltaMNPBudgetPolicy:
    """Named NP-budget policy for a neutral-meson mass splitting."""

    policy_id: str
    confidence_level: str
    doc_citation: str
    system: str
    budget: float
    construction: str
    experimental_delta_m: float
    sm_delta_m: float | None
    central_delta_m_residual: float | None
    experimental_sigma_delta_m: float | None
    sm_sigma_delta_m: float | None
    combined_sigma_delta_m: float | None
    legacy_full_delta_m_m12_budget: float
    note: str

    def as_diagnostics(self) -> dict[str, Any]:
        return {
            "budget_policy_id": self.policy_id,
            "confidence_level": self.confidence_level,
            "budget_doc_citation": self.doc_citation,
            "budget_construction": self.construction,
            "hard_veto_np_budget": float(self.budget),
            "central_np_budget": (
                None
                if self.central_delta_m_residual is None
                else float(self.central_delta_m_residual / 2.0)
            ),
            "experimental_delta_m_gev": float(self.experimental_delta_m),
            "sm_delta_m_gev": None if self.sm_delta_m is None else float(self.sm_delta_m),
            "central_delta_m_residual_gev": (
                None
                if self.central_delta_m_residual is None
                else float(self.central_delta_m_residual)
            ),
            "budget_experimental_delta_m_sigma_gev": (
                None
                if self.experimental_sigma_delta_m is None
                else float(self.experimental_sigma_delta_m)
            ),
            "budget_sm_delta_m_sigma_gev": (
                None if self.sm_sigma_delta_m is None else float(self.sm_sigma_delta_m)
            ),
            "budget_combined_delta_m_sigma_gev": (
                None
                if self.combined_sigma_delta_m is None
                else float(self.combined_sigma_delta_m)
            ),
            "legacy_full_delta_m_m12_budget": float(
                self.legacy_full_delta_m_m12_budget
            ),
            "budget_policy_note": self.note,
        }


@dataclass(frozen=True)
class DeltaF2ObservableSummary:
    """Compact per-system exclusion summary.

    In the live hadronic path, ``effective_amplitude`` is the evaluated
    observable used for pass/fail: ``abs(epsilon_K^NP)`` for kaons and
    ``abs(M12^NP)`` for the B/D systems. The dominant-operator weighted
    surrogate is retained only for the legacy ``use_hadronic=False`` fallback.
    """

    input: DeltaF2Input
    wilsons: DeltaF2WilsonCoefficients
    effective_amplitude: float
    coherent_amplitude: float
    ratio_to_bound: float
    passes: bool
    dominant_operator: str
    dominant_operator_size: float
    weighted_operator_sizes: Mapping[str, float]
    bound_override: float | None = None
    budget_policy_id: str | None = None
    confidence_level: str | None = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    @property
    def key(self) -> str:
        return self.input.key

    @property
    def system(self) -> str:
        return self.input.system

    @property
    def display_name(self) -> str:
        return self.input.display_name

    @property
    def column_name(self) -> str:
        return self.input.column_name

    @property
    def reject_reason(self) -> str:
        return self.input.reject_reason

    @property
    def bound(self) -> float:
        if self.bound_override is not None:
            return float(self.bound_override)
        return self.input.bound


@dataclass(frozen=True)
class DeltaF2ConstraintSummary:
    """Overall ``Delta F = 2`` evaluation for one fitted MFV point."""

    model_label: str
    input_bundle_label: str
    M_KK: float
    xi_KK: float
    observables: tuple[DeltaF2ObservableSummary, ...]
    lambda_ir_gev: float | None = None
    m_kk_physical_gev: float | None = None
    mass_convention_id: str | None = None
    coupling_policy_id: str = COUPLING_POLICY_PERTURBATIVE_4D_LEGACY
    g_s_4d: float | None = None
    g_eff: float | None = None
    g_s_multiplier: float | None = None

    @property
    def passes_all(self) -> bool:
        return all(item.passes for item in self.observables)

    @property
    def worst_ratio(self) -> float:
        if not self.observables:
            return 0.0
        return float(max(item.ratio_to_bound for item in self.observables))

    @property
    def max_ratio_to_bound(self) -> float:
        return self.worst_ratio

    @property
    def operator_convention(self) -> str:
        return self.model_label

    @property
    def input_bundle(self) -> str:
        return self.input_bundle_label

    @property
    def matching_scale(self) -> float:
        return self.M_KK

    @property
    def failing_reasons(self) -> tuple[str, ...]:
        return tuple(item.reject_reason for item in self.observables if not item.passes)

    @property
    def reject_reasons(self) -> tuple[str, ...]:
        return self.failing_reasons

    @property
    def by_system(self) -> dict[str, DeltaF2ObservableSummary]:
        aliases = {
            "epsilon_k": "K",
            "b_d": "B_d",
            "b_s": "B_s",
            "d": "D",
        }
        return {aliases.get(item.key, item.display_name): item for item in self.observables}

    def as_ratio_dict(self) -> dict[str, float]:
        return {item.column_name: float(item.ratio_to_bound) for item in self.observables}

    def get(self, key: str) -> DeltaF2ObservableSummary:
        for item in self.observables:
            if item.key == key or item.column_name == key or item.reject_reason == key:
                return item
        raise KeyError(key)


DEFAULT_DELTA_F2_INPUTS_V1: tuple[DeltaF2Input, ...] = (
    DeltaF2Input(
        key="epsilon_k",
        display_name="epsilon_K",
        column_name="epsilon_k_ratio",
        reject_reason="epsilon_k",
        sector="down",
        generations=(0, 1),
        bound=2.0e-8,
        ll_weight=1.0,
        rr_weight=1.0,
        lr1_weight=7.0,
        lr2_weight=2.0,
        note="Kaon mixing uses the strongest LR enhancement in the v1 input bundle.",
    ),
    DeltaF2Input(
        key="b_d",
        display_name="B_d mixing",
        column_name="b_d_ratio",
        reject_reason="b_d_mix",
        sector="down",
        generations=(0, 2),
        bound=DELTA_M_BD_BUDGET,
        ll_weight=1.0,
        rr_weight=1.0,
        lr1_weight=7.0,
        lr2_weight=2.0,
        note=(
            "B_d mixing with proper hadronic matrix elements (FLAG 2024 / PDG); "
            "budget follows the B001 one-sigma SM-vs-experiment room."
        ),
    ),
    DeltaF2Input(
        key="b_s",
        display_name="B_s mixing",
        column_name="b_s_ratio",
        reject_reason="b_s_mix",
        sector="down",
        generations=(1, 2),
        bound=DELTA_M_BS_BUDGET,
        ll_weight=1.0,
        rr_weight=1.0,
        lr1_weight=7.0,
        lr2_weight=2.0,
        note=(
            "B_s mixing with proper hadronic matrix elements (FLAG 2024 / PDG); "
            "budget follows the B003 one-sigma SM-vs-experiment room."
        ),
    ),
    DeltaF2Input(
        key="d",
        display_name="D mixing",
        column_name="d_ratio",
        reject_reason="d_mix",
        sector="up",
        generations=(0, 1),
        bound=DELTA_M_D_BUDGET,
        ll_weight=1.0,
        rr_weight=1.0,
        lr1_weight=7.0,
        lr2_weight=2.0,
        note=(
            "D0 mixing with proper hadronic matrix elements (FLAG 2024 / HFLAV); "
            "budget is a conservative long-distance envelope, not a Gaussian CL."
        ),
    ),
)


# Legacy operator-weight bounds (for backward compatibility with use_hadronic=False)
LEGACY_OPERATOR_WEIGHT_BOUNDS: dict[str, float] = {
    "epsilon_k": 2.0e-8,
    "b_d": 4.0e-7,
    "b_s": 5.5e-6,
    "d": 8.5e-9,
}


def default_delta_f2_inputs() -> tuple[DeltaF2Input, ...]:
    """Return the repo-owned v1 input bundle."""
    return DEFAULT_DELTA_F2_INPUTS_V1


def _coerce_couplings(
    source: QuarkFitResult | QuarkMassBasisCouplings,
    *,
    M_KK: float | None,
    xi_KK: float | None,
) -> QuarkMassBasisCouplings:
    if isinstance(source, QuarkMassBasisCouplings):
        return source
    required = ("M_KK", "left_down", "right_down", "left_up", "right_up")
    if all(hasattr(source, name) for name in required):
        return source  # type: ignore[return-value]
    # perturbative g_s (legacy repo_v1 behavior); callers that want g_s*
    # enhancement should pre-compute couplings with g_s_star=<value>.
    return compute_quark_kk_gluon_couplings(source, M_KK=M_KK, xi_KK=xi_KK, g_s_star=None)


def _pair_couplings(
    couplings: QuarkMassBasisCouplings,
    item: DeltaF2Input,
) -> tuple[complex, complex]:
    i, j = item.generations
    if item.sector == "down":
        left = complex(couplings.left_down[i, j])
        right = complex(couplings.right_down[i, j])
    else:
        left = complex(couplings.left_up[i, j])
        right = complex(couplings.right_up[i, j])
    return left, right


def compute_delta_f2_wilsons(
    source: QuarkFitResult | QuarkMassBasisCouplings,
    *,
    M_KK: float | None = None,
    xi_KK: float | None = None,
    inputs: Sequence[DeltaF2Input] | None = None,
) -> tuple[DeltaF2WilsonCoefficients, ...]:
    """Compute tree-level KK-gluon-inspired ``Delta F = 2`` Wilsons.

    The returned coefficients follow a compact v1 convention with one vector
    left-left operator, one vector right-right operator, and two left-right
    operators. The coefficients are matched directly at ``mu = M_KK``.
    """
    if xi_KK is not None and xi_KK <= 0.0:
        raise ValueError("xi_KK must be positive")
    couplings = _coerce_couplings(source, M_KK=M_KK, xi_KK=xi_KK)
    items = default_delta_f2_inputs() if inputs is None else tuple(inputs)
    out: list[DeltaF2WilsonCoefficients] = []
    prefactor = 1.0 / (couplings.M_KK**2)
    for item in items:
        left, right = _pair_couplings(couplings, item)
        out.append(
            DeltaF2WilsonCoefficients(
                input=item,
                M_KK=float(couplings.M_KK),
                matching_scale=float(couplings.M_KK),
                left_coupling=left,
                right_coupling=right,
                c1_vll=left * left * prefactor / 6.0,
                c1_vrr=right * right * prefactor / 6.0,
                c4_lr=-(left * right) * prefactor,
                c5_lr=(left * right) * prefactor / 3.0,
                lambda_ir_gev=getattr(couplings, "lambda_ir_gev", None),
                m_kk_physical_gev=getattr(couplings, "m_kk_physical_gev", couplings.M_KK),
                mass_convention_id=getattr(couplings, "mass_convention_id", None),
                g_s_4d=getattr(couplings, "g_s_4d", None),
                g_eff=getattr(couplings, "g_eff", getattr(couplings, "g_s", None)),
                g_s_multiplier=getattr(couplings, "g_s_multiplier", None),
                coupling_policy_id=getattr(
                    couplings,
                    "coupling_policy_id",
                    COUPLING_POLICY_PERTURBATIVE_4D_LEGACY,
                ),
                operator_convention_id=getattr(
                    couplings,
                    "operator_convention_id",
                    DELTA_F2_MODEL_V1,
                ),
            )
        )
    return tuple(out)


@dataclass(frozen=True)
class _HadronicEvaluation:
    ratio_to_bound: float
    effective_amplitude: float
    coherent_amplitude: float
    operator_sizes: Mapping[str, float]
    dominant_operator: str
    dominant_size: float
    bound: float | None = None
    budget_policy_id: str | None = None
    confidence_level: str | None = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


def _hadronic_eval_for_system(
    key: str,
    wilsons: DeltaF2WilsonCoefficients,
    *,
    epsilon_k_np_budget_override: float | None = None,
) -> _HadronicEvaluation | None:
    """Attempt proper hadronic evaluation for a known meson system.

    Returns an internal hadronic-evaluation payload or None if the system key
    is not recognized for hadronic evaluation.

    ``epsilon_k_np_budget_override`` is forwarded to :func:`evaluate_epsilon_k`
    only when ``key == "epsilon_k"``; it is silently ignored otherwise.
    """
    if key == "epsilon_k":
        eps_result = evaluate_epsilon_k(
            wilsons,
            epsilon_k_np_budget_override=epsilon_k_np_budget_override,
        )
        # Build per-operator breakdown for epsilon_K (imaginary parts)
        import math
        me = _kaon_matrix_elements()
        prefactor = abs(KAPPA_EPSILON / (math.sqrt(2.0) * DELTA_M_K))
        operator_sizes = {
            "C1_VLL": float(prefactor * abs(wilsons.c1_vll.imag * me["O1_VLL"])),
            "C1_VRR": float(prefactor * abs(wilsons.c1_vrr.imag * me["O1_VLL"])),
            "C4_LR": float(prefactor * abs(wilsons.c4_lr.imag * me["O4_LR"])),
            "C5_LR": float(prefactor * abs(wilsons.c5_lr.imag * me["O5_LR"])),
        }
        dominant_operator = max(operator_sizes, key=operator_sizes.get)
        dominant_size = float(operator_sizes[dominant_operator])
        effective_amplitude = eps_result.epsilon_k_np
        diagnostics = {
            "epsilon_k_np_signed": float(eps_result.epsilon_k_np_signed),
            "epsilon_k_np_abs": float(eps_result.epsilon_k_np),
            "epsilon_k_selected_budget_direction": (
                eps_result.selected_budget_direction
            ),
            "epsilon_k_selected_signed_budget": float(eps_result.epsilon_k_np_budget),
            "central_diagnostic_budget": float(eps_result.central_diagnostic_budget),
            "epsilon_k_np_is_absolute": True,
            "hadronic_scale_gev": float(wilsons.matching_scale),
            "running_order": eps_result.running_order,
            "running_bias_note": eps_result.running_bias_note,
            "m_kk_physical_gev": float(
                wilsons.m_kk_physical_gev
                if wilsons.m_kk_physical_gev is not None
                else wilsons.M_KK
            ),
            "lambda_ir_gev": (
                None if wilsons.lambda_ir_gev is None else float(wilsons.lambda_ir_gev)
            ),
            "mass_convention_id": wilsons.mass_convention_id,
            "coupling_policy_id": wilsons.coupling_policy_id,
            "operator_convention_id": wilsons.operator_convention_id,
            "g_s_4d": None if wilsons.g_s_4d is None else float(wilsons.g_s_4d),
            "g_eff": None if wilsons.g_eff is None else float(wilsons.g_eff),
            "g_s_multiplier": (
                None if wilsons.g_s_multiplier is None else float(wilsons.g_s_multiplier)
            ),
        }
        diagnostics.update(delta_f2_epsilon_k_budget_policy().as_diagnostics())
        return _HadronicEvaluation(
            ratio_to_bound=float(eps_result.ratio_to_budget),
            effective_amplitude=effective_amplitude,
            coherent_amplitude=float(eps_result.epsilon_k_np_signed),
            operator_sizes=operator_sizes,
            dominant_operator=dominant_operator,
            dominant_size=dominant_size,
            bound=float(eps_result.epsilon_k_np_budget),
            budget_policy_id=eps_result.budget_policy_id,
            confidence_level=eps_result.confidence_level,
            diagnostics=diagnostics,
        )
    elif key == "b_d":
        result = evaluate_bd_mixing(wilsons)
        me_vll, me_lr4, me_lr5 = _meson_matrix_elements(
            F_BD, M_BD, M_B_QUARK, M_D_QUARK_BD, B_1_BD, B_4_BD, B_5_BD
        )
        operator_sizes = {
            "C1_VLL": float(abs(wilsons.c1_vll * me_vll)),
            "C1_VRR": float(abs(wilsons.c1_vrr * me_vll)),
            "C4_LR": float(abs(wilsons.c4_lr * me_lr4)),
            "C5_LR": float(abs(wilsons.c5_lr * me_lr5)),
        }
        dominant_operator = max(operator_sizes, key=operator_sizes.get)
        dominant_size = float(operator_sizes[dominant_operator])
        diagnostics = result.diagnostics
        diagnostics = {
            **diagnostics,
            "abs_m12_np_gev": float(result.abs_m12_np),
            "hadronic_scale_gev": float(wilsons.matching_scale),
            "running_order": DELTA_F2_RUNNING_ORDER,
        }
        return _HadronicEvaluation(
            ratio_to_bound=result.ratio_to_budget,
            effective_amplitude=result.abs_m12_np,
            coherent_amplitude=result.abs_m12_np,
            operator_sizes=operator_sizes,
            dominant_operator=dominant_operator,
            dominant_size=dominant_size,
            bound=float(result.budget),
            budget_policy_id=result.budget_policy_id,
            confidence_level=result.confidence_level,
            diagnostics=diagnostics,
        )
    elif key == "b_s":
        result = evaluate_bs_mixing(wilsons)
        me_vll, me_lr4, me_lr5 = _meson_matrix_elements(
            F_BS, M_BS, M_B_QUARK, M_S_QUARK_BS, B_1_BS, B_4_BS, B_5_BS
        )
        operator_sizes = {
            "C1_VLL": float(abs(wilsons.c1_vll * me_vll)),
            "C1_VRR": float(abs(wilsons.c1_vrr * me_vll)),
            "C4_LR": float(abs(wilsons.c4_lr * me_lr4)),
            "C5_LR": float(abs(wilsons.c5_lr * me_lr5)),
        }
        dominant_operator = max(operator_sizes, key=operator_sizes.get)
        dominant_size = float(operator_sizes[dominant_operator])
        diagnostics = result.diagnostics
        diagnostics = {
            **diagnostics,
            "abs_m12_np_gev": float(result.abs_m12_np),
            "hadronic_scale_gev": float(wilsons.matching_scale),
            "running_order": DELTA_F2_RUNNING_ORDER,
        }
        return _HadronicEvaluation(
            ratio_to_bound=result.ratio_to_budget,
            effective_amplitude=result.abs_m12_np,
            coherent_amplitude=result.abs_m12_np,
            operator_sizes=operator_sizes,
            dominant_operator=dominant_operator,
            dominant_size=dominant_size,
            bound=float(result.budget),
            budget_policy_id=result.budget_policy_id,
            confidence_level=result.confidence_level,
            diagnostics=diagnostics,
        )
    elif key == "d":
        result = evaluate_d0_mixing(wilsons)
        me_vll, me_lr4, me_lr5 = _meson_matrix_elements(
            F_D, M_D0, M_C_QUARK, M_U_QUARK, B_1_D, B_4_D, B_5_D
        )
        operator_sizes = {
            "C1_VLL": float(abs(wilsons.c1_vll * me_vll)),
            "C1_VRR": float(abs(wilsons.c1_vrr * me_vll)),
            "C4_LR": float(abs(wilsons.c4_lr * me_lr4)),
            "C5_LR": float(abs(wilsons.c5_lr * me_lr5)),
        }
        dominant_operator = max(operator_sizes, key=operator_sizes.get)
        dominant_size = float(operator_sizes[dominant_operator])
        diagnostics = result.diagnostics
        diagnostics = {
            **diagnostics,
            "abs_m12_np_gev": float(result.abs_m12_np),
            "hadronic_scale_gev": float(wilsons.matching_scale),
            "running_order": DELTA_F2_RUNNING_ORDER,
        }
        return _HadronicEvaluation(
            ratio_to_bound=result.ratio_to_budget,
            effective_amplitude=result.abs_m12_np,
            coherent_amplitude=result.abs_m12_np,
            operator_sizes=operator_sizes,
            dominant_operator=dominant_operator,
            dominant_size=dominant_size,
            bound=float(result.budget),
            budget_policy_id=result.budget_policy_id,
            confidence_level=result.confidence_level,
            diagnostics=diagnostics,
        )
    return None


def evaluate_delta_f2_constraints(
    source: QuarkFitResult | QuarkMassBasisCouplings,
    *,
    M_KK: float | None = None,
    xi_KK: float | None = None,
    inputs: Sequence[DeltaF2Input] | None = None,
    apply_qcd_running: bool = True,
    mu_had: float | Mapping[str, float] | None = None,
    use_hadronic: bool = True,
    epsilon_k_np_budget_override: float | None = None,
) -> DeltaF2ConstraintSummary:
    """Evaluate the repo-owned ``Delta F = 2`` benchmark bundle.

    By default, Wilson coefficients are QCD-evolved from the matching scale
    (M_KK) down to a per-system hadronic scale using leading-log RG running
    before applying the exclusion bound: 3 GeV for kaons, 4.18 GeV for B_d/B_s,
    and 3 GeV for D0. This matches the dominant bag-parameter conventions more
    closely than the legacy global 2 GeV endpoint.

    Set ``apply_qcd_running=False`` to recover the previous behavior of using
    Wilson coefficients at the matching scale without running (backward
    compatible).

    When ``use_hadronic=True`` (the default), all meson systems are evaluated
    using proper hadronic matrix elements (BMU basis). Set
    ``use_hadronic=False`` to revert to the old operator-weight surrogate
    for B_d, B_s, and D (while epsilon_K always uses hadronic evaluation).

    Parameters
    ----------
    source : QuarkFitResult or QuarkMassBasisCouplings
        Fitted result or explicit mass-basis couplings.
    M_KK : float, optional
        KK mass scale in GeV.
    xi_KK : float
        KK convention factor.
    inputs : sequence of DeltaF2Input, optional
        Override the default input bundle.
    apply_qcd_running : bool
        If True (default), evolve Wilson coefficients from M_KK to mu_had.
    mu_had : float, mapping, or None
        Hadronic scale for RG evolution in GeV. ``None`` (default) selects the
        per-system endpoints; a scalar preserves a legacy common endpoint; a
        mapping may override selected input keys.
    use_hadronic : bool
        If True (default), use proper hadronic matrix elements for all systems.
    epsilon_k_np_budget_override : float, optional
        Override the default ε_K direction-aware one-sigma policy with an
        explicit scalar value. Used by the ``--epsilon-k-budget`` CLI sweep
        across central / low / high edges for the band-quote sensitivity study
        (see cleanup unit C02a-code and ``docs/phase_logs/phase2_h5_signoff.md``).
        Has no effect on B_d, B_s, or D systems.
    """
    couplings = _coerce_couplings(source, M_KK=M_KK, xi_KK=xi_KK)
    coefficients = compute_delta_f2_wilsons(
        couplings,
        M_KK=couplings.M_KK,
        xi_KK=couplings.xi_KK,
        inputs=inputs,
    )
    observables: list[DeltaF2ObservableSummary] = []
    for coeffs in coefficients:
        # Optionally evolve Wilson coefficients to the hadronic scale
        if apply_qcd_running:
            system_mu_had = _resolve_mu_had_for_system(coeffs.input.key, mu_had)
            evolved_coeffs = _evolve_wilsons(coeffs, mu_had=system_mu_had)
        else:
            evolved_coeffs = coeffs
        item = evolved_coeffs.input

        # Try proper hadronic evaluation
        hadronic_result = None
        if use_hadronic:
            hadronic_result = _hadronic_eval_for_system(
                item.key,
                evolved_coeffs,
                epsilon_k_np_budget_override=epsilon_k_np_budget_override,
            )

        if hadronic_result is not None:
            ratio_to_bound = hadronic_result.ratio_to_bound
            effective_amplitude = hadronic_result.effective_amplitude
            coherent_amplitude = hadronic_result.coherent_amplitude
            operator_sizes = hadronic_result.operator_sizes
            dominant_operator = hadronic_result.dominant_operator
            dominant_size = hadronic_result.dominant_size
            budget_policy_id = hadronic_result.budget_policy_id
            confidence_level = hadronic_result.confidence_level
            diagnostics = hadronic_result.diagnostics
            bound_override = hadronic_result.bound
        else:
            # Fallback: old operator-weight surrogate with legacy bounds
            weighted = {
                "C1_VLL": item.ll_weight * evolved_coeffs.c1_vll,
                "C1_VRR": item.rr_weight * evolved_coeffs.c1_vrr,
                "C4_LR": item.lr1_weight * evolved_coeffs.c4_lr,
                "C5_LR": item.lr2_weight * evolved_coeffs.c5_lr,
            }
            operator_sizes = {
                name: float(item.reference_scale**2 * abs(value))
                for name, value in weighted.items()
            }
            dominant_operator = max(operator_sizes, key=operator_sizes.get)
            coherent_amplitude = float(
                item.reference_scale**2 * abs(sum(weighted.values()))
            )
            effective_amplitude = float(operator_sizes[dominant_operator])
            dominant_size = effective_amplitude
            legacy_bound = LEGACY_OPERATOR_WEIGHT_BOUNDS.get(item.key, item.bound)
            ratio_to_bound = float(effective_amplitude / legacy_bound)
            budget_policy_id = None
            confidence_level = None
            diagnostics = {}
            bound_override = None

        observables.append(
            DeltaF2ObservableSummary(
                input=item,
                wilsons=evolved_coeffs,
                effective_amplitude=effective_amplitude,
                coherent_amplitude=coherent_amplitude,
                ratio_to_bound=ratio_to_bound,
                passes=ratio_to_bound <= 1.0,
                dominant_operator=dominant_operator,
                dominant_operator_size=dominant_size,
                weighted_operator_sizes=operator_sizes,
                bound_override=bound_override,
                budget_policy_id=budget_policy_id,
                confidence_level=confidence_level,
                diagnostics=diagnostics,
            )
        )
    return DeltaF2ConstraintSummary(
        model_label=str(getattr(couplings, "operator_convention_id", DELTA_F2_MODEL_V1)),
        input_bundle_label=DELTA_F2_INPUT_BUNDLE_V1,
        M_KK=float(couplings.M_KK),
        xi_KK=float(getattr(couplings, "xi_KK", xi_KK)),
        observables=tuple(observables),
        lambda_ir_gev=getattr(couplings, "lambda_ir_gev", None),
        m_kk_physical_gev=getattr(couplings, "m_kk_physical_gev", couplings.M_KK),
        mass_convention_id=getattr(couplings, "mass_convention_id", None),
        coupling_policy_id=getattr(
            couplings,
            "coupling_policy_id",
            COUPLING_POLICY_PERTURBATIVE_4D_LEGACY,
        ),
        g_s_4d=getattr(couplings, "g_s_4d", None),
        g_eff=getattr(couplings, "g_eff", getattr(couplings, "g_s", None)),
        g_s_multiplier=getattr(couplings, "g_s_multiplier", None),
    )


DeltaF2WilsonSet = DeltaF2WilsonCoefficients


def _evolve_wilsons(
    wilsons: DeltaF2WilsonCoefficients,
    mu_had: float = 2.0,
) -> DeltaF2WilsonCoefficients:
    """Return a new DeltaF2WilsonCoefficients with QCD-evolved coefficients.

    Uses leading-log QCD RG evolution from the matching scale (M_KK) down to
    mu_had.  VLL/VRR use the BMU current-current anomalous dimension, while
    the scalar LR coefficients C4_LR/C5_LR are evolved in the conventional
    scalar basis used by the B4/B5 matrix elements:

        O4_LR = (bar h^alpha P_L q^alpha)(bar h^beta P_R q^beta)
        O5_LR = (bar h^alpha P_L q^beta)(bar h^beta P_R q^alpha)

    With BMU's vector LR operator this means Q1_LR^BMU = -2 O5_LR and
    Q2_LR^BMU = O4_LR.  The FLAG-style B5 input is not sign-flipped.

    The ``qcd_running.evolve_deltaf2_wilsons`` symbol is imported at module
    scope (lifted in C03 cleanup, R04-I2, 2026-05-25).  ``qcd_running`` does
    not import ``deltaf2`` at module scope, so the lift is non-circular.

    Running-order caveat: Delta-F=2 running is LO; the LR (C4) enhancement is
    likely low by ~10-15% in M_KK vs NLO -- eps_K floors are therefore mildly
    CONSERVATIVE (understated) in this direction.
    """
    c_vll_low, c_vrr_low, c4_lr_low, c5_lr_low = evolve_deltaf2_wilsons(
        wilsons.c1_vll,
        wilsons.c1_vrr,
        wilsons.c4_lr,
        wilsons.c5_lr,
        mu_high=wilsons.matching_scale,
        mu_low=mu_had,
    )
    return DeltaF2WilsonCoefficients(
        input=wilsons.input,
        M_KK=wilsons.M_KK,
        matching_scale=mu_had,
        left_coupling=wilsons.left_coupling,
        right_coupling=wilsons.right_coupling,
        c1_vll=c_vll_low,
        c1_vrr=c_vrr_low,
        c4_lr=c4_lr_low,
        c5_lr=c5_lr_low,
        lambda_ir_gev=wilsons.lambda_ir_gev,
        m_kk_physical_gev=wilsons.m_kk_physical_gev,
        mass_convention_id=wilsons.mass_convention_id,
        g_s_4d=wilsons.g_s_4d,
        g_eff=wilsons.g_eff,
        g_s_multiplier=wilsons.g_s_multiplier,
        coupling_policy_id=wilsons.coupling_policy_id,
        operator_convention_id=wilsons.operator_convention_id,
    )


# ---------------------------------------------------------------------------
# Kaon hadronic parameters for proper epsilon_K and Delta m_K evaluation
# ---------------------------------------------------------------------------

F_K = 0.1557              # GeV, kaon decay constant (PDG)
M_K = 0.49761             # GeV, kaon mass (PDG)
DELTA_M_K = 3.484e-15     # GeV, K_L - K_S mass difference (PDG)
M_S_2GEV = 0.0934         # GeV, strange quark MS-bar mass at 2 GeV (FLAG)
M_D_2GEV = 0.00467        # GeV, down quark MS-bar mass at 2 GeV (FLAG)
B_1_K = 0.5503            # B_K MS-bar(2 GeV), FLAG 2024

# FLAG B4/B5 scale note:
# The LR kaon bag inputs below are FLAG 2024 MS-bar(3 GeV) values, so the
# default kaon Wilson endpoint is now 3 GeV.  The VLL B_K and quark-mass inputs
# retain the pre-existing 2 GeV literals; a full common-scheme kaon input
# refresh is a separate hadronic-input update.
B_4_K_3GEV = 0.903        # B_4^K MS-bar(3 GeV), FLAG 2024
B_5_K_3GEV = 0.691        # B_5^K MS-bar(3 GeV), FLAG 2024
B_4_K = B_4_K_3GEV        # compatibility alias; prefer B_4_K_3GEV in new code
B_5_K = B_5_K_3GEV        # compatibility alias; prefer B_5_K_3GEV in new code
KAPPA_EPSILON = 0.94       # multiplicative correction (Buras et al.)
EPSILON_K_EXP = 2.228e-3   # experimental value (PDG)
EPSILON_K_SM = 2.161e-3    # SM epsilon_K, Brod-Gorbahn-Stamou 2020
EPSILON_K_BUDGET_POLICY_ID = "epsilon_k_bgs2020_pdg2024_bgs_exp_one_sigma_v1"
EPSILON_K_BUDGET_CONFIDENCE_LEVEL = "68.27% one_sigma_sensitivity"
EPSILON_K_BUDGET_DOC_CITATION = "docs/audits/epsilon_k_sm_decision.md:34-37,71-100"
EPSILON_K_SIGMA_BGS = 0.18e-3
EPSILON_K_SIGMA_EXP = 0.011e-3
EPSILON_K_SM_CHOICE_SENSITIVITY = 0.15e-3
EPSILON_K_CENTRAL_RESIDUAL_SIGNED = EPSILON_K_EXP - EPSILON_K_SM
EPSILON_K_CENTRAL_BUDGET = abs(EPSILON_K_CENTRAL_RESIDUAL_SIGNED)
EPSILON_K_PRIMARY_COMBINED_SIGMA = math.sqrt(
    EPSILON_K_SIGMA_BGS**2 + EPSILON_K_SIGMA_EXP**2
)
EPSILON_K_SIGNED_LOWER_EDGE = (
    EPSILON_K_CENTRAL_RESIDUAL_SIGNED - EPSILON_K_PRIMARY_COMBINED_SIGMA
)
EPSILON_K_SIGNED_UPPER_EDGE = (
    EPSILON_K_CENTRAL_RESIDUAL_SIGNED + EPSILON_K_PRIMARY_COMBINED_SIGMA
)
EPSILON_K_BUDGET_LOWERS = abs(EPSILON_K_SIGNED_LOWER_EDGE)
EPSILON_K_BUDGET_RAISES = abs(EPSILON_K_SIGNED_UPPER_EDGE)
EPSILON_K_TIGHT_BUDGET = min(EPSILON_K_BUDGET_LOWERS, EPSILON_K_BUDGET_RAISES)
EPSILON_K_LOOSE_BUDGET = max(EPSILON_K_BUDGET_LOWERS, EPSILON_K_BUDGET_RAISES)

EPSILON_K_BUDGET_POLICY = EpsilonKBudgetPolicy(
    policy_id=EPSILON_K_BUDGET_POLICY_ID,
    confidence_level=EPSILON_K_BUDGET_CONFIDENCE_LEVEL,
    doc_citation=EPSILON_K_BUDGET_DOC_CITATION,
    experimental_value=EPSILON_K_EXP,
    sm_value=EPSILON_K_SM,
    central_residual_signed=EPSILON_K_CENTRAL_RESIDUAL_SIGNED,
    central_budget=EPSILON_K_CENTRAL_BUDGET,
    sigma_bgs=EPSILON_K_SIGMA_BGS,
    sigma_exp=EPSILON_K_SIGMA_EXP,
    primary_combined_sigma=EPSILON_K_PRIMARY_COMBINED_SIGMA,
    sm_choice_sensitivity=EPSILON_K_SM_CHOICE_SENSITIVITY,
    lower_signed_edge=EPSILON_K_SIGNED_LOWER_EDGE,
    upper_signed_edge=EPSILON_K_SIGNED_UPPER_EDGE,
    budget_lowers_epsilon_k=EPSILON_K_BUDGET_LOWERS,
    budget_raises_epsilon_k=EPSILON_K_BUDGET_RAISES,
    tight_budget=EPSILON_K_TIGHT_BUDGET,
    loose_budget=EPSILON_K_LOOSE_BUDGET,
)


def delta_f2_epsilon_k_budget_policy() -> EpsilonKBudgetPolicy:
    """Return the canonical shared epsilon_K budget policy."""

    return EPSILON_K_BUDGET_POLICY


B_D_MIXING_BUDGET_POLICY_ID = "b_d_delta_m_hpqcd2019_hflav2025_one_sigma_v1"
B_S_MIXING_BUDGET_POLICY_ID = "b_s_delta_m_flag2024_hflav2024_one_sigma_v1"
D0_MIXING_BUDGET_POLICY_ID = "d0_delta_m_exp_half_long_distance_envelope_v1"
DELTA_M_BUDGET_CONFIDENCE_LEVEL = "68.27% one_sigma_sensitivity"
D0_MIXING_BUDGET_CONFIDENCE_LEVEL = (
    "not_a_gaussian_cl_conservative_long_distance_envelope"
)

B_D_MIXING_BUDGET_POLICY = DeltaMNPBudgetPolicy(
    policy_id=B_D_MIXING_BUDGET_POLICY_ID,
    confidence_level=DELTA_M_BUDGET_CONFIDENCE_LEVEL,
    doc_citation=(
        "flavor_catalog/processes/beauty/B001.yaml:"
        "pdg_or_equivalent.canonical_experimental_average,"
        "standard_model_prediction; auxiliary_code_inputs.deltaf2_bd_constants"
    ),
    system="B_d",
    budget=DELTA_M_BD_BUDGET,
    construction="(|Delta m_exp - Delta m_SM| + sqrt(sigma_exp^2 + sigma_SM^2)) / 2",
    experimental_delta_m=DELTA_M_BD_EXP,
    sm_delta_m=DELTA_M_BD_SM,
    central_delta_m_residual=abs(DELTA_M_BD_EXP - DELTA_M_BD_SM),
    experimental_sigma_delta_m=DELTA_M_BD_EXP_SIGMA,
    sm_sigma_delta_m=DELTA_M_BD_SM_SIGMA,
    combined_sigma_delta_m=DELTA_M_BD_COMBINED_SIGMA,
    legacy_full_delta_m_m12_budget=DELTA_M_BD_EXP / 2.0,
    note=(
        "B001 catalog policy promoted into the core: HFLAV/PDG 2025 "
        "Delta m_d=0.5069(19) ps^-1 with the in-code GeV conversion, "
        "core SM central 3.6e-13 GeV, and HPQCD 2019 sigma_SM=0.062 ps^-1."
    ),
)

B_S_MIXING_BUDGET_POLICY = DeltaMNPBudgetPolicy(
    policy_id=B_S_MIXING_BUDGET_POLICY_ID,
    confidence_level=DELTA_M_BUDGET_CONFIDENCE_LEVEL,
    doc_citation=(
        "flavor_catalog/processes/beauty/B003.yaml:"
        "pdg_or_equivalent.canonical_hflav_recommended;"
        "auxiliary_theory_inputs.flag_2024_bmixing"
    ),
    system="B_s",
    budget=DELTA_M_BS_BUDGET,
    construction="(|Delta m_exp - Delta m_SM| + sqrt(sigma_exp^2 + sigma_SM^2)) / 2",
    experimental_delta_m=DELTA_M_BS_EXP,
    sm_delta_m=DELTA_M_BS_SM,
    central_delta_m_residual=abs(DELTA_M_BS_EXP - DELTA_M_BS_SM),
    experimental_sigma_delta_m=DELTA_M_BS_EXP_SIGMA,
    sm_sigma_delta_m=DELTA_M_BS_SM_SIGMA,
    combined_sigma_delta_m=DELTA_M_BS_COMBINED_SIGMA,
    legacy_full_delta_m_m12_budget=DELTA_M_BS_EXP / 2.0,
    note=(
        "B003 catalog policy promoted into the core: HFLAV Fall 2024 "
        "Delta m_s=17.766(6) ps^-1, hbar=6.582119569e-13 GeV ps, "
        "core SM central 1.17e-11 GeV, and FLAG 2024 "
        "f_Bs sqrt(Bhat_Bs)=256.1(5.7) MeV propagated as a 2x relative "
        "Delta m_s theory error."
    ),
)

D0_MIXING_BUDGET_POLICY = DeltaMNPBudgetPolicy(
    policy_id=D0_MIXING_BUDGET_POLICY_ID,
    confidence_level=D0_MIXING_BUDGET_CONFIDENCE_LEVEL,
    doc_citation="quarkConstraints/deltaf2.py:DELTA_M_D_EXP",
    system="D0",
    budget=DELTA_M_D_BUDGET,
    construction="Delta m_D_exp / 2",
    experimental_delta_m=DELTA_M_D_EXP,
    sm_delta_m=None,
    central_delta_m_residual=None,
    experimental_sigma_delta_m=None,
    sm_sigma_delta_m=None,
    combined_sigma_delta_m=None,
    legacy_full_delta_m_m12_budget=DELTA_M_D_BUDGET,
    note=(
        "D0 mixing remains a conservative long-distance envelope because the "
        "short-distance SM subtraction is not defensible; this is not a 68% "
        "Gaussian confidence level."
    ),
)


def delta_f2_bd_budget_policy() -> DeltaMNPBudgetPolicy:
    """Return the B_d Delta-m NP-budget policy."""

    return B_D_MIXING_BUDGET_POLICY


def delta_f2_bs_budget_policy() -> DeltaMNPBudgetPolicy:
    """Return the B_s Delta-m NP-budget policy."""

    return B_S_MIXING_BUDGET_POLICY


def delta_f2_d0_budget_policy() -> DeltaMNPBudgetPolicy:
    """Return the D0 conservative Delta-m envelope policy."""

    return D0_MIXING_BUDGET_POLICY

KAON_HADRONIC_PARAMS_V1 = "kaon_hadronic_params_bmu_2gev_v1"

# ---------------------------------------------------------------------------
# B_d meson hadronic parameters (FLAG 2024 / PDG)
# ---------------------------------------------------------------------------

F_BD = 0.1900               # GeV, decay constant (FLAG 2024)
M_BD = 5.27972              # GeV, meson mass (PDG)
M_B_QUARK = 4.18            # GeV, b quark MS-bar mass at m_b
M_D_QUARK_BD = 0.00467      # GeV, d quark MS-bar mass at 2 GeV
B_1_BD = 0.87               # VLL bag parameter (FLAG 2024, renormalized)
B_4_BD = 1.02               # LR bag parameter (FLAG 2024)
B_5_BD = 0.96               # LR bag parameter (FLAG 2024)

# ---------------------------------------------------------------------------
# B_s meson hadronic parameters (FLAG 2024 / PDG)
# ---------------------------------------------------------------------------

F_BS = 0.2303               # GeV, decay constant (FLAG 2024)
M_BS = 5.36692              # GeV, meson mass (PDG)
M_S_QUARK_BS = 0.0934       # GeV, s quark MS-bar mass at 2 GeV
B_1_BS = 0.87               # VLL bag parameter (FLAG 2024)
B_4_BS = 1.02               # LR bag parameter (FLAG 2024)
B_5_BS = 0.96               # LR bag parameter (FLAG 2024)

# ---------------------------------------------------------------------------
# D0 meson hadronic parameters (FLAG 2024 / PDG / HFLAV)
# ---------------------------------------------------------------------------

F_D = 0.2120                # GeV, decay constant (FLAG 2024)
M_D0 = 1.86484              # GeV, meson mass (PDG)
M_C_QUARK = 1.27            # GeV, c quark MS-bar mass at m_c
M_U_QUARK = 0.00216         # GeV, u quark MS-bar mass at 2 GeV
B_1_D = 0.75                # VLL bag parameter (less precise)
B_4_D = 1.0                 # LR bag parameter (estimated)
B_5_D = 1.0                 # LR bag parameter (estimated)
MESON_HADRONIC_PARAMS_V1 = "meson_hadronic_params_bmu_per_system_mu_v2"


def _kaon_matrix_elements() -> dict[str, float]:
    """Compute kaon matrix elements in the code's O1/O4/O5 basis.

    Returns a dict with keys 'O1_VLL', 'O4_LR', 'O5_LR' giving the real-valued
    hadronic matrix elements <K-bar|O_i|K> for each operator.  O1_VRR has the
    same matrix element as O1_VLL by parity.  O4_LR/O5_LR are the scalar LR
    operators paired with the conventional FLAG B4/B5 bag parameters:

        O4_LR = (bar s^alpha P_L d^alpha)(bar s^beta P_R d^beta)
        O5_LR = (bar s^alpha P_L d^beta)(bar s^beta P_R d^alpha)

    The O5_LR contraction below carries the conventional positive sign, so the
    BMU vector-LR Fierz map is Q1_LR^BMU = -2 O5_LR and B_5_K_3GEV stays
    positive.
    """
    fk2_mk = F_K**2 * M_K
    m_ratio_sq = (M_K / (M_S_2GEV + M_D_2GEV)) ** 2

    # GGMS (hep-ph/9604387) Eq. (8), M12-ready single-power-of-m_M normalization
    # (B=1).  The colour-SINGLET O4 carries the LARGE chiral coefficient
    # (R/4 + 1/24); the colour-CROSSED O5 carries the SMALL (R/12 + 1/8).  The
    # previous code had these swapped AND each x2 too large (missing 1/(2 m_M));
    # see PLAN §2.2.  O1 is (1/3), half the legacy (2/3), required so the SM box
    # reproduces the textbook epsilon_K master formula (PLAN §2.2 SM-box anchor).
    o1_vll = (1.0 / 3.0) * fk2_mk * B_1_K
    o4_lr = (m_ratio_sq * (1.0 / 4.0) + 1.0 / 24.0) * fk2_mk * B_4_K_3GEV
    o5_lr = (m_ratio_sq * (1.0 / 12.0) + 1.0 / 8.0) * fk2_mk * B_5_K_3GEV

    return {
        "O1_VLL": o1_vll,
        "O4_LR": o4_lr,
        "O5_LR": o5_lr,
    }


def _compute_m12_np(wilsons: DeltaF2WilsonCoefficients) -> complex:
    """Compute M_12^NP = sum_i C_i * <K-bar|O_i|K> for the kaon system."""
    me = _kaon_matrix_elements()
    return (
        wilsons.c1_vll * me["O1_VLL"]
        + wilsons.c1_vrr * me["O1_VLL"]   # VRR has same ME as VLL by parity
        + wilsons.c4_lr * me["O4_LR"]
        + wilsons.c5_lr * me["O5_LR"]
    )


@dataclass(frozen=True)
class EpsilonKResult:
    """NP contribution to epsilon_K from Delta F = 2 operators."""

    im_m12_np: float           # Im(M_12^NP) in GeV
    epsilon_k_np: float        # |epsilon_K^NP|
    epsilon_k_np_budget: float  # selected direction-aware allowed NP budget
    ratio_to_budget: float     # |epsilon_k_np_signed| / selected budget
    passes: bool               # ratio <= 1.0
    epsilon_k_np_signed: float = 0.0
    selected_budget_direction: str = "unknown"
    central_diagnostic_budget: float = EPSILON_K_CENTRAL_BUDGET
    budget_raises_epsilon_k: float = EPSILON_K_BUDGET_RAISES
    budget_lowers_epsilon_k: float = EPSILON_K_BUDGET_LOWERS
    primary_combined_sigma: float = EPSILON_K_PRIMARY_COMBINED_SIGMA
    sm_choice_sensitivity: float = EPSILON_K_SM_CHOICE_SENSITIVITY
    budget_policy_id: str = EPSILON_K_BUDGET_POLICY_ID
    confidence_level: str = EPSILON_K_BUDGET_CONFIDENCE_LEVEL
    running_order: str = DELTA_F2_RUNNING_ORDER
    running_bias_note: str = DELTA_F2_LO_RUNNING_BIAS_NOTE


@dataclass(frozen=True)
class DeltaMKResult:
    """NP contribution to Delta m_K."""

    abs_m12_np: float       # |M_12^NP| in GeV
    ratio_to_exp: float     # |M_12^NP| / (Delta_m_K / 2)
    passes: bool            # ratio <= 1.0 (NP shouldn't exceed experimental value)


def evaluate_epsilon_k(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    epsilon_k_np_budget_override: float | None = None,
) -> EpsilonKResult:
    """Evaluate the NP contribution to epsilon_K using proper hadronic matrix elements.

    This uses the physical formula:
        epsilon_K^NP = (kappa_epsilon / (sqrt(2) * Delta_m_K)) * Im(M_12^NP)
    and compares to the signed, direction-aware one-sigma NP budget from
    :func:`delta_f2_epsilon_k_budget_policy`.

    Parameters
    ----------
    wilsons : DeltaF2WilsonCoefficients
        Wilson coefficients (typically RG-evolved to mu_had).
    epsilon_k_np_budget_override : float, optional
        Replace the default direction-aware budget with an explicit numerical
        scalar value. Used by the legacy
        ``--epsilon-k-budget`` CLI sweep across central / low / high edges
        for the band-quote sensitivity study (cleanup unit C02a-code).
    """
    m12_np = _compute_m12_np(wilsons)
    im_m12_np = float(m12_np.imag)
    epsilon_k_np_signed = float(
        KAPPA_EPSILON / (math.sqrt(2.0) * DELTA_M_K) * im_m12_np
    )
    epsilon_k_np = abs(epsilon_k_np_signed)
    if epsilon_k_np_budget_override is None:
        policy = delta_f2_epsilon_k_budget_policy()
        budget, direction = policy.selected_signed_budget(epsilon_k_np_signed)
        budget_policy_id = policy.policy_id
        confidence_level = policy.confidence_level
    else:
        if epsilon_k_np_budget_override <= 0.0:
            raise ValueError(
                "epsilon_k_np_budget_override must be positive"
            )
        budget = float(epsilon_k_np_budget_override)
        direction = "manual_scalar_override"
        budget_policy_id = f"{EPSILON_K_BUDGET_POLICY_ID}:manual_scalar_override"
        confidence_level = "manual_scalar_override"
    ratio = abs(epsilon_k_np_signed) / budget if budget > 0.0 else float("inf")

    return EpsilonKResult(
        im_m12_np=im_m12_np,
        epsilon_k_np=epsilon_k_np,
        epsilon_k_np_budget=budget,
        ratio_to_budget=ratio,
        passes=ratio <= 1.0,
        epsilon_k_np_signed=epsilon_k_np_signed,
        selected_budget_direction=direction,
        central_diagnostic_budget=EPSILON_K_CENTRAL_BUDGET,
        budget_raises_epsilon_k=EPSILON_K_BUDGET_RAISES,
        budget_lowers_epsilon_k=EPSILON_K_BUDGET_LOWERS,
        primary_combined_sigma=EPSILON_K_PRIMARY_COMBINED_SIGMA,
        sm_choice_sensitivity=EPSILON_K_SM_CHOICE_SENSITIVITY,
        budget_policy_id=budget_policy_id,
        confidence_level=confidence_level,
        running_order=DELTA_F2_RUNNING_ORDER,
        running_bias_note=DELTA_F2_LO_RUNNING_BIAS_NOTE,
    )


def evaluate_delta_mk(
    wilsons: DeltaF2WilsonCoefficients,
) -> DeltaMKResult:
    """Evaluate the NP contribution to Delta m_K.

    The constraint is |M_12^NP| < Delta_m_K / 2 (NP should not exceed the
    experimental mass difference).
    """
    m12_np = _compute_m12_np(wilsons)
    abs_m12_np = abs(m12_np)
    half_dm = DELTA_M_K / 2.0
    ratio = abs_m12_np / half_dm if half_dm > 0.0 else float("inf")

    return DeltaMKResult(
        abs_m12_np=abs_m12_np,
        ratio_to_exp=ratio,
        passes=ratio <= 1.0,
    )


def evaluate_epsilon_k_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    mu_had: float = DELTA_F2_MU_HAD_KAON_GEV,
    *,
    epsilon_k_np_budget_override: float | None = None,
) -> EpsilonKResult:
    """Like ``evaluate_epsilon_k`` but with QCD RG evolution from M_KK to mu_had.

    The Wilson coefficients are evolved from the matching scale stored in
    ``wilsons.matching_scale`` (typically M_KK) down to ``mu_had`` (default 3 GeV)
    using leading-log QCD running before evaluating epsilon_K with the hadronic
    matrix elements at that scale.

    Parameters
    ----------
    wilsons : DeltaF2WilsonCoefficients
        Wilson coefficients at the matching scale.
    mu_had : float
        Hadronic scale in GeV (default 3.0 for the kaon LR path).
    epsilon_k_np_budget_override : float, optional
        Override for the NP budget; see :func:`evaluate_epsilon_k`.

    Returns
    -------
    EpsilonKResult
        The epsilon_K evaluation with RG-evolved Wilson coefficients.
    """
    evolved = _evolve_wilsons(wilsons, mu_had=mu_had)
    return evaluate_epsilon_k(
        evolved,
        epsilon_k_np_budget_override=epsilon_k_np_budget_override,
    )


def evaluate_delta_mk_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    mu_had: float = DELTA_F2_MU_HAD_KAON_GEV,
) -> DeltaMKResult:
    """Like ``evaluate_delta_mk`` but with QCD RG evolution from M_KK to mu_had.

    The Wilson coefficients are evolved from the matching scale stored in
    ``wilsons.matching_scale`` (typically M_KK) down to ``mu_had`` (default 3 GeV)
    using leading-log QCD running before evaluating Delta m_K.

    Parameters
    ----------
    wilsons : DeltaF2WilsonCoefficients
        Wilson coefficients at the matching scale.
    mu_had : float
        Hadronic scale in GeV (default 3.0 for the kaon LR path).

    Returns
    -------
    DeltaMKResult
        The Delta m_K evaluation with RG-evolved Wilson coefficients.
    """
    evolved = _evolve_wilsons(wilsons, mu_had=mu_had)
    return evaluate_delta_mk(evolved)


# ---------------------------------------------------------------------------
# Generic meson M_12^NP computation for ANY pseudoscalar meson
# ---------------------------------------------------------------------------


def _meson_matrix_elements(
    f_P: float,
    m_P: float,
    m_q1: float,
    m_q2: float,
    B_1: float,
    B_4: float,
    B_5: float,
) -> tuple[float, float, float]:
    """Hadronic matrix elements for a generic pseudoscalar meson.

    Returns (me_vll, me_lr4, me_lr5) in GeV^3.
    Same structure as ``_kaon_matrix_elements`` but with meson-specific inputs.
    """
    fp2_mp = f_P**2 * m_P
    r_chi = (m_P / (m_q1 + m_q2)) ** 2
    # GGMS Eq. (8) M12-ready normalization; colour-singlet O4 gets the LARGE
    # coefficient, colour-crossed O5 the SMALL one, O1 is (1/3).  See PLAN §2.2
    # and ``_kaon_matrix_elements``; all six ME sites must stay byte-identical.
    me_vll = (1.0 / 3.0) * fp2_mp * B_1
    me_lr4 = (r_chi / 4.0 + 1.0 / 24.0) * fp2_mp * B_4
    me_lr5 = (r_chi / 12.0 + 1.0 / 8.0) * fp2_mp * B_5
    return me_vll, me_lr4, me_lr5


def compute_m12_np(
    wilsons: DeltaF2WilsonCoefficients,
    f_P: float,
    m_P: float,
    m_q1: float,
    m_q2: float,
    B_1: float,
    B_4: float,
    B_5: float,
) -> complex:
    """Compute M_12^NP from Wilson coefficients and hadronic parameters."""
    me_vll, me_lr4, me_lr5 = _meson_matrix_elements(
        f_P, m_P, m_q1, m_q2, B_1, B_4, B_5
    )
    return (
        wilsons.c1_vll * me_vll
        + wilsons.c1_vrr * me_vll
        + wilsons.c4_lr * me_lr4
        + wilsons.c5_lr * me_lr5
    )


@dataclass(frozen=True)
class MesonMixingResult:
    """NP contribution to meson mixing for B_d, B_s, or D0."""

    system: str          # "B_d", "B_s", "D0"
    abs_m12_np: float    # |M_12^NP| in GeV
    budget: float        # allowed NP budget in GeV (half of Delta m)
    ratio_to_budget: float  # abs_m12_np / budget
    passes: bool         # ratio <= 1.0
    budget_policy_id: str = ""
    confidence_level: str = ""
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


def _bd_budget() -> float:
    """NP budget for B_d mixing under the one-sigma B001 policy."""
    return delta_f2_bd_budget_policy().budget


def _bs_budget() -> float:
    """NP budget for B_s mixing under the one-sigma B003 policy."""
    return delta_f2_bs_budget_policy().budget


def _d0_budget() -> float:
    """NP budget for D0 mixing: conservative long-distance envelope."""
    return delta_f2_d0_budget_policy().budget


def evaluate_bd_mixing(
    wilsons: DeltaF2WilsonCoefficients,
) -> MesonMixingResult:
    """Evaluate the NP contribution to B_d mixing using proper hadronic matrix elements."""
    m12 = compute_m12_np(
        wilsons, F_BD, M_BD, M_B_QUARK, M_D_QUARK_BD, B_1_BD, B_4_BD, B_5_BD
    )
    abs_m12 = abs(m12)
    budget = _bd_budget()
    ratio = abs_m12 / budget if budget > 0.0 else float("inf")
    policy = delta_f2_bd_budget_policy()
    return MesonMixingResult(
        system="B_d",
        abs_m12_np=abs_m12,
        budget=budget,
        ratio_to_budget=ratio,
        passes=ratio <= 1.0,
        budget_policy_id=policy.policy_id,
        confidence_level=policy.confidence_level,
        diagnostics=policy.as_diagnostics(),
    )


def evaluate_bs_mixing(
    wilsons: DeltaF2WilsonCoefficients,
) -> MesonMixingResult:
    """Evaluate the NP contribution to B_s mixing using proper hadronic matrix elements."""
    m12 = compute_m12_np(
        wilsons, F_BS, M_BS, M_B_QUARK, M_S_QUARK_BS, B_1_BS, B_4_BS, B_5_BS
    )
    abs_m12 = abs(m12)
    budget = _bs_budget()
    ratio = abs_m12 / budget if budget > 0.0 else float("inf")
    policy = delta_f2_bs_budget_policy()
    return MesonMixingResult(
        system="B_s",
        abs_m12_np=abs_m12,
        budget=budget,
        ratio_to_budget=ratio,
        passes=ratio <= 1.0,
        budget_policy_id=policy.policy_id,
        confidence_level=policy.confidence_level,
        diagnostics=policy.as_diagnostics(),
    )


def evaluate_d0_mixing(
    wilsons: DeltaF2WilsonCoefficients,
) -> MesonMixingResult:
    """Evaluate the NP contribution to D0 mixing using proper hadronic matrix elements."""
    m12 = compute_m12_np(
        wilsons, F_D, M_D0, M_C_QUARK, M_U_QUARK, B_1_D, B_4_D, B_5_D
    )
    abs_m12 = abs(m12)
    budget = _d0_budget()
    ratio = abs_m12 / budget if budget > 0.0 else float("inf")
    policy = delta_f2_d0_budget_policy()
    return MesonMixingResult(
        system="D0",
        abs_m12_np=abs_m12,
        budget=budget,
        ratio_to_budget=ratio,
        passes=ratio <= 1.0,
        budget_policy_id=policy.policy_id,
        confidence_level=policy.confidence_level,
        diagnostics=policy.as_diagnostics(),
    )


def evaluate_bd_mixing_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    mu_had: float = DELTA_F2_MU_HAD_B_GEV,
) -> MesonMixingResult:
    """Like ``evaluate_bd_mixing`` but with QCD RG evolution from M_KK to m_b."""
    evolved = _evolve_wilsons(wilsons, mu_had=mu_had)
    return evaluate_bd_mixing(evolved)


def evaluate_bs_mixing_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    mu_had: float = DELTA_F2_MU_HAD_B_GEV,
) -> MesonMixingResult:
    """Like ``evaluate_bs_mixing`` but with QCD RG evolution from M_KK to m_b."""
    evolved = _evolve_wilsons(wilsons, mu_had=mu_had)
    return evaluate_bs_mixing(evolved)


def evaluate_d0_mixing_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    mu_had: float = DELTA_F2_MU_HAD_D_GEV,
) -> MesonMixingResult:
    """Like ``evaluate_d0_mixing`` but with QCD RG evolution from M_KK to 3 GeV."""
    evolved = _evolve_wilsons(wilsons, mu_had=mu_had)
    return evaluate_d0_mixing(evolved)
