"""Tree-level KK-gluon-inspired ``Delta F = 2`` benchmark observables.

This module intentionally implements a compact, repo-owned v1 matching layer:

- exact fitted quark mass-basis couplings come from ``quarkConstraints.couplings``
- KK-gluon exchange is represented by a small fixed operator basis
- hadronic/bound inputs live in a deterministic in-repo bundle

The goal is to turn the current MFV fit machinery into a first exclusion slice,
not to provide a general-purpose EFT/RG package.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

from .couplings import QuarkMassBasisCouplings, compute_quark_kk_gluon_couplings
from .fit import QuarkFitResult
from .qcd_running import evolve_deltaf2_wilsons
from .scales import DEFAULT_QUARK_TARGET_SCALE_GEV, DEFAULT_QUARK_XI_KK

DELTA_F2_MODEL_V1 = "kk_gluon_tree_v1"
DELTA_F2_INPUT_BUNDLE_V1 = "deltaf2_inputs_mu3tev_v1"
DELTA_F2_OPERATOR_CONVENTION = DELTA_F2_MODEL_V1
DELTA_F2_INPUT_BUNDLE = DELTA_F2_INPUT_BUNDLE_V1


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

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {
            "C1_VLL": self.c1_vll,
            "C1_VRR": self.c1_vrr,
            "C4_LR": self.c4_lr,
            "C5_LR": self.c5_lr,
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
        return self.input.bound


@dataclass(frozen=True)
class DeltaF2ConstraintSummary:
    """Overall ``Delta F = 2`` evaluation for one fitted MFV point."""

    model_label: str
    input_bundle_label: str
    M_KK: float
    xi_KK: float
    observables: tuple[DeltaF2ObservableSummary, ...]

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
        bound=1.667e-13,
        ll_weight=1.0,
        rr_weight=1.0,
        lr1_weight=7.0,
        lr2_weight=2.0,
        note="B_d mixing with proper hadronic matrix elements (FLAG 2024 / PDG).",
    ),
    DeltaF2Input(
        key="b_s",
        display_name="B_s mixing",
        column_name="b_s_ratio",
        reject_reason="b_s_mix",
        sector="down",
        generations=(1, 2),
        bound=5.844e-12,
        ll_weight=1.0,
        rr_weight=1.0,
        lr1_weight=7.0,
        lr2_weight=2.0,
        note="B_s mixing with proper hadronic matrix elements (FLAG 2024 / PDG).",
    ),
    DeltaF2Input(
        key="d",
        display_name="D mixing",
        column_name="d_ratio",
        reject_reason="d_mix",
        sector="up",
        generations=(0, 1),
        bound=3.125e-15,
        ll_weight=1.0,
        rr_weight=1.0,
        lr1_weight=7.0,
        lr2_weight=2.0,
        note="D0 mixing with proper hadronic matrix elements (FLAG 2024 / HFLAV).",
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
    xi_KK: float,
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
    xi_KK: float = DEFAULT_QUARK_XI_KK,
    inputs: Sequence[DeltaF2Input] | None = None,
) -> tuple[DeltaF2WilsonCoefficients, ...]:
    """Compute tree-level KK-gluon-inspired ``Delta F = 2`` Wilsons.

    The returned coefficients follow a compact v1 convention with one vector
    left-left operator, one vector right-right operator, and two left-right
    operators. The coefficients are matched directly at ``mu = M_KK``.
    """
    if xi_KK <= 0.0:
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
            )
        )
    return tuple(out)


def _hadronic_eval_for_system(
    key: str,
    wilsons: DeltaF2WilsonCoefficients,
    *,
    epsilon_k_np_budget_override: float | None = None,
) -> tuple[float, float, float, dict[str, float], str, float] | None:
    """Attempt proper hadronic evaluation for a known meson system.

    Returns (ratio_to_bound, effective_amplitude, coherent_amplitude,
             operator_sizes, dominant_operator, dominant_size) or None if
    the system key is not recognized for hadronic evaluation.

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
        return (
            eps_result.ratio_to_budget,
            effective_amplitude,
            effective_amplitude,
            operator_sizes,
            dominant_operator,
            dominant_size,
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
        return (
            result.ratio_to_budget,
            result.abs_m12_np,
            result.abs_m12_np,
            operator_sizes,
            dominant_operator,
            dominant_size,
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
        return (
            result.ratio_to_budget,
            result.abs_m12_np,
            result.abs_m12_np,
            operator_sizes,
            dominant_operator,
            dominant_size,
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
        return (
            result.ratio_to_budget,
            result.abs_m12_np,
            result.abs_m12_np,
            operator_sizes,
            dominant_operator,
            dominant_size,
        )
    return None


def evaluate_delta_f2_constraints(
    source: QuarkFitResult | QuarkMassBasisCouplings,
    *,
    M_KK: float | None = None,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
    inputs: Sequence[DeltaF2Input] | None = None,
    apply_qcd_running: bool = True,
    mu_had: float = 2.0,
    use_hadronic: bool = True,
    epsilon_k_np_budget_override: float | None = None,
) -> DeltaF2ConstraintSummary:
    """Evaluate the repo-owned ``Delta F = 2`` benchmark bundle.

    By default, Wilson coefficients are QCD-evolved from the matching scale
    (M_KK) down to the hadronic scale ``mu_had`` (default 2 GeV) using
    leading-log RG running before applying the exclusion bound. This is the
    physically correct procedure since the hadronic matrix elements are
    evaluated at mu_had.

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
    mu_had : float
        Hadronic scale for RG evolution in GeV (default 2.0).
    use_hadronic : bool
        If True (default), use proper hadronic matrix elements for all systems.
    epsilon_k_np_budget_override : float, optional
        Override the default ε_K NP budget |EPSILON_K_EXP - EPSILON_K_SM|
        (~6.7e-5) with an explicit numerical value. Used by the
        ``--epsilon-k-budget`` CLI sweep across central / low / high edges
        for the band-quote sensitivity study (see cleanup unit
        C02a-code and ``docs/phase_logs/phase2_h5_signoff.md``).
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
            evolved_coeffs = _evolve_wilsons(coeffs, mu_had=mu_had)
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
            (
                ratio_to_bound,
                effective_amplitude,
                coherent_amplitude,
                operator_sizes,
                dominant_operator,
                dominant_size,
            ) = hadronic_result
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
            )
        )
    return DeltaF2ConstraintSummary(
        model_label=DELTA_F2_MODEL_V1,
        input_bundle_label=DELTA_F2_INPUT_BUNDLE_V1,
        M_KK=float(couplings.M_KK),
        xi_KK=float(getattr(couplings, "xi_KK", xi_KK)),
        observables=tuple(observables),
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

# FLAG B4/B5 scale caveat:
# The LR kaon bag inputs below are FLAG 2024 MS-bar(3 GeV) values, while
# ``mu_had`` defaults to 2 GeV and Wilson coefficients are normally evolved to
# that scale before evaluation.  The matrix-element path intentionally keeps
# using these existing numbers for now so current evaluations remain
# numerically unchanged, but this is a known scale-matching caveat.
# TODO: RG-run/convert B_4^K and B_5^K to 2 GeV in the matching scheme used
# here, then update the matrix elements, provenance, and compatibility aliases.
B_4_K_3GEV = 0.903        # B_4^K MS-bar(3 GeV), FLAG 2024
B_5_K_3GEV = 0.691        # B_5^K MS-bar(3 GeV), FLAG 2024
B_4_K = B_4_K_3GEV        # compatibility alias; prefer B_4_K_3GEV in new code
B_5_K = B_5_K_3GEV        # compatibility alias; prefer B_5_K_3GEV in new code
KAPPA_EPSILON = 0.94       # multiplicative correction (Buras et al.)
EPSILON_K_EXP = 2.228e-3   # experimental value (PDG)
EPSILON_K_SM = 2.161e-3    # SM epsilon_K, Brod-Gorbahn-Stamou 2020

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
DELTA_M_BD_EXP = 3.334e-13  # GeV, experimental (PDG)
DELTA_M_BD_SM = 3.6e-13     # GeV, SM prediction (CKMfitter)

# ---------------------------------------------------------------------------
# B_s meson hadronic parameters (FLAG 2024 / PDG)
# ---------------------------------------------------------------------------

F_BS = 0.2303               # GeV, decay constant (FLAG 2024)
M_BS = 5.36692              # GeV, meson mass (PDG)
M_S_QUARK_BS = 0.0934       # GeV, s quark MS-bar mass at 2 GeV
B_1_BS = 0.87               # VLL bag parameter (FLAG 2024)
B_4_BS = 1.02               # LR bag parameter (FLAG 2024)
B_5_BS = 0.96               # LR bag parameter (FLAG 2024)
DELTA_M_BS_EXP = 1.1688e-11 # GeV, experimental (PDG)
DELTA_M_BS_SM = 1.17e-11    # GeV, SM prediction (CKMfitter)

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
DELTA_M_D_EXP = 6.25e-15    # GeV, experimental (HFLAV)

MESON_HADRONIC_PARAMS_V1 = "meson_hadronic_params_bmu_2gev_v1"


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

    o1_vll = (2.0 / 3.0) * fk2_mk * B_1_K
    o4_lr = (m_ratio_sq * (1.0 / 6.0) + 1.0 / 4.0) * fk2_mk * B_4_K_3GEV
    o5_lr = (m_ratio_sq * (1.0 / 2.0) + 1.0 / 12.0) * fk2_mk * B_5_K_3GEV

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
    epsilon_k_np_budget: float  # allowed NP budget |epsilon_K^exp - epsilon_K^SM|
    ratio_to_budget: float     # epsilon_k_np / budget
    passes: bool               # ratio <= 1.0


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
    and compares to the experimental budget |epsilon_K^exp - epsilon_K^SM|.

    Parameters
    ----------
    wilsons : DeltaF2WilsonCoefficients
        Wilson coefficients (typically RG-evolved to mu_had).
    epsilon_k_np_budget_override : float, optional
        Replace the default budget |EPSILON_K_EXP - EPSILON_K_SM| (~6.7e-5)
        with an explicit numerical value. Used by the
        ``--epsilon-k-budget`` CLI sweep across central / low / high edges
        for the band-quote sensitivity study (cleanup unit C02a-code).
    """
    import math

    m12_np = _compute_m12_np(wilsons)
    im_m12_np = float(m12_np.imag)
    epsilon_k_np = abs(KAPPA_EPSILON / (math.sqrt(2.0) * DELTA_M_K) * im_m12_np)
    if epsilon_k_np_budget_override is None:
        budget = abs(EPSILON_K_EXP - EPSILON_K_SM)
    else:
        if epsilon_k_np_budget_override <= 0.0:
            raise ValueError(
                "epsilon_k_np_budget_override must be positive"
            )
        budget = float(epsilon_k_np_budget_override)
    ratio = epsilon_k_np / budget if budget > 0.0 else float("inf")

    return EpsilonKResult(
        im_m12_np=im_m12_np,
        epsilon_k_np=epsilon_k_np,
        epsilon_k_np_budget=budget,
        ratio_to_budget=ratio,
        passes=ratio <= 1.0,
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
    mu_had: float = 2.0,
    *,
    epsilon_k_np_budget_override: float | None = None,
) -> EpsilonKResult:
    """Like ``evaluate_epsilon_k`` but with QCD RG evolution from M_KK to mu_had.

    The Wilson coefficients are evolved from the matching scale stored in
    ``wilsons.matching_scale`` (typically M_KK) down to ``mu_had`` (default 2 GeV)
    using leading-log QCD running before evaluating epsilon_K with the hadronic
    matrix elements at that scale.

    Parameters
    ----------
    wilsons : DeltaF2WilsonCoefficients
        Wilson coefficients at the matching scale.
    mu_had : float
        Hadronic scale in GeV (default 2.0).
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
    mu_had: float = 2.0,
) -> DeltaMKResult:
    """Like ``evaluate_delta_mk`` but with QCD RG evolution from M_KK to mu_had.

    The Wilson coefficients are evolved from the matching scale stored in
    ``wilsons.matching_scale`` (typically M_KK) down to ``mu_had`` (default 2 GeV)
    using leading-log QCD running before evaluating Delta m_K.

    Parameters
    ----------
    wilsons : DeltaF2WilsonCoefficients
        Wilson coefficients at the matching scale.
    mu_had : float
        Hadronic scale in GeV (default 2.0).

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
    me_vll = (2.0 / 3.0) * fp2_mp * B_1
    me_lr4 = (r_chi / 6.0 + 0.25) * fp2_mp * B_4
    me_lr5 = (r_chi / 2.0 + 1.0 / 12.0) * fp2_mp * B_5
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


def _bd_budget() -> float:
    """NP budget for B_d mixing: max(DeltaM_exp/2, |DeltaM_exp - DeltaM_SM|/2)."""
    return max(DELTA_M_BD_EXP / 2.0, abs(DELTA_M_BD_EXP - DELTA_M_BD_SM) / 2.0)


def _bs_budget() -> float:
    """NP budget for B_s mixing: max(DeltaM_exp/2, |DeltaM_exp - DeltaM_SM|/2)."""
    return max(DELTA_M_BS_EXP / 2.0, abs(DELTA_M_BS_EXP - DELTA_M_BS_SM) / 2.0)


def _d0_budget() -> float:
    """NP budget for D0 mixing: DeltaM_exp/2 (long-distance SM dominates)."""
    return DELTA_M_D_EXP / 2.0


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
    return MesonMixingResult(
        system="B_d",
        abs_m12_np=abs_m12,
        budget=budget,
        ratio_to_budget=ratio,
        passes=ratio <= 1.0,
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
    return MesonMixingResult(
        system="B_s",
        abs_m12_np=abs_m12,
        budget=budget,
        ratio_to_budget=ratio,
        passes=ratio <= 1.0,
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
    return MesonMixingResult(
        system="D0",
        abs_m12_np=abs_m12,
        budget=budget,
        ratio_to_budget=ratio,
        passes=ratio <= 1.0,
    )


def evaluate_bd_mixing_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    mu_had: float = 2.0,
) -> MesonMixingResult:
    """Like ``evaluate_bd_mixing`` but with QCD RG evolution from M_KK to mu_had."""
    evolved = _evolve_wilsons(wilsons, mu_had=mu_had)
    return evaluate_bd_mixing(evolved)


def evaluate_bs_mixing_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    mu_had: float = 2.0,
) -> MesonMixingResult:
    """Like ``evaluate_bs_mixing`` but with QCD RG evolution from M_KK to mu_had."""
    evolved = _evolve_wilsons(wilsons, mu_had=mu_had)
    return evaluate_bs_mixing(evolved)


def evaluate_d0_mixing_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    mu_had: float = 2.0,
) -> MesonMixingResult:
    """Like ``evaluate_d0_mixing`` but with QCD RG evolution from M_KK to mu_had."""
    evolved = _evolve_wilsons(wilsons, mu_had=mu_had)
    return evaluate_d0_mixing(evolved)
