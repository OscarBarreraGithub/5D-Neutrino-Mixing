"""Lightweight ``B -> K*(892) l+ l-`` C9/C10 proxy machinery.

This module is append-only relative to :mod:`quarkConstraints.rare_b_dilepton`:
it reuses the shared ``b -> s l l`` CKM factors and documented C9/C10 RS
Wilson proxy, then adds only the vector-meson phase-space and form-factor
normalization needed for catalog LFU-ratio constraints such as B019.

Physics convention
------------------
The short-distance convention and Wilson matching are inherited from
``rare_b_dilepton``.  The exclusive K* rate is a massless-lepton, C9/C10-only
helicity-weight proxy,

    dBR/dq2 ~ tau_B G_F^2 alpha^2 |V_tb V_ts^*|^2 sqrt(lambda)
              beta_l(q2) H_K*(q2) (|C9|^2 + |C10|^2),

where ``H_K*(q2)`` is built from simple V, A1, A2 pole form factors for the
``B0 -> K*(892)0`` vector mode.  This is sufficient for LFU stress tests where
the SM electron denominator uses the same normalization and the dominant
observable is the C9/C10 response ratio.

NEEDS-HUMAN-PHYSICS: this is not a production ``B -> K* l l`` angular or
global-fit backend.  It omits C7, nonlocal charm, scalar/tensor operators,
right-handed transversity structure, lepton-mass/QED differences, S-wave
effects, form-factor covariance, and experimental bin covariance.  A complete
RS prediction needs electroweak KK/Z/Z' and lepton-sector inputs not carried by
``ParameterPoint``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping

from .couplings import QuarkMassBasisCouplings
from .rare_b_dilepton import (
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareBDileptonSMInputs,
    RareBDileptonWilsonCoefficients,
    ckm_factors,
    compute_rare_b_dilepton_wilsons,
)

RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_MODEL_V1 = (
    "rare_b_dilepton_b_to_kstar_mumu_c9_c10_form_factor_proxy_v1"
)
RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_INPUT_BUNDLE_V1 = (
    "rare_b_dilepton_b_to_kstar_mumu_vector_proxy_repo_ckm_v1"
)
RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_BUNDLE_V1 = (
    "b_to_kstar_va1a2_pole_proxy_v1"
)
RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_CITATION = (
    "Lightweight B0 -> K*(892)0 V/A1/A2 single-pole normalization inspired by "
    "standard B -> K* l+l- helicity form-factor decompositions; intended only "
    "as an LFU-ratio response proxy until full form-factor covariance is added."
)
RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: the V/A1/A2 pole parameters are a lightweight "
    "K* normalization proxy for LFU-ratio response tests, not a sourced "
    "precision form-factor/covariance bundle."
)
RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_LIMITATION_V1 = (
    "NEEDS-HUMAN-PHYSICS: C9/C10-only B -> K* mu mu proxy.  It omits C7, "
    "nonlocal charm, scalar/tensor operators, full right-handed transversity "
    "structure, lepton-mass/QED effects, S-wave pollution, form-factor "
    "covariance, and bin-covariance inputs."
)


@dataclass(frozen=True)
class RareBToKStarFormFactorInputs:
    """Masses, lifetime, and simple V/A1/A2 inputs for one ``B -> K*`` mode."""

    mode_key: str
    display_name: str
    parent_mass_gev: float
    daughter_mass_gev: float
    lifetime_ps: float
    v_0: float = 0.36
    a1_0: float = 0.27
    a2_0: float = 0.26
    v_pole_mass_gev: float = 5.42
    a1_pole_mass_gev: float = 5.83
    a2_pole_mass_gev: float = 5.83
    v_slope: float = 0.0
    a1_slope: float = 0.0
    a2_slope: float = 0.0
    source: str = RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_CITATION

    def __post_init__(self) -> None:
        for name in (
            "parent_mass_gev",
            "daughter_mass_gev",
            "lifetime_ps",
            "v_0",
            "a1_0",
            "a2_0",
            "v_pole_mass_gev",
            "a1_pole_mass_gev",
            "a2_pole_mass_gev",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        for name in ("v_slope", "a1_slope", "a2_slope"):
            if not math.isfinite(float(getattr(self, name))):
                raise ValueError(f"{name} must be finite")
        if self.parent_mass_gev <= self.daughter_mass_gev:
            raise ValueError("parent_mass_gev must be larger than daughter_mass_gev")

    @property
    def q2_max_gev2(self) -> float:
        """Return the kinematic endpoint ``(m_B - m_K*)^2``."""

        return float((self.parent_mass_gev - self.daughter_mass_gev) ** 2)


@dataclass(frozen=True)
class RareBToKStarDileptonInputs:
    """Reusable inputs for exclusive ``B -> K* mu+ mu-`` integration."""

    input_bundle: str = RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_INPUT_BUNDLE_V1
    c9_sm: float = 4.27
    integration_steps: int = 800
    short_distance_inputs: RareBDileptonSMInputs = field(
        default_factory=RareBDileptonSMInputs
    )
    bzero_kstarzero: RareBToKStarFormFactorInputs = field(
        default_factory=lambda: RareBToKStarFormFactorInputs(
            mode_key="bzero_kstarzero",
            display_name="B0 -> K*(892)0 mu+ mu-",
            parent_mass_gev=5.27966,
            daughter_mass_gev=0.89555,
            lifetime_ps=1.520,
            source=(
                "PDG-era B0/K*(892)0 masses and tau_B0; "
                + RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_CITATION
            ),
        )
    )
    constants_citation: str = (
        "Buras b->sll Hamiltonian convention inherited from rare_b_dilepton; "
        "massless-limit B -> K* ll helicity-weight proxy with simple V/A1/A2 "
        "form-factor normalization; repo CKM target "
        "quarkConstraints.modern.inputs.ModernDefaultCKMTarget"
    )

    def __post_init__(self) -> None:
        if not math.isfinite(float(self.c9_sm)) or float(self.c9_sm) == 0.0:
            raise ValueError("c9_sm must be finite and non-zero")
        if int(self.integration_steps) < 20:
            raise ValueError("integration_steps must be at least 20")

    def mode(self, mode_key: str) -> RareBToKStarFormFactorInputs:
        """Return the form-factor input block for ``mode_key``."""

        if mode_key == "bzero_kstarzero":
            return self.bzero_kstarzero
        raise ValueError(f"unsupported B -> K* mode {mode_key!r}")


@dataclass(frozen=True)
class RareBToKStarDileptonBranchingResult:
    """Partial branching-fraction prediction for ``B -> K* mu+ mu-``."""

    model_label: str
    input_bundle: str
    mode_key: str
    q2_min_gev2: float
    q2_max_gev2: float
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    ratio_to_sm: float
    c9_total: complex
    c10_total: complex
    c9_vector_np: complex
    c10_axial_np: complex
    lambda_t: complex
    form_factor_bundle: str
    wilsons: RareBDileptonWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str | bool] = field(default_factory=dict)


def default_b_to_kstar_dilepton_inputs() -> RareBToKStarDileptonInputs:
    """Return the repo-owned exclusive ``B -> K* mu+ mu-`` input bundle."""

    return RareBToKStarDileptonInputs()


def _positive_float(value: object, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _tau_ps_to_gev_inverse(tau_ps: float, hbar_gev_s: float) -> float:
    return float(tau_ps * 1.0e-12 / hbar_gev_s)


def _kallen(x: float, y: float, z: float) -> float:
    return float(x * x + y * y + z * z - 2.0 * (x * y + x * z + y * z))


def _simpson_integral(
    function: object,
    lower: float,
    upper: float,
    *,
    intervals: int,
) -> float:
    n = int(intervals)
    if n % 2:
        n += 1
    if n <= 0 or upper <= lower:
        raise ValueError("Simpson integration requires upper > lower and intervals > 0")
    h = (upper - lower) / n
    total = float(function(lower)) + float(function(upper))
    odd_sum = 0.0
    even_sum = 0.0
    for index in range(1, n):
        value = float(function(lower + index * h))
        if index % 2:
            odd_sum += value
        else:
            even_sum += value
    return float(h / 3.0 * (total + 4.0 * odd_sum + 2.0 * even_sum))


def _pole_form_factor(
    q2_gev2: float,
    *,
    f0: float,
    pole_mass_gev: float,
    slope: float,
    parent_mass_gev: float,
) -> float:
    q2 = float(q2_gev2)
    value = f0 * (1.0 + slope * q2 / parent_mass_gev**2) / (
        1.0 - q2 / pole_mass_gev**2
    )
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"B -> K* form factor at q2={q2!r} must be positive")
    return float(value)


def b_to_kstar_form_factors(
    q2_gev2: float,
    mode: RareBToKStarFormFactorInputs,
) -> Mapping[str, float]:
    """Return the proxy ``V, A1, A2`` form factors at ``q2``."""

    q2 = float(q2_gev2)
    if q2 < 0.0 or q2 >= (mode.parent_mass_gev + mode.daughter_mass_gev) ** 2:
        raise ValueError(f"q2_gev2={q2_gev2!r} outside B -> K* form-factor domain")
    return {
        "V": _pole_form_factor(
            q2,
            f0=mode.v_0,
            pole_mass_gev=mode.v_pole_mass_gev,
            slope=mode.v_slope,
            parent_mass_gev=mode.parent_mass_gev,
        ),
        "A1": _pole_form_factor(
            q2,
            f0=mode.a1_0,
            pole_mass_gev=mode.a1_pole_mass_gev,
            slope=mode.a1_slope,
            parent_mass_gev=mode.parent_mass_gev,
        ),
        "A2": _pole_form_factor(
            q2,
            f0=mode.a2_0,
            pole_mass_gev=mode.a2_pole_mass_gev,
            slope=mode.a2_slope,
            parent_mass_gev=mode.parent_mass_gev,
        ),
    }


def _kstar_helicity_weight(
    q2_gev2: float,
    *,
    mode: RareBToKStarFormFactorInputs,
) -> float:
    q2 = float(q2_gev2)
    if q2 <= 0.0 or q2 >= mode.q2_max_gev2:
        return 0.0
    lam = _kallen(mode.parent_mass_gev**2, mode.daughter_mass_gev**2, q2)
    if lam <= 0.0:
        return 0.0
    ff = b_to_kstar_form_factors(q2, mode)
    sqrt_lam = math.sqrt(lam)
    m_b = mode.parent_mass_gev
    m_v = mode.daughter_mass_gev
    h_perp = sqrt_lam * ff["V"] / (m_b + m_v)
    h_parallel = (m_b + m_v) * ff["A1"]
    h_long = (
        (m_b * m_b - m_v * m_v - q2) * (m_b + m_v) * ff["A1"]
        - lam * ff["A2"] / (m_b + m_v)
    ) / (2.0 * m_v * math.sqrt(q2))
    value = h_perp * h_perp + h_parallel * h_parallel + h_long * h_long
    if not math.isfinite(value):
        raise ValueError("B -> K* helicity weight must be finite")
    return float(max(0.0, value))


def _b_to_kstar_q2_bounds(
    mode: RareBToKStarFormFactorInputs,
    inputs: RareBToKStarDileptonInputs,
    q2_min_gev2: float | None,
    q2_max_gev2: float | None,
) -> tuple[float, float]:
    threshold = 4.0 * inputs.short_distance_inputs.muon_mass_gev**2
    lower = threshold if q2_min_gev2 is None else float(q2_min_gev2)
    upper = mode.q2_max_gev2 if q2_max_gev2 is None else float(q2_max_gev2)
    if lower < threshold - 1.0e-12:
        raise ValueError("q2_min_gev2 is below the dimuon threshold")
    if upper > mode.q2_max_gev2 + 1.0e-12:
        raise ValueError("q2_max_gev2 is above the B -> K* endpoint")
    if upper <= lower:
        raise ValueError("q2_max_gev2 must be larger than q2_min_gev2")
    return float(max(lower, threshold)), float(min(upper, mode.q2_max_gev2))


def _b_to_kstar_differential_branching_fraction(
    q2_gev2: float,
    *,
    mode: RareBToKStarFormFactorInputs,
    inputs: RareBToKStarDileptonInputs,
    c9_total: complex,
    c10_total: complex,
    lambda_t: complex,
) -> float:
    sd = inputs.short_distance_inputs
    q2 = float(q2_gev2)
    threshold = 4.0 * sd.muon_mass_gev**2
    if q2 <= threshold or q2 >= mode.q2_max_gev2:
        return 0.0
    lam = _kallen(mode.parent_mass_gev**2, mode.daughter_mass_gev**2, q2)
    if lam <= 0.0:
        return 0.0
    beta2 = 1.0 - threshold / q2
    if beta2 <= 0.0:
        return 0.0
    beta = math.sqrt(beta2)
    lepton_mass_factor = beta * (1.0 + 2.0 * sd.muon_mass_gev**2 / q2)
    helicity_weight = _kstar_helicity_weight(q2, mode=mode)
    tau = _tau_ps_to_gev_inverse(mode.lifetime_ps, sd.hbar_gev_s)
    prefactor = (
        tau
        * sd.gf_gev_minus2**2
        * sd.alpha_em_mz**2
        * abs(lambda_t) ** 2
        / (3072.0 * math.pi**5 * mode.parent_mass_gev**3)
    )
    amplitude_power = abs(complex(c9_total)) ** 2 + abs(complex(c10_total)) ** 2
    return float(
        prefactor
        * math.sqrt(lam)
        * lepton_mass_factor
        * helicity_weight
        * amplitude_power
    )


def _integrated_b_to_kstar_branching_fraction(
    *,
    mode: RareBToKStarFormFactorInputs,
    inputs: RareBToKStarDileptonInputs,
    q2_min_gev2: float,
    q2_max_gev2: float,
    c9_total: complex,
    c10_total: complex,
    lambda_t: complex,
) -> float:
    return _simpson_integral(
        lambda q2: _b_to_kstar_differential_branching_fraction(
            q2,
            mode=mode,
            inputs=inputs,
            c9_total=c9_total,
            c10_total=c10_total,
            lambda_t=lambda_t,
        ),
        q2_min_gev2,
        q2_max_gev2,
        intervals=inputs.integration_steps,
    )


def sm_b_to_kstar_mumu_branching_fraction(
    *,
    mode: str = "bzero_kstarzero",
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    inputs: RareBToKStarDileptonInputs | None = None,
) -> RareBToKStarDileptonBranchingResult:
    """Evaluate the SM-limit partial ``BR(B -> K* mu+ mu-)`` proxy."""

    return evaluate_b_to_kstar_mumu(
        None,
        mode=mode,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        inputs=inputs,
    )


def evaluate_b_to_kstar_mumu(
    source: QuarkMassBasisCouplings | RareBDileptonWilsonCoefficients | None = None,
    *,
    mode: str = "bzero_kstarzero",
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    m_kk_gev: float | None = None,
    inputs: RareBToKStarDileptonInputs | None = None,
) -> RareBToKStarDileptonBranchingResult:
    """Evaluate exclusive ``BR(B -> K* mu+ mu-)`` with the C9/C10 proxy."""

    p = default_b_to_kstar_dilepton_inputs() if inputs is None else inputs
    mode_inputs = p.mode(mode)
    q2_min, q2_max = _b_to_kstar_q2_bounds(mode_inputs, p, q2_min_gev2, q2_max_gev2)
    factors = ckm_factors("b_s", p.short_distance_inputs)

    wilsons: RareBDileptonWilsonCoefficients | None
    if source is None:
        wilsons = None
        c9_np = c9p_np = c10_np = c10p_np = 0.0j
    elif isinstance(source, RareBDileptonWilsonCoefficients):
        wilsons = source
        if wilsons.transition_key != "b_s":
            raise ValueError("B -> K* mu mu requires b_s Wilson coefficients")
        c9_np = wilsons.c9_np
        c9p_np = wilsons.c9p_np
        c10_np = wilsons.c10_np
        c10p_np = wilsons.c10p_np
    else:
        wilsons = compute_rare_b_dilepton_wilsons(
            source,
            transition="b_s",
            m_kk_gev=m_kk_gev,
            inputs=p.short_distance_inputs,
        )
        c9_np = wilsons.c9_np
        c9p_np = wilsons.c9p_np
        c10_np = wilsons.c10_np
        c10p_np = wilsons.c10p_np

    c9_vector_np = complex(c9_np + c9p_np)
    c10_axial_np = complex(c10_np + c10p_np)
    c9_total = complex(p.c9_sm + c9_vector_np)
    c10_total = complex(p.short_distance_inputs.c10_sm + c10_axial_np)
    c9_sm_total = complex(p.c9_sm)
    c10_sm_total = complex(p.short_distance_inputs.c10_sm)

    br = _integrated_b_to_kstar_branching_fraction(
        mode=mode_inputs,
        inputs=p,
        q2_min_gev2=q2_min,
        q2_max_gev2=q2_max,
        c9_total=c9_total,
        c10_total=c10_total,
        lambda_t=factors.lambda_t,
    )
    sm_br = _integrated_b_to_kstar_branching_fraction(
        mode=mode_inputs,
        inputs=p,
        q2_min_gev2=q2_min,
        q2_max_gev2=q2_max,
        c9_total=c9_sm_total,
        c10_total=c10_sm_total,
        lambda_t=factors.lambda_t,
    )
    if sm_br <= 0.0 or not math.isfinite(sm_br):
        raise ValueError("SM B -> K* ll proxy branching fraction must be positive")

    q2_mid = 0.5 * (q2_min + q2_max)
    ff_mid = b_to_kstar_form_factors(q2_mid, mode_inputs)
    diagnostics: dict[str, float | complex | str | bool] = {
        "mode_display_name": mode_inputs.display_name,
        "lambda_wolfenstein": float(factors.lambda_wolfenstein),
        "lambda_t": complex(factors.lambda_t),
        "c9_sm": float(p.c9_sm),
        "c10_sm": float(p.short_distance_inputs.c10_sm),
        "c9_total": complex(c9_total),
        "c10_total": complex(c10_total),
        "c9_vector_np": complex(c9_vector_np),
        "c10_axial_np": complex(c10_axial_np),
        "q2_min_gev2": float(q2_min),
        "q2_max_gev2": float(q2_max),
        "q2_bin_width_gev2": float(q2_max - q2_min),
        "differential_branching_fraction_at_bin_center": (
            _b_to_kstar_differential_branching_fraction(
                q2_mid,
                mode=mode_inputs,
                inputs=p,
                c9_total=c9_total,
                c10_total=c10_total,
                lambda_t=factors.lambda_t,
            )
        ),
        "average_differential_branching_fraction": float(br / (q2_max - q2_min)),
        "helicity_weight_at_bin_center": _kstar_helicity_weight(
            q2_mid,
            mode=mode_inputs,
        ),
        "form_factor_bundle": RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_BUNDLE_V1,
        "form_factor_source": mode_inputs.source,
        "form_factor_proxy_assumption": (
            RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_ASSUMPTION_V1
        ),
        "v_0": float(mode_inputs.v_0),
        "a1_0": float(mode_inputs.a1_0),
        "a2_0": float(mode_inputs.a2_0),
        "v_q2_mid": float(ff_mid["V"]),
        "a1_q2_mid": float(ff_mid["A1"]),
        "a2_q2_mid": float(ff_mid["A2"]),
        "integration_steps": float(p.integration_steps),
        "parent_mass_gev": float(mode_inputs.parent_mass_gev),
        "daughter_mass_gev": float(mode_inputs.daughter_mass_gev),
        "lifetime_ps": float(mode_inputs.lifetime_ps),
        "constants_citation": p.constants_citation,
        "matching_assumption": RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        "exclusive_kstar_limitations": RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_LIMITATION_V1,
        "c7_nonlocal_charm_omitted": True,
        "qed_lepton_mass_omitted": True,
    }
    if wilsons is not None:
        diagnostics.update(
            {
                "m_kk_gev": float(wilsons.M_KK),
                "matching_scale_gev": float(wilsons.matching_scale),
                "left_qb_coupling": complex(wilsons.left_qb_coupling),
                "right_qb_coupling": complex(wilsons.right_qb_coupling),
                "left_qb_overlap": complex(wilsons.left_qb_overlap),
                "right_qb_overlap": complex(wilsons.right_qb_overlap),
                "left_quark_delta": complex(wilsons.left_quark_delta),
                "right_quark_delta": complex(wilsons.right_quark_delta),
                "muon_left_delta": float(wilsons.muon_left_delta),
                "muon_right_delta": float(wilsons.muon_right_delta),
                "muon_vector_delta": float(wilsons.muon_vector_delta),
                "muon_axial_delta": float(wilsons.muon_axial_delta),
                "c9_np": complex(wilsons.c9_np),
                "c10_np": complex(wilsons.c10_np),
                "c9p_np": complex(wilsons.c9p_np),
                "c10p_np": complex(wilsons.c10p_np),
            }
        )

    return RareBToKStarDileptonBranchingResult(
        model_label=RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_MODEL_V1,
        input_bundle=p.input_bundle,
        mode_key=mode,
        q2_min_gev2=float(q2_min),
        q2_max_gev2=float(q2_max),
        branching_fraction=float(br),
        sm_branching_fraction=float(sm_br),
        np_shift_branching_fraction=float(br - sm_br),
        ratio_to_sm=float(br / sm_br),
        c9_total=complex(c9_total),
        c10_total=complex(c10_total),
        c9_vector_np=complex(c9_vector_np),
        c10_axial_np=complex(c10_axial_np),
        lambda_t=complex(factors.lambda_t),
        form_factor_bundle=RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_BUNDLE_V1,
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


__all__ = [
    "RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_MODEL_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_INPUT_BUNDLE_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_BUNDLE_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_CITATION",
    "RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_FORM_FACTOR_ASSUMPTION_V1",
    "RARE_B_DILEPTON_EXCLUSIVE_BKSTAR_LIMITATION_V1",
    "RareBToKStarFormFactorInputs",
    "RareBToKStarDileptonInputs",
    "RareBToKStarDileptonBranchingResult",
    "default_b_to_kstar_dilepton_inputs",
    "b_to_kstar_form_factors",
    "sm_b_to_kstar_mumu_branching_fraction",
    "evaluate_b_to_kstar_mumu",
]
