"""Baryonic ``b -> s mu+ mu-`` machinery for ``Lambda_b -> Lambda mu mu``.

This module is append-only relative to :mod:`quarkConstraints.rare_b_dilepton`:
it reuses that module's CKM inputs and documented C9/C10 RS Wilson proxy, then
adds a lightweight baryonic ``Lambda_b -> Lambda`` form-factor normalization.

Physics convention
------------------
The short-distance normalization follows the shared Buras ``b -> s l l``
Hamiltonian used for B005/B016.  The SM and NP rates keep only C9/C10 and use
the leading-HQET ``F_+`` and ``F_-`` dipole form factors from
Detmold-Lin-Meinel-Wingate, Phys. Rev. D 87 (2013) 074502,
arXiv:1212.4827.  The integrated rate is intended for the high-q2 bin where
those form factors are directly constrained.

NEEDS-HUMAN-PHYSICS
-------------------
The RS contribution inherits the rare_b_dilepton C9/C10 proxy.  A rigorous
baryonic likelihood needs the complete electroweak KK/Z/Z' and lepton-sector
matching, C7/nonlocal-charm terms, tensor/scalar operators, the full
Lambda_b -> Lambda ten-form-factor basis and covariance, and experimental bin
correlations.  The v1 implementation compresses chirality into the same
C9/C10 totals as B016 and treats the baryonic form-factor normalization as a
documented proxy uncertainty.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping

from .couplings import QuarkMassBasisCouplings
from .rare_b_dilepton import (
    RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_FRACTION,
    RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    RareBDileptonSMInputs,
    RareBDileptonWilsonCoefficients,
    ckm_factors,
    compute_rare_b_dilepton_wilsons,
)

RARE_B_BARYONIC_DILEPTON_MODEL_V1 = (
    "rare_b_dilepton_lambdab_lambda_mumu_c9_c10_hqet_ff_v1"
)
RARE_B_BARYONIC_DILEPTON_INPUT_BUNDLE_V1 = (
    "rare_b_dilepton_lambdab_lambda_mumu_detmold2013_repo_ckm_v1"
)
RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_BUNDLE_V1 = (
    "lambdab_lambda_fplus_fminus_detmold2013_hqet_dipole_v1"
)
RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_CITATION = (
    "Detmold, Lin, Meinel, Wingate, Phys. Rev. D 87 (2013) 074502, "
    "arXiv:1212.4827; leading-HQET F_+ and F_- dipole parameters "
    "N_+=3.188 GeV^2, X_+=1.852 GeV, N_-=4.124 GeV^2, X_-=1.634 GeV"
)
RARE_B_BARYONIC_DILEPTON_LIMITATION_V1 = (
    "NEEDS-HUMAN-PHYSICS: C9/C10-only Lambda_b -> Lambda mu mu evaluator. "
    "C7, nonlocal charm, scalar/tensor operators, the full ten baryonic form "
    "factors, lattice covariance, and q2-bin experimental correlations are "
    "not included in v1."
)
RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_FRACTION = (
    RARE_B_DILEPTON_EXCLUSIVE_BK_PROXY_THEORY_UNCERTAINTY_FRACTION
)
RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_RATIONALE = (
    "B021 reuses the B016 30% C9/C10-proxy theory envelope.  For the baryonic "
    "mode this also covers the leading-HQET Lambda_b -> Lambda form-factor "
    "normalization, missing covariance, omitted C7/nonlocal-charm terms, and "
    "the NEEDS-HUMAN-PHYSICS RS matching limitation."
)


@dataclass(frozen=True)
class RareBBaryonFormFactorInputs:
    """Masses, lifetime, and leading-HQET ``Lambda_b -> Lambda`` form factors."""

    mode_key: str = "lambdab_lambda"
    display_name: str = "Lambda_b0 -> Lambda mu+ mu-"
    parent_mass_gev: float = 5.6195
    daughter_mass_gev: float = 1.115683
    lifetime_ps: float = 1.425
    n_plus_gev2: float = 3.188
    x_plus_gev: float = 1.852
    n_minus_gev2: float = 4.124
    x_minus_gev: float = 1.634
    n_plus_sigma_gev2: float = 0.268
    x_plus_sigma_gev: float = 0.074
    n_minus_sigma_gev2: float = 0.750
    x_minus_sigma_gev: float = 0.144
    cov_n_plus_x_plus_gev3: float = 0.0198
    cov_n_minus_x_minus_gev3: float = 0.106
    hqet_c_gamma: float = 1.0
    hqet_c_v: float = 0.0
    source: str = RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_CITATION

    def __post_init__(self) -> None:
        for name in (
            "parent_mass_gev",
            "daughter_mass_gev",
            "lifetime_ps",
            "n_plus_gev2",
            "x_plus_gev",
            "n_minus_gev2",
            "x_minus_gev",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        for name in (
            "n_plus_sigma_gev2",
            "x_plus_sigma_gev",
            "n_minus_sigma_gev2",
            "x_minus_sigma_gev",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be non-negative and finite")
        for name in ("cov_n_plus_x_plus_gev3", "cov_n_minus_x_minus_gev3"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        for name in ("hqet_c_gamma", "hqet_c_v"):
            value = float(getattr(self, name))
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        if self.parent_mass_gev <= self.daughter_mass_gev:
            raise ValueError("parent_mass_gev must exceed daughter_mass_gev")

    @property
    def q2_max_gev2(self) -> float:
        """Return the kinematic endpoint ``(m_Lambda_b - m_Lambda)^2``."""

        return float((self.parent_mass_gev - self.daughter_mass_gev) ** 2)


@dataclass(frozen=True)
class RareBBaryonDileptonInputs:
    """Reusable inputs for baryonic ``Lambda_b -> Lambda mu+ mu-`` bins."""

    input_bundle: str = RARE_B_BARYONIC_DILEPTON_INPUT_BUNDLE_V1
    c9_sm: float = 4.27
    integration_steps: int = 800
    short_distance_inputs: RareBDileptonSMInputs = field(
        default_factory=RareBDileptonSMInputs
    )
    form_factor: RareBBaryonFormFactorInputs = field(
        default_factory=RareBBaryonFormFactorInputs
    )
    constants_citation: str = (
        "Buras b->sll Hamiltonian convention reused from rare_b_dilepton; "
        "Detmold-Lin-Meinel-Wingate leading-HQET Lambda_b -> Lambda F_+/- "
        "normalization for the high-q2 bin"
    )

    def __post_init__(self) -> None:
        if not math.isfinite(float(self.c9_sm)) or float(self.c9_sm) == 0.0:
            raise ValueError("c9_sm must be finite and non-zero")
        if int(self.integration_steps) < 20:
            raise ValueError("integration_steps must be at least 20")


@dataclass(frozen=True)
class RareBBaryonDileptonBranchingResult:
    """Partial branching-fraction prediction for ``Lambda_b -> Lambda mu mu``."""

    model_label: str
    input_bundle: str
    mode_key: str
    q2_min_gev2: float
    q2_max_gev2: float
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    ratio_to_sm: float
    average_differential_branching_fraction: float
    sm_average_differential_branching_fraction: float
    c9_total: complex
    c10_total: complex
    c9_vector_np: complex
    c10_axial_np: complex
    lambda_t: complex
    form_factor_bundle: str
    wilsons: RareBDileptonWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str | bool] = field(default_factory=dict)


def default_lambdab_to_lambda_dilepton_inputs() -> RareBBaryonDileptonInputs:
    """Return the repo-owned baryonic ``Lambda_b -> Lambda mu mu`` input bundle."""

    return RareBBaryonDileptonInputs()


def _positive_float(value: object, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _tau_ps_to_gev_inverse(tau_ps: float, hbar_gev_s: float) -> float:
    return float(tau_ps * 1.0e-12 / hbar_gev_s)


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


def _q2_bounds(
    form_factor: RareBBaryonFormFactorInputs,
    inputs: RareBBaryonDileptonInputs,
    q2_min_gev2: float | None,
    q2_max_gev2: float | None,
) -> tuple[float, float]:
    threshold = 4.0 * inputs.short_distance_inputs.muon_mass_gev**2
    lower = threshold if q2_min_gev2 is None else float(q2_min_gev2)
    upper = form_factor.q2_max_gev2 if q2_max_gev2 is None else float(q2_max_gev2)
    if lower < threshold - 1.0e-12:
        raise ValueError("q2_min_gev2 is below the dimuon threshold")
    if upper > form_factor.q2_max_gev2 + 1.0e-12:
        raise ValueError("q2_max_gev2 is above the Lambda_b -> Lambda endpoint")
    if upper <= lower:
        raise ValueError("q2_max_gev2 must be larger than q2_min_gev2")
    return float(max(lower, threshold)), float(min(upper, form_factor.q2_max_gev2))


def lambdab_to_lambda_fplus_fminus(
    q2_gev2: float,
    form_factor: RareBBaryonFormFactorInputs | None = None,
) -> tuple[float, float]:
    """Evaluate the Detmold et al. leading-HQET ``F_+`` and ``F_-`` form factors."""

    mode = RareBBaryonFormFactorInputs() if form_factor is None else form_factor
    q2 = float(q2_gev2)
    if q2 < 0.0 or q2 > mode.q2_max_gev2 + 1.0e-12:
        raise ValueError("q2_gev2 is outside the Lambda_b -> Lambda physical range")
    e_lambda = (
        mode.parent_mass_gev**2 + mode.daughter_mass_gev**2 - q2
    ) / (2.0 * mode.parent_mass_gev)
    delta_e = e_lambda - mode.daughter_mass_gev
    f_plus = mode.n_plus_gev2 / (mode.x_plus_gev + delta_e) ** 2
    f_minus = mode.n_minus_gev2 / (mode.x_minus_gev + delta_e) ** 2
    if not math.isfinite(f_plus) or f_plus <= 0.0:
        raise ValueError("F_+ must be positive and finite")
    if not math.isfinite(f_minus) or f_minus <= 0.0:
        raise ValueError("F_- must be positive and finite")
    return float(f_plus), float(f_minus)


def _baryonic_differential_branching_fraction(
    q2_gev2: float,
    *,
    inputs: RareBBaryonDileptonInputs,
    c9_total: complex,
    c10_total: complex,
    lambda_t: complex,
) -> float:
    sd = inputs.short_distance_inputs
    mode = inputs.form_factor
    q2 = float(q2_gev2)
    threshold = 4.0 * sd.muon_mass_gev**2
    if q2 <= threshold or q2 >= mode.q2_max_gev2:
        return 0.0

    m_parent = float(mode.parent_mass_gev)
    m_daughter = float(mode.daughter_mass_gev)
    phase_left = (m_parent - m_daughter) ** 2 - q2
    phase_right = (m_parent + m_daughter) ** 2 - q2
    beta2 = 1.0 - threshold / q2
    if phase_left <= 0.0 or phase_right <= 0.0 or beta2 <= 0.0:
        return 0.0

    f_plus, f_minus = lambdab_to_lambda_fplus_fminus(q2, mode)
    f_combined = phase_left * f_minus**2 + phase_right * f_plus**2
    g_combined = (
        m_parent**6
        - m_parent**4 * (3.0 * m_daughter**2 + q2)
        - m_parent**2
        * (q2 - m_daughter**2)
        * (3.0 * m_daughter**2 + q2)
        + (q2 - m_daughter**2) ** 3
    )
    c_gamma = float(mode.hqet_c_gamma)
    c_v = float(mode.hqet_c_v)
    muon_mass = float(sd.muon_mass_gev)

    a_10_10 = (
        (
            (2.0 * c_gamma**2 + 2.0 * c_gamma * c_v + c_v**2)
            * (2.0 * muon_mass**2 + q2)
            * (
                m_parent**4
                - 2.0 * m_parent**2 * m_daughter**2
                + (q2 - m_daughter**2) ** 2
            )
            + 2.0
            * m_parent**2
            * q2
            * (
                4.0 * c_gamma**2 * (q2 - 4.0 * muon_mass**2)
                - (2.0 * c_gamma * c_v + c_v**2)
                * (q2 - 10.0 * muon_mass**2)
            )
        )
        * f_combined
        + 4.0
        * c_gamma
        * (c_gamma + c_v)
        * (2.0 * muon_mass**2 + q2)
        * g_combined
        * f_plus
        * f_minus
    )
    a_9_9 = (
        (
            (2.0 * c_gamma**2 + 2.0 * c_gamma * c_v + c_v**2)
            * (m_parent**4 + (q2 - m_daughter**2) ** 2)
            - 2.0
            * m_parent**2
            * (
                2.0 * c_gamma**2 * (m_daughter**2 - 2.0 * q2)
                + (2.0 * c_gamma * c_v + c_v**2) * (m_daughter**2 + q2)
            )
        )
        * f_combined
        + 4.0 * c_gamma * (c_gamma + c_v) * g_combined * f_plus * f_minus
    )

    prefactor = (
        sd.alpha_em_mz**2
        * sd.gf_gev_minus2**2
        * abs(complex(lambda_t)) ** 2
        / (6144.0 * math.pi**5 * q2**2 * m_parent**5)
    )
    phase = math.sqrt(beta2) * math.sqrt(phase_left * phase_right)
    rate = prefactor * phase * (
        q2 * abs(complex(c10_total)) ** 2 * a_10_10
        + q2
        * (q2 + 2.0 * muon_mass**2)
        * abs(complex(c9_total)) ** 2
        * a_9_9
    )
    branching_density = (
        _tau_ps_to_gev_inverse(mode.lifetime_ps, sd.hbar_gev_s) * rate
    )
    if not math.isfinite(branching_density) or branching_density < 0.0:
        raise ValueError("Lambda_b -> Lambda mu mu dBR/dq2 must be finite")
    return float(branching_density)


def _integrated_branching_fraction(
    *,
    inputs: RareBBaryonDileptonInputs,
    q2_min_gev2: float,
    q2_max_gev2: float,
    c9_total: complex,
    c10_total: complex,
    lambda_t: complex,
) -> float:
    return _simpson_integral(
        lambda q2: _baryonic_differential_branching_fraction(
            q2,
            inputs=inputs,
            c9_total=c9_total,
            c10_total=c10_total,
            lambda_t=lambda_t,
        ),
        q2_min_gev2,
        q2_max_gev2,
        intervals=inputs.integration_steps,
    )


def sm_lambdab_to_lambda_mumu_branching_fraction(
    *,
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    inputs: RareBBaryonDileptonInputs | None = None,
) -> RareBBaryonDileptonBranchingResult:
    """Evaluate the SM-limit partial ``Lambda_b -> Lambda mu+ mu-`` BR."""

    return evaluate_lambdab_to_lambda_mumu(
        None,
        q2_min_gev2=q2_min_gev2,
        q2_max_gev2=q2_max_gev2,
        inputs=inputs,
    )


def evaluate_lambdab_to_lambda_mumu(
    source: QuarkMassBasisCouplings | RareBDileptonWilsonCoefficients | None = None,
    *,
    q2_min_gev2: float | None = None,
    q2_max_gev2: float | None = None,
    m_kk_gev: float | None = None,
    inputs: RareBBaryonDileptonInputs | None = None,
) -> RareBBaryonDileptonBranchingResult:
    """Evaluate ``BR(Lambda_b -> Lambda mu mu)`` with the C9/C10 proxy.

    NEEDS-HUMAN-PHYSICS: RS new physics is inherited from the documented
    ``b -> s mu mu`` C9/C10 proxy.  The baryonic rate is a high-q2 leading-HQET
    form-factor integral and does not include the complete baryonic
    form-factor/covariance or nonlocal-charm treatment.
    """

    p = default_lambdab_to_lambda_dilepton_inputs() if inputs is None else inputs
    q2_min, q2_max = _q2_bounds(p.form_factor, p, q2_min_gev2, q2_max_gev2)
    factors = ckm_factors("b_s", p.short_distance_inputs)

    wilsons: RareBDileptonWilsonCoefficients | None
    if source is None:
        wilsons = None
        c9_np = c9p_np = c10_np = c10p_np = 0.0j
    elif isinstance(source, RareBDileptonWilsonCoefficients):
        wilsons = source
        if wilsons.transition_key != "b_s":
            raise ValueError("Lambda_b -> Lambda mu mu requires b_s Wilsons")
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

    br = _integrated_branching_fraction(
        inputs=p,
        q2_min_gev2=q2_min,
        q2_max_gev2=q2_max,
        c9_total=c9_total,
        c10_total=c10_total,
        lambda_t=factors.lambda_t,
    )
    sm_br = _integrated_branching_fraction(
        inputs=p,
        q2_min_gev2=q2_min,
        q2_max_gev2=q2_max,
        c9_total=c9_sm_total,
        c10_total=c10_sm_total,
        lambda_t=factors.lambda_t,
    )
    q2_mid = 0.5 * (q2_min + q2_max)
    fplus_mid, fminus_mid = lambdab_to_lambda_fplus_fminus(q2_mid, p.form_factor)
    fplus_min, fminus_min = lambdab_to_lambda_fplus_fminus(q2_min, p.form_factor)
    fplus_max, fminus_max = lambdab_to_lambda_fplus_fminus(
        q2_max - min(1.0e-6, 0.5 * (q2_max - q2_min)),
        p.form_factor,
    )
    bin_width = q2_max - q2_min

    diagnostics: dict[str, float | complex | str | bool] = {
        "mode_display_name": p.form_factor.display_name,
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
        "q2_bin_width_gev2": float(bin_width),
        "differential_branching_fraction_at_bin_center": (
            _baryonic_differential_branching_fraction(
                q2_mid,
                inputs=p,
                c9_total=c9_total,
                c10_total=c10_total,
                lambda_t=factors.lambda_t,
            )
        ),
        "average_differential_branching_fraction": float(br / bin_width),
        "sm_average_differential_branching_fraction": float(sm_br / bin_width),
        "form_factor_bundle": RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_BUNDLE_V1,
        "form_factor_source": p.form_factor.source,
        "form_factor_uncertainty": RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_RATIONALE,
        "fplus_q2_min": float(fplus_min),
        "fminus_q2_min": float(fminus_min),
        "fplus_q2_mid": float(fplus_mid),
        "fminus_q2_mid": float(fminus_mid),
        "fplus_q2_max_minus_epsilon": float(fplus_max),
        "fminus_q2_max_minus_epsilon": float(fminus_max),
        "n_plus_gev2": float(p.form_factor.n_plus_gev2),
        "x_plus_gev": float(p.form_factor.x_plus_gev),
        "n_minus_gev2": float(p.form_factor.n_minus_gev2),
        "x_minus_gev": float(p.form_factor.x_minus_gev),
        "integration_steps": float(p.integration_steps),
        "parent_mass_gev": float(p.form_factor.parent_mass_gev),
        "daughter_mass_gev": float(p.form_factor.daughter_mass_gev),
        "lifetime_ps": float(p.form_factor.lifetime_ps),
        "hqet_c_gamma": float(p.form_factor.hqet_c_gamma),
        "hqet_c_v": float(p.form_factor.hqet_c_v),
        "constants_citation": p.constants_citation,
        "matching_assumption": RARE_B_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
        "baryonic_limitations": RARE_B_BARYONIC_DILEPTON_LIMITATION_V1,
        "c7_nonlocal_charm_omitted": True,
        "chirality_compressed_like_b_to_k": True,
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

    return RareBBaryonDileptonBranchingResult(
        model_label=RARE_B_BARYONIC_DILEPTON_MODEL_V1,
        input_bundle=p.input_bundle,
        mode_key=p.form_factor.mode_key,
        q2_min_gev2=float(q2_min),
        q2_max_gev2=float(q2_max),
        branching_fraction=float(br),
        sm_branching_fraction=float(sm_br),
        np_shift_branching_fraction=float(br - sm_br),
        ratio_to_sm=float(br / sm_br),
        average_differential_branching_fraction=float(br / bin_width),
        sm_average_differential_branching_fraction=float(sm_br / bin_width),
        c9_total=complex(c9_total),
        c10_total=complex(c10_total),
        c9_vector_np=complex(c9_vector_np),
        c10_axial_np=complex(c10_axial_np),
        lambda_t=complex(factors.lambda_t),
        form_factor_bundle=RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_BUNDLE_V1,
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


__all__ = [
    "RARE_B_BARYONIC_DILEPTON_MODEL_V1",
    "RARE_B_BARYONIC_DILEPTON_INPUT_BUNDLE_V1",
    "RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_BUNDLE_V1",
    "RARE_B_BARYONIC_DILEPTON_FORM_FACTOR_CITATION",
    "RARE_B_BARYONIC_DILEPTON_LIMITATION_V1",
    "RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_FRACTION",
    "RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_RATIONALE",
    "RareBBaryonFormFactorInputs",
    "RareBBaryonDileptonInputs",
    "RareBBaryonDileptonBranchingResult",
    "default_lambdab_to_lambda_dilepton_inputs",
    "lambdab_to_lambda_fplus_fminus",
    "sm_lambdab_to_lambda_mumu_branching_fraction",
    "evaluate_lambdab_to_lambda_mumu",
]
