"""K_S -> mu+mu- short-distance wrapper over rare-kaon dileptons.

This adapter is append-only relative to K006/K008/K009/K010.  It reuses the
K006 ``Y`` effective-input core with Phase-3a ``s -> d mu+mu-`` C10/C10p
Wilson shifts, but applies the K_S ell=0 CP projection rather than the K_L
real-amplitude projection:

    BR(K_S -> mu+mu-)_SD,ell=0
        = (tau_KS / tau_KL) kappa_mu [
            Im(lambda_t Y_t + Y_NP) / lambda^5
            - Im(lambda_c) P_c(Y) / lambda
        ]^2.

Equivalently, this is the imaginary part of the short-distance combination
``-lambda_c Y_c + lambda_t C10`` in the Dery-Ghosh-Grossman-Schacht
convention, up to the Buras/Isidori normalization used by K006.  The SM
short-distance value is O(1.9e-13), while the catalog's O(5e-12) SM number is
total-rate context and not an SD validation anchor.

The total measured K_S rate is long-distance dominated, so this is only the
constrained short-distance component requested for K012.

NEEDS-HUMAN-PHYSICS: a rigorous K_S dimuon prediction needs the time-dependent
K_S/K_L interference/extraction treatment, long-distance two-photon amplitude,
and non-light-Z or muon-sector effects outside Phase 3a.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping

from quarkConstraints.rs_semileptonic_wilsons import RSSemileptonicWilsonBundle

from .rare_kaon_dilepton import (
    RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION,
    RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1,
    KLongMuMuShortDistanceResult,
    QuarkMassBasisCouplings,
    RareKaonDileptonSMInputs,
    RareKaonDileptonWilsonCoefficients,
    _rs_semileptonic_coeff,
    _tag_rs_result,
    klong_mumu_short_distance_from_couplings,
    klong_mumu_short_distance_from_rs_semileptonic_wilsons,
    klong_mumu_short_distance_sm,
    rare_kaon_dilepton_ckm_factors,
    rare_kaon_dilepton_default_sm_inputs,
    rare_kaon_dilepton_kappa_mu,
    rare_kaon_dilepton_wilsons_from_couplings,
    rare_kaon_y_wilsons_from_rs_semileptonic,
)

RARE_KAON_KSHORT_MUMU_MODEL_V1 = "rare_kaon_kshort_mumu_sd_im_projection_v1"
KSHORT_MUMU_CP_PROJECTION_CITATION = (
    "Dery, Ghosh, Grossman, Schacht JHEP 07 (2021) 103, arXiv:2104.06427"
)
RARE_KAON_KSHORT_MUMU_PARAMETRIZATION_CITATION = (
    RARE_KAON_DILEPTON_PARAMETRIZATION_CITATION
    + "; K012 applies the K_S ell=0 imaginary CP projection "
    "Im[-lambda_c Y_c + lambda_t C10]^2 from "
    + KSHORT_MUMU_CP_PROJECTION_CITATION
)
RARE_KAON_KSHORT_MUMU_RS_MATCHING_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: K012 maps Phase-3a s -> d mu mu C10/C10p "
    "Wilson shifts into the K006 Y effective input and applies the K_S "
    "ell=0 imaginary CP projection. Full long-distance two-photon treatment, "
    "time-dependent extraction, and non-light-Z or muon-sector effects remain "
    "outside Phase 3a."
)
RARE_KAON_KSHORT_MUMU_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1 = (
    "No additional multiplicative QCD running is applied to the NP "
    "semileptonic axial-lepton Y/C10 shift between M_KK and the low scale: "
    "the color-singlet semileptonic current has LO QCD factor 1.0 in this "
    "restricted matching. Complete EW matching and operator mixing remain part "
    "of the NEEDS-HUMAN-PHYSICS limitation."
)
KSHORT_MUMU_LIFETIME_INPUTS_CITATION = (
    "PDG-era neutral-kaon lifetimes used only for the K006-to-K012 "
    "short-distance normalization conversion"
)


@dataclass(frozen=True)
class KShortMuMuLifetimeInputs:
    """Lifetime inputs for the K006-to-K012 short-distance normalization."""

    kshort_lifetime_ps: float = 89.54
    klong_lifetime_ps: float = 51160.0
    citation: str = KSHORT_MUMU_LIFETIME_INPUTS_CITATION

    def __post_init__(self) -> None:
        for name in ("kshort_lifetime_ps", "klong_lifetime_ps"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")

    @property
    def lifetime_ratio(self) -> float:
        """Return ``tau_KS / tau_KL``."""

        return float(self.kshort_lifetime_ps / self.klong_lifetime_ps)


@dataclass(frozen=True)
class KShortMuMuShortDistanceResult:
    """Short-distance branching-ratio proxy for ``K_S -> mu+ mu-``."""

    model_label: str
    input_bundle: str
    branching_fraction: float
    sm_branching_fraction: float
    np_shift_branching_fraction: float
    klong_branching_fraction: float
    klong_sm_branching_fraction: float
    klong_np_shift_branching_fraction: float
    lifetime_ratio: float
    kshort_lifetime_ps: float
    klong_lifetime_ps: float
    kappa_mu: float
    p_c_y: float
    y_t: float
    y_eff: complex
    y_np_total: complex
    short_distance_amplitude: float
    lambda_wolfenstein: float
    lambda_c: complex
    lambda_t: complex
    wilsons: RareKaonDileptonWilsonCoefficients | None = None
    diagnostics: Mapping[str, float | complex | str | bool] = field(
        default_factory=dict
    )


def kshort_mumu_lifetime_inputs_default() -> KShortMuMuLifetimeInputs:
    """Return the K_S/K_L lifetime inputs for the short-distance proxy."""

    return KShortMuMuLifetimeInputs()


def _kshort_mumu_sd_from_y_np(
    y_np_total: complex,
    *,
    inputs: RareKaonDileptonSMInputs,
    klong_result: KLongMuMuShortDistanceResult,
    lifetime_inputs: KShortMuMuLifetimeInputs,
    wilsons: RareKaonDileptonWilsonCoefficients | None = None,
) -> KShortMuMuShortDistanceResult:
    ratio = float(lifetime_inputs.lifetime_ratio)
    factors = rare_kaon_dilepton_ckm_factors(inputs)
    lam = float(factors.lambda_wolfenstein)
    y_eff = factors.lambda_t * inputs.y_t + complex(y_np_total)
    sm_y_eff = factors.lambda_t * inputs.y_t
    short_distance_amplitude = (
        y_eff.imag / lam**5 - factors.lambda_c.imag * inputs.p_c_y / lam
    )
    sm_short_distance_amplitude = (
        sm_y_eff.imag / lam**5 - factors.lambda_c.imag * inputs.p_c_y / lam
    )
    kappa = rare_kaon_dilepton_kappa_mu(inputs)
    branching_fraction = float(ratio * kappa * short_distance_amplitude**2)
    sm_branching_fraction = float(ratio * kappa * sm_short_distance_amplitude**2)
    np_shift = float(branching_fraction - sm_branching_fraction)

    diagnostics: dict[str, float | complex | str | bool] = dict(
        klong_result.diagnostics
    )
    if "sm_short_distance_amplitude" in diagnostics:
        diagnostics["klong_sm_real_short_distance_amplitude"] = diagnostics[
            "sm_short_distance_amplitude"
        ]
    diagnostics.update(
        {
            "matching_assumption": RARE_KAON_KSHORT_MUMU_RS_MATCHING_ASSUMPTION_V1,
            "k006_matching_assumption_reused": (
                RARE_KAON_DILEPTON_RS_MATCHING_ASSUMPTION_V1
            ),
            "parametrization_citation": (
                RARE_KAON_KSHORT_MUMU_PARAMETRIZATION_CITATION
            ),
            "cp_projection_citation": KSHORT_MUMU_CP_PROJECTION_CITATION,
            "semileptonic_running_diagnostic": (
                RARE_KAON_KSHORT_MUMU_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1
            ),
            "low_scale_qcd_running_factor": 1.0,
            "uses_k006_muon_sd_evaluator": False,
            "uses_k006_muon_sd_wilson_proxy": True,
            "uses_lifetime_rescaled_short_distance_proxy": False,
            "uses_imaginary_ks_ell0_projection": True,
            "kshort_ell0_cp_projection_implemented": True,
            "kshort_ell0_projection_formula": (
                "(tau_KS/tau_KL) * kappa_mu * "
                "[Im(lambda_t*Y_t + Y_NP)/lambda^5 "
                "- Im(lambda_c)*P_c(Y)/lambda]^2"
            ),
            "kshort_ell0_projection_source": (
                "Im[-lambda_c Y_c + lambda_t C10]^2"
            ),
            "short_distance_only": True,
            "long_distance_dominated_total_rate": True,
            "total_long_distance_not_modelled": True,
            "ks_ell0_time_dependent_clean_extraction_not_implemented": True,
            "sm_context_total_not_sd_anchor": True,
            "needs_human_physics": RARE_KAON_KSHORT_MUMU_RS_MATCHING_ASSUMPTION_V1,
            "lifetime_inputs_citation": lifetime_inputs.citation,
            "kshort_lifetime_ps": float(lifetime_inputs.kshort_lifetime_ps),
            "klong_lifetime_ps": float(lifetime_inputs.klong_lifetime_ps),
            "kshort_to_klong_lifetime_ratio": ratio,
            "klong_short_distance_branching_fraction": float(
                klong_result.branching_fraction
            ),
            "klong_sm_short_distance_branching_fraction": float(
                klong_result.sm_branching_fraction
            ),
            "klong_np_shift_branching_fraction": float(
                klong_result.np_shift_branching_fraction
            ),
            "klong_real_short_distance_amplitude": float(
                klong_result.short_distance_amplitude
            ),
            "sm_short_distance_amplitude": float(sm_short_distance_amplitude),
            "kshort_imaginary_short_distance_amplitude": float(
                short_distance_amplitude
            ),
            "sm_y_eff": complex(sm_y_eff),
            "y_eff": complex(y_eff),
            "y_np_total": complex(y_np_total),
        }
    )
    if wilsons is not None:
        diagnostics.update(
            {
                "m_kk_gev": float(wilsons.M_KK),
                "matching_scale_gev": float(wilsons.matching_scale),
                "left_sd_coupling": complex(wilsons.left_sd_coupling),
                "right_sd_coupling": complex(wilsons.right_sd_coupling),
                "left_sd_overlap": complex(wilsons.left_sd_overlap),
                "right_sd_overlap": complex(wilsons.right_sd_overlap),
                "left_quark_delta": complex(wilsons.left_quark_delta),
                "right_quark_delta": complex(wilsons.right_quark_delta),
                "muon_axial_delta": float(wilsons.muon_axial_delta),
                "y_np_left": complex(wilsons.y_np_left),
                "y_np_right": complex(wilsons.y_np_right),
                "y_np_total": complex(wilsons.y_np_total),
            }
        )

    return KShortMuMuShortDistanceResult(
        model_label=RARE_KAON_KSHORT_MUMU_MODEL_V1,
        input_bundle=klong_result.input_bundle,
        branching_fraction=branching_fraction,
        sm_branching_fraction=sm_branching_fraction,
        np_shift_branching_fraction=np_shift,
        klong_branching_fraction=float(klong_result.branching_fraction),
        klong_sm_branching_fraction=float(klong_result.sm_branching_fraction),
        klong_np_shift_branching_fraction=float(
            klong_result.np_shift_branching_fraction
        ),
        lifetime_ratio=ratio,
        kshort_lifetime_ps=float(lifetime_inputs.kshort_lifetime_ps),
        klong_lifetime_ps=float(lifetime_inputs.klong_lifetime_ps),
        kappa_mu=float(kappa),
        p_c_y=float(inputs.p_c_y),
        y_t=float(inputs.y_t),
        y_eff=complex(y_eff),
        y_np_total=complex(y_np_total),
        short_distance_amplitude=float(short_distance_amplitude),
        lambda_wolfenstein=lam,
        lambda_c=complex(factors.lambda_c),
        lambda_t=complex(factors.lambda_t),
        wilsons=wilsons,
        diagnostics=diagnostics,
    )


def kshort_mumu_short_distance_sm(
    inputs: RareKaonDileptonSMInputs | None = None,
    *,
    lifetime_inputs: KShortMuMuLifetimeInputs | None = None,
) -> KShortMuMuShortDistanceResult:
    """Evaluate the SM-limit ``K_S -> mu+ mu-`` short-distance proxy."""

    p = rare_kaon_dilepton_default_sm_inputs() if inputs is None else inputs
    lifetimes = (
        kshort_mumu_lifetime_inputs_default()
        if lifetime_inputs is None
        else lifetime_inputs
    )
    return _kshort_mumu_sd_from_y_np(
        0.0j,
        inputs=p,
        klong_result=klong_mumu_short_distance_sm(p),
        lifetime_inputs=lifetimes,
    )


def kshort_mumu_short_distance_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
    lifetime_inputs: KShortMuMuLifetimeInputs | None = None,
) -> KShortMuMuShortDistanceResult:
    """Evaluate ``BR(K_S -> mu+ mu-)_SD`` from mass-basis couplings."""

    p = rare_kaon_dilepton_default_sm_inputs() if inputs is None else inputs
    lifetimes = (
        kshort_mumu_lifetime_inputs_default()
        if lifetime_inputs is None
        else lifetime_inputs
    )
    wilsons = rare_kaon_dilepton_wilsons_from_couplings(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=p,
    )
    return _kshort_mumu_sd_from_y_np(
        wilsons.y_np_total,
        inputs=p,
        klong_result=klong_mumu_short_distance_from_couplings(
            couplings,
            m_kk_gev=m_kk_gev,
            inputs=p,
        ),
        lifetime_inputs=lifetimes,
        wilsons=wilsons,
    )


def kshort_mumu_short_distance_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    lepton: str = "mu",
    matching_scale_gev: float | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
    lifetime_inputs: KShortMuMuLifetimeInputs | None = None,
) -> KShortMuMuShortDistanceResult:
    """Evaluate ``BR(K_S -> mu+mu-)_SD`` from Phase-3a RS C10/C10p."""

    p = rare_kaon_dilepton_default_sm_inputs() if inputs is None else inputs
    lifetimes = (
        kshort_mumu_lifetime_inputs_default()
        if lifetime_inputs is None
        else lifetime_inputs
    )
    coeff = _rs_semileptonic_coeff(source, lepton=lepton)
    wilsons = rare_kaon_y_wilsons_from_rs_semileptonic(
        source,
        lepton=lepton,
        matching_scale_gev=matching_scale_gev,
        inputs=p,
    )
    result = _kshort_mumu_sd_from_y_np(
        wilsons.y_np_total,
        inputs=p,
        klong_result=klong_mumu_short_distance_from_rs_semileptonic_wilsons(
            source,
            lepton=lepton,
            matching_scale_gev=matching_scale_gev,
            inputs=p,
        ),
        lifetime_inputs=lifetimes,
        wilsons=wilsons,
    )
    return _tag_rs_result(result, coeff, inputs=p)


__all__ = [
    "RARE_KAON_KSHORT_MUMU_MODEL_V1",
    "KSHORT_MUMU_CP_PROJECTION_CITATION",
    "RARE_KAON_KSHORT_MUMU_PARAMETRIZATION_CITATION",
    "RARE_KAON_KSHORT_MUMU_RS_MATCHING_ASSUMPTION_V1",
    "RARE_KAON_KSHORT_MUMU_SEMILEPTONIC_RUNNING_DIAGNOSTIC_V1",
    "KSHORT_MUMU_LIFETIME_INPUTS_CITATION",
    "KShortMuMuLifetimeInputs",
    "KShortMuMuShortDistanceResult",
    "kshort_mumu_lifetime_inputs_default",
    "kshort_mumu_short_distance_sm",
    "kshort_mumu_short_distance_from_couplings",
    "kshort_mumu_short_distance_from_rs_semileptonic_wilsons",
]
