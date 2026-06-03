"""Catalog adapter for Phase-5a RS charged-current epsilon shifts."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

from quarkConstraints.rs_charged_current import (
    RSChargedCurrentCouplings,
    RS_CHARGED_CURRENT_MATCHING_ASSUMPTION_V1,
)

__all__ = [
    "CHARGED_CURRENT_MINIMAL_LH_STATUS",
    "CHARGED_CURRENT_NONMINIMAL_NEEDS_HUMAN",
    "B025_PARTIAL_MATCHING_STATUS",
    "EW003_CKM_TENSION_COVARIANCE_NEEDS_HUMAN",
    "EW003_CKM_TENSION_DIAGNOSTIC_STATUS",
    "ChargedCurrentEpsilon",
    "ChargedCurrentFirstRowShift",
    "RSChargedCurrentCouplings",
    "charged_current_epsilon",
    "charged_current_ckm_tension_diagnostics",
    "charged_current_light_epsilon",
    "charged_current_source_diagnostics",
    "ckm_first_row_delta_np",
    "shifted_abs_ckm",
    "shifted_branching_fraction",
    "shifted_lfu_ratio",
]

CHARGED_CURRENT_MINIMAL_LH_STATUS = (
    "RIGOROUS: minimal-RS left-handed W/W' charged-current epsilon from "
    "rs_charged_current, including the stored delta_G_F/G_F subtraction."
)
CHARGED_CURRENT_NONMINIMAL_NEEDS_HUMAN = (
    "NEEDS-HUMAN-PHYSICS: charged-Higgs, right-handed W, scalar/tensor WET, "
    "heavy-neutrino, radiative, and mode-integration effects are not included "
    "in the minimal left-handed rs_charged_current epsilon."
)
B025_PARTIAL_MATCHING_STATUS = (
    "PARTIAL: the minimal vector LFU-ratio tree is evaluated from epsilon_cb, "
    "but scalar/RH WET coefficients and B -> D form-factor integration are not "
    "built."
)
EW003_CKM_TENSION_COVARIANCE_NEEDS_HUMAN = (
    "NEEDS-HUMAN-PHYSICS: a rigorous inclusive-vs-exclusive |V_cb|/|V_ub| "
    "charged-current treatment needs a covariance/scheme input separating "
    "inclusive and exclusive CKM/form-factor uncertainties. The v1 "
    "rs_charged_current epsilon is common to both determinations, so it is "
    "reported only as a diagnostic and is not a naive NP veto."
)
EW003_CKM_TENSION_DIAGNOSTIC_STATUS = (
    "PARTIAL/data-level: EW003 keeps the data-level inclusive-vs-exclusive "
    "pull. Minimal-LH epsilon_cb and epsilon_ub are surfaced as diagnostics; "
    "a universal charged-current rescaling cancels in each inclusive/exclusive "
    "ratio."
)

_UP_INDEX = {"u": 0, "c": 1, "t": 2}
_DOWN_INDEX = {"d": 0, "s": 1, "b": 2}
_LEPTON_INDEX = {"e": 0, "mu": 1, "tau": 2}


@dataclass(frozen=True)
class ChargedCurrentEpsilon:
    """One charged-current amplitude shift and its observable multiplier."""

    up: str
    down: str
    lepton: str
    epsilon: complex

    @property
    def real(self) -> float:
        return float(self.epsilon.real)

    @property
    def abs_multiplier(self) -> float:
        return float(abs(1.0 + self.epsilon))

    @property
    def rate_multiplier(self) -> float:
        return float(self.abs_multiplier**2)

    @property
    def diagnostics(self) -> dict[str, Any]:
        return {
            "up_flavor": self.up,
            "down_flavor": self.down,
            "lepton_flavor": self.lepton,
            "epsilon": complex(self.epsilon),
            "epsilon_real": float(self.real),
            "epsilon_imag": float(self.epsilon.imag),
            "abs_1_plus_epsilon": float(self.abs_multiplier),
            "rate_multiplier": float(self.rate_multiplier),
        }


@dataclass(frozen=True)
class ChargedCurrentFirstRowShift:
    """First-row CKM unitarity shift from charged-current epsilons."""

    delta_ckm_np: float
    epsilon_ud_e: ChargedCurrentEpsilon
    epsilon_us_light: ChargedCurrentEpsilon
    epsilon_ub_light: ChargedCurrentEpsilon | None
    terms: Mapping[str, float]


def charged_current_epsilon(
    source: RSChargedCurrentCouplings,
    *,
    up: str,
    down: str,
    lepton: str,
) -> ChargedCurrentEpsilon:
    """Return ``epsilon_ij^a`` for a named quark/lepton channel."""

    i = _index(_UP_INDEX, up, "up")
    j = _index(_DOWN_INDEX, down, "down")
    a = _index(_LEPTON_INDEX, lepton, "lepton")
    if not isinstance(source, RSChargedCurrentCouplings):
        raise TypeError(
            "rs_charged_current must be an RSChargedCurrentCouplings instance"
        )
    value = _finite_complex(source.epsilon_ij_a(i, j, a), "epsilon")
    return ChargedCurrentEpsilon(up=up, down=down, lepton=lepton, epsilon=value)


def charged_current_light_epsilon(
    source: RSChargedCurrentCouplings,
    *,
    up: str,
    down: str,
) -> ChargedCurrentEpsilon:
    """Return the e/mu average epsilon used where YAML has no mode weights."""

    eps_e = charged_current_epsilon(source, up=up, down=down, lepton="e")
    eps_mu = charged_current_epsilon(source, up=up, down=down, lepton="mu")
    return ChargedCurrentEpsilon(
        up=up,
        down=down,
        lepton="light_average_e_mu",
        epsilon=_finite_complex(0.5 * (eps_e.epsilon + eps_mu.epsilon), "epsilon_light"),
    )


def ckm_first_row_delta_np(
    source: RSChargedCurrentCouplings,
    *,
    vud_abs: float,
    vus_abs: float,
    vub_abs: float | None = None,
) -> ChargedCurrentFirstRowShift:
    """Return ``Delta_CKM_NP ~= 2 sum |V_ij|^2 Re(epsilon_ij)``."""

    vud = _finite_real(vud_abs, "vud_abs")
    vus = _finite_real(vus_abs, "vus_abs")
    vub = None if vub_abs is None else _finite_real(vub_abs, "vub_abs")
    for name, value in (("vud_abs", vud), ("vus_abs", vus)):
        if value < 0.0:
            raise ValueError(f"{name} must be non-negative")
    if vub is not None and vub < 0.0:
        raise ValueError("vub_abs must be non-negative")

    eps_ud = charged_current_epsilon(source, up="u", down="d", lepton="e")
    eps_us = charged_current_light_epsilon(source, up="u", down="s")
    terms = {
        "ud_e": float(2.0 * vud * vud * eps_ud.real),
        "us_light_average_e_mu": float(2.0 * vus * vus * eps_us.real),
    }
    eps_ub = None
    if vub is not None:
        eps_ub = charged_current_light_epsilon(source, up="u", down="b")
        terms["ub_light_average_e_mu"] = float(2.0 * vub * vub * eps_ub.real)
    return ChargedCurrentFirstRowShift(
        delta_ckm_np=float(sum(terms.values())),
        epsilon_ud_e=eps_ud,
        epsilon_us_light=eps_us,
        epsilon_ub_light=eps_ub,
        terms=terms,
    )


def charged_current_ckm_tension_diagnostics(
    source: RSChargedCurrentCouplings,
) -> dict[str, Any]:
    """Return EW003 diagnostics for common ``|V_cb|``/``|V_ub|`` rescalings.

    EW003 compares inclusive and exclusive determinations of the same CKM
    element.  With only the v1 common charged-current epsilon available, both
    determinations receive the same ``|1 + epsilon|`` factor, so their ratio
    and Gaussian tension are unchanged.  A differential treatment needs an
    external covariance/scheme model and is intentionally not inferred here.
    """

    eps_cb = charged_current_light_epsilon(source, up="c", down="b")
    eps_ub = charged_current_light_epsilon(source, up="u", down="b")
    cb_ratio_multiplier = _common_rescaling_ratio_multiplier(eps_cb)
    ub_ratio_multiplier = _common_rescaling_ratio_multiplier(eps_ub)
    return {
        "ew003_matching_coverage": "PARTIAL",
        "ew003_charged_current_status": EW003_CKM_TENSION_DIAGNOSTIC_STATUS,
        "ew003_covariance_scheme_status": (
            EW003_CKM_TENSION_COVARIANCE_NEEDS_HUMAN
        ),
        "epsilon_cb": eps_cb.diagnostics,
        "epsilon_ub": eps_ub.diagnostics,
        "epsilon_cb_light_average_e_mu": eps_cb.epsilon,
        "epsilon_ub_light_average_e_mu": eps_ub.epsilon,
        "epsilon_cb_abs_1_plus": float(eps_cb.abs_multiplier),
        "epsilon_ub_abs_1_plus": float(eps_ub.abs_multiplier),
        "universal_cc_ratio_multiplier_vcb": float(cb_ratio_multiplier),
        "universal_cc_ratio_multiplier_vub": float(ub_ratio_multiplier),
        "universal_cc_pull_multiplier_vcb": float(cb_ratio_multiplier),
        "universal_cc_pull_multiplier_vub": float(ub_ratio_multiplier),
        "universal_cc_pull_shift_sigma_vcb": float(cb_ratio_multiplier - 1.0),
        "universal_cc_pull_shift_sigma_vub": float(ub_ratio_multiplier - 1.0),
        "universal_cc_cancellation_documented": True,
        "universal_cc_cancellation_applied_to_scalar": False,
    }


def shifted_abs_ckm(ckm_abs: float, epsilon: ChargedCurrentEpsilon) -> float:
    """Return ``|V|_app = |V| |1 + epsilon|``."""

    value = _finite_real(ckm_abs, "ckm_abs")
    if value < 0.0:
        raise ValueError("ckm_abs must be non-negative")
    return float(value * epsilon.abs_multiplier)


def shifted_lfu_ratio(
    sm_ratio: float,
    numerator_epsilon: ChargedCurrentEpsilon,
    denominator_epsilon: ChargedCurrentEpsilon,
) -> float:
    """Return ``R_SM |1+eps_num|^2 / |1+eps_den|^2``."""

    sm = _positive_real(sm_ratio, "sm_ratio")
    denominator = denominator_epsilon.rate_multiplier
    if denominator <= 0.0:
        raise ValueError("denominator |1+epsilon|^2 must be positive")
    return float(sm * numerator_epsilon.rate_multiplier / denominator)


def shifted_branching_fraction(
    sm_branching_fraction: float,
    epsilon: ChargedCurrentEpsilon,
) -> float:
    """Return ``BR_SM |1 + epsilon|^2``."""

    return float(_positive_real(sm_branching_fraction, "sm_branching_fraction") * epsilon.rate_multiplier)


def charged_current_source_diagnostics(
    source: RSChargedCurrentCouplings,
) -> dict[str, Any]:
    """Return reusable scalar diagnostics for a charged-current source."""

    if not isinstance(source, RSChargedCurrentCouplings):
        raise TypeError(
            "rs_charged_current must be an RSChargedCurrentCouplings instance"
        )
    return {
        "charged_current_model_label": source.model_label,
        "charged_current_input_bundle": source.input_bundle,
        "charged_current_matching_assumption": (
            RS_CHARGED_CURRENT_MATCHING_ASSUMPTION_V1
        ),
        "kk_ew_mass_gev": float(source.kk_ew_mass_gev),
        "m_w_gev": float(source.m_w_gev),
        "m_wprime_gev": float(source.m_wprime_gev),
        "eta_W": float(source.eta_W),
        "delta_G_F_over_G_F": float(source.delta_G_F_over_G_F),
        "minimal_lh_vector_matching_status": CHARGED_CURRENT_MINIMAL_LH_STATUS,
    }


def _index(mapping: Mapping[str, int], label: str, kind: str) -> int:
    key = str(label)
    try:
        return mapping[key]
    except KeyError as exc:
        raise ValueError(f"unknown {kind} flavor {label!r}") from exc


def _finite_complex(value: complex, name: str) -> complex:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    return number


def _finite_real(value: float, name: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _positive_real(value: float, name: str) -> float:
    number = _finite_real(value, name)
    if number <= 0.0:
        raise ValueError(f"{name} must be positive")
    return number


def _common_rescaling_ratio_multiplier(epsilon: ChargedCurrentEpsilon) -> float:
    multiplier = float(epsilon.abs_multiplier)
    if multiplier <= 0.0:
        raise ValueError("|1+epsilon| must be positive for ratio diagnostics")
    return float(multiplier / multiplier)
