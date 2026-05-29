"""Electroweak oblique-parameter likelihood and RS proxy.

The rigorous, model-independent part here is the correlated Gaussian
comparison of a predicted ``(S, T)`` point to a published electroweak global
fit with ``U`` fixed.  The Standard Model reference is ``S = T = U = 0`` by
construction.

RS matching status
------------------
NEEDS-HUMAN-PHYSICS.  A complete warped-model prediction requires the full
electroweak KK gauge sector, brane terms, custodial representation choices,
fermion embeddings, Higgs localization, and the mapping to vacuum-polarization
self energies or SMEFT bosonic operators.  The helper
``rs_minimal_oblique_proxy`` therefore implements only the documented classic
minimal-RS scaling

    Delta S = c_S v^2 / M_KK^2,
    Delta T = [pi L / (2 cos^2 theta_W)] v^2 / M_KK^2,
    Delta U = 0,

where ``L = k pi r_c`` is the warped volume.  ``c_S`` is supplied by the
catalog anchor, while the ``T`` coefficient is the standard volume-enhanced
minimal-RS parametric form.  This proxy is suitable for catalog screening and
must not be read as a full RS electroweak fit.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

OBLIQUE_STU_RS_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: EW001 uses a minimal-RS oblique proxy, "
    "Delta S = c_S v^2/M_KK^2 and Delta T = pi*L/(2*cW^2)*v^2/M_KK^2 "
    "with U=0. The full EW KK gauge sector, custodial structure, brane terms, "
    "fermion embeddings, Higgs localization, and EW global-fit matching are "
    "not available on ParameterPoint."
)

OBLIQUE_STU_LIKELIHOOD_V1 = "correlated_gaussian_ST_U_fixed_v1"
DEFAULT_HIGGS_VEV_GEV = 246.21965
DEFAULT_SIN2_THETA_W = 0.23122
DEFAULT_RS_VOLUME_LOG = 35.0
CHI2_2DOF_95 = 5.991464547107979


@dataclass(frozen=True)
class ObliqueSTFit:
    """Correlated ``(S, T)`` global-fit ellipse with ``U`` fixed."""

    s_central: float
    t_central: float
    sigma_s: float
    sigma_t: float
    rho_st: float
    u_fixed: float = 0.0
    chi2_budget: float = CHI2_2DOF_95
    confidence_level: float = 0.95
    source: str | None = None
    source_url: str | None = None
    covariance_source: str | None = None

    def __post_init__(self) -> None:
        _finite(self.s_central, "s_central")
        _finite(self.t_central, "t_central")
        _finite(self.u_fixed, "u_fixed")
        _positive(self.sigma_s, "sigma_s")
        _positive(self.sigma_t, "sigma_t")
        _positive(self.chi2_budget, "chi2_budget")
        _finite(self.confidence_level, "confidence_level")
        rho = float(self.rho_st)
        if not math.isfinite(rho) or not -1.0 < rho < 1.0:
            raise ValueError("rho_st must be finite and strictly between -1 and 1")
        cl = float(self.confidence_level)
        if not 0.0 < cl < 1.0:
            raise ValueError("confidence_level must be between zero and one")

    @property
    def covariance(self) -> tuple[tuple[float, float], tuple[float, float]]:
        cov_st = float(self.rho_st * self.sigma_s * self.sigma_t)
        return (
            (float(self.sigma_s * self.sigma_s), cov_st),
            (cov_st, float(self.sigma_t * self.sigma_t)),
        )

    @property
    def inverse_covariance(self) -> tuple[tuple[float, float], tuple[float, float]]:
        var_s = float(self.sigma_s * self.sigma_s)
        var_t = float(self.sigma_t * self.sigma_t)
        cov_st = float(self.rho_st * self.sigma_s * self.sigma_t)
        det = var_s * var_t - cov_st * cov_st
        if det <= 0.0:
            raise ValueError("ST covariance matrix must be positive definite")
        return ((var_t / det, -cov_st / det), (-cov_st / det, var_s / det))


@dataclass(frozen=True)
class ObliqueSTUShift:
    """Predicted oblique-parameter shift."""

    s: float
    t: float
    u: float
    m_kk_gev: float
    s_coefficient: float
    t_coefficient: float
    higgs_vev_gev: float
    sin2_theta_w: float
    rs_volume_log: float
    matching_assumption: str = OBLIQUE_STU_RS_PROXY_V1

    def __post_init__(self) -> None:
        for name in (
            "s",
            "t",
            "u",
            "s_coefficient",
            "t_coefficient",
            "higgs_vev_gev",
            "sin2_theta_w",
            "rs_volume_log",
        ):
            _finite(getattr(self, name), name)
        _positive(self.m_kk_gev, "m_kk_gev")


@dataclass(frozen=True)
class ObliqueSTComparison:
    """Comparison of one ``(S, T)`` prediction to the fit ellipse."""

    prediction: ObliqueSTUShift
    fit: ObliqueSTFit
    delta_s_from_fit_center: float
    delta_t_from_fit_center: float
    chi2: float
    chi2_budget: float
    ratio_to_budget: float
    passes: bool
    diagnostics: Mapping[str, float | str | tuple[tuple[float, float], tuple[float, float]]]


def minimal_rs_t_coefficient(
    *,
    rs_volume_log: float = DEFAULT_RS_VOLUME_LOG,
    sin2_theta_w: float = DEFAULT_SIN2_THETA_W,
) -> float:
    """Return the volume-enhanced minimal-RS ``Delta T`` coefficient."""
    volume = _positive(rs_volume_log, "rs_volume_log")
    sin2 = _finite(sin2_theta_w, "sin2_theta_w")
    cos2 = 1.0 - sin2
    if cos2 <= 0.0:
        raise ValueError("sin2_theta_w must be less than one")
    return float(math.pi * volume / (2.0 * cos2))


def rs_minimal_oblique_proxy(
    *,
    m_kk_gev: float,
    s_coefficient: float,
    higgs_vev_gev: float = DEFAULT_HIGGS_VEV_GEV,
    sin2_theta_w: float = DEFAULT_SIN2_THETA_W,
    rs_volume_log: float = DEFAULT_RS_VOLUME_LOG,
) -> ObliqueSTUShift:
    """Return the documented minimal-RS ``S,T,U`` proxy at ``M_KK``."""
    mass = _positive(m_kk_gev, "m_kk_gev")
    vev = _positive(higgs_vev_gev, "higgs_vev_gev")
    coeff_s = _finite(s_coefficient, "s_coefficient")
    coeff_t = minimal_rs_t_coefficient(
        rs_volume_log=rs_volume_log,
        sin2_theta_w=sin2_theta_w,
    )
    scale = float((vev / mass) ** 2)
    return ObliqueSTUShift(
        s=float(coeff_s * scale),
        t=float(coeff_t * scale),
        u=0.0,
        m_kk_gev=float(mass),
        s_coefficient=float(coeff_s),
        t_coefficient=float(coeff_t),
        higgs_vev_gev=float(vev),
        sin2_theta_w=float(sin2_theta_w),
        rs_volume_log=float(rs_volume_log),
    )


def compare_oblique_st_to_fit(
    *,
    s: float,
    t: float,
    fit: ObliqueSTFit,
) -> tuple[float, float, float, bool]:
    """Return ``(chi2, ratio, budget, passes)`` for a point in the ST plane."""
    s_pred = _finite(s, "s")
    t_pred = _finite(t, "t")
    ds = float(s_pred - fit.s_central)
    dt = float(t_pred - fit.t_central)
    inv = fit.inverse_covariance
    chi2 = float(ds * (inv[0][0] * ds + inv[0][1] * dt) + dt * (inv[1][0] * ds + inv[1][1] * dt))
    budget = float(fit.chi2_budget)
    ratio = float(chi2 / budget)
    return chi2, ratio, budget, bool(ratio <= 1.0)


def evaluate_rs_oblique_proxy(
    *,
    m_kk_gev: float,
    fit: ObliqueSTFit,
    s_coefficient: float,
    higgs_vev_gev: float = DEFAULT_HIGGS_VEV_GEV,
    sin2_theta_w: float = DEFAULT_SIN2_THETA_W,
    rs_volume_log: float = DEFAULT_RS_VOLUME_LOG,
) -> ObliqueSTComparison:
    """Evaluate the documented RS proxy against an ``S,T`` fit ellipse."""
    prediction = rs_minimal_oblique_proxy(
        m_kk_gev=m_kk_gev,
        s_coefficient=s_coefficient,
        higgs_vev_gev=higgs_vev_gev,
        sin2_theta_w=sin2_theta_w,
        rs_volume_log=rs_volume_log,
    )
    chi2, ratio, budget, passes = compare_oblique_st_to_fit(
        s=prediction.s,
        t=prediction.t,
        fit=fit,
    )
    diagnostics = {
        "likelihood_model": OBLIQUE_STU_LIKELIHOOD_V1,
        "matching_assumption": OBLIQUE_STU_RS_PROXY_V1,
        "covariance": fit.covariance,
        "inverse_covariance": fit.inverse_covariance,
        "fit_s_central": float(fit.s_central),
        "fit_t_central": float(fit.t_central),
        "fit_u_fixed": float(fit.u_fixed),
        "sigma_s": float(fit.sigma_s),
        "sigma_t": float(fit.sigma_t),
        "rho_st": float(fit.rho_st),
        "confidence_level": float(fit.confidence_level),
    }
    return ObliqueSTComparison(
        prediction=prediction,
        fit=fit,
        delta_s_from_fit_center=float(prediction.s - fit.s_central),
        delta_t_from_fit_center=float(prediction.t - fit.t_central),
        chi2=float(chi2),
        chi2_budget=float(budget),
        ratio_to_budget=float(ratio),
        passes=passes,
        diagnostics=diagnostics,
    )


def _finite(value: float, name: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _positive(value: float, name: str) -> float:
    number = _finite(value, name)
    if number <= 0.0:
        raise ValueError(f"{name} must be positive")
    return number


__all__ = [
    "CHI2_2DOF_95",
    "DEFAULT_HIGGS_VEV_GEV",
    "DEFAULT_RS_VOLUME_LOG",
    "DEFAULT_SIN2_THETA_W",
    "OBLIQUE_STU_LIKELIHOOD_V1",
    "OBLIQUE_STU_RS_PROXY_V1",
    "ObliqueSTComparison",
    "ObliqueSTFit",
    "ObliqueSTUShift",
    "compare_oblique_st_to_fit",
    "evaluate_rs_oblique_proxy",
    "minimal_rs_t_coefficient",
    "rs_minimal_oblique_proxy",
]
