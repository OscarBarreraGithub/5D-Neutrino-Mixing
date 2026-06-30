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
minimal-RS scaling, with a tree-level custodial proxy coefficient selected
when ``ew_model="custodial_rs_plr"``:

    Delta S = c_S v^2 / M_KK^2,
    Delta T_minimal = x_1^2 [pi L / (2 cos^2 theta_W)] v^2 / M_KK^2,
    Delta T_custodial = -x_1^2 [pi / (4 cos^2 theta_W L)] v^2 / M_KK^2,
    Delta U = 0,

where ``L = k pi r_c`` is the warped volume and ``x_1 ~ 2.4487`` is the first
gauge-KK Bessel root (``M_KK = x_1 Lambda_IR``).  ``c_S`` is supplied by the
catalog anchor and is calibrated in the PHYSICAL-M_KK convention.  The
volume-enhanced ``T`` coefficient ``pi L/(2 cos^2)`` is derived in the
GEOMETRIC-Lambda_IR convention, so it is multiplied by ``x_1^2`` here to put
Delta T in the same physical-M_KK convention as Delta S (PLAN §5/M2); without
this factor Delta T was ~x_1^2 ~ 6 too small (anti-conservative).  This proxy
is suitable for catalog screening and must not be read as a full RS
electroweak fit.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping

OBLIQUE_STU_RS_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: EW001 uses an RS oblique proxy, "
    "Delta S = c_S v^2/M_KK^2 with U=0, and model-selected Delta T: "
    "minimal_rs uses x1^2*pi*L/(2*cW^2)*v^2/M_KK^2 while custodial_rs_plr uses "
    "-x1^2*pi/(4*cW^2*L)*v^2/M_KK^2 (x1 = first gauge-KK root, M_KK = x1*Lambda_IR; "
    "the x1^2 factor puts the geometric-Lambda_IR Delta T coefficient in the same "
    "physical-M_KK convention as the c_S calibration). "
    "The full EW KK gauge sector, custodial structure, brane terms, "
    "fermion embeddings, Higgs localization, and EW global-fit matching are "
    "not available on ParameterPoint."
)

OBLIQUE_STU_LIKELIHOOD_V1 = "correlated_gaussian_ST_U_fixed_v1"
DEFAULT_HIGGS_VEV_GEV = 246.21965
DEFAULT_SIN2_THETA_W = 0.23122
DEFAULT_RS_VOLUME_LOG = 35.0
CHI2_2DOF_95 = 5.991464547107979

# First electroweak gauge-KK Neumann-Neumann Bessel root, x_1 ~ 2.4487, so the
# physical KK mass is M_KK = x_1 * Lambda_IR (mirrors quarkConstraints.scales
# .GAUGE_KK_ROOT_NN).  The CGHNP volume-enhanced Delta T coefficient
# pi*L/(2 cos^2) is derived in the GEOMETRIC-Lambda_IR convention (its natural
# scale is Lambda_IR), but the proxy applies the SAME physical (v/M_KK)^2 scale
# as the Delta S coefficient (which IS calibrated in the physical-M_KK
# convention).  To evaluate Delta T in a single, documented convention that is
# consistent with the Delta S calibration we express the Delta T coefficient in
# the physical-M_KK convention by multiplying it by x_1^2: this is the
# "rescale the coefficient" option of PLAN §5/M2, equivalent to using
# (v/Lambda_IR)^2 = x_1^2 (v/M_KK)^2 for the Delta T term.  Without it Delta T
# was ~x_1^2 ~ 6.0 too small (anti-conservative).
GAUGE_KK_ROOT_NN = 2.448687135269161
PHYSICAL_MKK_OVER_LAMBDA_IR_SQ = GAUGE_KK_ROOT_NN * GAUGE_KK_ROOT_NN
MINIMAL_RS_EW_MODEL = "minimal_rs"
# SU(2)_R only: the custodial SU(2)_R protects the oblique T parameter (CGHNP
# Eq. 153) but WITHOUT the discrete P_LR, so the Z b_L vertex is left unprotected
# (it stays at its minimal value in rs_ew_couplings.py).  This isolates what the
# SU(2)_R custodial gauge structure buys (the T cure) from what the extra P_LR
# discrete symmetry buys (the Z b_L cure).  For the oblique sector S,T,U this is
# identical to custodial_rs_plr, because P_LR does not act on the oblique params.
CUSTODIAL_RS_SU2R_EW_MODEL = "custodial_rs_su2r"
CUSTODIAL_RS_PLR_EW_MODEL = "custodial_rs_plr"
SUPPORTED_RS_EW_MODELS = (
    MINIMAL_RS_EW_MODEL,
    CUSTODIAL_RS_SU2R_EW_MODEL,
    CUSTODIAL_RS_PLR_EW_MODEL,
)


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
    ew_model: str = MINIMAL_RS_EW_MODEL
    t_tree: float | None = None
    delta_t_loop: float = 0.0
    loop_metadata: Mapping[str, Any] = field(default_factory=dict)

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
        _finite(self.delta_t_loop, "delta_t_loop")
        if self.t_tree is not None:
            _finite(self.t_tree, "t_tree")
        _positive(self.m_kk_gev, "m_kk_gev")
        _validate_ew_model(self.ew_model)
        object.__setattr__(self, "loop_metadata", dict(self.loop_metadata))


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
    diagnostics: Mapping[str, Any]


def minimal_rs_t_coefficient(
    *,
    rs_volume_log: float = DEFAULT_RS_VOLUME_LOG,
    sin2_theta_w: float = DEFAULT_SIN2_THETA_W,
) -> float:
    """Return the volume-enhanced minimal-RS ``Delta T`` coefficient.

    SOURCE: CGHNP (Casagrande-Goertz-Haisch-Neubert-Pfoh, arXiv:0807.4937)
    Eq. (147), T = (pi v^2 / (2 cos^2 theta_W M_KK^2)) (L - 1/(2L)) in the
    GEOMETRIC-Lambda_IR convention.  This implementation keeps the leading L
    (a deliberate, conservative ~0.04% over-estimate that drops the -1/(2L)
    subleading term -> tighter floor) and applies the x_1^2 geometric->physical
    conversion below.  See docs/CUSTODIAL_PROVENANCE.md (item 1).
    """
    volume = _positive(rs_volume_log, "rs_volume_log")
    sin2 = _finite(sin2_theta_w, "sin2_theta_w")
    cos2 = 1.0 - sin2
    if cos2 <= 0.0:
        raise ValueError("sin2_theta_w must be less than one")
    # x_1^2 converts the geometric-Lambda_IR coefficient to the physical-M_KK
    # convention used by the proxy's shared (v/M_KK)^2 scale (PLAN §5/M2).
    return float(PHYSICAL_MKK_OVER_LAMBDA_IR_SQ * math.pi * volume / (2.0 * cos2))


def custodial_rs_plr_t_coefficient(
    *,
    rs_volume_log: float = DEFAULT_RS_VOLUME_LOG,
    sin2_theta_w: float = DEFAULT_SIN2_THETA_W,
) -> float:
    """Return the custodial-P_LR tree proxy ``Delta T`` coefficient.

    SOURCE: CGHNP (arXiv:0807.4937) Eq. (153), the custodial SU(2)_R result
    T = -(pi v^2 / (4 cos^2 theta_W M_KK^2)) (1/L) in the GEOMETRIC-Lambda_IR
    convention, which CGHNP attribute to Agashe-Delgado-May-Sundrum
    (hep-ph/0308036, [33] therein).  The 1/L suppression AND the sign flip
    (relative to minimal_rs_t_coefficient) are both genuinely in CGHNP (153);
    the x_1^2 geometric->physical conversion below matches the minimal case.
    See docs/CUSTODIAL_PROVENANCE.md (item 2).
    """
    volume = _positive(rs_volume_log, "rs_volume_log")
    sin2 = _finite(sin2_theta_w, "sin2_theta_w")
    cos2 = 1.0 - sin2
    if cos2 <= 0.0:
        raise ValueError("sin2_theta_w must be less than one")
    # x_1^2 converts the geometric-Lambda_IR coefficient to the physical-M_KK
    # convention used by the proxy's shared (v/M_KK)^2 scale (PLAN §5/M2).
    return float(-PHYSICAL_MKK_OVER_LAMBDA_IR_SQ * math.pi / (4.0 * cos2 * volume))


def rs_oblique_t_coefficient(
    *,
    ew_model: str = MINIMAL_RS_EW_MODEL,
    rs_volume_log: float = DEFAULT_RS_VOLUME_LOG,
    sin2_theta_w: float = DEFAULT_SIN2_THETA_W,
) -> float:
    """Return the EW-model-selected RS oblique ``Delta T`` coefficient."""
    model = _validate_ew_model(ew_model)
    # Both custodial models protect T via SU(2)_R (CGHNP Eq. 153); the difference
    # between them (P_LR protection of Z b_L) lives in rs_ew_couplings.py, not in
    # the oblique sector, so they share the same Delta T coefficient here.
    if model in (CUSTODIAL_RS_PLR_EW_MODEL, CUSTODIAL_RS_SU2R_EW_MODEL):
        return custodial_rs_plr_t_coefficient(
            rs_volume_log=rs_volume_log,
            sin2_theta_w=sin2_theta_w,
        )
    return minimal_rs_t_coefficient(
        rs_volume_log=rs_volume_log,
        sin2_theta_w=sin2_theta_w,
    )


def rs_minimal_oblique_proxy(
    *,
    m_kk_gev: float,
    s_coefficient: float,
    ew_model: str = MINIMAL_RS_EW_MODEL,
    higgs_vev_gev: float = DEFAULT_HIGGS_VEV_GEV,
    sin2_theta_w: float = DEFAULT_SIN2_THETA_W,
    rs_volume_log: float = DEFAULT_RS_VOLUME_LOG,
    delta_t_loop: float = 0.0,
    loop_metadata: Mapping[str, Any] | None = None,
) -> ObliqueSTUShift:
    """Return the documented RS ``S,T,U`` proxy at ``M_KK``."""
    model = _validate_ew_model(ew_model)
    mass = _positive(m_kk_gev, "m_kk_gev")
    vev = _positive(higgs_vev_gev, "higgs_vev_gev")
    coeff_s = _finite(s_coefficient, "s_coefficient")
    coeff_t = rs_oblique_t_coefficient(
        ew_model=model,
        rs_volume_log=rs_volume_log,
        sin2_theta_w=sin2_theta_w,
    )
    try:
        loop_delta = _finite(delta_t_loop, "delta_t_loop")
    except (TypeError, ValueError) as exc:
        raise ValueError("delta_t_loop must be finite") from exc
    loop_meta = {} if loop_metadata is None else dict(loop_metadata)
    if bool(loop_meta.get("top_partner_loop_numerics_included", False)) and not math.isfinite(loop_delta):
        raise ValueError("finite delta_t_loop is required when loop metadata says it is included")
    scale = float((vev / mass) ** 2)
    t_tree = float(coeff_t * scale)
    return ObliqueSTUShift(
        s=float(coeff_s * scale),
        t=float(t_tree + loop_delta),
        u=0.0,
        m_kk_gev=float(mass),
        s_coefficient=float(coeff_s),
        t_coefficient=float(coeff_t),
        higgs_vev_gev=float(vev),
        sin2_theta_w=float(sin2_theta_w),
        rs_volume_log=float(rs_volume_log),
        ew_model=model,
        t_tree=t_tree,
        delta_t_loop=loop_delta,
        loop_metadata=loop_meta,
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
    ew_model: str = MINIMAL_RS_EW_MODEL,
    higgs_vev_gev: float = DEFAULT_HIGGS_VEV_GEV,
    sin2_theta_w: float = DEFAULT_SIN2_THETA_W,
    rs_volume_log: float = DEFAULT_RS_VOLUME_LOG,
    delta_t_loop: float = 0.0,
    loop_metadata: Mapping[str, Any] | None = None,
) -> ObliqueSTComparison:
    """Evaluate the documented RS proxy against an ``S,T`` fit ellipse."""
    prediction = rs_minimal_oblique_proxy(
        m_kk_gev=m_kk_gev,
        s_coefficient=s_coefficient,
        ew_model=ew_model,
        higgs_vev_gev=higgs_vev_gev,
        sin2_theta_w=sin2_theta_w,
        rs_volume_log=rs_volume_log,
        delta_t_loop=delta_t_loop,
        loop_metadata=loop_metadata,
    )
    chi2, ratio, budget, passes = compare_oblique_st_to_fit(
        s=prediction.s,
        t=prediction.t,
        fit=fit,
    )
    diagnostics = {
        "likelihood_model": OBLIQUE_STU_LIKELIHOOD_V1,
        "matching_assumption": OBLIQUE_STU_RS_PROXY_V1,
        "ew_model": prediction.ew_model,
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
    if float(prediction.delta_t_loop) != 0.0 or prediction.loop_metadata:
        diagnostics.update(
            {
                "t_tree_prediction": float(
                    prediction.t if prediction.t_tree is None else prediction.t_tree
                ),
                "t_loop_prediction": float(prediction.delta_t_loop),
            }
        )
        diagnostics.update(dict(prediction.loop_metadata))
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


def _validate_ew_model(ew_model: str) -> str:
    model = str(ew_model)
    if model not in SUPPORTED_RS_EW_MODELS:
        raise ValueError(
            f"unsupported ew_model {model!r}; supported models are "
            f"{SUPPORTED_RS_EW_MODELS}"
        )
    return model


__all__ = [
    "CHI2_2DOF_95",
    "DEFAULT_HIGGS_VEV_GEV",
    "DEFAULT_RS_VOLUME_LOG",
    "DEFAULT_SIN2_THETA_W",
    "OBLIQUE_STU_LIKELIHOOD_V1",
    "OBLIQUE_STU_RS_PROXY_V1",
    "CUSTODIAL_RS_PLR_EW_MODEL",
    "ObliqueSTComparison",
    "ObliqueSTFit",
    "ObliqueSTUShift",
    "compare_oblique_st_to_fit",
    "custodial_rs_plr_t_coefficient",
    "evaluate_rs_oblique_proxy",
    "minimal_rs_t_coefficient",
    "rs_oblique_t_coefficient",
    "rs_minimal_oblique_proxy",
]
