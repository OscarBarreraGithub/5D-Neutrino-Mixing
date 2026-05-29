"""Neutrino-trident ``sigma/sigma_SM`` adapter.

This adapter owns the small model-independent trident-response formula used by
L023.  It does not attempt a nuclear/flux calculation.  The observable is the
ratio to the SM prediction, so the rigorous SM limit is exactly one in this
normalization.

NEEDS-HUMAN-PHYSICS
-------------------
The current ``ParameterPoint`` does not carry the RS EW KK/Z/Z' spectrum or the
mass-basis ``nu_mu nu_mu mu mu`` neutral-current couplings needed for a full
trident recast.  The proxy path therefore accepts caller-supplied effective
coupling shifts, or a simple heavy-Z'-like contact estimate, and reports that
assumption in diagnostics.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping

__all__ = [
    "DEFAULT_TRIDENT_SIN2_THETA_W",
    "DEFAULT_TRIDENT_VEV_GEV",
    "TRIDENT_PARAMETRIZATION_CITATION",
    "TRIDENT_PROXY_ASSUMPTION_V1",
    "TridentProxyInput",
    "TridentRatioResult",
    "trident_proxy_input",
    "trident_zprime_proxy_input",
    "trident_sigma_ratio_from_effective_couplings",
    "trident_sigma_ratio_from_lepton_input",
]


DEFAULT_TRIDENT_SIN2_THETA_W = 0.23126
DEFAULT_TRIDENT_VEV_GEV = 246.0

TRIDENT_PARAMETRIZATION_CITATION = (
    "Altmannshofer, Gori, Pospelov, Yavin, Phys. Rev. Lett. 113 (2014) "
    "091801, arXiv:1406.2332; ratio approximation "
    "(C_V^2 + C_A^2)/(C_V_SM^2 + C_A_SM^2)"
)

TRIDENT_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: full RS EW KK/Z/Z' and lepton/neutrino "
    "neutral-current matching to nu_mu nu_mu mu mu is not available on "
    "ParameterPoint; v1 uses caller-supplied effective shifts in the "
    "CC-normalized trident C_V/C_A coefficients, or a heavy-Z'-like contact "
    "proxy Delta C_{V,A}=2 v^2 g_nu g_{V,A}/M^2."
)


@dataclass(frozen=True)
class TridentProxyInput:
    """Explicit effective-coupling proxy for muon-neutrino trident production.

    ``delta_c_vector`` and ``delta_c_axial`` shift the CC-normalized
    coefficients

        C_V^SM = 1 + 4 sin^2(theta_W),   C_A^SM = 1.

    A heavy vector mediator with ``M`` well above the trident momentum transfer
    may be mapped to these shifts with ``trident_zprime_proxy_input``.  That
    mapping is only a documented proxy for L023, not a full RS calculation.
    """

    delta_c_vector: float
    delta_c_axial: float = 0.0
    source: str = "caller-supplied neutrino-trident effective-coupling proxy"
    input_kind: str = "effective_shift"
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class TridentRatioResult:
    """Prediction for ``sigma(nu_mu N -> nu_mu N mu+ mu-) / sigma_SM``."""

    sigma_ratio: float
    sm_sigma_ratio: float
    c_vector_sm: float
    c_axial_sm: float
    c_vector_total: float
    c_axial_total: float
    delta_c_vector: float
    delta_c_axial: float
    sin2_theta_w: float
    diagnostics: Mapping[str, Any]


def _finite_float(value: Any, *, name: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return number


def _positive_finite(value: Any, *, name: str) -> float:
    number = _finite_float(value, name=name)
    if number <= 0.0:
        raise ValueError(f"{name} must be positive, got {value!r}")
    return number


def trident_proxy_input(
    delta_c_vector: float,
    delta_c_axial: float = 0.0,
    *,
    source: str = "caller-supplied neutrino-trident effective-coupling proxy",
) -> TridentProxyInput:
    """Build a shape-checked direct effective-shift proxy."""

    return TridentProxyInput(
        delta_c_vector=_finite_float(delta_c_vector, name="delta_c_vector"),
        delta_c_axial=_finite_float(delta_c_axial, name="delta_c_axial"),
        source=str(source),
        input_kind="effective_shift",
        diagnostics={"proxy_source": str(source)},
    )


def trident_zprime_proxy_input(
    *,
    g_nu_mu: float,
    mediator_mass_gev: float,
    g_mu_vector: float | None = None,
    g_mu_axial: float | None = None,
    g_mu_left: float | None = None,
    g_mu_right: float | None = None,
    vev_gev: float = DEFAULT_TRIDENT_VEV_GEV,
    source: str = "caller-supplied heavy-Zprime neutrino-trident proxy",
) -> TridentProxyInput:
    """Map a heavy-Z'-like contact interaction to trident coefficient shifts.

    ``g_mu_vector`` and ``g_mu_axial`` are defined as
    ``(g_mu_left + g_mu_right)/2`` and ``(g_mu_left - g_mu_right)/2`` when the
    chiral couplings are supplied.  The normalization follows the common
    ``L_mu-L_tau`` trident estimate, where a vectorlike coupling gives
    ``Delta C_V = 2 v^2 g'^2 / M^2``.
    """

    if g_mu_vector is None:
        if g_mu_left is None or g_mu_right is None:
            raise ValueError(
                "provide g_mu_vector or both g_mu_left and g_mu_right"
            )
        g_mu_vector = 0.5 * (
            _finite_float(g_mu_left, name="g_mu_left")
            + _finite_float(g_mu_right, name="g_mu_right")
        )
    if g_mu_axial is None:
        if g_mu_left is None or g_mu_right is None:
            g_mu_axial = 0.0
        else:
            g_mu_axial = 0.5 * (
                _finite_float(g_mu_left, name="g_mu_left")
                - _finite_float(g_mu_right, name="g_mu_right")
            )

    g_nu = _finite_float(g_nu_mu, name="g_nu_mu")
    g_v = _finite_float(g_mu_vector, name="g_mu_vector")
    g_a = _finite_float(g_mu_axial, name="g_mu_axial")
    mass = _positive_finite(mediator_mass_gev, name="mediator_mass_gev")
    vev = _positive_finite(vev_gev, name="vev_gev")
    prefactor = 2.0 * vev * vev / (mass * mass)

    return TridentProxyInput(
        delta_c_vector=float(prefactor * g_nu * g_v),
        delta_c_axial=float(prefactor * g_nu * g_a),
        source=str(source),
        input_kind="heavy_zprime_contact",
        diagnostics={
            "proxy_source": str(source),
            "g_nu_mu": float(g_nu),
            "g_mu_vector": float(g_v),
            "g_mu_axial": float(g_a),
            "mediator_mass_gev": float(mass),
            "vev_gev": float(vev),
            "zprime_contact_prefactor": float(prefactor),
        },
    )


def trident_sigma_ratio_from_effective_couplings(
    delta_c_vector: float = 0.0,
    delta_c_axial: float = 0.0,
    *,
    sin2_theta_w: float = DEFAULT_TRIDENT_SIN2_THETA_W,
    extra_diagnostics: Mapping[str, Any] | None = None,
) -> TridentRatioResult:
    """Return the trident cross-section ratio from effective coupling shifts."""

    sin2 = _finite_float(sin2_theta_w, name="sin2_theta_w")
    if not 0.0 < sin2 < 1.0:
        raise ValueError("sin2_theta_w must lie between 0 and 1")
    delta_v = _finite_float(delta_c_vector, name="delta_c_vector")
    delta_a = _finite_float(delta_c_axial, name="delta_c_axial")

    c_vector_sm = 1.0 + 4.0 * sin2
    c_axial_sm = 1.0
    denominator = c_vector_sm * c_vector_sm + c_axial_sm * c_axial_sm
    c_vector_total = c_vector_sm + delta_v
    c_axial_total = c_axial_sm + delta_a
    numerator = c_vector_total * c_vector_total + c_axial_total * c_axial_total
    ratio = numerator / denominator

    diagnostics = {
        "model_label": "neutrino_trident_effective_cv_ca_proxy_v1",
        "ratio_formula": (
            "sigma/sigma_SM = ((C_V_SM + Delta C_V)^2 + "
            "(C_A_SM + Delta C_A)^2) / (C_V_SM^2 + C_A_SM^2)"
        ),
        "parametrization_citation": TRIDENT_PARAMETRIZATION_CITATION,
        "needs_human_physics": TRIDENT_PROXY_ASSUMPTION_V1,
    }
    if extra_diagnostics is not None:
        diagnostics.update(dict(extra_diagnostics))

    return TridentRatioResult(
        sigma_ratio=float(ratio),
        sm_sigma_ratio=1.0,
        c_vector_sm=float(c_vector_sm),
        c_axial_sm=float(c_axial_sm),
        c_vector_total=float(c_vector_total),
        c_axial_total=float(c_axial_total),
        delta_c_vector=float(delta_v),
        delta_c_axial=float(delta_a),
        sin2_theta_w=float(sin2),
        diagnostics=diagnostics,
    )


def _mapping_value(mapping: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    raise KeyError(f"none of keys {keys!r} present")


def _has_any(mapping: Mapping[str, Any], keys: tuple[str, ...]) -> bool:
    return any(key in mapping for key in keys)


def _proxy_from_mapping(
    mapping: Mapping[str, Any],
    *,
    m_kk_gev: float | None,
) -> TridentProxyInput:
    if "trident" in mapping:
        nested = mapping["trident"]
        if not isinstance(nested, Mapping):
            raise TypeError("lepton input 'trident' entry must be a mapping")
        mapping = nested

    direct_keys = ("trident_delta_c_vector", "delta_c_vector", "delta_cv")
    if _has_any(mapping, direct_keys):
        return trident_proxy_input(
            _mapping_value(mapping, *direct_keys),
            _mapping_value(
                mapping,
                "trident_delta_c_axial",
                "delta_c_axial",
                "delta_ca",
            )
            if _has_any(
                mapping,
                ("trident_delta_c_axial", "delta_c_axial", "delta_ca"),
            )
            else 0.0,
            source=str(
                mapping.get(
                    "source",
                    "mapping neutrino-trident effective-coupling proxy",
                )
            ),
        )

    mass = (
        _positive_finite(m_kk_gev, name="m_kk_gev")
        if m_kk_gev is not None
        else _positive_finite(
            _mapping_value(mapping, "mediator_mass_gev", "m_kk_gev", "M_KK"),
            name="mediator_mass_gev",
        )
    )
    return trident_zprime_proxy_input(
        g_nu_mu=_mapping_value(mapping, "g_nu_mu", "g_nu", "gprime_nu_mu"),
        g_mu_vector=mapping.get("g_mu_vector", mapping.get("g_mu_v")),
        g_mu_axial=mapping.get("g_mu_axial", mapping.get("g_mu_a")),
        g_mu_left=mapping.get("g_mu_left", mapping.get("g_mu_l")),
        g_mu_right=mapping.get("g_mu_right", mapping.get("g_mu_r")),
        mediator_mass_gev=mass,
        source=str(mapping.get("source", "mapping heavy-Zprime trident proxy")),
    )


def _coerce_proxy(
    lepton_input: Any,
    *,
    m_kk_gev: float | None,
) -> TridentProxyInput:
    if isinstance(lepton_input, TridentProxyInput):
        return lepton_input
    if isinstance(lepton_input, Mapping):
        return _proxy_from_mapping(lepton_input, m_kk_gev=m_kk_gev)
    raise TypeError(
        "lepton input must be a TridentProxyInput or a mapping with direct "
        "trident_delta_c_* shifts, or heavy-Zprime g_nu/g_mu/M entries"
    )


def trident_sigma_ratio_from_lepton_input(
    lepton_input: Any,
    *,
    m_kk_gev: float | None = None,
    sin2_theta_w: float = DEFAULT_TRIDENT_SIN2_THETA_W,
) -> TridentRatioResult:
    """Evaluate ``sigma/sigma_SM`` from a caller-supplied trident proxy."""

    proxy = _coerce_proxy(lepton_input, m_kk_gev=m_kk_gev)
    diagnostics = {
        "input_kind": proxy.input_kind,
        "proxy_source": proxy.source,
        "kk_ew_mass_override_used": m_kk_gev is not None,
        **dict(proxy.diagnostics),
    }
    return trident_sigma_ratio_from_effective_couplings(
        proxy.delta_c_vector,
        proxy.delta_c_axial,
        sin2_theta_w=sin2_theta_w,
        extra_diagnostics=diagnostics,
    )
