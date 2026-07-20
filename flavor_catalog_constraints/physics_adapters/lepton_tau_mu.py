"""Tau-muon radiative charged-LFV adapter.

This adapter reuses the L001 lepton-dipole machinery for ``mu -> e gamma`` and
pins it to ``tau -> mu gamma`` by permuting the charged-lepton rows before
calling :func:`flavorConstraints.muToEGamma.check_mu_to_e_gamma_raw`.  The core
still computes its documented ``(0, 1)`` entry; after the row permutation this
is the canonical ``(mu, tau)`` dipole entry of the original matrix.

NEEDS-HUMAN-PHYSICS: a rigorous RS prediction for ``tau -> mu gamma`` requires
lepton-sector mass-basis couplings and loop-level radiative matching that are
not present on ``ParameterPoint``.  The accepted inputs are therefore explicit
proxies in the same Perez-Randall NDA normalization used by L001.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from flavorConstraints.muToEGamma import (
    check_mu_to_e_gamma as _check_mu_to_e_gamma,
)
from flavorConstraints.muToEGamma import (
    check_mu_to_e_gamma_raw as _check_mu_to_e_gamma_raw,
)
from quarkConstraints.lfv_three_body import (
    PDG_TAU_LEPTONIC_BRANCHING_FRACTION_SOURCE,
    TAU_TO_MU_NUNU_BRANCHING_FRACTION,
)

from .lepton import (
    LEPTON_DIPOLE_PROXY_ASSUMPTION_V1,
    MuToEGammaProxyInput,
    mu_to_e_gamma_coefficient_from_limit,
    mu_to_e_gamma_proxy_input,
)

__all__ = [
    "TAU_TO_MU_GAMMA_PROXY_V1",
    "TauToMuGammaProxyInput",
    "TauToMuGammaBranchingResult",
    "tau_to_mu_gamma_proxy_input",
    "tau_to_mu_gamma_from_lepton_input",
]

TAU_TO_MU_GAMMA_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: tau->mu gamma is evaluated as a documented proxy "
    "using the same Perez-Randall NDA lepton-dipole normalization as L001. "
    "The adapter maps the L001 core's (e,mu) matrix element to the tau->mu "
    "(mu,tau) element by a charged-lepton row permutation; full tau radiative "
    "RS loop matching is not implemented."
)

_ROW_PERMUTATION_MU_TAU_E = (1, 2, 0)
_FINAL_FLAVOR = "mu"
_INITIAL_FLAVOR = "tau"
_FINAL_INDEX = 1
_INITIAL_INDEX = 2
_INITIAL_FLAVOR_KEYS = ("initial_flavor", "parent_flavor")
_FINAL_FLAVOR_KEYS = ("final_flavor", "daughter_flavor")
_M_KK_KEYS = ("m_kk_gev", "M_KK_gev", "M_KK")
_CHARGED_LEPTON_ALIASES = {
    "electron": "e",
    "e-": "e",
    "e+": "e",
    "muon": "mu",
    "mu-": "mu",
    "mu+": "mu",
    "tau-": "tau",
    "tau+": "tau",
}
_CHARGED_LEPTONS = frozenset({"e", "mu", "tau"})


@dataclass(frozen=True)
class TauToMuGammaProxyInput:
    """Explicit proxy inputs for the tau->mu radiative LFV dipole.

    ``y_n_bar`` are rescaled neutrino Yukawa eigenvalues and ``pmns`` is the
    charged-lepton mass-basis mixing matrix in the conventional row order
    ``(e, mu, tau)``.
    """

    y_n_bar: tuple[complex, ...]
    pmns: tuple[tuple[complex, ...], ...]
    m_kk_gev: float
    source: str = "caller-supplied tau->mu gamma dipole proxy"


@dataclass(frozen=True)
class TauToMuGammaBranchingResult:
    """Branching-fraction view of the tau->mu dipole-bound proxy."""

    branching_fraction: float
    sm_branching_fraction: float
    br_limit: float
    passes: bool
    ratio_to_limit: float
    dipole_lhs: float
    dipole_rhs: float
    dipole_ratio_to_bound: float
    c_lfv: float
    prefactor_br: float
    m_kk_gev: float
    reference_scale_gev: float
    off_diagonal_23: complex
    off_diagonal_32: complex
    product_matrix: tuple[tuple[complex, ...], ...]
    row_permutation: tuple[int, int, int]
    input_kind: str
    used_proxy: bool
    diagnostics: Mapping[str, Any]


@dataclass(frozen=True)
class _YukawaResultRowView:
    """Minimal row-permuted view accepted by the mu->e gamma core."""

    Y_N_matrix: np.ndarray
    params: Mapping[str, Any]


def tau_to_mu_gamma_proxy_input(
    y_n_bar: Any,
    pmns: Any,
    m_kk_gev: float,
    *,
    source: str = "caller-supplied tau->mu gamma dipole proxy",
) -> TauToMuGammaProxyInput:
    """Build a shape-checked proxy input for ``tau_to_mu_gamma_from_lepton_input``."""

    validated = mu_to_e_gamma_proxy_input(
        y_n_bar,
        pmns,
        m_kk_gev,
        source=source,
    )
    return TauToMuGammaProxyInput(
        y_n_bar=validated.y_n_bar,
        pmns=validated.pmns,
        m_kk_gev=validated.m_kk_gev,
        source=str(source),
    )


def tau_to_mu_gamma_from_lepton_input(
    lepton_input: Any,
    *,
    br_limit: float,
    prefactor_br: float,
    reference_scale_gev: float = 3000.0,
    m_kk_gev: float | None = None,
) -> TauToMuGammaBranchingResult:
    """Evaluate proxy ``BR(tau -> mu gamma)`` from a YukawaResult or proxy input."""

    _assert_tau_to_mu_flavors(lepton_input)
    limit = _positive_finite(br_limit, name="br_limit")
    raw_prefactor = _positive_finite(prefactor_br, name="prefactor_br")
    # PDG B(tau->mu nu nubar) converts the reused muon-normalized NDA width to a tau BR.
    prefactor = float(raw_prefactor * TAU_TO_MU_NUNU_BRANCHING_FRACTION)
    reference_scale = _positive_finite(
        reference_scale_gev,
        name="reference_scale_gev",
    )
    c_lfv = mu_to_e_gamma_coefficient_from_limit(limit, prefactor)

    if isinstance(lepton_input, Mapping) and "yukawa_result" in lepton_input:
        lepton_input = lepton_input["yukawa_result"]

    if _is_yukawa_result_like(lepton_input):
        m_kk = _m_kk_from_yukawa_result(lepton_input, m_kk_gev)
        y_matrix = np.asarray(lepton_input.Y_N_matrix, dtype=complex)
        if y_matrix.shape[0] != 3:
            raise ValueError("yukawa_result.Y_N_matrix must have three charged-lepton rows")
        params = lepton_input.params
        k = float(params["k"])
        product = (2.0 * k * y_matrix) @ (2.0 * k * y_matrix).conj().T
        core = _check_mu_to_e_gamma(
            _YukawaResultRowView(
                Y_N_matrix=y_matrix[list(_ROW_PERMUTATION_MU_TAU_E), :],
                params=params,
            ),
            C=c_lfv,
            reference_scale=reference_scale,
            M_KK_override=m_kk,
        )
        return _result_from_core(
            core,
            product_matrix=product,
            br_limit=limit,
            prefactor_br=prefactor,
            raw_prefactor_br=raw_prefactor,
            c_lfv=c_lfv,
            reference_scale_gev=reference_scale,
            m_kk_gev=m_kk,
            input_kind="YukawaResult",
            used_proxy=False,
            extra_diagnostics={"proxy_source": None},
        )

    proxy = _coerce_proxy(lepton_input)
    m_kk = _positive_finite(
        proxy.m_kk_gev if m_kk_gev is None else m_kk_gev,
        name="m_kk_gev",
    )
    y_n_bar = np.asarray(proxy.y_n_bar, dtype=complex)
    pmns = np.asarray(proxy.pmns, dtype=complex)
    y_matrix = pmns @ np.diag(y_n_bar)
    product = y_matrix @ y_matrix.conj().T
    core = _check_mu_to_e_gamma_raw(
        y_n_bar,
        pmns[list(_ROW_PERMUTATION_MU_TAU_E), :],
        M_KK=m_kk,
        C=c_lfv,
        reference_scale=reference_scale,
    )
    return _result_from_core(
        core,
        product_matrix=product,
        br_limit=limit,
        prefactor_br=prefactor,
        raw_prefactor_br=raw_prefactor,
        c_lfv=c_lfv,
        reference_scale_gev=reference_scale,
        m_kk_gev=m_kk,
        input_kind="TauToMuGammaProxyInput",
        used_proxy=True,
        extra_diagnostics={
            "proxy_source": proxy.source,
            "needs_human_physics": TAU_TO_MU_GAMMA_PROXY_V1,
            "shared_lepton_proxy_assumption": LEPTON_DIPOLE_PROXY_ASSUMPTION_V1,
        },
    )


def _positive_finite(value: Any, *, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be a positive finite number")
    return number


def _matrix_tuple(matrix: Any) -> tuple[tuple[complex, ...], ...]:
    array = np.asarray(matrix, dtype=complex)
    return tuple(tuple(complex(value) for value in row) for row in array)


def _branching_fraction(
    *,
    dipole_lhs: float,
    prefactor_br: float,
    reference_scale_gev: float,
    m_kk_gev: float,
) -> float:
    return float(prefactor_br * dipole_lhs**2 * (reference_scale_gev / m_kk_gev) ** 4)


def _result_from_core(
    core: Mapping[str, Any],
    *,
    product_matrix: np.ndarray,
    br_limit: float,
    prefactor_br: float,
    raw_prefactor_br: float,
    c_lfv: float,
    reference_scale_gev: float,
    m_kk_gev: float,
    input_kind: str,
    used_proxy: bool,
    extra_diagnostics: Mapping[str, Any],
) -> TauToMuGammaBranchingResult:
    lhs = float(core["lhs"])
    rhs = float(core["rhs"])
    dipole_ratio = float(core["ratio"])
    branching_fraction = _branching_fraction(
        dipole_lhs=lhs,
        prefactor_br=prefactor_br,
        reference_scale_gev=reference_scale_gev,
        m_kk_gev=m_kk_gev,
    )
    ratio_to_limit = float(branching_fraction / br_limit)
    off_23 = complex(product_matrix[_FINAL_INDEX, _INITIAL_INDEX])
    off_32 = complex(product_matrix[_INITIAL_INDEX, _FINAL_INDEX])
    diagnostics = {
        "core_passes": bool(core["passes"]),
        "reused_physics_module": "flavorConstraints.muToEGamma",
        "reused_adapter_validation": "flavor_catalog_constraints.physics_adapters.lepton",
        "branching_formula": (
            "BR_NP = muon_normalized_prefactor_br * B(tau->mu nu nubar) * "
            "|(Y_N_bar Y_N_bar^dagger)_{mu tau}|^2 "
            "* (reference_scale_gev / M_KK)^4"
        ),
        "muon_normalized_prefactor_br": float(raw_prefactor_br),
        "tau_leptonic_branching_fraction": float(
            TAU_TO_MU_NUNU_BRANCHING_FRACTION
        ),
        "tau_leptonic_branching_fraction_source": (
            PDG_TAU_LEPTONIC_BRANCHING_FRACTION_SOURCE
        ),
        "flavor_basis": "(e, mu, tau)",
        "final_flavor": _FINAL_FLAVOR,
        "initial_flavor": _INITIAL_FLAVOR,
        "final_initial_indices_zero_based": (_FINAL_INDEX, _INITIAL_INDEX),
        "row_permutation_for_mu_to_e_core": _ROW_PERMUTATION_MU_TAU_E,
        "core_off_diagonal_12_after_permutation": complex(core["off_diagonal_12"]),
        **dict(extra_diagnostics),
    }
    return TauToMuGammaBranchingResult(
        branching_fraction=branching_fraction,
        sm_branching_fraction=0.0,
        br_limit=float(br_limit),
        passes=bool(ratio_to_limit <= 1.0),
        ratio_to_limit=ratio_to_limit,
        dipole_lhs=lhs,
        dipole_rhs=rhs,
        dipole_ratio_to_bound=dipole_ratio,
        c_lfv=float(c_lfv),
        prefactor_br=float(prefactor_br),
        m_kk_gev=float(m_kk_gev),
        reference_scale_gev=float(reference_scale_gev),
        off_diagonal_23=off_23,
        off_diagonal_32=off_32,
        product_matrix=_matrix_tuple(product_matrix),
        row_permutation=_ROW_PERMUTATION_MU_TAU_E,
        input_kind=input_kind,
        used_proxy=used_proxy,
        diagnostics=diagnostics,
    )


def _is_yukawa_result_like(value: Any) -> bool:
    return hasattr(value, "Y_N_matrix") and hasattr(value, "params")


def _m_kk_from_yukawa_result(yukawa_result: Any, override: float | None) -> float:
    if override is not None:
        return _positive_finite(override, name="m_kk_gev")
    try:
        params = yukawa_result.params
        return _positive_finite(
            params.get("M_KK", params["Lambda_IR"]),
            name="yukawa_result M_KK/Lambda_IR",
        )
    except Exception as exc:
        raise ValueError(
            "yukawa_result must include params with 'Lambda_IR' or 'M_KK'"
        ) from exc


def _mapping_value(mapping: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    raise KeyError(f"none of keys {keys!r} present")


def _proxy_from_mapping(mapping: Mapping[str, Any]) -> TauToMuGammaProxyInput:
    return tau_to_mu_gamma_proxy_input(
        _mapping_value(mapping, "y_n_bar", "Y_N_bar"),
        _mapping_value(mapping, "pmns", "PMNS", "pmns_matrix"),
        _mapping_value(mapping, *_M_KK_KEYS),
        source=str(mapping.get("source", "mapping tau->mu gamma dipole proxy")),
    )


def _coerce_proxy(value: Any) -> TauToMuGammaProxyInput:
    if isinstance(value, TauToMuGammaProxyInput):
        return value
    if isinstance(value, MuToEGammaProxyInput):
        return TauToMuGammaProxyInput(
            y_n_bar=value.y_n_bar,
            pmns=value.pmns,
            m_kk_gev=value.m_kk_gev,
            source=value.source,
        )
    if isinstance(value, Mapping):
        return _proxy_from_mapping(value)
    raise TypeError(
        "lepton input must be a YukawaResult-like object, a "
        "TauToMuGammaProxyInput, or a mapping with y_n_bar, pmns, and m_kk_gev"
    )


def _assert_tau_to_mu_flavors(value: Any) -> None:
    initial = _first_explicit_flavor(value, _INITIAL_FLAVOR_KEYS)
    final = _first_explicit_flavor(value, _FINAL_FLAVOR_KEYS)
    if initial is not None and _canonical_flavor(initial) != _INITIAL_FLAVOR:
        raise ValueError(
            "tau->mu gamma adapter is pinned to initial_flavor='tau'; "
            f"got {initial!r}"
        )
    if final is not None and _canonical_flavor(final) != _FINAL_FLAVOR:
        raise ValueError(
            "tau->mu gamma adapter is pinned to final_flavor='mu'; "
            f"got {final!r}"
        )


def _first_explicit_flavor(value: Any, keys: tuple[str, ...]) -> Any:
    if isinstance(value, Mapping):
        for key in keys:
            if key in value:
                return value[key]
    else:
        for key in keys:
            if hasattr(value, key):
                return getattr(value, key)
    return None


def _canonical_flavor(value: Any) -> str:
    flavor = str(value).strip().lower()
    flavor = _CHARGED_LEPTON_ALIASES.get(flavor, flavor)
    if flavor not in _CHARGED_LEPTONS:
        raise ValueError(f"unknown charged-lepton flavor {value!r}")
    return flavor
