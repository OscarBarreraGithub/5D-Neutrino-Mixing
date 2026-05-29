"""Adapter over :mod:`flavorConstraints.muToEGamma`.

This is the catalog boundary for charged-lepton dipole LFV constraints.
Constraint modules import this adapter only; the underlying implementation
remains isolated in ``flavorConstraints``.

NEEDS-HUMAN-PHYSICS
-------------------
The current ``ParameterPoint`` scaffold is quark-sector first and does not
carry a full lepton mass-basis RS coupling object.  Until that exists, this
adapter accepts either the repo's ``YukawaResult`` or an explicit
``MuToEGammaProxyInput`` made from caller-supplied lepton spurions
``Y_N_bar``, PMNS, and ``M_KK``.  The proxy path is documented and flagged in
diagnostics; it is not a substitute for a loop-level RS lepton-dipole match.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

import numpy as np

from flavorConstraints.muToEGamma import (
    check_mu_to_e_gamma as _check_mu_to_e_gamma,
    check_mu_to_e_gamma_raw as _check_mu_to_e_gamma_raw,
    coefficient_from_br_limit as _coefficient_from_br_limit,
)

__all__ = [
    "LEPTON_DIPOLE_PROXY_ASSUMPTION_V1",
    "MuToEGammaProxyInput",
    "MuToEGammaBranchingResult",
    "mu_to_e_gamma_proxy_input",
    "mu_to_e_gamma_coefficient_from_limit",
    "mu_to_e_gamma_from_lepton_input",
]

LEPTON_DIPOLE_PROXY_ASSUMPTION_V1 = (
    "NEEDS-HUMAN-PHYSICS: lepton-sector RS mass-basis couplings are not "
    "available on the quark-only ParameterPoint; proxy inputs must be supplied "
    "explicitly as Y_N_bar, PMNS, and M_KK, and the dipole normalization uses "
    "the repo's Perez-Randall NDA mu->e gamma convention."
)


@dataclass(frozen=True)
class MuToEGammaProxyInput:
    """Explicit proxy inputs for the repo's mu->e gamma NDA checker.

    ``y_n_bar`` are the rescaled neutrino Yukawa eigenvalues and ``pmns`` is the
    charged-lepton mass-basis mixing matrix.  This is intentionally thin: it
    supplies only the data needed by ``check_mu_to_e_gamma_raw`` and carries a
    source label for diagnostics.
    """

    y_n_bar: tuple[complex, ...]
    pmns: tuple[tuple[complex, ...], ...]
    m_kk_gev: float
    source: str = "caller-supplied lepton dipole proxy"


@dataclass(frozen=True)
class MuToEGammaBranchingResult:
    """Branching-fraction view of the lower-level dipole-bound result."""

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
    off_diagonal_12: complex
    product_matrix: tuple[tuple[complex, ...], ...]
    input_kind: str
    used_proxy: bool
    diagnostics: Mapping[str, Any]


def _positive_finite(value: Any, *, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be a positive finite number")
    return number


def _matrix_tuple(matrix: Any) -> tuple[tuple[complex, ...], ...]:
    array = np.asarray(matrix, dtype=complex)
    return tuple(tuple(complex(value) for value in row) for row in array)


def mu_to_e_gamma_proxy_input(
    y_n_bar: Any,
    pmns: Any,
    m_kk_gev: float,
    *,
    source: str = "caller-supplied lepton dipole proxy",
) -> MuToEGammaProxyInput:
    """Build a shape-checked proxy input for ``mu_to_e_gamma_from_lepton_input``."""

    y_array = np.asarray(y_n_bar, dtype=complex).reshape(-1)
    if y_array.shape != (3,):
        raise ValueError("y_n_bar must contain exactly three entries")
    pmns_array = np.asarray(pmns, dtype=complex)
    if pmns_array.shape != (3, 3):
        raise ValueError("pmns must be a 3x3 matrix")
    return MuToEGammaProxyInput(
        y_n_bar=tuple(complex(value) for value in y_array),
        pmns=_matrix_tuple(pmns_array),
        m_kk_gev=_positive_finite(m_kk_gev, name="m_kk_gev"),
        source=str(source),
    )


def mu_to_e_gamma_coefficient_from_limit(
    br_limit: float,
    prefactor_br: float,
) -> float:
    """Return the LFV coefficient ``C = sqrt(BR_limit / prefactor_br)``."""

    return float(
        _coefficient_from_br_limit(
            _positive_finite(br_limit, name="br_limit"),
            prefactor=_positive_finite(prefactor_br, name="prefactor_br"),
        )
    )


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
    br_limit: float,
    prefactor_br: float,
    c_lfv: float,
    reference_scale_gev: float,
    m_kk_gev: float,
    input_kind: str,
    used_proxy: bool,
    extra_diagnostics: Mapping[str, Any],
) -> MuToEGammaBranchingResult:
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
    diagnostics = {
        "core_passes": bool(core["passes"]),
        "branching_formula": (
            "BR_NP = prefactor_br * |(Y_N_bar Y_N_bar^dagger)_{e mu}|^2 "
            "* (reference_scale_gev / M_KK)^4"
        ),
        **dict(extra_diagnostics),
    }
    return MuToEGammaBranchingResult(
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
        off_diagonal_12=complex(core["off_diagonal_12"]),
        product_matrix=_matrix_tuple(core["product_matrix"]),
        input_kind=input_kind,
        used_proxy=used_proxy,
        diagnostics=diagnostics,
    )


def _is_yukawa_result_like(value: Any) -> bool:
    return hasattr(value, "Y_N_matrix") and hasattr(value, "params")


def _mapping_value(mapping: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    raise KeyError(f"none of keys {keys!r} present")


def _proxy_from_mapping(mapping: Mapping[str, Any]) -> MuToEGammaProxyInput:
    return mu_to_e_gamma_proxy_input(
        _mapping_value(mapping, "y_n_bar", "Y_N_bar"),
        _mapping_value(mapping, "pmns", "PMNS", "pmns_matrix"),
        _mapping_value(mapping, "m_kk_gev", "M_KK", "M_KK_gev"),
        source=str(mapping.get("source", "mapping lepton dipole proxy")),
    )


def _coerce_proxy(value: Any) -> MuToEGammaProxyInput:
    if isinstance(value, MuToEGammaProxyInput):
        return value
    if isinstance(value, Mapping):
        return _proxy_from_mapping(value)
    raise TypeError(
        "lepton input must be a YukawaResult-like object, a MuToEGammaProxyInput, "
        "or a mapping with y_n_bar, pmns, and m_kk_gev"
    )


def mu_to_e_gamma_from_lepton_input(
    lepton_input: Any,
    *,
    br_limit: float,
    prefactor_br: float,
    reference_scale_gev: float = 3000.0,
    m_kk_gev: float | None = None,
) -> MuToEGammaBranchingResult:
    """Evaluate ``BR(mu -> e gamma)`` from a YukawaResult or proxy input."""

    limit = _positive_finite(br_limit, name="br_limit")
    prefactor = _positive_finite(prefactor_br, name="prefactor_br")
    reference_scale = _positive_finite(
        reference_scale_gev,
        name="reference_scale_gev",
    )
    c_lfv = mu_to_e_gamma_coefficient_from_limit(limit, prefactor)

    if isinstance(lepton_input, Mapping) and "yukawa_result" in lepton_input:
        lepton_input = lepton_input["yukawa_result"]

    if _is_yukawa_result_like(lepton_input):
        m_kk = _m_kk_from_yukawa_result(lepton_input, m_kk_gev)
        core = _check_mu_to_e_gamma(
            lepton_input,
            C=c_lfv,
            reference_scale=reference_scale,
            M_KK_override=m_kk,
        )
        return _result_from_core(
            core,
            br_limit=limit,
            prefactor_br=prefactor,
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
    core = _check_mu_to_e_gamma_raw(
        np.asarray(proxy.y_n_bar, dtype=complex),
        np.asarray(proxy.pmns, dtype=complex),
        M_KK=m_kk,
        C=c_lfv,
        reference_scale=reference_scale,
    )
    return _result_from_core(
        core,
        br_limit=limit,
        prefactor_br=prefactor,
        c_lfv=c_lfv,
        reference_scale_gev=reference_scale,
        m_kk_gev=m_kk,
        input_kind="MuToEGammaProxyInput",
        used_proxy=True,
        extra_diagnostics={
            "proxy_source": proxy.source,
            "needs_human_physics": LEPTON_DIPOLE_PROXY_ASSUMPTION_V1,
        },
    )
