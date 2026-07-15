"""Adapter over :mod:`flavorConstraints.muToEGamma`.

This is the catalog boundary for charged-lepton dipole LFV constraints.
Constraint modules import this adapter only; the underlying implementation
remains isolated in ``flavorConstraints``.

The locked L001 path uses :class:`LMFVLeptonParameters`, a read-only carrier for
the Perez-Randall LMFV NDA spurion.  Legacy caller-supplied proxy inputs remain
supported and are still flagged in diagnostics.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

import numpy as np

from flavorConstraints.muToEGamma import (
    PEREZ_RANDALL_LFV_M_KK_CONVENTION,
    check_mu_to_e_gamma as _check_mu_to_e_gamma,
    check_mu_to_e_gamma_raw as _check_mu_to_e_gamma_raw,
    coefficient_from_br_limit as _coefficient_from_br_limit,
)

__all__ = [
    "LEPTON_DIPOLE_PROXY_ASSUMPTION_V1",
    "LMFVLeptonParameters",
    "MuToEGammaProxyInput",
    "MuToEGammaBranchingResult",
    "lmfv_lepton_parameters_from_yukawa_result",
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


LMFV_LEPTON_PARAMETERS_SOURCE_V1 = "yukawa.compute_yukawas.YukawaResult"
LMFV_LEPTON_PARAMETERS_MATCHING_STATUS_V1 = "locked_lmfv_nda_carrier"
_LMFV_FORBIDDEN_STATUS_TERMS = ("partial", "deferred", "proxy", "recast")


def _finite_array(array: np.ndarray, *, name: str) -> None:
    if not np.all(np.isfinite(array.real)) or not np.all(np.isfinite(array.imag)):
        raise ValueError(f"{name} must contain only finite values")


def _readonly_complex_array(value: Any, *, name: str, shape: tuple[int, ...]) -> np.ndarray:
    array = np.array(value, dtype=np.complex128, copy=True)
    if array.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    _finite_array(array, name=name)
    array.setflags(write=False)
    return array


def _readonly_real_array(value: Any, *, name: str, shape: tuple[int, ...]) -> np.ndarray:
    array = np.array(value, dtype=np.float64, copy=True)
    if array.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    array.setflags(write=False)
    return array


def _broadcast_real_triplet(value: Any, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if array.shape == ():
        return np.full(3, float(array), dtype=np.float64)
    if array.shape != (3,):
        raise ValueError(f"{name} must be scalar or have shape (3,)")
    return np.array(array, dtype=np.float64, copy=True)


def _finite_float(value: Any, *, name: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


@dataclass(frozen=True)
class LMFVLeptonParameters:
    """Read-only Perez-Randall LMFV NDA carrier for ``mu -> e gamma``."""

    Y_N: np.ndarray
    Y_N_bar: np.ndarray
    Y_N_matrix: np.ndarray
    Y_N_bar_matrix: np.ndarray
    pmns: np.ndarray
    lmfv_spurion: np.ndarray
    M_KK_gev: float
    M_N_gev: float
    c_L: np.ndarray
    c_E: np.ndarray
    c_N: np.ndarray
    v_gev: float
    k_gev: float
    epsilon: float
    Lambda_IR_gev: float
    ordering: str
    majorana_alpha: float
    majorana_beta: float
    max_abs_ybar: float
    source: str = LMFV_LEPTON_PARAMETERS_SOURCE_V1
    matching_status: str = LMFV_LEPTON_PARAMETERS_MATCHING_STATUS_V1

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "Y_N",
            _readonly_complex_array(self.Y_N, name="Y_N", shape=(3,)),
        )
        object.__setattr__(
            self,
            "Y_N_bar",
            _readonly_complex_array(self.Y_N_bar, name="Y_N_bar", shape=(3,)),
        )
        for name in ("Y_N_matrix", "Y_N_bar_matrix", "pmns", "lmfv_spurion"):
            object.__setattr__(
                self,
                name,
                _readonly_complex_array(getattr(self, name), name=name, shape=(3, 3)),
            )
        for name in ("c_L", "c_E", "c_N"):
            object.__setattr__(
                self,
                name,
                _readonly_real_array(getattr(self, name), name=name, shape=(3,)),
            )

        for name in ("M_KK_gev", "M_N_gev", "v_gev", "k_gev", "Lambda_IR_gev"):
            value = _positive_finite(getattr(self, name), name=name)
            object.__setattr__(self, name, value)
        object.__setattr__(self, "epsilon", _finite_float(self.epsilon, name="epsilon"))
        object.__setattr__(
            self,
            "majorana_alpha",
            _finite_float(self.majorana_alpha, name="majorana_alpha"),
        )
        object.__setattr__(
            self,
            "majorana_beta",
            _finite_float(self.majorana_beta, name="majorana_beta"),
        )
        max_abs_ybar = _finite_float(self.max_abs_ybar, name="max_abs_ybar")
        if max_abs_ybar < 0.0:
            raise ValueError("max_abs_ybar must be non-negative")
        object.__setattr__(self, "max_abs_ybar", max_abs_ybar)

        ordering = str(self.ordering)
        if not ordering:
            raise ValueError("ordering must be non-empty")
        object.__setattr__(self, "ordering", ordering)
        object.__setattr__(self, "source", str(self.source))
        status = str(self.matching_status)
        lowered_status = status.lower()
        if any(term in lowered_status for term in _LMFV_FORBIDDEN_STATUS_TERMS):
            raise ValueError("matching_status contains a proxy/partial/deferred marker")
        object.__setattr__(self, "matching_status", status)

        identity = np.eye(3, dtype=np.complex128)
        if not np.allclose(
            self.pmns.conjugate().T @ self.pmns,
            identity,
            rtol=0.0,
            atol=1.0e-8,
        ):
            raise ValueError("pmns must be unitary")
        if not np.allclose(
            self.Y_N_bar,
            2.0 * self.k_gev * self.Y_N,
            rtol=1.0e-9,
            atol=1.0e-24,
        ):
            raise ValueError("Y_N_bar must equal 2*k_gev*Y_N")
        expected_ybar_matrix = 2.0 * self.k_gev * self.Y_N_matrix
        expected_pmns_matrix = self.pmns @ np.diag(self.Y_N_bar)
        if not np.allclose(
            self.Y_N_bar_matrix,
            expected_ybar_matrix,
            rtol=1.0e-9,
            atol=1.0e-12,
        ) or not np.allclose(
            self.Y_N_bar_matrix,
            expected_pmns_matrix,
            rtol=1.0e-9,
            atol=1.0e-12,
        ):
            raise ValueError(
                "Y_N_bar_matrix must equal both 2*k_gev*Y_N_matrix and "
                "pmns@diag(Y_N_bar)"
            )
        expected_spurion = self.Y_N_bar_matrix @ self.Y_N_bar_matrix.conjugate().T
        if not np.allclose(
            self.lmfv_spurion,
            expected_spurion,
            rtol=1.0e-9,
            atol=1.0e-12,
        ):
            raise ValueError("lmfv_spurion must equal Y_N_bar_matrix@Y_N_bar_matrix^dagger")

    @property
    def y_n_bar(self) -> np.ndarray:
        """Alias matching legacy proxy naming."""

        return self.Y_N_bar

    @property
    def pmns_matrix(self) -> np.ndarray:
        """Alias matching legacy proxy naming."""

        return self.pmns


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


def lmfv_lepton_parameters_from_yukawa_result(
    yukawa_result: Any,
    *,
    m_kk_gev: float,
    source: str = LMFV_LEPTON_PARAMETERS_SOURCE_V1,
) -> LMFVLeptonParameters:
    """Build the locked LMFV carrier from a repo ``YukawaResult``."""

    from neutrinos.neutrinoValues import get_pmns

    try:
        params = dict(yukawa_result.params)
    except Exception as exc:
        raise ValueError("yukawa_result must expose a params mapping") from exc

    k_gev = _positive_finite(params["k"], name="k")
    y_n = _readonly_complex_array(getattr(yukawa_result, "Y_N"), name="Y_N", shape=(3,))
    y_n_bar = _readonly_complex_array(
        getattr(yukawa_result, "Y_N_bar"),
        name="Y_N_bar",
        shape=(3,),
    )
    y_n_matrix = _readonly_complex_array(
        getattr(yukawa_result, "Y_N_matrix"),
        name="Y_N_matrix",
        shape=(3, 3),
    )
    ordering = str(params["ordering"])
    pmns = get_pmns(
        ordering,
        float(params["majorana_alpha"]),
        float(params["majorana_beta"]),
    )
    y_n_bar_matrix = 2.0 * k_gev * y_n_matrix
    max_abs_ybar = max(
        float(np.max(np.abs(getattr(yukawa_result, "Y_N_bar")))),
        float(np.max(np.abs(getattr(yukawa_result, "Y_E_bar", [0.0])))),
    )
    return LMFVLeptonParameters(
        Y_N=y_n,
        Y_N_bar=y_n_bar,
        Y_N_matrix=y_n_matrix,
        Y_N_bar_matrix=y_n_bar_matrix,
        pmns=pmns,
        lmfv_spurion=y_n_bar_matrix @ y_n_bar_matrix.conjugate().T,
        M_KK_gev=_positive_finite(m_kk_gev, name="m_kk_gev"),
        M_N_gev=_positive_finite(params["M_N"], name="M_N"),
        c_L=_broadcast_real_triplet(params["c_L"], name="c_L"),
        c_E=_readonly_real_array(params["c_E"], name="c_E", shape=(3,)),
        c_N=_broadcast_real_triplet(params["c_N"], name="c_N"),
        v_gev=_positive_finite(params["v"], name="v"),
        k_gev=k_gev,
        epsilon=_finite_float(getattr(yukawa_result, "epsilon"), name="epsilon"),
        Lambda_IR_gev=_positive_finite(params["Lambda_IR"], name="Lambda_IR"),
        ordering=ordering,
        majorana_alpha=_finite_float(
            params["majorana_alpha"],
            name="majorana_alpha",
        ),
        majorana_beta=_finite_float(params["majorana_beta"], name="majorana_beta"),
        max_abs_ybar=max_abs_ybar,
        source=source,
        matching_status=LMFV_LEPTON_PARAMETERS_MATCHING_STATUS_V1,
    )


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


def _is_lmfv_lepton_parameters(value: Any) -> bool:
    return isinstance(value, LMFVLeptonParameters)


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
    c_lfv: float | None = None,
    reference_scale_gev: float = 3000.0,
    m_kk_gev: float | None = None,
) -> MuToEGammaBranchingResult:
    """Evaluate ``BR(mu -> e gamma)`` from an LMFV carrier, YukawaResult, or proxy."""

    limit = _positive_finite(br_limit, name="br_limit")
    prefactor = _positive_finite(prefactor_br, name="prefactor_br")
    reference_scale = _positive_finite(
        reference_scale_gev,
        name="reference_scale_gev",
    )
    coefficient = (
        mu_to_e_gamma_coefficient_from_limit(limit, prefactor)
        if c_lfv is None
        else _positive_finite(c_lfv, name="c_lfv")
    )

    if isinstance(lepton_input, Mapping) and "yukawa_result" in lepton_input:
        lepton_input = lepton_input["yukawa_result"]

    if _is_lmfv_lepton_parameters(lepton_input):
        carrier = LMFVLeptonParameters(
            Y_N=lepton_input.Y_N,
            Y_N_bar=lepton_input.Y_N_bar,
            Y_N_matrix=lepton_input.Y_N_matrix,
            Y_N_bar_matrix=lepton_input.Y_N_bar_matrix,
            pmns=lepton_input.pmns,
            lmfv_spurion=lepton_input.lmfv_spurion,
            M_KK_gev=lepton_input.M_KK_gev,
            M_N_gev=lepton_input.M_N_gev,
            c_L=lepton_input.c_L,
            c_E=lepton_input.c_E,
            c_N=lepton_input.c_N,
            v_gev=lepton_input.v_gev,
            k_gev=lepton_input.k_gev,
            epsilon=lepton_input.epsilon,
            Lambda_IR_gev=lepton_input.Lambda_IR_gev,
            ordering=lepton_input.ordering,
            majorana_alpha=lepton_input.majorana_alpha,
            majorana_beta=lepton_input.majorana_beta,
            max_abs_ybar=lepton_input.max_abs_ybar,
            source=lepton_input.source,
            matching_status=lepton_input.matching_status,
        )
        m_kk = _positive_finite(carrier.M_KK_gev, name="LMFVLeptonParameters.M_KK_gev")
        if m_kk_gev is not None and not math.isclose(
            _positive_finite(m_kk_gev, name="m_kk_gev"),
            m_kk,
            rel_tol=1.0e-12,
            abs_tol=1.0e-9,
        ):
            raise ValueError("m_kk_gev override must match LMFVLeptonParameters.M_KK_gev")
        core = _check_mu_to_e_gamma_raw(
            carrier.Y_N_bar,
            carrier.pmns,
            M_KK=m_kk,
            C=coefficient,
            reference_scale=reference_scale,
        )
        return _result_from_core(
            core,
            br_limit=limit,
            prefactor_br=prefactor,
            c_lfv=coefficient,
            reference_scale_gev=reference_scale,
            m_kk_gev=m_kk,
            input_kind="LMFVLeptonParameters",
            used_proxy=False,
            extra_diagnostics={
                "extra_used": "lepton_lmfv_parameters",
                "source": carrier.source,
                "matching_status": carrier.matching_status,
                "m_kk_convention": PEREZ_RANDALL_LFV_M_KK_CONVENTION,
                "m_kk_source": "LMFVLeptonParameters.M_KK_gev",
                "max_abs_ybar": float(carrier.max_abs_ybar),
            },
        )

    if _is_yukawa_result_like(lepton_input):
        m_kk = _m_kk_from_yukawa_result(lepton_input, m_kk_gev)
        core = _check_mu_to_e_gamma(
            lepton_input,
            C=coefficient,
            reference_scale=reference_scale,
            M_KK_override=m_kk,
        )
        return _result_from_core(
            core,
            br_limit=limit,
            prefactor_br=prefactor,
            c_lfv=coefficient,
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
        C=coefficient,
        reference_scale=reference_scale,
    )
    return _result_from_core(
        core,
        br_limit=limit,
        prefactor_br=prefactor,
        c_lfv=coefficient,
        reference_scale_gev=reference_scale,
        m_kk_gev=m_kk,
        input_kind="MuToEGammaProxyInput",
        used_proxy=True,
        extra_diagnostics={
            "proxy_source": proxy.source,
            "needs_human_physics": LEPTON_DIPOLE_PROXY_ASSUMPTION_V1,
        },
    )
