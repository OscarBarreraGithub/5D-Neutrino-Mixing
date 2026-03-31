"""Mass-basis KK-gluon couplings for the quark-sector MFV module.

This layer rotates the non-universal zero-mode overlap profiles into the exact
quark mass basis returned by :mod:`quarkConstraints.fit`. The resulting
Hermitian matrices are the building blocks for downstream ``Delta F = 2``
matching.

The repo default remains the bookkeeping convention ``M_KK = Lambda_IR``.
Callers that want a different physical KK-mass convention should pass an
explicit ``xi_KK`` or ``M_KK``.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt

import numpy as np

from qcd import alpha_s

from .fit import QuarkFitResult


def _default_quark_m_kk_from_lambda_ir(Lambda_IR: float, xi_KK: float = 1.0) -> float:
    """Resolve the KK scale using the public helper when it is available."""
    try:
        from .scales import default_quark_m_kk_from_lambda_ir

        return default_quark_m_kk_from_lambda_ir(Lambda_IR, xi_KK=xi_KK)
    except Exception:
        if Lambda_IR <= 0.0:
            raise ValueError("Lambda_IR must be positive")
        if xi_KK <= 0.0:
            raise ValueError("xi_KK must be positive")
        return float(xi_KK * Lambda_IR)


def _hermitian(matrix: np.ndarray) -> np.ndarray:
    """Project small numerical noise back onto the Hermitian subspace."""
    arr = np.asarray(matrix, dtype=np.complex128)
    return 0.5 * (arr + arr.conjugate().T)


def _mass_basis_overlap(rotation: np.ndarray, profile_values: np.ndarray) -> np.ndarray:
    """Rotate a diagonal overlap profile into the corresponding mass basis."""
    profile = np.diag(np.asarray(profile_values, dtype=float) ** 2)
    return _hermitian(rotation.conjugate().T @ profile @ rotation)


def _off_diagonal_norm(matrix: np.ndarray) -> float:
    """Return the Frobenius norm of off-diagonal entries only."""
    arr = np.asarray(matrix, dtype=np.complex128).copy()
    np.fill_diagonal(arr, 0.0)
    return float(np.linalg.norm(arr, ord="fro"))


@dataclass(frozen=True)
class QuarkMassBasisCouplings:
    """KK-gluon couplings in the exact quark mass basis.

    The ``*_overlap`` matrices are the dimensionless overlap structures. The
    ``left_*`` / ``right_*`` matrices multiply them by the running QCD gauge
    coupling ``g_s(M_KK)``.
    """

    M_KK: float
    xi_KK: float
    alpha_s: float
    g_s: float
    left_overlap: np.ndarray
    right_up_overlap: np.ndarray
    right_down_overlap: np.ndarray
    left_up: np.ndarray
    left_down: np.ndarray
    right_up: np.ndarray
    right_down: np.ndarray

    @property
    def left_down_offdiag_norm(self) -> float:
        return _off_diagonal_norm(self.left_down)

    @property
    def left_up_offdiag_norm(self) -> float:
        return _off_diagonal_norm(self.left_up)

    @property
    def right_down_offdiag_norm(self) -> float:
        return _off_diagonal_norm(self.right_down)

    @property
    def right_up_offdiag_norm(self) -> float:
        return _off_diagonal_norm(self.right_up)


def compute_quark_kk_gluon_couplings(
    fit_result: QuarkFitResult,
    *,
    M_KK: float | None = None,
    xi_KK: float = 1.0,
) -> QuarkMassBasisCouplings:
    """Return mass-basis KK-gluon couplings for a fitted MFV point.

    Parameters
    ----------
    fit_result
        Exact quark fit result from :func:`quarkConstraints.fit.fit_quark_sector`.
    M_KK
        Explicit KK scale in GeV. When omitted, the helper uses the repo's
        default quark convention ``M_KK = xi_KK * Lambda_IR``.
    xi_KK
        Conversion factor between the geometric IR scale and the KK scale when
        ``M_KK`` is not supplied. The repo default is ``1.0``.
    """
    if xi_KK <= 0.0:
        raise ValueError("xi_KK must be positive")

    state = fit_result.state
    resolved_mkk = (
        float(M_KK)
        if M_KK is not None
        else _default_quark_m_kk_from_lambda_ir(state.point.Lambda_IR, xi_KK=xi_KK)
    )
    if resolved_mkk <= 0.0:
        raise ValueError("M_KK must be positive")
    resolved_xi_kk = (
        float(resolved_mkk / state.point.Lambda_IR)
        if M_KK is not None
        else float(xi_KK)
    )

    running_alpha_s = float(alpha_s(resolved_mkk, precision="high"))
    g_s = float(sqrt(4.0 * pi * running_alpha_s))

    left_overlap = _mass_basis_overlap(fit_result.U_L_d, state.F_Q)
    right_down_overlap = _mass_basis_overlap(fit_result.U_R_d, state.F_d)
    right_up_overlap = _mass_basis_overlap(fit_result.U_R_u, state.F_u)
    left_up_overlap = _mass_basis_overlap(fit_result.U_L_u, state.F_Q)

    return QuarkMassBasisCouplings(
        M_KK=resolved_mkk,
        xi_KK=resolved_xi_kk,
        alpha_s=running_alpha_s,
        g_s=g_s,
        left_overlap=left_overlap,
        right_up_overlap=right_up_overlap,
        right_down_overlap=right_down_overlap,
        left_up=_hermitian(g_s * left_up_overlap),
        left_down=_hermitian(g_s * left_overlap),
        right_up=_hermitian(g_s * right_up_overlap),
        right_down=_hermitian(g_s * right_down_overlap),
    )
