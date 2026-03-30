"""Scale conventions for the quark-sector MFV module."""

from __future__ import annotations

DEFAULT_QUARK_XI_KK = 1.0
DEFAULT_QUARK_TARGET_SCALE_GEV = 3000.0


def default_quark_m_kk_from_lambda_ir(
    Lambda_IR: float,
    xi_KK: float = DEFAULT_QUARK_XI_KK,
) -> float:
    """Map the geometric IR scale to the quark-sector KK convention.

    Parameters
    ----------
    Lambda_IR : float
        Geometric IR scale, ``Lambda_IR = 1 / z_v``.
    xi_KK : float, optional
        Explicit KK-mass convention. The repo default keeps ``xi_KK = 1.0``,
        i.e. the bookkeeping choice ``M_KK = Lambda_IR``.
    """
    if Lambda_IR <= 0:
        raise ValueError("Lambda_IR must be positive")
    if xi_KK <= 0:
        raise ValueError("xi_KK must be positive")
    return float(xi_KK * Lambda_IR)
