"""Shared CGHNP profile brackets translated to repo conventions."""

from __future__ import annotations

import math


def casagrande_cghnp_B_profile(c: float, F: float, *, name: str) -> float:
    """Diagonal CGHNP fermion-KK bracket in repo variables.

    CGHNP (0807.4937) uses the convention dictionary
    ``c_CGHNP = -c_repo`` and ``F_CGHNP^2 = 2 f_IR,repo^2``.  In repo variables:

        B(c, F) = 1/(1 + 2c) * (1/(2 F^2) - 1 + 2 F^2/(3 - 2c)).

    The raw CGHNP-looking ``1/(1 - 2c)`` form is wrong when fed repo ``c`` and
    repo ``f_IR``; this helper is shared by the B1 Zbb fix and Higgs-Yukawa
    matching so the convention transformation stays single-sourced.
    """

    if not math.isfinite(c):
        raise ValueError(f"{name} c must be finite")
    if not math.isfinite(F) or F <= 0.0:
        raise ValueError(f"{name} F must be positive and finite")
    denom_left = 1.0 + 2.0 * c
    denom_right = 3.0 - 2.0 * c
    if denom_left == 0.0 or denom_right == 0.0:
        raise ValueError(f"{name} Casagrande B(c) denominator is singular")
    f_sq = F * F
    value = (1.0 / denom_left) * (
        1.0 / (2.0 * f_sq) - 1.0 + (2.0 * f_sq) / denom_right
    )
    if not math.isfinite(value):
        raise ValueError(f"{name} Casagrande B(c) is non-finite")
    return float(value)


__all__ = ["casagrande_cghnp_B_profile"]
