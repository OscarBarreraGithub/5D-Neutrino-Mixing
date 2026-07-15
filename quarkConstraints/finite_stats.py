"""Finite-sample statistical helpers for scan acceptance rates."""

from __future__ import annotations

import math


def wilson_upper_limit(k: int, n: int, z: float = 1.92) -> float:
    """Return the Wilson-score upper limit for ``k`` successes in ``n`` trials.

    The default ``z=1.92`` follows the convention used in the quark-scan
    zero-pass audit, where a zero-success observation gives
    ``z**2 / (n + z**2) ~= 3.69 / n`` for large ``n``.
    This is a local audit convention, not a standard 95 percent one-sided
    upper limit; as a Gaussian score it corresponds to about 97.3 percent
    one-sided coverage.
    """
    if n <= 0:
        raise ValueError("n must be positive")
    if k < 0:
        raise ValueError("k must be non-negative")
    if k > n:
        raise ValueError("k cannot exceed n")
    if z <= 0 or not math.isfinite(z):
        raise ValueError("z must be positive and finite")

    z2 = z * z
    root = math.sqrt(z2 + 4.0 * k * (1.0 - k / n))
    return (z2 + 2.0 * k + z * root) / (2.0 * (n + z2))
