"""Session-level fixtures for ``tests/constraints/``.

Two responsibilities:

1. Trigger lazy ``ConstraintRegistry.discover()`` exactly once per
   pytest session (R1 REC-1 in the plan).
2. Provide ``sm_point`` and ``excluded_point`` fixtures that produce
   :class:`ParameterPoint` instances pre-populated with the
   ``quark_mass_basis_couplings`` extra. Constraints can then be
   exercised against a deterministic SM-like / NP-excluded benchmark
   without each test reinventing the helper geometry.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

# Make sure the repo root is on sys.path so ``import
# flavor_catalog_constraints`` works inside the test session, even if
# the package was not pip-installed in editable mode.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from flavor_catalog_constraints import ConstraintRegistry  # noqa: E402
from flavor_catalog_constraints.base import (  # noqa: E402
    ParameterPoint,
    ParameterPointExtras,
)
from quarkConstraints.couplings import QuarkMassBasisCouplings  # noqa: E402


@pytest.fixture(scope="session", autouse=True)
def _discover_constraints():
    """Trigger lazy discovery once per pytest session (R1 REC-1)."""
    ConstraintRegistry.discover()


def _zero_couplings(M_KK: float = 3000.0) -> QuarkMassBasisCouplings:
    """Return a ``QuarkMassBasisCouplings`` with all matrices zeroed.

    epsilon_K^NP is identically zero when all KK-gluon off-diagonals
    vanish (no Im(M_12^NP) source). This is the SM-baseline point —
    every Delta F=2 constraint must report ``passes=True`` on it.
    """
    zeros = np.zeros((3, 3), dtype=np.complex128)
    return QuarkMassBasisCouplings(
        M_KK=float(M_KK),
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=zeros.copy(),
        right_up_overlap=zeros.copy(),
        right_down_overlap=zeros.copy(),
        left_up=zeros.copy(),
        left_down=zeros.copy(),
        right_up=zeros.copy(),
        right_down=zeros.copy(),
    )


def _np_excluded_couplings(
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Return couplings with a complex (s,d) off-diagonal in ``left_down`` and
    a slightly different complex entry in ``right_down``.

    ``c4_lr = -(left * right) / M_KK^2`` carries an imaginary part as
    long as ``Im(left * right) != 0``; that part feeds
    ``Im(M_12^NP)`` and hence ``epsilon_K^NP``. At ``M_KK = 3 TeV``
    couplings ~ 5e-2 yield ratios many orders of magnitude above the
    experimental budget (~6.7e-5), which is the desired "obviously
    excluded" regime.
    """
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    sd_left = 0.05 + 0.02j
    sd_right = 0.05 - 0.03j
    # Hermitian off-diagonals.
    left_down[0, 1] = sd_left
    left_down[1, 0] = np.conj(sd_left)
    right_down[0, 1] = sd_right
    right_down[1, 0] = np.conj(sd_right)
    return QuarkMassBasisCouplings(
        M_KK=float(M_KK),
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=zeros.copy(),
        right_up_overlap=zeros.copy(),
        right_down_overlap=zeros.copy(),
        left_up=zeros.copy(),
        left_down=left_down,
        right_up=zeros.copy(),
        right_down=right_down,
    )


@pytest.fixture
def sm_point() -> ParameterPoint:
    """SM-baseline ``ParameterPoint``: all NP couplings vanish.

    Every Delta F = 2 constraint must report ``passes=True`` here.
    """
    extras: ParameterPointExtras = {
        "quark_mass_basis_couplings": _zero_couplings(),
        "kk_gluon_mass_gev": 3000.0,
    }
    return ParameterPoint(raw=None, extras=extras)


@pytest.fixture
def excluded_point() -> ParameterPoint:
    """Clearly excluded ``ParameterPoint``: large imaginary (s,d) coupling.

    Every kaon-mixing constraint must report ``passes=False`` here.
    """
    extras: ParameterPointExtras = {
        "quark_mass_basis_couplings": _np_excluded_couplings(),
        "kk_gluon_mass_gev": 3000.0,
    }
    return ParameterPoint(raw=None, extras=extras)
