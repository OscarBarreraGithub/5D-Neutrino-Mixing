"""Adapter for :mod:`quarkConstraints.deltaf2`.

This is the import boundary that catalogued constraints (K001, …) use
to reach the Delta F = 2 physics core. Per the append-only adapter
convention (plan §E, R1 M-6), constraints may **add** new wrapper
functions to this module but may not modify the signature of an
existing one.

The adapter passes the physics-module dataclasses
(``EpsilonKResult``, ``DeltaMKResult``, …) through unchanged so that
a signature change in ``quarkConstraints.deltaf2`` touches *one* file
(this one), not every constraint.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.deltaf2 import (
    DEFAULT_DELTA_F2_INPUTS_V1,
    DeltaF2WilsonCoefficients,
    DeltaMKResult,
    EpsilonKResult,
    compute_delta_f2_wilsons,
    evaluate_delta_mk as _evaluate_delta_mk,
    evaluate_epsilon_k as _evaluate_epsilon_k,
)

__all__ = [
    "DEFAULT_DELTA_F2_INPUTS_V1",
    "DeltaF2WilsonCoefficients",
    "DeltaMKResult",
    "EpsilonKResult",
    "QuarkMassBasisCouplings",
    "evaluate_delta_mk_from_couplings",
    "evaluate_delta_mk_from_wilsons",
    "evaluate_epsilon_k_from_couplings",
    "evaluate_epsilon_k_from_wilsons",
]


def evaluate_epsilon_k_from_wilsons(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    epsilon_k_np_budget_override: float | None = None,
) -> EpsilonKResult:
    """Thin pass-through to :func:`quarkConstraints.deltaf2.evaluate_epsilon_k`.

    Returns the underlying :class:`EpsilonKResult` dataclass unchanged.
    """
    return _evaluate_epsilon_k(
        wilsons,
        epsilon_k_np_budget_override=epsilon_k_np_budget_override,
    )


def evaluate_epsilon_k_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    epsilon_k_np_budget_override: float | None = None,
) -> EpsilonKResult:
    """Compute Delta F=2 Wilsons then evaluate epsilon_K.

    Convenience wrapper for constraints that hold a
    :class:`QuarkMassBasisCouplings` extra and want to skip the
    explicit Wilson-coefficient step. Picks out the ``epsilon_k``
    entry from the default Delta F=2 input bundle.
    """
    wilsons_tuple = compute_delta_f2_wilsons(
        couplings,
        inputs=DEFAULT_DELTA_F2_INPUTS_V1,
    )
    eps_wilsons = None
    for w in wilsons_tuple:
        if w.input.key == "epsilon_k":
            eps_wilsons = w
            break
    if eps_wilsons is None:
        raise RuntimeError(
            "epsilon_k entry missing from DEFAULT_DELTA_F2_INPUTS_V1"
        )
    return _evaluate_epsilon_k(
        eps_wilsons,
        epsilon_k_np_budget_override=epsilon_k_np_budget_override,
    )


def evaluate_delta_mk_from_wilsons(
    wilsons: DeltaF2WilsonCoefficients,
) -> DeltaMKResult:
    """Thin pass-through to :func:`quarkConstraints.deltaf2.evaluate_delta_mk`.

    Returns the underlying :class:`DeltaMKResult` dataclass unchanged.
    """
    return _evaluate_delta_mk(wilsons)


def evaluate_delta_mk_from_couplings(
    couplings: QuarkMassBasisCouplings,
) -> DeltaMKResult:
    """Compute Delta F=2 Wilsons then evaluate Delta m_K.

    Convenience wrapper for constraints that hold a
    :class:`QuarkMassBasisCouplings` extra and want to skip the
    explicit Wilson-coefficient step. Picks out the kaon ``epsilon_k``
    entry from the default Delta F=2 input bundle; the same kaon Wilsons
    feed both ``epsilon_K`` and ``Delta m_K``.
    """
    wilsons_tuple = compute_delta_f2_wilsons(
        couplings,
        inputs=DEFAULT_DELTA_F2_INPUTS_V1,
    )
    kaon_wilsons = None
    for w in wilsons_tuple:
        if w.input.key == "epsilon_k":
            kaon_wilsons = w
            break
    if kaon_wilsons is None:
        raise RuntimeError(
            "epsilon_k entry missing from DEFAULT_DELTA_F2_INPUTS_V1"
        )
    return _evaluate_delta_mk(kaon_wilsons)
