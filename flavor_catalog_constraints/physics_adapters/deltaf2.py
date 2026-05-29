"""Adapter over :mod:`quarkConstraints.deltaf2` (the Delta F = 2 core).

This is the import boundary the kaon / beauty / charm mixing constraints
use to reach the Delta F = 2 physics. Per the append-only convention,
constraints may *add* wrappers here but must not change the signature of
an existing one; an upstream signature change is absorbed inside the
wrapper body so callers see a stable surface.

The upstream result dataclasses (``EpsilonKResult``, ``DeltaMKResult``)
are re-exported unchanged.
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
    evaluate_epsilon_k_with_running as _evaluate_epsilon_k_with_running,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "DeltaF2WilsonCoefficients",
    "EpsilonKResult",
    "DeltaMKResult",
    "epsilon_k_wilsons_from_couplings",
    "epsilon_k_from_wilsons",
    "epsilon_k_from_wilsons_with_running",
    "epsilon_k_from_couplings",
    "delta_mk_from_couplings",
]


def _kaon_wilsons(couplings: QuarkMassBasisCouplings) -> DeltaF2WilsonCoefficients:
    """Return the kaon Delta F=2 Wilson coefficients for ``couplings``.

    Centralizes the "pick the epsilon_k entry out of the default input
    bundle" step so both wrappers below share it. If upstream renames
    the bundle key, this is the one place to update.
    """
    wilsons = compute_delta_f2_wilsons(couplings, inputs=DEFAULT_DELTA_F2_INPUTS_V1)
    for w in wilsons:
        if w.input.key == "epsilon_k":
            return w
    raise RuntimeError("epsilon_k entry missing from DEFAULT_DELTA_F2_INPUTS_V1")


def epsilon_k_from_couplings(couplings: QuarkMassBasisCouplings) -> EpsilonKResult:
    """Compute the NP contribution to epsilon_K from mass-basis couplings."""
    return _evaluate_epsilon_k(_kaon_wilsons(couplings))


def delta_mk_from_couplings(couplings: QuarkMassBasisCouplings) -> DeltaMKResult:
    """Compute the NP contribution to Delta m_K from mass-basis couplings."""
    return _evaluate_delta_mk(_kaon_wilsons(couplings))


def epsilon_k_wilsons_from_couplings(
    couplings: QuarkMassBasisCouplings,
) -> DeltaF2WilsonCoefficients:
    """Return the kaon Delta F=2 Wilson coefficients used for epsilon_K."""
    return _kaon_wilsons(couplings)


def epsilon_k_from_wilsons(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    epsilon_k_np_budget: float | None = None,
) -> EpsilonKResult:
    """Compute epsilon_K from kaon Wilsons, optionally with a catalog budget."""
    return _evaluate_epsilon_k(
        wilsons,
        epsilon_k_np_budget_override=epsilon_k_np_budget,
    )


def epsilon_k_from_wilsons_with_running(
    wilsons: DeltaF2WilsonCoefficients,
    *,
    mu_had: float = 2.0,
    epsilon_k_np_budget: float | None = None,
) -> EpsilonKResult:
    """Compute epsilon_K after QCD-evolving kaon Wilsons to ``mu_had``."""
    return _evaluate_epsilon_k_with_running(
        wilsons,
        mu_had=mu_had,
        epsilon_k_np_budget_override=epsilon_k_np_budget,
    )
