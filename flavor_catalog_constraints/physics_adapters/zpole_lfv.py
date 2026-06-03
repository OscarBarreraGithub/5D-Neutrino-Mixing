"""Adapter over the LFV extension of the shared Z-pole machinery.

Constraint modules import this adapter only.  The effective-coupling width
weights are reused from :mod:`quarkConstraints.zpole`; the new LFV glue lives
in :mod:`quarkConstraints.zpole_lfv`.

The Phase-4b constraint path uses the direct
``z_lfv_branching_fraction_from_couplings`` export with explicit
``rs_ew_couplings.z_delta_g_L/R_e`` matrix elements.  Legacy proxy helpers are
kept for older callers, but T015-T017 do not use them.
"""

from __future__ import annotations

from quarkConstraints.zpole_lfv import (
    ZPOLE_LFV_MODEL_V1,
    ZPOLE_LFV_PROXY_V1,
    ZPoleLFVBranchingResult,
    ZPoleLFVCouplingProxy,
    ZPoleLFVProxyInput,
    sm_total_width_weight as _sm_total_width_weight,
    z_lfv_branching_fraction_from_couplings as _branching_from_couplings,
    z_lfv_branching_fraction_with_proxy as _branching_with_proxy,
    z_lfv_coupling_proxy as _coupling_proxy,
    z_lfv_effective_coupling_limit as _effective_coupling_limit,
    z_lfv_proxy_input as _proxy_input,
)
from quarkConstraints.zpole import ZPoleSMInputs

__all__ = [
    "ZPOLE_LFV_MODEL_V1",
    "ZPOLE_LFV_PROXY_V1",
    "ZPoleSMInputs",
    "ZPoleLFVBranchingResult",
    "ZPoleLFVCouplingProxy",
    "ZPoleLFVProxyInput",
    "zpole_lfv_proxy_input",
    "zpole_lfv_sm_total_width_weight",
    "zpole_lfv_branching_fraction_from_couplings",
    "z_lfv_branching_fraction_from_couplings",
    "zpole_lfv_effective_coupling_limit",
    "zpole_lfv_coupling_proxy",
    "zpole_lfv_branching_fraction_with_proxy",
]


def zpole_lfv_proxy_input(
    left_emu_overlap: complex,
    right_emu_overlap: complex,
    m_kk_gev: float,
    *,
    source: str = "caller-supplied lepton neutral-current proxy",
) -> ZPoleLFVProxyInput:
    """Build a proxy input accepted by the LFV Z-pole adapter."""

    return _proxy_input(
        left_emu_overlap,
        right_emu_overlap,
        m_kk_gev,
        source=source,
    )


def zpole_lfv_sm_total_width_weight(
    inputs: ZPoleSMInputs | None = None,
) -> dict[str, float]:
    """Return SM total Z-width weights from the shared Z-pole convention."""

    return _sm_total_width_weight(inputs)


def zpole_lfv_branching_fraction_from_couplings(
    *,
    delta_g_left: complex,
    delta_g_right: complex,
    initial_flavor: str = "e",
    final_flavor: str = "mu",
    br_limit: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleLFVBranchingResult:
    """Return ``BR(Z -> l_i l_j)`` from off-diagonal effective couplings."""

    return _branching_from_couplings(
        delta_g_left=delta_g_left,
        delta_g_right=delta_g_right,
        initial_flavor=initial_flavor,
        final_flavor=final_flavor,
        br_limit=br_limit,
        inputs=inputs,
    )


def z_lfv_branching_fraction_from_couplings(
    *,
    delta_g_left: complex,
    delta_g_right: complex,
    initial_flavor: str = "e",
    final_flavor: str = "mu",
    br_limit: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleLFVBranchingResult:
    """Adapter-exported LFV Z branching fraction helper used by constraints."""

    return zpole_lfv_branching_fraction_from_couplings(
        delta_g_left=delta_g_left,
        delta_g_right=delta_g_right,
        initial_flavor=initial_flavor,
        final_flavor=final_flavor,
        br_limit=br_limit,
        inputs=inputs,
    )


def zpole_lfv_effective_coupling_limit(
    br_limit: float,
    *,
    inputs: ZPoleSMInputs | None = None,
) -> float:
    """Return the limit on ``sqrt(|delta_g_L|^2 + |delta_g_R|^2)``."""

    return _effective_coupling_limit(br_limit, inputs=inputs)


def zpole_lfv_coupling_proxy(
    source: object,
    *,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleLFVCouplingProxy:
    """Return the documented lepton-overlap proxy for ``Z e mu`` couplings."""

    return _coupling_proxy(source, m_kk_gev=m_kk_gev, inputs=inputs)


def zpole_lfv_branching_fraction_with_proxy(
    source: object,
    *,
    br_limit: float | None = None,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> tuple[ZPoleLFVBranchingResult, ZPoleLFVCouplingProxy]:
    """Evaluate ``BR(Z -> e mu)`` with the documented proxy."""

    return _branching_with_proxy(
        source,
        br_limit=br_limit,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
