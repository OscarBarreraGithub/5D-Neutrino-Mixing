"""Adapter over :mod:`quarkConstraints.zpole`.

Constraint modules import this adapter only.  The underlying Z-pole
pseudo-observable arithmetic and the documented RS ``Zbb`` coupling-shift
proxy live in ``quarkConstraints.zpole``.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.zpole import (
    ZPOLE_DOWN_FCNC_PROXY_V1,
    ZPOLE_INPUT_BUNDLE_V1,
    ZPOLE_MODEL_V1,
    ZPOLE_RS_ZBB_PROXY_V1,
    ZPoleCouplings,
    ZPoleDownFCNCBranchingResult,
    ZPoleDownFCNCCouplingProxy,
    ZbbCouplingShiftProxy,
    ZPoleQuarkObservables,
    ZPoleSMInputs,
    asymmetry_parameter as _asymmetry_parameter,
    default_sm_inputs as _default_sm_inputs,
    down_fcnc_branching_fraction_from_couplings as _down_fcnc_branching_from_couplings,
    down_fcnc_branching_fraction_with_proxy as _down_fcnc_branching_with_proxy,
    down_fcnc_coupling_proxy as _down_fcnc_coupling_proxy,
    down_fcnc_effective_coupling_limit as _down_fcnc_effective_coupling_limit,
    down_fcnc_sm_hadronic_width_weight as _down_fcnc_sm_hadronic_width_weight,
    down_fcnc_sm_total_width_weight as _down_fcnc_sm_total_width_weight,
    evaluate_quark_pseudo_observables as _evaluate_quark_pseudo_observables,
    evaluate_zbb_with_proxy as _evaluate_zbb_with_proxy,
    forward_backward_asymmetry as _forward_backward_asymmetry,
    inputs_with_bottom_radiator as _inputs_with_bottom_radiator,
    partial_width_weight as _partial_width_weight,
    r_quark as _r_quark,
    shifted_couplings as _shifted_couplings,
    sm_couplings as _sm_couplings,
    zbb_coupling_shift_proxy as _zbb_coupling_shift_proxy,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "ZPOLE_MODEL_V1",
    "ZPOLE_INPUT_BUNDLE_V1",
    "ZPOLE_RS_ZBB_PROXY_V1",
    "ZPoleSMInputs",
    "ZPoleCouplings",
    "ZPoleQuarkObservables",
    "ZbbCouplingShiftProxy",
    "ZPoleDownFCNCBranchingResult",
    "ZPoleDownFCNCCouplingProxy",
    "ZPOLE_DOWN_FCNC_PROXY_V1",
    "zpole_default_sm_inputs",
    "zpole_inputs_with_bottom_radiator",
    "zpole_sm_couplings",
    "zpole_shifted_couplings",
    "zpole_asymmetry_parameter",
    "zpole_forward_backward_asymmetry",
    "zpole_partial_width_weight",
    "zpole_r_quark",
    "zpole_evaluate_quark",
    "zpole_zbb_coupling_shift_proxy",
    "zpole_evaluate_zbb_with_proxy",
    "zpole_down_fcnc_sm_hadronic_width_weight",
    "zpole_down_fcnc_sm_total_width_weight",
    "zpole_down_fcnc_branching_fraction_from_couplings",
    "zpole_down_fcnc_effective_coupling_limit",
    "zpole_down_fcnc_coupling_proxy",
    "zpole_down_fcnc_branching_fraction_with_proxy",
]


def zpole_default_sm_inputs() -> ZPoleSMInputs:
    """Return the default Z-pole pseudo-observable input bundle."""
    return _default_sm_inputs()


def zpole_inputs_with_bottom_radiator(
    target_r_b: float,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleSMInputs:
    """Return inputs whose SM-limit ``R_b`` equals ``target_r_b``."""
    return _inputs_with_bottom_radiator(target_r_b, inputs)


def zpole_sm_couplings(
    flavor: str,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleCouplings:
    """Return SM effective chiral Z couplings for ``flavor``."""
    return _sm_couplings(flavor, inputs)


def zpole_shifted_couplings(
    base: ZPoleCouplings,
    *,
    delta_g_left: complex = 0.0j,
    delta_g_right: complex = 0.0j,
) -> ZPoleCouplings:
    """Return shifted effective chiral Z couplings."""
    return _shifted_couplings(
        base,
        delta_g_left=delta_g_left,
        delta_g_right=delta_g_right,
    )


def zpole_asymmetry_parameter(couplings: ZPoleCouplings) -> float:
    """Return ``A_f`` from effective chiral couplings."""
    return _asymmetry_parameter(couplings)


def zpole_forward_backward_asymmetry(
    final_state: ZPoleCouplings,
    *,
    initial_state: ZPoleCouplings | None = None,
) -> float:
    """Return ``A_FB^0,f`` from initial and final asymmetry parameters."""
    return _forward_backward_asymmetry(final_state, initial_state=initial_state)


def zpole_partial_width_weight(
    couplings: ZPoleCouplings,
    *,
    radiator: float = 1.0,
    n_color: int | None = None,
) -> float:
    """Return a relative partial-width weight."""
    return _partial_width_weight(couplings, radiator=radiator, n_color=n_color)


def zpole_r_quark(
    target_quark: str,
    couplings_by_flavor: dict[str, ZPoleCouplings] | None = None,
    *,
    inputs: ZPoleSMInputs | None = None,
) -> float:
    """Return ``R_q`` for a target quark."""
    return _r_quark(target_quark, couplings_by_flavor, inputs=inputs)


def zpole_evaluate_quark(
    target_quark: str,
    couplings_by_flavor: dict[str, ZPoleCouplings] | None = None,
    *,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleQuarkObservables:
    """Evaluate ``R_q``, ``A_q``, and ``A_FB^0,q``."""
    return _evaluate_quark_pseudo_observables(
        target_quark,
        couplings_by_flavor,
        inputs=inputs,
    )


def zpole_zbb_coupling_shift_proxy(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> ZbbCouplingShiftProxy:
    """Return the documented RS ``Zbb`` coupling-shift proxy."""
    return _zbb_coupling_shift_proxy(couplings, m_kk_gev=m_kk_gev, inputs=inputs)


def zpole_evaluate_zbb_with_proxy(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> tuple[ZPoleQuarkObservables, ZbbCouplingShiftProxy]:
    """Evaluate ``Z -> b bbar`` pseudo-observables with the RS proxy."""
    return _evaluate_zbb_with_proxy(couplings, m_kk_gev=m_kk_gev, inputs=inputs)


def zpole_down_fcnc_sm_hadronic_width_weight(
    inputs: ZPoleSMInputs | None = None,
) -> dict[str, float]:
    """Return SM hadronic Z-width weights from the shared Z-pole convention."""

    return _down_fcnc_sm_hadronic_width_weight(inputs)


def zpole_down_fcnc_sm_total_width_weight(
    inputs: ZPoleSMInputs | None = None,
) -> dict[str, float]:
    """Return SM total Z-width weights from the shared Z-pole convention."""

    return _down_fcnc_sm_total_width_weight(inputs)


def zpole_down_fcnc_branching_fraction_from_couplings(
    *,
    flavor_i: str,
    flavor_j: str,
    delta_g_left: complex,
    delta_g_right: complex,
    br_limit: float | None = None,
    inputs: ZPoleSMInputs | None = None,
    charge_state_factor: float = 2.0,
) -> ZPoleDownFCNCBranchingResult:
    """Return ``BR(Z -> q_i qbar_j + q_j qbar_i)`` from FCNC couplings."""

    return _down_fcnc_branching_from_couplings(
        flavor_i=flavor_i,
        flavor_j=flavor_j,
        delta_g_left=delta_g_left,
        delta_g_right=delta_g_right,
        br_limit=br_limit,
        inputs=inputs,
        charge_state_factor=charge_state_factor,
    )


def zpole_down_fcnc_effective_coupling_limit(
    br_limit: float,
    *,
    flavor_i: str,
    flavor_j: str,
    inputs: ZPoleSMInputs | None = None,
    charge_state_factor: float = 2.0,
) -> float:
    """Return the limit on ``sqrt(|delta_g_L|^2 + |delta_g_R|^2)``."""

    return _down_fcnc_effective_coupling_limit(
        br_limit,
        flavor_i=flavor_i,
        flavor_j=flavor_j,
        inputs=inputs,
        charge_state_factor=charge_state_factor,
    )


def zpole_down_fcnc_coupling_proxy(
    source: QuarkMassBasisCouplings,
    *,
    flavor_i: str,
    flavor_j: str,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
) -> ZPoleDownFCNCCouplingProxy:
    """Map quark mass-basis overlaps onto an off-diagonal down-sector Z proxy."""

    return _down_fcnc_coupling_proxy(
        source,
        flavor_i=flavor_i,
        flavor_j=flavor_j,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def zpole_down_fcnc_branching_fraction_with_proxy(
    source: QuarkMassBasisCouplings,
    *,
    flavor_i: str,
    flavor_j: str,
    br_limit: float | None = None,
    m_kk_gev: float | None = None,
    inputs: ZPoleSMInputs | None = None,
    charge_state_factor: float = 2.0,
) -> tuple[ZPoleDownFCNCBranchingResult, ZPoleDownFCNCCouplingProxy]:
    """Evaluate an off-diagonal down-sector Z decay with the documented proxy."""

    return _down_fcnc_branching_with_proxy(
        source,
        flavor_i=flavor_i,
        flavor_j=flavor_j,
        br_limit=br_limit,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
        charge_state_factor=charge_state_factor,
    )
