"""Adapter over :mod:`quarkConstraints.zpole`.

Constraint modules import this adapter only.  The underlying Z-pole
pseudo-observable arithmetic and the documented RS ``Zbb`` coupling-shift
proxy live in ``quarkConstraints.zpole``.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.zpole import (
    ZPOLE_INPUT_BUNDLE_V1,
    ZPOLE_MODEL_V1,
    ZPOLE_RS_ZBB_PROXY_V1,
    ZPoleCouplings,
    ZbbCouplingShiftProxy,
    ZPoleQuarkObservables,
    ZPoleSMInputs,
    asymmetry_parameter as _asymmetry_parameter,
    default_sm_inputs as _default_sm_inputs,
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
