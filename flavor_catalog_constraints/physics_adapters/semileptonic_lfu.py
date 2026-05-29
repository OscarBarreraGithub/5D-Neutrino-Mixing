"""Adapter over :mod:`quarkConstraints.semileptonic_lfu`.

This is the catalog boundary for charged-current semileptonic LFU-ratio
machinery.  Constraint modules import this adapter only; the underlying
SM-ratio normalization hook and documented ``b -> c tau nu`` proxy live in
``quarkConstraints`` for reuse by B025 and later B026.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.semileptonic_lfu import (
    SEMILEPTONIC_LFU_INPUT_BUNDLE_V1,
    SEMILEPTONIC_LFU_MODEL_V1,
    SEMILEPTONIC_LFU_OPERATOR_CONVENTION,
    SEMILEPTONIC_LFU_RS_MATCHING_ASSUMPTION_V1,
    SemileptonicLFUInputs,
    SemileptonicLFUResult,
    SemileptonicLFUWilsonProxy,
    compute_semileptonic_lfu_wilson_proxy as _compute_semileptonic_lfu_wilson_proxy,
    default_inputs as _default_inputs,
    evaluate_rd_lfu_ratio as _evaluate_rd_lfu_ratio,
    inputs_with_sm_ratio as _inputs_with_sm_ratio,
    ratio_from_scalar_shift as _ratio_from_scalar_shift,
    sm_lfu_ratio as _sm_lfu_ratio,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "SEMILEPTONIC_LFU_MODEL_V1",
    "SEMILEPTONIC_LFU_OPERATOR_CONVENTION",
    "SEMILEPTONIC_LFU_INPUT_BUNDLE_V1",
    "SEMILEPTONIC_LFU_RS_MATCHING_ASSUMPTION_V1",
    "SemileptonicLFUInputs",
    "SemileptonicLFUWilsonProxy",
    "SemileptonicLFUResult",
    "semileptonic_lfu_default_inputs",
    "semileptonic_lfu_inputs_with_sm_ratio",
    "semileptonic_lfu_wilson_proxy_from_couplings",
    "semileptonic_lfu_ratio_from_scalar_shift",
    "semileptonic_lfu_sm_ratio",
    "rd_lfu_ratio_from_couplings",
]


def semileptonic_lfu_default_inputs() -> SemileptonicLFUInputs:
    """Return the default semileptonic LFU proxy input bundle."""

    return _default_inputs()


def semileptonic_lfu_inputs_with_sm_ratio(
    sm_lfu_ratio: float,
    *,
    mode: str = "B->D",
    inputs: SemileptonicLFUInputs | None = None,
) -> SemileptonicLFUInputs:
    """Return inputs with a catalog-loaded SM LFU ratio installed."""

    return _inputs_with_sm_ratio(sm_lfu_ratio, mode=mode, inputs=inputs)


def semileptonic_lfu_wilson_proxy_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: SemileptonicLFUInputs | None = None,
) -> SemileptonicLFUWilsonProxy:
    """Return the v1 charged-current LFU Wilson proxy."""

    return _compute_semileptonic_lfu_wilson_proxy(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def semileptonic_lfu_ratio_from_scalar_shift(
    scalar_amplitude_shift: complex = 0.0j,
    *,
    inputs: SemileptonicLFUInputs | None = None,
) -> SemileptonicLFUResult:
    """Evaluate the LFU ratio from an explicit scalar amplitude shift."""

    return _ratio_from_scalar_shift(scalar_amplitude_shift, inputs=inputs)


def semileptonic_lfu_sm_ratio(
    inputs: SemileptonicLFUInputs | None = None,
) -> SemileptonicLFUResult:
    """Evaluate the SM-limit LFU ratio from the supplied input bundle."""

    return _sm_lfu_ratio(inputs)


def rd_lfu_ratio_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: SemileptonicLFUInputs | None = None,
) -> SemileptonicLFUResult:
    """Evaluate ``R(D)`` from mass-basis couplings with the v1 proxy."""

    return _evaluate_rd_lfu_ratio(couplings, m_kk_gev=m_kk_gev, inputs=inputs)
