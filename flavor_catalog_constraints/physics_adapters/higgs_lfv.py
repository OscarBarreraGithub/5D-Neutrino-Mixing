"""Adapter over :mod:`quarkConstraints.higgs_lfv`.

Constraint modules import this adapter only.  The reusable effective-Yukawa
width formula and the documented Higgs LFV proxy live in
``quarkConstraints.higgs_lfv``.

NEEDS-HUMAN-PHYSICS: ``ParameterPoint`` does not carry a rigorous
off-diagonal charged-lepton Higgs-Yukawa matrix.  The proxy path accepts
caller-supplied effective Yukawa entries and is diagnostic-only until a full
RS Higgs/lepton matching object exists.
"""

from __future__ import annotations

from quarkConstraints.higgs_lfv import (
    HIGGS_LFV_INPUT_BUNDLE_V1,
    HIGGS_LFV_MODEL_V1,
    HIGGS_LFV_RS_PROXY_V1,
    HIGGS_LFV_YUKAWA_CONVENTION,
    HiggsLFVBranchingResult,
    HiggsLFVInputs,
    HiggsLFVYukawaProxy,
    HiggsLFVYukawaProxyInput,
)
from quarkConstraints.higgs_lfv import (
    default_higgs_lfv_inputs as _default_inputs,
)
from quarkConstraints.higgs_lfv import (
    h_lfv_branching_fraction_from_yukawas as _branching_from_yukawas,
)
from quarkConstraints.higgs_lfv import (
    h_lfv_branching_fraction_with_proxy as _branching_with_proxy,
)
from quarkConstraints.higgs_lfv import (
    h_lfv_effective_yukawa_limit as _effective_yukawa_limit,
)
from quarkConstraints.higgs_lfv import (
    h_lfv_partial_width as _partial_width,
)
from quarkConstraints.higgs_lfv import (
    h_lfv_yukawa_proxy as _yukawa_proxy,
)
from quarkConstraints.higgs_lfv import (
    h_lfv_yukawa_proxy_input as _proxy_input,
)

__all__ = [
    "HIGGS_LFV_MODEL_V1",
    "HIGGS_LFV_INPUT_BUNDLE_V1",
    "HIGGS_LFV_YUKAWA_CONVENTION",
    "HIGGS_LFV_RS_PROXY_V1",
    "HiggsLFVInputs",
    "HiggsLFVYukawaProxyInput",
    "HiggsLFVYukawaProxy",
    "HiggsLFVBranchingResult",
    "higgs_lfv_default_inputs",
    "higgs_lfv_yukawa_proxy_input",
    "higgs_lfv_partial_width",
    "higgs_lfv_branching_fraction_from_yukawas",
    "higgs_lfv_effective_yukawa_limit",
    "higgs_lfv_yukawa_proxy",
    "higgs_lfv_branching_fraction_with_proxy",
]


def higgs_lfv_default_inputs() -> HiggsLFVInputs:
    """Return the default Higgs LFV input bundle."""

    return _default_inputs()


def higgs_lfv_yukawa_proxy_input(
    initial_flavor: str,
    final_flavor: str,
    yukawa_ij: complex,
    yukawa_ji: complex = 0.0j,
    *,
    source: str = "caller-supplied Higgs LFV Yukawa proxy",
) -> HiggsLFVYukawaProxyInput:
    """Build a proxy input accepted by the Higgs LFV adapter."""

    return _proxy_input(
        initial_flavor,
        final_flavor,
        yukawa_ij,
        yukawa_ji,
        source=source,
    )


def higgs_lfv_partial_width(
    yukawa_ij: complex = 0.0j,
    yukawa_ji: complex = 0.0j,
    *,
    inputs: HiggsLFVInputs | None = None,
) -> float:
    """Return ``Gamma(h -> l_i l_j)`` from effective LFV Yukawas."""

    return _partial_width(yukawa_ij, yukawa_ji, inputs=inputs)


def higgs_lfv_branching_fraction_from_yukawas(
    *,
    initial_flavor: str,
    final_flavor: str,
    yukawa_ij: complex = 0.0j,
    yukawa_ji: complex = 0.0j,
    br_limit: float | None = None,
    inputs: HiggsLFVInputs | None = None,
) -> HiggsLFVBranchingResult:
    """Return ``BR(h -> l_i l_j)`` from effective LFV Higgs Yukawas."""

    return _branching_from_yukawas(
        initial_flavor=initial_flavor,
        final_flavor=final_flavor,
        yukawa_ij=yukawa_ij,
        yukawa_ji=yukawa_ji,
        br_limit=br_limit,
        inputs=inputs,
    )


def higgs_lfv_effective_yukawa_limit(
    br_limit: float,
    *,
    inputs: HiggsLFVInputs | None = None,
) -> float:
    """Return the limit on ``sqrt(|Y_ij|^2 + |Y_ji|^2)``."""

    return _effective_yukawa_limit(br_limit, inputs=inputs)


def higgs_lfv_yukawa_proxy(
    source: object,
    *,
    initial_flavor: str,
    final_flavor: str,
) -> HiggsLFVYukawaProxy:
    """Return the documented proxy for LFV Higgs-Yukawa entries."""

    return _yukawa_proxy(
        source,
        initial_flavor=initial_flavor,
        final_flavor=final_flavor,
    )


def higgs_lfv_branching_fraction_with_proxy(
    source: object,
    *,
    initial_flavor: str,
    final_flavor: str,
    br_limit: float | None = None,
    inputs: HiggsLFVInputs | None = None,
) -> tuple[HiggsLFVBranchingResult, HiggsLFVYukawaProxy]:
    """Evaluate ``BR(h -> l_i l_j)`` with the documented Yukawa proxy."""

    return _branching_with_proxy(
        source,
        initial_flavor=initial_flavor,
        final_flavor=final_flavor,
        br_limit=br_limit,
        inputs=inputs,
    )
