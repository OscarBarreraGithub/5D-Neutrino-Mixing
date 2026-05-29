"""Adapter for LFV ``D+ -> pi+ e mu`` rare-charm semileptonic machinery.

Constraint modules import this adapter only.  The implementation composes the
repo-owned C007 ``D -> pi`` form-factor machinery with the C006 documented
``e mu`` LFV Wilson proxy in :mod:`quarkConstraints.rare_charm_lfv_semileptonic`.

NEEDS-HUMAN-PHYSICS: ``ParameterPoint`` does not carry the full off-diagonal
charged-lepton neutral-current couplings or the scalar/tensor/resonance
information needed for a rigorous RS recast.  The v1 observable is a
short-distance full-q2 proxy.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_charm_lfv_dilepton import (
    RARE_CHARM_LFV_DILEPTON_PROXY_V1,
    RareCharmLFVLeptonCouplingProxy,
    RareCharmLFVLeptonProxyInput,
    RareCharmLFVWilsonCoefficients,
    rare_charm_lfv_lepton_coupling_proxy as _lepton_proxy,
    rare_charm_lfv_proxy_input as _proxy_input,
)
from quarkConstraints.rare_charm_lfv_semileptonic import (
    RARE_CHARM_DTOPI_EMU_MODEL_V1,
    RARE_CHARM_DTOPI_EMU_OPERATOR_CONVENTION,
    RARE_CHARM_DTOPI_EMU_PARAMETRIZATION_CITATION,
    RARE_CHARM_DTOPI_EMU_PROXY_V1,
    RARE_CHARM_DTOPI_EMU_Q2_TREATMENT_V1,
    RareCharmDToPiLFVBranchingResult,
    RareCharmDToPiMuMuInputs,
    default_dtopi_emu_inputs as _default_dtopi_emu_inputs,
    dtopi_emu_differential_branching_fraction as _differential_branching,
    dtopi_emu_q2_range as _q2_range,
    dtopi_emu_sm as _dtopi_emu_sm,
    evaluate_dplus_to_piplus_emu as _evaluate_dplus_to_piplus_emu,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_CHARM_LFV_DILEPTON_PROXY_V1",
    "RARE_CHARM_DTOPI_EMU_MODEL_V1",
    "RARE_CHARM_DTOPI_EMU_OPERATOR_CONVENTION",
    "RARE_CHARM_DTOPI_EMU_PARAMETRIZATION_CITATION",
    "RARE_CHARM_DTOPI_EMU_PROXY_V1",
    "RARE_CHARM_DTOPI_EMU_Q2_TREATMENT_V1",
    "RareCharmLFVLeptonProxyInput",
    "RareCharmLFVLeptonCouplingProxy",
    "RareCharmLFVWilsonCoefficients",
    "RareCharmDToPiMuMuInputs",
    "RareCharmDToPiLFVBranchingResult",
    "rare_charm_lfv_proxy_input",
    "rare_charm_lfv_lepton_coupling_proxy",
    "rare_charm_dtopi_emu_default_inputs",
    "rare_charm_dtopi_emu_q2_range",
    "rare_charm_dtopi_emu_differential_branching_fraction",
    "dplus_piplus_emu_sm",
    "dplus_piplus_emu_from_couplings",
]


def rare_charm_lfv_proxy_input(
    left_emu_overlap: complex,
    right_emu_overlap: complex,
    m_kk_gev: float,
    *,
    source: str = "caller-supplied rare-charm e-mu lepton proxy",
) -> RareCharmLFVLeptonProxyInput:
    """Build a proxy input accepted by the LFV rare-charm adapters."""

    return _proxy_input(
        left_emu_overlap,
        right_emu_overlap,
        m_kk_gev,
        source=source,
    )


def rare_charm_lfv_lepton_coupling_proxy(
    source: object,
    *,
    m_kk_gev: float | None = None,
) -> RareCharmLFVLeptonCouplingProxy:
    """Return the documented e-mu lepton-coupling proxy."""

    return _lepton_proxy(source, m_kk_gev=m_kk_gev)


def rare_charm_dtopi_emu_default_inputs() -> RareCharmDToPiMuMuInputs:
    """Return the shared ``D+ -> pi+`` short-distance input bundle."""

    return _default_dtopi_emu_inputs()


def rare_charm_dtopi_emu_q2_range(
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> tuple[float, float]:
    """Return the kinematic ``q2`` range for ``D+ -> pi+ e mu``."""

    return _q2_range(inputs)


def rare_charm_dtopi_emu_differential_branching_fraction(
    q2_gev2: float,
    *,
    c9_lfv_semileptonic: complex,
    c10_lfv_semileptonic: complex,
    inputs: RareCharmDToPiMuMuInputs | None = None,
    charge_mode: str = "eplus_muminus",
) -> float:
    """Evaluate ``dBR_SD(D+ -> pi+ e mu) / dq2``."""

    return _differential_branching(
        q2_gev2,
        c9_lfv_semileptonic=c9_lfv_semileptonic,
        c10_lfv_semileptonic=c10_lfv_semileptonic,
        inputs=inputs,
        charge_mode=charge_mode,
    )


def dplus_piplus_emu_sm(
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> RareCharmDToPiLFVBranchingResult:
    """Return the catalog SM-limit ``D+ -> pi+ e mu`` branching fraction."""

    return _dtopi_emu_sm(inputs)


def dplus_piplus_emu_from_couplings(
    quark_couplings: QuarkMassBasisCouplings,
    lepton_couplings: object,
    *,
    m_kk_gev: float | None = None,
    inputs: RareCharmDToPiMuMuInputs | None = None,
    charge_mode: str = "eplus_muminus",
) -> RareCharmDToPiLFVBranchingResult:
    """Evaluate smooth ``BR_SD(D+ -> pi+ e mu)`` from proxy couplings."""

    return _evaluate_dplus_to_piplus_emu(
        quark_couplings,
        lepton_couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
        charge_mode=charge_mode,
    )
