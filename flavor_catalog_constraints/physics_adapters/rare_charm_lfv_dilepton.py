"""Adapter for LFV ``D0 -> e mu`` rare-charm machinery.

Constraint modules import this adapter only.  The implementation reuses the
shared :mod:`quarkConstraints.rare_charm_dilepton` inputs and quark-side
Wilson matching through :mod:`quarkConstraints.rare_charm_lfv_dilepton`.

NEEDS-HUMAN-PHYSICS: ``ParameterPoint`` does not carry a rigorous off-diagonal
``e mu`` lepton neutral-current coupling.  The proxy path accepts an explicit
lepton overlap spurion and maps it to a Z-like LFV lepton coupling while the
quark ``c -> u`` matching is reused from the rare-charm dilepton core.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_charm_dilepton import RareCharmDileptonSMInputs
from quarkConstraints.rare_charm_lfv_dilepton import (
    RARE_CHARM_LFV_DILEPTON_MODEL_V1,
    RARE_CHARM_LFV_DILEPTON_OPERATOR_CONVENTION,
    RARE_CHARM_LFV_DILEPTON_PROXY_V1,
    RareCharmLFVBranchingResult,
    RareCharmLFVLeptonCouplingProxy,
    RareCharmLFVLeptonProxyInput,
    RareCharmLFVWilsonCoefficients,
    compute_rare_charm_lfv_wilsons as _compute_lfv_wilsons,
    evaluate_d0_to_emu as _evaluate_d0_to_emu,
    rare_charm_lfv_lepton_coupling_proxy as _lepton_proxy,
    rare_charm_lfv_proxy_input as _proxy_input,
    rare_charm_lfv_sm_branching_fraction as _sm_branching_fraction,
)
from quarkConstraints.rare_charm_dilepton import default_sm_inputs as _default_sm_inputs

__all__ = [
    "QuarkMassBasisCouplings",
    "RareCharmDileptonSMInputs",
    "RARE_CHARM_LFV_DILEPTON_MODEL_V1",
    "RARE_CHARM_LFV_DILEPTON_OPERATOR_CONVENTION",
    "RARE_CHARM_LFV_DILEPTON_PROXY_V1",
    "RareCharmLFVLeptonProxyInput",
    "RareCharmLFVLeptonCouplingProxy",
    "RareCharmLFVWilsonCoefficients",
    "RareCharmLFVBranchingResult",
    "rare_charm_lfv_default_sm_inputs",
    "rare_charm_lfv_proxy_input",
    "rare_charm_lfv_lepton_coupling_proxy",
    "rare_charm_lfv_wilsons_from_couplings",
    "rare_charm_lfv_sm_branching_fraction",
    "d0_emu_from_couplings",
]


def rare_charm_lfv_default_sm_inputs() -> RareCharmDileptonSMInputs:
    """Return the shared rare-charm dilepton SM input bundle."""

    return _default_sm_inputs()


def rare_charm_lfv_proxy_input(
    left_emu_overlap: complex,
    right_emu_overlap: complex,
    m_kk_gev: float,
    *,
    source: str = "caller-supplied rare-charm e-mu lepton proxy",
) -> RareCharmLFVLeptonProxyInput:
    """Build a proxy input accepted by the rare-charm LFV adapter."""

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
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLFVLeptonCouplingProxy:
    """Return the documented e-mu lepton-coupling proxy."""

    return _lepton_proxy(source, m_kk_gev=m_kk_gev, inputs=inputs)


def rare_charm_lfv_wilsons_from_couplings(
    quark_couplings: QuarkMassBasisCouplings,
    lepton_couplings: object,
    *,
    transition: str = "c_u",
    m_kk_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLFVWilsonCoefficients:
    """Return the v1 ``c -> u e mu`` Wilson proxy."""

    return _compute_lfv_wilsons(
        quark_couplings,
        lepton_couplings,
        transition=transition,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def rare_charm_lfv_sm_branching_fraction(
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLFVBranchingResult:
    """Return the catalog SM-limit ``D0 -> e mu`` branching fraction."""

    return _sm_branching_fraction(inputs)


def d0_emu_from_couplings(
    quark_couplings: QuarkMassBasisCouplings,
    lepton_couplings: object,
    *,
    m_kk_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
) -> RareCharmLFVBranchingResult:
    """Evaluate short-distance ``BR(D0 -> e+- mu-+)`` from proxy couplings."""

    return _evaluate_d0_to_emu(
        quark_couplings,
        lepton_couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
