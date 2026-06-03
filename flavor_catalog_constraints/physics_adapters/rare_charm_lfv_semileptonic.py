"""Adapter for LFV ``D+ -> pi+ e mu`` rare-charm semileptonic machinery.

Constraint modules import this adapter only.  The production Phase-4c path
composes the repo-owned C007 ``D -> pi`` form-factor machinery with the
rigorous off-diagonal ``e mu`` block from
``rs_semileptonic_wilsons.lfv_llqq``.  The older explicit-spurion functions
remain as legacy compatibility entry points only.

With the current diagonal charged-lepton fit, the Phase-4a LFV block is
rigorously zero at tree level, so the evaluated LFV rate is zero and
non-vetoing.  Nonzero tree-level LFV requires non-diagonal lepton structure;
loop-induced LFV is deferred.  Scalar/tensor and resonance effects remain
outside this tree-level rewire.
"""

from __future__ import annotations

from dataclasses import replace

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_charm_lfv_dilepton import (
    RARE_CHARM_LFV_DILEPTON_PROXY_V1,
    RareCharmLFVLeptonCouplingProxy,
    RareCharmLFVLeptonProxyInput,
    RareCharmLFVWilsonCoefficients,
    rare_charm_lfv_lepton_coupling_proxy as _lepton_proxy,
    rare_charm_lfv_proxy_input as _proxy_input,
)
from quarkConstraints.rs_semileptonic_wilsons import RSSemileptonicWilsonBundle
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
from .rare_charm_lfv_dilepton import (
    RARE_CHARM_LFV_TREE_LEVEL_NOTE_V1,
    rare_charm_lfv_coeff_from_rs_semileptonic,
    rare_charm_lfv_rs_semileptonic_diagnostics,
    rare_charm_lfv_wilsons_from_rs_semileptonic,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_CHARM_LFV_DILEPTON_PROXY_V1",
    "RARE_CHARM_LFV_TREE_LEVEL_NOTE_V1",
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
    "rare_charm_lfv_wilsons_from_rs_semileptonic",
    "rare_charm_dtopi_emu_default_inputs",
    "rare_charm_dtopi_emu_q2_range",
    "rare_charm_dtopi_emu_differential_branching_fraction",
    "dplus_piplus_emu_sm",
    "dplus_piplus_emu_from_couplings",
    "dplus_piplus_emu_from_rs_semileptonic_wilsons",
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


def dplus_piplus_emu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareCharmDToPiMuMuInputs | None = None,
    charge_mode: str = "eplus_muminus",
    lepton_pair: str = "e_mu",
) -> RareCharmDToPiLFVBranchingResult:
    """Evaluate ``D+ -> pi+ e mu`` from Phase-4a LFV llqq Wilsons."""

    coeff = rare_charm_lfv_coeff_from_rs_semileptonic(
        source,
        transition="c_to_u",
        lepton_pair=lepton_pair,
    )
    wilsons = rare_charm_lfv_wilsons_from_rs_semileptonic(
        source,
        matching_scale_gev=matching_scale_gev,
        lepton_pair=lepton_pair,
    )
    result = _evaluate_dplus_to_piplus_emu(
        wilsons,
        inputs=inputs,
        charge_mode=charge_mode,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        rare_charm_lfv_rs_semileptonic_diagnostics(coeff, source=source)
    )
    diagnostics["base_matching_assumption"] = coeff.matching_assumption
    diagnostics["lfv_lepton_matching_assumption"] = coeff.matching_assumption
    diagnostics["matching_assumption"] = coeff.matching_assumption
    diagnostics["c_to_u_lfv_llqq_rs_semileptonic_rewired"] = True
    diagnostics["rare_charm_lfv_proxy_reused"] = False
    diagnostics["lfv_tree_level_note"] = RARE_CHARM_LFV_TREE_LEVEL_NOTE_V1
    return replace(result, diagnostics=diagnostics)
