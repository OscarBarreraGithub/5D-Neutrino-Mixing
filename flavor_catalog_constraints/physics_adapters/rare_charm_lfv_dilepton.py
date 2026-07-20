"""Adapter for LFV ``D0 -> e mu`` rare-charm machinery.

Constraint modules import this adapter only.  The production Phase-4c path
reads the rigorous off-diagonal ``e mu`` block from
``rs_semileptonic_wilsons.lfv_llqq`` and feeds its C9/C10 inputs directly to
the rare-charm LFV core.  The older explicit-spurion functions remain as legacy
compatibility entry points only.

With the current diagonal charged-lepton fit, the Phase-4a LFV block is
rigorously zero at tree level, so the evaluated LFV rate is zero and
non-vetoing.  Nonzero tree-level LFV requires non-diagonal lepton structure;
loop-induced LFV is deferred.
"""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Any

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_charm_dilepton import RareCharmDileptonSMInputs
from quarkConstraints.rare_charm_dilepton import default_sm_inputs as _default_sm_inputs
from quarkConstraints.rare_charm_lfv_dilepton import (
    RARE_CHARM_LFV_DILEPTON_MODEL_V1,
    RARE_CHARM_LFV_DILEPTON_OPERATOR_CONVENTION,
    RARE_CHARM_LFV_DILEPTON_PROXY_V1,
    RareCharmLFVBranchingResult,
    RareCharmLFVLeptonCouplingProxy,
    RareCharmLFVLeptonProxyInput,
    RareCharmLFVWilsonCoefficients,
)
from quarkConstraints.rare_charm_lfv_dilepton import (
    compute_rare_charm_lfv_wilsons as _compute_lfv_wilsons,
)
from quarkConstraints.rare_charm_lfv_dilepton import (
    evaluate_d0_to_emu as _evaluate_d0_to_emu,
)
from quarkConstraints.rare_charm_lfv_dilepton import (
    rare_charm_lfv_lepton_coupling_proxy as _lepton_proxy,
)
from quarkConstraints.rare_charm_lfv_dilepton import (
    rare_charm_lfv_proxy_input as _proxy_input,
)
from quarkConstraints.rare_charm_lfv_dilepton import (
    rare_charm_lfv_sm_branching_fraction as _sm_branching_fraction,
)
from quarkConstraints.rs_semileptonic_wilsons import (
    RSLFVSemileptonicWilsonCoefficients,
    RSSemileptonicWilsonBundle,
)

RARE_CHARM_LFV_TREE_LEVEL_NOTE_V1 = (
    "tree-level LFV rigorous from Phase-4a lfv_llqq light-Z contacts "
    "(=0 for the diagonal charged-lepton fit); nonzero only with "
    "non-diagonal lepton structure / loop-induced LFV deferred"
)
RARE_CHARM_LFV_RS_MATCHING_STATUS_V1 = (
    "rs_semileptonic_lfv_llqq_additive_no_wilson_prefactor_reuse_no_second_"
    "1_over_M_KK_squared"
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RareCharmDileptonSMInputs",
    "RARE_CHARM_LFV_DILEPTON_MODEL_V1",
    "RARE_CHARM_LFV_DILEPTON_OPERATOR_CONVENTION",
    "RARE_CHARM_LFV_DILEPTON_PROXY_V1",
    "RARE_CHARM_LFV_TREE_LEVEL_NOTE_V1",
    "RARE_CHARM_LFV_RS_MATCHING_STATUS_V1",
    "RareCharmLFVLeptonProxyInput",
    "RareCharmLFVLeptonCouplingProxy",
    "RareCharmLFVWilsonCoefficients",
    "RareCharmLFVBranchingResult",
    "rare_charm_lfv_default_sm_inputs",
    "rare_charm_lfv_proxy_input",
    "rare_charm_lfv_lepton_coupling_proxy",
    "rare_charm_lfv_wilsons_from_couplings",
    "rare_charm_lfv_coeff_from_rs_semileptonic",
    "rare_charm_lfv_rs_semileptonic_diagnostics",
    "rare_charm_lfv_wilsons_from_rs_semileptonic",
    "rare_charm_lfv_sm_branching_fraction",
    "d0_emu_from_couplings",
    "d0_emu_from_rs_semileptonic_wilsons",
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


def rare_charm_lfv_coeff_from_rs_semileptonic(
    source: RSSemileptonicWilsonBundle,
    *,
    transition: str = "c_to_u",
    lepton_pair: str = "e_mu",
) -> RSLFVSemileptonicWilsonCoefficients:
    """Return the Phase-4a LFV ``c -> u e mu`` Wilson block."""

    try:
        coeff = source.lfv_llqq[transition][lepton_pair]
    except (AttributeError, KeyError, TypeError) as exc:
        raise ValueError(
            "rs_semileptonic_wilsons.lfv_llqq"
            f"[{transition!r}][{lepton_pair!r}] is not available"
        ) from exc
    if coeff.transition_key != "c_u":
        raise ValueError(
            f"lfv_llqq[{transition!r}][{lepton_pair!r}] "
            f"transition_key={coeff.transition_key!r}, expected 'c_u'"
        )
    if coeff.lepton_pair_key != lepton_pair:
        raise ValueError(
            f"LFV lepton pair {coeff.lepton_pair_key!r} does not match "
            f"{lepton_pair!r}"
        )
    return coeff


def rare_charm_lfv_rs_semileptonic_diagnostics(
    coeff: RSLFVSemileptonicWilsonCoefficients,
    *,
    source: RSSemileptonicWilsonBundle | None = None,
) -> dict[str, Any]:
    """Diagnostics common to the rigorous LFV rare-charm rewire."""

    return {
        "rs_semileptonic_wilsons_present": True,
        "rs_semileptonic_model_label": coeff.model_label,
        "rs_semileptonic_operator_convention": coeff.operator_convention,
        "rs_semileptonic_matching_assumption": coeff.matching_assumption,
        "rs_semileptonic_matching_status": RARE_CHARM_LFV_RS_MATCHING_STATUS_V1,
        "tree_level_matching_status": (
            "rigorous_tree_light_z_lfv_llqq_from_rs_semileptonic_wilsons"
        ),
        "lfv_tree_level_note": RARE_CHARM_LFV_TREE_LEVEL_NOTE_V1,
        "loop_lfv_status": "loop_induced_lfv_deferred",
        "rs_semileptonic_transition_key": coeff.transition_key,
        "rs_semileptonic_lepton_pair_key": coeff.lepton_pair_key,
        "rs_semileptonic_quark_sector": coeff.quark_sector,
        "rs_semileptonic_final_quark_index": int(coeff.final_quark_index),
        "rs_semileptonic_initial_quark_index": int(coeff.initial_quark_index),
        "rs_semileptonic_final_lepton_index": int(coeff.final_lepton_index),
        "rs_semileptonic_initial_lepton_index": int(coeff.initial_lepton_index),
        "rs_semileptonic_lambda_ckm_name": coeff.lambda_ckm_name,
        "rs_semileptonic_lambda_ckm": complex(coeff.lambda_ckm),
        "rs_semileptonic_contact_units": coeff.contact_units,
        "rs_semileptonic_contacts": {
            key: complex(value) for key, value in coeff.contacts.items()
        },
        "rs_semileptonic_wilson_coefficients": {
            key: complex(value) for key, value in coeff.wilsons.items()
        },
        "lfv_contact_factorized": False,
        "wilson_prefactor_reused": False,
        "second_mkk_suppression_applied": False,
        "includes_heavy_neutral_exchange": (
            None if source is None else bool(source.includes_heavy_neutral_exchange)
        ),
        "includes_heavy_neutral_lepton": (
            None if source is None else bool(source.includes_heavy_neutral_lepton)
        ),
    }


def rare_charm_lfv_wilsons_from_rs_semileptonic(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    lepton_pair: str = "e_mu",
) -> RareCharmLFVWilsonCoefficients:
    """Return rare-charm LFV core Wilson inputs from Phase-4a llqq Wilsons."""

    coeff = rare_charm_lfv_coeff_from_rs_semileptonic(
        source,
        transition="c_to_u",
        lepton_pair=lepton_pair,
    )
    scale = _diagnostic_matching_scale(matching_scale_gev)
    left_vector_contact = complex(coeff.contact_LL + coeff.contact_LR)
    right_vector_contact = complex(coeff.contact_RL + coeff.contact_RR)
    return RareCharmLFVWilsonCoefficients(
        model_label=coeff.model_label,
        base_model_label=RARE_CHARM_LFV_DILEPTON_MODEL_V1,
        operator_convention=coeff.operator_convention,
        matching_assumption=coeff.matching_assumption,
        transition_key="c_u",
        M_KK=scale,
        matching_scale=scale,
        lambda_b=complex(coeff.lambda_ckm),
        left_uc_coupling=left_vector_contact,
        right_uc_coupling=right_vector_contact,
        left_uc_overlap=0.0j,
        right_uc_overlap=0.0j,
        left_quark_delta=left_vector_contact,
        right_quark_delta=right_vector_contact,
        lepton_left_delta_emu=complex(coeff.contact_LL + coeff.contact_RL),
        lepton_right_delta_emu=complex(coeff.contact_LR + coeff.contact_RR),
        lepton_vector_delta_emu=complex(
            coeff.contact_LL + coeff.contact_LR + coeff.contact_RL + coeff.contact_RR
        ),
        lepton_axial_delta_emu=complex(
            coeff.contact_LR + coeff.contact_RR - coeff.contact_LL - coeff.contact_RL
        ),
        c9_lfv_np=complex(coeff.c9_lfv_np),
        c10_lfv_np=complex(coeff.c10_lfv_np),
        c9p_lfv_np=complex(coeff.c9p_lfv_np),
        c10p_lfv_np=complex(coeff.c10p_lfv_np),
        base_same_flavor_wilsons=None,
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


def d0_emu_from_rs_semileptonic_wilsons(
    source: RSSemileptonicWilsonBundle,
    *,
    matching_scale_gev: float | None = None,
    inputs: RareCharmDileptonSMInputs | None = None,
    lepton_pair: str = "e_mu",
) -> RareCharmLFVBranchingResult:
    """Evaluate ``D0 -> e mu`` from Phase-4a LFV llqq Wilsons."""

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
    result = _evaluate_d0_to_emu(wilsons, inputs=inputs)
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        rare_charm_lfv_rs_semileptonic_diagnostics(coeff, source=source)
    )
    diagnostics["base_matching_assumption"] = coeff.matching_assumption
    diagnostics["matching_assumption"] = coeff.matching_assumption
    diagnostics["c_to_u_lfv_llqq_rs_semileptonic_rewired"] = True
    diagnostics["rare_charm_lfv_proxy_reused"] = False
    return replace(result, diagnostics=diagnostics)


def _diagnostic_matching_scale(matching_scale_gev: float | None) -> float:
    if matching_scale_gev is None:
        return 0.0
    number = float(matching_scale_gev)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError("matching_scale_gev must be positive and finite")
    return number
