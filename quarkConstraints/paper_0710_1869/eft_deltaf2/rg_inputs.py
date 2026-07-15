"""Frozen LO RG defaults for the paper-owned ``Delta F = 2`` EFT slice."""

from __future__ import annotations

from qcd.constants import THRESHOLD_LIST

PAPER_0710_1869_DELTAF2_RG_INPUTS_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.rg_inputs.v1"
)
PAPER_0710_1869_DELTAF2_RG_ORDER_ID = "lo"
PAPER_0710_1869_DELTAF2_RG_SCHEME_ID = "bmu.hep-ph-0005183.ndr-ms.lo.v1"
PAPER_0710_1869_DELTAF2_RG_ADM_SOURCE_ID = (
    "bmu.hep-ph-0005183.eq2.21-plus-eq2.23.lo.v1"
)
PAPER_0710_1869_DELTAF2_RG_ALPHA_S_POLICY_ID = (
    "qcd.alpha_s.lo.beta0.thresholds_continuous.v1"
)
PAPER_0710_1869_DELTAF2_RG_THRESHOLD_POLICY_ID = (
    "qcd.thresholds.pdg2024_msbar.c-b-t.mu-eq-mh.continuous.v1"
)
PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.lr_basis_contract.v1"
)
PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID = (
    "paper_q4q5_to_bmu_lr_basis.map_frozen.audit_ready.v3"
)
PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_INACTIVE_ID = (
    "paper_q4q5_to_bmu_lr_basis.map_frozen.lr_running_inactive.v3"
)
PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ACTIVE_ID = (
    "paper_q4q5_to_bmu_lr_basis.map_frozen.lr_rg_active.custom_lr_hadronic_active.custom_lr_only_observable_active.custom_combined_observable_active.default_export_lr_capable.v7"
)
PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID = (
    PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ACTIVE_ID
)
PAPER_0710_1869_DELTAF2_RG_LR_BASIS_DIRECTION_ID = (
    "paper_q4q5_operator_order.to_bmu_q1lr_q2lr_operator_order.map_frozen.v3"
)
PAPER_0710_1869_DELTAF2_RG_LR_FIERZ_VALIDITY_STATEMENT_ID = (
    "bmu.ndr_ms.lr_sector.fierz_validity.required_for_paper_o4o5_to_bmu_q1lrq2lr_map.v1"
)
PAPER_0710_1869_DELTAF2_RG_LR_RUNNING_BRIDGE_ID = (
    "paper_o4o5_to_bmu_q1lrq2lr.lo_conjugated_wilson_running.v2"
)
PAPER_0710_1869_DELTAF2_RG_LR_BMU_BASIS_ID = (
    "bmu.lr_basis.q1lr_q2lr.ndr_ms.lo.v1"
)
PAPER_0710_1869_DELTAF2_RG_LR_BMU_OPERATOR_ORDER = (
    "Q1_LR_BMU",
    "Q2_LR_BMU",
)
PAPER_0710_1869_DELTAF2_RG_SCALE_SEMANTICS_ID = (
    "evolved_wilsons_reuse_matching_scale_field_for_current_scale.v1"
)
PAPER_0710_1869_DELTAF2_RG_SUPPORTED_OPERATOR_SUBSET_ID = (
    "paper_basis.q1_vll_q1_vrr_q4_lr_q5_lr.rg_supported.lo.v3"
)
PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV = 2.0


def default_paper_0710_1869_rg_thresholds() -> tuple[tuple[float, int, int], ...]:
    """Return the frozen PDG 2024 threshold tuple used by the LO paper RG slice."""
    return tuple(
        (float(mass), int(nf_below), int(nf_above))
        for mass, nf_below, nf_above in THRESHOLD_LIST
    )


__all__ = [
    "PAPER_0710_1869_DELTAF2_DEFAULT_MU_HAD_GEV",
    "PAPER_0710_1869_DELTAF2_RG_ADM_SOURCE_ID",
    "PAPER_0710_1869_DELTAF2_RG_ALPHA_S_POLICY_ID",
    "PAPER_0710_1869_DELTAF2_RG_INPUTS_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID",
    "PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_RG_LR_BASIS_DIRECTION_ID",
    "PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ACTIVE_ID",
    "PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_INACTIVE_ID",
    "PAPER_0710_1869_DELTAF2_RG_LR_FIERZ_VALIDITY_STATEMENT_ID",
    "PAPER_0710_1869_DELTAF2_RG_LR_BMU_BASIS_ID",
    "PAPER_0710_1869_DELTAF2_RG_LR_BMU_OPERATOR_ORDER",
    "PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID",
    "PAPER_0710_1869_DELTAF2_RG_LR_RUNNING_BRIDGE_ID",
    "PAPER_0710_1869_DELTAF2_RG_ORDER_ID",
    "PAPER_0710_1869_DELTAF2_RG_SCALE_SEMANTICS_ID",
    "PAPER_0710_1869_DELTAF2_RG_SCHEME_ID",
    "PAPER_0710_1869_DELTAF2_RG_SUPPORTED_OPERATOR_SUBSET_ID",
    "PAPER_0710_1869_DELTAF2_RG_THRESHOLD_POLICY_ID",
    "default_paper_0710_1869_rg_thresholds",
]
