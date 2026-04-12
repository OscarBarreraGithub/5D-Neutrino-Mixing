"""Paper-owned observables for the 0710.1869 ``Delta F = 2`` slice."""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

from ..conventions import PAPER_0710_1869_MODE_ID, PAPER_0710_1869_PAPER_ID
from ..validation import require_known_schema_id, require_member, require_nonempty_identifier
from .hadronic import (
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID,
    PAPER_0710_1869_DELTAF2_BD_MATRIX_ELEMENT_FORMULA_ID,
    PAPER_0710_1869_DELTAF2_BD_PARITY_RELATION_ID,
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID,
    PAPER_0710_1869_DELTAF2_BS_MATRIX_ELEMENT_FORMULA_ID,
    PAPER_0710_1869_DELTAF2_BS_PARITY_RELATION_ID,
    PAPER_0710_1869_DELTAF2_D0_HADRONIC_SYSTEM_ID,
    PAPER_0710_1869_DELTAF2_D0_MATRIX_ELEMENT_FORMULA_ID,
    PAPER_0710_1869_DELTAF2_D0_PARITY_RELATION_ID,
    PAPER_0710_1869_DELTAF2_HADRONIC_CUSTOM_PROVENANCE_MODE_ID,
    PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID,
    PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID,
    PAPER_0710_1869_DELTAF2_KAON_LR_CHIRAL_RATIO_DEFINITION_ID,
    PAPER_0710_1869_DELTAF2_KAON_LR_Q4_MATRIX_ELEMENT_FORMULA_ID,
    PAPER_0710_1869_DELTAF2_KAON_LR_Q5_MATRIX_ELEMENT_FORMULA_ID,
    PAPER_0710_1869_DELTAF2_KAON_MATRIX_ELEMENT_FORMULA_ID,
    PAPER_0710_1869_DELTAF2_KAON_PARITY_RELATION_ID,
    PAPER_0710_1869_DELTAF2_LR_HADRONIC_CUSTOM_PROVENANCE_MODE_ID,
    Paper07101869BMesonHadronicBundle,
    Paper07101869D0HadronicBundle,
    Paper07101869KaonHadronicBundle,
    Paper07101869KaonLRHadronicInputs,
    default_paper_0710_1869_kaon_hadronic_bundle,
)
from .matching_kkgluon import default_paper_0710_1869_kaon_matching
from .operators import (
    PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID,
    PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID,
    PAPER_0710_1869_DELTAF2_Q4_LR,
    PAPER_0710_1869_DELTAF2_Q5_LR,
)
from .rg import (
    PAPER_0710_1869_DELTAF2_RG_WILSON_SCHEMA_ID,
    Paper07101869DeltaF2RGResult,
    Paper07101869DeltaF2RGWilsonSnapshot,
    run_deltaf2_wilsons_lo,
)
from .rg_inputs import (
    PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID,
    PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID,
    PAPER_0710_1869_DELTAF2_RG_SCHEME_ID,
)

PAPER_0710_1869_DELTAF2_KAON_NP_OBSERVABLE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_np_observable_result.v1"
)
PAPER_0710_1869_DELTAF2_KAON_NP_OBSERVABLE_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_np_observable_summary.v1"
)
PAPER_0710_1869_DELTAF2_KAON_NP_OBSERVABLE_SCOPE_ID = "kaon.np_only.m12.v1"
PAPER_0710_1869_DELTAF2_KAON_NP_INTERPRETATION_ID = "kaon.np_only.v1"
PAPER_0710_1869_DELTAF2_KAON_M12_OBSERVABLE_ID = "M12_K_NP"
PAPER_0710_1869_DELTAF2_KAON_DELTA_M_OBSERVABLE_ID = "Delta_m_K_NP"
PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_OBSERVABLE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_lr_only_observable_result.v1"
)
PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_OBSERVABLE_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_lr_only_observable_summary.v1"
)
PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_OBSERVABLE_SCOPE_ID = "kaon.lr_only.custom.m12.v1"
PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_INTERPRETATION_ID = "kaon.lr_only.custom.v1"
PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_M12_OBSERVABLE_ID = "M12_K_LR_NP"
PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_DELTA_M_OBSERVABLE_ID = "Delta_m_K_LR_NP"
PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_OBSERVABLE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_custom_total_observable_result.v1"
)
PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_OBSERVABLE_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_custom_total_observable_summary.v1"
)
PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_OBSERVABLE_SCOPE_ID = (
    "kaon.np_only.custom_total.q1_plus_lr.v1"
)
PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_INTERPRETATION_ID = (
    "kaon.np_only.custom_total.q1_plus_lr.v1"
)
PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_M12_OBSERVABLE_ID = (
    "M12_K_NP_CUSTOM_TOTAL"
)
PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_DELTA_M_OBSERVABLE_ID = (
    "Delta_m_K_NP_CUSTOM_TOTAL"
)
PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_OBSERVABLE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.bd_custom_q1_observable_result.v1"
)
PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_OBSERVABLE_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.bd_custom_q1_observable_summary.v1"
)
PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_OBSERVABLE_SCOPE_ID = "bd.np_only.custom_q1.m12.v1"
PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_INTERPRETATION_ID = "bd.np_only.custom_q1.v1"
PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_M12_OBSERVABLE_ID = "M12_Bd_NP_CUSTOM_Q1"
PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_DELTA_M_OBSERVABLE_ID = (
    "Delta_m_Bd_NP_CUSTOM_Q1"
)
PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_OBSERVABLE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.bs_custom_q1_observable_result.v1"
)
PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_OBSERVABLE_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.bs_custom_q1_observable_summary.v1"
)
PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_OBSERVABLE_SCOPE_ID = "bs.np_only.custom_q1.m12.v1"
PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_INTERPRETATION_ID = "bs.np_only.custom_q1.v1"
PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_M12_OBSERVABLE_ID = "M12_Bs_NP_CUSTOM_Q1"
PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_DELTA_M_OBSERVABLE_ID = (
    "Delta_m_Bs_NP_CUSTOM_Q1"
)
PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_OBSERVABLE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.d0_custom_q1_observable_result.v1"
)
PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_OBSERVABLE_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.d0_custom_q1_observable_summary.v1"
)
PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_OBSERVABLE_SCOPE_ID = "d0.np_only.custom_q1.m12.v1"
PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_INTERPRETATION_ID = "d0.np_only.custom_q1.v1"
PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_M12_OBSERVABLE_ID = "M12_D0_NP"
PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_DELTA_M_OBSERVABLE_ID = "Delta_m_D0_NP"

_ZERO_TOLERANCE = 1e-30


def _require_finite_complex(name: str, value: complex) -> complex:
    complex_value = complex(value)
    if not math.isfinite(complex_value.real) or not math.isfinite(complex_value.imag):
        raise ValueError(f"{name} must be finite")
    real_part = 0.0 if complex_value.real == 0.0 else float(complex_value.real)
    imag_part = 0.0 if complex_value.imag == 0.0 else float(complex_value.imag)
    return complex(real_part, imag_part)


def _complex_as_dict(value: complex) -> dict[str, float]:
    complex_value = _require_finite_complex("value", value)
    return {
        "real": float(complex_value.real),
        "imag": float(complex_value.imag),
    }


def _resolve_rg_wilsons(
    rg_result_or_wilsons: Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot,
) -> Paper07101869DeltaF2RGWilsonSnapshot:
    if isinstance(rg_result_or_wilsons, Paper07101869DeltaF2RGResult):
        return rg_result_or_wilsons.wilsons
    if isinstance(rg_result_or_wilsons, Paper07101869DeltaF2RGWilsonSnapshot):
        return rg_result_or_wilsons
    raise ValueError(
        "rg_result_or_wilsons must be a Paper07101869DeltaF2RGResult or "
        "Paper07101869DeltaF2RGWilsonSnapshot"
    )


def _require_hadronic_compatibility(
    *,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot,
    hadronic_bundle: Paper07101869KaonHadronicBundle,
) -> None:
    hadronic_contract = hadronic_bundle.contract
    if hadronic_contract.system_id != hadronic_bundle.system_id:
        raise ValueError("hadronic contract system_id must match the kaon bundle system_id")
    if hadronic_contract.operator_basis_id != hadronic_bundle.operator_basis_id:
        raise ValueError(
            "hadronic contract operator_basis_id must match the frozen kaon bundle basis"
        )
    if (
        hadronic_contract.operator_normalization_id
        != hadronic_bundle.operator_normalization_id
    ):
        raise ValueError(
            "hadronic contract operator_normalization_id must match the frozen kaon bundle "
            "normalization"
        )
    if (
        hadronic_contract.renormalization_scheme_id
        != hadronic_bundle.renormalization_scheme_id
    ):
        raise ValueError(
            "hadronic contract renormalization_scheme_id must match the frozen kaon bundle "
            "scheme"
        )
    if (
        hadronic_contract.matrix_element_formula_id
        != hadronic_bundle.matrix_element_formula_id
    ):
        raise ValueError(
            "hadronic contract matrix_element_formula_id must match the frozen kaon formula"
        )
    if (
        hadronic_contract.hamiltonian_convention_id
        != hadronic_bundle.hamiltonian_convention_id
    ):
        raise ValueError(
            "hadronic contract hamiltonian_convention_id must match the frozen paper convention"
        )
    if hadronic_contract.parity_relation_id != hadronic_bundle.parity_relation_id:
        raise ValueError(
            "hadronic contract parity_relation_id must match the frozen VLL/VRR parity relation"
        )
    if not math.isclose(
        hadronic_contract.mu_had_GeV,
        hadronic_bundle.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError("hadronic contract mu_had_GeV must match the kaon bundle mu_had_GeV")
    if not math.isclose(
        hadronic_contract.evaluation_scale_GeV,
        hadronic_bundle.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "hadronic contract evaluation_scale_GeV must match the kaon bundle mu_had_GeV"
        )
    if hadronic_bundle.system_id != PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID:
        raise ValueError("hadronic bundle must be kaon-only for the PR5a observable path")
    if hadronic_bundle.operator_basis_id != PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID:
        raise ValueError("hadronic bundle operator_basis_id must match the frozen paper basis")
    if (
        hadronic_bundle.operator_normalization_id
        != PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    ):
        raise ValueError(
            "hadronic bundle operator_normalization_id must match the frozen paper normalization"
        )
    if hadronic_bundle.renormalization_scheme_id != PAPER_0710_1869_DELTAF2_RG_SCHEME_ID:
        raise ValueError(
            "hadronic bundle renormalization_scheme_id must match the frozen RG scheme"
        )
    if (
        hadronic_bundle.hamiltonian_convention_id
        != PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID
    ):
        raise ValueError(
            "hadronic bundle hamiltonian_convention_id must match the frozen paper convention"
        )
    if (
        hadronic_bundle.matrix_element_formula_id
        != PAPER_0710_1869_DELTAF2_KAON_MATRIX_ELEMENT_FORMULA_ID
    ):
        raise ValueError(
            "hadronic bundle matrix_element_formula_id must match the frozen kaon formula"
        )
    if hadronic_bundle.parity_relation_id != PAPER_0710_1869_DELTAF2_KAON_PARITY_RELATION_ID:
        raise ValueError(
            "hadronic bundle parity_relation_id must match the frozen VLL/VRR parity relation "
            "before a shared Q1 matrix element can be used"
        )
    if wilsons.system_id != PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID:
        raise ValueError("kaon observable evaluation requires kaon RG-tagged wilsons")
    if wilsons.sector_id != "down":
        raise ValueError("kaon observable evaluation requires the down-sector kaon system")
    if wilsons.operator_basis_id != hadronic_bundle.operator_basis_id:
        raise ValueError(
            "hadronic bundle operator_basis_id is incompatible with the RG Wilson snapshot"
        )
    if wilsons.operator_normalization_id != hadronic_bundle.operator_normalization_id:
        raise ValueError(
            "hadronic bundle operator_normalization_id is incompatible with the RG Wilson snapshot"
        )
    if wilsons.renormalization_scheme_id != hadronic_bundle.renormalization_scheme_id:
        raise ValueError(
            "hadronic bundle renormalization_scheme_id is incompatible with the RG Wilson snapshot"
        )
    if not math.isclose(
        wilsons.matching_scale_GeV,
        hadronic_bundle.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "hadronic bundle mu_had_GeV must match the evolved Wilson evaluation scale"
        )
    if abs(wilsons.q4_lr) > _ZERO_TOLERANCE or abs(wilsons.q5_lr) > _ZERO_TOLERANCE:
        raise ValueError(
            "kaon NP observables currently support only Q1_VLL/Q1_VRR on the "
            "default/exported path; "
            "non-zero Q4_LR/Q5_LR remain blocked by LR contract "
            f"{PAPER_0710_1869_DELTAF2_RG_LR_BASIS_CONTRACT_ID!r} under status "
            f"{PAPER_0710_1869_DELTAF2_RG_LR_BASIS_STATUS_ID!r}. "
            "Use the custom LR-only observable surface with exact LR hadronic alignment "
            "instead."
        )


def _require_lr_hadronic_compatibility(
    *,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot,
    hadronic_inputs: Paper07101869KaonLRHadronicInputs,
) -> None:
    hadronic_contract = hadronic_inputs.contract
    if hadronic_contract.system_id != hadronic_inputs.system_id:
        raise ValueError("hadronic contract system_id must match the LR hadronic system_id")
    if hadronic_contract.operator_basis_id != hadronic_inputs.operator_basis_id:
        raise ValueError(
            "hadronic contract operator_basis_id must match the LR hadronic paper basis"
        )
    if (
        hadronic_contract.operator_normalization_id
        != hadronic_inputs.operator_normalization_id
    ):
        raise ValueError(
            "hadronic contract operator_normalization_id must match the LR hadronic "
            "normalization"
        )
    if (
        hadronic_contract.renormalization_scheme_id
        != hadronic_inputs.renormalization_scheme_id
    ):
        raise ValueError(
            "hadronic contract renormalization_scheme_id must match the LR hadronic "
            "renormalization_scheme_id"
        )
    if (
        hadronic_contract.hamiltonian_convention_id
        != hadronic_inputs.hamiltonian_convention_id
    ):
        raise ValueError(
            "hadronic contract hamiltonian_convention_id must match the LR hadronic "
            "hamiltonian convention"
        )
    if (
        hadronic_contract.q4_matrix_element_formula_id
        != hadronic_inputs.q4_matrix_element_formula_id
    ):
        raise ValueError(
            "hadronic contract q4_matrix_element_formula_id must match the LR hadronic "
            "Q4 matrix-element formula"
        )
    if (
        hadronic_contract.q5_matrix_element_formula_id
        != hadronic_inputs.q5_matrix_element_formula_id
    ):
        raise ValueError(
            "hadronic contract q5_matrix_element_formula_id must match the LR hadronic "
            "Q5 matrix-element formula"
        )
    if (
        hadronic_contract.chiral_ratio_definition_id
        != hadronic_inputs.chiral_ratio_definition_id
    ):
        raise ValueError(
            "hadronic contract chiral_ratio_definition_id must match the LR hadronic "
            "custom R_chi definition"
        )
    if not math.isclose(
        hadronic_contract.mu_had_GeV,
        hadronic_inputs.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError("hadronic contract mu_had_GeV must match the LR hadronic mu_had_GeV")
    if not math.isclose(
        hadronic_contract.evaluation_scale_GeV,
        hadronic_inputs.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "hadronic contract evaluation_scale_GeV must match the LR hadronic mu_had_GeV"
        )
    if hadronic_inputs.system_id != PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID:
        raise ValueError("LR-only observable evaluation requires kaon hadronic inputs")
    if hadronic_inputs.operator_basis_id != PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID:
        raise ValueError(
            "LR hadronic operator_basis_id must match the frozen paper DeltaF=2 basis"
        )
    if (
        hadronic_inputs.operator_normalization_id
        != PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    ):
        raise ValueError(
            "LR hadronic operator_normalization_id must match the frozen paper "
            "normalization"
        )
    if (
        hadronic_inputs.hamiltonian_convention_id
        != PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID
    ):
        raise ValueError(
            "LR hadronic hamiltonian_convention_id must match the frozen paper "
            "convention"
        )
    require_member(
        "hadronic_inputs.input_provenance_mode_id",
        hadronic_inputs.input_provenance_mode_id,
        (PAPER_0710_1869_DELTAF2_LR_HADRONIC_CUSTOM_PROVENANCE_MODE_ID,),
    )
    if (
        hadronic_inputs.q4_matrix_element_formula_id
        != PAPER_0710_1869_DELTAF2_KAON_LR_Q4_MATRIX_ELEMENT_FORMULA_ID
    ):
        raise ValueError(
            "LR hadronic Q4 matrix-element formula must match the frozen BV 2004 contract"
        )
    if (
        hadronic_inputs.q5_matrix_element_formula_id
        != PAPER_0710_1869_DELTAF2_KAON_LR_Q5_MATRIX_ELEMENT_FORMULA_ID
    ):
        raise ValueError(
            "LR hadronic Q5 matrix-element formula must match the frozen BV 2004 contract"
        )
    if (
        hadronic_inputs.chiral_ratio_definition_id
        != PAPER_0710_1869_DELTAF2_KAON_LR_CHIRAL_RATIO_DEFINITION_ID
    ):
        raise ValueError(
            "LR hadronic chiral_ratio_definition_id must match the frozen custom R_chi "
            "definition"
        )
    if wilsons.system_id != PAPER_0710_1869_DELTAF2_KAON_HADRONIC_SYSTEM_ID:
        raise ValueError("LR-only observable evaluation requires kaon RG-tagged wilsons")
    if wilsons.sector_id != "down":
        raise ValueError("LR-only observable evaluation requires the down-sector kaon system")
    if wilsons.operator_basis_id != hadronic_inputs.operator_basis_id:
        raise ValueError(
            "LR hadronic operator_basis_id is incompatible with the RG Wilson snapshot"
        )
    if wilsons.operator_normalization_id != hadronic_inputs.operator_normalization_id:
        raise ValueError(
            "LR hadronic operator_normalization_id is incompatible with the RG Wilson snapshot"
        )
    if wilsons.renormalization_scheme_id != hadronic_inputs.renormalization_scheme_id:
        raise ValueError(
            "LR hadronic renormalization_scheme_id is incompatible with the RG Wilson snapshot"
        )
    if not math.isclose(
        wilsons.matching_scale_GeV,
        hadronic_inputs.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "LR hadronic mu_had_GeV must match the evolved Wilson evaluation scale"
        )
    if not math.isclose(
        wilsons.matching_scale_GeV,
        hadronic_contract.evaluation_scale_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "LR hadronic evaluation_scale_GeV must match the evolved Wilson evaluation scale"
        )


def _require_custom_q1_hadronic_compatibility(
    *,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot,
    hadronic_bundle: Paper07101869KaonHadronicBundle,
) -> None:
    q1_only_wilsons = replace(wilsons, q4_lr=0.0 + 0.0j, q5_lr=0.0 + 0.0j)
    _require_hadronic_compatibility(
        wilsons=q1_only_wilsons,
        hadronic_bundle=hadronic_bundle,
    )
    require_member(
        "hadronic_bundle.input_provenance_mode_id",
        hadronic_bundle.input_provenance_mode_id,
        (PAPER_0710_1869_DELTAF2_HADRONIC_CUSTOM_PROVENANCE_MODE_ID,),
    )


def _require_custom_combined_alignment(
    *,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot,
    q1_hadronic_bundle: Paper07101869KaonHadronicBundle,
    lr_hadronic_inputs: Paper07101869KaonLRHadronicInputs,
) -> None:
    alignment_fields = (
        "system_id",
        "operator_basis_id",
        "operator_normalization_id",
        "renormalization_scheme_id",
    )
    for field_name in alignment_fields:
        q1_value = getattr(q1_hadronic_bundle, field_name)
        lr_value = getattr(lr_hadronic_inputs, field_name)
        if q1_value != lr_value:
            raise ValueError(
                f"custom combined kaon observable requires Q1 and LR hadronic alignment on "
                f"{field_name}"
            )
        if getattr(wilsons, field_name) != q1_value:
            raise ValueError(
                f"custom combined kaon observable requires RG Wilson alignment on "
                f"{field_name}"
            )
    if (
        q1_hadronic_bundle.hamiltonian_convention_id
        != lr_hadronic_inputs.hamiltonian_convention_id
    ):
        raise ValueError(
            "custom combined kaon observable requires Q1 and LR hadronic alignment on "
            "hamiltonian_convention_id"
        )
    if (
        q1_hadronic_bundle.hamiltonian_convention_id
        != PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID
    ):
        raise ValueError(
            "custom combined kaon observable requires the frozen paper "
            "hamiltonian_convention_id"
        )
    if not math.isclose(
        q1_hadronic_bundle.mu_had_GeV,
        lr_hadronic_inputs.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "custom combined kaon observable requires Q1 and LR hadronic alignment on mu_had_GeV"
        )
    if not math.isclose(
        q1_hadronic_bundle.contract.evaluation_scale_GeV,
        lr_hadronic_inputs.contract.evaluation_scale_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "custom combined kaon observable requires Q1 and LR hadronic alignment on "
            "evaluation_scale_GeV"
        )
    if not math.isclose(
        wilsons.matching_scale_GeV,
        q1_hadronic_bundle.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "custom combined kaon observable requires the evolved Wilson evaluation scale "
            "to match the hadronic mu_had_GeV"
        )
    if not math.isclose(
        wilsons.matching_scale_GeV,
        q1_hadronic_bundle.contract.evaluation_scale_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "custom combined kaon observable requires the evolved Wilson evaluation scale "
            "to match the hadronic evaluation_scale_GeV"
        )


def _compute_q1_only_m12_value(
    *,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot,
    hadronic_bundle: Paper07101869KaonHadronicBundle,
) -> complex:
    matrix_element = hadronic_bundle.q1_matrix_element_GeV4
    return _require_finite_complex(
        "M12_K_NP_GeV",
        (wilsons.q1_vll + wilsons.q1_vrr) * matrix_element / (2.0 * hadronic_bundle.m_K0_GeV),
    )


_CUSTOM_B_Q1_FORMULA_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BD_MATRIX_ELEMENT_FORMULA_ID
    ),
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BS_MATRIX_ELEMENT_FORMULA_ID
    ),
}
_CUSTOM_B_Q1_PARITY_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: PAPER_0710_1869_DELTAF2_BD_PARITY_RELATION_ID,
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: PAPER_0710_1869_DELTAF2_BS_PARITY_RELATION_ID,
}
_CUSTOM_B_Q1_SCHEMA_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_OBSERVABLE_SCHEMA_ID
    ),
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_OBSERVABLE_SCHEMA_ID
    ),
}
_CUSTOM_B_Q1_SUMMARY_SCHEMA_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_OBSERVABLE_SUMMARY_SCHEMA_ID
    ),
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_OBSERVABLE_SUMMARY_SCHEMA_ID
    ),
}
_CUSTOM_B_Q1_SCOPE_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_OBSERVABLE_SCOPE_ID
    ),
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_OBSERVABLE_SCOPE_ID
    ),
}
_CUSTOM_B_Q1_INTERPRETATION_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_INTERPRETATION_ID
    ),
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_INTERPRETATION_ID
    ),
}
_CUSTOM_B_Q1_M12_ID_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_M12_OBSERVABLE_ID
    ),
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_M12_OBSERVABLE_ID
    ),
}
_CUSTOM_B_Q1_DELTA_M_ID_BY_SYSTEM = {
    PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_DELTA_M_OBSERVABLE_ID
    ),
    PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID: (
        PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_DELTA_M_OBSERVABLE_ID
    ),
}


def _require_d0_q1_hadronic_compatibility(
    *,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot,
    hadronic_bundle: Paper07101869D0HadronicBundle,
) -> None:
    hadronic_contract = hadronic_bundle.contract
    if hadronic_contract.system_id != hadronic_bundle.system_id:
        raise ValueError("hadronic contract system_id must match the D0 bundle system_id")
    if hadronic_contract.operator_basis_id != hadronic_bundle.operator_basis_id:
        raise ValueError("hadronic contract operator_basis_id must match the D0 bundle basis")
    if hadronic_contract.operator_normalization_id != hadronic_bundle.operator_normalization_id:
        raise ValueError(
            "hadronic contract operator_normalization_id must match the D0 bundle normalization"
        )
    if (
        hadronic_contract.renormalization_scheme_id
        != hadronic_bundle.renormalization_scheme_id
    ):
        raise ValueError(
            "hadronic contract renormalization_scheme_id must match the D0 bundle scheme"
        )
    if hadronic_contract.hamiltonian_convention_id != hadronic_bundle.hamiltonian_convention_id:
        raise ValueError(
            "hadronic contract hamiltonian_convention_id must match the frozen paper convention"
        )
    if (
        hadronic_contract.matrix_element_formula_id
        != hadronic_bundle.matrix_element_formula_id
    ):
        raise ValueError(
            "hadronic contract matrix_element_formula_id must match "
            "the D0 Q1 matrix-element formula"
        )
    if hadronic_contract.parity_relation_id != hadronic_bundle.parity_relation_id:
        raise ValueError(
            "hadronic contract parity_relation_id must match the D0 VLL/VRR parity relation"
        )
    if not math.isclose(
        hadronic_contract.mu_had_GeV,
        hadronic_bundle.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError("hadronic contract mu_had_GeV must match the D0 bundle mu_had_GeV")
    if not math.isclose(
        hadronic_contract.evaluation_scale_GeV,
        hadronic_bundle.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "hadronic contract evaluation_scale_GeV must match the D0 bundle mu_had_GeV"
        )
    if hadronic_bundle.system_id != PAPER_0710_1869_DELTAF2_D0_HADRONIC_SYSTEM_ID:
        raise ValueError("hadronic bundle must be D0-only for the custom D0 observable path")
    if hadronic_bundle.operator_basis_id != PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID:
        raise ValueError("hadronic bundle operator_basis_id must match the frozen paper basis")
    if (
        hadronic_bundle.operator_normalization_id
        != PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    ):
        raise ValueError(
            "hadronic bundle operator_normalization_id must match the frozen paper normalization"
        )
    if (
        hadronic_bundle.hamiltonian_convention_id
        != PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID
    ):
        raise ValueError(
            "hadronic bundle hamiltonian_convention_id must match the frozen paper convention"
        )
    if (
        hadronic_bundle.matrix_element_formula_id
        != PAPER_0710_1869_DELTAF2_D0_MATRIX_ELEMENT_FORMULA_ID
    ):
        raise ValueError(
            "hadronic bundle matrix_element_formula_id must match the frozen D0 Q1 formula"
        )
    if hadronic_bundle.parity_relation_id != PAPER_0710_1869_DELTAF2_D0_PARITY_RELATION_ID:
        raise ValueError(
            "hadronic bundle parity_relation_id must match the frozen D0 VLL/VRR parity relation"
        )
    require_member(
        "hadronic_bundle.input_provenance_mode_id",
        hadronic_bundle.input_provenance_mode_id,
        (PAPER_0710_1869_DELTAF2_HADRONIC_CUSTOM_PROVENANCE_MODE_ID,),
    )
    if wilsons.system_id != PAPER_0710_1869_DELTAF2_D0_HADRONIC_SYSTEM_ID:
        raise ValueError("custom D0 observable evaluation requires D0 RG-tagged wilsons")
    if wilsons.sector_id != "up":
        raise ValueError("custom D0 observable evaluation requires the up-sector system")
    if tuple(wilsons.generations) != (0, 1):
        raise ValueError(
            "custom D0 observable evaluation requires generations=(0, 1) "
            "for the D0 slice"
        )
    if wilsons.operator_basis_id != hadronic_bundle.operator_basis_id:
        raise ValueError(
            "hadronic bundle operator_basis_id is incompatible with the RG Wilson snapshot"
        )
    if wilsons.operator_normalization_id != hadronic_bundle.operator_normalization_id:
        raise ValueError(
            "hadronic bundle operator_normalization_id is incompatible with the RG Wilson snapshot"
        )
    if wilsons.renormalization_scheme_id != hadronic_bundle.renormalization_scheme_id:
        raise ValueError(
            "hadronic bundle renormalization_scheme_id is incompatible with the RG Wilson snapshot"
        )
    if not math.isclose(
        wilsons.matching_scale_GeV,
        hadronic_bundle.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "hadronic bundle mu_had_GeV must match the evolved Wilson evaluation scale"
        )
    if abs(wilsons.q4_lr) > _ZERO_TOLERANCE or abs(wilsons.q5_lr) > _ZERO_TOLERANCE:
        raise ValueError(
            "custom D0 Q1 observables currently support only Q1_VLL/Q1_VRR; non-zero "
            "Q4_LR/Q5_LR remain out of scope for the D0 slice."
        )


def _require_custom_b_q1_hadronic_compatibility(
    *,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot,
    hadronic_bundle: Paper07101869BMesonHadronicBundle,
    expected_system_id: str,
) -> None:
    require_member("expected_system_id", expected_system_id, tuple(_CUSTOM_B_Q1_SCOPE_BY_SYSTEM))
    if hadronic_bundle.system_id != expected_system_id:
        raise ValueError(
            f"custom {expected_system_id} observable evaluation requires "
            f"{expected_system_id} hadronic inputs"
        )
    hadronic_contract = hadronic_bundle.contract
    if hadronic_contract.system_id != hadronic_bundle.system_id:
        raise ValueError("hadronic contract system_id must match the B-system bundle system_id")
    if hadronic_contract.operator_basis_id != hadronic_bundle.operator_basis_id:
        raise ValueError(
            "hadronic contract operator_basis_id must match the B-system bundle basis"
        )
    if hadronic_contract.operator_normalization_id != hadronic_bundle.operator_normalization_id:
        raise ValueError(
            "hadronic contract operator_normalization_id must match the B-system "
            "bundle normalization"
        )
    if hadronic_contract.renormalization_scheme_id != hadronic_bundle.renormalization_scheme_id:
        raise ValueError(
            "hadronic contract renormalization_scheme_id must match the B-system bundle scheme"
        )
    if hadronic_contract.hamiltonian_convention_id != hadronic_bundle.hamiltonian_convention_id:
        raise ValueError(
            "hadronic contract hamiltonian_convention_id must match the frozen paper convention"
        )
    if hadronic_contract.matrix_element_formula_id != hadronic_bundle.matrix_element_formula_id:
        raise ValueError(
            "hadronic contract matrix_element_formula_id must match the B-system "
            "Q1 matrix-element formula"
        )
    if hadronic_contract.parity_relation_id != hadronic_bundle.parity_relation_id:
        raise ValueError(
            "hadronic contract parity_relation_id must match the B-system VLL/VRR parity relation"
        )
    if not math.isclose(
        hadronic_contract.mu_had_GeV,
        hadronic_bundle.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError("hadronic contract mu_had_GeV must match the B-system bundle mu_had_GeV")
    if not math.isclose(
        hadronic_contract.evaluation_scale_GeV,
        hadronic_bundle.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "hadronic contract evaluation_scale_GeV must match the B-system bundle mu_had_GeV"
        )
    if hadronic_bundle.operator_basis_id != PAPER_0710_1869_DELTAF2_OPERATOR_BASIS_ID:
        raise ValueError("hadronic bundle operator_basis_id must match the frozen paper basis")
    if (
        hadronic_bundle.operator_normalization_id
        != PAPER_0710_1869_DELTAF2_OPERATOR_NORMALIZATION_ID
    ):
        raise ValueError(
            "hadronic bundle operator_normalization_id must match the frozen paper normalization"
        )
    if (
        hadronic_bundle.hamiltonian_convention_id
        != PAPER_0710_1869_DELTAF2_HAMILTONIAN_CONVENTION_ID
    ):
        raise ValueError(
            "hadronic bundle hamiltonian_convention_id must match the frozen paper convention"
        )
    if (
        hadronic_bundle.matrix_element_formula_id
        != _CUSTOM_B_Q1_FORMULA_BY_SYSTEM[expected_system_id]
    ):
        raise ValueError(
            "hadronic bundle matrix_element_formula_id must match the frozen B-system Q1 formula"
        )
    if hadronic_bundle.parity_relation_id != _CUSTOM_B_Q1_PARITY_BY_SYSTEM[expected_system_id]:
        raise ValueError(
            "hadronic bundle parity_relation_id must match the frozen B-system "
            "VLL/VRR parity relation"
        )
    require_member(
        "hadronic_bundle.input_provenance_mode_id",
        hadronic_bundle.input_provenance_mode_id,
        (PAPER_0710_1869_DELTAF2_HADRONIC_CUSTOM_PROVENANCE_MODE_ID,),
    )
    if wilsons.system_id != expected_system_id:
        raise ValueError(
            f"custom {expected_system_id} observable evaluation requires "
            f"{expected_system_id} RG-tagged wilsons"
        )
    if wilsons.sector_id != "down":
        raise ValueError(
            f"custom {expected_system_id} observable evaluation requires the down-sector system"
        )
    if wilsons.operator_basis_id != hadronic_bundle.operator_basis_id:
        raise ValueError(
            "hadronic bundle operator_basis_id is incompatible with the RG Wilson snapshot"
        )
    if wilsons.operator_normalization_id != hadronic_bundle.operator_normalization_id:
        raise ValueError(
            "hadronic bundle operator_normalization_id is incompatible with the RG Wilson snapshot"
        )
    if wilsons.renormalization_scheme_id != hadronic_bundle.renormalization_scheme_id:
        raise ValueError(
            "hadronic bundle renormalization_scheme_id is incompatible with the RG Wilson snapshot"
        )
    if not math.isclose(
        wilsons.matching_scale_GeV,
        hadronic_bundle.mu_had_GeV,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "hadronic bundle mu_had_GeV must match the evolved Wilson evaluation scale"
        )
    if abs(wilsons.q4_lr) > _ZERO_TOLERANCE or abs(wilsons.q5_lr) > _ZERO_TOLERANCE:
        raise ValueError(
            f"custom {expected_system_id} Q1 observables currently support only Q1_VLL/Q1_VRR; "
            "non-zero Q4_LR/Q5_LR remain out of scope for the B-system slice."
        )


def _compute_custom_b_q1_m12_value(
    *,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot,
    hadronic_bundle: Paper07101869BMesonHadronicBundle,
    m12_name: str,
) -> complex:
    return _require_finite_complex(
        m12_name,
        (
            (wilsons.q1_vll + wilsons.q1_vrr)
            * hadronic_bundle.q1_matrix_element_GeV4
            / (2.0 * hadronic_bundle.meson_mass_GeV)
        ),
    )


def _compute_d0_q1_m12_value(
    *,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot,
    hadronic_bundle: Paper07101869D0HadronicBundle,
) -> complex:
    return _require_finite_complex(
        "M12_D0_NP_GeV",
        (
            (wilsons.q1_vll + wilsons.q1_vrr)
            * hadronic_bundle.q1_matrix_element_GeV4
            / (2.0 * hadronic_bundle.meson_mass_GeV)
        ),
    )


def _compute_lr_only_m12_value(
    *,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot,
    hadronic_inputs: Paper07101869KaonLRHadronicInputs,
) -> complex:
    return _require_finite_complex(
        "M12_K_LR_NP_GeV",
        (
            wilsons.q4_lr * hadronic_inputs.q4_matrix_element_GeV4
            + wilsons.q5_lr * hadronic_inputs.q5_matrix_element_GeV4
        )
        / (2.0 * hadronic_inputs.m_K0_GeV),
    )


@dataclass(frozen=True)
class Paper07101869KaonNPObservableResult:
    """Self-describing kaon NP-only observable result at the hadronic scale."""

    hadronic_bundle: Paper07101869KaonHadronicBundle
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot
    M12_K_NP_GeV: complex
    delta_m_K_NP_GeV: float
    observable_scope_id: str = PAPER_0710_1869_DELTAF2_KAON_NP_OBSERVABLE_SCOPE_ID
    interpretation: str = PAPER_0710_1869_DELTAF2_KAON_NP_INTERPRETATION_ID
    m12_observable_id: str = PAPER_0710_1869_DELTAF2_KAON_M12_OBSERVABLE_ID
    delta_m_observable_id: str = PAPER_0710_1869_DELTAF2_KAON_DELTA_M_OBSERVABLE_ID
    schema_id: str = PAPER_0710_1869_DELTAF2_KAON_NP_OBSERVABLE_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    notes: str = (
        "PR5a kaon NP-only observable result: M12_K^NP from Q1_VLL/Q1_VRR only at mu_had. "
        "epsilon_K and LR contributions remain outside the honest supported scope under "
        "the frozen LR audit contract."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_KAON_NP_OBSERVABLE_SCHEMA_ID,
            ),
        )
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        if not isinstance(self.hadronic_bundle, Paper07101869KaonHadronicBundle):
            raise ValueError("hadronic_bundle must be a Paper07101869KaonHadronicBundle")
        if not isinstance(self.wilsons, Paper07101869DeltaF2RGWilsonSnapshot):
            raise ValueError("wilsons must be a Paper07101869DeltaF2RGWilsonSnapshot")
        require_known_schema_id(
            "wilsons.schema_id",
            self.wilsons.schema_id,
            expected=PAPER_0710_1869_DELTAF2_RG_WILSON_SCHEMA_ID,
        )
        _require_hadronic_compatibility(wilsons=self.wilsons, hadronic_bundle=self.hadronic_bundle)
        for field_name in (
            "observable_scope_id",
            "interpretation",
            "m12_observable_id",
            "delta_m_observable_id",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        object.__setattr__(
            self,
            "M12_K_NP_GeV",
            _require_finite_complex("M12_K_NP_GeV", self.M12_K_NP_GeV),
        )
        numeric_delta_m = float(self.delta_m_K_NP_GeV)
        if not math.isfinite(numeric_delta_m):
            raise ValueError("delta_m_K_NP_GeV must be finite")
        object.__setattr__(self, "delta_m_K_NP_GeV", numeric_delta_m)
        expected_delta_m = 2.0 * float(self.M12_K_NP_GeV.real)
        if not math.isclose(
            self.delta_m_K_NP_GeV,
            expected_delta_m,
            rel_tol=0.0,
            abs_tol=1e-18,
        ):
            raise ValueError("delta_m_K_NP_GeV must equal 2 * Re(M12_K_NP_GeV)")

    @property
    def benchmark_id(self) -> str:
        return self.wilsons.benchmark_id

    @property
    def system_id(self) -> str:
        return self.hadronic_bundle.system_id

    @property
    def scale_label(self) -> str:
        return self.wilsons.scale_label

    @property
    def mu_had_GeV(self) -> float:
        return self.hadronic_bundle.mu_had_GeV

    @property
    def operator_basis_id(self) -> str:
        return self.hadronic_bundle.operator_basis_id

    @property
    def operator_normalization_id(self) -> str:
        return self.hadronic_bundle.operator_normalization_id

    @property
    def renormalization_scheme_id(self) -> str:
        return self.hadronic_bundle.renormalization_scheme_id

    @property
    def re_M12_K_NP_GeV(self) -> float:
        return float(self.M12_K_NP_GeV.real)

    @property
    def im_M12_K_NP_GeV(self) -> float:
        return float(self.M12_K_NP_GeV.imag)

    @property
    def tags(self) -> dict[str, object]:
        return {
            **self.wilsons.tags,
            **self.hadronic_bundle.tags,
            "observable_scope_id": self.observable_scope_id,
            "interpretation": self.interpretation,
            "m12_observable_id": self.m12_observable_id,
            "delta_m_observable_id": self.delta_m_observable_id,
        }

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "benchmark_id": self.benchmark_id,
            "system_id": self.system_id,
            "scale_label": self.scale_label,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "rg_scheme_id": self.wilsons.contract.rg_scheme_id,
            "matrix_element_formula_id": self.hadronic_bundle.matrix_element_formula_id,
            "hamiltonian_convention_id": self.hadronic_bundle.hamiltonian_convention_id,
            "observable_scope_id": self.observable_scope_id,
            "interpretation": self.interpretation,
            "m12_observable_id": self.m12_observable_id,
            "delta_m_observable_id": self.delta_m_observable_id,
            "mu_had_GeV": self.mu_had_GeV,
            "matching_scale_GeV": self.wilsons.matching_scale_GeV,
            "source_matching_scale_GeV": self.wilsons.source_matching_scale_GeV,
            "propagator_mass_GeV": self.wilsons.propagator_mass_GeV,
            "q1_matrix_element_GeV4": self.hadronic_bundle.q1_matrix_element_GeV4,
            "M12_K_NP_GeV": _complex_as_dict(self.M12_K_NP_GeV),
            "re_M12_K_NP_GeV": self.re_M12_K_NP_GeV,
            "im_M12_K_NP_GeV": self.im_M12_K_NP_GeV,
            "delta_m_K_NP_GeV": self.delta_m_K_NP_GeV,
            "observables": {
                self.m12_observable_id: {
                    "re": self.re_M12_K_NP_GeV,
                    "im": self.im_M12_K_NP_GeV,
                },
                self.delta_m_observable_id: self.delta_m_K_NP_GeV,
            },
            "coefficients": {
                "Q1_VLL": _complex_as_dict(self.wilsons.q1_vll),
                "Q1_VRR": _complex_as_dict(self.wilsons.q1_vrr),
                PAPER_0710_1869_DELTAF2_Q4_LR: _complex_as_dict(self.wilsons.q4_lr),
                PAPER_0710_1869_DELTAF2_Q5_LR: _complex_as_dict(self.wilsons.q5_lr),
            },
            "hadronic_bundle": self.hadronic_bundle.as_dict(),
            "wilsons": self.wilsons.as_dict(),
            "tags": self.tags,
            "notes": self.notes,
        }

    def summary(self) -> dict[str, object]:
        return {
            **self.as_dict(),
            "schema_id": PAPER_0710_1869_DELTAF2_KAON_NP_OBSERVABLE_SUMMARY_SCHEMA_ID,
        }


@dataclass(frozen=True)
class Paper07101869BMesonCustomQ1ObservableResult:
    """Self-describing custom Q1-only B_d/B_s observable result at ``mu_had``."""

    hadronic_bundle: Paper07101869BMesonHadronicBundle
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot
    M12_NP_GeV: complex
    delta_m_NP_GeV: float
    observable_scope_id: str | None = None
    interpretation: str | None = None
    m12_observable_id: str | None = None
    delta_m_observable_id: str | None = None
    schema_id: str | None = None
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    notes: str | None = None

    def __post_init__(self) -> None:
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        if not isinstance(self.hadronic_bundle, Paper07101869BMesonHadronicBundle):
            raise ValueError("hadronic_bundle must be a Paper07101869BMesonHadronicBundle")
        if not isinstance(self.wilsons, Paper07101869DeltaF2RGWilsonSnapshot):
            raise ValueError("wilsons must be a Paper07101869DeltaF2RGWilsonSnapshot")
        require_known_schema_id(
            "wilsons.schema_id",
            self.wilsons.schema_id,
            expected=PAPER_0710_1869_DELTAF2_RG_WILSON_SCHEMA_ID,
        )
        system_id = self.hadronic_bundle.system_id
        expected_schema_id = _CUSTOM_B_Q1_SCHEMA_BY_SYSTEM[system_id]
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                expected_schema_id if self.schema_id is None else self.schema_id,
                expected=expected_schema_id,
            ),
        )
        _require_custom_b_q1_hadronic_compatibility(
            wilsons=self.wilsons,
            hadronic_bundle=self.hadronic_bundle,
            expected_system_id=system_id,
        )
        resolved_scope_id = (
            _CUSTOM_B_Q1_SCOPE_BY_SYSTEM[system_id]
            if self.observable_scope_id is None
            else self.observable_scope_id
        )
        resolved_interpretation_id = (
            _CUSTOM_B_Q1_INTERPRETATION_BY_SYSTEM[system_id]
            if self.interpretation is None
            else self.interpretation
        )
        resolved_m12_id = (
            _CUSTOM_B_Q1_M12_ID_BY_SYSTEM[system_id]
            if self.m12_observable_id is None
            else self.m12_observable_id
        )
        resolved_delta_m_id = (
            _CUSTOM_B_Q1_DELTA_M_ID_BY_SYSTEM[system_id]
            if self.delta_m_observable_id is None
            else self.delta_m_observable_id
        )
        for field_name, value in (
            ("observable_scope_id", resolved_scope_id),
            ("interpretation", resolved_interpretation_id),
            ("m12_observable_id", resolved_m12_id),
            ("delta_m_observable_id", resolved_delta_m_id),
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, value),
            )
        resolved_notes = (
            f"Custom {system_id} NP-only observable result: M12 and Delta_m from Q1_VLL/Q1_VRR "
            "only under exact RG/hadronic alignment. This does not widen the default/exported "
            "kaon Q1-only observable or artifact surface."
            if self.notes is None
            else self.notes
        )
        object.__setattr__(self, "notes", require_nonempty_identifier("notes", resolved_notes))
        object.__setattr__(
            self,
            "M12_NP_GeV",
            _require_finite_complex("M12_NP_GeV", self.M12_NP_GeV),
        )
        numeric_delta_m = float(self.delta_m_NP_GeV)
        if not math.isfinite(numeric_delta_m):
            raise ValueError("delta_m_NP_GeV must be finite")
        object.__setattr__(self, "delta_m_NP_GeV", numeric_delta_m)
        expected_delta_m = 2.0 * float(self.M12_NP_GeV.real)
        if not math.isclose(self.delta_m_NP_GeV, expected_delta_m, rel_tol=0.0, abs_tol=1e-18):
            raise ValueError("delta_m_NP_GeV must equal 2 * Re(M12_NP_GeV)")

    @property
    def benchmark_id(self) -> str:
        return self.wilsons.benchmark_id

    @property
    def system_id(self) -> str:
        return self.hadronic_bundle.system_id

    @property
    def scale_label(self) -> str:
        return self.wilsons.scale_label

    @property
    def mu_had_GeV(self) -> float:
        return self.hadronic_bundle.mu_had_GeV

    @property
    def operator_basis_id(self) -> str:
        return self.hadronic_bundle.operator_basis_id

    @property
    def operator_normalization_id(self) -> str:
        return self.hadronic_bundle.operator_normalization_id

    @property
    def renormalization_scheme_id(self) -> str:
        return self.hadronic_bundle.renormalization_scheme_id

    @property
    def re_M12_NP_GeV(self) -> float:
        return float(self.M12_NP_GeV.real)

    @property
    def im_M12_NP_GeV(self) -> float:
        return float(self.M12_NP_GeV.imag)

    @property
    def tags(self) -> dict[str, object]:
        return {
            **self.wilsons.tags,
            **self.hadronic_bundle.tags,
            "observable_scope_id": self.observable_scope_id,
            "interpretation": self.interpretation,
            "m12_observable_id": self.m12_observable_id,
            "delta_m_observable_id": self.delta_m_observable_id,
        }

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "benchmark_id": self.benchmark_id,
            "system_id": self.system_id,
            "scale_label": self.scale_label,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "rg_scheme_id": self.wilsons.contract.rg_scheme_id,
            "matrix_element_formula_id": self.hadronic_bundle.matrix_element_formula_id,
            "hamiltonian_convention_id": self.hadronic_bundle.hamiltonian_convention_id,
            "observable_scope_id": self.observable_scope_id,
            "interpretation": self.interpretation,
            "m12_observable_id": self.m12_observable_id,
            "delta_m_observable_id": self.delta_m_observable_id,
            "mu_had_GeV": self.mu_had_GeV,
            "matching_scale_GeV": self.wilsons.matching_scale_GeV,
            "source_matching_scale_GeV": self.wilsons.source_matching_scale_GeV,
            "propagator_mass_GeV": self.wilsons.propagator_mass_GeV,
            "q1_matrix_element_GeV4": self.hadronic_bundle.q1_matrix_element_GeV4,
            "M12_NP_GeV": _complex_as_dict(self.M12_NP_GeV),
            "re_M12_NP_GeV": self.re_M12_NP_GeV,
            "im_M12_NP_GeV": self.im_M12_NP_GeV,
            "delta_m_NP_GeV": self.delta_m_NP_GeV,
            "observables": {
                self.m12_observable_id: {
                    "re": self.re_M12_NP_GeV,
                    "im": self.im_M12_NP_GeV,
                },
                self.delta_m_observable_id: self.delta_m_NP_GeV,
            },
            "coefficients": {
                "Q1_VLL": _complex_as_dict(self.wilsons.q1_vll),
                "Q1_VRR": _complex_as_dict(self.wilsons.q1_vrr),
                PAPER_0710_1869_DELTAF2_Q4_LR: _complex_as_dict(self.wilsons.q4_lr),
                PAPER_0710_1869_DELTAF2_Q5_LR: _complex_as_dict(self.wilsons.q5_lr),
            },
            "hadronic_bundle": self.hadronic_bundle.as_dict(),
            "wilsons": self.wilsons.as_dict(),
            "tags": self.tags,
            "notes": self.notes,
        }

    def summary(self) -> dict[str, object]:
        return {
            **self.as_dict(),
            "schema_id": _CUSTOM_B_Q1_SUMMARY_SCHEMA_BY_SYSTEM[self.system_id],
        }


@dataclass(frozen=True)
class Paper07101869D0CustomQ1ObservableResult:
    """Self-describing custom Q1-only D0 observable result at ``mu_had``."""

    hadronic_bundle: Paper07101869D0HadronicBundle
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot
    M12_D0_NP_GeV: complex
    delta_m_D0_NP_GeV: float
    observable_scope_id: str = PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_OBSERVABLE_SCOPE_ID
    interpretation: str = PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_INTERPRETATION_ID
    m12_observable_id: str = PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_M12_OBSERVABLE_ID
    delta_m_observable_id: str = PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_DELTA_M_OBSERVABLE_ID
    schema_id: str = PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_OBSERVABLE_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    notes: str = (
        "Custom D0 NP-only observable result: M12_D0^NP and Delta_m_D0^NP from "
        "Q1_VLL/Q1_VRR only under exact RG/hadronic alignment. This does not widen the "
        "default/exported kaon Q1-only observable or artifact surface."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_OBSERVABLE_SCHEMA_ID,
            ),
        )
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        if not isinstance(self.hadronic_bundle, Paper07101869D0HadronicBundle):
            raise ValueError("hadronic_bundle must be a Paper07101869D0HadronicBundle")
        if not isinstance(self.wilsons, Paper07101869DeltaF2RGWilsonSnapshot):
            raise ValueError("wilsons must be a Paper07101869DeltaF2RGWilsonSnapshot")
        require_known_schema_id(
            "wilsons.schema_id",
            self.wilsons.schema_id,
            expected=PAPER_0710_1869_DELTAF2_RG_WILSON_SCHEMA_ID,
        )
        _require_d0_q1_hadronic_compatibility(
            wilsons=self.wilsons,
            hadronic_bundle=self.hadronic_bundle,
        )
        for field_name in (
            "observable_scope_id",
            "interpretation",
            "m12_observable_id",
            "delta_m_observable_id",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        object.__setattr__(
            self,
            "M12_D0_NP_GeV",
            _require_finite_complex("M12_D0_NP_GeV", self.M12_D0_NP_GeV),
        )
        numeric_delta_m = float(self.delta_m_D0_NP_GeV)
        if not math.isfinite(numeric_delta_m):
            raise ValueError("delta_m_D0_NP_GeV must be finite")
        object.__setattr__(self, "delta_m_D0_NP_GeV", numeric_delta_m)
        expected_delta_m = 2.0 * float(self.M12_D0_NP_GeV.real)
        if not math.isclose(
            self.delta_m_D0_NP_GeV,
            expected_delta_m,
            rel_tol=0.0,
            abs_tol=1e-18,
        ):
            raise ValueError("delta_m_D0_NP_GeV must equal 2 * Re(M12_D0_NP_GeV)")

    @property
    def benchmark_id(self) -> str:
        return self.wilsons.benchmark_id

    @property
    def system_id(self) -> str:
        return self.hadronic_bundle.system_id

    @property
    def scale_label(self) -> str:
        return self.wilsons.scale_label

    @property
    def mu_had_GeV(self) -> float:
        return self.hadronic_bundle.mu_had_GeV

    @property
    def operator_basis_id(self) -> str:
        return self.hadronic_bundle.operator_basis_id

    @property
    def operator_normalization_id(self) -> str:
        return self.hadronic_bundle.operator_normalization_id

    @property
    def renormalization_scheme_id(self) -> str:
        return self.hadronic_bundle.renormalization_scheme_id

    @property
    def re_M12_D0_NP_GeV(self) -> float:
        return float(self.M12_D0_NP_GeV.real)

    @property
    def im_M12_D0_NP_GeV(self) -> float:
        return float(self.M12_D0_NP_GeV.imag)

    @property
    def tags(self) -> dict[str, object]:
        return {
            **self.wilsons.tags,
            **self.hadronic_bundle.tags,
            "observable_scope_id": self.observable_scope_id,
            "interpretation": self.interpretation,
            "m12_observable_id": self.m12_observable_id,
            "delta_m_observable_id": self.delta_m_observable_id,
        }

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "benchmark_id": self.benchmark_id,
            "system_id": self.system_id,
            "scale_label": self.scale_label,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "rg_scheme_id": self.wilsons.contract.rg_scheme_id,
            "matrix_element_formula_id": self.hadronic_bundle.matrix_element_formula_id,
            "hamiltonian_convention_id": self.hadronic_bundle.hamiltonian_convention_id,
            "observable_scope_id": self.observable_scope_id,
            "interpretation": self.interpretation,
            "m12_observable_id": self.m12_observable_id,
            "delta_m_observable_id": self.delta_m_observable_id,
            "mu_had_GeV": self.mu_had_GeV,
            "matching_scale_GeV": self.wilsons.matching_scale_GeV,
            "source_matching_scale_GeV": self.wilsons.source_matching_scale_GeV,
            "propagator_mass_GeV": self.wilsons.propagator_mass_GeV,
            "q1_matrix_element_GeV4": self.hadronic_bundle.q1_matrix_element_GeV4,
            "M12_NP_GeV": _complex_as_dict(self.M12_D0_NP_GeV),
            "M12_D0_NP_GeV": _complex_as_dict(self.M12_D0_NP_GeV),
            "re_M12_NP_GeV": self.re_M12_D0_NP_GeV,
            "re_M12_D0_NP_GeV": self.re_M12_D0_NP_GeV,
            "im_M12_NP_GeV": self.im_M12_D0_NP_GeV,
            "im_M12_D0_NP_GeV": self.im_M12_D0_NP_GeV,
            "delta_m_NP_GeV": self.delta_m_D0_NP_GeV,
            "delta_m_D0_NP_GeV": self.delta_m_D0_NP_GeV,
            "observables": {
                self.m12_observable_id: {
                    "re": self.re_M12_D0_NP_GeV,
                    "im": self.im_M12_D0_NP_GeV,
                },
                self.delta_m_observable_id: self.delta_m_D0_NP_GeV,
            },
            "coefficients": {
                "Q1_VLL": _complex_as_dict(self.wilsons.q1_vll),
                "Q1_VRR": _complex_as_dict(self.wilsons.q1_vrr),
                PAPER_0710_1869_DELTAF2_Q4_LR: _complex_as_dict(self.wilsons.q4_lr),
                PAPER_0710_1869_DELTAF2_Q5_LR: _complex_as_dict(self.wilsons.q5_lr),
            },
            "hadronic_bundle": self.hadronic_bundle.as_dict(),
            "wilsons": self.wilsons.as_dict(),
            "tags": self.tags,
            "notes": self.notes,
        }

    def summary(self) -> dict[str, object]:
        return {
            **self.as_dict(),
            "schema_id": PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_OBSERVABLE_SUMMARY_SCHEMA_ID,
        }


@dataclass(frozen=True)
class Paper07101869KaonLROnlyObservableResult:
    """Self-describing custom LR-only kaon observable result at ``mu_had``."""

    hadronic_inputs: Paper07101869KaonLRHadronicInputs
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot
    M12_K_LR_NP_GeV: complex
    delta_m_K_LR_NP_GeV: float
    observable_scope_id: str = PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_OBSERVABLE_SCOPE_ID
    interpretation: str = PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_INTERPRETATION_ID
    m12_observable_id: str = PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_M12_OBSERVABLE_ID
    delta_m_observable_id: str = (
        PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_DELTA_M_OBSERVABLE_ID
    )
    schema_id: str = PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_OBSERVABLE_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    notes: str = (
        "LR-OBS-1 custom LR-only kaon observable result: M12_K^NP from Q4_LR/Q5_LR only "
        "using explicit custom LR hadronic inputs at mu_had under exact RG/hadronic "
        "scheme and scale alignment. This does not widen the default Q1-only observable "
        "or artifact surface."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_OBSERVABLE_SCHEMA_ID,
            ),
        )
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        if not isinstance(self.hadronic_inputs, Paper07101869KaonLRHadronicInputs):
            raise ValueError("hadronic_inputs must be a Paper07101869KaonLRHadronicInputs")
        if not isinstance(self.wilsons, Paper07101869DeltaF2RGWilsonSnapshot):
            raise ValueError("wilsons must be a Paper07101869DeltaF2RGWilsonSnapshot")
        require_known_schema_id(
            "wilsons.schema_id",
            self.wilsons.schema_id,
            expected=PAPER_0710_1869_DELTAF2_RG_WILSON_SCHEMA_ID,
        )
        _require_lr_hadronic_compatibility(
            wilsons=self.wilsons,
            hadronic_inputs=self.hadronic_inputs,
        )
        for field_name in (
            "observable_scope_id",
            "interpretation",
            "m12_observable_id",
            "delta_m_observable_id",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        object.__setattr__(
            self,
            "M12_K_LR_NP_GeV",
            _require_finite_complex("M12_K_LR_NP_GeV", self.M12_K_LR_NP_GeV),
        )
        numeric_delta_m = float(self.delta_m_K_LR_NP_GeV)
        if not math.isfinite(numeric_delta_m):
            raise ValueError("delta_m_K_LR_NP_GeV must be finite")
        object.__setattr__(self, "delta_m_K_LR_NP_GeV", numeric_delta_m)
        expected_delta_m = 2.0 * float(self.M12_K_LR_NP_GeV.real)
        if not math.isclose(
            self.delta_m_K_LR_NP_GeV,
            expected_delta_m,
            rel_tol=0.0,
            abs_tol=1e-18,
        ):
            raise ValueError("delta_m_K_LR_NP_GeV must equal 2 * Re(M12_K_LR_NP_GeV)")

    @property
    def benchmark_id(self) -> str:
        return self.wilsons.benchmark_id

    @property
    def system_id(self) -> str:
        return self.hadronic_inputs.system_id

    @property
    def scale_label(self) -> str:
        return self.wilsons.scale_label

    @property
    def mu_had_GeV(self) -> float:
        return self.hadronic_inputs.mu_had_GeV

    @property
    def operator_basis_id(self) -> str:
        return self.hadronic_inputs.operator_basis_id

    @property
    def operator_normalization_id(self) -> str:
        return self.hadronic_inputs.operator_normalization_id

    @property
    def renormalization_scheme_id(self) -> str:
        return self.hadronic_inputs.renormalization_scheme_id

    @property
    def re_M12_K_LR_NP_GeV(self) -> float:
        return float(self.M12_K_LR_NP_GeV.real)

    @property
    def im_M12_K_LR_NP_GeV(self) -> float:
        return float(self.M12_K_LR_NP_GeV.imag)

    @property
    def tags(self) -> dict[str, object]:
        return {
            **self.wilsons.tags,
            **self.hadronic_inputs.tags,
            "observable_scope_id": self.observable_scope_id,
            "interpretation": self.interpretation,
            "m12_observable_id": self.m12_observable_id,
            "delta_m_observable_id": self.delta_m_observable_id,
        }

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "benchmark_id": self.benchmark_id,
            "system_id": self.system_id,
            "scale_label": self.scale_label,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "rg_scheme_id": self.wilsons.contract.rg_scheme_id,
            "lr_basis_contract_id": self.wilsons.contract.lr_basis_contract_id,
            "lr_basis_status_id": self.wilsons.contract.lr_basis_status_id,
            "lr_running_bridge_id": self.wilsons.contract.lr_running_bridge_id,
            "hamiltonian_convention_id": self.hadronic_inputs.hamiltonian_convention_id,
            "q4_matrix_element_formula_id": self.hadronic_inputs.q4_matrix_element_formula_id,
            "q5_matrix_element_formula_id": self.hadronic_inputs.q5_matrix_element_formula_id,
            "observable_scope_id": self.observable_scope_id,
            "interpretation": self.interpretation,
            "m12_observable_id": self.m12_observable_id,
            "delta_m_observable_id": self.delta_m_observable_id,
            "mu_had_GeV": self.mu_had_GeV,
            "matching_scale_GeV": self.wilsons.matching_scale_GeV,
            "source_matching_scale_GeV": self.wilsons.source_matching_scale_GeV,
            "propagator_mass_GeV": self.wilsons.propagator_mass_GeV,
            "q4_matrix_element_GeV4": self.hadronic_inputs.q4_matrix_element_GeV4,
            "q5_matrix_element_GeV4": self.hadronic_inputs.q5_matrix_element_GeV4,
            "M12_K_LR_NP_GeV": _complex_as_dict(self.M12_K_LR_NP_GeV),
            "re_M12_K_LR_NP_GeV": self.re_M12_K_LR_NP_GeV,
            "im_M12_K_LR_NP_GeV": self.im_M12_K_LR_NP_GeV,
            "delta_m_K_LR_NP_GeV": self.delta_m_K_LR_NP_GeV,
            "observables": {
                self.m12_observable_id: {
                    "re": self.re_M12_K_LR_NP_GeV,
                    "im": self.im_M12_K_LR_NP_GeV,
                },
                self.delta_m_observable_id: self.delta_m_K_LR_NP_GeV,
            },
            "coefficients": {
                PAPER_0710_1869_DELTAF2_Q4_LR: _complex_as_dict(self.wilsons.q4_lr),
                PAPER_0710_1869_DELTAF2_Q5_LR: _complex_as_dict(self.wilsons.q5_lr),
            },
            "hadronic_inputs": self.hadronic_inputs.as_dict(),
            "wilsons": self.wilsons.as_dict(),
            "tags": self.tags,
            "notes": self.notes,
        }

    def summary(self) -> dict[str, object]:
        return {
            **self.as_dict(),
            "schema_id": PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_OBSERVABLE_SUMMARY_SCHEMA_ID,
        }


@dataclass(frozen=True)
class Paper07101869KaonCustomTotalObservableResult:
    """Self-describing custom combined kaon observable result at ``mu_had``."""

    q1_hadronic_bundle: Paper07101869KaonHadronicBundle
    lr_hadronic_inputs: Paper07101869KaonLRHadronicInputs
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot
    M12_K_NP_Q1_GeV: complex
    M12_K_LR_NP_GeV: complex
    M12_K_NP_TOTAL_GeV: complex
    delta_m_K_NP_TOTAL_GeV: float
    observable_scope_id: str = PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_OBSERVABLE_SCOPE_ID
    interpretation: str = PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_INTERPRETATION_ID
    m12_observable_id: str = PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_M12_OBSERVABLE_ID
    delta_m_observable_id: str = (
        PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_DELTA_M_OBSERVABLE_ID
    )
    schema_id: str = PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_OBSERVABLE_SCHEMA_ID
    mode_id: str = PAPER_0710_1869_MODE_ID
    paper_id: str = PAPER_0710_1869_PAPER_ID
    notes: str = (
        "LR-TOTAL-1 custom combined kaon observable result: M12_K^NP from the "
        "Q1_VLL/Q1_VRR and Q4_LR/Q5_LR contributions summed under exact RG/Q1/LR "
        "hadronic alignment. This leaves the default/exported observable and artifact "
        "surface unchanged and keeps epsilon_K blocked."
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "schema_id",
            require_known_schema_id(
                "schema_id",
                self.schema_id,
                expected=PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_OBSERVABLE_SCHEMA_ID,
            ),
        )
        require_member("mode_id", self.mode_id, (PAPER_0710_1869_MODE_ID,))
        require_member("paper_id", self.paper_id, (PAPER_0710_1869_PAPER_ID,))
        if not isinstance(self.q1_hadronic_bundle, Paper07101869KaonHadronicBundle):
            raise ValueError("q1_hadronic_bundle must be a Paper07101869KaonHadronicBundle")
        if not isinstance(self.lr_hadronic_inputs, Paper07101869KaonLRHadronicInputs):
            raise ValueError("lr_hadronic_inputs must be a Paper07101869KaonLRHadronicInputs")
        if not isinstance(self.wilsons, Paper07101869DeltaF2RGWilsonSnapshot):
            raise ValueError("wilsons must be a Paper07101869DeltaF2RGWilsonSnapshot")
        require_known_schema_id(
            "wilsons.schema_id",
            self.wilsons.schema_id,
            expected=PAPER_0710_1869_DELTAF2_RG_WILSON_SCHEMA_ID,
        )
        _require_custom_q1_hadronic_compatibility(
            wilsons=self.wilsons,
            hadronic_bundle=self.q1_hadronic_bundle,
        )
        _require_lr_hadronic_compatibility(
            wilsons=self.wilsons,
            hadronic_inputs=self.lr_hadronic_inputs,
        )
        _require_custom_combined_alignment(
            wilsons=self.wilsons,
            q1_hadronic_bundle=self.q1_hadronic_bundle,
            lr_hadronic_inputs=self.lr_hadronic_inputs,
        )
        for field_name in (
            "observable_scope_id",
            "interpretation",
            "m12_observable_id",
            "delta_m_observable_id",
            "notes",
        ):
            object.__setattr__(
                self,
                field_name,
                require_nonempty_identifier(field_name, getattr(self, field_name)),
            )
        for field_name in (
            "M12_K_NP_Q1_GeV",
            "M12_K_LR_NP_GeV",
            "M12_K_NP_TOTAL_GeV",
        ):
            object.__setattr__(
                self,
                field_name,
                _require_finite_complex(field_name, getattr(self, field_name)),
            )
        expected_total = self.M12_K_NP_Q1_GeV + self.M12_K_LR_NP_GeV
        if abs(self.M12_K_NP_TOTAL_GeV - expected_total) > 1e-24:
            raise ValueError("M12_K_NP_TOTAL_GeV must equal M12_K_NP_Q1_GeV + M12_K_LR_NP_GeV")
        numeric_delta_m = float(self.delta_m_K_NP_TOTAL_GeV)
        if not math.isfinite(numeric_delta_m):
            raise ValueError("delta_m_K_NP_TOTAL_GeV must be finite")
        object.__setattr__(self, "delta_m_K_NP_TOTAL_GeV", numeric_delta_m)
        expected_delta_m = 2.0 * float(self.M12_K_NP_TOTAL_GeV.real)
        if not math.isclose(
            self.delta_m_K_NP_TOTAL_GeV,
            expected_delta_m,
            rel_tol=0.0,
            abs_tol=1e-18,
        ):
            raise ValueError("delta_m_K_NP_TOTAL_GeV must equal 2 * Re(M12_K_NP_TOTAL_GeV)")

    @property
    def benchmark_id(self) -> str:
        return self.wilsons.benchmark_id

    @property
    def system_id(self) -> str:
        return self.q1_hadronic_bundle.system_id

    @property
    def scale_label(self) -> str:
        return self.wilsons.scale_label

    @property
    def mu_had_GeV(self) -> float:
        return self.q1_hadronic_bundle.mu_had_GeV

    @property
    def operator_basis_id(self) -> str:
        return self.q1_hadronic_bundle.operator_basis_id

    @property
    def operator_normalization_id(self) -> str:
        return self.q1_hadronic_bundle.operator_normalization_id

    @property
    def renormalization_scheme_id(self) -> str:
        return self.q1_hadronic_bundle.renormalization_scheme_id

    @property
    def re_M12_K_NP_TOTAL_GeV(self) -> float:
        return float(self.M12_K_NP_TOTAL_GeV.real)

    @property
    def im_M12_K_NP_TOTAL_GeV(self) -> float:
        return float(self.M12_K_NP_TOTAL_GeV.imag)

    @property
    def tags(self) -> dict[str, object]:
        return {
            **self.wilsons.tags,
            **self.q1_hadronic_bundle.tags,
            **self.lr_hadronic_inputs.tags,
            "observable_scope_id": self.observable_scope_id,
            "interpretation": self.interpretation,
            "m12_observable_id": self.m12_observable_id,
            "delta_m_observable_id": self.delta_m_observable_id,
        }

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_id": self.schema_id,
            "mode_id": self.mode_id,
            "paper_id": self.paper_id,
            "benchmark_id": self.benchmark_id,
            "system_id": self.system_id,
            "scale_label": self.scale_label,
            "operator_basis_id": self.operator_basis_id,
            "operator_normalization_id": self.operator_normalization_id,
            "renormalization_scheme_id": self.renormalization_scheme_id,
            "rg_scheme_id": self.wilsons.contract.rg_scheme_id,
            "lr_basis_contract_id": self.wilsons.contract.lr_basis_contract_id,
            "lr_basis_status_id": self.wilsons.contract.lr_basis_status_id,
            "lr_running_bridge_id": self.wilsons.contract.lr_running_bridge_id,
            "hamiltonian_convention_id": self.q1_hadronic_bundle.hamiltonian_convention_id,
            "observable_scope_id": self.observable_scope_id,
            "interpretation": self.interpretation,
            "m12_observable_id": self.m12_observable_id,
            "delta_m_observable_id": self.delta_m_observable_id,
            "mu_had_GeV": self.mu_had_GeV,
            "matching_scale_GeV": self.wilsons.matching_scale_GeV,
            "source_matching_scale_GeV": self.wilsons.source_matching_scale_GeV,
            "propagator_mass_GeV": self.wilsons.propagator_mass_GeV,
            "q1_matrix_element_GeV4": self.q1_hadronic_bundle.q1_matrix_element_GeV4,
            "q4_matrix_element_GeV4": self.lr_hadronic_inputs.q4_matrix_element_GeV4,
            "q5_matrix_element_GeV4": self.lr_hadronic_inputs.q5_matrix_element_GeV4,
            "M12_K_NP_Q1_GeV": _complex_as_dict(self.M12_K_NP_Q1_GeV),
            "M12_K_LR_NP_GeV": _complex_as_dict(self.M12_K_LR_NP_GeV),
            "M12_K_NP_TOTAL_GeV": _complex_as_dict(self.M12_K_NP_TOTAL_GeV),
            "re_M12_K_NP_TOTAL_GeV": self.re_M12_K_NP_TOTAL_GeV,
            "im_M12_K_NP_TOTAL_GeV": self.im_M12_K_NP_TOTAL_GeV,
            "delta_m_K_NP_TOTAL_GeV": self.delta_m_K_NP_TOTAL_GeV,
            "observables": {
                self.m12_observable_id: {
                    "re": self.re_M12_K_NP_TOTAL_GeV,
                    "im": self.im_M12_K_NP_TOTAL_GeV,
                },
                self.delta_m_observable_id: self.delta_m_K_NP_TOTAL_GeV,
            },
            "component_observables": {
                PAPER_0710_1869_DELTAF2_KAON_M12_OBSERVABLE_ID: _complex_as_dict(
                    self.M12_K_NP_Q1_GeV
                ),
                PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_M12_OBSERVABLE_ID: _complex_as_dict(
                    self.M12_K_LR_NP_GeV
                ),
            },
            "coefficients": {
                "Q1_VLL": _complex_as_dict(self.wilsons.q1_vll),
                "Q1_VRR": _complex_as_dict(self.wilsons.q1_vrr),
                PAPER_0710_1869_DELTAF2_Q4_LR: _complex_as_dict(self.wilsons.q4_lr),
                PAPER_0710_1869_DELTAF2_Q5_LR: _complex_as_dict(self.wilsons.q5_lr),
            },
            "q1_hadronic_bundle": self.q1_hadronic_bundle.as_dict(),
            "lr_hadronic_inputs": self.lr_hadronic_inputs.as_dict(),
            "wilsons": self.wilsons.as_dict(),
            "tags": self.tags,
            "notes": self.notes,
        }

    def summary(self) -> dict[str, object]:
        return {
            **self.as_dict(),
            "schema_id": PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_OBSERVABLE_SUMMARY_SCHEMA_ID,
        }


def compute_kaon_np_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
) -> Paper07101869KaonNPObservableResult:
    """Compute the kaon NP-only observable subset from evolved paper-mode Wilsons."""

    resolved_hadronic_bundle = (
        default_paper_0710_1869_kaon_hadronic_bundle()
        if hadronic_bundle is None
        else hadronic_bundle
    )
    if not isinstance(resolved_hadronic_bundle, Paper07101869KaonHadronicBundle):
        raise ValueError("hadronic_bundle must be a Paper07101869KaonHadronicBundle")
    resolved_input = (
        run_deltaf2_wilsons_lo(default_paper_0710_1869_kaon_matching().wilsons)
        if rg_result_or_wilsons is None
        else rg_result_or_wilsons
    )
    resolved_wilsons = _resolve_rg_wilsons(resolved_input)
    _require_hadronic_compatibility(
        wilsons=resolved_wilsons,
        hadronic_bundle=resolved_hadronic_bundle,
    )
    m12_value = _compute_q1_only_m12_value(
        wilsons=resolved_wilsons,
        hadronic_bundle=resolved_hadronic_bundle,
    )
    return Paper07101869KaonNPObservableResult(
        hadronic_bundle=resolved_hadronic_bundle,
        wilsons=resolved_wilsons,
        M12_K_NP_GeV=m12_value,
        delta_m_K_NP_GeV=2.0 * float(m12_value.real),
    )


def compute_kaon_lr_only_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_inputs: Paper07101869KaonLRHadronicInputs | None = None,
) -> Paper07101869KaonLROnlyObservableResult:
    """Compute the custom LR-only kaon observable subset from evolved paper-mode Wilsons."""

    if hadronic_inputs is None:
        raise ValueError(
            "hadronic_inputs must be provided for LR-only observable evaluation"
        )
    if rg_result_or_wilsons is None:
        raise ValueError(
            "rg_result_or_wilsons must be provided for LR-only observable evaluation"
        )
    if not isinstance(hadronic_inputs, Paper07101869KaonLRHadronicInputs):
        raise ValueError("hadronic_inputs must be a Paper07101869KaonLRHadronicInputs")
    resolved_wilsons = _resolve_rg_wilsons(rg_result_or_wilsons)
    _require_lr_hadronic_compatibility(
        wilsons=resolved_wilsons,
        hadronic_inputs=hadronic_inputs,
    )
    m12_value = _compute_lr_only_m12_value(
        wilsons=resolved_wilsons,
        hadronic_inputs=hadronic_inputs,
    )
    return Paper07101869KaonLROnlyObservableResult(
        hadronic_inputs=hadronic_inputs,
        wilsons=resolved_wilsons,
        M12_K_LR_NP_GeV=m12_value,
        delta_m_K_LR_NP_GeV=2.0 * float(m12_value.real),
    )


def compute_kaon_custom_total_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    q1_hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
    lr_hadronic_inputs: Paper07101869KaonLRHadronicInputs | None = None,
) -> Paper07101869KaonCustomTotalObservableResult:
    """Compute the custom combined kaon observable subset from evolved paper-mode Wilsons."""

    if q1_hadronic_bundle is None:
        raise ValueError(
            "q1_hadronic_bundle must be provided for custom combined observable evaluation"
        )
    if lr_hadronic_inputs is None:
        raise ValueError(
            "lr_hadronic_inputs must be provided for custom combined observable evaluation"
        )
    if rg_result_or_wilsons is None:
        raise ValueError(
            "rg_result_or_wilsons must be provided for custom combined observable evaluation"
        )
    if not isinstance(q1_hadronic_bundle, Paper07101869KaonHadronicBundle):
        raise ValueError("q1_hadronic_bundle must be a Paper07101869KaonHadronicBundle")
    if not isinstance(lr_hadronic_inputs, Paper07101869KaonLRHadronicInputs):
        raise ValueError("lr_hadronic_inputs must be a Paper07101869KaonLRHadronicInputs")
    resolved_wilsons = _resolve_rg_wilsons(rg_result_or_wilsons)
    _require_custom_q1_hadronic_compatibility(
        wilsons=resolved_wilsons,
        hadronic_bundle=q1_hadronic_bundle,
    )
    _require_lr_hadronic_compatibility(
        wilsons=resolved_wilsons,
        hadronic_inputs=lr_hadronic_inputs,
    )
    _require_custom_combined_alignment(
        wilsons=resolved_wilsons,
        q1_hadronic_bundle=q1_hadronic_bundle,
        lr_hadronic_inputs=lr_hadronic_inputs,
    )
    q1_m12_value = _compute_q1_only_m12_value(
        wilsons=resolved_wilsons,
        hadronic_bundle=q1_hadronic_bundle,
    )
    lr_m12_value = _compute_lr_only_m12_value(
        wilsons=resolved_wilsons,
        hadronic_inputs=lr_hadronic_inputs,
    )
    total_m12_value = _require_finite_complex(
        "M12_K_NP_TOTAL_GeV",
        q1_m12_value + lr_m12_value,
    )
    return Paper07101869KaonCustomTotalObservableResult(
        q1_hadronic_bundle=q1_hadronic_bundle,
        lr_hadronic_inputs=lr_hadronic_inputs,
        wilsons=resolved_wilsons,
        M12_K_NP_Q1_GeV=q1_m12_value,
        M12_K_LR_NP_GeV=lr_m12_value,
        M12_K_NP_TOTAL_GeV=total_m12_value,
        delta_m_K_NP_TOTAL_GeV=2.0 * float(total_m12_value.real),
    )


def _compute_custom_b_q1_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869BMesonHadronicBundle | None = None,
    expected_system_id: str,
) -> Paper07101869BMesonCustomQ1ObservableResult:
    if hadronic_bundle is None:
        raise ValueError(
            f"hadronic_bundle must be provided for custom {expected_system_id} "
            "observable evaluation"
        )
    if rg_result_or_wilsons is None:
        raise ValueError(
            f"rg_result_or_wilsons must be provided for custom {expected_system_id} "
            "observable evaluation"
        )
    if not isinstance(hadronic_bundle, Paper07101869BMesonHadronicBundle):
        raise ValueError("hadronic_bundle must be a Paper07101869BMesonHadronicBundle")
    resolved_wilsons = _resolve_rg_wilsons(rg_result_or_wilsons)
    _require_custom_b_q1_hadronic_compatibility(
        wilsons=resolved_wilsons,
        hadronic_bundle=hadronic_bundle,
        expected_system_id=expected_system_id,
    )
    m12_value = _compute_custom_b_q1_m12_value(
        wilsons=resolved_wilsons,
        hadronic_bundle=hadronic_bundle,
        m12_name=f"M12_{expected_system_id}_NP_GeV",
    )
    return Paper07101869BMesonCustomQ1ObservableResult(
        hadronic_bundle=hadronic_bundle,
        wilsons=resolved_wilsons,
        M12_NP_GeV=m12_value,
        delta_m_NP_GeV=2.0 * float(m12_value.real),
    )


def compute_bd_custom_q1_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869BMesonHadronicBundle | None = None,
) -> Paper07101869BMesonCustomQ1ObservableResult:
    """Compute the custom Q1-only B_d observable subset from evolved paper-mode Wilsons."""

    return _compute_custom_b_q1_observables(
        rg_result_or_wilsons,
        hadronic_bundle=hadronic_bundle,
        expected_system_id=PAPER_0710_1869_DELTAF2_BD_HADRONIC_SYSTEM_ID,
    )


def compute_bs_custom_q1_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869BMesonHadronicBundle | None = None,
) -> Paper07101869BMesonCustomQ1ObservableResult:
    """Compute the custom Q1-only B_s observable subset from evolved paper-mode Wilsons."""

    return _compute_custom_b_q1_observables(
        rg_result_or_wilsons,
        hadronic_bundle=hadronic_bundle,
        expected_system_id=PAPER_0710_1869_DELTAF2_BS_HADRONIC_SYSTEM_ID,
    )


def compute_d0_custom_q1_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869D0HadronicBundle | None = None,
) -> Paper07101869D0CustomQ1ObservableResult:
    """Compute the custom Q1-only D0 observable subset from evolved paper-mode Wilsons."""

    if hadronic_bundle is None:
        raise ValueError("hadronic_bundle must be provided for custom D0 observable evaluation")
    if rg_result_or_wilsons is None:
        raise ValueError(
            "rg_result_or_wilsons must be provided for custom D0 observable evaluation"
        )
    if not isinstance(hadronic_bundle, Paper07101869D0HadronicBundle):
        raise ValueError("hadronic_bundle must be a Paper07101869D0HadronicBundle")
    resolved_wilsons = _resolve_rg_wilsons(rg_result_or_wilsons)
    _require_d0_q1_hadronic_compatibility(
        wilsons=resolved_wilsons,
        hadronic_bundle=hadronic_bundle,
    )
    m12_value = _compute_d0_q1_m12_value(
        wilsons=resolved_wilsons,
        hadronic_bundle=hadronic_bundle,
    )
    return Paper07101869D0CustomQ1ObservableResult(
        hadronic_bundle=hadronic_bundle,
        wilsons=resolved_wilsons,
        M12_D0_NP_GeV=m12_value,
        delta_m_D0_NP_GeV=2.0 * float(m12_value.real),
    )


def evaluate_paper_0710_1869_bd_custom_q1_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869BMesonHadronicBundle | None = None,
) -> Paper07101869BMesonCustomQ1ObservableResult:
    """Compatibility alias for the custom Q1-only B_d observable computation."""

    return compute_bd_custom_q1_observables(
        rg_result_or_wilsons,
        hadronic_bundle=hadronic_bundle,
    )


def evaluate_bd_custom_q1_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869BMesonHadronicBundle | None = None,
) -> Paper07101869BMesonCustomQ1ObservableResult:
    """Compatibility alias for the custom Q1-only B_d observable computation."""

    return evaluate_paper_0710_1869_bd_custom_q1_observables(
        rg_result_or_wilsons,
        hadronic_bundle=hadronic_bundle,
    )


def evaluate_paper_0710_1869_bs_custom_q1_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869BMesonHadronicBundle | None = None,
) -> Paper07101869BMesonCustomQ1ObservableResult:
    """Compatibility alias for the custom Q1-only B_s observable computation."""

    return compute_bs_custom_q1_observables(
        rg_result_or_wilsons,
        hadronic_bundle=hadronic_bundle,
    )


def evaluate_bs_custom_q1_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869BMesonHadronicBundle | None = None,
) -> Paper07101869BMesonCustomQ1ObservableResult:
    """Compatibility alias for the custom Q1-only B_s observable computation."""

    return evaluate_paper_0710_1869_bs_custom_q1_observables(
        rg_result_or_wilsons,
        hadronic_bundle=hadronic_bundle,
    )


def evaluate_paper_0710_1869_d0_custom_q1_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869D0HadronicBundle | None = None,
) -> Paper07101869D0CustomQ1ObservableResult:
    """Compatibility alias for the custom Q1-only D0 observable computation."""

    return compute_d0_custom_q1_observables(
        rg_result_or_wilsons,
        hadronic_bundle=hadronic_bundle,
    )


def evaluate_d0_custom_q1_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869D0HadronicBundle | None = None,
) -> Paper07101869D0CustomQ1ObservableResult:
    """Compatibility alias for the custom Q1-only D0 observable computation."""

    return evaluate_paper_0710_1869_d0_custom_q1_observables(
        rg_result_or_wilsons,
        hadronic_bundle=hadronic_bundle,
    )


def evaluate_paper_0710_1869_kaon_lr_only_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_inputs: Paper07101869KaonLRHadronicInputs | None = None,
) -> Paper07101869KaonLROnlyObservableResult:
    """Compatibility alias for the custom LR-only kaon observable computation."""

    return compute_kaon_lr_only_observables(
        rg_result_or_wilsons,
        hadronic_inputs=hadronic_inputs,
    )


def evaluate_kaon_lr_only_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_inputs: Paper07101869KaonLRHadronicInputs | None = None,
) -> Paper07101869KaonLROnlyObservableResult:
    """Compatibility alias for the custom LR-only kaon observable computation."""

    return evaluate_paper_0710_1869_kaon_lr_only_observables(
        rg_result_or_wilsons,
        hadronic_inputs=hadronic_inputs,
    )


def evaluate_paper_0710_1869_kaon_custom_total_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    q1_hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
    lr_hadronic_inputs: Paper07101869KaonLRHadronicInputs | None = None,
) -> Paper07101869KaonCustomTotalObservableResult:
    """Compatibility alias for the custom combined kaon observable computation."""

    return compute_kaon_custom_total_observables(
        rg_result_or_wilsons,
        q1_hadronic_bundle=q1_hadronic_bundle,
        lr_hadronic_inputs=lr_hadronic_inputs,
    )


def evaluate_kaon_custom_total_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    q1_hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
    lr_hadronic_inputs: Paper07101869KaonLRHadronicInputs | None = None,
) -> Paper07101869KaonCustomTotalObservableResult:
    """Compatibility alias for the custom combined kaon observable computation."""

    return evaluate_paper_0710_1869_kaon_custom_total_observables(
        rg_result_or_wilsons,
        q1_hadronic_bundle=q1_hadronic_bundle,
        lr_hadronic_inputs=lr_hadronic_inputs,
    )


def build_paper_0710_1869_kaon_lr_only_observable_result(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    hadronic_inputs: Paper07101869KaonLRHadronicInputs | None = None,
) -> Paper07101869KaonLROnlyObservableResult:
    """Return the deterministic custom LR-only kaon observable result."""

    if rg_result is not None and wilsons is not None:
        raise ValueError("provide either rg_result or wilsons, not both")
    resolved_input = rg_result if rg_result is not None else wilsons
    return compute_kaon_lr_only_observables(
        resolved_input,
        hadronic_inputs=hadronic_inputs,
    )


def build_paper_0710_1869_kaon_lr_only_observable_summary(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    hadronic_inputs: Paper07101869KaonLRHadronicInputs | None = None,
) -> dict[str, object]:
    """Return the deterministic summary mapping for the custom LR-only observable result."""

    return build_paper_0710_1869_kaon_lr_only_observable_result(
        rg_result=rg_result,
        wilsons=wilsons,
        hadronic_inputs=hadronic_inputs,
    ).summary()


def build_paper_0710_1869_kaon_custom_total_observable_result(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    q1_hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
    lr_hadronic_inputs: Paper07101869KaonLRHadronicInputs | None = None,
) -> Paper07101869KaonCustomTotalObservableResult:
    """Return the deterministic custom combined kaon observable result."""

    if rg_result is not None and wilsons is not None:
        raise ValueError("provide either rg_result or wilsons, not both")
    resolved_input = rg_result if rg_result is not None else wilsons
    return compute_kaon_custom_total_observables(
        resolved_input,
        q1_hadronic_bundle=q1_hadronic_bundle,
        lr_hadronic_inputs=lr_hadronic_inputs,
    )


def build_paper_0710_1869_kaon_custom_total_observable_summary(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    q1_hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
    lr_hadronic_inputs: Paper07101869KaonLRHadronicInputs | None = None,
) -> dict[str, object]:
    """Return the deterministic summary mapping for the custom combined kaon observable result."""

    return build_paper_0710_1869_kaon_custom_total_observable_result(
        rg_result=rg_result,
        wilsons=wilsons,
        q1_hadronic_bundle=q1_hadronic_bundle,
        lr_hadronic_inputs=lr_hadronic_inputs,
    ).summary()


def evaluate_paper_0710_1869_kaon_np_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
) -> Paper07101869KaonNPObservableResult:
    """Compatibility alias for the kaon NP-only observable computation."""

    return compute_kaon_np_observables(
        rg_result_or_wilsons,
        hadronic_bundle=hadronic_bundle,
    )


def evaluate_kaon_np_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
) -> Paper07101869KaonNPObservableResult:
    """Compatibility alias for the kaon NP-only observable computation."""

    return evaluate_paper_0710_1869_kaon_np_observables(
        rg_result_or_wilsons,
        hadronic_bundle=hadronic_bundle,
    )


def evaluate_paper_0710_1869_kaon_observables(
    rg_result_or_wilsons: (
        Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot | None
    ) = None,
    *,
    hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
) -> Paper07101869KaonNPObservableResult:
    """Compatibility alias matching the acceptance benchmark observable naming."""

    return compute_kaon_np_observables(
        rg_result_or_wilsons,
        hadronic_bundle=hadronic_bundle,
    )


def build_paper_0710_1869_kaon_np_observable_result(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
) -> Paper07101869KaonNPObservableResult:
    """Return the deterministic kaon NP-only observable result."""

    if rg_result is not None and wilsons is not None:
        raise ValueError("provide either rg_result or wilsons, not both")
    resolved_input: Paper07101869DeltaF2RGResult | Paper07101869DeltaF2RGWilsonSnapshot
    if rg_result is None and wilsons is None:
        resolved_input = run_deltaf2_wilsons_lo(default_paper_0710_1869_kaon_matching().wilsons)
    else:
        resolved_input = rg_result if rg_result is not None else wilsons  # type: ignore[assignment]
    return compute_kaon_np_observables(
        resolved_input,
        hadronic_bundle=hadronic_bundle,
    )


def default_paper_0710_1869_kaon_np_observable_result() -> Paper07101869KaonNPObservableResult:
    """Return the default kaon NP-only observable result."""

    return build_paper_0710_1869_kaon_np_observable_result()


def build_paper_0710_1869_kaon_observables(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
) -> Paper07101869KaonNPObservableResult:
    """Compatibility alias for the deterministic kaon observable result."""

    return build_paper_0710_1869_kaon_np_observable_result(
        rg_result=rg_result,
        wilsons=wilsons,
        hadronic_bundle=hadronic_bundle,
    )


def default_paper_0710_1869_kaon_observables() -> Paper07101869KaonNPObservableResult:
    """Compatibility alias for the default kaon observable result."""

    return default_paper_0710_1869_kaon_np_observable_result()


def build_paper_0710_1869_kaon_np_observable_summary(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
) -> dict[str, object]:
    """Return the deterministic summary mapping for the kaon NP-only observable result."""

    return build_paper_0710_1869_kaon_np_observable_result(
        rg_result=rg_result,
        wilsons=wilsons,
        hadronic_bundle=hadronic_bundle,
    ).summary()


def default_paper_0710_1869_kaon_np_observable_summary() -> dict[str, object]:
    """Return the summary mapping for the default kaon NP-only observable result."""

    return default_paper_0710_1869_kaon_np_observable_result().summary()


def build_paper_0710_1869_kaon_observables_summary(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    hadronic_bundle: Paper07101869KaonHadronicBundle | None = None,
) -> dict[str, object]:
    """Compatibility alias for the deterministic kaon observable summary."""

    return build_paper_0710_1869_kaon_np_observable_summary(
        rg_result=rg_result,
        wilsons=wilsons,
        hadronic_bundle=hadronic_bundle,
    )


def default_paper_0710_1869_kaon_observables_summary() -> dict[str, object]:
    """Compatibility alias for the default kaon observable summary."""

    return default_paper_0710_1869_kaon_np_observable_summary()


def build_paper_0710_1869_bd_custom_q1_observable_result(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    hadronic_bundle: Paper07101869BMesonHadronicBundle | None = None,
) -> Paper07101869BMesonCustomQ1ObservableResult:
    """Return the deterministic custom B_d Q1-only observable result."""

    if rg_result is not None and wilsons is not None:
        raise ValueError("provide either rg_result or wilsons, not both")
    resolved_input = rg_result if rg_result is not None else wilsons
    return compute_bd_custom_q1_observables(
        resolved_input,
        hadronic_bundle=hadronic_bundle,
    )


def build_paper_0710_1869_bd_custom_q1_observable_summary(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    hadronic_bundle: Paper07101869BMesonHadronicBundle | None = None,
) -> dict[str, object]:
    """Return the deterministic summary mapping for the custom B_d Q1-only observable result."""

    return build_paper_0710_1869_bd_custom_q1_observable_result(
        rg_result=rg_result,
        wilsons=wilsons,
        hadronic_bundle=hadronic_bundle,
    ).summary()


def build_paper_0710_1869_bs_custom_q1_observable_result(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    hadronic_bundle: Paper07101869BMesonHadronicBundle | None = None,
) -> Paper07101869BMesonCustomQ1ObservableResult:
    """Return the deterministic custom B_s Q1-only observable result."""

    if rg_result is not None and wilsons is not None:
        raise ValueError("provide either rg_result or wilsons, not both")
    resolved_input = rg_result if rg_result is not None else wilsons
    return compute_bs_custom_q1_observables(
        resolved_input,
        hadronic_bundle=hadronic_bundle,
    )


def build_paper_0710_1869_bs_custom_q1_observable_summary(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    hadronic_bundle: Paper07101869BMesonHadronicBundle | None = None,
) -> dict[str, object]:
    """Return the deterministic summary mapping for the custom B_s Q1-only observable result."""

    return build_paper_0710_1869_bs_custom_q1_observable_result(
        rg_result=rg_result,
        wilsons=wilsons,
        hadronic_bundle=hadronic_bundle,
    ).summary()


def build_paper_0710_1869_d0_custom_q1_observable_result(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    hadronic_bundle: Paper07101869D0HadronicBundle | None = None,
) -> Paper07101869D0CustomQ1ObservableResult:
    """Return the deterministic custom D0 Q1-only observable result."""

    if rg_result is not None and wilsons is not None:
        raise ValueError("provide either rg_result or wilsons, not both")
    resolved_input = rg_result if rg_result is not None else wilsons
    return compute_d0_custom_q1_observables(
        resolved_input,
        hadronic_bundle=hadronic_bundle,
    )


def build_paper_0710_1869_d0_custom_q1_observable_summary(
    *,
    rg_result: Paper07101869DeltaF2RGResult | None = None,
    wilsons: Paper07101869DeltaF2RGWilsonSnapshot | None = None,
    hadronic_bundle: Paper07101869D0HadronicBundle | None = None,
) -> dict[str, object]:
    """Return the deterministic summary mapping for the custom D0 Q1-only observable result."""

    return build_paper_0710_1869_d0_custom_q1_observable_result(
        rg_result=rg_result,
        wilsons=wilsons,
        hadronic_bundle=hadronic_bundle,
    ).summary()


__all__ = [
    "PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_DELTA_M_OBSERVABLE_ID",
    "PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_INTERPRETATION_ID",
    "PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_M12_OBSERVABLE_ID",
    "PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_OBSERVABLE_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_OBSERVABLE_SCOPE_ID",
    "PAPER_0710_1869_DELTAF2_BD_CUSTOM_Q1_OBSERVABLE_SUMMARY_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_DELTA_M_OBSERVABLE_ID",
    "PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_INTERPRETATION_ID",
    "PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_M12_OBSERVABLE_ID",
    "PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_OBSERVABLE_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_OBSERVABLE_SCOPE_ID",
    "PAPER_0710_1869_DELTAF2_BS_CUSTOM_Q1_OBSERVABLE_SUMMARY_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_DELTA_M_OBSERVABLE_ID",
    "PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_INTERPRETATION_ID",
    "PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_M12_OBSERVABLE_ID",
    "PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_OBSERVABLE_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_OBSERVABLE_SCOPE_ID",
    "PAPER_0710_1869_DELTAF2_D0_CUSTOM_Q1_OBSERVABLE_SUMMARY_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_DELTA_M_OBSERVABLE_ID",
    "PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_INTERPRETATION_ID",
    "PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_M12_OBSERVABLE_ID",
    "PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_OBSERVABLE_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_OBSERVABLE_SCOPE_ID",
    "PAPER_0710_1869_DELTAF2_KAON_CUSTOM_TOTAL_OBSERVABLE_SUMMARY_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_DELTA_M_OBSERVABLE_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_DELTA_M_OBSERVABLE_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_INTERPRETATION_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_M12_OBSERVABLE_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_OBSERVABLE_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_OBSERVABLE_SCOPE_ID",
    "PAPER_0710_1869_DELTAF2_KAON_LR_ONLY_OBSERVABLE_SUMMARY_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_M12_OBSERVABLE_ID",
    "PAPER_0710_1869_DELTAF2_KAON_NP_OBSERVABLE_SCHEMA_ID",
    "PAPER_0710_1869_DELTAF2_KAON_NP_INTERPRETATION_ID",
    "PAPER_0710_1869_DELTAF2_KAON_NP_OBSERVABLE_SCOPE_ID",
    "PAPER_0710_1869_DELTAF2_KAON_NP_OBSERVABLE_SUMMARY_SCHEMA_ID",
    "Paper07101869BMesonCustomQ1ObservableResult",
    "Paper07101869D0CustomQ1ObservableResult",
    "Paper07101869KaonCustomTotalObservableResult",
    "Paper07101869KaonLROnlyObservableResult",
    "Paper07101869KaonNPObservableResult",
    "build_paper_0710_1869_bd_custom_q1_observable_result",
    "build_paper_0710_1869_bd_custom_q1_observable_summary",
    "build_paper_0710_1869_bs_custom_q1_observable_result",
    "build_paper_0710_1869_bs_custom_q1_observable_summary",
    "build_paper_0710_1869_d0_custom_q1_observable_result",
    "build_paper_0710_1869_d0_custom_q1_observable_summary",
    "build_paper_0710_1869_kaon_custom_total_observable_result",
    "build_paper_0710_1869_kaon_custom_total_observable_summary",
    "build_paper_0710_1869_kaon_lr_only_observable_result",
    "build_paper_0710_1869_kaon_lr_only_observable_summary",
    "build_paper_0710_1869_kaon_observables",
    "build_paper_0710_1869_kaon_observables_summary",
    "build_paper_0710_1869_kaon_np_observable_result",
    "build_paper_0710_1869_kaon_np_observable_summary",
    "compute_bd_custom_q1_observables",
    "compute_bs_custom_q1_observables",
    "compute_d0_custom_q1_observables",
    "compute_kaon_custom_total_observables",
    "compute_kaon_lr_only_observables",
    "compute_kaon_np_observables",
    "default_paper_0710_1869_kaon_observables",
    "default_paper_0710_1869_kaon_observables_summary",
    "default_paper_0710_1869_kaon_np_observable_result",
    "default_paper_0710_1869_kaon_np_observable_summary",
    "evaluate_bd_custom_q1_observables",
    "evaluate_bs_custom_q1_observables",
    "evaluate_d0_custom_q1_observables",
    "evaluate_kaon_custom_total_observables",
    "evaluate_kaon_lr_only_observables",
    "evaluate_kaon_np_observables",
    "evaluate_paper_0710_1869_bd_custom_q1_observables",
    "evaluate_paper_0710_1869_bs_custom_q1_observables",
    "evaluate_paper_0710_1869_d0_custom_q1_observables",
    "evaluate_paper_0710_1869_kaon_custom_total_observables",
    "evaluate_paper_0710_1869_kaon_lr_only_observables",
    "evaluate_paper_0710_1869_kaon_observables",
    "evaluate_paper_0710_1869_kaon_np_observables",
]
