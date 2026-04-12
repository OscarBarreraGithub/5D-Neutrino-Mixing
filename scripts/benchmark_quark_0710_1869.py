#!/usr/bin/env python3
"""Contract-first acceptance benchmark for the paper_0710_1869 path.

The script imports the paper package through the real
``quarkConstraints.paper_0710_1869`` path, checks for repo-v1 leakage, and
builds a deterministic structural summary from the paper-only modules.
"""

from __future__ import annotations

import argparse
import ast
import dataclasses
import hashlib
import importlib
import inspect
import json
import math
import subprocess
import sys
import tempfile
import textwrap
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
PAPER_PACKAGE_DIR = REPO_ROOT / "quarkConstraints" / "paper_0710_1869"
EXPECTED_MODE_ID = "paper_0710_1869"
REQUIRED_MODULE_IMPORT_ORDER = (
    "validation",
    "conventions",
    "scales",
    "scan",
)
OPTIONAL_MODULE_IMPORT_ORDER = (
    "inputs",
    "benchmarks",
    "model",
    "couplings",
    "kkgluon",
    "artifacts",
    "verifier",
    "eft_deltaf2.operators",
    "eft_deltaf2.matching_kkgluon",
    "eft_deltaf2.rg",
    "eft_deltaf2.hadronic",
    "eft_deltaf2.observables",
)
IMPLEMENTATION_DEFAULT_MATCHING_EXPORT_NAMES = (
    "match_default_paper_0710_1869_kaon_kk_gluon_deltaf2",
)
MATCHING_ALIAS_STEMS = ("kaon_matching", "deltaf2_matching", "kkgluon_matching")
DEFAULT_MATCHING_SUMMARY_NAMES = tuple(
    name
    for stem in MATCHING_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}_summary",
        f"build_paper_0710_1869_{stem}_summary",
    )
)
DEFAULT_MATCHING_EXPORT_NAMES = IMPLEMENTATION_DEFAULT_MATCHING_EXPORT_NAMES + tuple(
    name
    for stem in MATCHING_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}",
        f"build_paper_0710_1869_{stem}",
    )
)
RG_ALIAS_STEMS = ("kaon_rg", "deltaf2_rg", "lo_rg", "rg_lo")
IMPLEMENTATION_DEFAULT_RG_SUMMARY_NAMES = (
    "default_paper_0710_1869_default_kaon_rg_summary",
    "build_paper_0710_1869_default_kaon_rg_summary",
)
DEFAULT_RG_SUMMARY_NAMES = tuple(
    name
    for stem in RG_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}_summary",
        f"build_paper_0710_1869_{stem}_summary",
    )
) + IMPLEMENTATION_DEFAULT_RG_SUMMARY_NAMES
DEFAULT_RG_EXPORT_NAMES = tuple(
    name
    for stem in RG_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}",
        f"build_paper_0710_1869_{stem}",
    )
)
PARAMETERIZED_RG_EXPORT_NAMES = (
    "evolve_deltaf2_wilsons_lo",
    "run_deltaf2_wilsons_lo",
    "run_paper_0710_1869_deltaf2_wilsons_lo",
    "evolve_paper_0710_1869_deltaf2_wilsons_lo",
)
HADRONIC_ALIAS_STEMS = (
    "default_kaon_hadronic",
    "kaon_hadronic",
    "hadronic_inputs",
    "hadronic",
)
DEFAULT_HADRONIC_SUMMARY_NAMES = tuple(
    name
    for stem in HADRONIC_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}_summary",
        f"build_paper_0710_1869_{stem}_summary",
    )
)
DEFAULT_HADRONIC_EXPORT_NAMES = tuple(
    name
    for stem in HADRONIC_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}",
        f"build_paper_0710_1869_{stem}",
    )
)
HADRONIC_BUILDER_NAMES = (
    "build_paper_0710_1869_kaon_hadronic_bundle",
    "build_paper_0710_1869_kaon_hadronic",
)
CUSTOM_LR_HADRONIC_BUILDER_NAMES = (
    "build_paper_0710_1869_custom_kaon_lr_hadronic_bundle",
    "build_paper_0710_1869_kaon_lr_hadronic_bundle",
    "build_paper_0710_1869_custom_kaon_lr_hadronic",
    "build_paper_0710_1869_kaon_lr_hadronic",
)
DEFAULT_LR_HADRONIC_EXPORT_NAMES = (
    "default_paper_0710_1869_kaon_lr_hadronic_inputs",
    "default_paper_0710_1869_kaon_lr_hadronic_bundle",
    "default_paper_0710_1869_kaon_lr_hadronic",
    "default_paper_0710_1869_kaon_lr_hadronic_summary",
)
DEFAULT_LR_HADRONIC_VALUE_EXPORT_NAMES = DEFAULT_LR_HADRONIC_EXPORT_NAMES[:-1]
KAON_LR_R_CHI_FREEZE_EXPORT_NAMES = (
    "default_paper_0710_1869_kaon_lr_r_chi_freeze",
    "build_paper_0710_1869_kaon_lr_r_chi_freeze",
)
KAON_LR_R_CHI_SUMMARY_EXPORT_NAMES = (
    "default_paper_0710_1869_kaon_lr_r_chi_summary",
    "build_paper_0710_1869_kaon_lr_r_chi_summary",
)
HADRONIC_SOURCE_REF_CLASS_NAMES = (
    "Paper07101869LRHadronicSourceRef",
    "Paper07101869HadronicSourceRef",
)
LR_HADRONIC_B4_PARAM_NAMES = ("B4_mu_had", "b4_mu_had", "B4", "b4")
LR_HADRONIC_B5_PARAM_NAMES = ("B5_mu_had", "b5_mu_had", "B5", "b5")
LR_HADRONIC_R_CHI_PARAM_NAMES = ("R_chi_mu_had", "r_chi_mu_had", "R_chi", "r_chi")
LR_HADRONIC_B4_SOURCE_PARAM_NAMES = (
    "B4_source",
    "b4_source",
    "B4_source_ref",
    "b4_source_ref",
)
LR_HADRONIC_B5_SOURCE_PARAM_NAMES = (
    "B5_source",
    "b5_source",
    "B5_source_ref",
    "b5_source_ref",
)
LR_HADRONIC_R_CHI_SOURCE_PARAM_NAMES = (
    "R_chi_source",
    "r_chi_source",
    "R_chi_source_ref",
    "r_chi_source_ref",
)
DEFAULT_HADRONIC_FORBIDDEN_LR_FIELDS = (
    "q4_matrix_element_GeV4",
    "q5_matrix_element_GeV4",
    "B4_mu_had",
    "B5_mu_had",
    "R_chi_mu_had",
)
LR_HADRONIC_PROBE_B4 = 0.92
LR_HADRONIC_PROBE_B5 = 0.71
LR_HADRONIC_PROBE_R_CHI = 31.5
LR_HADRONIC_PROBE_SCHEME_ID = "paper_0710_1869.deltaf2.kaon_lr_hadronic.custom_probe.v1"
LR_HADRONIC_PROBE_MU_HAD_GEV = 3.25
LR_OBSERVABLE_PROBE_Q4 = complex(1.25e-12, -4.0e-13)
LR_OBSERVABLE_PROBE_Q5 = complex(-7.5e-13, 3.0e-13)
LR_HADRONIC_PROBE_BUNDLE_ID = "hadronic.kaon.lr_custom_probe.v1"
LR_HADRONIC_PROBE_SOURCE_ID = "hadronic.kaon.lr_custom_probe.sources.v1"
LR_HADRONIC_PROBE_B4_SOURCE_ID = "hadronic.kaon.lr_custom_probe.b4.v1"
LR_HADRONIC_PROBE_B5_SOURCE_ID = "hadronic.kaon.lr_custom_probe.b5.v1"
LR_HADRONIC_PROBE_R_CHI_SOURCE_ID = "hadronic.kaon.lr_custom_probe.rchi.v1"
OBSERVABLE_ALIAS_STEMS = (
    "default_kaon_observables",
    "kaon_np_observable",
    "kaon_np_observables",
    "kaon_np_observable_result",
    "kaon_observables",
    "kaon_observable",
    "observables",
    "observable",
)
DEFAULT_OBSERVABLE_SUMMARY_NAMES = tuple(
    name
    for stem in OBSERVABLE_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}_summary",
        f"build_paper_0710_1869_{stem}_summary",
    )
)
DEFAULT_OBSERVABLE_EXPORT_NAMES = tuple(
    name
    for stem in OBSERVABLE_ALIAS_STEMS
    for name in (
        f"default_paper_0710_1869_{stem}",
        f"build_paper_0710_1869_{stem}",
    )
)
PARAMETERIZED_OBSERVABLE_EXPORT_NAMES = (
    "evaluate_paper_0710_1869_kaon_np_observables",
    "evaluate_paper_0710_1869_kaon_observables",
    "evaluate_deltaf2_kaon_observables",
    "evaluate_paper_0710_1869_observables",
    "build_paper_0710_1869_kaon_np_observable_result",
    "compute_paper_0710_1869_kaon_observables",
    "compute_kaon_np_observables",
    "evaluate_kaon_np_observables",
)
CUSTOM_LR_OBSERVABLE_EXPORT_NAMES = (
    "evaluate_paper_0710_1869_kaon_lr_only_observables",
    "evaluate_kaon_lr_only_observables",
    "build_paper_0710_1869_kaon_lr_only_observable_result",
    "compute_kaon_lr_only_observables",
)
CUSTOM_TOTAL_OBSERVABLE_EXPORT_NAMES = (
    "evaluate_paper_0710_1869_kaon_custom_total_observables",
    "evaluate_kaon_custom_total_observables",
    "build_paper_0710_1869_kaon_custom_total_observable_result",
    "compute_kaon_custom_total_observables",
)
RG_RESULT_PARAM_NAMES = (
    "rg_result_or_wilsons",
    "rg_result",
    "deltaf2_rg_result",
    "deltaf2_result",
    "result",
    "rg_output",
)
WILSON_PARAM_NAMES = (
    "wilsons",
    "deltaf2_wilsons",
    "evolved_wilsons",
    "coefficients",
)
HADRONIC_PARAM_NAMES = (
    "hadronic",
    "hadronic_inputs",
    "hadronic_bundle",
    "input_bundle",
    "matrix_elements",
)
SUPPORTED_OBSERVABLE_ROWS = {
    "M12_K_NP.re",
    "M12_K_NP.im",
    "Delta_m_K_NP",
}
SUPPORTED_OBSERVABLE_CONTAINER_KEYS = {"M12_K_NP", "Delta_m_K_NP"}
SUPPORTED_OBSERVABLE_FLAT_FIELDS = {
    "re_M12_K_NP_GeV": "M12_K_NP.re",
    "im_M12_K_NP_GeV": "M12_K_NP.im",
    "delta_m_K_NP_GeV": "Delta_m_K_NP",
}
EXPECTED_OBSERVABLE_SCOPE_ID = "kaon.np_only.m12.v1"
EXPECTED_INTERPRETATION_ID = "kaon.np_only.v1"
EXPECTED_M12_OBSERVABLE_ID = "M12_K_NP"
EXPECTED_DELTA_M_OBSERVABLE_ID = "Delta_m_K_NP"
EXPECTED_KAON_LR_R_CHI_MASS_SCHEME_ID = "pdg.2024.msbar.running_masses.at_2gev.v1"
EXPECTED_KAON_LR_R_CHI_ACTIVE_FLAVOR_POLICY_ID = (
    "pdg.2024.quark_masses.n_l_4.at_2gev.explicit.v1"
)
EXPECTED_KAON_LR_R_CHI_NO_HIDDEN_CONVERSION_POLICY_ID = (
    "none.freeze_pdg2024_msbar_nl4_inputs_at_2gev.v1"
)
EXPECTED_KAON_LR_R_CHI_M_S_2GEV_GEV = 0.09274
EXPECTED_KAON_LR_R_CHI_M_D_2GEV_GEV = 0.00469
EXPECTED_KAON_LR_R_CHI_M_K0_GEV = 0.497611
EXPECTED_KAON_LR_R_CHI_EXACT_VALUE = 26.085222120747908
EXPECTED_KAON_LR_R_CHI_FREEZE_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_lr_r_chi_freeze.v1"
)
EXPECTED_KAON_LR_R_CHI_SUMMARY_SCHEMA_ID = (
    "quarkConstraints.paper_0710_1869.eft_deltaf2.kaon_lr_r_chi_summary.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_BUNDLE_ID = "hadronic.kaon.lr.default.etm2013_ms_2gev.v1"
EXPECTED_DEFAULT_LR_HADRONIC_SOURCE_ID = (
    "hadronic.kaon.lr.default.etm2013_ms_2gev.aggregate.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_INPUT_POLICY_ID = (
    "default_source.etm2013.table1.ms_2gev.no_hidden_conversion.v1"
)
EXPECTED_CUSTOM_LR_HADRONIC_INPUT_POLICY_ID = (
    "custom_input_only.default_lr_bundle_frozen_separately.no_auto_consumption.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_B4_SOURCE_ID = (
    "hadronic.kaon.lr.b4.etm2013.table1.ms_2gev.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_B5_SOURCE_ID = (
    "hadronic.kaon.lr.b5.etm2013.table1.ms_2gev.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_R_CHI_SOURCE_ID = "hadronic.kaon.lr.r_chi.pdg2024_msbar_nl4.v1"
EXPECTED_DEFAULT_LR_HADRONIC_MASS_SOURCE_ID = "pdg.2024.k0.mass.v1"
EXPECTED_DEFAULT_LR_HADRONIC_DECAY_CONSTANT_SOURCE_ID = "pdg.2024.fkplus.eq72.14.v1"
EXPECTED_DEFAULT_LR_OPERATOR_BASIS_ID = "kk_gluon_tree_np_only.v1"
EXPECTED_DEFAULT_LR_OPERATOR_NORMALIZATION_ID = (
    "paper_0710_1869.deltaf2.kk_gluon_tree_color_normalization.v1"
)
EXPECTED_DEFAULT_LR_HAMILTONIAN_CONVENTION_ID = "heff.sum_ci_qi.no_hc_factor.v1"
EXPECTED_DEFAULT_LR_Q4_FORMULA_ID = (
    "kaon.q4_lr.o4_scalar_lr.bv2004.eq5.matrix_element.mu_had.v1"
)
EXPECTED_DEFAULT_LR_Q5_FORMULA_ID = (
    "kaon.q5_lr.o5_scalar_lr.bv2004.eq5.matrix_element.mu_had.v1"
)
EXPECTED_DEFAULT_LR_HADRONIC_ETM_CITATION = (
    "ETM Collaboration, JHEP 03 (2013) 089, arXiv:1207.1287"
)
EXPECTED_DEFAULT_LR_HADRONIC_B4_VALUE = 0.78
EXPECTED_DEFAULT_LR_HADRONIC_B5_VALUE = 0.57
EXPECTED_DEFAULT_LR_HADRONIC_PROVENANCE_IDS = [
    EXPECTED_DEFAULT_LR_HADRONIC_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_B4_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_B5_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_R_CHI_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_MASS_SOURCE_ID,
    EXPECTED_DEFAULT_LR_HADRONIC_DECAY_CONSTANT_SOURCE_ID,
]
EXPECTED_DEFAULT_LR_SUPPORTED_OPERATORS = {"Q4_LR", "Q5_LR"}
EXPECTED_DEFAULT_LR_UNSUPPORTED_OPERATORS = {"Q1_VLL", "Q1_VRR"}
EXPECTED_CUSTOM_TOTAL_SCOPE_ID = "kaon.np_only.custom_total.q1_plus_lr.v1"
EXPECTED_CUSTOM_TOTAL_INTERPRETATION_ID = "kaon.np_only.custom_total.q1_plus_lr.v1"
EXPECTED_CUSTOM_TOTAL_M12_OBSERVABLE_ID = "M12_K_NP_CUSTOM_TOTAL"
EXPECTED_CUSTOM_TOTAL_DELTA_M_OBSERVABLE_ID = "Delta_m_K_NP_CUSTOM_TOTAL"
B_CUSTOM_Q1_PROBE_MU_HAD_GEV = 4.2
B_CUSTOM_Q1_PROBE_CONFIG = {
    "B_d": {
        "matching_names": (
            "default_paper_0710_1869_bd_matching",
            "build_paper_0710_1869_bd_matching",
            "match_default_paper_0710_1869_bd_kk_gluon_deltaf2",
        ),
        "hadronic_builder_names": (
            "build_paper_0710_1869_bd_hadronic_bundle",
            "build_paper_0710_1869_bd_hadronic",
        ),
        "observable_names": (
            "evaluate_paper_0710_1869_bd_custom_q1_observables",
            "evaluate_bd_custom_q1_observables",
            "build_paper_0710_1869_bd_custom_q1_observable_result",
            "compute_bd_custom_q1_observables",
        ),
        "scope_id": "bd.np_only.custom_q1.m12.v1",
        "interpretation_id": "bd.np_only.custom_q1.v1",
        "m12_id": "M12_Bd_NP_CUSTOM_Q1",
        "delta_m_id": "Delta_m_Bd_NP_CUSTOM_Q1",
        "mass_param": "m_Bd_GeV",
        "decay_param": "f_Bd_GeV",
        "bag_param": "B_Bd_mu_had",
        "meson_mass_GeV": 5.27965,
        "meson_decay_constant_GeV": 0.190,
        "bag_parameter_mu_had": 0.92,
        "source_prefix": "hadronic.bd.custom_q1_probe",
    },
    "B_s": {
        "matching_names": (
            "default_paper_0710_1869_bs_matching",
            "build_paper_0710_1869_bs_matching",
            "match_default_paper_0710_1869_bs_kk_gluon_deltaf2",
        ),
        "hadronic_builder_names": (
            "build_paper_0710_1869_bs_hadronic_bundle",
            "build_paper_0710_1869_bs_hadronic",
        ),
        "observable_names": (
            "evaluate_paper_0710_1869_bs_custom_q1_observables",
            "evaluate_bs_custom_q1_observables",
            "build_paper_0710_1869_bs_custom_q1_observable_result",
            "compute_bs_custom_q1_observables",
        ),
        "scope_id": "bs.np_only.custom_q1.m12.v1",
        "interpretation_id": "bs.np_only.custom_q1.v1",
        "m12_id": "M12_Bs_NP_CUSTOM_Q1",
        "delta_m_id": "Delta_m_Bs_NP_CUSTOM_Q1",
        "mass_param": "m_Bs_GeV",
        "decay_param": "f_Bs_GeV",
        "bag_param": "B_Bs_mu_had",
        "meson_mass_GeV": 5.36688,
        "meson_decay_constant_GeV": 0.230,
        "bag_parameter_mu_had": 0.91,
        "source_prefix": "hadronic.bs.custom_q1_probe",
    },
    "D0": {
        "matching_names": (
            "default_paper_0710_1869_d0_matching",
            "build_paper_0710_1869_d0_matching",
            "match_default_paper_0710_1869_d0_kk_gluon_deltaf2",
        ),
        "hadronic_builder_names": (
            "build_paper_0710_1869_d0_hadronic_bundle",
            "build_paper_0710_1869_d0_hadronic",
        ),
        "observable_names": (
            "evaluate_paper_0710_1869_d0_custom_q1_observables",
            "evaluate_d0_custom_q1_observables",
            "build_paper_0710_1869_d0_custom_q1_observable_result",
            "compute_d0_custom_q1_observables",
        ),
        "scope_id": "d0.np_only.custom_q1.m12.v1",
        "interpretation_id": "d0.np_only.custom_q1.v1",
        "m12_id": "M12_D0_NP",
        "delta_m_id": "Delta_m_D0_NP",
        "mass_param": "m_D0_GeV",
        "decay_param": "f_D0_GeV",
        "bag_param": "B_D0_mu_had",
        "meson_mass_GeV": 1.86484,
        "meson_decay_constant_GeV": 0.212,
        "bag_parameter_mu_had": 0.80,
        "probe_q1_vll": complex(3.25e-13, -1.10e-13),
        "probe_q1_vrr": complex(-4.50e-14, 2.50e-14),
        "require_nonzero_probe": True,
        "source_prefix": "hadronic.d0.custom_q1_probe",
    },
}
EXPECTED_SUPPORTED_OPERATORS = {"Q1_VLL", "Q1_VRR"}
EXPECTED_UNSUPPORTED_OPERATORS = {"Q4_LR", "Q5_LR"}
EXPECTED_SUPPORTED_OPERATOR_ORDER = ("Q1_VLL", "Q1_VRR")
EXPECTED_GUARDED_LR_OPERATOR_ORDER = ("Q4_LR", "Q5_LR")
EXPECTED_BMU_LR_OPERATOR_ORDER = ("Q1_LR_BMU", "Q2_LR_BMU")
EXPECTED_PAPER_TO_BMU_OPERATOR_MAP = ((0.0, 2.0), (1.0, 0.0))
EXPECTED_BMU_TO_PAPER_OPERATOR_MAP = ((0.0, 1.0), (0.5, 0.0))
EXPECTED_PAPER_TO_BMU_WILSON_MAP = ((0.0, 0.5), (1.0, 0.0))
EXPECTED_BMU_TO_PAPER_WILSON_MAP = ((0.0, 1.0), (2.0, 0.0))
EXPECTED_ARTIFACT_OBSERVABLE_ROWS = {
    "M12_K_NP.re",
    "M12_K_NP.im",
    "Delta_m_K_NP",
}
EXPECTED_M12_K_NP = {
    "real": -1.3495753042583394e-14,
    "imag": -9.322420653906486e-15,
}
EXPECTED_DELTA_M_K_NP_GEV = -2.6991506085166787e-14
ARTIFACT_EXPORT_FILENAMES = {
    "wilson": "wilsons.json",
    "hadronic": "hadronic.json",
    "observable": "observables.json",
    "provenance": "provenance.json",
}
DEFAULT_MATCHING_OBJECT_NAMES = (
    "default_paper_0710_1869_kaon_matching",
    "build_paper_0710_1869_kaon_matching",
    "default_paper_0710_1869_deltaf2_matching",
    "build_paper_0710_1869_deltaf2_matching",
    "default_paper_0710_1869_kkgluon_matching",
    "build_paper_0710_1869_kkgluon_matching",
    "match_default_paper_0710_1869_kaon_kk_gluon_deltaf2",
)
DEFAULT_RG_LOW_SCALE_GEV = 2.0
CUSTOM_Q1_HADRONIC_PROBE_BAG_SOURCE_ID = "hadronic.kaon.custom_total_probe.bk.v1"
CUSTOM_Q1_HADRONIC_PROBE_BAG_SOURCE_CITATION = "Custom caller-supplied B_K(mu_had) source"
CUSTOM_Q1_HADRONIC_PROBE_BAG_SOURCE_LOCATOR = "B_K(mu_had) external input"
CUSTOM_Q1_HADRONIC_PROBE_YEAR = 2026
CUSTOM_Q1_HADRONIC_PROBE_NOTES = "Custom Q1 hadronic source for combined-observable acceptance."
REQUIRED_SCALE_FIELDS = (
    "Lambda_IR_GeV",
    "m_g1_GeV",
    "xi_g",
    "mu_match_GeV",
    "mu_gs_GeV",
    "m_KK_eff_GeV",
)
FORBIDDEN_MODULES = {
    "quarkConstraints.benchmarks",
    "quarkConstraints.couplings",
    "quarkConstraints.deltaf2",
    "quarkConstraints.fit",
    "quarkConstraints.model",
    "quarkConstraints.proxies",
    "quarkConstraints.scan",
    "quarkConstraints.scales",
    "quarkConstraints.validation",
}
FORBIDDEN_ROOT_EXPORTS = {
    "BulkMassMap",
    "DeltaF2ConstraintSummary",
    "DeltaF2Input",
    "DeltaF2ObservableSummary",
    "DeltaF2WilsonSet",
    "compute_delta_f2_wilsons",
    "compute_quark_kk_gluon_couplings",
    "default_delta_f2_inputs",
    "evaluate_delta_f2_constraints",
}

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--require-package",
        action="store_true",
        help="Fail if quarkConstraints/paper_0710_1869 is not present.",
    )
    parser.add_argument(
        "--emit-json",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--emit-default-matching-json",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--emit-default-rg-json",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--emit-default-hadronic-json",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--emit-default-lr-hadronic-json",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--emit-default-observable-json",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--emit-custom-bd-observable-probe-json",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--emit-custom-bs-observable-probe-json",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--emit-custom-d0-observable-probe-json",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--emit-lr-r-chi-freeze-json",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--export-artifacts-dir",
        type=Path,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--verify-artifacts-dir",
        type=Path,
        help=argparse.SUPPRESS,
    )
    return parser.parse_args()


def _py_files_under(path: Path) -> list[Path]:
    return sorted(file for file in path.rglob("*.py") if file.is_file())


def _read_ast(path: Path) -> ast.AST:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _canonicalize(value: Any) -> Any:
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _canonicalize(dataclasses.asdict(value))
    if isinstance(value, Mapping):
        items = sorted(value.items(), key=lambda item: str(item[0]))
        return {str(key): _canonicalize(val) for key, val in items}
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    if isinstance(value, float):
        return float(value)
    if hasattr(value, "tolist"):
        try:
            return _canonicalize(value.tolist())
        except Exception:
            pass
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_canonicalize(item) for item in value]
    if hasattr(value, "__dict__") and value.__dict__:
        return _canonicalize(vars(value))
    return repr(value)


def _dump_canonical_json(payload: Any) -> str:
    return json.dumps(_canonicalize(payload), indent=2, sort_keys=True) + "\n"


def _canonical_payload_json_and_sha256(payload: Any) -> tuple[str, str]:
    canonical_json = _dump_canonical_json(payload)
    return canonical_json, hashlib.sha256(canonical_json.encode("utf-8")).hexdigest()


def _write_canonical_json(path: Path, payload: Any) -> Path:
    path.write_text(_dump_canonical_json(payload), encoding="utf-8")
    return path


def _sha256_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _build_default_artifact_exports(modules: Mapping[str, Any]) -> dict[str, Any]:
    artifact_module = modules["artifacts"]
    export_set = artifact_module.build_default_paper_0710_1869_kaon_artifact_export_set()
    return {
        "point_id": export_set.wilson_bundle.metadata.point_id,
        "bundle_ids": {
            "wilsons": export_set.wilson_bundle.metadata.bundle_id,
            "hadronic": export_set.hadronic_bundle.metadata.bundle_id,
            "observables": export_set.observable_bundle.metadata.bundle_id,
            "provenance": export_set.provenance_bundle.metadata.bundle_id,
        },
        "wilson_bundle": export_set.wilson_bundle,
        "hadronic_bundle": export_set.hadronic_bundle,
        "observable_bundle": export_set.observable_bundle,
        "provenance_bundle": export_set.provenance_bundle,
    }


def _default_artifact_file_paths(
    artifact_module: Any,
    destination: Path,
) -> dict[str, Path]:
    path_factory = getattr(
        artifact_module,
        "default_paper_0710_1869_kaon_artifact_export_paths",
        None,
    )
    if callable(path_factory):
        resolved = path_factory(destination)
        return {
            "wilson": resolved.wilson_path,
            "hadronic": resolved.hadronic_path,
            "observable": resolved.observable_path,
            "provenance": resolved.provenance_path,
        }
    return {
        name: destination / filename
        for name, filename in ARTIFACT_EXPORT_FILENAMES.items()
    }


def _export_default_artifacts(
    modules: Mapping[str, Any],
    destination: Path,
) -> dict[str, Any]:
    destination.mkdir(parents=True, exist_ok=True)
    exports = _build_default_artifact_exports(modules)
    artifact_module = modules["artifacts"]
    file_paths = _default_artifact_file_paths(artifact_module, destination)
    exports["wilson_bundle"].write_json(file_paths["wilson"])
    exports["hadronic_bundle"].write_json(file_paths["hadronic"])
    exports["observable_bundle"].write_json(file_paths["observable"])
    exports["provenance_bundle"].write_json(file_paths["provenance"])
    return {
        **exports,
        "file_paths": file_paths,
        "file_sha256": {
            name: _sha256_digest(path)
            for name, path in file_paths.items()
        },
    }


def _run_independent_verifier(
    *,
    artifact_dir: Path | None = None,
    wilson_path: Path | None = None,
    hadronic_path: Path | None = None,
    observable_path: Path | None = None,
    provenance_path: Path | None = None,
) -> dict[str, Any]:
    artifact_root_literal = repr(str(artifact_dir) if artifact_dir is not None else None)
    wilson_path_literal = repr(str(wilson_path) if wilson_path is not None else None)
    hadronic_path_literal = repr(str(hadronic_path) if hadronic_path is not None else None)
    observable_path_literal = repr(
        str(observable_path) if observable_path is not None else None
    )
    provenance_path_literal = repr(
        str(provenance_path) if provenance_path is not None else None
    )
    script = textwrap.dedent(
        f"""
        import ast
        import dataclasses
        import importlib.util
        import json
        import pathlib
        import sys
        import types

        REPO_ROOT = pathlib.Path({str(REPO_ROOT)!r})
        PACKAGE_DIR = REPO_ROOT / "quarkConstraints" / "paper_0710_1869"
        PACKAGE_NAME = "_paper07101869_standalone"
        EXPECTED_OBSERVABLE_ROWS = {sorted(EXPECTED_ARTIFACT_OBSERVABLE_ROWS)!r}
        EXPECTED_SUPPORTED_OPERATORS = {sorted(EXPECTED_SUPPORTED_OPERATORS)!r}
        EXPECTED_UNSUPPORTED_OPERATORS = {sorted(EXPECTED_UNSUPPORTED_OPERATORS)!r}
        EXPECTED_FILENAMES = {ARTIFACT_EXPORT_FILENAMES!r}
        ALLOWED_RELATIVE_IMPORTS = {{"artifacts", "conventions"}}
        ALLOWED_ABSOLUTE_LOCAL_IMPORTS = {{
            "quarkConstraints.paper_0710_1869.artifacts",
            "quarkConstraints.paper_0710_1869.conventions",
        }}

        package = types.ModuleType(PACKAGE_NAME)
        package.__path__ = [str(PACKAGE_DIR)]
        sys.modules[PACKAGE_NAME] = package

        def load_module(module_basename: str):
            module_path = PACKAGE_DIR / f"{{module_basename}}.py"
            spec = importlib.util.spec_from_file_location(
                f"{{PACKAGE_NAME}}.{{module_basename}}",
                module_path,
            )
            if spec is None or spec.loader is None:
                raise RuntimeError(f"failed to load {{module_path}}")
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
            return module

        def _scan_verifier_imports() -> list[str]:
            source_path = PACKAGE_DIR / "verifier.py"
            tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
            unexpected: list[str] = []
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        module_name = alias.name
                        if module_name == "quarkConstraints" or module_name.startswith(
                            "quarkConstraints."
                        ):
                            if module_name not in ALLOWED_ABSOLUTE_LOCAL_IMPORTS:
                                unexpected.append(f"import {{module_name}}")
                elif isinstance(node, ast.ImportFrom):
                    module_name = node.module or ""
                    if node.level == 1:
                        if module_name not in ALLOWED_RELATIVE_IMPORTS:
                            unexpected.append(f"from .{{module_name}} import ...")
                    elif node.level > 1:
                        prefix = "." * node.level
                        unexpected.append(f"from {{prefix}}{{module_name}} import ...")
                    elif module_name == "quarkConstraints" or module_name.startswith(
                        "quarkConstraints."
                    ):
                        if module_name not in ALLOWED_ABSOLUTE_LOCAL_IMPORTS:
                            unexpected.append(f"from {{module_name}} import ...")
            return sorted(set(unexpected))

        artifacts = load_module("artifacts")
        verifier = load_module("verifier")

        forbidden_roots = {{
            "quarkConstraints.deltaf2",
            "quarkConstraints.couplings",
            "quarkConstraints.model",
            "quarkConstraints.paper_0710_1869.couplings",
            "quarkConstraints.paper_0710_1869.kkgluon",
            "quarkConstraints.paper_0710_1869.model",
            "quarkConstraints.paper_0710_1869.eft_deltaf2",
        }}

        loaded_forbidden = sorted(
            name
            for name in sys.modules
            if any(
                name == root or name.startswith(f"{{root}}.")
                for root in forbidden_roots
            )
        )
        unexpected_import_targets = _scan_verifier_imports()
        import_isolation_ok = not loaded_forbidden and not unexpected_import_targets

        def _emit_payload(payload: dict[str, object]) -> None:
            print(json.dumps(payload, sort_keys=True))

        try:
            artifact_root = {artifact_root_literal}
            if artifact_root is not None:
                resolved = artifacts.default_paper_0710_1869_kaon_artifact_export_paths(
                    pathlib.Path(artifact_root)
                )
                wilson_path = pathlib.Path(resolved.wilson_path)
                hadronic_path = pathlib.Path(resolved.hadronic_path)
                observable_path = pathlib.Path(resolved.observable_path)
                provenance_path = pathlib.Path(resolved.provenance_path)
            else:
                wilson_path = pathlib.Path({wilson_path_literal})
                hadronic_path = pathlib.Path({hadronic_path_literal})
                observable_path = pathlib.Path({observable_path_literal})
                provenance_path = pathlib.Path({provenance_path_literal})

            inputs = verifier.read_verifier_inputs(
                wilson_path,
                hadronic_path,
                observable_path,
                provenance_path,
            )
            wilson_bundle = inputs.wilson_bundle
            hadronic_bundle = inputs.hadronic_bundle
            observable_bundle = inputs.observable_bundle
            provenance_bundle = inputs.provenance_bundle

            report = verifier.verify_inputs(inputs)
            issue_codes = [issue.code for issue in report.issues]
            issues = [
                {{
                    "code": issue.code,
                    "message": issue.message,
                    "subject": getattr(issue, "subject", None),
                }}
                for issue in report.issues
            ]

            if unexpected_import_targets:
                issue_codes.append("import_isolation_static_violation")
                issues.append(
                    {{
                        "code": "import_isolation_static_violation",
                        "message": (
                            "verifier source imports modules outside the allowed "
                            "artifact-only boundary"
                        ),
                        "subject": None,
                    }}
                )
            if loaded_forbidden:
                issue_codes.append("import_isolation_runtime_violation")
                issues.append(
                    {{
                        "code": "import_isolation_runtime_violation",
                        "message": "standalone verifier loaded forbidden paper or repo-v1 modules",
                        "subject": None,
                    }}
                )

            observable_records = {{
                record.name: record
                for record in observable_bundle.observables
            }}
            observable_values = {{
                name: float(record.value)
                for name, record in observable_records.items()
            }}
            def _serialize_numeric_check(check):
                return {{
                    "expected": float(check.expected),
                    "actual": (
                        float(check.actual)
                        if check.actual is not None
                        else None
                    ),
                    "abs_diff": (
                        float(check.abs_diff)
                        if check.abs_diff is not None
                        else None
                    ),
                    "rel_diff": (
                        float(check.rel_diff)
                        if check.rel_diff is not None
                        else None
                    ),
                    "abs_tol": float(check.abs_tol),
                    "rel_tol": float(check.rel_tol),
                    "expected_units": check.expected_units,
                    "actual_units": check.actual_units,
                    "present": bool(check.present),
                    "issue_code": check.issue_code,
                    "ok": bool(check.ok),
                }}

            numeric_checks = {{
                check.name: _serialize_numeric_check(check)
                for check in report.numeric_checks
            }}
            tolerance_policy = dataclasses.asdict(report.tolerance_policy)
            recomputed_m12 = report.recomputed_m12_GeV.to_complex()

            numeric_issue_codes = {{
                "q1_matrix_element_mismatch",
                "m12_reconstruction_mismatch",
                "m12_imag_reconstruction_mismatch",
                "delta_m_reconstruction_mismatch",
            }}
            scope_issue_codes = {{
                "supported_operator_subset_mismatch",
                "unsupported_operator_subset_mismatch",
                "observable_supported_operator_subset_mismatch",
                "observable_unsupported_operator_subset_mismatch",
                "lr_coefficients_nonzero",
                "observable_row_missing",
                "unexpected_observable_row",
                "required_wilson_operator_missing",
                "unexpected_wilson_operator",
            }}
            import_issue_codes = {{
                "import_isolation_static_violation",
                "import_isolation_runtime_violation",
            }}
            structural_issue_codes = sorted(
                set(issue_codes) - numeric_issue_codes - scope_issue_codes - import_issue_codes
            )
            scope_ok = not any(code in scope_issue_codes for code in issue_codes)
            numeric_match_ok = all(bool(entry.get("ok")) for entry in numeric_checks.values())
            schema_ok = not structural_issue_codes

            _emit_payload(
                {{
                    "ok": schema_ok and scope_ok and numeric_match_ok and import_isolation_ok,
                    "import_isolation_ok": import_isolation_ok,
                    "schema_ok": schema_ok,
                    "numeric_match_ok": numeric_match_ok,
                    "scope_ok": scope_ok,
                    "tolerance_policy": tolerance_policy,
                    "tolerances": tolerance_policy,
                    "bundle_ids": {{
                        "wilsons": wilson_bundle.metadata.bundle_id,
                        "hadronic": hadronic_bundle.metadata.bundle_id,
                        "observables": observable_bundle.metadata.bundle_id,
                        "provenance": provenance_bundle.metadata.bundle_id,
                    }},
                    "point_id": wilson_bundle.metadata.point_id,
                    "issue_codes": issue_codes,
                    "issues": issues,
                    "numeric_checks": numeric_checks,
                    "numeric_diffs": {{
                        name: {{
                            "abs_diff": entry["abs_diff"],
                            "rel_diff": entry["rel_diff"],
                            "abs_tol": entry["abs_tol"],
                            "rel_tol": entry["rel_tol"],
                            "issue_code": entry["issue_code"],
                            "ok": entry["ok"],
                        }}
                        for name, entry in numeric_checks.items()
                    }},
                    "compatibility_shims": {{}},
                    "loaded_forbidden_modules": loaded_forbidden,
                    "unexpected_import_targets": unexpected_import_targets,
                    "observable_values": observable_values,
                    "observable_rows": sorted(observable_values),
                    "recomputed": {{
                        "q1_matrix_element_GeV4": float(report.recomputed_q1_matrix_element_GeV4),
                        "M12_K_NP": {{
                            "real": float(recomputed_m12.real),
                            "imag": float(recomputed_m12.imag),
                        }},
                        "Delta_m_K_NP": float(report.recomputed_delta_m_GeV),
                    }},
                    "paths": {{
                        "wilson": str(wilson_path),
                        "hadronic": str(hadronic_path),
                        "observable": str(observable_path),
                        "provenance": str(provenance_path),
                    }},
                    "expected": {{
                        "observable_rows": EXPECTED_OBSERVABLE_ROWS,
                        "supported_operators": EXPECTED_SUPPORTED_OPERATORS,
                        "unsupported_operators": EXPECTED_UNSUPPORTED_OPERATORS,
                        "filenames": EXPECTED_FILENAMES,
                    }},
                }}
            )
        except Exception as exc:
            tolerance_policy = dataclasses.asdict(verifier.TolerancePolicy())
            _emit_payload(
                {{
                    "ok": False,
                    "import_isolation_ok": import_isolation_ok,
                    "schema_ok": False,
                    "numeric_match_ok": False,
                    "scope_ok": False,
                    "tolerance_policy": tolerance_policy,
                    "tolerances": tolerance_policy,
                    "bundle_ids": {{}},
                    "point_id": None,
                    "issue_codes": ["artifact_schema_error"],
                    "issues": [
                        {{
                            "code": "artifact_schema_error",
                            "message": str(exc),
                            "subject": None,
                        }}
                    ],
                    "numeric_checks": {{}},
                    "numeric_diffs": {{}},
                    "compatibility_shims": {{}},
                    "loaded_forbidden_modules": loaded_forbidden,
                    "unexpected_import_targets": unexpected_import_targets,
                    "observable_values": {{}},
                    "observable_rows": [],
                    "recomputed": {{}},
                    "paths": {{}},
                    "expected": {{
                        "observable_rows": EXPECTED_OBSERVABLE_ROWS,
                        "supported_operators": EXPECTED_SUPPORTED_OPERATORS,
                        "unsupported_operators": EXPECTED_UNSUPPORTED_OPERATORS,
                        "filenames": EXPECTED_FILENAMES,
                    }},
                    "load_error": {{
                        "type": type(exc).__name__,
                        "message": str(exc),
                    }},
                }}
            )
        """
    )
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=False,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "independent verifier subprocess failed:\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    return json.loads(completed.stdout)


def _summarize_static_contracts(package_dir: Path) -> tuple[list[str], dict[str, Any]]:
    py_files = _py_files_under(package_dir)
    failures: list[str] = []

    for path in py_files:
        tree = _read_ast(path)
        relpath = path.relative_to(REPO_ROOT)

        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    module_name = alias.name
                    if module_name in FORBIDDEN_MODULES or any(
                        module_name.startswith(f"{forbidden}.")
                        for forbidden in FORBIDDEN_MODULES
                    ):
                        failures.append(f"{relpath}: forbidden import `{module_name}`")
                    if module_name == "quarkConstraints":
                        failures.append(
                            f"{relpath}: absolute import `quarkConstraints` risks "
                            "repo-v1 leakage"
                        )
                    if module_name.startswith("quarkConstraints.") and not module_name.startswith(
                        "quarkConstraints.paper_0710_1869"
                    ):
                        failures.append(
                            f"{relpath}: absolute import `{module_name}` leaves "
                            "the canonical paper package boundary"
                        )

            elif isinstance(node, ast.ImportFrom):
                module_name = node.module or ""
                imported_names = {alias.name for alias in node.names}
                if module_name in FORBIDDEN_MODULES or any(
                    module_name.startswith(f"{forbidden}.")
                    for forbidden in FORBIDDEN_MODULES
                ):
                    failures.append(f"{relpath}: forbidden import-from `{module_name}`")
                if module_name.startswith("quarkConstraints.") and not module_name.startswith(
                    "quarkConstraints.paper_0710_1869"
                ):
                    failures.append(
                        f"{relpath}: absolute import-from `{module_name}` leaves "
                        "the canonical paper package boundary"
                    )
                if module_name == "quarkConstraints":
                    bad_names = sorted(imported_names & FORBIDDEN_ROOT_EXPORTS)
                    if bad_names:
                        failures.append(
                            f"{relpath}: forbidden root exports imported from "
                            f"quarkConstraints: {', '.join(bad_names)}"
                        )
                if (
                    path.parent == package_dir
                    and node.level >= 2
                    and module_name in {"couplings", "deltaf2", "model"}
                ):
                    failures.append(
                        f"{relpath}: relative import reaches repo-v1 `{module_name}`"
                    )

            elif isinstance(node, ast.Name) and node.id == "BulkMassMap":
                failures.append(f"{relpath}: forbidden symbol `BulkMassMap`")

    summary = {
        "python_files": [str(path.relative_to(REPO_ROOT)) for path in py_files],
        "forbidden_references": sorted(set(failures)),
    }
    return sorted(set(failures)), summary


def _import_paper_modules() -> tuple[dict[str, Any], dict[str, str], list[str]]:
    imported_modules: dict[str, Any] = {}
    import_summary: dict[str, str] = {}
    failures: list[str] = []

    for module_name in REQUIRED_MODULE_IMPORT_ORDER:
        module_path = PAPER_PACKAGE_DIR.joinpath(*module_name.split(".")).with_suffix(".py")
        if not module_path.exists():
            failures.append(
                f"{module_path.relative_to(REPO_ROOT)}: required paper module is missing"
            )
            continue
        try:
            imported_modules[module_name] = importlib.import_module(
                f"quarkConstraints.paper_0710_1869.{module_name}"
            )
            import_summary[module_name] = "ok"
        except Exception as exc:
            import_summary[module_name] = f"{type(exc).__name__}: {exc}"
            failures.append(
                f"{module_path.relative_to(REPO_ROOT)}: import failed: "
                f"{type(exc).__name__}: {exc}"
            )

    for module_name in OPTIONAL_MODULE_IMPORT_ORDER:
        module_path = PAPER_PACKAGE_DIR.joinpath(*module_name.split(".")).with_suffix(".py")
        if not module_path.exists():
            import_summary[module_name] = "missing"
            continue
        try:
            imported_modules[module_name] = importlib.import_module(
                f"quarkConstraints.paper_0710_1869.{module_name}"
            )
            import_summary[module_name] = "ok"
        except Exception as exc:
            import_summary[module_name] = f"{type(exc).__name__}: {exc}"
            failures.append(
                f"{module_path.relative_to(REPO_ROOT)}: import failed: "
                f"{type(exc).__name__}: {exc}"
            )

    return imported_modules, import_summary, failures


def _build_coupling_summary(modules: Mapping[str, Any]) -> tuple[list[str], dict[str, Any]]:
    module = modules.get("couplings")
    if module is None:
        return [], {"status": "missing"}

    summary_fn = None
    for name in (
        "coupling_contract_summary",
        "build_coupling_contract_summary",
        "paper_0710_1869_coupling_contract_summary",
    ):
        candidate = getattr(module, name, None)
        if callable(candidate):
            summary_fn = candidate
            break

    if summary_fn is None:
        return (
            [
                "couplings module is present but exposes no contract summary callable; "
                "expected one of: coupling_contract_summary / build_coupling_contract_summary / "
                "paper_0710_1869_coupling_contract_summary"
            ],
            {"status": "missing_contract_summary_api"},
        )

    payload = _canonicalize(summary_fn())
    second = _canonicalize(summary_fn())
    deterministic = json.dumps(payload, sort_keys=True) == json.dumps(second, sort_keys=True)
    failures = []
    if not deterministic:
        failures.append("coupling contract summary is not deterministic within one process")
    summary: dict[str, Any] = {"status": "ok", "deterministic": deterministic, "payload": payload}

    kkgluon_module = modules.get("kkgluon")
    if kkgluon_module is not None:
        kk_summary_fn = getattr(
            kkgluon_module,
            "default_paper_0710_1869_kk_gluon_benchmark_summary",
            None,
        )
        if not callable(kk_summary_fn):
            failures.append(
                "kkgluon module is present but exposes no default benchmark summary callable"
            )
        else:
            kk_payload = _canonicalize(kk_summary_fn())
            kk_second = _canonicalize(kk_summary_fn())
            kk_deterministic = json.dumps(kk_payload, sort_keys=True) == json.dumps(
                kk_second,
                sort_keys=True,
            )
            if not kk_deterministic:
                failures.append("kkgluon benchmark summary is not deterministic within one process")
            summary["kkgluon"] = {
                "deterministic": kk_deterministic,
                "payload": kk_payload,
            }
    return failures, summary


def _matching_payload_from_value(value: Any) -> dict[str, Any]:
    if hasattr(value, "summary") and callable(value.summary):
        try:
            return _matching_payload_from_value(value.summary())
        except TypeError:
            pass
    if hasattr(value, "as_dict") and callable(value.as_dict):
        payload = _canonicalize(value.as_dict())
    else:
        payload = _canonicalize(value)
    if not isinstance(payload, Mapping):
        raise AssertionError("default matching export does not canonicalize to a mapping")
    return dict(payload)


def _payload_from_value(value: Any) -> dict[str, Any]:
    return _matching_payload_from_value(value)


def _get_callable(module: Any, names: Sequence[str]) -> Any:
    for name in names:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    return None


def _get_class(module: Any, names: Sequence[str]) -> type[Any] | None:
    for name in names:
        candidate = getattr(module, name, None)
        if isinstance(candidate, type):
            return candidate
    return None


def _default_matching_export_callable(module: Any) -> Any:
    for name in DEFAULT_MATCHING_SUMMARY_NAMES + DEFAULT_MATCHING_EXPORT_NAMES:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    return None


def _default_matching_object_callable(module: Any) -> Any:
    for name in DEFAULT_MATCHING_OBJECT_NAMES:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    return None


def _default_matching_payload(module: Any) -> dict[str, Any]:
    summary_fn = _default_matching_export_callable(module)
    if summary_fn is None:
        raise AssertionError(
            "matching module is present but exposes no default benchmark export; expected one of "
            + ", ".join(DEFAULT_MATCHING_SUMMARY_NAMES + DEFAULT_MATCHING_EXPORT_NAMES)
        )
    return _matching_payload_from_value(summary_fn())


def _nested_value(mapping: Mapping[str, Any], key_paths: Sequence[Sequence[str]]) -> Any | None:
    for path in key_paths:
        current: Any = mapping
        for key in path:
            if not isinstance(current, Mapping) or key not in current:
                current = None
                break
            current = current[key]
        if current is not None:
            return current
    return None


def _mapping_contains_key_deep(value: Any, target_key: str) -> bool:
    if isinstance(value, Mapping):
        if target_key in value:
            return True
        return any(_mapping_contains_key_deep(item, target_key) for item in value.values())
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return any(_mapping_contains_key_deep(item, target_key) for item in value)
    return False


def _string_list(value: Any) -> list[str]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [str(item) for item in value]
    return []


def _real_matrix_2x2(value: Any) -> tuple[tuple[float, float], tuple[float, float]] | None:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        return None
    rows = list(value)
    if len(rows) != 2:
        return None
    normalized: list[tuple[float, float]] = []
    for row in rows:
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes, bytearray)):
            return None
        entries = list(row)
        if len(entries) != 2:
            return None
        try:
            normalized.append(tuple(float(entry) for entry in entries))
        except (TypeError, ValueError):
            return None
    return (normalized[0], normalized[1])


def _real_matrix_multiply(
    left: tuple[tuple[float, float], tuple[float, float]],
    right: tuple[tuple[float, float], tuple[float, float]],
) -> tuple[tuple[float, float], tuple[float, float]]:
    return (
        (
            (left[0][0] * right[0][0]) + (left[0][1] * right[1][0]),
            (left[0][0] * right[0][1]) + (left[0][1] * right[1][1]),
        ),
        (
            (left[1][0] * right[0][0]) + (left[1][1] * right[1][0]),
            (left[1][0] * right[0][1]) + (left[1][1] * right[1][1]),
        ),
    )


def _real_matrix_vector_multiply(
    matrix: tuple[tuple[float, float], tuple[float, float]],
    vector: tuple[float, float],
) -> tuple[float, float]:
    return (
        (matrix[0][0] * vector[0]) + (matrix[0][1] * vector[1]),
        (matrix[1][0] * vector[0]) + (matrix[1][1] * vector[1]),
    )


def _complex_as_dict(value: complex) -> dict[str, float]:
    number = complex(value)
    return {
        "real": float(number.real),
        "imag": float(number.imag),
    }


def _complex_from_payload(value: Any) -> complex:
    if not isinstance(value, Mapping):
        raise TypeError("complex payload must be a mapping")
    return complex(float(value["real"]), float(value["imag"]))


def _complex_vector_as_dict_list(
    values: Sequence[complex],
) -> list[dict[str, float]]:
    return [_complex_as_dict(complex(value)) for value in values]


def _complex_matrix_as_dict(
    matrix: Sequence[Sequence[complex]],
) -> list[list[dict[str, float]]]:
    return [
        [_complex_as_dict(complex(entry)) for entry in row]
        for row in matrix
    ]


def _complex_matrix_multiply(
    left: tuple[tuple[complex, complex], tuple[complex, complex]],
    right: tuple[tuple[complex, complex], tuple[complex, complex]],
) -> tuple[tuple[complex, complex], tuple[complex, complex]]:
    return (
        (
            (left[0][0] * right[0][0]) + (left[0][1] * right[1][0]),
            (left[0][0] * right[0][1]) + (left[0][1] * right[1][1]),
        ),
        (
            (left[1][0] * right[0][0]) + (left[1][1] * right[1][0]),
            (left[1][0] * right[0][1]) + (left[1][1] * right[1][1]),
        ),
    )


def _complex_matrix_vector_multiply(
    matrix: tuple[tuple[complex, complex], tuple[complex, complex]],
    vector: tuple[complex, complex],
) -> tuple[complex, complex]:
    return (
        (matrix[0][0] * vector[0]) + (matrix[0][1] * vector[1]),
        (matrix[1][0] * vector[0]) + (matrix[1][1] * vector[1]),
    )


def _complex_vector_close(
    left: tuple[complex, complex],
    right: tuple[complex, complex],
    *,
    abs_tol: float = 1.0e-12,
) -> bool:
    return all(abs(lhs - rhs) <= abs_tol for lhs, rhs in zip(left, right, strict=True))


def _complex_matrix_close(
    left: tuple[tuple[complex, complex], tuple[complex, complex]],
    right: tuple[tuple[complex, complex], tuple[complex, complex]],
    *,
    abs_tol: float = 1.0e-12,
) -> bool:
    return all(
        _complex_vector_close(lhs, rhs, abs_tol=abs_tol)
        for lhs, rhs in zip(left, right, strict=True)
    )


def _lr_block_from_paper_basis_matrix(
    matrix: Sequence[Sequence[complex]],
) -> tuple[tuple[complex, complex], tuple[complex, complex]] | None:
    rows = list(matrix)
    if len(rows) != 4:
        return None
    if any(len(row) != 4 for row in rows):
        return None
    return (
        (complex(rows[2][2]), complex(rows[2][3])),
        (complex(rows[3][2]), complex(rows[3][3])),
    )


def _rg_contract_wilson_scheme_id(payload: Mapping[str, Any]) -> Any | None:
    return _nested_value(
        payload,
        (("contract", "renormalization_scheme_id"),),
    )


def _rg_contract_rg_scheme_id(payload: Mapping[str, Any]) -> Any | None:
    return _nested_value(
        payload,
        (("contract", "rg_scheme_id"),),
    )


def _rg_contract_lr_basis_contract_id(payload: Mapping[str, Any]) -> Any | None:
    return _nested_value(
        payload,
        (("contract", "lr_basis_contract_id"),),
    )


def _rg_contract_lr_basis_status_id(payload: Mapping[str, Any]) -> Any | None:
    return _nested_value(
        payload,
        (
            ("contract", "lr_basis_status_id"),
            ("lr_basis_status_id",),
            ("evolution", "lr_basis_status_id"),
        ),
    )


def _output_wilson_scheme_id(payload: Mapping[str, Any]) -> Any | None:
    return _nested_value(
        payload,
        (
            ("scheme_id",),
            ("renormalization_scheme_id",),
            ("matching_summary", "scheme_id"),
            ("matching_summary", "renormalization_scheme_id"),
            ("matching_summary", "tags", "renormalization_scheme_id"),
            ("tags", "renormalization_scheme_id"),
            ("wilsons", "tags", "renormalization_scheme_id"),
        ),
    )


def _output_rg_scheme_id(payload: Mapping[str, Any]) -> Any | None:
    return _nested_value(
        payload,
        (
            ("rg_scheme_id",),
            ("evolution", "rg_scheme_id"),
        ),
    )


def _rg_contract_matching_input_scheme_id(payload: Mapping[str, Any]) -> Any | None:
    return _nested_value(payload, (("contract", "matching_input_scheme_id"),))


def _rg_contract_matching_input_matching_id(payload: Mapping[str, Any]) -> Any | None:
    return _nested_value(payload, (("contract", "matching_input_matching_id"),))


def _rg_contract_matching_input_wilson_schema_id(payload: Mapping[str, Any]) -> Any | None:
    return _nested_value(payload, (("contract", "matching_input_wilson_schema_id"),))


def _matching_summary_scheme_id(payload: Mapping[str, Any]) -> Any | None:
    return _nested_value(
        payload,
        (
            ("matching_summary", "scheme_id"),
            ("matching_summary", "renormalization_scheme_id"),
            ("matching_summary", "tags", "renormalization_scheme_id"),
            ("input_matching_summary", "scheme_id"),
            ("input_matching_summary", "renormalization_scheme_id"),
            ("input_matching_summary", "tags", "renormalization_scheme_id"),
            ("source_renormalization_scheme_id",),
            ("tags", "source_renormalization_scheme_id"),
            ("wilsons", "source_renormalization_scheme_id"),
            ("wilsons", "tags", "source_renormalization_scheme_id"),
        ),
    )


def _matching_summary_matching_id(payload: Mapping[str, Any]) -> Any | None:
    return _nested_value(
        payload,
        (
            ("matching_summary", "matching_id"),
            ("input_matching_summary", "matching_id"),
            ("tags", "matching_input_matching_id"),
            ("wilsons", "tags", "matching_input_matching_id"),
        ),
    )


def _matching_summary_wilson_schema_id(payload: Mapping[str, Any]) -> Any | None:
    return _nested_value(
        payload,
        (
            ("matching_summary", "wilson_schema_id"),
            ("input_matching_summary", "wilson_schema_id"),
            ("source_wilson_schema_id",),
            ("tags", "source_wilson_schema_id"),
            ("wilsons", "source_wilson_schema_id"),
            ("wilsons", "tags", "source_wilson_schema_id"),
        ),
    )


def _default_rg_export_callable(module: Any) -> Any:
    for name in DEFAULT_RG_SUMMARY_NAMES + DEFAULT_RG_EXPORT_NAMES:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    return None


def _invoke_rg_evolution(callable_obj: Any, *, matching_or_wilsons: Any, mu_low_GeV: float) -> Any:
    import inspect

    parameters = inspect.signature(callable_obj).parameters
    kwargs: dict[str, Any] = {}
    wilsons = getattr(matching_or_wilsons, "wilsons", matching_or_wilsons)

    for name in ("wilsons", "deltaf2_wilsons", "wilsons_in", "coefficients"):
        if name in parameters:
            kwargs[name] = wilsons
            break
    for name in ("matching", "deltaf2_matching", "match"):
        if name in parameters:
            kwargs[name] = matching_or_wilsons
            break
    for name in ("mu_low_GeV", "mu_eval_GeV", "evaluation_scale_GeV", "scale_GeV"):
        if name in parameters:
            kwargs[name] = mu_low_GeV
            break

    missing = [
        name
        for name, parameter in parameters.items()
        if parameter.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
        and parameter.default is inspect._empty
        and name not in kwargs
    ]
    if missing:
        raise AssertionError(
            "RG evolution callable requires unsupported arguments for acceptance benchmark: "
            + ", ".join(missing)
        )
    return callable_obj(**kwargs)


def _default_rg_payload(modules: Mapping[str, Any]) -> dict[str, Any]:
    module = modules["eft_deltaf2.rg"]
    summary_fn = _default_rg_export_callable(module)
    if callable(summary_fn):
        return _payload_from_value(summary_fn())

    evolve_fn = None
    for name in PARAMETERIZED_RG_EXPORT_NAMES:
        candidate = getattr(module, name, None)
        if callable(candidate):
            evolve_fn = candidate
            break
    if not callable(evolve_fn):
        raise AssertionError(
            "RG module is present but exposes no default export; expected one of "
            + ", ".join(
                DEFAULT_RG_SUMMARY_NAMES
                + DEFAULT_RG_EXPORT_NAMES
                + PARAMETERIZED_RG_EXPORT_NAMES
            )
        )

    matching_module = modules.get("eft_deltaf2.matching_kkgluon")
    if matching_module is None:
        raise AssertionError("RG module is present but eft_deltaf2.matching_kkgluon is missing")
    matching_fn = _default_matching_object_callable(matching_module)
    if not callable(matching_fn):
        raise AssertionError(
            "matching module is present but exposes no default object export required by RG"
        )

    return _payload_from_value(
        _invoke_rg_evolution(
            evolve_fn,
            matching_or_wilsons=matching_fn(),
            mu_low_GeV=DEFAULT_RG_LOW_SCALE_GEV,
        )
    )


def _default_hadronic_payload(modules: Mapping[str, Any]) -> dict[str, Any]:
    module = modules["eft_deltaf2.hadronic"]
    callable_obj = _get_callable(
        module,
        DEFAULT_HADRONIC_SUMMARY_NAMES + DEFAULT_HADRONIC_EXPORT_NAMES,
    )
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic module is present but exposes no default export; expected one of "
            + ", ".join(DEFAULT_HADRONIC_SUMMARY_NAMES + DEFAULT_HADRONIC_EXPORT_NAMES)
        )
    return _payload_from_value(callable_obj())


def _default_hadronic_builder_payload(modules: Mapping[str, Any]) -> dict[str, Any]:
    module = modules["eft_deltaf2.hadronic"]
    builder = _get_callable(
        module,
        HADRONIC_BUILDER_NAMES,
    )
    if not callable(builder):
        raise AssertionError(
            "hadronic module is present but exposes no default hadronic builder; expected one "
            "of " + ", ".join(HADRONIC_BUILDER_NAMES)
        )
    return _payload_from_value(builder())


def _default_hadronic_value(modules: Mapping[str, Any]) -> Any:
    module = modules["eft_deltaf2.hadronic"]
    callable_obj = _get_callable(
        module,
        DEFAULT_HADRONIC_EXPORT_NAMES + DEFAULT_HADRONIC_SUMMARY_NAMES,
    )
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic module is present but exposes no default kaon export; expected one of "
            + ", ".join(DEFAULT_HADRONIC_EXPORT_NAMES + DEFAULT_HADRONIC_SUMMARY_NAMES)
        )
    return callable_obj()


def _default_kaon_lr_hadronic_payload(modules: Mapping[str, Any]) -> dict[str, Any]:
    module = modules["eft_deltaf2.hadronic"]
    callable_obj = _get_callable(
        module,
        DEFAULT_LR_HADRONIC_EXPORT_NAMES,
    )
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic module is present but exposes no default LR hadronic export; expected one "
            "of " + ", ".join(DEFAULT_LR_HADRONIC_EXPORT_NAMES)
        )
    return _payload_from_value(callable_obj())


def _default_kaon_lr_hadronic_value(modules: Mapping[str, Any]) -> Any:
    module = modules["eft_deltaf2.hadronic"]
    callable_obj = _get_callable(
        module,
        DEFAULT_LR_HADRONIC_VALUE_EXPORT_NAMES,
    )
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic module is present but exposes no default LR hadronic object export; "
            "expected one of " + ", ".join(DEFAULT_LR_HADRONIC_VALUE_EXPORT_NAMES)
        )
    return callable_obj()


def _default_kaon_lr_hadronic_builder_payload(
    modules: Mapping[str, Any],
) -> dict[str, Any]:
    module = modules["eft_deltaf2.hadronic"]
    builder = getattr(module, "build_paper_0710_1869_default_kaon_lr_hadronic_inputs", None)
    if not callable(builder):
        raise AssertionError(
            "hadronic module is present but exposes no default LR hadronic builder; "
            "expected build_paper_0710_1869_default_kaon_lr_hadronic_inputs"
        )
    return _payload_from_value(builder())


def _default_kaon_lr_r_chi_freeze_payload(modules: Mapping[str, Any]) -> dict[str, Any]:
    module = modules["eft_deltaf2.hadronic"]
    callable_obj = _get_callable(
        module,
        KAON_LR_R_CHI_SUMMARY_EXPORT_NAMES + KAON_LR_R_CHI_FREEZE_EXPORT_NAMES,
    )
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic module is present but exposes no frozen LR R_chi export; expected one of "
            + ", ".join(KAON_LR_R_CHI_SUMMARY_EXPORT_NAMES + KAON_LR_R_CHI_FREEZE_EXPORT_NAMES)
        )
    return _payload_from_value(callable_obj())


def _default_kaon_lr_r_chi_freeze_value(modules: Mapping[str, Any]) -> Any:
    module = modules["eft_deltaf2.hadronic"]
    callable_obj = _get_callable(
        module,
        KAON_LR_R_CHI_FREEZE_EXPORT_NAMES + KAON_LR_R_CHI_SUMMARY_EXPORT_NAMES,
    )
    if not callable(callable_obj):
        raise AssertionError(
            "hadronic module is present but exposes no frozen LR R_chi value; expected one of "
            + ", ".join(KAON_LR_R_CHI_FREEZE_EXPORT_NAMES + KAON_LR_R_CHI_SUMMARY_EXPORT_NAMES)
        )
    return callable_obj()


def _custom_lr_hadronic_builder(module: Any) -> Any:
    for name in CUSTOM_LR_HADRONIC_BUILDER_NAMES:
        candidate = getattr(module, name, None)
        if callable(candidate):
            return candidate
    for name in HADRONIC_BUILDER_NAMES:
        candidate = getattr(module, name, None)
        if not callable(candidate):
            continue
        parameters = inspect.signature(candidate).parameters
        required_groups = (
            LR_HADRONIC_B4_PARAM_NAMES,
            LR_HADRONIC_B5_PARAM_NAMES,
            LR_HADRONIC_R_CHI_PARAM_NAMES,
            LR_HADRONIC_B4_SOURCE_PARAM_NAMES,
            LR_HADRONIC_B5_SOURCE_PARAM_NAMES,
            LR_HADRONIC_R_CHI_SOURCE_PARAM_NAMES,
        )
        if all(any(name in parameters for name in group) for group in required_groups):
            return candidate
    raise AssertionError(
        "hadronic module is present but exposes no custom LR builder; expected one of "
        + ", ".join(CUSTOM_LR_HADRONIC_BUILDER_NAMES)
        + " or an LR-extended "
        + ", ".join(HADRONIC_BUILDER_NAMES)
    )


def _build_hadronic_source_ref(
    module: Any,
    *,
    source_id: str,
    citation: str,
    locator_label: str,
    year: int,
    scheme_id: str,
    scale_GeV: float,
    notes: str,
) -> Any:
    source_ref_cls = _get_class(module, HADRONIC_SOURCE_REF_CLASS_NAMES)
    if source_ref_cls is None:
        raise AssertionError(
            "hadronic module is present but exposes no source-ref class; expected one of "
            + ", ".join(HADRONIC_SOURCE_REF_CLASS_NAMES)
        )

    parameters = inspect.signature(source_ref_cls).parameters
    kwargs: dict[str, Any] = {}
    if "source_id" in parameters:
        kwargs["source_id"] = source_id
    if "source_kind" in parameters:
        kwargs["source_kind"] = "review-average"
    if "citation" in parameters:
        kwargs["citation"] = citation
    if "locator_label" in parameters:
        kwargs["locator_label"] = locator_label
    if "year" in parameters:
        kwargs["year"] = year
    if "renormalization_scheme_id" in parameters:
        kwargs["renormalization_scheme_id"] = scheme_id
    if "scale_GeV" in parameters:
        kwargs["scale_GeV"] = scale_GeV
    if "transformation_id" in parameters:
        kwargs["transformation_id"] = "none"
    if "notes" in parameters:
        kwargs["notes"] = notes

    missing = [
        name
        for name, parameter in parameters.items()
        if parameter.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
        and parameter.default is inspect._empty
        and name not in kwargs
    ]
    if missing:
        raise AssertionError(
            "hadronic source-ref constructor requires unsupported arguments for LR-HAD-1 "
            "acceptance benchmarking: "
            + ", ".join(missing)
        )
    return source_ref_cls(**kwargs)


def _invoke_custom_lr_hadronic_builder(
    callable_obj: Any,
    *,
    default_bundle: Any,
    B4_mu_had: float,
    B5_mu_had: float,
    R_chi_mu_had: float,
    B4_source: Any,
    B5_source: Any,
    R_chi_source: Any,
    renormalization_scheme_id: str | None = None,
    mu_had_GeV: float | None = None,
    source_prefix: str = "hadronic.kaon.lr_custom_probe",
) -> Any:
    parameters = inspect.signature(callable_obj).parameters
    kwargs: dict[str, Any] = {}

    def assign_first(candidates: Sequence[str], value: Any) -> None:
        for name in candidates:
            if name in parameters:
                kwargs[name] = value
                return

    assign_first(LR_HADRONIC_B4_PARAM_NAMES, B4_mu_had)
    assign_first(LR_HADRONIC_B5_PARAM_NAMES, B5_mu_had)
    assign_first(LR_HADRONIC_R_CHI_PARAM_NAMES, R_chi_mu_had)
    assign_first(LR_HADRONIC_B4_SOURCE_PARAM_NAMES, B4_source)
    assign_first(LR_HADRONIC_B5_SOURCE_PARAM_NAMES, B5_source)
    assign_first(LR_HADRONIC_R_CHI_SOURCE_PARAM_NAMES, R_chi_source)

    resolved_mu_had_GeV = (
        float(default_bundle.mu_had_GeV) if mu_had_GeV is None else float(mu_had_GeV)
    )
    resolved_renormalization_scheme_id = (
        str(default_bundle.renormalization_scheme_id)
        if renormalization_scheme_id is None
        else str(renormalization_scheme_id)
    )
    if "mu_had_GeV" in parameters:
        kwargs["mu_had_GeV"] = resolved_mu_had_GeV
    if "m_K0_GeV" in parameters:
        kwargs["m_K0_GeV"] = float(default_bundle.m_K0_GeV)
    if "f_K_GeV" in parameters:
        kwargs["f_K_GeV"] = float(default_bundle.f_K_GeV)
    if "renormalization_scheme_id" in parameters:
        kwargs["renormalization_scheme_id"] = resolved_renormalization_scheme_id
    if "scheme_id" in parameters:
        kwargs["scheme_id"] = resolved_renormalization_scheme_id
    if "mass_source" in parameters and hasattr(default_bundle, "mass_source"):
        kwargs["mass_source"] = default_bundle.mass_source
    if "decay_constant_source" in parameters and hasattr(default_bundle, "decay_constant_source"):
        kwargs["decay_constant_source"] = default_bundle.decay_constant_source
    if "bundle_id" in parameters:
        kwargs["bundle_id"] = f"{source_prefix}.v1"
    if "source_id" in parameters:
        kwargs["source_id"] = f"{source_prefix}.sources.v1"
    if "provenance_ids" in parameters:
        provenance_ids = [f"{source_prefix}.sources.v1"]
        for attr_name in ("mass_source", "decay_constant_source", "bag_parameter_source"):
            source = getattr(default_bundle, attr_name, None)
            source_id = getattr(source, "source_id", None)
            if source_id:
                provenance_ids.append(str(source_id))
        provenance_ids.extend(
            (
                f"{source_prefix}.b4.v1",
                f"{source_prefix}.b5.v1",
                f"{source_prefix}.rchi.v1",
            )
        )
        kwargs["provenance_ids"] = tuple(dict.fromkeys(provenance_ids))

    missing = [
        name
        for name, parameter in parameters.items()
        if parameter.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
        and parameter.default is inspect._empty
        and name not in kwargs
    ]
    if missing:
        raise AssertionError(
            "custom LR hadronic builder requires unsupported arguments for acceptance "
            "benchmarking: "
            + ", ".join(missing)
        )
    return callable_obj(**kwargs)


def _extract_lr_matrix_element(
    payload: Mapping[str, Any],
    operator_name: str,
) -> float | None:
    short_name = operator_name.lower().split("_", 1)[0]
    candidate = _nested_value(
        payload,
        (
            (f"{short_name}_matrix_element_GeV4",),
            (f"{short_name}_matrix_element_gev4",),
            (f"{operator_name.lower()}_matrix_element_GeV4",),
            (f"{operator_name.lower()}_matrix_element_gev4",),
            ("matrix_elements_GeV4", operator_name),
            ("matrix_elements", operator_name, "GeV4"),
            ("matrix_elements", operator_name, "value_GeV4"),
            ("matrix_elements", operator_name),
            ("lr_matrix_elements_GeV4", operator_name),
            ("lr_matrix_elements", operator_name, "GeV4"),
            ("lr_matrix_elements", operator_name),
        ),
    )
    if candidate is None:
        return None
    try:
        return float(candidate)
    except (TypeError, ValueError):
        return None


def _source_metadata_text(payload: Mapping[str, Any], source_key: str) -> str:
    source_payload = _nested_value(payload, ((source_key,),))
    if not isinstance(source_payload, Mapping):
        return ""
    return " ".join(
        str(source_payload.get(field_name, ""))
        for field_name in ("citation", "locator_label", "notes")
    )


def _normalized_metadata_text(value: object) -> str:
    lowered = str(value).strip().lower().replace("-", " ").replace("_", " ")
    return " ".join(lowered.split())


def _is_current_custom_lr_input_policy(policy_id: object) -> bool:
    lowered = str(policy_id).strip().lower()
    return bool(lowered) and lowered == EXPECTED_CUSTOM_LR_HADRONIC_INPUT_POLICY_ID.lower()


def _has_current_custom_lr_note_core(notes: object) -> bool:
    lowered = _normalized_metadata_text(notes)
    return bool(lowered) and "defaults not frozen" not in lowered


def _is_current_custom_lr_contract_notes(notes: object) -> bool:
    lowered = _normalized_metadata_text(notes)
    return (
        _has_current_custom_lr_note_core(notes)
        and "frozen default lr hadronic bundle remains separate" in lowered
        and "not auto consumed on this custom surface" in lowered
    )


def _is_current_custom_lr_bundle_notes(notes: object) -> bool:
    lowered = _normalized_metadata_text(notes)
    return (
        _has_current_custom_lr_note_core(notes)
        and "frozen default lr hadronic bundle remains separate" in lowered
        and (
            "lr only and custom combined observable surfaces still "
            "require explicit custom lr inputs"
            in lowered
        )
    )


def _is_custom_input_only_guard_error(error: str | None) -> bool:
    if not error:
        return False
    lowered = error.lower()
    return "input_provenance_mode_id" in lowered or "custom" in lowered


def _replace_rg_lr_wilsons(
    rg_value: Any,
    *,
    q4_lr: complex,
    q5_lr: complex,
    q1_vll: complex | None = None,
    q1_vrr: complex | None = None,
    renormalization_scheme_id: str | None = None,
    matching_scale_GeV: float | None = None,
) -> Any:
    if not dataclasses.is_dataclass(rg_value) or isinstance(rg_value, type):
        raise AssertionError("RG result is not a dataclass and cannot be patched for LR probing")
    wilsons = getattr(rg_value, "wilsons", None)
    if not dataclasses.is_dataclass(wilsons) or isinstance(wilsons, type):
        raise AssertionError("RG result exposes no dataclass wilsons payload for LR probing")
    resolved_contract = wilsons.contract
    if renormalization_scheme_id is not None:
        resolved_contract = dataclasses.replace(
            resolved_contract,
            renormalization_scheme_id=renormalization_scheme_id,
        )
    return dataclasses.replace(
        rg_value,
        wilsons=dataclasses.replace(
            wilsons,
            contract=resolved_contract,
            q1_vll=wilsons.q1_vll if q1_vll is None else q1_vll,
            q1_vrr=wilsons.q1_vrr if q1_vrr is None else q1_vrr,
            q4_lr=q4_lr,
            q5_lr=q5_lr,
            matching_scale_GeV=(
                wilsons.matching_scale_GeV
                if matching_scale_GeV is None
                else matching_scale_GeV
            ),
        ),
    )


def _default_rg_value(modules: Mapping[str, Any]) -> Any:
    module = modules["eft_deltaf2.rg"]
    evolve_fn = _get_callable(module, PARAMETERIZED_RG_EXPORT_NAMES)
    if not callable(evolve_fn):
        raise AssertionError(
            "RG module is present but exposes no callable evolution API; expected one of "
            + ", ".join(PARAMETERIZED_RG_EXPORT_NAMES)
        )

    matching_module = modules.get("eft_deltaf2.matching_kkgluon")
    if matching_module is None:
        raise AssertionError("RG module is present but eft_deltaf2.matching_kkgluon is missing")
    matching_fn = _default_matching_object_callable(matching_module)
    if not callable(matching_fn):
        raise AssertionError(
            "matching module is present but exposes no default object export required by RG"
        )

    return _invoke_rg_evolution(
        evolve_fn,
        matching_or_wilsons=matching_fn(),
        mu_low_GeV=DEFAULT_RG_LOW_SCALE_GEV,
    )


def _invoke_observable_evaluation(callable_obj: Any, *, rg_value: Any, hadronic_value: Any) -> Any:
    import inspect

    parameters = inspect.signature(callable_obj).parameters
    kwargs: dict[str, Any] = {}

    for name in RG_RESULT_PARAM_NAMES:
        if name in parameters:
            kwargs[name] = rg_value
            break
    else:
        for name in WILSON_PARAM_NAMES:
            if name in parameters:
                kwargs[name] = getattr(rg_value, "wilsons", rg_value)
                break

    for name in HADRONIC_PARAM_NAMES:
        if name in parameters:
            kwargs[name] = hadronic_value
            break

    missing = [
        name
        for name, parameter in parameters.items()
        if parameter.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
        and parameter.default is inspect._empty
        and name not in kwargs
    ]
    if missing:
        raise AssertionError(
            "observable evaluator requires unsupported arguments for acceptance benchmark: "
            + ", ".join(missing)
        )
    return callable_obj(**kwargs)


def _invoke_custom_total_observable_evaluation(
    callable_obj: Any,
    *,
    rg_value: Any,
    q1_hadronic_bundle: Any,
    lr_hadronic_inputs: Any,
) -> Any:
    parameters = inspect.signature(callable_obj).parameters
    kwargs: dict[str, Any] = {}

    for name in RG_RESULT_PARAM_NAMES:
        if name in parameters:
            kwargs[name] = rg_value
            break
    else:
        for name in WILSON_PARAM_NAMES:
            if name in parameters:
                kwargs[name] = getattr(rg_value, "wilsons", rg_value)
                break

    for name in ("q1_hadronic_bundle", "q1_bundle", "hadronic_bundle"):
        if name in parameters:
            kwargs[name] = q1_hadronic_bundle
            break
    for name in ("lr_hadronic_inputs", "lr_inputs", "hadronic_inputs"):
        if name in parameters:
            kwargs[name] = lr_hadronic_inputs
            break

    missing = [
        name
        for name, parameter in parameters.items()
        if parameter.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
        and parameter.default is inspect._empty
        and name not in kwargs
    ]
    if missing:
        raise AssertionError(
            "custom combined observable evaluator requires unsupported arguments for "
            "acceptance benchmark: " + ", ".join(missing)
        )
    return callable_obj(**kwargs)


def _build_custom_q1_hadronic_bundle(
    modules: Mapping[str, Any],
    *,
    mu_had_GeV: float = DEFAULT_RG_LOW_SCALE_GEV,
) -> Any:
    hadronic_module = modules["eft_deltaf2.hadronic"]
    builder = _get_callable(hadronic_module, HADRONIC_BUILDER_NAMES)
    if not callable(builder):
        raise AssertionError(
            "hadronic module is present but exposes no Q1 hadronic builder; expected one of "
            + ", ".join(HADRONIC_BUILDER_NAMES)
        )
    default_bundle = _default_hadronic_value(modules)
    custom_bag_source = dataclasses.replace(
        default_bundle.bag_parameter_source,
        source_id=CUSTOM_Q1_HADRONIC_PROBE_BAG_SOURCE_ID,
        citation=CUSTOM_Q1_HADRONIC_PROBE_BAG_SOURCE_CITATION,
        locator_label=CUSTOM_Q1_HADRONIC_PROBE_BAG_SOURCE_LOCATOR,
        year=CUSTOM_Q1_HADRONIC_PROBE_YEAR,
        notes=CUSTOM_Q1_HADRONIC_PROBE_NOTES,
    )
    return builder(
        mu_had_GeV=float(mu_had_GeV),
        m_K0_GeV=float(default_bundle.m_K0_GeV),
        f_K_GeV=float(default_bundle.f_K_GeV),
        hat_B_K_rgi_source_value=float(default_bundle.hat_B_K_rgi_source_value),
        bag_parameter_source=custom_bag_source,
    )


def _get_custom_b_q1_observable_evaluator(module: Any, system_id: str) -> Any:
    config = B_CUSTOM_Q1_PROBE_CONFIG[system_id]
    callable_obj = _get_callable(module, config["observable_names"])
    if not callable(callable_obj):
        raise AssertionError(
            f"observable module exposes no callable {system_id} custom Q1 evaluator; "
            "expected one of "
            + ", ".join(config["observable_names"])
        )
    return callable_obj


def _default_b_matching_object(modules: Mapping[str, Any], system_id: str) -> Any:
    module = modules["eft_deltaf2.matching_kkgluon"]
    config = B_CUSTOM_Q1_PROBE_CONFIG[system_id]
    callable_obj = _get_callable(module, config["matching_names"])
    if not callable(callable_obj):
        raise AssertionError(
            f"matching module exposes no callable {system_id} default object; expected one of "
            + ", ".join(config["matching_names"])
        )
    return callable_obj()


def _build_custom_b_rg_result(
    modules: Mapping[str, Any],
    system_id: str,
    *,
    mu_low_GeV: float = B_CUSTOM_Q1_PROBE_MU_HAD_GEV,
) -> Any:
    module = modules["eft_deltaf2.rg"]
    callable_obj = _get_callable(module, PARAMETERIZED_RG_EXPORT_NAMES)
    if not callable(callable_obj):
        raise AssertionError(
            "RG module is present but exposes no evolution API; expected one of "
            + ", ".join(PARAMETERIZED_RG_EXPORT_NAMES)
        )
    rg_result = _invoke_rg_evolution(
        callable_obj,
        matching_or_wilsons=_default_b_matching_object(modules, system_id),
        mu_low_GeV=float(mu_low_GeV),
    )
    config = B_CUSTOM_Q1_PROBE_CONFIG[system_id]
    wilson_updates: dict[str, complex] = {}
    if "probe_q1_vll" in config:
        wilson_updates["q1_vll"] = complex(config["probe_q1_vll"])
    if "probe_q1_vrr" in config:
        wilson_updates["q1_vrr"] = complex(config["probe_q1_vrr"])
    if not wilson_updates:
        return rg_result
    updated_wilsons = dataclasses.replace(rg_result.wilsons, **wilson_updates)
    return dataclasses.replace(rg_result, wilsons=updated_wilsons)


def _build_custom_b_hadronic_bundle(
    modules: Mapping[str, Any],
    system_id: str,
    *,
    renormalization_scheme_id: str,
    mu_had_GeV: float,
) -> Any:
    hadronic_module = modules["eft_deltaf2.hadronic"]
    config = B_CUSTOM_Q1_PROBE_CONFIG[system_id]
    builder = _get_callable(hadronic_module, config["hadronic_builder_names"])
    if not callable(builder):
        raise AssertionError(
            f"hadronic module exposes no callable {system_id} custom bundle builder; "
            "expected one of "
            + ", ".join(config["hadronic_builder_names"])
        )
    source_prefix = str(config["source_prefix"])
    mass_source = _build_hadronic_source_ref(
        hadronic_module,
        source_id=f"{source_prefix}.mass.v1",
        citation=f"Custom {system_id} meson mass input",
        locator_label=f"{system_id} mass external input",
        year=2026,
        scheme_id=renormalization_scheme_id,
        scale_GeV=mu_had_GeV,
        notes=f"Custom {system_id} mass source for Q1 observable acceptance.",
    )
    decay_constant_source = _build_hadronic_source_ref(
        hadronic_module,
        source_id=f"{source_prefix}.decay_constant.v1",
        citation=f"Custom {system_id} decay constant input",
        locator_label=f"{system_id} decay constant external input",
        year=2026,
        scheme_id=renormalization_scheme_id,
        scale_GeV=mu_had_GeV,
        notes=f"Custom {system_id} decay-constant source for Q1 observable acceptance.",
    )
    bag_parameter_source = _build_hadronic_source_ref(
        hadronic_module,
        source_id=f"{source_prefix}.bag_parameter.v1",
        citation=f"Custom {system_id} bag-parameter input",
        locator_label=f"{system_id} bag parameter external input",
        year=2026,
        scheme_id=renormalization_scheme_id,
        scale_GeV=mu_had_GeV,
        notes=f"Custom {system_id} bag-parameter source for Q1 observable acceptance.",
    )
    parameters = inspect.signature(builder).parameters
    kwargs: dict[str, Any] = {
        "mu_had_GeV": float(mu_had_GeV),
        "renormalization_scheme_id": str(renormalization_scheme_id),
        str(config["mass_param"]): float(config["meson_mass_GeV"]),
        str(config["decay_param"]): float(config["meson_decay_constant_GeV"]),
        str(config["bag_param"]): float(config["bag_parameter_mu_had"]),
        "mass_source": mass_source,
        "decay_constant_source": decay_constant_source,
        "bag_parameter_source": bag_parameter_source,
    }
    if "bundle_id" in parameters:
        kwargs["bundle_id"] = f"{source_prefix}.bundle.v1"
    if "source_id" in parameters:
        kwargs["source_id"] = f"{source_prefix}.sources.v1"
    if "provenance_ids" in parameters:
        kwargs["provenance_ids"] = (
            f"{source_prefix}.sources.v1",
            f"{source_prefix}.mass.v1",
            f"{source_prefix}.decay_constant.v1",
            f"{source_prefix}.bag_parameter.v1",
        )
    return builder(**kwargs)


def _default_observable_payload(modules: Mapping[str, Any]) -> dict[str, Any]:
    module = modules["eft_deltaf2.observables"]
    for name in DEFAULT_OBSERVABLE_SUMMARY_NAMES + DEFAULT_OBSERVABLE_EXPORT_NAMES:
        candidate = getattr(module, name, None)
        if not callable(candidate):
            continue
        try:
            return _payload_from_value(candidate())
        except TypeError:
            continue

    evaluate_fn = _get_callable(module, PARAMETERIZED_OBSERVABLE_EXPORT_NAMES)
    if not callable(evaluate_fn):
        raise AssertionError(
            "observable module is present but exposes no default export; expected one of "
            + ", ".join(
                DEFAULT_OBSERVABLE_SUMMARY_NAMES
                + DEFAULT_OBSERVABLE_EXPORT_NAMES
                + PARAMETERIZED_OBSERVABLE_EXPORT_NAMES
            )
        )

    if "eft_deltaf2.hadronic" not in modules:
        raise AssertionError("observable module is present but eft_deltaf2.hadronic is missing")
    if "eft_deltaf2.rg" not in modules:
        raise AssertionError("observable module is present but eft_deltaf2.rg is missing")

    return _payload_from_value(
        _invoke_observable_evaluation(
            evaluate_fn,
            rg_value=_default_rg_value(modules),
            hadronic_value=_default_hadronic_value(modules),
        )
    )


def _positive_finite(value: Any) -> bool:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return False
    return numeric > 0.0 and math.isfinite(numeric)


def _is_kaon_system(value: Any) -> bool:
    if value is None:
        return False
    lowered = str(value).strip().lower()
    return lowered in {"k", "k0", "kaon"} or "kaon" in lowered


def _is_np_only_interpretation(value: Any) -> bool:
    if value is None:
        return False
    lowered = str(value).strip().lower().replace("_", " ").replace("-", " ")
    return "np only" in lowered or ("np" in lowered and "sm" not in lowered)


def _float_from_value(value: Any) -> float:
    if isinstance(value, bool):
        raise TypeError("boolean values are not observable numbers")
    if isinstance(value, (int, float)):
        numeric = float(value)
        if not math.isfinite(numeric):
            raise ValueError("observable value must be finite")
        return numeric
    if isinstance(value, Mapping):
        for key in ("value", "central_value", "observable_value"):
            if key in value:
                return _float_from_value(value[key])
    raise TypeError(f"unsupported observable value encoding: {value!r}")


def _extract_named_observables(payload: Mapping[str, Any]) -> dict[str, float]:
    container = _nested_value(
        payload,
        (
            ("observables",),
            ("results",),
            ("observable_values",),
            ("quantities",),
        ),
    )
    rows: dict[str, float] = {}
    if isinstance(container, Mapping):
        unexpected = sorted(
            str(key) for key in container if str(key) not in SUPPORTED_OBSERVABLE_CONTAINER_KEYS
        )
        if unexpected:
            raise AssertionError(
                "observable payload exposes unsupported named outputs: "
                + ", ".join(unexpected)
            )
        if "M12_K_NP" not in container or "Delta_m_K_NP" not in container:
            raise AssertionError(
                "observable payload must expose both M12_K_NP and Delta_m_K_NP"
            )
        m12_value = container["M12_K_NP"]
        if not isinstance(m12_value, Mapping) or not {"re", "im"} <= set(m12_value):
            raise AssertionError("M12_K_NP output must expose both re and im components")
        rows["M12_K_NP.re"] = _float_from_value(m12_value["re"])
        rows["M12_K_NP.im"] = _float_from_value(m12_value["im"])
        rows["Delta_m_K_NP"] = _float_from_value(container["Delta_m_K_NP"])
        return rows
    if isinstance(container, Sequence) and not isinstance(container, (str, bytes, bytearray)):
        unexpected_names: list[str] = []
        for index, item in enumerate(container):
            if not isinstance(item, Mapping):
                raise AssertionError("observable rows must be mappings")
            name = item.get("name") or item.get("observable_id") or item.get("quantity")
            if not name:
                name = f"item_{index}"
            if name not in SUPPORTED_OBSERVABLE_CONTAINER_KEYS:
                unexpected_names.append(str(name))
                continue
            if name == "M12_K_NP":
                if not {"re", "im"} <= set(item):
                    raise AssertionError("M12_K_NP output must expose both re and im components")
                rows["M12_K_NP.re"] = _float_from_value(item["re"])
                rows["M12_K_NP.im"] = _float_from_value(item["im"])
            else:
                rows[str(name)] = _float_from_value(item)
        if unexpected_names:
            raise AssertionError(
                "observable payload exposes unsupported named outputs: "
                + ", ".join(sorted(unexpected_names))
            )
        return rows

    for key, row_name in SUPPORTED_OBSERVABLE_FLAT_FIELDS.items():
        if key in payload:
            rows[row_name] = _float_from_value(payload[key])
    if rows:
        unexpected_top_level: list[str] = []
        for key, value in payload.items():
            lowered = str(key).lower()
            if key in SUPPORTED_OBSERVABLE_FLAT_FIELDS or key == "M12_K_NP_GeV":
                continue
            if "epsilon" in lowered:
                unexpected_top_level.append(str(key))
                continue
            if key.startswith("re_M12") or key.startswith("im_M12") or key.startswith("delta_m"):
                try:
                    _float_from_value(value)
                except (TypeError, ValueError):
                    continue
                unexpected_top_level.append(str(key))
        if unexpected_top_level:
            raise AssertionError(
                "observable payload exposes unsupported top-level outputs: "
                + ", ".join(sorted(unexpected_top_level))
            )
        if set(rows) != SUPPORTED_OBSERVABLE_ROWS:
            raise AssertionError(
                "observable payload must expose exactly the PR5a kaon NP-only outputs"
            )
        return rows

    for key in payload:
        lowered = str(key).lower()
        if "epsilon" in lowered or "m12" in lowered or "delta_m" in lowered:
            raise AssertionError(
                "observable payload exposes unsupported observable-like fields outside the "
                "frozen PR5a kaon NP-only surface"
            )
    raise AssertionError("observable payload does not expose the PR5a kaon NP-only outputs")


def _first_mapping_value(mapping: Mapping[str, Any], keys: Sequence[str]) -> Any | None:
    for key in keys:
        if key in mapping:
            return mapping[key]
    return None


def _extract_matching_coefficients(payload: Mapping[str, Any]) -> Mapping[str, Any] | Sequence[Any]:
    coefficients = _first_mapping_value(
        payload,
        ("coefficients", "wilson_coefficients", "C_i", "C"),
    )
    if coefficients is None:
        raise AssertionError("matching payload does not expose C_i(mu_match)")
    if not isinstance(coefficients, (Mapping, Sequence)) or isinstance(
        coefficients,
        (str, bytes, bytearray),
    ):
        raise AssertionError("matching C_i(mu_match) export has an unsupported shape")
    return coefficients


def _build_eft_matching_summary(modules: Mapping[str, Any]) -> tuple[list[str], dict[str, Any]]:
    module = modules.get("eft_deltaf2.matching_kkgluon")
    if module is None:
        return [], {"status": "missing"}

    if "eft_deltaf2.operators" not in modules:
        return (
            ["matching module is present but eft_deltaf2.operators is missing"],
            {"status": "missing_operators_module"},
        )

    try:
        payload = _default_matching_payload(module)
    except AssertionError as exc:
        return ([str(exc)], {"status": "missing_default_matching_export"})

    second = _default_matching_payload(module)
    deterministic = json.dumps(payload, sort_keys=True) == json.dumps(second, sort_keys=True)
    failures = []
    if not deterministic:
        failures.append("default matching export is not deterministic within one process")
    cross_process_payload = _cross_process_default_matching_payload()
    cross_process_deterministic = json.dumps(payload, sort_keys=True) == json.dumps(
        cross_process_payload,
        sort_keys=True,
    )
    if not cross_process_deterministic:
        failures.append("default matching export is not deterministic across processes")
    if not isinstance(payload, Mapping):
        failures.append("default matching export does not canonicalize to a mapping")
        return failures, {"status": "invalid_payload", "deterministic": deterministic}

    basis_id = _first_mapping_value(payload, ("basis_id", "operator_basis_id"))
    scheme_id = _first_mapping_value(payload, ("scheme_id", "renormalization_scheme_id"))
    mu_match_GeV = _first_mapping_value(payload, ("mu_match_GeV", "matching_scale_GeV"))
    propagator_mass_GeV = _first_mapping_value(payload, ("propagator_mass_GeV",))
    if not basis_id:
        failures.append("default matching export is missing a basis tag")
    if not scheme_id:
        failures.append("default matching export is missing a renormalization scheme tag")
    if mu_match_GeV is None:
        failures.append("default matching export is missing a mu_match tag")
    if propagator_mass_GeV is None:
        failures.append("default matching export is missing a propagator-mass tag")

    try:
        coefficients = _extract_matching_coefficients(payload)
    except AssertionError as exc:
        failures.append(str(exc))
        coefficients_summary: dict[str, Any] = {"status": "missing"}
    else:
        coefficients_summary = {
            "status": "ok",
            "count": len(coefficients),
            "keys": (
                sorted(str(key) for key in coefficients.keys())
                if isinstance(coefficients, Mapping)
                else [str(index) for index, _ in enumerate(coefficients)]
            ),
        }

    return failures, {
        "status": "ok",
        "deterministic": deterministic,
        "cross_process_deterministic": cross_process_deterministic,
        "basis_id": basis_id,
        "scheme_id": scheme_id,
        "mu_match_GeV": float(mu_match_GeV) if mu_match_GeV is not None else None,
        "propagator_mass_GeV": (
            float(propagator_mass_GeV) if propagator_mass_GeV is not None else None
        ),
        "default_value_collision": (
            float(mu_match_GeV) == float(propagator_mass_GeV)
            if mu_match_GeV is not None and propagator_mass_GeV is not None
            else None
        ),
        "coefficients": coefficients_summary,
        "payload": payload,
    }


def _build_eft_rg_lo_summary(modules: Mapping[str, Any]) -> tuple[list[str], dict[str, Any]]:
    module = modules.get("eft_deltaf2.rg")
    if module is None:
        return [], {"status": "missing"}

    try:
        payload = _default_rg_payload(modules)
    except AssertionError as exc:
        return ([str(exc)], {"status": "missing_default_rg_export"})

    second = _default_rg_payload(modules)
    deterministic = json.dumps(payload, sort_keys=True) == json.dumps(second, sort_keys=True)
    failures = []
    if not deterministic:
        failures.append("default LO RG export is not deterministic within one process")
    cross_process_payload = _cross_process_default_rg_payload()
    cross_process_deterministic = json.dumps(payload, sort_keys=True) == json.dumps(
        cross_process_payload,
        sort_keys=True,
    )
    if not cross_process_deterministic:
        failures.append("default LO RG export is not deterministic across processes")

    basis_id = _nested_value(
        payload,
        (
            ("basis_id",),
            ("operator_basis_id",),
            ("contract", "operator_basis_id"),
            ("tags", "operator_basis_id"),
            ("wilsons", "tags", "operator_basis_id"),
        ),
    )
    scheme_id = _output_wilson_scheme_id(payload)
    contract_wilson_scheme_id = _rg_contract_wilson_scheme_id(payload)
    rg_scheme_id = _output_rg_scheme_id(payload)
    contract_rg_scheme_id = _rg_contract_rg_scheme_id(payload)
    contract_lr_basis_contract_id = _rg_contract_lr_basis_contract_id(payload)
    contract_lr_basis_status_id = _rg_contract_lr_basis_status_id(payload)
    matching_input_scheme_id = _rg_contract_matching_input_scheme_id(payload)
    matching_input_matching_id = _rg_contract_matching_input_matching_id(payload)
    matching_input_wilson_schema_id = _rg_contract_matching_input_wilson_schema_id(payload)
    matching_summary_scheme_id = _matching_summary_scheme_id(payload)
    matching_summary_matching_id = _matching_summary_matching_id(payload)
    matching_summary_wilson_schema_id = _matching_summary_wilson_schema_id(payload)
    lr_contract_callable = getattr(
        module,
        "default_paper_0710_1869_deltaf2_lr_basis_contract",
        None,
    )
    lr_contract_payload = (
        _payload_from_value(lr_contract_callable()) if callable(lr_contract_callable) else {}
    )
    lr_contract_contract_id = _nested_value(lr_contract_payload, (("contract_id",),))
    lr_contract_status_id = _nested_value(lr_contract_payload, (("status_id",),))
    lr_contract_definition_ids = _nested_value(
        lr_contract_payload,
        (("paper_operator_definition_ids",),),
    )
    paper_lr_operator_order = _string_list(
        _nested_value(lr_contract_payload, (("paper_operator_order",),))
    )
    bmu_lr_operator_order = _string_list(
        _nested_value(lr_contract_payload, (("bmu_lr_operator_order",),))
    )
    paper_to_bmu_operator_map = _real_matrix_2x2(
        _nested_value(lr_contract_payload, (("paper_to_bmu_operator_map_matrix",),))
    )
    bmu_to_paper_operator_map = _real_matrix_2x2(
        _nested_value(lr_contract_payload, (("bmu_to_paper_operator_map_matrix",),))
    )
    paper_to_bmu_wilson_map = _real_matrix_2x2(
        _nested_value(lr_contract_payload, (("paper_to_bmu_wilson_map_matrix",),))
    )
    bmu_to_paper_wilson_map = _real_matrix_2x2(
        _nested_value(lr_contract_payload, (("bmu_to_paper_wilson_map_matrix",),))
    )
    lr_definitions_frozen = _nested_value(lr_contract_payload, (("lr_definitions_frozen",),))
    mapping_matrix_frozen = _nested_value(
        lr_contract_payload,
        (("mapping_matrix_frozen",),),
    )
    lr_running_activated = _nested_value(lr_contract_payload, (("lr_running_activated",),))
    supported_operator_subset_id = _nested_value(
        payload,
        (("contract", "supported_operator_subset_id"),),
    )
    supported_operator_names = _string_list(
        _nested_value(
            payload,
            (
                ("contract", "supported_operator_names"),
                ("evolution", "supported_operator_names"),
            ),
        )
    )
    unsupported_operator_names = _string_list(
        _nested_value(
            payload,
            (
                ("contract", "unsupported_operator_names"),
                ("evolution", "unsupported_operator_names"),
            ),
        )
    )
    lr_basis_map_supported = _nested_value(
        payload,
        (
            ("lr_basis_map_supported",),
            ("evolution", "lr_basis_map_supported"),
        ),
    )
    rg_order = _nested_value(
        payload,
        (
            ("rg_order_id",),
            ("rg_order",),
            ("contract", "rg_order_id"),
            ("evolution", "rg_order_id"),
        ),
    )
    alpha_s_policy_id = _nested_value(
        payload,
        (
            ("alpha_s_policy_id",),
            ("contract", "alpha_s_policy_id"),
            ("evolution", "alpha_s_policy_id"),
        ),
    )
    threshold_policy_id = _nested_value(
        payload,
        (
            ("threshold_policy_id",),
            ("contract", "threshold_policy_id"),
            ("evolution", "threshold_policy_id"),
        ),
    )
    mu_low_GeV = _nested_value(
        payload,
        (
            ("mu_low_GeV",),
            ("evaluation_scale_GeV",),
            ("running_scale_GeV",),
            ("scale_GeV",),
            ("wilsons", "matching_scale_GeV"),
            ("matching_scale_GeV",),
        ),
    )
    segments = _nested_value(
        payload,
        (
            ("segments",),
            ("evolution", "segments"),
        ),
    )
    coefficients = _nested_value(
        payload,
        (
            ("coefficients",),
            ("evolved_coefficients",),
            ("wilsons", "coefficients"),
            ("output_wilsons", "coefficients"),
        ),
    )

    checks = {
        "basis_id_present": bool(basis_id),
        "scheme_id_present": bool(scheme_id),
        "contract_wilson_scheme_present": bool(contract_wilson_scheme_id),
        "contract_rg_scheme_present": bool(contract_rg_scheme_id),
        "lr_contract_id_present": bool(contract_lr_basis_contract_id),
        "lr_contract_status_present": bool(contract_lr_basis_status_id),
        "matching_input_scheme_present": bool(matching_input_scheme_id),
        "matching_input_matching_id_present": bool(matching_input_matching_id),
        "matching_input_wilson_schema_present": bool(matching_input_wilson_schema_id),
        "wilson_scheme_matches_contract": (
            scheme_id == contract_wilson_scheme_id
            if scheme_id is not None and contract_wilson_scheme_id is not None
            else False
        ),
        "rg_scheme_matches_contract": (
            rg_scheme_id == contract_rg_scheme_id
            if rg_scheme_id is not None and contract_rg_scheme_id is not None
            else False
        ),
        "lr_contract_ids_match_freeze_contract": (
            contract_lr_basis_contract_id == lr_contract_contract_id
            if contract_lr_basis_contract_id is not None and lr_contract_contract_id is not None
            else False
        ),
        "lr_contract_status_matches_freeze_contract": (
            contract_lr_basis_status_id == lr_contract_status_id
            if contract_lr_basis_status_id is not None and lr_contract_status_id is not None
            else False
        ),
        "lr_contract_definitions_frozen": lr_definitions_frozen is True,
        "lr_contract_mapping_matrix_is_frozen": mapping_matrix_frozen is True,
        "lr_contract_running_is_active": lr_running_activated is True,
        "public_lr_rg_support_is_active": lr_basis_map_supported is True,
        "supported_operator_surface_is_full_lr_rg": (
            tuple(supported_operator_names)
            == EXPECTED_SUPPORTED_OPERATOR_ORDER + EXPECTED_GUARDED_LR_OPERATOR_ORDER
        ),
        "unsupported_operator_surface_is_empty": tuple(unsupported_operator_names) == tuple(),
        "matching_input_scheme_matches_summary": (
            matching_summary_scheme_id == matching_input_scheme_id
            if matching_summary_scheme_id is not None and matching_input_scheme_id is not None
            else False
        ),
        "matching_input_matching_id_matches_summary": (
            matching_summary_matching_id == matching_input_matching_id
            if matching_summary_matching_id is not None
            and matching_input_matching_id is not None
            else False
        ),
        "matching_input_wilson_schema_matches_summary": (
            matching_summary_wilson_schema_id == matching_input_wilson_schema_id
            if matching_summary_wilson_schema_id is not None
            and matching_input_wilson_schema_id is not None
            else False
        ),
        "rg_order_is_lo": str(rg_order).lower() == "lo",
        "alpha_s_policy_present": bool(alpha_s_policy_id),
        "threshold_policy_present": bool(threshold_policy_id),
        "low_scale_matches_acceptance_target": (
            float(mu_low_GeV) == float(DEFAULT_RG_LOW_SCALE_GEV)
            if mu_low_GeV is not None
            else False
        ),
        "segments_present": isinstance(segments, Sequence)
        and not isinstance(segments, (str, bytes, bytearray))
        and len(segments) > 0,
        "coefficients_present": isinstance(coefficients, Mapping) and len(coefficients) > 0,
    }
    failures.extend(
        f"LO RG acceptance check failed: {name}" for name, ok in checks.items() if not ok
    )

    return failures, {
        "status": "ok",
        "deterministic": deterministic,
        "cross_process_deterministic": cross_process_deterministic,
        "basis_id": basis_id,
        "scheme_id": scheme_id,
        "contract_wilson_scheme_id": contract_wilson_scheme_id,
        "rg_scheme_id": rg_scheme_id,
        "contract_rg_scheme_id": contract_rg_scheme_id,
        "matching_input_scheme_id": matching_input_scheme_id,
        "matching_input_matching_id": matching_input_matching_id,
        "matching_input_wilson_schema_id": matching_input_wilson_schema_id,
        "rg_order": rg_order,
        "alpha_s_policy_id": alpha_s_policy_id,
        "threshold_policy_id": threshold_policy_id,
        "mu_low_GeV": float(mu_low_GeV) if mu_low_GeV is not None else None,
        "lr_running_state": {
            "contract_id": contract_lr_basis_contract_id,
            "status_id": contract_lr_basis_status_id,
            "supported_operator_subset_id": supported_operator_subset_id,
            "supported_operator_names": supported_operator_names,
            "unsupported_operator_names": unsupported_operator_names,
            "lr_running_activated": (
                bool(lr_running_activated)
                if isinstance(lr_running_activated, bool)
                else None
            ),
            "lr_basis_map_supported": (
                bool(lr_basis_map_supported)
                if isinstance(lr_basis_map_supported, bool)
                else None
            ),
        },
        "lr_contract_freeze": {
            "contract_id": lr_contract_contract_id,
            "status_id": lr_contract_status_id,
            "contract_present": bool(lr_contract_contract_id),
            "status_present": bool(lr_contract_status_id),
            "definitions_frozen": lr_definitions_frozen is True,
            "mapping_matrix_frozen": mapping_matrix_frozen,
            "running_active": lr_running_activated is True,
            "paper_operator_order": paper_lr_operator_order,
            "bmu_lr_operator_order": bmu_lr_operator_order,
            "paper_to_bmu_operator_map_matrix": (
                [list(row) for row in paper_to_bmu_operator_map]
                if paper_to_bmu_operator_map is not None
                else None
            ),
            "bmu_to_paper_operator_map_matrix": (
                [list(row) for row in bmu_to_paper_operator_map]
                if bmu_to_paper_operator_map is not None
                else None
            ),
            "paper_to_bmu_wilson_map_matrix": (
                [list(row) for row in paper_to_bmu_wilson_map]
                if paper_to_bmu_wilson_map is not None
                else None
            ),
            "bmu_to_paper_wilson_map_matrix": (
                [list(row) for row in bmu_to_paper_wilson_map]
                if bmu_to_paper_wilson_map is not None
                else None
            ),
            "operator_map_round_trip_is_identity": (
                paper_to_bmu_operator_map is not None
                and bmu_to_paper_operator_map is not None
                and _real_matrix_multiply(bmu_to_paper_operator_map, paper_to_bmu_operator_map)
                == ((1.0, 0.0), (0.0, 1.0))
            ),
            "wilson_map_round_trip_is_identity": (
                paper_to_bmu_wilson_map is not None
                and bmu_to_paper_wilson_map is not None
                and _real_matrix_multiply(bmu_to_paper_wilson_map, paper_to_bmu_wilson_map)
                == ((1.0, 0.0), (0.0, 1.0))
            ),
            "operator_basis_vectors_map_as_frozen": (
                paper_to_bmu_operator_map is not None
                and _real_matrix_vector_multiply(paper_to_bmu_operator_map, (1.0, 0.0))
                == (0.0, 1.0)
                and _real_matrix_vector_multiply(paper_to_bmu_operator_map, (0.0, 1.0))
                == (2.0, 0.0)
            ),
            "wilson_basis_vectors_map_as_frozen": (
                paper_to_bmu_wilson_map is not None
                and _real_matrix_vector_multiply(paper_to_bmu_wilson_map, (1.0, 0.0))
                == (0.0, 1.0)
                and _real_matrix_vector_multiply(paper_to_bmu_wilson_map, (0.0, 1.0))
                == (0.5, 0.0)
            ),
            "contract_matches_rg_export": (
                lr_contract_contract_id == contract_lr_basis_contract_id
                if lr_contract_contract_id is not None and contract_lr_basis_contract_id is not None
                else False
            ),
            "status_matches_rg_export": (
                lr_contract_status_id == contract_lr_basis_status_id
                if lr_contract_status_id is not None and contract_lr_basis_status_id is not None
                else False
            ),
            "paper_operator_definition_ids": (
                list(lr_contract_definition_ids)
                if isinstance(lr_contract_definition_ids, Sequence)
                and not isinstance(lr_contract_definition_ids, (str, bytes, bytearray))
                else []
            ),
        },
        "segments": segments,
        "checks": checks,
        "payload": payload,
    }


def _build_lr_contract_freeze_summary(
    modules: Mapping[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    required_modules = (
        "eft_deltaf2.operators",
        "eft_deltaf2.matching_kkgluon",
        "eft_deltaf2.rg",
    )
    missing_modules = [name for name in required_modules if name not in modules]
    if missing_modules:
        return [], {
            "status": "missing",
            "missing_modules": ", ".join(missing_modules),
        }

    operators_module = modules["eft_deltaf2.operators"]
    matching_module = modules["eft_deltaf2.matching_kkgluon"]
    rg_module = modules["eft_deltaf2.rg"]

    operator_basis_fn = _get_callable(
        operators_module,
        ("default_paper_0710_1869_deltaf2_operator_basis",),
    )
    if not callable(operator_basis_fn):
        return (
            [
                "operators module is present but exposes no default paper operator basis; "
                "expected default_paper_0710_1869_deltaf2_operator_basis"
            ],
            {"status": "missing_operator_basis_api"},
        )

    lr_basis_contract_fn = _get_callable(
        rg_module,
        ("default_paper_0710_1869_deltaf2_lr_basis_contract",),
    )
    if not callable(lr_basis_contract_fn):
        return (
            [
                "RG module is present but exposes no default LR basis contract; "
                "expected default_paper_0710_1869_deltaf2_lr_basis_contract"
            ],
            {"status": "missing_lr_basis_contract_api"},
        )

    rg_contract_fn = _get_callable(
        rg_module,
        ("default_paper_0710_1869_deltaf2_rg_contract",),
    )
    if not callable(rg_contract_fn):
        return (
            [
                "RG module is present but exposes no default RG contract; "
                "expected default_paper_0710_1869_deltaf2_rg_contract"
            ],
            {"status": "missing_rg_contract_api"},
        )

    try:
        operator_basis_payload = _canonicalize(operator_basis_fn())
        matching_payload = _default_matching_payload(matching_module)
        lr_basis_contract_payload = _canonicalize(lr_basis_contract_fn())
        rg_contract_payload = _canonicalize(rg_contract_fn())
    except AssertionError as exc:
        return ([str(exc)], {"status": "invalid_payload"})

    payloads = (
        operator_basis_payload,
        matching_payload,
        lr_basis_contract_payload,
        rg_contract_payload,
    )
    if not all(isinstance(payload, Mapping) for payload in payloads):
        return (
            [
                "LR-freeze acceptance summary requires mapping-shaped "
                "operator, matching, and RG payloads"
            ],
            {"status": "invalid_payload"},
        )

    operator_basis_payload = dict(operator_basis_payload)
    matching_payload = dict(matching_payload)
    lr_basis_contract_payload = dict(lr_basis_contract_payload)
    rg_contract_payload = dict(rg_contract_payload)

    deterministic = (
        json.dumps(operator_basis_payload, sort_keys=True)
        == json.dumps(_canonicalize(operator_basis_fn()), sort_keys=True)
        and json.dumps(matching_payload, sort_keys=True)
        == json.dumps(_default_matching_payload(matching_module), sort_keys=True)
        and json.dumps(lr_basis_contract_payload, sort_keys=True)
        == json.dumps(_canonicalize(lr_basis_contract_fn()), sort_keys=True)
        and json.dumps(rg_contract_payload, sort_keys=True)
        == json.dumps(_canonicalize(rg_contract_fn()), sort_keys=True)
    )

    lr_operator_entries: list[dict[str, Any]] = []
    operators = operator_basis_payload.get("operators")
    if isinstance(operators, Sequence) and not isinstance(operators, (str, bytes, bytearray)):
        lr_operator_entries = [
            dict(entry)
            for entry in operators
            if isinstance(entry, Mapping)
            and str(entry.get("name")) in EXPECTED_UNSUPPORTED_OPERATORS
        ]

    paper_lr_operator_names = [str(entry.get("name")) for entry in lr_operator_entries]
    paper_lr_definition_ids = [
        str(entry.get("definition_id"))
        for entry in lr_operator_entries
        if entry.get("definition_id") is not None
    ]
    projector_normalization_id = operator_basis_payload.get("projector_normalization_id")
    projector_normalization_note = operator_basis_payload.get("projector_normalization_note")
    paper_operator_order = _string_list(lr_basis_contract_payload.get("paper_operator_order"))
    bmu_operator_order = _string_list(lr_basis_contract_payload.get("bmu_lr_operator_order"))
    paper_to_bmu_operator_map = _real_matrix_2x2(
        lr_basis_contract_payload.get("paper_to_bmu_operator_map_matrix")
    )
    bmu_to_paper_operator_map = _real_matrix_2x2(
        lr_basis_contract_payload.get("bmu_to_paper_operator_map_matrix")
    )
    paper_to_bmu_wilson_map = _real_matrix_2x2(
        lr_basis_contract_payload.get("paper_to_bmu_wilson_map_matrix")
    )
    bmu_to_paper_wilson_map = _real_matrix_2x2(
        lr_basis_contract_payload.get("bmu_to_paper_wilson_map_matrix")
    )
    mapping_matrix_frozen = lr_basis_contract_payload.get("mapping_matrix_frozen")
    matching_supported_operator_names = _string_list(
        matching_payload.get("supported_observable_operator_names")
    )
    matching_unsupported_operator_names = _string_list(
        matching_payload.get("unsupported_observable_operator_names")
    )
    rg_supported_operator_names = _string_list(rg_contract_payload.get("supported_operator_names"))
    rg_unsupported_operator_names = _string_list(
        rg_contract_payload.get("unsupported_operator_names")
    )
    rg_lr_basis_map_supported = _nested_value(
        _default_rg_payload(modules),
        (
            ("lr_basis_map_supported",),
            ("evolution", "lr_basis_map_supported"),
        ),
    )

    checks = {
        "source_payloads_are_deterministic": deterministic,
        "paper_lr_operator_order_is_frozen": (
            tuple(paper_lr_operator_names) == EXPECTED_GUARDED_LR_OPERATOR_ORDER
        ),
        "paper_lr_definition_ids_are_explicit": (
            len(paper_lr_definition_ids) == len(EXPECTED_GUARDED_LR_OPERATOR_ORDER)
            and all(paper_lr_definition_ids)
        ),
        "paper_operator_order_matches_frozen_map": (
            tuple(paper_operator_order) == EXPECTED_GUARDED_LR_OPERATOR_ORDER
        ),
        "bmu_operator_order_matches_frozen_map": (
            tuple(bmu_operator_order) == EXPECTED_BMU_LR_OPERATOR_ORDER
        ),
        "projector_normalization_id_present": bool(projector_normalization_id),
        "projector_normalization_note_mentions_projector_convention": (
            isinstance(projector_normalization_note, str)
            and "P_L = (1-gamma5)/2" in projector_normalization_note
            and "P_R = (1+gamma5)/2" in projector_normalization_note
        ),
        "lr_contract_references_paper_basis": (
            lr_basis_contract_payload.get("paper_operator_basis_id")
            == operator_basis_payload.get("basis_id")
        ),
        "lr_contract_references_paper_lr_definitions": (
            paper_lr_definition_ids
            == _string_list(lr_basis_contract_payload.get("paper_operator_definition_ids"))
        ),
        "lr_definitions_are_frozen": lr_basis_contract_payload.get("lr_definitions_frozen") is True,
        "mapping_matrix_is_frozen": mapping_matrix_frozen is True,
        "paper_to_bmu_operator_map_is_exact": (
            paper_to_bmu_operator_map == EXPECTED_PAPER_TO_BMU_OPERATOR_MAP
        ),
        "bmu_to_paper_operator_map_is_exact": (
            bmu_to_paper_operator_map == EXPECTED_BMU_TO_PAPER_OPERATOR_MAP
        ),
        "paper_to_bmu_wilson_map_is_exact": (
            paper_to_bmu_wilson_map == EXPECTED_PAPER_TO_BMU_WILSON_MAP
        ),
        "bmu_to_paper_wilson_map_is_exact": (
            bmu_to_paper_wilson_map == EXPECTED_BMU_TO_PAPER_WILSON_MAP
        ),
        "operator_map_round_trip_is_identity": (
            paper_to_bmu_operator_map is not None
            and bmu_to_paper_operator_map is not None
            and _real_matrix_multiply(bmu_to_paper_operator_map, paper_to_bmu_operator_map)
            == ((1.0, 0.0), (0.0, 1.0))
            and _real_matrix_multiply(paper_to_bmu_operator_map, bmu_to_paper_operator_map)
            == ((1.0, 0.0), (0.0, 1.0))
        ),
        "wilson_map_round_trip_is_identity": (
            paper_to_bmu_wilson_map is not None
            and bmu_to_paper_wilson_map is not None
            and _real_matrix_multiply(bmu_to_paper_wilson_map, paper_to_bmu_wilson_map)
            == ((1.0, 0.0), (0.0, 1.0))
            and _real_matrix_multiply(paper_to_bmu_wilson_map, bmu_to_paper_wilson_map)
            == ((1.0, 0.0), (0.0, 1.0))
        ),
        "operator_basis_vectors_map_as_frozen": (
            paper_to_bmu_operator_map is not None
            and _real_matrix_vector_multiply(paper_to_bmu_operator_map, (1.0, 0.0))
            == (0.0, 1.0)
            and _real_matrix_vector_multiply(paper_to_bmu_operator_map, (0.0, 1.0))
            == (2.0, 0.0)
        ),
        "wilson_basis_vectors_map_as_frozen": (
            paper_to_bmu_wilson_map is not None
            and _real_matrix_vector_multiply(paper_to_bmu_wilson_map, (1.0, 0.0))
                == (0.0, 1.0)
                and _real_matrix_vector_multiply(paper_to_bmu_wilson_map, (0.0, 1.0))
                == (0.5, 0.0)
            ),
        "lr_running_is_active": (
            lr_basis_contract_payload.get("lr_running_activated") is True
        ),
        "matching_references_lr_contract": (
            matching_payload.get("lr_observable_support_contract_id")
            == lr_basis_contract_payload.get("contract_id")
        ),
        "matching_references_lr_contract_status": (
            matching_payload.get("lr_basis_status_id") == lr_basis_contract_payload.get("status_id")
        ),
        "matching_supported_subset_stays_q1_only": (
            tuple(matching_supported_operator_names) == EXPECTED_SUPPORTED_OPERATOR_ORDER
        ),
        "matching_unsupported_subset_is_lr_only": (
            tuple(matching_unsupported_operator_names) == EXPECTED_GUARDED_LR_OPERATOR_ORDER
        ),
        "rg_contract_references_lr_contract": (
            rg_contract_payload.get("lr_basis_contract_id")
            == lr_basis_contract_payload.get("contract_id")
        ),
        "rg_contract_references_lr_contract_status": (
            rg_contract_payload.get("lr_basis_status_id")
            == lr_basis_contract_payload.get("status_id")
        ),
        "rg_supported_subset_is_lr_running_surface": (
            tuple(rg_supported_operator_names)
            == EXPECTED_SUPPORTED_OPERATOR_ORDER + EXPECTED_GUARDED_LR_OPERATOR_ORDER
        ),
        "rg_unsupported_subset_is_empty": (
            tuple(rg_unsupported_operator_names) == tuple()
        ),
    }
    failures = [
        f"LR-freeze acceptance check failed: {name}"
        for name, ok in checks.items()
        if not ok
    ]

    return failures, {
        "status": "ok",
        "deterministic": deterministic,
        "paper_operator_basis_id": operator_basis_payload.get("basis_id"),
        "paper_lr_operator_names": paper_lr_operator_names,
        "paper_lr_definition_ids": paper_lr_definition_ids,
        "projector_normalization_id": projector_normalization_id,
        "projector_normalization_note": projector_normalization_note,
        "lr_basis_contract_id": lr_basis_contract_payload.get("contract_id"),
        "lr_basis_status_id": lr_basis_contract_payload.get("status_id"),
        "paper_operator_order": paper_operator_order,
        "bmu_lr_basis_id": lr_basis_contract_payload.get("bmu_lr_basis_id"),
        "bmu_lr_operator_order": bmu_operator_order,
        "paper_to_bmu_operator_map_matrix": (
            [list(row) for row in paper_to_bmu_operator_map]
            if paper_to_bmu_operator_map is not None
            else None
        ),
        "bmu_to_paper_operator_map_matrix": (
            [list(row) for row in bmu_to_paper_operator_map]
            if bmu_to_paper_operator_map is not None
            else None
        ),
        "paper_to_bmu_wilson_map_matrix": (
            [list(row) for row in paper_to_bmu_wilson_map]
            if paper_to_bmu_wilson_map is not None
            else None
        ),
        "bmu_to_paper_wilson_map_matrix": (
            [list(row) for row in bmu_to_paper_wilson_map]
            if bmu_to_paper_wilson_map is not None
            else None
        ),
        "mapping_matrix_frozen": mapping_matrix_frozen,
        "matching_lr_support_contract_id": matching_payload.get(
            "lr_observable_support_contract_id"
        ),
        "matching_lr_basis_status_id": matching_payload.get("lr_basis_status_id"),
        "matching_lr_support_status_id": matching_payload.get(
            "lr_observable_support_status_id"
        ),
        "lr_running_activated": lr_basis_contract_payload.get("lr_running_activated"),
        "rg_lr_basis_map_supported": rg_lr_basis_map_supported,
        "rg_supported_operator_subset_id": rg_contract_payload.get(
            "supported_operator_subset_id"
        ),
        "matching_supported_observable_operator_names": matching_supported_operator_names,
        "matching_unsupported_observable_operator_names": matching_unsupported_operator_names,
        "rg_supported_operator_names": rg_supported_operator_names,
        "rg_unsupported_operator_names": rg_unsupported_operator_names,
        "checks": checks,
    }


def _build_eft_lr_running_probe_summary(
    modules: Mapping[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    rg_module = modules.get("eft_deltaf2.rg")
    operators_module = modules.get("eft_deltaf2.operators")
    if rg_module is None or operators_module is None:
        return [], {"status": "missing"}

    required_callables = {
        "contract": getattr(rg_module, "default_paper_0710_1869_deltaf2_rg_contract", None),
        "evolution": getattr(rg_module, "compute_deltaf2_lo_evolution_matrix", None),
        "paper_to_bmu": getattr(rg_module, "map_paper_lr_wilsons_to_bmu_lr", None),
        "bmu_to_paper": getattr(rg_module, "map_bmu_lr_wilsons_to_paper_lr", None),
        "paper_to_bmu_w": getattr(rg_module, "paper_lr_to_bmu_lr_wilson_map_matrix", None),
        "bmu_to_paper_w": getattr(rg_module, "bmu_lr_to_paper_lr_wilson_map_matrix", None),
        "evolve": _get_callable(rg_module, PARAMETERIZED_RG_EXPORT_NAMES),
        "wilson_contract": getattr(
            operators_module,
            "default_paper_0710_1869_deltaf2_wilson_contract",
            None,
        ),
        "wilson_type": getattr(
            operators_module,
            "Paper07101869DeltaF2WilsonCoefficients",
            None,
        ),
    }
    missing = [
        name
        for name, callable_obj in required_callables.items()
        if not callable(callable_obj)
    ]
    if missing:
        return [], {
            "status": "missing_api",
            "missing": sorted(missing),
        }

    probe_input = (complex(1.25, -0.5), complex(-0.75, 0.25))
    try:
        contract = required_callables["contract"]()
        lr_basis_contract = rg_module.default_paper_0710_1869_deltaf2_lr_basis_contract()
        evolution = required_callables["evolution"](
            3000.0,
            DEFAULT_RG_LOW_SCALE_GEV,
            contract=contract,
        )
        paper_to_bmu_w = required_callables["paper_to_bmu_w"]()
        bmu_to_paper_w = required_callables["bmu_to_paper_w"]()
        mapped_bmu_input = required_callables["paper_to_bmu"](probe_input)
        evolved_bmu_vector = _complex_matrix_vector_multiply(
            evolution.bmu_lr_evolution_matrix,
            mapped_bmu_input,
        )
        mapped_back_paper_vector = required_callables["bmu_to_paper"](evolved_bmu_vector)
        segment_product = (
            (1.0 + 0.0j, 0.0 + 0.0j),
            (0.0 + 0.0j, 1.0 + 0.0j),
        )
        for segment in evolution.segments:
            segment_product = _complex_matrix_multiply(
                segment.bmu_lr_matrix,
                segment_product,
            )
        expected_paper_lr_block = _complex_matrix_multiply(
            bmu_to_paper_w,
            _complex_matrix_multiply(evolution.bmu_lr_evolution_matrix, paper_to_bmu_w),
        )
        paper_lr_block = _lr_block_from_paper_basis_matrix(evolution.paper_basis_evolution_matrix)
        evolved_paper_from_block = (
            _complex_matrix_vector_multiply(paper_lr_block, probe_input)
            if paper_lr_block is not None
            else None
        )

        toy_wilsons = required_callables["wilson_type"](
            contract=required_callables["wilson_contract"](),
            benchmark_id="benchmark.synthetic_lr_probe.v1",
            scale_label="synthetic_lr_probe",
            system_id="kaon",
            sector_id="down",
            generations=(0, 1),
            matching_scale_GeV=3000.0,
            propagator_mass_GeV=3000.0,
            left_coupling=0.0 + 0.0j,
            right_coupling=0.0 + 0.0j,
            q1_vll=0.0 + 0.0j,
            q1_vrr=0.0 + 0.0j,
            q4_lr=probe_input[0],
            q5_lr=probe_input[1],
        )
        public_run_error = None
        public_rg_accepts_nonzero_lr = False
        public_evolved_paper_vector = None
        public_rg_matches_expected = None
        try:
            public_result = required_callables["evolve"](
                toy_wilsons,
                mu_low_GeV=DEFAULT_RG_LOW_SCALE_GEV,
            )
            public_coefficients = getattr(public_result, "coefficients", None)
            if isinstance(public_coefficients, Mapping):
                public_evolved_paper_vector = (
                    complex(public_coefficients["Q4_LR"]),
                    complex(public_coefficients["Q5_LR"]),
                )
                public_rg_accepts_nonzero_lr = True
                public_rg_matches_expected = _complex_vector_close(
                    public_evolved_paper_vector,
                    mapped_back_paper_vector,
                )
        except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
            public_run_error = f"{type(exc).__name__}: {exc}"

        return [], {
            "status": "ok",
            "probe_input_paper_lr_wilsons": _complex_vector_as_dict_list(probe_input),
            "mapped_input_bmu_lr_wilsons": _complex_vector_as_dict_list(mapped_bmu_input),
            "bmu_lr_evolution_matrix": _complex_matrix_as_dict(
                evolution.bmu_lr_evolution_matrix
            ),
            "segment_product_bmu_lr_matrix": _complex_matrix_as_dict(segment_product),
            "expected_paper_lr_block": _complex_matrix_as_dict(expected_paper_lr_block),
            "paper_lr_block_from_public_evolution": (
                _complex_matrix_as_dict(paper_lr_block)
                if paper_lr_block is not None
                else None
            ),
            "evolved_bmu_lr_wilsons": _complex_vector_as_dict_list(evolved_bmu_vector),
            "mapped_back_paper_lr_wilsons": _complex_vector_as_dict_list(
                mapped_back_paper_vector
            ),
            "paper_lr_block_times_input": (
                _complex_vector_as_dict_list(evolved_paper_from_block)
                if evolved_paper_from_block is not None
                else None
            ),
            "public_evolved_paper_lr_wilsons": (
                _complex_vector_as_dict_list(public_evolved_paper_vector)
                if public_evolved_paper_vector is not None
                else None
            ),
            "public_run_error": public_run_error,
            "lr_running_activated": bool(
                getattr(lr_basis_contract, "lr_running_activated", False)
            ),
            "lr_basis_map_supported": bool(
                getattr(evolution, "lr_basis_map_supported", False)
            ),
            "checks": {
                "segment_product_matches_bmu_total": _complex_matrix_close(
                    segment_product,
                    evolution.bmu_lr_evolution_matrix,
                ),
                "paper_block_matches_winv_u_bmu_w": (
                    paper_lr_block is not None
                    and _complex_matrix_close(paper_lr_block, expected_paper_lr_block)
                ),
                "mapped_back_vector_matches_conjugated_block": (
                    evolved_paper_from_block is not None
                    and _complex_vector_close(
                        evolved_paper_from_block,
                        mapped_back_paper_vector,
                    )
                ),
                "public_rg_accepts_nonzero_lr": public_rg_accepts_nonzero_lr,
                "public_rg_matches_winv_u_bmu_w": public_rg_matches_expected is True,
            },
        }
    except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
        return [], {
            "status": "error",
            "error": f"{type(exc).__name__}: {exc}",
        }


def _build_hadronic_summary(modules: Mapping[str, Any]) -> tuple[list[str], dict[str, Any]]:
    module = modules.get("eft_deltaf2.hadronic")
    if module is None:
        return [], {"status": "missing"}

    try:
        payload = _default_hadronic_payload(modules)
    except AssertionError as exc:
        return ([str(exc)], {"status": "missing_default_hadronic_export"})

    second = _default_hadronic_payload(modules)
    builder_default_payload = _default_hadronic_builder_payload(modules)
    deterministic = json.dumps(payload, sort_keys=True) == json.dumps(second, sort_keys=True)
    failures = []
    if not deterministic:
        failures.append("default hadronic export is not deterministic within one process")
    if builder_default_payload != payload:
        failures.append("default hadronic export does not match the builder default payload")
    cross_process_payload = _cross_process_default_hadronic_payload()
    cross_process_deterministic = json.dumps(payload, sort_keys=True) == json.dumps(
        cross_process_payload,
        sort_keys=True,
    )
    if not cross_process_deterministic:
        failures.append("default hadronic export is not deterministic across processes")

    schema_tag = _nested_value(
        payload,
        (
            ("schema_id",),
            ("metadata", "schema_name"),
            ("contract", "schema_id"),
        ),
    )
    mass_source_id = _nested_value(
        payload,
        (
            ("mass_source", "source_id"),
        ),
    )
    system_id = _nested_value(
        payload,
        (
            ("system_id",),
            ("system",),
            ("tags", "system_id"),
            ("contract", "system_id"),
        ),
    )
    scheme_id = _nested_value(
        payload,
        (
            ("renormalization_scheme_id",),
            ("scheme_id",),
            ("tags", "renormalization_scheme_id"),
            ("contract", "renormalization_scheme_id"),
        ),
    )
    contract_system_id = _nested_value(
        payload,
        (
            ("contract", "system_id"),
        ),
    )
    contract_scheme_id = _nested_value(
        payload,
        (
            ("contract", "renormalization_scheme_id"),
            ("contract", "scheme_id"),
        ),
    )
    normalization_id = _nested_value(
        payload,
        (
            ("operator_normalization_id",),
            ("normalization_id",),
            ("tags", "operator_normalization_id"),
            ("contract", "operator_normalization_id"),
        ),
    )
    contract_normalization_id = _nested_value(
        payload,
        (
            ("contract", "operator_normalization_id"),
        ),
    )
    basis_id = _nested_value(
        payload,
        (
            ("operator_basis_id",),
            ("basis_id",),
            ("tags", "operator_basis_id"),
            ("contract", "operator_basis_id"),
        ),
    )
    contract_basis_id = _nested_value(
        payload,
        (
            ("contract", "operator_basis_id"),
        ),
    )
    mu_had_GeV = _nested_value(
        payload,
        (
            ("mu_had_GeV",),
            ("evaluation_scale_GeV",),
            ("hadronic_scale_GeV",),
            ("tags", "mu_had_GeV"),
            ("tags", "evaluation_scale_GeV"),
            ("contract", "mu_had_GeV"),
        ),
    )
    contract_mu_had_GeV = _nested_value(
        payload,
        (
            ("contract", "mu_had_GeV"),
            ("contract", "evaluation_scale_GeV"),
        ),
    )
    decay_constant_source_id = _nested_value(
        payload,
        (
            ("decay_constant_source", "source_id"),
        ),
    )
    bag_parameter_source_id = _nested_value(
        payload,
        (
            ("bag_parameter_source", "source_id"),
        ),
    )
    bundle_source_id = _nested_value(
        payload,
        (
            ("source_id",),
            ("tags", "source_id"),
        ),
    )
    provenance_ids = _nested_value(
        payload,
        (
            ("provenance_ids",),
            ("tags", "provenance_ids"),
        ),
    )
    parity_relation_id = _nested_value(
        payload,
        (
            ("parity_relation_id",),
            ("tags", "parity_relation_id"),
            ("contract", "parity_relation_id"),
        ),
    )
    contract_parity_relation_id = _nested_value(
        payload,
        (
            ("contract", "parity_relation_id"),
        ),
    )
    bag_parameter_source_scheme_id = _nested_value(
        payload,
        (
            ("bag_parameter_source_scheme_id",),
        ),
    )
    bag_parameter_source_scheme = _nested_value(
        payload,
        (
            ("bag_parameter_source", "renormalization_scheme_id"),
        ),
    )
    bag_parameter_transformation_id = _nested_value(
        payload,
        (
            ("bag_parameter_transformation_id",),
        ),
    )
    bag_parameter_source_transformation_id = _nested_value(
        payload,
        (
            ("bag_parameter_source", "transformation_id"),
        ),
    )
    operator_names = _nested_value(
        payload,
        (
            ("supported_operator_names",),
            ("contract", "supported_operator_names"),
        ),
    )
    unsupported_operator_names = _nested_value(
        payload,
        (
            ("unsupported_operator_names",),
            ("contract", "unsupported_operator_names"),
        ),
    )

    checks = {
        "schema_tag_present": bool(schema_tag),
        "schema_tag_mentions_hadronic": schema_tag is not None
        and "hadronic" in str(schema_tag).lower(),
        "system_is_kaon": _is_kaon_system(system_id),
        "scheme_id_present": bool(scheme_id),
        "operator_basis_id_present": bool(basis_id),
        "operator_normalization_id_present": bool(normalization_id),
        "contract_system_matches_top_level": contract_system_id == system_id,
        "contract_basis_matches_top_level": contract_basis_id == basis_id,
        "contract_normalization_matches_top_level": (
            contract_normalization_id == normalization_id
        ),
        "contract_scheme_matches_top_level": contract_scheme_id == scheme_id,
        "contract_mu_had_matches_top_level": (
            _positive_finite(contract_mu_had_GeV)
            and float(contract_mu_had_GeV) == float(mu_had_GeV)
        ),
        "mu_had_present_and_positive": _positive_finite(mu_had_GeV),
        "source_ids_present": all(
            bool(source_id)
            for source_id in (
                mass_source_id,
                decay_constant_source_id,
                bag_parameter_source_id,
            )
        ),
        "bundle_source_present": bool(bundle_source_id),
        "bundle_source_in_provenance_ids": bool(bundle_source_id)
        and isinstance(provenance_ids, Sequence)
        and not isinstance(provenance_ids, (str, bytes, bytearray))
        and bundle_source_id in {str(item) for item in provenance_ids},
        "parity_relation_tag_present": bool(parity_relation_id),
        "parity_relation_matches_supported_relation": bool(parity_relation_id)
        and "q1_vll_equals_q1_vrr" in str(parity_relation_id).lower(),
        "contract_parity_relation_matches_top_level": (
            contract_parity_relation_id == parity_relation_id
        ),
        "bag_source_scheme_matches_bundle": (
            bag_parameter_source_scheme == bag_parameter_source_scheme_id
        ),
        "bag_source_transformation_matches_bundle": (
            bag_parameter_source_transformation_id == bag_parameter_transformation_id
        ),
        "default_export_matches_builder_default": builder_default_payload == payload,
        "supported_operator_subset_is_pr5a": isinstance(operator_names, Sequence)
        and not isinstance(operator_names, (str, bytes, bytearray))
        and {str(item) for item in operator_names} == EXPECTED_SUPPORTED_OPERATORS,
        "unsupported_operator_subset_is_lr_only": isinstance(unsupported_operator_names, Sequence)
        and not isinstance(unsupported_operator_names, (str, bytes, bytearray))
        and {str(item) for item in unsupported_operator_names}
        == EXPECTED_UNSUPPORTED_OPERATORS,
    }
    failures.extend(
        f"hadronic acceptance check failed: {name}" for name, ok in checks.items() if not ok
    )

    return failures, {
        "status": "ok",
        "deterministic": deterministic,
        "cross_process_deterministic": cross_process_deterministic,
        "schema_tag": schema_tag,
        "system_id": system_id,
        "scheme_id": scheme_id,
        "operator_basis_id": basis_id,
        "operator_normalization_id": normalization_id,
        "parity_relation_id": parity_relation_id,
        "bag_parameter_source_scheme_id": bag_parameter_source_scheme_id,
        "bag_parameter_transformation_id": bag_parameter_transformation_id,
        "mu_had_GeV": float(mu_had_GeV) if mu_had_GeV is not None else None,
        "checks": checks,
        "payload": payload,
    }


def _build_default_lr_hadronic_probe_summary(
    modules: Mapping[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    hadronic_module = modules.get("eft_deltaf2.hadronic")
    if hadronic_module is None:
        return [], {"status": "missing"}

    try:
        payload = _default_kaon_lr_hadronic_payload(modules)
        second_payload = _default_kaon_lr_hadronic_payload(modules)
        default_lr_value = _default_kaon_lr_hadronic_value(modules)
        builder_payload = _default_kaon_lr_hadronic_builder_payload(modules)
        cross_process_payload = _cross_process_default_lr_hadronic_payload()
    except AssertionError as exc:
        return ([str(exc)], {"status": "missing_default_lr_hadronic_export"})
    except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
        return ([f"{type(exc).__name__}: {exc}"], {"status": "error", "error": str(exc)})

    payload_canonical_json, payload_sha256 = _canonical_payload_json_and_sha256(payload)
    second_payload_canonical_json, second_payload_sha256 = (
        _canonical_payload_json_and_sha256(second_payload)
    )
    cross_process_payload_canonical_json, cross_process_payload_sha256 = (
        _canonical_payload_json_and_sha256(cross_process_payload)
    )

    system_id = _nested_value(payload, (("system_id",),))
    scheme_id = _nested_value(
        payload,
        (
            ("scheme_id",),
            ("renormalization_scheme_id",),
            ("contract", "renormalization_scheme_id"),
        ),
    )
    mu_had_GeV = _nested_value(
        payload,
        (
            ("mu_had_GeV",),
            ("evaluation_scale_GeV",),
            ("contract", "mu_had_GeV"),
        ),
    )
    bundle_id = _nested_value(payload, (("bundle_id",),))
    source_id = _nested_value(payload, (("source_id",),))
    basis_id = _nested_value(
        payload,
        (
            ("operator_basis_id",),
            ("contract", "operator_basis_id"),
        ),
    )
    normalization_id = _nested_value(
        payload,
        (
            ("operator_normalization_id",),
            ("contract", "operator_normalization_id"),
        ),
    )
    hamiltonian_convention_id = _nested_value(
        payload,
        (
            ("hamiltonian_convention_id",),
            ("contract", "hamiltonian_convention_id"),
        ),
    )
    contract_basis_id = _nested_value(payload, (("contract", "operator_basis_id"),))
    contract_normalization_id = _nested_value(
        payload,
        (("contract", "operator_normalization_id"),),
    )
    contract_hamiltonian_convention_id = _nested_value(
        payload,
        (("contract", "hamiltonian_convention_id"),),
    )
    provenance_ids = _nested_value(
        payload,
        (
            ("provenance_ids",),
            ("tags", "provenance_ids"),
        ),
    )
    input_policy_id = _nested_value(
        payload,
        (
            ("input_policy_id",),
            ("contract", "input_policy_id"),
            ("tags", "input_policy_id"),
        ),
    )
    b4_mu_had = _nested_value(payload, (("B4_mu_had",),))
    b5_mu_had = _nested_value(payload, (("B5_mu_had",),))
    r_chi_mu_had = _nested_value(payload, (("R_chi_mu_had",),))
    m_K0_GeV = _nested_value(payload, (("m_K0_GeV",),))
    f_K_GeV = _nested_value(payload, (("f_K_GeV",),))
    q4_formula_id = _nested_value(
        payload,
        (
            ("q4_matrix_element_formula_id",),
            ("contract", "q4_matrix_element_formula_id"),
        ),
    )
    q5_formula_id = _nested_value(
        payload,
        (
            ("q5_matrix_element_formula_id",),
            ("contract", "q5_matrix_element_formula_id"),
        ),
    )
    q4_matrix_element_GeV4 = _extract_lr_matrix_element(payload, "Q4_LR")
    q5_matrix_element_GeV4 = _extract_lr_matrix_element(payload, "Q5_LR")
    supported_operator_names = _string_list(
        _nested_value(
            payload,
            (
                ("supported_operator_names",),
                ("contract", "supported_operator_names"),
            ),
        )
    )
    unsupported_operator_names = _string_list(
        _nested_value(
            payload,
            (
                ("unsupported_operator_names",),
                ("contract", "unsupported_operator_names"),
            ),
        )
    )

    b4_source_id = _nested_value(payload, (("b4_source", "source_id"),))
    b5_source_id = _nested_value(payload, (("b5_source", "source_id"),))
    b4_source_transformation_id = _nested_value(
        payload,
        (("b4_source", "transformation_id"),),
    )
    b5_source_transformation_id = _nested_value(
        payload,
        (("b5_source", "transformation_id"),),
    )
    r_chi_source_id = _nested_value(payload, (("r_chi_source", "source_id"),))
    mass_source_id = _nested_value(payload, (("mass_source", "source_id"),))
    decay_constant_source_id = _nested_value(
        payload,
        (("decay_constant_source", "source_id"),),
    )
    b4_source_scheme_id = _nested_value(payload, (("b4_source", "renormalization_scheme_id"),))
    b5_source_scheme_id = _nested_value(payload, (("b5_source", "renormalization_scheme_id"),))
    r_chi_source_scheme_id = _nested_value(
        payload,
        (("r_chi_source", "renormalization_scheme_id"),),
    )
    b4_source_scale_GeV = _nested_value(payload, (("b4_source", "scale_GeV"),))
    b5_source_scale_GeV = _nested_value(payload, (("b5_source", "scale_GeV"),))
    r_chi_source_scale_GeV = _nested_value(payload, (("r_chi_source", "scale_GeV"),))
    b4_source_metadata_text = _source_metadata_text(payload, "b4_source")
    b5_source_metadata_text = _source_metadata_text(payload, "b5_source")

    expected_q4 = None
    expected_q5 = None
    if all(
        _positive_finite(value)
        for value in (b4_mu_had, r_chi_mu_had, m_K0_GeV, f_K_GeV)
    ):
        expected_q4 = (
            2.0
            * float(r_chi_mu_had)
            * (float(m_K0_GeV) ** 2)
            * (float(f_K_GeV) ** 2)
            * float(b4_mu_had)
        )
    if all(
        _positive_finite(value)
        for value in (b5_mu_had, r_chi_mu_had, m_K0_GeV, f_K_GeV)
    ):
        expected_q5 = (
            (2.0 / 3.0)
            * float(r_chi_mu_had)
            * (float(m_K0_GeV) ** 2)
            * (float(f_K_GeV) ** 2)
            * float(b5_mu_had)
        )

    custom_lr_only_rejection_error = None
    custom_combined_rejection_error = None
    custom_lr_only_rejects_default_lr_bundle = False
    custom_combined_rejects_default_lr_bundle = False
    observable_module = modules.get("eft_deltaf2.observables")
    if (
        observable_module is not None
        and bool(scheme_id)
        and _positive_finite(mu_had_GeV)
    ):
        lr_evaluator = _get_callable(observable_module, CUSTOM_LR_OBSERVABLE_EXPORT_NAMES)
        if callable(lr_evaluator):
            try:
                probe_rg_value = _replace_rg_lr_wilsons(
                    _default_rg_value(modules),
                    q4_lr=LR_OBSERVABLE_PROBE_Q4,
                    q5_lr=LR_OBSERVABLE_PROBE_Q5,
                    renormalization_scheme_id=str(scheme_id),
                    matching_scale_GeV=float(mu_had_GeV),
                )
                _invoke_observable_evaluation(
                    lr_evaluator,
                    rg_value=probe_rg_value,
                    hadronic_value=default_lr_value,
                )
            except ValueError as exc:
                custom_lr_only_rejection_error = f"{type(exc).__name__}: {exc}"
                custom_lr_only_rejects_default_lr_bundle = _is_custom_input_only_guard_error(
                    custom_lr_only_rejection_error
                )
            except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
                custom_lr_only_rejection_error = f"{type(exc).__name__}: {exc}"

        combined_evaluator = _get_callable(
            observable_module,
            CUSTOM_TOTAL_OBSERVABLE_EXPORT_NAMES,
        )
        if callable(combined_evaluator):
            try:
                custom_q1_hadronic = _build_custom_q1_hadronic_bundle(
                    modules,
                    mu_had_GeV=float(mu_had_GeV),
                )
                combined_rg_value = _replace_rg_lr_wilsons(
                    _default_rg_value(modules),
                    q4_lr=LR_OBSERVABLE_PROBE_Q4,
                    q5_lr=LR_OBSERVABLE_PROBE_Q5,
                    renormalization_scheme_id=str(scheme_id),
                    matching_scale_GeV=float(mu_had_GeV),
                )
                _invoke_custom_total_observable_evaluation(
                    combined_evaluator,
                    rg_value=combined_rg_value,
                    q1_hadronic_bundle=custom_q1_hadronic,
                    lr_hadronic_inputs=default_lr_value,
                )
            except ValueError as exc:
                custom_combined_rejection_error = f"{type(exc).__name__}: {exc}"
                custom_combined_rejects_default_lr_bundle = _is_custom_input_only_guard_error(
                    custom_combined_rejection_error
                )
            except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
                custom_combined_rejection_error = f"{type(exc).__name__}: {exc}"

    default_lr_rejection_checks: dict[str, bool] = {}
    default_lr_rejection_errors: dict[str, str] = {}
    for check_name, mutation in (
        (
            "drifted_b4_value_rejected",
            lambda bundle: {"B4_mu_had": float(bundle.B4_mu_had) + 0.01},
        ),
        (
            "drifted_b5_value_rejected",
            lambda bundle: {"B5_mu_had": float(bundle.B5_mu_had) + 0.01},
        ),
        (
            "drifted_source_id_rejected",
            lambda bundle: {"source_id": f"{bundle.source_id}.drift"},
        ),
        (
            "drifted_b4_source_id_rejected",
            lambda bundle: {
                "b4_source": dataclasses.replace(
                    bundle.b4_source,
                    source_id=f"{bundle.b4_source.source_id}.drift",
                )
            },
        ),
        (
            "drifted_b5_source_id_rejected",
            lambda bundle: {
                "b5_source": dataclasses.replace(
                    bundle.b5_source,
                    source_id=f"{bundle.b5_source.source_id}.drift",
                )
            },
        ),
        (
            "drifted_r_chi_source_id_rejected",
            lambda bundle: {
                "r_chi_source": dataclasses.replace(
                    bundle.r_chi_source,
                    source_id=f"{bundle.r_chi_source.source_id}.drift",
                )
            },
        ),
        (
            "drifted_bundle_scheme_rejected",
            lambda bundle: {
                "renormalization_scheme_id": f"{bundle.renormalization_scheme_id}.drift"
            },
        ),
        (
            "forced_r_chi_source_bundle_scheme_rejected",
            lambda bundle: {
                "r_chi_source": dataclasses.replace(
                    bundle.r_chi_source,
                    renormalization_scheme_id=str(bundle.renormalization_scheme_id),
                )
            },
        ),
        (
            "drifted_mu_had_rejected",
            lambda bundle: {"mu_had_GeV": float(bundle.mu_had_GeV) + 0.25},
        ),
        (
            "hidden_conversion_policy_drift_rejected",
            lambda _bundle: {
                "input_policy_id": "default_source.etm2013.table1.ms_2gev.hidden_conversion.v1"
            },
        ),
        (
            "hidden_r_chi_conversion_drift_rejected",
            lambda bundle: {
                "r_chi_source": dataclasses.replace(
                    bundle.r_chi_source,
                    transformation_id="ri_mom_to_ms.hidden_conversion_probe.v1",
                )
            },
        ),
    ):
        try:
            dataclasses.replace(default_lr_value, **mutation(default_lr_value))
        except ValueError as exc:
            default_lr_rejection_checks[check_name] = True
            default_lr_rejection_errors[check_name] = f"{type(exc).__name__}: {exc}"
        except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
            default_lr_rejection_checks[check_name] = False
            default_lr_rejection_errors[check_name] = f"{type(exc).__name__}: {exc}"
        else:  # pragma: no cover - surfaced in acceptance payload
            default_lr_rejection_checks[check_name] = False
            default_lr_rejection_errors[check_name] = "mutation unexpectedly accepted"

    export_matches_builder_default = builder_payload == payload
    checks = {
        "probe_payload_is_deterministic": payload_canonical_json == second_payload_canonical_json,
        "probe_payload_is_cross_process_deterministic": (
            payload_canonical_json == cross_process_payload_canonical_json
        ),
        "probe_payload_sha256_matches_repeat": payload_sha256 == second_payload_sha256,
        "probe_payload_sha256_matches_cross_process": (
            payload_sha256 == cross_process_payload_sha256
        ),
        "system_is_kaon": _is_kaon_system(system_id),
        "mu_had_is_2gev": _positive_finite(mu_had_GeV)
        and math.isclose(float(mu_had_GeV), 2.0, rel_tol=0.0, abs_tol=1.0e-12),
        "default_export_matches_builder_default": export_matches_builder_default,
        "bundle_id_is_exact": bundle_id == EXPECTED_DEFAULT_LR_HADRONIC_BUNDLE_ID,
        "source_id_is_exact": source_id == EXPECTED_DEFAULT_LR_HADRONIC_SOURCE_ID,
        "operator_basis_id_is_exact": basis_id == EXPECTED_DEFAULT_LR_OPERATOR_BASIS_ID,
        "operator_normalization_id_is_exact": (
            normalization_id == EXPECTED_DEFAULT_LR_OPERATOR_NORMALIZATION_ID
        ),
        "hamiltonian_convention_id_is_exact": (
            hamiltonian_convention_id == EXPECTED_DEFAULT_LR_HAMILTONIAN_CONVENTION_ID
        ),
        "input_policy_id_is_exact": (
            input_policy_id == EXPECTED_DEFAULT_LR_HADRONIC_INPUT_POLICY_ID
        ),
        "provenance_ids_match_frozen_source_package": isinstance(provenance_ids, Sequence)
        and not isinstance(provenance_ids, (str, bytes, bytearray))
        and list(provenance_ids) == EXPECTED_DEFAULT_LR_HADRONIC_PROVENANCE_IDS,
        "supported_operator_subset_is_lr_only": (
            set(supported_operator_names) == EXPECTED_DEFAULT_LR_SUPPORTED_OPERATORS
        ),
        "unsupported_operator_subset_is_q1_only": (
            set(unsupported_operator_names) == EXPECTED_DEFAULT_LR_UNSUPPORTED_OPERATORS
        ),
        "b4_value_matches_etm2013_table1": _positive_finite(b4_mu_had)
        and math.isclose(
            float(b4_mu_had),
            EXPECTED_DEFAULT_LR_HADRONIC_B4_VALUE,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "b5_value_matches_etm2013_table1": _positive_finite(b5_mu_had)
        and math.isclose(
            float(b5_mu_had),
            EXPECTED_DEFAULT_LR_HADRONIC_B5_VALUE,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "r_chi_matches_frozen_default": _positive_finite(r_chi_mu_had)
        and math.isclose(
            float(r_chi_mu_had),
            EXPECTED_KAON_LR_R_CHI_EXACT_VALUE,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "q4_formula_id_is_exact": q4_formula_id == EXPECTED_DEFAULT_LR_Q4_FORMULA_ID,
        "q5_formula_id_is_exact": q5_formula_id == EXPECTED_DEFAULT_LR_Q5_FORMULA_ID,
        "b4_source_id_is_exact": b4_source_id == EXPECTED_DEFAULT_LR_HADRONIC_B4_SOURCE_ID,
        "b5_source_id_is_exact": b5_source_id == EXPECTED_DEFAULT_LR_HADRONIC_B5_SOURCE_ID,
        "b4_source_transformation_id_is_none": b4_source_transformation_id == "none",
        "b5_source_transformation_id_is_none": b5_source_transformation_id == "none",
        "r_chi_source_id_is_exact": (
            r_chi_source_id == EXPECTED_DEFAULT_LR_HADRONIC_R_CHI_SOURCE_ID
        ),
        "mass_source_id_is_exact": mass_source_id == EXPECTED_DEFAULT_LR_HADRONIC_MASS_SOURCE_ID,
        "decay_constant_source_id_is_exact": (
            decay_constant_source_id == EXPECTED_DEFAULT_LR_HADRONIC_DECAY_CONSTANT_SOURCE_ID
        ),
        "contract_operator_basis_id_is_exact": (
            contract_basis_id == EXPECTED_DEFAULT_LR_OPERATOR_BASIS_ID
        ),
        "contract_operator_normalization_id_is_exact": (
            contract_normalization_id == EXPECTED_DEFAULT_LR_OPERATOR_NORMALIZATION_ID
        ),
        "contract_hamiltonian_convention_id_is_exact": (
            contract_hamiltonian_convention_id
            == EXPECTED_DEFAULT_LR_HAMILTONIAN_CONVENTION_ID
        ),
        "b4_source_scheme_matches_bundle": bool(scheme_id)
        and b4_source_scheme_id == scheme_id,
        "b5_source_scheme_matches_bundle": bool(scheme_id)
        and b5_source_scheme_id == scheme_id,
        "r_chi_source_scheme_matches_frozen_mass_scheme": (
            r_chi_source_scheme_id == EXPECTED_KAON_LR_R_CHI_MASS_SCHEME_ID
        ),
        "r_chi_source_scheme_stays_distinct_from_bundle_scheme": bool(scheme_id)
        and bool(r_chi_source_scheme_id)
        and r_chi_source_scheme_id != scheme_id,
        "b4_source_scale_matches_bundle": _positive_finite(mu_had_GeV)
        and _positive_finite(b4_source_scale_GeV)
        and math.isclose(
            float(b4_source_scale_GeV),
            float(mu_had_GeV),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ),
        "b5_source_scale_matches_bundle": _positive_finite(mu_had_GeV)
        and _positive_finite(b5_source_scale_GeV)
        and math.isclose(
            float(b5_source_scale_GeV),
            float(mu_had_GeV),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ),
        "r_chi_source_scale_matches_bundle": _positive_finite(mu_had_GeV)
        and _positive_finite(r_chi_source_scale_GeV)
        and math.isclose(
            float(r_chi_source_scale_GeV),
            float(mu_had_GeV),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ),
        "b4_source_mentions_etm2013_table1_ms_2gev": (
            EXPECTED_DEFAULT_LR_HADRONIC_ETM_CITATION in b4_source_metadata_text
            and "Table 1" in b4_source_metadata_text
            and "Buras" in b4_source_metadata_text
            and "2 GeV" in b4_source_metadata_text
            and "B4" in b4_source_metadata_text
            and f"{EXPECTED_DEFAULT_LR_HADRONIC_B4_VALUE:.2f}" in b4_source_metadata_text
        ),
        "b5_source_mentions_etm2013_table1_ms_2gev": (
            EXPECTED_DEFAULT_LR_HADRONIC_ETM_CITATION in b5_source_metadata_text
            and "Table 1" in b5_source_metadata_text
            and "Buras" in b5_source_metadata_text
            and "2 GeV" in b5_source_metadata_text
            and "B5" in b5_source_metadata_text
            and f"{EXPECTED_DEFAULT_LR_HADRONIC_B5_VALUE:.2f}" in b5_source_metadata_text
        ),
        "q4_matches_bv2004_eq5": expected_q4 is not None
        and q4_matrix_element_GeV4 is not None
        and math.isclose(q4_matrix_element_GeV4, expected_q4, rel_tol=0.0, abs_tol=1.0e-15),
        "q5_matches_bv2004_eq5": expected_q5 is not None
        and q5_matrix_element_GeV4 is not None
        and math.isclose(q5_matrix_element_GeV4, expected_q5, rel_tol=0.0, abs_tol=1.0e-15),
        "custom_lr_only_rejects_default_lr_bundle": custom_lr_only_rejects_default_lr_bundle,
        "custom_combined_rejects_default_lr_bundle": (
            custom_combined_rejects_default_lr_bundle
        ),
        "custom_lr_surfaces_still_require_custom_inputs": (
            custom_lr_only_rejects_default_lr_bundle
            and custom_combined_rejects_default_lr_bundle
        ),
        **default_lr_rejection_checks,
    }
    failures = [
        f"default LR hadronic acceptance check failed: {name}"
        for name, ok in checks.items()
        if not ok
    ]
    if not export_matches_builder_default:
        failures.append(
            "default LR hadronic acceptance check failed: builder payload drifts from export"
        )
    return failures, {
        "status": "ok",
        "deterministic": checks["probe_payload_is_deterministic"],
        "cross_process_deterministic": checks["probe_payload_is_cross_process_deterministic"],
        "payload_sha256": payload_sha256,
        "same_process_repeat_payload_sha256": second_payload_sha256,
        "cross_process_payload_sha256": cross_process_payload_sha256,
        "system_id": system_id,
        "scheme_id": scheme_id,
        "mu_had_GeV": float(mu_had_GeV) if mu_had_GeV is not None else None,
        "bundle_id": bundle_id,
        "source_id": source_id,
        "provenance_ids": (
            list(provenance_ids)
            if isinstance(provenance_ids, Sequence)
            and not isinstance(provenance_ids, (str, bytes, bytearray))
            else None
        ),
        "input_policy_id": input_policy_id,
        "B4_mu_had": float(b4_mu_had) if b4_mu_had is not None else None,
        "B5_mu_had": float(b5_mu_had) if b5_mu_had is not None else None,
        "R_chi_mu_had": float(r_chi_mu_had) if r_chi_mu_had is not None else None,
        "m_K0_GeV": float(m_K0_GeV) if m_K0_GeV is not None else None,
        "f_K_GeV": float(f_K_GeV) if f_K_GeV is not None else None,
        "q4_matrix_element_GeV4": q4_matrix_element_GeV4,
        "q5_matrix_element_GeV4": q5_matrix_element_GeV4,
        "b4_source_id": b4_source_id,
        "b5_source_id": b5_source_id,
        "r_chi_source_id": r_chi_source_id,
        "mass_source_id": mass_source_id,
        "decay_constant_source_id": decay_constant_source_id,
        "b4_source_scheme_id": b4_source_scheme_id,
        "b5_source_scheme_id": b5_source_scheme_id,
        "r_chi_source_scheme_id": r_chi_source_scheme_id,
        "b4_source_scale_GeV": (
            float(b4_source_scale_GeV) if b4_source_scale_GeV is not None else None
        ),
        "b5_source_scale_GeV": (
            float(b5_source_scale_GeV) if b5_source_scale_GeV is not None else None
        ),
        "r_chi_source_scale_GeV": (
            float(r_chi_source_scale_GeV) if r_chi_source_scale_GeV is not None else None
        ),
        "b4_source_metadata_text": b4_source_metadata_text,
        "b5_source_metadata_text": b5_source_metadata_text,
        "custom_lr_only_rejection_error": custom_lr_only_rejection_error,
        "custom_combined_rejection_error": custom_combined_rejection_error,
        "default_lr_rejection_errors": default_lr_rejection_errors,
        "checks": checks,
        "payload": payload,
    }


def _build_lr_r_chi_freeze_probe_summary(
    modules: Mapping[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    hadronic_module = modules.get("eft_deltaf2.hadronic")
    if hadronic_module is None:
        return [], {"status": "missing"}

    try:
        payload = _default_kaon_lr_r_chi_freeze_payload(modules)
        second_payload = _default_kaon_lr_r_chi_freeze_payload(modules)
        freeze_value = _default_kaon_lr_r_chi_freeze_value(modules)
        freeze_object_payload = _payload_from_value(freeze_value.as_dict())
        cross_process_payload = _cross_process_lr_r_chi_freeze_payload()
    except AssertionError as exc:
        return ([str(exc)], {"status": "missing_lr_r_chi_freeze_helper"})
    except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
        return ([f"{type(exc).__name__}: {exc}"], {"status": "error", "error": str(exc)})

    payload_canonical_json, payload_sha256 = _canonical_payload_json_and_sha256(payload)
    second_payload_canonical_json, second_payload_sha256 = (
        _canonical_payload_json_and_sha256(second_payload)
    )
    cross_process_payload_canonical_json, cross_process_payload_sha256 = (
        _canonical_payload_json_and_sha256(cross_process_payload)
    )
    system_id = _nested_value(payload, (("system_id",),))
    object_schema_id = _nested_value(freeze_object_payload, (("schema_id",),))
    summary_schema_id = _nested_value(payload, (("schema_id",),))
    mu_had_GeV = _nested_value(payload, (("mu_had_GeV",),))
    operator_scheme_id = _nested_value(
        payload,
        (
            ("operator_scheme_id",),
            ("operator_renormalization_scheme_id",),
        ),
    )
    mass_scheme_id = _nested_value(
        payload,
        (
            ("mass_scheme_id",),
            ("mass_renormalization_scheme_id",),
        ),
    )
    active_flavor_policy_id = _nested_value(
        payload,
        (
            ("mass_active_flavor_policy_id",),
            ("active_flavor_policy_id",),
        ),
    )
    input_provenance_mode_id = _nested_value(payload, (("input_provenance_mode_id",),))
    input_policy_id = _nested_value(payload, (("input_policy_id",),))
    freeze_id = _nested_value(payload, (("freeze_id",),))
    source_id = _nested_value(payload, (("source_id",),))
    provenance_ids = _nested_value(payload, (("provenance_ids",),))
    derivation_formula_id = _nested_value(payload, (("derivation_formula_id",),))
    derivation_formula_source_id = _nested_value(
        payload,
        (("derivation_formula_source_id",),),
    )
    no_hidden_conversion_policy_id = _nested_value(
        payload,
        (("no_hidden_conversion_policy_id",),),
    )
    m_K0_GeV = _nested_value(payload, (("m_K0_GeV",),))
    m_s_mu_had_GeV = _nested_value(payload, (("m_s_mu_had_GeV",),))
    m_d_mu_had_GeV = _nested_value(payload, (("m_d_mu_had_GeV",),))
    r_chi_mu_had = _nested_value(payload, (("R_chi_mu_had",),))
    kaon_mass_source_id = _nested_value(payload, (("kaon_mass_source", "source_id"),))
    strange_mass_source_id = _nested_value(
        payload,
        (("strange_mass_source", "source_id"),),
    )
    down_mass_source_id = _nested_value(payload, (("down_mass_source", "source_id"),))
    strange_mass_source_scale_GeV = _nested_value(
        payload,
        (("strange_mass_source", "scale_GeV"),),
    )
    down_mass_source_scale_GeV = _nested_value(
        payload,
        (("down_mass_source", "scale_GeV"),),
    )

    expected_r_chi = None
    if all(
        _positive_finite(value)
        for value in (m_K0_GeV, m_s_mu_had_GeV, m_d_mu_had_GeV, r_chi_mu_had)
    ):
        expected_r_chi = (
            float(m_K0_GeV) / (float(m_s_mu_had_GeV) + float(m_d_mu_had_GeV))
        ) ** 2

    default_lr_hadronic_available = any(
        callable(getattr(hadronic_module, export_name, None))
        for export_name in DEFAULT_LR_HADRONIC_EXPORT_NAMES
    )

    custom_lr_only_rejection_error = None
    custom_combined_rejection_error = None
    custom_lr_only_rejects_r_chi_freeze = False
    custom_combined_rejects_r_chi_freeze = False
    observable_module = modules.get("eft_deltaf2.observables")
    if observable_module is not None:
        lr_evaluator = _get_callable(observable_module, CUSTOM_LR_OBSERVABLE_EXPORT_NAMES)
        if callable(lr_evaluator):
            try:
                probe_rg_value = _replace_rg_lr_wilsons(
                    _default_rg_value(modules),
                    q4_lr=LR_OBSERVABLE_PROBE_Q4,
                    q5_lr=LR_OBSERVABLE_PROBE_Q5,
                    renormalization_scheme_id=str(operator_scheme_id),
                    matching_scale_GeV=float(mu_had_GeV),
                )
                _invoke_observable_evaluation(
                    lr_evaluator,
                    rg_value=probe_rg_value,
                    hadronic_value=freeze_value,
                )
            except ValueError as exc:
                custom_lr_only_rejection_error = f"{type(exc).__name__}: {exc}"
                custom_lr_only_rejects_r_chi_freeze = (
                    custom_lr_only_rejection_error
                    == "ValueError: hadronic_inputs must be a Paper07101869KaonLRHadronicInputs"
                )
            except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
                custom_lr_only_rejection_error = f"{type(exc).__name__}: {exc}"

        combined_evaluator = _get_callable(
            observable_module,
            CUSTOM_TOTAL_OBSERVABLE_EXPORT_NAMES,
        )
        if callable(combined_evaluator):
            try:
                custom_q1_hadronic = _build_custom_q1_hadronic_bundle(
                    modules,
                    mu_had_GeV=float(mu_had_GeV),
                )
                combined_rg_value = _replace_rg_lr_wilsons(
                    _default_rg_value(modules),
                    q4_lr=LR_OBSERVABLE_PROBE_Q4,
                    q5_lr=LR_OBSERVABLE_PROBE_Q5,
                    renormalization_scheme_id=str(operator_scheme_id),
                    matching_scale_GeV=float(mu_had_GeV),
                )
                _invoke_custom_total_observable_evaluation(
                    combined_evaluator,
                    rg_value=combined_rg_value,
                    q1_hadronic_bundle=custom_q1_hadronic,
                    lr_hadronic_inputs=freeze_value,
                )
            except ValueError as exc:
                custom_combined_rejection_error = f"{type(exc).__name__}: {exc}"
                custom_combined_rejects_r_chi_freeze = (
                    custom_combined_rejection_error
                    == (
                        "ValueError: lr_hadronic_inputs must be a "
                        "Paper07101869KaonLRHadronicInputs"
                    )
                )
            except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
                custom_combined_rejection_error = f"{type(exc).__name__}: {exc}"

    checks = {
        "probe_payload_is_deterministic": payload_canonical_json == second_payload_canonical_json,
        "probe_payload_is_cross_process_deterministic": (
            payload_canonical_json == cross_process_payload_canonical_json
        ),
        "probe_payload_sha256_matches_repeat": payload_sha256 == second_payload_sha256,
        "probe_payload_sha256_matches_cross_process": (
            payload_sha256 == cross_process_payload_sha256
        ),
        "object_schema_is_exact": (
            object_schema_id == EXPECTED_KAON_LR_R_CHI_FREEZE_SCHEMA_ID
        ),
        "summary_schema_is_exact": (
            summary_schema_id == EXPECTED_KAON_LR_R_CHI_SUMMARY_SCHEMA_ID
        ),
        "object_and_summary_schema_are_distinct": (
            object_schema_id is not None
            and summary_schema_id is not None
            and object_schema_id != summary_schema_id
        ),
        "system_is_kaon": _is_kaon_system(system_id),
        "mass_scale_is_2gev": _positive_finite(mu_had_GeV)
        and math.isclose(float(mu_had_GeV), 2.0, rel_tol=0.0, abs_tol=1.0e-12),
        "formula_matches_bv2004": expected_r_chi is not None
        and math.isclose(
            float(r_chi_mu_had),
            expected_r_chi,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "mass_sources_frozen": isinstance(provenance_ids, Sequence)
        and not isinstance(provenance_ids, (str, bytes, bytearray))
        and list(provenance_ids)
        == [source_id, kaon_mass_source_id, strange_mass_source_id, down_mass_source_id],
        "mass_sources_at_declared_scale": all(
            _positive_finite(value)
            and math.isclose(float(value), float(mu_had_GeV), rel_tol=0.0, abs_tol=1.0e-12)
            for value in (strange_mass_source_scale_GeV, down_mass_source_scale_GeV)
        ),
        "mass_scheme_is_explicit": bool(mass_scheme_id),
        "mass_scheme_id_is_exact": (
            mass_scheme_id == EXPECTED_KAON_LR_R_CHI_MASS_SCHEME_ID
        ),
        "active_flavor_policy_is_frozen": (
            active_flavor_policy_id
            == EXPECTED_KAON_LR_R_CHI_ACTIVE_FLAVOR_POLICY_ID
        ),
        "active_flavor_policy_is_explicit_n_l_4": (
            active_flavor_policy_id
            == EXPECTED_KAON_LR_R_CHI_ACTIVE_FLAVOR_POLICY_ID
        ),
        "derivation_formula_frozen": bool(derivation_formula_id)
        and bool(derivation_formula_source_id),
        "provenance_mode_is_derived_default": bool(input_provenance_mode_id)
        and "derived" in str(input_provenance_mode_id).lower(),
        "input_policy_is_r_chi_freeze_only": bool(input_policy_id)
        and "freeze_only" in str(input_policy_id).lower(),
        "no_hidden_conversion_policy": (
            no_hidden_conversion_policy_id
            == EXPECTED_KAON_LR_R_CHI_NO_HIDDEN_CONVERSION_POLICY_ID
        ),
        "no_hidden_conversion_policy_is_explicit": (
            no_hidden_conversion_policy_id
            == EXPECTED_KAON_LR_R_CHI_NO_HIDDEN_CONVERSION_POLICY_ID
        ),
        "m_s_mu_had_matches_exact_freeze": _positive_finite(m_s_mu_had_GeV)
        and math.isclose(
            float(m_s_mu_had_GeV),
            EXPECTED_KAON_LR_R_CHI_M_S_2GEV_GEV,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "m_d_mu_had_matches_exact_freeze": _positive_finite(m_d_mu_had_GeV)
        and math.isclose(
            float(m_d_mu_had_GeV),
            EXPECTED_KAON_LR_R_CHI_M_D_2GEV_GEV,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "m_k0_gev_matches_exact_freeze": _positive_finite(m_K0_GeV)
        and math.isclose(
            float(m_K0_GeV),
            EXPECTED_KAON_LR_R_CHI_M_K0_GEV,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "r_chi_mu_had_matches_exact_freeze": _positive_finite(r_chi_mu_had)
        and math.isclose(
            float(r_chi_mu_had),
            EXPECTED_KAON_LR_R_CHI_EXACT_VALUE,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "operator_vs_mass_scheme_not_aliased": bool(operator_scheme_id)
        and bool(mass_scheme_id)
        and operator_scheme_id != mass_scheme_id,
        "custom_lr_only_rejects_r_chi_freeze": custom_lr_only_rejects_r_chi_freeze,
        "custom_combined_rejects_r_chi_freeze": custom_combined_rejects_r_chi_freeze,
        "custom_lr_only_rejection_message_is_explicit_type_guard": (
            custom_lr_only_rejection_error
            == "ValueError: hadronic_inputs must be a Paper07101869KaonLRHadronicInputs"
        ),
        "custom_combined_rejection_message_is_explicit_type_guard": (
            custom_combined_rejection_error
            == (
                "ValueError: lr_hadronic_inputs must be a "
                "Paper07101869KaonLRHadronicInputs"
            )
        ),
        "custom_lr_surfaces_still_require_custom_inputs": (
            custom_lr_only_rejects_r_chi_freeze and custom_combined_rejects_r_chi_freeze
        ),
    }
    failures = [
        f"LR R_chi freeze acceptance check failed: {name}"
        for name, ok in checks.items()
        if not ok
    ]
    return failures, {
        "status": "ok",
        "deterministic": checks["probe_payload_is_deterministic"],
        "cross_process_deterministic": checks["probe_payload_is_cross_process_deterministic"],
        "payload_sha256": payload_sha256,
        "same_process_repeat_payload_sha256": second_payload_sha256,
        "cross_process_payload_sha256": cross_process_payload_sha256,
        "object_schema_id": object_schema_id,
        "summary_schema_id": summary_schema_id,
        "system_id": system_id,
        "mu_had_GeV": float(mu_had_GeV) if mu_had_GeV is not None else None,
        "freeze_id": freeze_id,
        "source_id": source_id,
        "provenance_ids": list(provenance_ids) if isinstance(provenance_ids, Sequence) else None,
        "input_provenance_mode_id": input_provenance_mode_id,
        "input_policy_id": input_policy_id,
        "operator_renormalization_scheme_id": operator_scheme_id,
        "operator_scheme_id": operator_scheme_id,
        "mass_renormalization_scheme_id": mass_scheme_id,
        "mass_scheme_id": mass_scheme_id,
        "mass_active_flavor_policy_id": active_flavor_policy_id,
        "derivation_formula_id": derivation_formula_id,
        "derivation_formula_source_id": derivation_formula_source_id,
        "no_hidden_conversion_policy_id": no_hidden_conversion_policy_id,
        "m_K0_GeV": float(m_K0_GeV) if m_K0_GeV is not None else None,
        "m_s_mu_had_GeV": float(m_s_mu_had_GeV) if m_s_mu_had_GeV is not None else None,
        "m_d_mu_had_GeV": float(m_d_mu_had_GeV) if m_d_mu_had_GeV is not None else None,
        "R_chi_mu_had": float(r_chi_mu_had) if r_chi_mu_had is not None else None,
        "kaon_mass_source_id": kaon_mass_source_id,
        "strange_mass_source_id": strange_mass_source_id,
        "down_mass_source_id": down_mass_source_id,
        "strange_mass_source_scale_GeV": (
            float(strange_mass_source_scale_GeV)
            if strange_mass_source_scale_GeV is not None
            else None
        ),
        "down_mass_source_scale_GeV": (
            float(down_mass_source_scale_GeV)
            if down_mass_source_scale_GeV is not None
            else None
        ),
        "default_lr_hadronic_available": default_lr_hadronic_available,
        "default_lr_hadronic_exports_present": {
            export_name: callable(getattr(hadronic_module, export_name, None))
            for export_name in DEFAULT_LR_HADRONIC_EXPORT_NAMES
        },
        "custom_lr_only_rejection_error": custom_lr_only_rejection_error,
        "custom_combined_rejection_error": custom_combined_rejection_error,
        "checks": checks,
        "payload": payload,
    }


def _build_custom_lr_hadronic_probe_summary(
    modules: Mapping[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    module = modules.get("eft_deltaf2.hadronic")
    if module is None:
        return [], {"status": "missing"}

    try:
        default_payload = _default_hadronic_payload(modules)
        default_bundle = _default_hadronic_value(modules)
    except AssertionError as exc:
        return ([str(exc)], {"status": "missing_default_hadronic_export"})

    default_lr_fields_absent = all(
        not _mapping_contains_key_deep(default_payload, field_name)
        for field_name in DEFAULT_HADRONIC_FORBIDDEN_LR_FIELDS
    )
    default_supported_operator_names = _string_list(
        _nested_value(
            default_payload,
            (
                ("supported_operator_names",),
                ("contract", "supported_operator_names"),
            ),
        )
    )
    default_unsupported_operator_names = _string_list(
        _nested_value(
            default_payload,
            (
                ("unsupported_operator_names",),
                ("contract", "unsupported_operator_names"),
            ),
        )
    )
    default_bundle_unchanged = (
        set(default_supported_operator_names) == EXPECTED_SUPPORTED_OPERATORS
        and set(default_unsupported_operator_names) == EXPECTED_UNSUPPORTED_OPERATORS
        and default_lr_fields_absent
    )

    try:
        builder = _custom_lr_hadronic_builder(module)
        declared_scheme_id = LR_HADRONIC_PROBE_SCHEME_ID
        declared_mu_had_GeV = LR_HADRONIC_PROBE_MU_HAD_GEV
        m_K0_GeV = float(default_bundle.m_K0_GeV)
        f_K_GeV = float(default_bundle.f_K_GeV)

        b4_source = _build_hadronic_source_ref(
            module,
            source_id=LR_HADRONIC_PROBE_B4_SOURCE_ID,
            citation="Becirevic and Villadoro, hep-lat/0408029",
            locator_label="eq. (5), O4 scalar LR matrix element",
            year=2004,
            scheme_id=declared_scheme_id,
            scale_GeV=declared_mu_had_GeV,
            notes="Custom LR hadronic probe source for B4(mu_had).",
        )
        b5_source = _build_hadronic_source_ref(
            module,
            source_id=LR_HADRONIC_PROBE_B5_SOURCE_ID,
            citation="Becirevic and Villadoro, hep-lat/0408029",
            locator_label="eq. (5), O5 scalar LR matrix element",
            year=2004,
            scheme_id=declared_scheme_id,
            scale_GeV=declared_mu_had_GeV,
            notes="Custom LR hadronic probe source for B5(mu_had).",
        )
        r_chi_source = _build_hadronic_source_ref(
            module,
            source_id=LR_HADRONIC_PROBE_R_CHI_SOURCE_ID,
            citation="Custom caller-supplied chiral ratio input",
            locator_label="R_chi(mu_had) external input",
            year=2026,
            scheme_id=declared_scheme_id,
            scale_GeV=declared_mu_had_GeV,
            notes="Custom LR hadronic probe source for R_chi(mu_had).",
        )
        built = _invoke_custom_lr_hadronic_builder(
            builder,
            default_bundle=default_bundle,
            B4_mu_had=LR_HADRONIC_PROBE_B4,
            B5_mu_had=LR_HADRONIC_PROBE_B5,
            R_chi_mu_had=LR_HADRONIC_PROBE_R_CHI,
            B4_source=b4_source,
            B5_source=b5_source,
            R_chi_source=r_chi_source,
            renormalization_scheme_id=declared_scheme_id,
            mu_had_GeV=declared_mu_had_GeV,
        )
        payload = _payload_from_value(built)
    except AssertionError as exc:
        return (
            [str(exc)],
            {
                "status": "missing_custom_lr_hadronic_builder",
                "default_bundle_unchanged": default_bundle_unchanged,
                "default_lr_fields_absent": default_lr_fields_absent,
            },
        )
    except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
        return (
            [f"{type(exc).__name__}: {exc}"],
            {
                "status": "error",
                "default_bundle_unchanged": default_bundle_unchanged,
                "default_lr_fields_absent": default_lr_fields_absent,
                "error": f"{type(exc).__name__}: {exc}",
            },
        )

    scheme_id = _nested_value(
        payload,
        (
            ("scheme_id",),
            ("renormalization_scheme_id",),
            ("tags", "renormalization_scheme_id"),
            ("contract", "renormalization_scheme_id"),
        ),
    )
    mu_had_GeV = _nested_value(
        payload,
        (
            ("mu_had_GeV",),
            ("tags", "mu_had_GeV"),
            ("contract", "mu_had_GeV"),
        ),
    )
    input_provenance_mode_id = _nested_value(
        payload,
        (
            ("input_provenance_mode_id",),
            ("tags", "input_provenance_mode_id"),
        ),
    )
    input_policy_id = _nested_value(payload, (("input_policy_id",),))
    contract_input_policy_id = _nested_value(payload, (("contract", "input_policy_id"),))
    bundle_notes = _nested_value(payload, (("notes",),))
    contract_notes = _nested_value(payload, (("contract", "notes"),))
    contract_id = _nested_value(
        payload,
        (
            ("lr_hadronic_contract_id",),
            ("contract", "lr_hadronic_contract_id"),
            ("contract", "contract_id"),
            ("contract", "schema_id"),
            ("contract", "input_policy_id"),
            ("contract_id",),
        ),
    )
    q4_formula_id = _nested_value(
        payload,
        (
            ("q4_matrix_element_formula_id",),
            ("formula_ids", "Q4_LR"),
            ("formula_ids", "q4_lr"),
            ("contract", "q4_matrix_element_formula_id"),
            ("contract", "formula_ids", "Q4_LR"),
            ("contract", "formula_id"),
            ("matrix_element_formula_id",),
        ),
    )
    q5_formula_id = _nested_value(
        payload,
        (
            ("q5_matrix_element_formula_id",),
            ("formula_ids", "Q5_LR"),
            ("formula_ids", "q5_lr"),
            ("contract", "q5_matrix_element_formula_id"),
            ("contract", "formula_ids", "Q5_LR"),
            ("contract", "formula_id"),
            ("matrix_element_formula_id",),
        ),
    )
    q4_matrix_element_GeV4 = _extract_lr_matrix_element(payload, "Q4_LR")
    q5_matrix_element_GeV4 = _extract_lr_matrix_element(payload, "Q5_LR")
    expected_q4 = (
        2.0
        * LR_HADRONIC_PROBE_R_CHI
        * (m_K0_GeV**2)
        * (f_K_GeV**2)
        * LR_HADRONIC_PROBE_B4
    )
    expected_q5 = (
        (2.0 / 3.0)
        * LR_HADRONIC_PROBE_R_CHI
        * (m_K0_GeV**2)
        * (f_K_GeV**2)
        * LR_HADRONIC_PROBE_B5
    )

    default_q1_only_observable_surface_rejects_custom_lr_bundle = False
    observable_error = None
    observable_module = modules.get("eft_deltaf2.observables")
    if observable_module is not None:
        try:
            evaluator = _get_callable(
                observable_module,
                PARAMETERIZED_OBSERVABLE_EXPORT_NAMES,
            )
            if not callable(evaluator):
                raise AssertionError(
                    "observable module exposes no callable evaluator for LR-HAD acceptance probe"
                )
            probed_rg_value = _replace_rg_lr_wilsons(
                _default_rg_value(modules),
                q4_lr=complex(1.0e-12, 0.0),
                q5_lr=complex(-5.0e-13, 0.0),
            )
            _invoke_observable_evaluation(
                evaluator,
                rg_value=probed_rg_value,
                hadronic_value=built,
            )
        except ValueError as exc:
            observable_error = str(exc)
            default_q1_only_observable_surface_rejects_custom_lr_bundle = (
                (
                    "Q4_LR" in observable_error
                    and "Q5_LR" in observable_error
                    and "blocked" in observable_error
                )
                or "hadronic_bundle must be a Paper07101869KaonHadronicBundle"
                in observable_error
            )
        except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
            observable_error = f"{type(exc).__name__}: {exc}"

    checks = {
        "default_bundle_unchanged": default_bundle_unchanged,
        "default_lr_fields_absent": default_lr_fields_absent,
        "custom_mode_only": bool(input_provenance_mode_id)
        and "custom" in str(input_provenance_mode_id).lower(),
        "input_policy_id_present": bool(input_policy_id),
        "contract_input_policy_id_present": bool(contract_input_policy_id),
        "input_policy_id_matches_contract": input_policy_id == contract_input_policy_id,
        "input_policy_id_is_current_custom_frozen_default_bundle_no_auto_consumption": (
            _is_current_custom_lr_input_policy(input_policy_id)
            and _is_current_custom_lr_input_policy(contract_input_policy_id)
        ),
        "bundle_notes_match_current_custom_lr_contract": _is_current_custom_lr_bundle_notes(
            bundle_notes
        ),
        "contract_notes_match_current_custom_lr_contract": (
            _is_current_custom_lr_contract_notes(contract_notes)
        ),
        "contract_id_present": bool(contract_id),
        "q4_formula_id_present": bool(q4_formula_id),
        "q5_formula_id_present": bool(q5_formula_id),
        "declared_scheme_matches_payload": scheme_id == declared_scheme_id,
        "declared_mu_had_matches_payload": _positive_finite(mu_had_GeV)
        and math.isclose(
            float(mu_had_GeV),
            declared_mu_had_GeV,
            rel_tol=0.0,
            abs_tol=1e-12,
        ),
        "source_scheme_matches_declared_scheme": all(
            bool(source is not None)
            and str(source.renormalization_scheme_id) == declared_scheme_id
            for source in (b4_source, b5_source, r_chi_source)
        ),
        "source_scale_matches_declared_mu_had": all(
            math.isclose(
                float(source.scale_GeV),
                declared_mu_had_GeV,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            for source in (b4_source, b5_source, r_chi_source)
        ),
        "q4_matches_bv2004_eq5": q4_matrix_element_GeV4 is not None
        and math.isclose(q4_matrix_element_GeV4, expected_q4, rel_tol=0.0, abs_tol=1e-15),
        "q5_matches_bv2004_eq5": q5_matrix_element_GeV4 is not None
        and math.isclose(q5_matrix_element_GeV4, expected_q5, rel_tol=0.0, abs_tol=1e-15),
        "default_q1_only_observable_surface_rejects_custom_lr_bundle": (
            default_q1_only_observable_surface_rejects_custom_lr_bundle
        ),
    }
    failures = [
        f"custom LR hadronic acceptance check failed: {name}"
        for name, ok in checks.items()
        if not ok
    ]
    return failures, {
        "status": "ok",
        "default_bundle_unchanged": default_bundle_unchanged,
        "default_lr_fields_absent": default_lr_fields_absent,
        "input_provenance_mode_id": input_provenance_mode_id,
        "input_policy_id": input_policy_id,
        "contract_input_policy_id": contract_input_policy_id,
        "bundle_notes": bundle_notes,
        "contract_notes": contract_notes,
        "contract_id": contract_id,
        "q4_formula_id": q4_formula_id,
        "q5_formula_id": q5_formula_id,
        "declared_scheme_id": declared_scheme_id,
        "declared_mu_had_GeV": declared_mu_had_GeV,
        "scheme_id": scheme_id,
        "mu_had_GeV": float(mu_had_GeV) if mu_had_GeV is not None else None,
        "m_K0_GeV": m_K0_GeV,
        "f_K_GeV": f_K_GeV,
        "B4_mu_had": LR_HADRONIC_PROBE_B4,
        "B5_mu_had": LR_HADRONIC_PROBE_B5,
        "R_chi_mu_had": LR_HADRONIC_PROBE_R_CHI,
        "b4_source_id": getattr(b4_source, "source_id", None),
        "b5_source_id": getattr(b5_source, "source_id", None),
        "r_chi_source_id": getattr(r_chi_source, "source_id", None),
        "b4_source_scheme_id": getattr(b4_source, "renormalization_scheme_id", None),
        "b5_source_scheme_id": getattr(b5_source, "renormalization_scheme_id", None),
        "r_chi_source_scheme_id": getattr(r_chi_source, "renormalization_scheme_id", None),
        "b4_source_scale_GeV": float(getattr(b4_source, "scale_GeV", 0.0)),
        "b5_source_scale_GeV": float(getattr(b5_source, "scale_GeV", 0.0)),
        "r_chi_source_scale_GeV": float(getattr(r_chi_source, "scale_GeV", 0.0)),
        "q4_matrix_element_GeV4": q4_matrix_element_GeV4,
        "q5_matrix_element_GeV4": q5_matrix_element_GeV4,
        "observable_error": observable_error,
        "checks": checks,
        "payload": payload,
    }


def _build_custom_lr_observable_probe_summary(
    modules: Mapping[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    observable_module = modules.get("eft_deltaf2.observables")
    hadronic_module = modules.get("eft_deltaf2.hadronic")
    if observable_module is None or hadronic_module is None:
        return [], {"status": "missing"}

    try:
        evaluator = _get_callable(observable_module, CUSTOM_LR_OBSERVABLE_EXPORT_NAMES)
        if not callable(evaluator):
            raise AssertionError(
                "observable module exposes no callable custom LR-only evaluator; expected one "
                "of " + ", ".join(CUSTOM_LR_OBSERVABLE_EXPORT_NAMES)
            )
        default_bundle = _default_hadronic_value(modules)
        builder = _custom_lr_hadronic_builder(hadronic_module)
        declared_scheme_id = LR_HADRONIC_PROBE_SCHEME_ID
        declared_mu_had_GeV = LR_HADRONIC_PROBE_MU_HAD_GEV
        b4_source = _build_hadronic_source_ref(
            hadronic_module,
            source_id=LR_HADRONIC_PROBE_B4_SOURCE_ID,
            citation="Becirevic and Villadoro, hep-lat/0408029",
            locator_label="eq. (5), O4 scalar LR matrix element",
            year=2004,
            scheme_id=declared_scheme_id,
            scale_GeV=declared_mu_had_GeV,
            notes="Custom LR observable probe source for B4(mu_had).",
        )
        b5_source = _build_hadronic_source_ref(
            hadronic_module,
            source_id=LR_HADRONIC_PROBE_B5_SOURCE_ID,
            citation="Becirevic and Villadoro, hep-lat/0408029",
            locator_label="eq. (5), O5 scalar LR matrix element",
            year=2004,
            scheme_id=declared_scheme_id,
            scale_GeV=declared_mu_had_GeV,
            notes="Custom LR observable probe source for B5(mu_had).",
        )
        r_chi_source = _build_hadronic_source_ref(
            hadronic_module,
            source_id=LR_HADRONIC_PROBE_R_CHI_SOURCE_ID,
            citation="Custom caller-supplied chiral ratio input",
            locator_label="R_chi(mu_had) external input",
            year=2026,
            scheme_id=declared_scheme_id,
            scale_GeV=declared_mu_had_GeV,
            notes="Custom LR observable probe source for R_chi(mu_had).",
        )
        custom_hadronic = _invoke_custom_lr_hadronic_builder(
            builder,
            default_bundle=default_bundle,
            B4_mu_had=LR_HADRONIC_PROBE_B4,
            B5_mu_had=LR_HADRONIC_PROBE_B5,
            R_chi_mu_had=LR_HADRONIC_PROBE_R_CHI,
            B4_source=b4_source,
            B5_source=b5_source,
            R_chi_source=r_chi_source,
            renormalization_scheme_id=declared_scheme_id,
            mu_had_GeV=declared_mu_had_GeV,
        )
        default_rg_wilsons = _default_rg_value(modules).wilsons
        probe_rg_value = dataclasses.replace(
            default_rg_wilsons,
            contract=dataclasses.replace(
                default_rg_wilsons.contract,
                renormalization_scheme_id=declared_scheme_id,
            ),
            q1_vll=complex(0.0, 0.0),
            q1_vrr=complex(0.0, 0.0),
            q4_lr=LR_OBSERVABLE_PROBE_Q4,
            q5_lr=LR_OBSERVABLE_PROBE_Q5,
            matching_scale_GeV=declared_mu_had_GeV,
        )
        result = _invoke_observable_evaluation(
            evaluator,
            rg_value=probe_rg_value,
            hadronic_value=custom_hadronic,
        )
        payload = _payload_from_value(result)
    except AssertionError as exc:
        return ([str(exc)], {"status": "missing_custom_lr_observable_evaluator"})
    except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
        return ([f"{type(exc).__name__}: {exc}"], {"status": "error", "error": str(exc)})

    scheme_id = _nested_value(
        payload,
        (
            ("scheme_id",),
            ("renormalization_scheme_id",),
            ("tags", "renormalization_scheme_id"),
        ),
    )
    mu_had_GeV = _nested_value(
        payload,
        (
            ("mu_had_GeV",),
            ("tags", "mu_had_GeV"),
        ),
    )
    m12_payload = _nested_value(
        payload,
        (
            ("M12_K_LR_NP_GeV",),
            ("observables", "M12_K_LR_NP"),
        ),
    )
    delta_m_value = _nested_value(
        payload,
        (
            ("delta_m_K_LR_NP_GeV",),
            ("observables", "Delta_m_K_LR_NP"),
        ),
    )
    q4_matrix_element_GeV4 = _extract_lr_matrix_element(payload, "Q4_LR")
    q5_matrix_element_GeV4 = _extract_lr_matrix_element(payload, "Q5_LR")
    expected_m12 = (
        LR_OBSERVABLE_PROBE_Q4 * float(custom_hadronic.q4_matrix_element_GeV4)
        + LR_OBSERVABLE_PROBE_Q5 * float(custom_hadronic.q5_matrix_element_GeV4)
    ) / (2.0 * float(custom_hadronic.m_K0_GeV))
    expected_delta_m = 2.0 * float(expected_m12.real)
    observed_m12 = None
    if isinstance(m12_payload, Mapping):
        try:
            observed_m12 = complex(
                float(m12_payload["real"]),
                float(m12_payload["imag"]),
            )
        except (KeyError, TypeError, ValueError):
            observed_m12 = None
    observed_delta_m = None
    if delta_m_value is not None:
        try:
            observed_delta_m = float(delta_m_value)
        except (TypeError, ValueError):
            observed_delta_m = None

    checks = {
        "declared_scheme_matches_payload": scheme_id == declared_scheme_id,
        "declared_mu_had_matches_payload": _positive_finite(mu_had_GeV)
        and math.isclose(
            float(mu_had_GeV),
            declared_mu_had_GeV,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ),
        "m12_matches_hand_calculation": observed_m12 is not None
        and math.isclose(observed_m12.real, expected_m12.real, rel_tol=0.0, abs_tol=1.0e-24)
        and math.isclose(observed_m12.imag, expected_m12.imag, rel_tol=0.0, abs_tol=1.0e-24),
        "delta_m_matches_hand_calculation": observed_delta_m is not None
        and math.isclose(observed_delta_m, expected_delta_m, rel_tol=0.0, abs_tol=1.0e-24),
        "q4_matrix_element_present": q4_matrix_element_GeV4 is not None,
        "q5_matrix_element_present": q5_matrix_element_GeV4 is not None,
    }
    failures = [
        f"custom LR observable acceptance check failed: {name}"
        for name, ok in checks.items()
        if not ok
    ]
    return failures, {
        "status": "ok",
        "declared_scheme_id": declared_scheme_id,
        "declared_mu_had_GeV": declared_mu_had_GeV,
        "scheme_id": scheme_id,
        "mu_had_GeV": float(mu_had_GeV) if mu_had_GeV is not None else None,
        "probe_q4_lr": _complex_as_dict(LR_OBSERVABLE_PROBE_Q4),
        "probe_q5_lr": _complex_as_dict(LR_OBSERVABLE_PROBE_Q5),
        "q4_matrix_element_GeV4": q4_matrix_element_GeV4,
        "q5_matrix_element_GeV4": q5_matrix_element_GeV4,
        "M12_K_LR_NP_GeV": (
            _complex_as_dict(observed_m12) if observed_m12 is not None else None
        ),
        "delta_m_K_LR_NP_GeV": observed_delta_m,
        "checks": checks,
        "payload": payload,
    }


def _build_custom_combined_observable_probe_summary(
    modules: Mapping[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    observable_module = modules.get("eft_deltaf2.observables")
    hadronic_module = modules.get("eft_deltaf2.hadronic")
    if observable_module is None or hadronic_module is None:
        return [], {"status": "missing"}

    try:
        combined_evaluator = _get_callable(
            observable_module,
            CUSTOM_TOTAL_OBSERVABLE_EXPORT_NAMES,
        )
        if not callable(combined_evaluator):
            raise AssertionError(
                "observable module exposes no callable custom combined evaluator; expected one "
                "of " + ", ".join(CUSTOM_TOTAL_OBSERVABLE_EXPORT_NAMES)
            )
        q1_evaluator = _get_callable(observable_module, PARAMETERIZED_OBSERVABLE_EXPORT_NAMES)
        if not callable(q1_evaluator):
            raise AssertionError(
                "observable module exposes no callable Q1 observable evaluator for "
                "combined-surface benchmarking"
            )
        lr_evaluator = _get_callable(observable_module, CUSTOM_LR_OBSERVABLE_EXPORT_NAMES)
        if not callable(lr_evaluator):
            raise AssertionError(
                "observable module exposes no callable LR-only evaluator for "
                "combined-surface benchmarking"
            )
        custom_q1_hadronic = _build_custom_q1_hadronic_bundle(modules)
        default_bundle = _default_hadronic_value(modules)
        lr_builder = _custom_lr_hadronic_builder(hadronic_module)
        source_prefix = "hadronic.kaon.custom_total_probe"
        custom_lr_hadronic = _invoke_custom_lr_hadronic_builder(
            lr_builder,
            default_bundle=default_bundle,
            B4_mu_had=LR_HADRONIC_PROBE_B4,
            B5_mu_had=LR_HADRONIC_PROBE_B5,
            R_chi_mu_had=LR_HADRONIC_PROBE_R_CHI,
            B4_source=_build_hadronic_source_ref(
                hadronic_module,
                source_id=f"{source_prefix}.b4.v1",
                citation="Becirevic and Villadoro, hep-lat/0408029",
                locator_label="eq. (5), O4 scalar LR matrix element",
                year=2004,
                scheme_id=str(custom_q1_hadronic.renormalization_scheme_id),
                scale_GeV=float(custom_q1_hadronic.mu_had_GeV),
                notes="Custom combined observable probe source for B4(mu_had).",
            ),
            B5_source=_build_hadronic_source_ref(
                hadronic_module,
                source_id=f"{source_prefix}.b5.v1",
                citation="Becirevic and Villadoro, hep-lat/0408029",
                locator_label="eq. (5), O5 scalar LR matrix element",
                year=2004,
                scheme_id=str(custom_q1_hadronic.renormalization_scheme_id),
                scale_GeV=float(custom_q1_hadronic.mu_had_GeV),
                notes="Custom combined observable probe source for B5(mu_had).",
            ),
            R_chi_source=_build_hadronic_source_ref(
                hadronic_module,
                source_id=f"{source_prefix}.rchi.v1",
                citation="Custom caller-supplied chiral ratio input",
                locator_label="R_chi(mu_had) external input",
                year=2026,
                scheme_id=str(custom_q1_hadronic.renormalization_scheme_id),
                scale_GeV=float(custom_q1_hadronic.mu_had_GeV),
                notes="Custom combined observable probe source for R_chi(mu_had).",
            ),
            renormalization_scheme_id=str(custom_q1_hadronic.renormalization_scheme_id),
            mu_had_GeV=float(custom_q1_hadronic.mu_had_GeV),
            source_prefix=source_prefix,
        )
        probe_rg_value = _replace_rg_lr_wilsons(
            _default_rg_value(modules),
            q4_lr=LR_OBSERVABLE_PROBE_Q4,
            q5_lr=LR_OBSERVABLE_PROBE_Q5,
            renormalization_scheme_id=str(custom_q1_hadronic.renormalization_scheme_id),
            matching_scale_GeV=float(custom_q1_hadronic.mu_had_GeV),
        )
        result = _invoke_custom_total_observable_evaluation(
            combined_evaluator,
            rg_value=probe_rg_value,
            q1_hadronic_bundle=custom_q1_hadronic,
            lr_hadronic_inputs=custom_lr_hadronic,
        )
        payload = _payload_from_value(result)
        q1_payload = _payload_from_value(
            _invoke_observable_evaluation(
                q1_evaluator,
                rg_value=_replace_rg_lr_wilsons(
                    probe_rg_value,
                    q1_vll=probe_rg_value.wilsons.q1_vll,
                    q1_vrr=probe_rg_value.wilsons.q1_vrr,
                    q4_lr=0.0 + 0.0j,
                    q5_lr=0.0 + 0.0j,
                ),
                hadronic_value=custom_q1_hadronic,
            )
        )
        lr_payload = _payload_from_value(
            _invoke_observable_evaluation(
                lr_evaluator,
                rg_value=_replace_rg_lr_wilsons(
                    probe_rg_value,
                    q1_vll=0.0 + 0.0j,
                    q1_vrr=0.0 + 0.0j,
                    q4_lr=LR_OBSERVABLE_PROBE_Q4,
                    q5_lr=LR_OBSERVABLE_PROBE_Q5,
                ),
                hadronic_value=custom_lr_hadronic,
            )
        )
    except AssertionError as exc:
        return ([str(exc)], {"status": "missing_custom_combined_observable_evaluator"})
    except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
        return ([f"{type(exc).__name__}: {exc}"], {"status": "error", "error": str(exc)})

    declared_scheme_id = str(custom_q1_hadronic.renormalization_scheme_id)
    declared_mu_had_GeV = float(custom_q1_hadronic.mu_had_GeV)
    q1_matrix_element = float(custom_q1_hadronic.q1_matrix_element_GeV4)
    q4_matrix_element = float(custom_lr_hadronic.q4_matrix_element_GeV4)
    q5_matrix_element = float(custom_lr_hadronic.q5_matrix_element_GeV4)
    hand_q1 = (
        (probe_rg_value.wilsons.q1_vll + probe_rg_value.wilsons.q1_vrr)
        * q1_matrix_element
        / (2.0 * float(custom_q1_hadronic.m_K0_GeV))
    )
    hand_lr = (
        LR_OBSERVABLE_PROBE_Q4 * q4_matrix_element
        + LR_OBSERVABLE_PROBE_Q5 * q5_matrix_element
    ) / (2.0 * float(custom_lr_hadronic.m_K0_GeV))
    hand_total = hand_q1 + hand_lr
    hand_delta_m = 2.0 * float(hand_total.real)
    observed_q1 = _complex_from_payload(payload["M12_K_NP_Q1_GeV"])
    observed_lr = _complex_from_payload(payload["M12_K_LR_NP_GeV"])
    observed_total = _complex_from_payload(payload["M12_K_NP_TOTAL_GeV"])
    observed_delta_m = float(payload["delta_m_K_NP_TOTAL_GeV"])
    existing_q1 = _complex_from_payload(q1_payload["M12_K_NP_GeV"])
    existing_lr = _complex_from_payload(lr_payload["M12_K_LR_NP_GeV"])
    checks = {
        "scope_id_is_custom_total": (
            payload.get("observable_scope_id") == EXPECTED_CUSTOM_TOTAL_SCOPE_ID
        ),
        "interpretation_is_custom_total": (
            payload.get("interpretation") == EXPECTED_CUSTOM_TOTAL_INTERPRETATION_ID
        ),
        "observable_ids_are_custom_total": (
            payload.get("m12_observable_id") == EXPECTED_CUSTOM_TOTAL_M12_OBSERVABLE_ID
            and payload.get("delta_m_observable_id")
            == EXPECTED_CUSTOM_TOTAL_DELTA_M_OBSERVABLE_ID
        ),
        "declared_scheme_matches_payload": (
            payload.get("renormalization_scheme_id") == declared_scheme_id
        ),
        "declared_mu_had_matches_payload": math.isclose(
            float(payload["mu_had_GeV"]),
            declared_mu_had_GeV,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ),
        "q1_component_matches_hand_calculation": (
            math.isclose(observed_q1.real, hand_q1.real, rel_tol=0.0, abs_tol=1.0e-24)
            and math.isclose(observed_q1.imag, hand_q1.imag, rel_tol=0.0, abs_tol=1.0e-24)
        ),
        "lr_component_matches_hand_calculation": (
            math.isclose(observed_lr.real, hand_lr.real, rel_tol=0.0, abs_tol=1.0e-24)
            and math.isclose(observed_lr.imag, hand_lr.imag, rel_tol=0.0, abs_tol=1.0e-24)
        ),
        "total_matches_hand_calculation": (
            math.isclose(observed_total.real, hand_total.real, rel_tol=0.0, abs_tol=1.0e-24)
            and math.isclose(observed_total.imag, hand_total.imag, rel_tol=0.0, abs_tol=1.0e-24)
        ),
        "delta_m_matches_hand_calculation": math.isclose(
            observed_delta_m,
            hand_delta_m,
            rel_tol=0.0,
            abs_tol=1.0e-24,
        ),
        "total_equals_component_sum": abs(observed_total - (observed_q1 + observed_lr))
        <= 1.0e-24,
        "total_matches_existing_q1_plus_lr_surfaces": (
            abs(observed_total - (existing_q1 + existing_lr)) <= 1.0e-24
        ),
        "custom_q1_bundle_active": "custom"
        in str(custom_q1_hadronic.input_provenance_mode_id).lower(),
        "epsilon_k_blocked": "epsilon_K_NP" not in payload.get("observables", {}),
    }
    failures = [
        f"custom combined observable acceptance check failed: {name}"
        for name, ok in checks.items()
        if not ok
    ]
    return failures, {
        "status": "ok",
        "declared_scheme_id": declared_scheme_id,
        "declared_mu_had_GeV": declared_mu_had_GeV,
        "scheme_id": payload.get("renormalization_scheme_id"),
        "mu_had_GeV": float(payload["mu_had_GeV"]),
        "probe_q1_vll": _complex_as_dict(probe_rg_value.wilsons.q1_vll),
        "probe_q1_vrr": _complex_as_dict(probe_rg_value.wilsons.q1_vrr),
        "probe_q4_lr": _complex_as_dict(LR_OBSERVABLE_PROBE_Q4),
        "probe_q5_lr": _complex_as_dict(LR_OBSERVABLE_PROBE_Q5),
        "q1_matrix_element_GeV4": q1_matrix_element,
        "q4_matrix_element_GeV4": q4_matrix_element,
        "q5_matrix_element_GeV4": q5_matrix_element,
        "M12_K_NP_Q1_GeV": _complex_as_dict(observed_q1),
        "M12_K_LR_NP_GeV": _complex_as_dict(observed_lr),
        "M12_K_NP_TOTAL_GeV": _complex_as_dict(observed_total),
        "delta_m_K_NP_TOTAL_GeV": observed_delta_m,
        "checks": checks,
        "payload": payload,
    }


def _build_custom_b_observable_probe_payload(
    modules: Mapping[str, Any],
    system_id: str,
) -> tuple[dict[str, Any], Any, Any]:
    observable_module = modules.get("eft_deltaf2.observables")
    assert observable_module is not None, "eft_deltaf2.observables module is required"

    evaluator = _get_custom_b_q1_observable_evaluator(observable_module, system_id)
    rg_result = _build_custom_b_rg_result(modules, system_id)
    declared_scheme_id = str(rg_result.wilsons.renormalization_scheme_id)
    declared_mu_had_GeV = float(rg_result.wilsons.matching_scale_GeV)
    hadronic_bundle = _build_custom_b_hadronic_bundle(
        modules,
        system_id,
        renormalization_scheme_id=declared_scheme_id,
        mu_had_GeV=declared_mu_had_GeV,
    )
    result = _invoke_observable_evaluation(
        evaluator,
        rg_value=rg_result,
        hadronic_value=hadronic_bundle,
    )
    return _payload_from_value(result), rg_result, hadronic_bundle


def _build_custom_b_observable_probe_summary(
    modules: Mapping[str, Any],
    system_id: str,
) -> tuple[list[str], dict[str, Any]]:
    observable_module = modules.get("eft_deltaf2.observables")
    if observable_module is None:
        return [], {"status": "missing"}

    try:
        payload, rg_result, hadronic_bundle = _build_custom_b_observable_probe_payload(
            modules,
            system_id,
        )
        declared_scheme_id = str(rg_result.wilsons.renormalization_scheme_id)
        declared_mu_had_GeV = float(rg_result.wilsons.matching_scale_GeV)
        second_payload, _, _ = _build_custom_b_observable_probe_payload(
            modules,
            system_id,
        )
    except AssertionError as exc:
        return ([str(exc)], {"status": f"missing_custom_{system_id.lower()}_observable_evaluator"})
    except Exception as exc:  # pragma: no cover - surfaced in acceptance payload
        return ([f"{type(exc).__name__}: {exc}"], {"status": "error", "error": str(exc)})

    config = B_CUSTOM_Q1_PROBE_CONFIG[system_id]
    payload_canonical_json, payload_sha256 = _canonical_payload_json_and_sha256(payload)
    second_payload_canonical_json, same_process_repeat_payload_sha256 = (
        _canonical_payload_json_and_sha256(second_payload)
    )
    deterministic = payload_canonical_json == second_payload_canonical_json
    cross_process_payload = _cross_process_custom_b_observable_probe_payload(system_id)
    cross_process_payload_canonical_json, cross_process_payload_sha256 = (
        _canonical_payload_json_and_sha256(cross_process_payload)
    )
    cross_process_deterministic = payload_canonical_json == cross_process_payload_canonical_json
    expected_m12 = (
        (rg_result.wilsons.q1_vll + rg_result.wilsons.q1_vrr)
        * float(hadronic_bundle.q1_matrix_element_GeV4)
        / (2.0 * float(hadronic_bundle.meson_mass_GeV))
    )
    expected_delta_m = 2.0 * float(expected_m12.real)
    observed_m12 = _complex_from_payload(payload["M12_NP_GeV"])
    observed_delta_m = float(payload["delta_m_NP_GeV"])
    checks = {
        "system_id_matches_probe": payload.get("system_id") == system_id,
        "scope_id_is_custom_b_surface": payload.get("observable_scope_id")
        == config["scope_id"],
        "interpretation_is_custom_b_surface": payload.get("interpretation")
        == config["interpretation_id"],
        "observable_ids_are_custom_b_surface": (
            payload.get("m12_observable_id") == config["m12_id"]
            and payload.get("delta_m_observable_id") == config["delta_m_id"]
        ),
        "declared_scheme_matches_payload": (
            payload.get("renormalization_scheme_id") == declared_scheme_id
        ),
        "declared_mu_had_matches_payload": math.isclose(
            float(payload["mu_had_GeV"]),
            declared_mu_had_GeV,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ),
        "matching_scale_matches_payload": math.isclose(
            float(payload["matching_scale_GeV"]),
            declared_mu_had_GeV,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ),
        "m12_matches_hand_calculation": (
            math.isclose(observed_m12.real, expected_m12.real, rel_tol=0.0, abs_tol=1.0e-24)
            and math.isclose(observed_m12.imag, expected_m12.imag, rel_tol=0.0, abs_tol=1.0e-24)
        ),
        "delta_m_matches_hand_calculation": math.isclose(
            observed_delta_m,
            expected_delta_m,
            rel_tol=0.0,
            abs_tol=1.0e-24,
        ),
        "probe_payload_is_deterministic": deterministic,
        "probe_payload_is_cross_process_deterministic": cross_process_deterministic,
        "probe_payload_sha256_matches_repeat": (
            payload_sha256 == same_process_repeat_payload_sha256
        ),
        "probe_payload_sha256_matches_cross_process": (
            payload_sha256 == cross_process_payload_sha256
        ),
        "q1_matrix_element_present": _positive_finite(payload.get("q1_matrix_element_GeV4")),
        "lr_coefficients_remain_zero": (
            payload.get("coefficients", {}).get("Q4_LR", {}).get("real") == 0.0
            and payload.get("coefficients", {}).get("Q4_LR", {}).get("imag") == 0.0
            and payload.get("coefficients", {}).get("Q5_LR", {}).get("real") == 0.0
            and payload.get("coefficients", {}).get("Q5_LR", {}).get("imag") == 0.0
        ),
    }
    if bool(config.get("require_nonzero_probe", False)):
        checks["probe_payload_is_nonzero"] = abs(observed_m12) > 0.0 and abs(observed_delta_m) > 0.0
    failures = [
        f"custom {system_id} observable acceptance check failed: {name}"
        for name, ok in checks.items()
        if not ok
    ]
    return failures, {
        "status": "ok",
        "deterministic": deterministic,
        "cross_process_deterministic": cross_process_deterministic,
        "system_id": system_id,
        "declared_scheme_id": declared_scheme_id,
        "declared_mu_had_GeV": declared_mu_had_GeV,
        "scheme_id": payload.get("renormalization_scheme_id"),
        "mu_had_GeV": float(payload["mu_had_GeV"]),
        "matching_scale_GeV": float(payload["matching_scale_GeV"]),
        "probe_q1_vll": _complex_as_dict(rg_result.wilsons.q1_vll),
        "probe_q1_vrr": _complex_as_dict(rg_result.wilsons.q1_vrr),
        "payload_sha256": payload_sha256,
        "same_process_repeat_payload_sha256": same_process_repeat_payload_sha256,
        "cross_process_payload_sha256": cross_process_payload_sha256,
        "q1_matrix_element_GeV4": float(payload["q1_matrix_element_GeV4"]),
        "meson_mass_GeV": float(hadronic_bundle.meson_mass_GeV),
        "M12_NP_GeV": _complex_as_dict(observed_m12),
        "delta_m_NP_GeV": observed_delta_m,
        "checks": checks,
        "payload": payload,
    }


def _build_observable_summary(modules: Mapping[str, Any]) -> tuple[list[str], dict[str, Any]]:
    module = modules.get("eft_deltaf2.observables")
    if module is None:
        return [], {"status": "missing"}

    try:
        payload = _default_observable_payload(modules)
    except AssertionError as exc:
        return ([str(exc)], {"status": "missing_default_observable_export"})

    second = _default_observable_payload(modules)
    deterministic = json.dumps(payload, sort_keys=True) == json.dumps(second, sort_keys=True)
    failures = []
    if not deterministic:
        failures.append("default observable export is not deterministic within one process")
    cross_process_payload = _cross_process_default_observable_payload()
    cross_process_deterministic = json.dumps(payload, sort_keys=True) == json.dumps(
        cross_process_payload,
        sort_keys=True,
    )
    if not cross_process_deterministic:
        failures.append("default observable export is not deterministic across processes")

    interpretation = _nested_value(
        payload,
        (
            ("interpretation",),
            ("observable_scope_id",),
            ("contract", "interpretation"),
            ("contract", "observable_scope_id"),
            ("tags", "interpretation"),
            ("tags", "observable_scope_id"),
        ),
    )
    observable_scope_id = _nested_value(
        payload,
        (
            ("observable_scope_id",),
            ("contract", "observable_scope_id"),
            ("tags", "observable_scope_id"),
        ),
    )
    m12_observable_id = _nested_value(
        payload,
        (
            ("m12_observable_id",),
            ("tags", "m12_observable_id"),
        ),
    )
    delta_m_observable_id = _nested_value(
        payload,
        (
            ("delta_m_observable_id",),
            ("tags", "delta_m_observable_id"),
        ),
    )
    system_id = _nested_value(
        payload,
        (
            ("system_id",),
            ("system",),
            ("tags", "system_id"),
            ("contract", "system_id"),
        ),
    )
    scheme_id = _nested_value(
        payload,
        (
            ("renormalization_scheme_id",),
            ("scheme_id",),
            ("tags", "renormalization_scheme_id"),
            ("contract", "renormalization_scheme_id"),
        ),
    )
    mu_had_GeV = _nested_value(
        payload,
        (
            ("mu_had_GeV",),
            ("evaluation_scale_GeV",),
            ("hadronic_scale_GeV",),
            ("tags", "mu_had_GeV"),
            ("tags", "evaluation_scale_GeV"),
            ("contract", "mu_had_GeV"),
        ),
    )
    parity_relation_id = _nested_value(
        payload,
        (
            ("hadronic_bundle", "parity_relation_id"),
            ("tags", "parity_relation_id"),
        ),
    )
    hadronic_mu_had_GeV = None
    hadronic_parity_relation_id = None
    if "eft_deltaf2.hadronic" in modules:
        hadronic_payload = _default_hadronic_payload(modules)
        hadronic_mu_had_GeV = _nested_value(
            hadronic_payload,
            (
                ("mu_had_GeV",),
                ("evaluation_scale_GeV",),
                ("hadronic_scale_GeV",),
                ("tags", "mu_had_GeV"),
                ("tags", "evaluation_scale_GeV"),
                ("contract", "mu_had_GeV"),
            ),
        )
        hadronic_parity_relation_id = _nested_value(
            hadronic_payload,
            (
                ("parity_relation_id",),
                ("tags", "parity_relation_id"),
            ),
        )
    try:
        observables = _extract_named_observables(payload)
    except AssertionError as exc:
        failures.append(str(exc))
        observables = {}

    checks = {
        "interpretation_is_exact_np_only": interpretation == EXPECTED_INTERPRETATION_ID,
        "observable_scope_is_pr5a_surface": observable_scope_id == EXPECTED_OBSERVABLE_SCOPE_ID,
        "m12_observable_id_is_supported": m12_observable_id == EXPECTED_M12_OBSERVABLE_ID,
        "delta_m_observable_id_is_supported": (
            delta_m_observable_id == EXPECTED_DELTA_M_OBSERVABLE_ID
        ),
        "system_is_kaon": _is_kaon_system(system_id),
        "scheme_id_present": bool(scheme_id),
        "mu_had_present_and_positive": _positive_finite(mu_had_GeV),
        "observable_names_match_supported_surface": set(observables) == SUPPORTED_OBSERVABLE_ROWS,
        "observable_values_finite": all(
            math.isfinite(float(item)) for item in observables.values()
        ),
        "m12_re_matches_frozen_reference": math.isclose(
            observables.get("M12_K_NP.re", float("nan")),
            EXPECTED_M12_K_NP["real"],
            rel_tol=0.0,
            abs_tol=1.0e-24,
        ),
        "m12_im_matches_frozen_reference": math.isclose(
            observables.get("M12_K_NP.im", float("nan")),
            EXPECTED_M12_K_NP["imag"],
            rel_tol=0.0,
            abs_tol=1.0e-24,
        ),
        "delta_m_matches_frozen_reference": math.isclose(
            observables.get("Delta_m_K_NP", float("nan")),
            EXPECTED_DELTA_M_K_NP_GEV,
            rel_tol=0.0,
            abs_tol=1.0e-24,
        ),
        "mu_had_matches_hadronic_bundle": (
            float(mu_had_GeV) == float(hadronic_mu_had_GeV)
            if mu_had_GeV is not None and hadronic_mu_had_GeV is not None
            else False
        ),
        "parity_relation_tag_present": bool(parity_relation_id),
        "parity_relation_matches_hadronic_bundle": (
            parity_relation_id == hadronic_parity_relation_id
            if parity_relation_id is not None and hadronic_parity_relation_id is not None
            else False
        ),
    }
    failures.extend(
        f"observable acceptance check failed: {name}" for name, ok in checks.items() if not ok
    )

    return failures, {
        "status": "ok" if not failures else "failed",
        "deterministic": deterministic,
        "cross_process_deterministic": cross_process_deterministic,
        "interpretation": interpretation,
        "observable_scope_id": observable_scope_id,
        "m12_observable_id": m12_observable_id,
        "delta_m_observable_id": delta_m_observable_id,
        "system_id": system_id,
        "scheme_id": scheme_id,
        "parity_relation_id": parity_relation_id,
        "mu_had_GeV": float(mu_had_GeV) if mu_had_GeV is not None else None,
        "hadronic_mu_had_GeV": (
            float(hadronic_mu_had_GeV) if hadronic_mu_had_GeV is not None else None
        ),
        "observable_names": sorted(observables),
        "checks": checks,
        "payload": payload,
    }


def _materialize_artifact_summary(
    modules: Mapping[str, Any],
    export_dir: Path,
) -> tuple[list[str], dict[str, Any]]:
    exported = _export_default_artifacts(modules, export_dir)
    artifact_module = modules["artifacts"]
    file_paths = exported["file_paths"]
    verifier_payload = _run_independent_verifier(
        wilson_path=file_paths["wilson"],
        hadronic_path=file_paths["hadronic"],
        observable_path=file_paths["observable"],
        provenance_path=file_paths["provenance"],
    )
    with tempfile.TemporaryDirectory(prefix="paper_0710_1869_writer_") as writer_tmpdir:
        writer_paths = artifact_module.write_default_paper_0710_1869_kaon_artifact_exports(
            writer_tmpdir
        )
        writer_file_sha256 = {
            "wilson": _sha256_digest(writer_paths.wilson_path),
            "hadronic": _sha256_digest(writer_paths.hadronic_path),
            "observable": _sha256_digest(writer_paths.observable_path),
            "provenance": _sha256_digest(writer_paths.provenance_path),
        }
    canonical_paths = artifact_module.default_paper_0710_1869_kaon_artifact_export_paths()
    canonical_file_paths = {
        "wilson": canonical_paths.wilson_path,
        "hadronic": canonical_paths.hadronic_path,
        "observable": canonical_paths.observable_path,
        "provenance": canonical_paths.provenance_path,
    }
    canonical_export_files_present = all(path.exists() for path in canonical_file_paths.values())
    canonical_file_sha256 = (
        {name: _sha256_digest(path) for name, path in canonical_file_paths.items()}
        if canonical_export_files_present
        else {}
    )
    hadronic_bundle = exported["hadronic_bundle"]
    observable_bundle = exported["observable_bundle"]
    wilson_bundle = exported["wilson_bundle"]
    wilson_scale_map = {scale.name: float(scale.value_gev) for scale in wilson_bundle.scales}

    observable_names = [record.name for record in observable_bundle.observables]
    observable_values = {
        record.name: float(record.value)
        for record in observable_bundle.observables
    }
    coefficient_summary = {
        record.operator: {
            "real": float(record.value.real),
            "imag": float(record.value.imag),
        }
        for record in wilson_bundle.coefficients
    }
    checks = {
        "exported_files_present": all(path.exists() for path in file_paths.values()),
        "schema_ok": bool(verifier_payload.get("schema_ok")),
        "import_isolation_ok": bool(verifier_payload.get("import_isolation_ok")),
        "numeric_match_ok": bool(verifier_payload.get("numeric_match_ok")),
        "scope_ok": bool(verifier_payload.get("scope_ok")),
        "point_id_present": bool(verifier_payload.get("point_id")),
        "bundle_ids_present": all(
            bool(bundle_id) for bundle_id in verifier_payload.get("bundle_ids", {}).values()
        ),
        "observable_source_wilson_bundle_linked": (
            observable_bundle.source_wilson_bundle_id == wilson_bundle.metadata.bundle_id
        ),
        "observable_source_hadronic_bundle_linked": (
            observable_bundle.source_hadronic_bundle_id == hadronic_bundle.metadata.bundle_id
        ),
        "hadronic_mu_had_matches_wilson_scale": (
            wilson_scale_map.get("mu_had") == float(hadronic_bundle.mu_had_GeV)
        ),
        "hadronic_operator_normalization_matches_wilson": (
            hadronic_bundle.operator_normalization == wilson_bundle.operator_normalization
        ),
        "observable_operator_normalization_matches_wilson": (
            observable_bundle.operator_normalization == wilson_bundle.operator_normalization
        ),
        "observable_rows_match_honest_scope": set(observable_names)
        == EXPECTED_ARTIFACT_OBSERVABLE_ROWS,
        "observable_values_match_frozen_reference": (
            math.isclose(
                observable_values.get("M12_K_NP.re", float("nan")),
                EXPECTED_M12_K_NP["real"],
                rel_tol=0.0,
                abs_tol=1.0e-24,
            )
            and math.isclose(
                observable_values.get("M12_K_NP.im", float("nan")),
                EXPECTED_M12_K_NP["imag"],
                rel_tol=0.0,
                abs_tol=1.0e-24,
            )
            and math.isclose(
                observable_values.get("Delta_m_K_NP", float("nan")),
                EXPECTED_DELTA_M_K_NP_GEV,
                rel_tol=0.0,
                abs_tol=1.0e-24,
            )
        ),
        "epsilon_k_not_exported": "epsilon_K_NP" not in observable_names,
        "lr_coefficients_absent_or_zero": all(
            name not in coefficient_summary
            or (
                coefficient_summary[name]["real"] == 0.0
                and coefficient_summary[name]["imag"] == 0.0
            )
            for name in sorted(EXPECTED_UNSUPPORTED_OPERATORS)
        ),
        "canonical_export_files_present": canonical_export_files_present,
        "tracked_default_exports_match_current_export": (
            canonical_export_files_present and exported["file_sha256"] == canonical_file_sha256
        ),
        "writer_outputs_are_deterministic": exported["file_sha256"] == writer_file_sha256,
    }
    failures = [
        f"artifact acceptance check failed: {name}"
        for name, ok in checks.items()
        if not ok
    ]
    failures.extend(
        f"artifact verifier issue: {code}"
        for code in verifier_payload.get("issue_codes", [])
    )
    return failures, {
        "status": "ok" if not failures else "failed",
        "point_id": verifier_payload.get("point_id"),
        "bundle_ids": verifier_payload.get("bundle_ids", {}),
        "tolerance_policy": verifier_payload.get("tolerance_policy", {}),
        "tolerances": verifier_payload.get("tolerances", {}),
        "schema_ok": bool(verifier_payload.get("schema_ok")),
        "import_isolation_ok": bool(verifier_payload.get("import_isolation_ok")),
        "numeric_match_ok": bool(verifier_payload.get("numeric_match_ok")),
        "scope_ok": bool(verifier_payload.get("scope_ok")),
        "observable_rows": observable_names,
        "observable_values": observable_values,
        "coefficients": coefficient_summary,
        "hadronic_bundle_id": hadronic_bundle.metadata.bundle_id,
        "hadronic_schema_id": hadronic_bundle.metadata.schema_name,
        "file_sha256": exported["file_sha256"],
        "canonical_file_sha256": canonical_file_sha256,
        "writer_file_sha256": writer_file_sha256,
        "issue_codes": verifier_payload.get("issue_codes", []),
        "issues": verifier_payload.get("issues", []),
        "unexpected_import_targets": verifier_payload.get("unexpected_import_targets", []),
        "loaded_forbidden_modules": verifier_payload.get("loaded_forbidden_modules", []),
        "numeric_checks": verifier_payload.get("numeric_checks", {}),
        "numeric_diffs": verifier_payload.get("numeric_diffs", {}),
        "recomputed": verifier_payload.get("recomputed", {}),
        "load_error": verifier_payload.get("load_error"),
        "checks": checks,
    }


def _build_artifact_summary(
    modules: Mapping[str, Any],
    export_dir: Path | None = None,
) -> tuple[list[str], dict[str, Any]]:
    required_modules = {
        "artifacts",
        "verifier",
        "eft_deltaf2.hadronic",
        "eft_deltaf2.matching_kkgluon",
        "eft_deltaf2.rg",
        "eft_deltaf2.observables",
    }
    missing = sorted(name for name in required_modules if name not in modules)
    if missing:
        return [], {"status": "missing", "missing_modules": missing}

    if export_dir is not None:
        return _materialize_artifact_summary(modules, export_dir)

    with tempfile.TemporaryDirectory(prefix="paper_0710_1869_artifacts_") as tmpdir:
        return _materialize_artifact_summary(modules, Path(tmpdir))


def _forbidden_modules_loaded() -> list[str]:
    loaded = []
    for name in sorted(sys.modules):
        if name in FORBIDDEN_MODULES:
            loaded.append(name)
            continue
        if name == "deltaf2" or any(
            name.startswith(f"{forbidden}.") for forbidden in FORBIDDEN_MODULES
        ):
            loaded.append(name)
    return loaded


def _build_convention_summary(module: Any) -> tuple[list[str], dict[str, Any]]:
    factory = getattr(module, "default_paper_0710_1869_conventions", None)
    if not callable(factory):
        return ["conventions.default_paper_0710_1869_conventions is missing"], {}

    conventions = factory()
    payload = _canonicalize(conventions)
    checks = {
        "mode_id": payload.get("mode_id") == EXPECTED_MODE_ID,
        "observable_scope": payload.get("observable_scope") == "np_only",
        "rg_order": payload.get("rg_order") == "lo",
        "kk_gluon_normalization_id": "mu_gs"
        in str(payload.get("kk_gluon_normalization_id", "")),
        "provenance_policy_id": "deterministic_bundle.required"
        in str(payload.get("provenance_policy_id", "")),
        "verifier_policy_id": "independent_verifier.required"
        in str(payload.get("verifier_policy_id", "")),
    }
    failures = [
        f"conventions contract check failed: {name}"
        for name, ok in checks.items()
        if not ok
    ]
    return failures, {"payload": payload, "checks": checks}


def _build_scale_summary(module: Any) -> tuple[list[str], dict[str, Any]]:
    factory = getattr(module, "default_paper_0710_1869_scales", None)
    if not callable(factory):
        return ["scales.default_paper_0710_1869_scales is missing"], {}

    scale_point = factory()
    payload = _canonicalize(scale_point)
    missing_fields = [name for name in REQUIRED_SCALE_FIELDS if name not in payload]
    expected_propagator_mass_GeV = (
        float(payload["m_g1_GeV"])
        if payload.get("m_KK_eff_GeV") is None
        else float(payload["m_KK_eff_GeV"])
    )
    checks = {
        "mode_id": payload.get("mode_id") == EXPECTED_MODE_ID,
        "required_fields_present": not missing_fields,
        "xi_g_matches_mass_ratio": float(payload["xi_g"])
        == float(payload["m_g1_GeV"]) / float(payload["Lambda_IR_GeV"]),
        "propagator_mass_matches_scale_contract": (
            float(scale_point.propagator_mass_GeV) == expected_propagator_mass_GeV
        ),
        "effective_kk_flag_matches_payload": (
            bool(scale_point.has_explicit_effective_kk_scale)
            == (payload.get("m_KK_eff_GeV") is not None)
        ),
    }
    failures = [
        f"scale contract check failed: {name}" for name, ok in checks.items() if not ok
    ]
    failures.extend(f"scale contract missing field: {field}" for field in missing_fields)
    summary = {
        "payload": payload,
        "checks": checks,
        "propagator_mass_GeV": float(scale_point.propagator_mass_GeV),
        "has_explicit_effective_kk_scale": bool(
            scale_point.has_explicit_effective_kk_scale
        ),
    }
    return failures, summary


def _structural_payload(modules: Mapping[str, Any]) -> dict[str, Any]:
    conventions = modules["conventions"].default_paper_0710_1869_conventions()
    scale_point = modules["scales"].default_paper_0710_1869_scales()
    request = modules["scan"].Paper07101869ScanRequest(
        conventions=conventions,
        scale_points=(scale_point,),
    )
    rows = modules["scan"].build_structural_scan_rows(request)
    return {
        "conventions": conventions,
        "request": request,
        "rows": rows,
    }


def _build_structural_summary(modules: Mapping[str, Any]) -> tuple[list[str], dict[str, Any]]:
    missing = [name for name in ("conventions", "scales", "scan") if name not in modules]
    if missing:
        missing_list = ", ".join(missing)
        return [f"structural benchmark is missing modules: {missing_list}"], {}

    first = _canonicalize(_structural_payload(modules))
    second = _canonicalize(_structural_payload(modules))
    deterministic = json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)

    rows = first.get("rows", [])
    request_scale_points = first.get("request", {}).get("scale_points", [])
    scale_by_label = {
        str(scale.get("label")): scale
        for scale in request_scale_points
        if isinstance(scale, Mapping) and scale.get("label") is not None
    }

    def _row_expected_propagator_mass_GeV(row: Mapping[str, Any]) -> float:
        point_id = str(row.get("point_id"))
        scale = scale_by_label.get(point_id, {})
        effective_scale = scale.get("m_KK_eff_GeV")
        if effective_scale is not None:
            return float(effective_scale)
        return float(scale.get("m_g1_GeV"))

    checks = {
        "deterministic": deterministic,
        "row_count_positive": len(rows) > 0,
        "row_mode_ids_match": all(row.get("mode_id") == EXPECTED_MODE_ID for row in rows),
        "row_statuses_are_structural_only": all(
            row.get("status") == "structural_only" for row in rows
        ),
        "propagator_mass_matches_scale_contract": all(
            float(row.get("propagator_mass_GeV")) == _row_expected_propagator_mass_GeV(row)
            for row in rows
        ),
    }
    failures = [
        f"structural benchmark check failed: {name}"
        for name, ok in checks.items()
        if not ok
    ]
    return failures, {"checks": checks, "payload": first}


def _build_pr1_benchmark_summary(modules: Mapping[str, Any]) -> tuple[list[str], dict[str, Any]]:
    missing = [name for name in ("benchmarks", "model") if name not in modules]
    if missing:
        missing_list = ", ".join(missing)
        return [], {"status": "missing", "missing_modules": missing_list}

    benchmark = modules["benchmarks"].default_paper_0710_1869_pr1_benchmark()
    eq3_summary = modules["model"].evaluate_default_paper_0710_1869_eq3_consistency()
    benchmark_payload = _canonicalize(benchmark)
    eq3_payload = _canonicalize(eq3_summary)
    checks = {
        "benchmark_status_is_structural": benchmark_payload.get("status")
        == "sourced_structural_only",
        "eq3_is_not_faked_exact": eq3_payload.get("exact_agreement") is False,
        "eq3_residual_is_finite": float(eq3_payload.get("max_abs_residual", float("nan"))) >= 0.0,
        "eq3_residual_is_bounded": float(eq3_payload.get("max_abs_residual", 1.0)) < 0.25,
    }
    failures = [f"PR1 benchmark check failed: {name}" for name, ok in checks.items() if not ok]
    return failures, {
        "status": "ok",
        "checks": checks,
        "benchmark": benchmark_payload,
        "eq3": eq3_payload,
    }


def _collect_summary(export_dir: Path | None = None) -> tuple[dict[str, Any], list[str]]:
    failures: list[str] = []

    static_failures, static_summary = _summarize_static_contracts(PAPER_PACKAGE_DIR)
    failures.extend(static_failures)

    imported_modules, import_summary, import_failures = _import_paper_modules()
    failures.extend(import_failures)

    leaked_modules = _forbidden_modules_loaded()
    if leaked_modules:
        failures.append(
            "paper acceptance imported repo-v1 modules: " + ", ".join(leaked_modules)
        )

    convention_failures: list[str] = []
    convention_summary: dict[str, Any] = {}
    if "conventions" in imported_modules:
        convention_failures, convention_summary = _build_convention_summary(
            imported_modules["conventions"]
        )
        failures.extend(convention_failures)

    scale_failures: list[str] = []
    scale_summary: dict[str, Any] = {}
    if "scales" in imported_modules:
        scale_failures, scale_summary = _build_scale_summary(imported_modules["scales"])
        failures.extend(scale_failures)

    structural_failures, structural_summary = _build_structural_summary(imported_modules)
    failures.extend(structural_failures)
    pr1_failures, pr1_summary = _build_pr1_benchmark_summary(imported_modules)
    failures.extend(pr1_failures)
    coupling_failures, coupling_summary = _build_coupling_summary(imported_modules)
    failures.extend(coupling_failures)
    matching_failures, matching_summary = _build_eft_matching_summary(imported_modules)
    failures.extend(matching_failures)
    rg_failures, rg_summary = _build_eft_rg_lo_summary(imported_modules)
    failures.extend(rg_failures)
    lr_running_probe_failures, lr_running_probe_summary = _build_eft_lr_running_probe_summary(
        imported_modules
    )
    failures.extend(lr_running_probe_failures)
    lr_freeze_failures, lr_freeze_summary = _build_lr_contract_freeze_summary(imported_modules)
    failures.extend(lr_freeze_failures)
    hadronic_failures, hadronic_summary = _build_hadronic_summary(imported_modules)
    failures.extend(hadronic_failures)
    default_lr_hadronic_failures, default_lr_hadronic_summary = (
        _build_default_lr_hadronic_probe_summary(imported_modules)
    )
    failures.extend(default_lr_hadronic_failures)
    lr_r_chi_freeze_failures, lr_r_chi_freeze_summary = _build_lr_r_chi_freeze_probe_summary(
        imported_modules
    )
    failures.extend(lr_r_chi_freeze_failures)
    custom_lr_hadronic_failures, custom_lr_hadronic_summary = (
        _build_custom_lr_hadronic_probe_summary(imported_modules)
    )
    failures.extend(custom_lr_hadronic_failures)
    custom_lr_observable_failures, custom_lr_observable_summary = (
        _build_custom_lr_observable_probe_summary(imported_modules)
    )
    failures.extend(custom_lr_observable_failures)
    custom_combined_observable_failures, custom_combined_observable_summary = (
        _build_custom_combined_observable_probe_summary(imported_modules)
    )
    failures.extend(custom_combined_observable_failures)
    custom_bd_observable_failures, custom_bd_observable_summary = (
        _build_custom_b_observable_probe_summary(imported_modules, "B_d")
    )
    failures.extend(custom_bd_observable_failures)
    custom_bs_observable_failures, custom_bs_observable_summary = (
        _build_custom_b_observable_probe_summary(imported_modules, "B_s")
    )
    failures.extend(custom_bs_observable_failures)
    custom_d0_observable_failures, custom_d0_observable_summary = (
        _build_custom_b_observable_probe_summary(imported_modules, "D0")
    )
    failures.extend(custom_d0_observable_failures)
    observable_failures, observable_summary = _build_observable_summary(imported_modules)
    failures.extend(observable_failures)
    artifact_failures, artifact_summary = _build_artifact_summary(
        imported_modules,
        export_dir=export_dir,
    )
    failures.extend(artifact_failures)

    def _append_post_aggregation_failure(
        summary_name: str,
        check_name: str,
        passed: bool,
    ) -> None:
        if not passed:
            failures.append(f"{summary_name}: {check_name} is false")

    if custom_lr_observable_summary.get("status") == "ok":
        custom_checks = custom_lr_observable_summary.setdefault("checks", {})
        default_observables_unchanged = (
            observable_summary.get("status") == "ok"
            and set(observable_summary.get("observable_names", [])) == SUPPORTED_OBSERVABLE_ROWS
            and observable_summary.get("checks", {}).get("m12_re_matches_frozen_reference")
            is True
            and observable_summary.get("checks", {}).get("m12_im_matches_frozen_reference")
            is True
            and observable_summary.get("checks", {}).get("delta_m_matches_frozen_reference")
            is True
        )
        default_artifacts_unchanged = (
            artifact_summary.get("status") == "ok"
            and set(artifact_summary.get("observable_rows", []))
            == EXPECTED_ARTIFACT_OBSERVABLE_ROWS
            and artifact_summary.get("checks", {}).get("observable_values_match_frozen_reference")
            is True
            and artifact_summary.get("checks", {}).get(
                "tracked_default_exports_match_current_export"
            )
            is True
        )
        custom_lr_observable_summary["default_observables_unchanged"] = (
            default_observables_unchanged
        )
        custom_lr_observable_summary["default_artifacts_unchanged"] = (
            default_artifacts_unchanged
        )
        custom_checks["default_observables_unchanged"] = default_observables_unchanged
        custom_checks["default_artifacts_unchanged"] = default_artifacts_unchanged
        _append_post_aggregation_failure(
            "custom_lr_observable_probe",
            "default_observables_unchanged",
            default_observables_unchanged,
        )
        _append_post_aggregation_failure(
            "custom_lr_observable_probe",
            "default_artifacts_unchanged",
            default_artifacts_unchanged,
        )

    if default_lr_hadronic_summary.get("status") == "ok":
        custom_checks = default_lr_hadronic_summary.setdefault("checks", {})
        default_q1_hadronic_bundle_stays_q1_only = (
            hadronic_summary.get("status") == "ok"
            and hadronic_summary.get("checks", {}).get("supported_operator_subset_is_pr5a")
            is True
            and hadronic_summary.get("checks", {}).get("unsupported_operator_subset_is_lr_only")
            is True
        )
        default_observables_unchanged = (
            observable_summary.get("status") == "ok"
            and set(observable_summary.get("observable_names", [])) == SUPPORTED_OBSERVABLE_ROWS
            and observable_summary.get("checks", {}).get("m12_re_matches_frozen_reference")
            is True
            and observable_summary.get("checks", {}).get("m12_im_matches_frozen_reference")
            is True
            and observable_summary.get("checks", {}).get("delta_m_matches_frozen_reference")
            is True
        )
        default_artifacts_unchanged = (
            artifact_summary.get("status") == "ok"
            and set(artifact_summary.get("observable_rows", []))
            == EXPECTED_ARTIFACT_OBSERVABLE_ROWS
            and artifact_summary.get("checks", {}).get("observable_values_match_frozen_reference")
            is True
            and artifact_summary.get("checks", {}).get(
                "tracked_default_exports_match_current_export"
            )
            is True
        )
        default_lr_hadronic_summary["default_q1_hadronic_bundle_stays_q1_only"] = (
            default_q1_hadronic_bundle_stays_q1_only
        )
        default_lr_hadronic_summary["default_observables_unchanged"] = (
            default_observables_unchanged
        )
        default_lr_hadronic_summary["default_artifacts_unchanged"] = (
            default_artifacts_unchanged
        )
        custom_checks["default_q1_hadronic_bundle_stays_q1_only"] = (
            default_q1_hadronic_bundle_stays_q1_only
        )
        custom_checks["default_observables_unchanged"] = default_observables_unchanged
        custom_checks["default_artifacts_unchanged"] = default_artifacts_unchanged
        _append_post_aggregation_failure(
            "default_lr_hadronic_probe",
            "default_q1_hadronic_bundle_stays_q1_only",
            default_q1_hadronic_bundle_stays_q1_only,
        )
        _append_post_aggregation_failure(
            "default_lr_hadronic_probe",
            "default_observables_unchanged",
            default_observables_unchanged,
        )
        _append_post_aggregation_failure(
            "default_lr_hadronic_probe",
            "default_artifacts_unchanged",
            default_artifacts_unchanged,
        )

    if lr_r_chi_freeze_summary.get("status") == "ok":
        custom_checks = lr_r_chi_freeze_summary.setdefault("checks", {})
        default_hadronic_bundle_stays_q1_only = (
            hadronic_summary.get("status") == "ok"
            and hadronic_summary.get("checks", {}).get("supported_operator_subset_is_pr5a")
            is True
            and hadronic_summary.get("checks", {}).get("unsupported_operator_subset_is_lr_only")
            is True
        )
        default_observables_unchanged = (
            observable_summary.get("status") == "ok"
            and set(observable_summary.get("observable_names", [])) == SUPPORTED_OBSERVABLE_ROWS
        )
        default_artifacts_unchanged = (
            artifact_summary.get("status") == "ok"
            and set(artifact_summary.get("observable_rows", []))
            == EXPECTED_ARTIFACT_OBSERVABLE_ROWS
        )
        lr_r_chi_freeze_summary["default_hadronic_bundle_stays_q1_only"] = (
            default_hadronic_bundle_stays_q1_only
        )
        lr_r_chi_freeze_summary["default_observables_unchanged"] = (
            default_observables_unchanged
        )
        lr_r_chi_freeze_summary["default_artifacts_unchanged"] = (
            default_artifacts_unchanged
        )
        custom_checks["default_hadronic_bundle_stays_q1_only"] = (
            default_hadronic_bundle_stays_q1_only
        )
        custom_checks["default_observables_unchanged"] = default_observables_unchanged
        custom_checks["default_artifacts_unchanged"] = default_artifacts_unchanged
        _append_post_aggregation_failure(
            "lr_r_chi_freeze_probe",
            "default_hadronic_bundle_stays_q1_only",
            default_hadronic_bundle_stays_q1_only,
        )
        _append_post_aggregation_failure(
            "lr_r_chi_freeze_probe",
            "default_observables_unchanged",
            default_observables_unchanged,
        )
        _append_post_aggregation_failure(
            "lr_r_chi_freeze_probe",
            "default_artifacts_unchanged",
            default_artifacts_unchanged,
        )
        default_observable_values_match_frozen_reference = (
            observable_summary.get("checks", {}).get("m12_re_matches_frozen_reference") is True
            and observable_summary.get("checks", {}).get("m12_im_matches_frozen_reference")
            is True
            and observable_summary.get("checks", {}).get("delta_m_matches_frozen_reference")
            is True
        )
        custom_checks["default_observable_values_match_frozen_reference"] = (
            default_observable_values_match_frozen_reference
        )
        _append_post_aggregation_failure(
            "lr_r_chi_freeze_probe",
            "default_observable_values_match_frozen_reference",
            default_observable_values_match_frozen_reference,
        )
        default_artifact_exports_match_tracked_files = (
            artifact_summary.get("checks", {}).get(
                "tracked_default_exports_match_current_export"
            )
            is True
        )
        custom_checks["default_artifact_exports_match_tracked_files"] = (
            default_artifact_exports_match_tracked_files
        )
        _append_post_aggregation_failure(
            "lr_r_chi_freeze_probe",
            "default_artifact_exports_match_tracked_files",
            default_artifact_exports_match_tracked_files,
        )

    if custom_combined_observable_summary.get("status") == "ok":
        custom_checks = custom_combined_observable_summary.setdefault("checks", {})
        default_observables_unchanged = (
            observable_summary.get("status") == "ok"
            and set(observable_summary.get("observable_names", [])) == SUPPORTED_OBSERVABLE_ROWS
        )
        default_artifacts_unchanged = (
            artifact_summary.get("status") == "ok"
            and set(artifact_summary.get("observable_rows", []))
            == EXPECTED_ARTIFACT_OBSERVABLE_ROWS
        )
        custom_combined_observable_summary["default_observables_unchanged"] = (
            default_observables_unchanged
        )
        custom_combined_observable_summary["default_artifacts_unchanged"] = (
            default_artifacts_unchanged
        )
        custom_checks["default_observables_unchanged"] = default_observables_unchanged
        custom_checks["default_artifacts_unchanged"] = default_artifacts_unchanged
        _append_post_aggregation_failure(
            "custom_combined_observable_probe",
            "default_observables_unchanged",
            default_observables_unchanged,
        )
        _append_post_aggregation_failure(
            "custom_combined_observable_probe",
            "default_artifacts_unchanged",
            default_artifacts_unchanged,
        )

    for custom_b_summary in (
        custom_bd_observable_summary,
        custom_bs_observable_summary,
        custom_d0_observable_summary,
    ):
        if custom_b_summary.get("status") != "ok":
            continue
        custom_checks = custom_b_summary.setdefault("checks", {})
        default_observables_unchanged = (
            observable_summary.get("status") == "ok"
            and set(observable_summary.get("observable_names", [])) == SUPPORTED_OBSERVABLE_ROWS
        )
        default_artifacts_unchanged = (
            artifact_summary.get("status") == "ok"
            and set(artifact_summary.get("observable_rows", []))
            == EXPECTED_ARTIFACT_OBSERVABLE_ROWS
        )
        custom_b_summary["default_observables_unchanged"] = default_observables_unchanged
        custom_b_summary["default_artifacts_unchanged"] = default_artifacts_unchanged
        custom_checks["default_observables_unchanged"] = default_observables_unchanged
        custom_checks["default_artifacts_unchanged"] = default_artifacts_unchanged
        _append_post_aggregation_failure(
            f"custom_{custom_b_summary.get('system_id', 'unknown')}_observable_probe",
            "default_observables_unchanged",
            default_observables_unchanged,
        )
        _append_post_aggregation_failure(
            f"custom_{custom_b_summary.get('system_id', 'unknown')}_observable_probe",
            "default_artifacts_unchanged",
            default_artifacts_unchanged,
        )

    summary = {
        "static": static_summary,
        "imports": import_summary,
        "conventions": convention_summary,
        "scales": scale_summary,
        "structural_benchmark": structural_summary,
        "pr1_benchmark": pr1_summary,
        "couplings": coupling_summary,
        "eft_matching": matching_summary,
        "eft_rg_lo": rg_summary,
        "eft_rg_lr_probe": lr_running_probe_summary,
        "lr_contract_freeze": lr_freeze_summary,
        "hadronic_inputs": hadronic_summary,
        "default_lr_hadronic_probe": default_lr_hadronic_summary,
        "lr_r_chi_freeze_probe": lr_r_chi_freeze_summary,
        "custom_lr_hadronic_probe": custom_lr_hadronic_summary,
        "custom_lr_observable_probe": custom_lr_observable_summary,
        "custom_combined_observable_probe": custom_combined_observable_summary,
        "custom_bd_observable_probe": custom_bd_observable_summary,
        "custom_bs_observable_probe": custom_bs_observable_summary,
        "custom_d0_observable_probe": custom_d0_observable_summary,
        "observables": observable_summary,
        "artifacts": artifact_summary,
    }
    return summary, failures


def _cross_process_summary(require_package: bool) -> dict[str, Any]:
    command = [sys.executable, str(Path(__file__).resolve())]
    if require_package:
        command.append("--require-package")
    command.append("--emit-json")
    completed = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def _cross_process_default_matching_payload() -> dict[str, Any]:
    command = [sys.executable, str(Path(__file__).resolve()), "--emit-default-matching-json"]
    completed = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def _cross_process_default_rg_payload() -> dict[str, Any]:
    command = [sys.executable, str(Path(__file__).resolve()), "--emit-default-rg-json"]
    completed = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def _cross_process_default_hadronic_payload() -> dict[str, Any]:
    command = [sys.executable, str(Path(__file__).resolve()), "--emit-default-hadronic-json"]
    completed = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def _cross_process_default_lr_hadronic_payload() -> dict[str, Any]:
    command = [sys.executable, str(Path(__file__).resolve()), "--emit-default-lr-hadronic-json"]
    completed = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def _cross_process_lr_r_chi_freeze_payload() -> dict[str, Any]:
    command = [sys.executable, str(Path(__file__).resolve()), "--emit-lr-r-chi-freeze-json"]
    completed = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def _cross_process_default_observable_payload() -> dict[str, Any]:
    command = [sys.executable, str(Path(__file__).resolve()), "--emit-default-observable-json"]
    completed = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def _cross_process_custom_b_observable_probe_payload(system_id: str) -> dict[str, Any]:
    emit_flag = {
        "B_d": "--emit-custom-bd-observable-probe-json",
        "B_s": "--emit-custom-bs-observable-probe-json",
        "D0": "--emit-custom-d0-observable-probe-json",
    }[system_id]
    command = [sys.executable, str(Path(__file__).resolve()), emit_flag]
    completed = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    return json.loads(completed.stdout)


def main() -> int:
    args = _parse_args()

    if not PAPER_PACKAGE_DIR.exists():
        if args.require_package:
            return 1
        if (
            args.emit_json
            or args.emit_default_matching_json
            or args.emit_default_rg_json
            or args.emit_default_hadronic_json
            or args.emit_default_lr_hadronic_json
            or args.emit_default_observable_json
            or args.emit_lr_r_chi_freeze_json
            or args.emit_custom_bd_observable_probe_json
            or args.emit_custom_bs_observable_probe_json
            or args.emit_custom_d0_observable_probe_json
        ):
            print(json.dumps({"status": "package_missing"}, sort_keys=True))
        else:
            print("paper_0710_1869 acceptance benchmark")
            print(f"  repo_root    = {REPO_ROOT}")
            print(f"  package_dir  = {PAPER_PACKAGE_DIR.relative_to(REPO_ROOT)}")
            print("  status       = package missing in this checkout")
            print("SKIP: acceptance benchmark deferred until the paper package exists.")
        return 0

    if args.emit_default_matching_json:
        imported_modules, _, _ = _import_paper_modules()
        module = imported_modules.get("eft_deltaf2.matching_kkgluon")
        if module is None:
            print(json.dumps({"status": "missing"}, sort_keys=True))
            return 0
        print(json.dumps(_default_matching_payload(module), sort_keys=True))
        return 0

    if args.emit_default_rg_json:
        imported_modules, _, _ = _import_paper_modules()
        module = imported_modules.get("eft_deltaf2.rg")
        if module is None:
            print(json.dumps({"status": "missing"}, sort_keys=True))
            return 0
        print(json.dumps(_default_rg_payload(imported_modules), sort_keys=True))
        return 0

    if args.emit_default_hadronic_json:
        imported_modules, _, _ = _import_paper_modules()
        module = imported_modules.get("eft_deltaf2.hadronic")
        if module is None:
            print(json.dumps({"status": "missing"}, sort_keys=True))
            return 0
        print(json.dumps(_default_hadronic_payload(imported_modules), sort_keys=True))
        return 0

    if args.emit_default_lr_hadronic_json:
        imported_modules, _, _ = _import_paper_modules()
        module = imported_modules.get("eft_deltaf2.hadronic")
        if module is None:
            print(json.dumps({"status": "missing"}, sort_keys=True))
            return 0
        try:
            print(json.dumps(_default_kaon_lr_hadronic_payload(imported_modules), sort_keys=True))
        except AssertionError:
            print(json.dumps({"status": "missing"}, sort_keys=True))
        return 0

    if args.emit_lr_r_chi_freeze_json:
        imported_modules, _, _ = _import_paper_modules()
        module = imported_modules.get("eft_deltaf2.hadronic")
        if module is None:
            print(json.dumps({"status": "missing"}, sort_keys=True))
            return 0
        try:
            print(
                json.dumps(
                    _default_kaon_lr_r_chi_freeze_payload(imported_modules),
                    sort_keys=True,
                )
            )
        except AssertionError:
            print(json.dumps({"status": "missing"}, sort_keys=True))
        return 0

    if args.emit_default_observable_json:
        imported_modules, _, _ = _import_paper_modules()
        module = imported_modules.get("eft_deltaf2.observables")
        if module is None:
            print(json.dumps({"status": "missing"}, sort_keys=True))
            return 0
        print(json.dumps(_default_observable_payload(imported_modules), sort_keys=True))
        return 0

    if (
        args.emit_custom_bd_observable_probe_json
        or args.emit_custom_bs_observable_probe_json
        or args.emit_custom_d0_observable_probe_json
    ):
        imported_modules, _, _ = _import_paper_modules()
        if args.emit_custom_bd_observable_probe_json:
            system_id = "B_d"
        elif args.emit_custom_bs_observable_probe_json:
            system_id = "B_s"
        else:
            system_id = "D0"
        try:
            payload, _, _ = _build_custom_b_observable_probe_payload(imported_modules, system_id)
        except AssertionError:
            print(json.dumps({"status": "missing"}, sort_keys=True))
            return 0
        print(json.dumps(payload, sort_keys=True))
        return 0

    if args.verify_artifacts_dir is not None:
        payload = _run_independent_verifier(artifact_dir=args.verify_artifacts_dir)
        print(json.dumps(payload, sort_keys=True))
        return 0 if payload.get("ok") else 1

    summary, failures = _collect_summary(export_dir=args.export_artifacts_dir)
    if not failures and not args.emit_json:
        cross_process_summary = _cross_process_summary(args.require_package)
        if json.dumps(summary, sort_keys=True) != json.dumps(cross_process_summary, sort_keys=True):
            failures.append("cross-process summary mismatch indicates nondeterministic output")

    if args.emit_json:
        print(json.dumps(summary, sort_keys=True))
        return 0 if not failures else 1

    print("paper_0710_1869 acceptance benchmark")
    print(f"  repo_root    = {REPO_ROOT}")
    print(f"  package_dir  = {PAPER_PACKAGE_DIR.relative_to(REPO_ROOT)}")
    print("Contract summary:")
    print(json.dumps(summary, indent=2, sort_keys=True))

    if failures:
        print("FAILURES:")
        for item in failures:
            print(f"  - {item}")
        return 1

    print("PASS: paper_0710_1869 contracts are consistent with this benchmark probe.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
