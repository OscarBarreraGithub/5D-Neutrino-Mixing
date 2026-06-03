"""Production tests for K009 (K_L -> pi0 mu+ mu- direct CP)."""

from __future__ import annotations

import math
from pathlib import Path
import re

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.kaon import K009 as k009_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.model import RotationParameters, ckm_like_unitary
from quarkConstraints.rare_kaon_dilepton import (
    KLONG_PI0EE_SM_Y7A_BAR,
    KLONG_PI0EE_SM_Y7V_BAR,
)
from tests.rare_kaon_phase3d_helpers import (
    core_y7_wilsons_from_rs_coeff,
    rs_coeff,
    sample_rare_kaon_point,
    scaled_rare_kaon_point,
    sm_limit_rare_kaon_point,
)

_PID = "K009"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "kaon" / "K009.yaml"
_LIMIT_RE = re.compile(r"^\s*<\s*(?P<value>[0-9.eE+-]+)\s*$")
_LEADING_FLOAT_RE = re.compile(
    r"^\s*(?P<value>[+-]?[0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)(?:\s+.*)?$"
)
_VALUE_WITH_UNCERT_RE = re.compile(
    r"^\s*\(?\s*(?P<value>[0-9.eE+-]+)\s*\+/-\s*"
    r"(?P<uncertainty>[0-9.eE+-]+)\s*\)?(?:\s*x\s*10\^(?P<exponent>[+-]?[0-9]+))?\s*$"
)


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _value_entry(observable: str):
    for entry in _yaml()["pdg_or_equivalent"]["values"]:
        if entry["observable"] == observable:
            return entry
    raise AssertionError(f"no K009 value entry for {observable!r}")


def _limit_value(value: str) -> float:
    match = _LIMIT_RE.match(value)
    assert match is not None
    return float(match.group("value"))


def _leading_float(value: str) -> float:
    match = _LEADING_FLOAT_RE.match(value)
    assert match is not None
    return float(match.group("value"))


def _value_with_uncertainty(value: str) -> tuple[float, float]:
    match = _VALUE_WITH_UNCERT_RE.match(value)
    assert match is not None
    exponent = int(match.group("exponent") or "0")
    scale = 10.0**exponent
    return (
        float(match.group("value")) * scale,
        float(match.group("uncertainty")) * scale,
    )


def _sd_couplings(
    left: complex,
    right: complex = 0.0j,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with only the s-d slot populated."""

    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[0, 1] = left
    left_down[1, 0] = np.conj(left)
    right_down[0, 1] = right
    right_down[1, 0] = np.conj(right)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=zeros,
        left_down=left_down,
        right_up=zeros,
        right_down=right_down,
    )


def _independent_pi0mumu_components(
    inputs,
    chpt,
    lambda_y7v_np_proxy: complex = 0.0j,
    lambda_y7a_np_proxy: complex = 0.0j,
):
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=inputs.theta12,
            theta13=inputs.theta13,
            theta23=inputs.theta23,
            delta=inputs.delta,
        )
    )
    lambda_t = complex(np.conjugate(matrix[2, 1]) * matrix[2, 0])
    lambda_y7v_eff = lambda_t * KLONG_PI0EE_SM_Y7V_BAR + lambda_y7v_np_proxy
    lambda_y7a_eff = lambda_t * KLONG_PI0EE_SM_Y7A_BAR + lambda_y7a_np_proxy
    vector_ratio = lambda_y7v_eff.imag / 1.0e-4
    axial_ratio = lambda_y7a_eff.imag / 1.0e-4
    y7_norm = KLONG_PI0EE_SM_Y7V_BAR**2 + KLONG_PI0EE_SM_Y7A_BAR**2
    c_dir_y7va = chpt.c_dir / y7_norm
    c_int_y7v = chpt.c_int / abs(KLONG_PI0EE_SM_Y7V_BAR)
    ratio = math.sqrt(vector_ratio * vector_ratio + axial_ratio * axial_ratio)
    direct = chpt.normalization * c_dir_y7va * (
        vector_ratio * vector_ratio + axial_ratio * axial_ratio
    )
    indirect = chpt.normalization * chpt.c_mix * chpt.a_s_abs * chpt.a_s_abs
    interference_abs = (
        chpt.normalization * c_int_y7v * chpt.a_s_abs * abs(vector_ratio)
    )
    cpc = chpt.normalization * chpt.c_cpc
    return {
        "lambda_t": lambda_t,
        "lambda_y7v_eff": lambda_y7v_eff,
        "lambda_y7a_eff": lambda_y7a_eff,
        "ratio": ratio,
        "vector_ratio": vector_ratio,
        "axial_ratio": axial_ratio,
        "c_dir_y7va": c_dir_y7va,
        "c_int_y7v": c_int_y7v,
        "direct": direct,
        "indirect": indirect,
        "interference_abs": interference_abs,
        "cpc": cpc,
        "constructive": indirect + direct + cpc + interference_abs,
        "destructive": max(indirect + direct + cpc - interference_abs, 0.0),
    }


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "kaon"
    assert constraint.observable == "BR(K_L -> pi0 mu+ mu-)_direct CP"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml()["pdg_or_equivalent"]
    ktev_limit = _value_entry(k009_module._KTEV_LIMIT_OBS)
    ktev_observed = _value_entry(k009_module._KTEV_OBSERVED_OBS)
    ktev_background = _value_entry(k009_module._KTEV_BACKGROUND_OBS)
    c_mix = _value_entry(k009_module._C_MIX_OBS)
    c_int = _value_entry(k009_module._C_INT_OBS)
    c_dir = _value_entry(k009_module._C_DIR_OBS)
    c_cpc = _value_entry(k009_module._C_CPC_OBS)
    a_s = _value_entry(k009_module._A_S_OBS)
    sm_constructive = _value_entry(k009_module._SM_CONSTRUCTIVE_OBS)
    theory_coeffs = _yaml()["theory_decomposition"]["muon_coefficients"]
    sm_destructive, sm_destructive_unc = _value_with_uncertainty(
        theory_coeffs["SM_destructive_interference"]
    )
    sm_im, sm_im_unc = _value_with_uncertainty(
        theory_coeffs["Im_lambda_t_over_1e_minus4_SM"]
    )

    assert constraint.anchor.experimental_limit.value == pytest.approx(
        _limit_value(pdg["value"])
    )
    assert constraint.anchor.experimental_limit.source_url == pdg["source_url"]
    assert constraint.anchor.experimental_limit.is_upper_limit is True
    assert constraint.anchor.ktev_limit.value == pytest.approx(
        _limit_value(ktev_limit["value"])
    )
    assert constraint.anchor.ktev_observed_candidates.value == pytest.approx(
        float(ktev_observed["value"])
    )
    assert constraint.anchor.ktev_expected_background.value == pytest.approx(
        float(ktev_background["value"])
    )
    assert constraint.anchor.ktev_expected_background.uncertainty == pytest.approx(
        float(ktev_background["uncertainty"])
    )
    assert constraint.anchor.c_mix.value == pytest.approx(
        _leading_float(c_mix["value"])
    )
    assert constraint.anchor.c_mix.uncertainty == pytest.approx(
        _leading_float(c_mix["uncertainty"])
    )
    assert constraint.anchor.c_int.value == pytest.approx(
        _leading_float(c_int["value"])
    )
    assert constraint.anchor.c_dir.value == pytest.approx(float(c_dir["value"]))
    assert constraint.anchor.c_cpc.value == pytest.approx(float(c_cpc["value"]))
    assert constraint.anchor.a_s_abs.value == pytest.approx(float(a_s["value"]))
    assert constraint.anchor.sm_im_lambda_t_over_1e_minus4.value == pytest.approx(
        sm_im
    )
    assert constraint.anchor.sm_im_lambda_t_over_1e_minus4.uncertainty == (
        pytest.approx(sm_im_unc)
    )
    assert constraint.anchor.sm_constructive_total.value == pytest.approx(
        float(sm_constructive["value"])
    )
    assert constraint.anchor.sm_constructive_total.uncertainty == pytest.approx(
        float(sm_constructive["uncertainty"])
    )
    assert constraint.anchor.sm_destructive_total.value == pytest.approx(
        sm_destructive
    )
    assert constraint.anchor.sm_destructive_total.uncertainty == pytest.approx(
        sm_destructive_unc
    )
    assert constraint.anchor.budget == pytest.approx(_limit_value(pdg["value"]))
    assert constraint.anchor.budget_band.sm_subtracted is False
    assert abs(
        constraint.sm_result.constructive_total_branching_fraction
        - constraint.anchor.sm_constructive_total.value
    ) < constraint.anchor.sm_constructive_total.uncertainty
    assert abs(
        constraint.sm_result.destructive_total_branching_fraction
        - constraint.anchor.sm_destructive_total.value
    ) < constraint.anchor.sm_destructive_total.uncertainty


def test_k009_anchor_loud_fail_probe():
    with pytest.raises(k009_module.AnchorError):
        k009_module._load_value_anchor(_PID, "not a K009 observable")


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())
    old_style_result = constraint.evaluate(
        point_builder.build_from_quark_couplings(_sd_couplings(left=1.0e-2j))
    )

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.sm_prediction == pytest.approx(
        fcc.get(_PID).sm_result.direct_cp_branching_fraction
    )
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert old_style_result.passes is True
    assert old_style_result.predicted is None
    assert old_style_result.diagnostics["evaluated"] is False
    assert old_style_result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert old_style_result.diagnostics["legacy_quark_mass_basis_couplings_present"] is True


def test_sm_limit_branching_fraction_matches_independent_recomputation():
    constraint = fcc.get(_PID)
    point = sm_limit_rare_kaon_point()
    result = constraint.evaluate(point)
    independent = _independent_pi0mumu_components(
        constraint.sm_inputs,
        constraint.chpt_inputs,
    )

    assert result.predicted == pytest.approx(independent["direct"])
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.diagnostics["constructive_total_branching_fraction"] == (
        pytest.approx(independent["constructive"])
    )
    assert result.diagnostics["destructive_total_branching_fraction"] == (
        pytest.approx(independent["destructive"])
    )
    assert result.diagnostics["direct_cp_amplitude_ratio"] == pytest.approx(
        independent["ratio"]
    )
    assert result.diagnostics["direct_cp_vector_amplitude_ratio"] == pytest.approx(
        independent["vector_ratio"]
    )
    assert result.diagnostics["direct_cp_axial_amplitude_ratio"] == pytest.approx(
        independent["axial_ratio"]
    )
    assert result.diagnostics["cpc_branching_fraction"] == pytest.approx(
        5.2e-12
    )
    assert result.predicted < constraint.anchor.budget
    assert result.diagnostics["sm_constructive_anchor_pull"] == pytest.approx(
        (
            independent["constructive"]
            - constraint.anchor.sm_constructive_total.value
        )
        / constraint.anchor.sm_constructive_total.uncertainty
    )
    assert abs(result.diagnostics["sm_constructive_anchor_pull"]) < 1.0
    assert abs(result.diagnostics["sm_destructive_anchor_pull"]) < 1.0
    assert result.diagnostics["lambda_y7v_np_proxy"] == pytest.approx(0.0j)
    assert result.diagnostics["lambda_y7a_np_proxy"] == pytest.approx(0.0j)
    assert result.diagnostics["wilson_prefactor_reused"] is False
    assert result.diagnostics["second_mkk_suppression_applied"] is False
    assert result.passes is True


def test_y7_structure_separates_direct_cp_and_vector_interference():
    constraint = fcc.get(_PID)
    sm = _independent_pi0mumu_components(constraint.sm_inputs, constraint.chpt_inputs)
    axial_only = _independent_pi0mumu_components(
        constraint.sm_inputs,
        constraint.chpt_inputs,
        lambda_y7a_np_proxy=2.0e-4j,
    )
    vector_only = _independent_pi0mumu_components(
        constraint.sm_inputs,
        constraint.chpt_inputs,
        lambda_y7v_np_proxy=2.0e-4j,
    )

    assert axial_only["direct"] > sm["direct"]
    assert axial_only["interference_abs"] == pytest.approx(
        sm["interference_abs"]
    )
    assert vector_only["interference_abs"] != pytest.approx(
        sm["interference_abs"]
    )


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    point = sample_rare_kaon_point()
    result = fcc.get(_PID).evaluate(point)

    assert result.process_id == _PID
    for value in (
        result.predicted,
        result.ratio,
        result.budget,
        result.sm_prediction,
        result.experimental,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    for key in (
        "left_sd_coupling",
        "right_sd_coupling",
        "lambda_t",
        "lambda_y7v_eff",
        "lambda_y7a_eff",
        "lambda_y7v_np_proxy",
        "lambda_y7a_np_proxy",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "direct_cp_branching_fraction",
        "constructive_total_branching_fraction",
        "destructive_total_branching_fraction",
        "direct_cp_amplitude_ratio",
        "direct_cp_vector_amplitude_ratio",
        "direct_cp_axial_amplitude_ratio",
        "np_shift_branching_fraction",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["short_distance_direct_cp_only"] is True
    assert result.diagnostics["interference_sign_fixed"] is False
    assert result.diagnostics["uses_y7v_y7a_wilson_structure"] is True
    assert result.diagnostics["uses_k008_y7v_y7a_rs_proxy"] is False
    assert result.diagnostics["uses_k008_y7v_y7a_core_inputs"] is True
    assert result.diagnostics["interference_uses_vector_y7v_only"] is True
    assert result.diagnostics["semileptonic_qcd_running_applied"] is False
    assert result.diagnostics["semileptonic_qcd_running_effect_fraction"] == 0.0
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert result.diagnostics["rs_semileptonic_matching_status"].endswith(
        "no_second_1_over_M_KK_squared"
    )
    assert result.diagnostics["rare_kaon_proxy_reused"] is False


@pytest.mark.parametrize(
    ("scale", "expected_pass"),
    [
        (1.0, True),
        (1.0e5, False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    scale: float,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = scaled_rare_kaon_point(scale=scale)
    result = constraint.evaluate(point)
    wilsons = core_y7_wilsons_from_rs_coeff(
        rs_coeff(point, lepton="mu"),
        inputs=constraint.sm_inputs,
        matching_scale_gev=point.extras["kk_ew_mass_gev"],
    )
    direct = _independent_pi0mumu_components(
        constraint.sm_inputs,
        constraint.chpt_inputs,
        lambda_y7v_np_proxy=wilsons.lambda_y7v_np_proxy,
        lambda_y7a_np_proxy=wilsons.lambda_y7a_np_proxy,
    )

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(direct["direct"])
    assert result.sm_prediction == pytest.approx(
        constraint.sm_result.direct_cp_branching_fraction
    )
    assert result.ratio == pytest.approx(direct["direct"] / constraint.anchor.budget)
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_is_diagnostic_only_for_rigorous_wilsons():
    default_point = scaled_rare_kaon_point(scale=20.0)
    ew_point = scaled_rare_kaon_point(scale=20.0, kk_ew_mass_gev=6000.0)
    default_result = fcc.get(_PID).evaluate(default_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert ew_result.predicted == pytest.approx(default_result.predicted)
    assert ew_result.diagnostics["lambda_y7a_np_proxy"] == pytest.approx(
        default_result.diagnostics["lambda_y7a_np_proxy"]
    )
    assert ew_result.diagnostics["second_mkk_suppression_applied"] is False


def test_evaluate_is_pure_and_deterministic():
    point = sample_rare_kaon_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
