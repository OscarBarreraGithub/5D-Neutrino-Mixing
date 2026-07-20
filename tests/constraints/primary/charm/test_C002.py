"""Production tests for C002 (neutral-D indirect CP violation)."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.deltaf2 import (
    B_1_D,
    B_4_D,
    B_5_D,
    DEFAULT_DELTA_F2_INPUTS_V1,
    F_D,
    M_C_QUARK,
    M_D0,
    M_U_QUARK,
    _evolve_wilsons,
    compute_delta_f2_wilsons,
    compute_m12_np,
    evaluate_d0_mixing_with_running,
)

_PID = "C002"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charm" / "C002.yaml"
_C001_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charm" / "C001.yaml"
_MU_HAD_GEV = 3.0
_FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
_PLUS_MINUS_RE = re.compile(
    rf"(?P<value>{_FLOAT_RE})\s+\+(?P<plus>{_FLOAT_RE})\/-(?P<minus>{_FLOAT_RE})"
)
_CI_RE = re.compile(rf"\[(?P<low>{_FLOAT_RE}),\s*(?P<high>{_FLOAT_RE})\]")


def _yaml_pdg_block(path: Path = _SIDECAR):
    with open(path) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _snapshot_phi_m_values(snapshot_path: str):
    path = _REPO_ROOT / snapshot_path
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped.startswith("phiM") or "(degrees)" not in stripped:
            continue
        all_cpv = list(_PLUS_MINUS_RE.finditer(stripped))[-1]
        ci = _CI_RE.search(stripped)
        assert ci is not None
        return {
            "value": float(all_cpv.group("value")),
            "uncertainty_plus": float(all_cpv.group("plus")),
            "uncertainty_minus": float(all_cpv.group("minus")),
            "confidence_interval_95_percent": (
                float(ci.group("low")),
                float(ci.group("high")),
            ),
        }
    raise AssertionError("phiM row missing from C002 HFLAV snapshot")


def _uc_couplings(
    left: complex,
    right: complex,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with only the u-c slot populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_up = zeros.copy()
    right_up = zeros.copy()
    left_up[0, 1] = left
    left_up[1, 0] = np.conj(left)
    right_up[0, 1] = right
    right_up[1, 0] = np.conj(right)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=left_up,
        left_down=zeros.copy(),
        right_up=right_up,
        right_down=zeros.copy(),
    )


def _d0_wilsons(couplings: QuarkMassBasisCouplings):
    wilsons = compute_delta_f2_wilsons(
        couplings,
        inputs=DEFAULT_DELTA_F2_INPUTS_V1,
    )
    return next(w for w in wilsons if w.input.key == "d")


def _expected_cpv_numbers(couplings: QuarkMassBasisCouplings):
    constraint = fcc.get(_PID)
    wilsons = _d0_wilsons(couplings)
    evolved = _evolve_wilsons(wilsons, mu_had=_MU_HAD_GEV)
    core_mixing = evaluate_d0_mixing_with_running(wilsons, mu_had=_MU_HAD_GEV)
    m12_np = complex(
        compute_m12_np(
            evolved,
            F_D,
            M_D0,
            M_C_QUARK,
            M_U_QUARK,
            B_1_D,
            B_4_D,
            B_5_D,
        )
    )
    predicted = abs(float(m12_np.imag)) / constraint.anchor.m12_budget_gev
    ratio = predicted / constraint.anchor.budget
    return core_mixing, m12_np, predicted, ratio


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charm"
    assert constraint.observable == "D0-D0bar indirect CP violation"


def test_anchor_matches_yaml():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    c001_pdg = _yaml_pdg_block(_C001_SIDECAR)
    q_over_p = pdg["q_over_p"]
    phi_d = pdg["phi_D"]
    phi_m = _snapshot_phi_m_values(phi_d["snapshot_path"])
    no_cpv = pdg["no_indirect_cpv_test"]
    delta_m_d = c001_pdg["Delta_m_D"]

    assert constraint.anchor.q_over_p.value == pytest.approx(q_over_p["value"])
    assert constraint.anchor.q_over_p.uncertainty_plus == pytest.approx(
        q_over_p["uncertainty_plus"]
    )
    assert constraint.anchor.q_over_p.uncertainty_minus == pytest.approx(
        q_over_p["uncertainty_minus"]
    )
    assert constraint.anchor.q_over_p.confidence_interval_95 == pytest.approx(
        tuple(q_over_p["confidence_interval_95_percent"])
    )
    assert constraint.anchor.q_over_p.source_url == q_over_p["source_url"]
    assert constraint.anchor.phi_d.value == pytest.approx(phi_d["value"])
    assert constraint.anchor.phi_d.uncertainty_plus == pytest.approx(
        phi_d["uncertainty_plus"]
    )
    assert constraint.anchor.phi_d.uncertainty_minus == pytest.approx(
        phi_d["uncertainty_minus"]
    )
    assert constraint.anchor.phi_d.confidence_interval_95 == pytest.approx(
        tuple(phi_d["confidence_interval_95_percent"])
    )
    assert constraint.anchor.phi_d.source_url == phi_d["source_url"]
    assert constraint.anchor.phi_m.value == pytest.approx(phi_m["value"])
    assert constraint.anchor.phi_m.uncertainty_plus == pytest.approx(
        phi_m["uncertainty_plus"]
    )
    assert constraint.anchor.phi_m.uncertainty_minus == pytest.approx(
        phi_m["uncertainty_minus"]
    )
    assert constraint.anchor.phi_m.confidence_interval_95 == pytest.approx(
        phi_m["confidence_interval_95_percent"]
    )
    assert constraint.anchor.phi_m.source_url == phi_d["source_url"]
    assert constraint.anchor.no_indirect_cpv.delta_chi2 == pytest.approx(
        no_cpv["delta_chi2"]
    )
    assert constraint.anchor.no_indirect_cpv.significance_sigma == pytest.approx(
        no_cpv["significance_sigma"]
    )
    assert constraint.anchor.delta_m_d_gev.value == pytest.approx(
        delta_m_d["value_GeV"]
    )
    assert constraint.anchor.m12_budget_gev == pytest.approx(
        delta_m_d["value_GeV"] / 2.0
    )

    q_budget = max(
        abs(q_over_p["confidence_interval_95_percent"][0] - 1.0),
        abs(q_over_p["confidence_interval_95_percent"][1] - 1.0),
    )
    phi_budget_degrees = max(
        abs(phi_d["confidence_interval_95_percent"][0]),
        abs(phi_d["confidence_interval_95_percent"][1]),
    )
    phi_budget_sine = abs(math.sin(math.radians(phi_budget_degrees)))
    phi_m_budget_degrees = max(
        abs(phi_m["confidence_interval_95_percent"][0]),
        abs(phi_m["confidence_interval_95_percent"][1]),
    )
    phi_m_budget_sine = abs(math.sin(math.radians(phi_m_budget_degrees)))
    central_proxy = abs(math.sin(math.radians(phi_m["value"])))

    assert constraint.anchor.budget_band.q_over_p_95_deviation_budget == pytest.approx(
        q_budget
    )
    assert constraint.anchor.budget_band.phi_d_95_deviation_degrees == pytest.approx(
        phi_budget_degrees
    )
    assert constraint.anchor.budget_band.phi_d_95_deviation_sine_budget == pytest.approx(
        phi_budget_sine
    )
    assert constraint.anchor.budget_band.phi_m_95_deviation_degrees == pytest.approx(
        phi_m_budget_degrees
    )
    assert constraint.anchor.budget_band.phi_m_95_deviation_sine_budget == pytest.approx(
        phi_m_budget_sine
    )
    assert constraint.anchor.budget == pytest.approx(phi_m_budget_sine)
    assert constraint.anchor.budget == pytest.approx(0.0258280005)
    assert constraint.anchor.value == pytest.approx(central_proxy)
    assert constraint.anchor.value == pytest.approx(0.0005235988)
    assert constraint.anchor.m12_budget_gev == pytest.approx(3.281e-15)
    assert constraint.anchor.cp_odd_m12_budget_gev == pytest.approx(8.47417e-17)


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"


def test_evaluate_runs_end_to_end_with_real_couplings_and_real_finite_fields():
    couplings = _uc_couplings(
        left=1.0e-4 + 0.5e-4j,
        right=0.7e-4 + 0.2e-4j,
    )
    point = point_builder.build_from_quark_couplings(couplings)
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
    assert isinstance(result.diagnostics["m12_np_gev"], complex)
    assert isinstance(result.diagnostics["left_uc_coupling"], complex)
    assert isinstance(result.diagnostics["wilson_coefficients"]["C4_LR"], complex)
    for key in (
        "re_m12_np_gev",
        "im_m12_np_gev",
        "abs_m12_np_gev",
        "m12_np_phase_degrees",
        "m12_np_raw_argument_degrees",
        "m12_np_phase_sine_abs",
        "amplitude_ratio_to_delta_m_d_half",
        "cp_odd_fraction",
        "q_over_p_deviation_proxy",
        "phi_d_sine_deviation_proxy",
        "hadronic_scale_gev",
        "matching_scale_gev",
        "m_kk_gev",
        "m12_normalization_budget_gev",
        "cp_odd_m12_budget_gev",
        "delta_m_d_experimental_gev",
        "delta_m_d_uncertainty_gev",
        "q_over_p_value",
        "q_over_p_95_deviation_budget",
        "phi_d_95_sine_budget",
        "phi_m_value_degrees",
        "phi_m_95_sine_budget",
        "no_indirect_cpv_delta_chi2",
        "no_indirect_cpv_significance_sigma",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert "q_over_p_phase_proxy_degrees" not in result.diagnostics
    assert "not a q/p or phi_D prediction" in result.diagnostics[
        "m12_np_raw_argument_note"
    ]
    assert result.diagnostics["qcd_running_applied"] is True
    assert result.diagnostics["hadronic_scale_gev"] == pytest.approx(_MU_HAD_GEV)
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["sm_long_distance_phase_grounded"] is False
    assert result.diagnostics["full_hflav_likelihood_used"] is False


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_uc_couplings(left=1.0e-4, right=0.7e-4), True),
        (_uc_couplings(left=2.0e-3 + 1.0e-3j, right=1.4e-3 + 0.4e-3j), False),
    ],
)
def test_pass_fail_and_numbers_match_running_evaluator(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    core_mixing, m12_np, expected_predicted, expected_ratio = _expected_cpv_numbers(
        couplings
    )

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(expected_predicted)
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["m12_np_gev"] == pytest.approx(m12_np)
    assert abs(m12_np) == pytest.approx(core_mixing.abs_m12_np)
    assert result.diagnostics["abs_m12_np_from_magnitude_evaluator_gev"] == pytest.approx(
        core_mixing.abs_m12_np
    )
    assert result.diagnostics["m12_normalization_budget_gev"] == pytest.approx(
        constraint.anchor.m12_budget_gev
    )
    if expected_pass:
        assert result.ratio <= 1.0
    else:
        assert result.ratio > 10.0


def test_real_negative_m12_is_not_cp_odd():
    couplings = _uc_couplings(left=1.0e-4, right=0.7e-4)
    result = fcc.get(_PID).evaluate(point_builder.build_from_quark_couplings(couplings))

    assert result.diagnostics["m12_np_gev"].real < 0.0
    assert result.diagnostics["im_m12_np_gev"] == pytest.approx(0.0, abs=1e-30)
    assert abs(abs(result.diagnostics["m12_np_raw_argument_degrees"]) - 180.0) < 1e-12
    assert result.predicted == pytest.approx(0.0, abs=1e-30)
    assert result.ratio == pytest.approx(0.0, abs=1e-30)


def test_evaluate_is_pure_and_deterministic():
    couplings = _uc_couplings(
        left=1.0e-4 + 0.5e-4j,
        right=0.7e-4 + 0.2e-4j,
    )
    before_left_up = couplings.left_up.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_up, before_left_up)
