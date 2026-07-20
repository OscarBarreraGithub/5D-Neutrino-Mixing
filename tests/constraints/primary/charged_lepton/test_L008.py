"""Production tests for L008 (tau -> e gamma)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.lepton_tau_e import (
    tau_to_e_gamma_from_lepton_input,
    tau_to_e_gamma_proxy_input,
)
from flavorConstraints.muToEGamma import check_mu_to_e_gamma_raw
from quarkConstraints.lfv_three_body import TAU_TO_E_NUNU_BRANCHING_FRACTION

_PID = "L008"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charged_lepton" / "L008.yaml"
_L001_SIDECAR = (
    _REPO_ROOT / "flavor_catalog" / "processes" / "charged_lepton" / "L001.yaml"
)


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _l001_yaml_pdg_block():
    with open(_L001_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _tau_e_rotation_pmns() -> np.ndarray:
    theta = 0.5
    c = math.cos(theta)
    s = math.sin(theta)
    return np.asarray(
        [
            [c, 0.0, s],
            [0.0, 1.0, 0.0],
            [-s, 0.0, c],
        ],
        dtype=complex,
    )


def _proxy_point(
    y_n_bar: tuple[float, float, float],
    *,
    m_kk_gev: float = 3000.0,
):
    proxy = tau_to_e_gamma_proxy_input(
        y_n_bar,
        _tau_e_rotation_pmns(),
        m_kk_gev,
        source="L008 test proxy",
    )
    return point_builder.make_point(lepton_mass_basis_couplings=proxy)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charged_lepton"
    assert constraint.observable == "BR(tau -> e gamma)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    l001 = _l001_yaml_pdg_block()
    exp = pdg["primary_current_limit"]
    babar = pdg["primary_experiment"]
    belle = pdg["recent_belle_cross_check"]
    belle_ii = pdg["belle_ii_projection"]
    prefactor = l001["repo_default"]["prefac_br"]

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.primary_experiment.value == pytest.approx(babar["value"])
    assert constraint.anchor.primary_experiment.source_url == babar["source_url"]
    assert constraint.anchor.belle_cross_check.value == pytest.approx(belle["value"])
    assert constraint.anchor.belle_cross_check.source_url == belle["source_url"]
    assert constraint.anchor.belle_ii_projection.value == pytest.approx(
        belle_ii["projected_limit"]
    )
    assert constraint.anchor.belle_ii_projection.source_url == belle_ii["source_url"]
    assert constraint.anchor.dipole_prefactor_br.value == pytest.approx(
        prefactor["value"]
    )
    assert constraint.anchor.budget == pytest.approx(exp["value"])
    assert constraint.anchor.budget == pytest.approx(3.3e-8)

    with pytest.raises(anchors.AnchorError):
        anchors.load_anchor(
            _PID,
            family="charged_lepton",
            candidates=("no_such_block",),
        )


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.sm_prediction is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["unevaluated_reason"] == (
        "no tau->e gamma dipole prediction available "
        "(lepton-sector RS couplings not on ParameterPoint)"
    )
    assert "non-vetoing only" in result.diagnostics["passes_semantics"]
    assert result.diagnostics["missing_extra"] == "lepton_mass_basis_couplings"
    assert "needs_human_physics" not in result.diagnostics


def test_invalid_lepton_input_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(lepton_mass_basis_couplings={"y_n_bar": [1.0, 2.0]})
    )

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.sm_prediction is None
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["invalid_extra"] == "lepton_mass_basis_couplings"
    assert result.diagnostics["exception_type"] in {"KeyError", "ValueError"}
    assert "needs_human_physics" not in result.diagnostics


def test_tau_to_e_gamma_rejects_caller_flavor_override():
    constraint = fcc.get(_PID)
    bad = {
        "initial_flavor": "mu",
        "final_flavor": "e",
        "y_n_bar": [0.1, 0.2, 0.3],
        "pmns": _tau_e_rotation_pmns(),
        "m_kk_gev": 3000.0,
        "source": "bad mu->e probe",
    }

    with pytest.raises(ValueError, match="pinned to initial_flavor='tau'"):
        tau_to_e_gamma_from_lepton_input(
            bad,
            br_limit=constraint.anchor.budget,
            prefactor_br=constraint.anchor.dipole_prefactor_br.value,
        )

    result = constraint.evaluate(point_builder.make_point(lepton_mass_basis_couplings=bad))
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["exception_type"] == "ValueError"
    assert "pinned to initial_flavor='tau'" in result.diagnostics["exception"]


def test_proxy_numerics_match_mu_to_e_core_with_tau_e_rotation():
    constraint = fcc.get(_PID)
    y_n_bar = (0.10, 0.20, 0.30)
    m_kk_gev = 3000.0
    reference_scale_gev = 3000.0
    result = constraint.evaluate(_proxy_point(y_n_bar))
    effective_prefactor = (
        constraint.anchor.dipole_prefactor_br.value
        * TAU_TO_E_NUNU_BRANCHING_FRACTION
    )

    c_lfv = math.sqrt(constraint.anchor.budget / effective_prefactor)
    core = check_mu_to_e_gamma_raw(
        np.asarray(y_n_bar, dtype=complex),
        _tau_e_rotation_pmns()[[0, 2, 1], :],
        M_KK=m_kk_gev,
        C=c_lfv,
        reference_scale=reference_scale_gev,
    )
    expected_br = float(
        effective_prefactor
        * float(core["lhs"]) ** 2
        * (reference_scale_gev / m_kk_gev) ** 4
    )
    expected_ratio = expected_br / constraint.anchor.budget
    expected_off_13 = complex(core["off_diagonal_12"])
    expected_off_12 = complex(core["product_matrix"][0, 2])

    assert result.predicted == pytest.approx(expected_br)
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.diagnostics["dipole_lhs"] == pytest.approx(core["lhs"])
    assert result.diagnostics["dipole_rhs"] == pytest.approx(core["rhs"])
    assert result.diagnostics["dipole_ratio_to_bound"] == pytest.approx(
        core["ratio"]
    )
    assert result.diagnostics["off_diagonal_13"] == pytest.approx(expected_off_13)
    assert result.diagnostics["product_matrix"][0][2] == pytest.approx(expected_off_13)
    assert result.diagnostics["product_matrix"][0][1] == pytest.approx(expected_off_12)
    assert result.diagnostics["core_off_diagonal_12_after_permutation"] == (
        pytest.approx(expected_off_13)
    )
    assert abs(expected_off_13) > 0.0
    assert expected_off_12 == pytest.approx(0.0j)
    assert result.diagnostics["lfv_coefficient"] == pytest.approx(c_lfv)
    assert result.diagnostics["prefactor_br"] == pytest.approx(effective_prefactor)
    assert result.diagnostics["muon_normalized_prefactor_br"] == pytest.approx(
        constraint.anchor.dipole_prefactor_br.value
    )
    assert result.diagnostics["tau_leptonic_branching_fraction"] == pytest.approx(
        TAU_TO_E_NUNU_BRANCHING_FRACTION
    )
    assert result.diagnostics["evaluated"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["row_permutation"] == (0, 2, 1)
    assert result.diagnostics["final_initial_indices_zero_based"] == (0, 2)


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(_proxy_point((0.03, 0.07, 0.11)))

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
    assert isinstance(result.diagnostics["off_diagonal_13"], complex)
    assert isinstance(result.diagnostics["off_diagonal_31"], complex)
    assert isinstance(result.diagnostics["product_matrix"][0][2], complex)
    for key in (
        "dipole_lhs",
        "dipole_rhs",
        "dipole_ratio_to_bound",
        "lfv_coefficient",
        "prefactor_br",
        "m_kk_gev",
        "reference_scale_gev",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["used_proxy"] is True
    assert result.diagnostics["input_kind"] == "TauToEGammaProxyInput"


@pytest.mark.parametrize(
    ("point", "expected_pass"),
    [
        (_proxy_point((0.01, 0.02, 0.03)), True),
        (_proxy_point((0.10, 1.00, 6.00)), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(point, expected_pass: bool):
    result = fcc.get(_PID).evaluate(point)

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_changes_dipole_scaling():
    proxy = tau_to_e_gamma_proxy_input(
        (0.10, 0.20, 0.30),
        _tau_e_rotation_pmns(),
        3000.0,
        source="L008 test proxy",
    )
    default_point = point_builder.make_point(lepton_mass_basis_couplings=proxy)
    heavy_point = point_builder.make_point(
        lepton_mass_basis_couplings=proxy,
        kk_ew_mass_gev=6000.0,
    )

    default_result = fcc.get(_PID).evaluate(default_point)
    heavy_result = fcc.get(_PID).evaluate(heavy_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert heavy_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert heavy_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert heavy_result.predicted == pytest.approx(default_result.predicted / 16.0)


def test_evaluate_is_pure_and_deterministic():
    proxy = tau_to_e_gamma_proxy_input(
        (0.03, 0.07, 0.11),
        _tau_e_rotation_pmns(),
        3000.0,
        source="L008 test proxy",
    )
    point = point_builder.make_point(lepton_mass_basis_couplings=proxy)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.get_extra("lepton_mass_basis_couplings") == proxy
