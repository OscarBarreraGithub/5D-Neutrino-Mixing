"""Production tests for L001 (mu -> e gamma)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.lepton import (
    mu_to_e_gamma_proxy_input,
)

_PID = "L001"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charged_lepton" / "L001.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _rotation_pmns() -> np.ndarray:
    theta = 0.5
    c = math.cos(theta)
    s = math.sin(theta)
    return np.asarray(
        [
            [c, s, 0.0],
            [-s, c, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=complex,
    )


def _proxy_point(
    y_n_bar: tuple[float, float, float],
    *,
    m_kk_gev: float = 3000.0,
):
    proxy = mu_to_e_gamma_proxy_input(
        y_n_bar,
        _rotation_pmns(),
        m_kk_gev,
        source="L001 test proxy",
    )
    return point_builder.make_point(lepton_mass_basis_couplings=proxy)


def _manual_branching_fraction(
    y_n_bar: tuple[float, float, float],
    *,
    prefactor_br: float,
    m_kk_gev: float = 3000.0,
    reference_scale_gev: float = 3000.0,
) -> tuple[float, float, complex]:
    y = np.asarray(y_n_bar, dtype=complex)
    pmns = _rotation_pmns()
    y_matrix = pmns @ np.diag(y)
    product = y_matrix @ y_matrix.conj().T
    off_diagonal = complex(product[0, 1])
    lhs = float(abs(off_diagonal))
    branching_fraction = float(
        prefactor_br * lhs * lhs * (reference_scale_gev / m_kk_gev) ** 4
    )
    return branching_fraction, lhs, off_diagonal


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charged_lepton"
    assert constraint.observable == "BR(mu -> e gamma)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["primary_current_limit"]
    repo = pdg["repo_default"]

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.budget == pytest.approx(1.5e-13)
    assert constraint.anchor.repo_default_br_limit.value == pytest.approx(
        repo["br_limit"]["value"]
    )
    assert constraint.anchor.prefactor_br.value == pytest.approx(
        repo["prefac_br"]["value"]
    )
    assert constraint.anchor.lfv_c.value == pytest.approx(repo["lfv_C"]["value"])
    assert constraint.anchor.lfv_c.value == pytest.approx(
        math.sqrt(exp["value"] / repo["prefac_br"]["value"])
    )
    assert constraint.anchor.c_paper.value == pytest.approx(
        repo["c_paper"]["value"]
    )

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
    assert result.notes == (
        "NOT EVALUATED — no lepton dipole prediction available "
        "(lepton-sector RS couplings not on ParameterPoint)"
    )
    assert "pass" not in result.notes.lower()
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["unevaluated_reason"] == (
        "no lepton dipole prediction available "
        "(lepton-sector RS couplings not on ParameterPoint)"
    )
    assert "non-vetoing only" in result.diagnostics["passes_semantics"]
    assert result.diagnostics["missing_extra"] == "lepton_mass_basis_couplings"
    assert "needs_human_physics" not in result.diagnostics


def test_invalid_lepton_input_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(lepton_mass_basis_couplings={"y_n_bar": [1.0, 2.0]})
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.sm_prediction is None
    assert result.notes.startswith("NOT EVALUATED —")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["invalid_extra"] == "lepton_mass_basis_couplings"
    assert result.diagnostics["exception_type"] in {"KeyError", "ValueError"}
    assert "needs_human_physics" not in result.diagnostics


def test_proxy_numerics_match_independent_dipole_recomputation():
    constraint = fcc.get(_PID)
    y_n_bar = (0.10, 0.20, 0.30)
    result = constraint.evaluate(_proxy_point(y_n_bar))
    expected_br, expected_lhs, expected_off = _manual_branching_fraction(
        y_n_bar,
        prefactor_br=constraint.anchor.prefactor_br.value,
    )

    assert result.predicted == pytest.approx(expected_br)
    assert result.ratio == pytest.approx(expected_br / constraint.anchor.budget)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.diagnostics["dipole_lhs"] == pytest.approx(expected_lhs)
    assert result.diagnostics["off_diagonal_12"] == pytest.approx(expected_off)
    assert result.diagnostics["lfv_coefficient"] == pytest.approx(
        constraint.anchor.lfv_c.value
    )
    assert result.diagnostics["prefactor_br"] == pytest.approx(
        constraint.anchor.prefactor_br.value
    )
    assert result.diagnostics["evaluated"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


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
    assert isinstance(result.diagnostics["off_diagonal_12"], complex)
    assert isinstance(result.diagnostics["product_matrix"][0][1], complex)
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
    assert result.diagnostics["input_kind"] == "MuToEGammaProxyInput"


@pytest.mark.parametrize(
    ("point", "expected_pass"),
    [
        (_proxy_point((0.01, 0.02, 0.03)), True),
        (_proxy_point((0.10, 0.20, 0.30)), False),
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
    proxy = mu_to_e_gamma_proxy_input(
        (0.10, 0.20, 0.30),
        _rotation_pmns(),
        3000.0,
        source="L001 test proxy",
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
    proxy = mu_to_e_gamma_proxy_input(
        (0.03, 0.07, 0.11),
        _rotation_pmns(),
        3000.0,
        source="L001 test proxy",
    )
    point = point_builder.make_point(lepton_mass_basis_couplings=proxy)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.get_extra("lepton_mass_basis_couplings") == proxy
