"""Production tests for L005 (mu->e conversion in titanium)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from quarkConstraints.mu_e_conversion import (
    mu_e_conversion_from_components,
    titanium_nuclear_inputs,
    zero_mu_e_conversion_coefficients,
)

_PID = "L005"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charged_lepton" / "L005.yaml"
_L001_SIDECAR = (
    _REPO_ROOT / "flavor_catalog" / "processes" / "charged_lepton" / "L001.yaml"
)

_G_F = 1.1663787e-5
_HBAR_GEV_S = 6.582119569e-25
_M_MU_GEV = 0.1056583755
_D_TI = 0.0864
_VP_TI = 0.0396
_VN_TI = 0.0468
_SP_TI = 0.0368
_SN_TI = 0.0435
_CAPTURE_TI_S_INV = 2.59e6


def _yaml_pdg_block(path: Path = _SIDECAR):
    with open(path) as handle:
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


def _point_from_mapping(mapping):
    return point_builder.make_point(lepton_mass_basis_couplings=dict(mapping))


def _manual_l001_dipole_br(
    y_n_bar: tuple[float, float, float],
    *,
    prefactor_br: float,
    m_kk_gev: float = 3000.0,
    reference_scale_gev: float = 3000.0,
) -> tuple[float, complex]:
    y = np.asarray(y_n_bar, dtype=complex)
    pmns = _rotation_pmns()
    y_matrix = pmns @ np.diag(y)
    product = y_matrix @ y_matrix.conj().T
    off_diagonal = complex(product[0, 1])
    lhs = float(abs(off_diagonal))
    return (
        float(prefactor_br * lhs * lhs * (reference_scale_gev / m_kk_gev) ** 4),
        off_diagonal,
    )


def _manual_conversion_rate(
    *,
    dipole_parent_br: float,
    g_lv_p: complex = 0.0j,
    g_lv_n: complex = 0.0j,
    g_rv_p: complex = 0.0j,
    g_rv_n: complex = 0.0j,
    g_ls_p: complex = 0.0j,
    g_ls_n: complex = 0.0j,
    g_rs_p: complex = 0.0j,
    g_rs_n: complex = 0.0j,
) -> tuple[float, float, complex, complex]:
    scalar_left = g_ls_p * _SP_TI + g_ls_n * _SN_TI
    scalar_right = g_rs_p * _SP_TI + g_rs_n * _SN_TI
    vector_left = g_lv_p * _VP_TI + g_lv_n * _VN_TI
    vector_right = g_rv_p * _VP_TI + g_rv_n * _VN_TI
    contact_left = complex(scalar_left + vector_left)
    contact_right = complex(scalar_right + vector_right)
    contact_inner = float(abs(contact_left) ** 2 + abs(contact_right) ** 2)
    dipole_norm = math.sqrt(dipole_parent_br / (384.0 * math.pi**2))
    dipole_inner = float(_D_TI**2 * dipole_norm**2)
    interference_inner = float(2.0 * _D_TI * dipole_norm * math.sqrt(contact_inner))
    lower_inner = max(0.0, dipole_inner + contact_inner - interference_inner)
    upper_inner = dipole_inner + contact_inner + interference_inner
    capture_gev = _CAPTURE_TI_S_INV * _HBAR_GEV_S
    kko_dimension_factor = _M_MU_GEV**5
    return (
        float(2.0 * _G_F**2 * kko_dimension_factor * lower_inner / capture_gev),
        float(2.0 * _G_F**2 * kko_dimension_factor * upper_inner / capture_gev),
        contact_left,
        contact_right,
    )


def test_pure_dipole_kko_normalization_benchmark():
    result = mu_e_conversion_from_components(
        dipole_parent_branching_fraction=1.0,
        coefficients=zero_mu_e_conversion_coefficients(),
        nuclear_inputs=titanium_nuclear_inputs(),
    )

    assert result.conversion_rate == pytest.approx(0.004139613411251543)
    assert result.dipole_component == pytest.approx(result.conversion_rate)
    assert result.diagnostics["kko_overlap_dimension_factor_gev5"] == pytest.approx(
        _M_MU_GEV**5
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charged_lepton"
    assert constraint.observable == "CR(mu Ti -> e Ti)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    primary = pdg["primary_current_limit"]
    publication = pdg["experimental_publication"]
    l001 = _yaml_pdg_block(_L001_SIDECAR)

    assert constraint.anchor.value == pytest.approx(primary["value"])
    assert constraint.anchor.budget == pytest.approx(primary["value"])
    assert constraint.anchor.budget == pytest.approx(6.1e-13)
    assert constraint.anchor.target_material == primary["target_material"] == "Ti"
    assert constraint.anchor.relation == primary["relation"]
    assert constraint.anchor.confidence_level == primary["confidence_level"]
    assert constraint.anchor.experiment == "SINDRUM II"
    assert constraint.anchor.source_url == primary["source_url"]
    assert constraint.anchor.pdg_note == primary["pdg_note"]
    assert constraint.anchor.experimental_publication.doi == publication["doi"]
    assert constraint.anchor.dipole_br_limit.value == pytest.approx(
        l001["primary_current_limit"]["value"]
    )
    assert constraint.anchor.dipole_prefactor_br.value == pytest.approx(
        l001["repo_default"]["prefac_br"]["value"]
    )
    assert constraint.nuclear_inputs.target == "Ti"
    assert constraint.nuclear_inputs.isotope == "48Ti"
    assert constraint.nuclear_inputs.D == pytest.approx(_D_TI)
    assert constraint.nuclear_inputs.V_p == pytest.approx(_VP_TI)
    assert constraint.nuclear_inputs.V_n == pytest.approx(_VN_TI)
    assert constraint.nuclear_inputs.S_p == pytest.approx(_SP_TI)
    assert constraint.nuclear_inputs.S_n == pytest.approx(_SN_TI)
    assert constraint.nuclear_inputs.capture_rate_s_inv == pytest.approx(
        _CAPTURE_TI_S_INV
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
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.notes.startswith("NOT EVALUATED --")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "lepton_mass_basis_couplings"
    assert result.diagnostics["budget_limit_status"] == "observed_experimental_bound"
    assert "needs_human_physics" not in result.diagnostics


def test_invalid_lepton_input_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(lepton_mass_basis_couplings={"source": "empty"})
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.notes.startswith("NOT EVALUATED --")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["invalid_extra"] == "lepton_mass_basis_couplings"
    assert result.diagnostics["exception_type"] == "TypeError"


def test_absent_vector_ignores_legacy_vector_but_preserves_scalar_dipole():
    constraint = fcc.get(_PID)
    y_n_bar = (1.0e-4, 2.0e-4, 3.0e-4)
    lepton = {
        "y_n_bar": y_n_bar,
        "pmns": _rotation_pmns(),
        "m_kk_gev": 3000.0,
        "g_lv_p": 2.0e-12 + 0.5e-12j,
        "g_lv_n": -1.0e-12j,
        "g_ls_p": 1.5e-12,
        "g_rs_n": -0.4e-12j,
        "source": "L005 test low-energy proxy",
    }
    result = constraint.evaluate(_point_from_mapping(lepton))
    expected_dipole_br, expected_off = _manual_l001_dipole_br(
        y_n_bar,
        prefactor_br=constraint.anchor.dipole_prefactor_br.value,
    )
    expected_cr_lower, expected_cr_upper, expected_left, expected_right = (
        _manual_conversion_rate(
            dipole_parent_br=expected_dipole_br,
            g_ls_p=lepton["g_ls_p"],
            g_rs_n=lepton["g_rs_n"],
        )
    )

    assert result.predicted == pytest.approx(expected_cr_lower)
    assert result.ratio == pytest.approx(expected_cr_lower / constraint.anchor.budget)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.diagnostics["conversion_rate_lower"] == pytest.approx(
        expected_cr_lower
    )
    assert result.diagnostics["conversion_rate_upper"] == pytest.approx(
        expected_cr_upper
    )
    assert result.diagnostics["conversion_rate_interval"] == pytest.approx(
        (expected_cr_lower, expected_cr_upper)
    )
    assert result.diagnostics["dipole_parent_branching_fraction"] == pytest.approx(
        expected_dipole_br
    )
    assert result.diagnostics["dipole_off_diagonal_12"] == pytest.approx(
        expected_off
    )
    assert result.diagnostics["contact_left_nuclear_amplitude"] == pytest.approx(
        expected_left
    )
    assert result.diagnostics["contact_right_nuclear_amplitude"] == pytest.approx(
        expected_right
    )
    assert result.diagnostics["vector_component"] == pytest.approx(0.0)
    assert result.diagnostics["legacy_vector_proxy_ignored"] is True
    assert result.diagnostics["vector_tree_missing_extra"] == "rs_ew_couplings"
    assert result.diagnostics["target"] == "Ti"
    assert result.diagnostics["dipole_input_present"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["dipole_contact_interference_treatment"] == (
        "unknown_relative_phase_interval_NEEDS-HUMAN-PHYSICS"
    )
    assert result.diagnostics["dipole_contact_relative_phase_status"] == (
        "NEEDS-HUMAN-PHYSICS"
    )
    assert result.diagnostics["conversion_rate_verdict_branch"] == (
        "lower_envelope_unknown_relative_phase"
    )


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(
        _point_from_mapping(
            {
                "g_lv_p": 1.0e-12 + 0.5e-12j,
                "g_rv_n": -0.3e-12j,
                "g_ls_p": 0.2e-12,
                "g_rs_n": -0.1e-12,
                "source": "L005 finite field probe",
            }
        )
    )

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
        "g_lv_p",
        "g_rv_n",
        "contact_left_nuclear_amplitude",
        "contact_right_nuclear_amplitude",
        "left_nuclear_amplitude",
        "right_nuclear_amplitude",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "dipole_component",
        "scalar_component",
        "vector_component",
        "contact_component",
        "conversion_rate_lower",
        "conversion_rate_upper",
        "capture_rate_gev",
    ):
        value = result.diagnostics[key]
        assert isinstance(value, float)
        assert math.isfinite(value)
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["dipole_input_present"] is False
    assert result.diagnostics["vector_component"] == pytest.approx(0.0)
    assert result.diagnostics["legacy_vector_proxy_ignored"] is True
    assert result.diagnostics["vector_tree_missing_extra"] == "rs_ew_couplings"
    assert result.diagnostics["target"] == "Ti"
    assert result.diagnostics["target_material"] == "Ti"
    assert result.diagnostics["budget_limit_status"] == "observed_experimental_bound"
    assert result.diagnostics["budget_verdict_role"] == "HARD observed upper-limit veto"


@pytest.mark.parametrize(
    "lepton",
    [
        {"g_lv_p": 1.0e-13, "source": "safe stale L005 vector proxy"},
        {"g_lv_p": 1.0e-6, "source": "excluded stale L005 vector proxy"},
    ],
)
def test_legacy_vector_only_proxy_is_unevaluated_not_real_pass(lepton):
    result = fcc.get(_PID).evaluate(_point_from_mapping(lepton))

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["exception_type"] == "TypeError"


def test_evaluate_is_pure_and_deterministic():
    lepton = {
        "g_lv_p": 2.0e-12,
        "g_ls_n": -0.5e-12j,
        "source": "L005 purity proxy",
    }
    point = _point_from_mapping(lepton)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert first.diagnostics["vector_component"] == pytest.approx(0.0)
    assert first.diagnostics["legacy_vector_proxy_ignored"] is True
    assert point.get_extra("lepton_mass_basis_couplings") == lepton
