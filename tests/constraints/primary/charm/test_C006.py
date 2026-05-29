"""Production tests for C006 (D0 -> e+- mu-+ LFV)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.rare_charm_lfv_dilepton import (
    rare_charm_lfv_proxy_input,
)
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.model import RotationParameters, ckm_like_unitary

_PID = "C006"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charm" / "C006.yaml"


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _uc_couplings(
    left: complex,
    right: complex = 0.0j,
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


def _lfv_proxy(
    left_emu_overlap: complex,
    right_emu_overlap: complex = 0.0j,
    *,
    m_kk_gev: float = 3000.0,
):
    return rare_charm_lfv_proxy_input(
        left_emu_overlap,
        right_emu_overlap,
        m_kk_gev,
        source="C006 test lepton proxy",
    )


def _point(
    quark: QuarkMassBasisCouplings,
    lepton,
    *,
    kk_ew_mass_gev: float | None = None,
):
    extras = {
        "quark_mass_basis_couplings": quark,
        "lepton_mass_basis_couplings": lepton,
    }
    if kk_ew_mass_gev is not None:
        extras["kk_ew_mass_gev"] = kk_ew_mass_gev
    return point_builder.make_point(**extras)


def _manual_lfv_prediction(
    couplings: QuarkMassBasisCouplings,
    lepton_proxy,
    inputs,
    *,
    m_kk_gev: float | None = None,
    charge_state_factor: float = 2.0,
) -> float:
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=inputs.theta12,
            theta13=inputs.theta13,
            theta23=inputs.theta23,
            delta=inputs.delta,
        )
    )
    lambda_b = complex(np.conjugate(matrix[1, 2]) * matrix[0, 2])
    resolved_mkk = couplings.M_KK if m_kk_gev is None else float(m_kk_gev)
    g_weak = math.sqrt(4.0 * math.pi * inputs.alpha_em_mz / inputs.sin2_theta_w)
    g_z = g_weak / math.sqrt(1.0 - inputs.sin2_theta_w)
    neutral_delta = g_z / 2.0
    lepton_left = g_z * lepton_proxy.left_emu_overlap
    lepton_right = g_z * lepton_proxy.right_emu_overlap
    lepton_vector = lepton_left + lepton_right
    lepton_axial = lepton_right - lepton_left
    left_overlap = couplings.left_up[0, 1] / couplings.g_s
    right_overlap = couplings.right_up[0, 1] / couplings.g_s
    prefactor_wilson = -math.pi / (
        math.sqrt(2.0)
        * inputs.gf_gev_minus2
        * inputs.alpha_em_mz
        * lambda_b
        * resolved_mkk**2
    )
    c9 = prefactor_wilson * neutral_delta * left_overlap * lepton_vector
    c10 = prefactor_wilson * neutral_delta * left_overlap * lepton_axial
    c9p = prefactor_wilson * neutral_delta * right_overlap * lepton_vector
    c10p = prefactor_wilson * neutral_delta * right_overlap * lepton_axial
    c9_combo = c9 - c9p
    c10_combo = c10 - c10p

    meson = inputs.d0
    electron = inputs.electron
    muon = inputs.muon
    r_e = electron.mass_gev / meson.meson_mass_gev
    r_mu = muon.mass_gev / meson.meson_mass_gev
    x = r_e * r_e
    y = r_mu * r_mu
    phase_lambda = 1.0 + x * x + y * y - 2.0 * (x + y + x * y)
    mass_sum = r_e + r_mu
    mass_diff = r_mu - r_e
    vector_term = (1.0 - mass_sum**2) * abs(mass_diff * c9_combo) ** 2
    axial_term = (1.0 - mass_diff**2) * abs(mass_sum * c10_combo) ** 2
    tau_gev_inv = meson.lifetime_ps * 1.0e-12 / inputs.hbar_gev_s
    return float(
        charge_state_factor
        * tau_gev_inv
        * inputs.gf_gev_minus2**2
        * inputs.alpha_em_mz**2
        / (64.0 * math.pi**3)
        * meson.decay_constant_gev**2
        * meson.meson_mass_gev**3
        * abs(lambda_b) ** 2
        * math.sqrt(phase_lambda)
        * (vector_term + axial_term)
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charm"
    assert constraint.observable == "BR(D0 -> e+- mu-+) LFV"


def test_anchor_matches_yaml_and_budget():
    constraint = fcc.get(_PID)
    data = _yaml()
    pdg = data["pdg_or_equivalent"]
    current = pdg["canonical_current_limit"]
    lhcb = pdg["lhcb_primary_result"]
    belle = pdg["previous_limits"]["belle_2010"]
    babar = pdg["previous_limits"]["babar_2012"]
    rs = data["paper_era_reference"]["rs_baseline"]

    assert constraint.anchor.current_limit.value == pytest.approx(current["value"])
    assert constraint.anchor.current_limit.confidence_level == pytest.approx(
        current["confidence_level"]
    )
    assert constraint.anchor.current_limit.source_url == current["source_url"]
    assert constraint.anchor.lhcb_primary_limit.value == pytest.approx(lhcb["value"])
    assert constraint.anchor.lhcb_primary_limit.integrated_luminosity_fb_inv == (
        pytest.approx(lhcb["integrated_luminosity_fb_inv"])
    )
    assert constraint.anchor.lhcb_primary_limit.collision_energy_tev == tuple(
        float(item) for item in lhcb["collision_energy_TeV"]
    )
    assert constraint.anchor.belle_previous_limit.value == pytest.approx(belle["value"])
    assert constraint.anchor.babar_previous_limit.value == pytest.approx(babar["value"])
    assert constraint.anchor.rs_baseline.generic_rs_kk_gluon_scale_tev == pytest.approx(
        rs["generic_rs_kk_gluon_scale_TeV"]
    )
    assert (
        constraint.anchor.rs_baseline.composite_pseudo_goldstone_kk_gluon_scale_tev
        == pytest.approx(rs["composite_pseudo_goldstone_kk_gluon_scale_TeV"])
    )
    assert constraint.anchor.budget == pytest.approx(current["value"])
    assert constraint.anchor.budget == pytest.approx(1.3e-8)


def test_anchor_loading_fails_loudly_for_missing_candidate_or_value():
    with pytest.raises(fcc.AnchorError, match="none of the expected anchor keys"):
        fcc.load_anchor(_PID, family="charm", candidates=("missing_c006_anchor",))
    with pytest.raises(fcc.AnchorError, match="has no 'missing_limit' field"):
        fcc.load_anchor(
            _PID,
            family="charm",
            candidates=("canonical_current_limit",),
            value_key="missing_limit",
        )


def test_evaluate_without_required_inputs_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extras"] == (
        "quark_mass_basis_couplings",
        "lepton_mass_basis_couplings",
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_invalid_lepton_proxy_is_unevaluated_not_real_pass():
    quark = _uc_couplings(left=1.0)
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(
            quark_mass_basis_couplings=quark,
            lepton_mass_basis_couplings={},
        )
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["invalid_extra"] == (
        "quark_mass_basis_couplings",
        "lepton_mass_basis_couplings",
    )
    assert result.diagnostics["exception_type"] == "KeyError"


def test_zero_np_prediction_and_sm_lfv_policy():
    constraint = fcc.get(_PID)
    quark = _uc_couplings(left=0.0j, right=0.0j)
    lepton = _lfv_proxy(left_emu_overlap=0.5, right_emu_overlap=0.2j)
    result = constraint.evaluate(_point(quark, lepton))

    assert constraint.sm_result.branching_fraction == pytest.approx(0.0)
    assert result.predicted == pytest.approx(0.0)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.diagnostics["sm_lfv_policy"].startswith("D0 -> e mu is charged-LFV")
    assert result.diagnostics["charge_conjugate_modes_included"] is True
    assert result.passes is True


def test_lfv_proxy_prediction_matches_independent_manual_formula():
    constraint = fcc.get(_PID)
    quark = _uc_couplings(left=2.0e-2 + 0.3e-2j, right=0.4e-2j)
    lepton = _lfv_proxy(left_emu_overlap=0.40 + 0.10j, right_emu_overlap=0.05j)
    result = constraint.evaluate(_point(quark, lepton))
    manual = _manual_lfv_prediction(quark, lepton, constraint.sm_inputs)

    assert result.predicted == pytest.approx(manual, rel=1e-12, abs=0.0)
    assert result.predicted == pytest.approx(
        3.752990413820220023e-13,
        rel=1e-12,
        abs=0.0,
    )
    assert result.diagnostics["branching_fraction_one_charge"] == pytest.approx(
        manual / 2.0,
        rel=1e-12,
        abs=0.0,
    )
    assert result.diagnostics["charge_state_factor"] == pytest.approx(2.0)
    assert result.ratio == pytest.approx(result.predicted / constraint.anchor.budget)
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["rs_matching_assumption"]


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    quark = _uc_couplings(left=1.0 + 0.2j, right=0.5j)
    lepton = _lfv_proxy(left_emu_overlap=0.3 + 0.1j, right_emu_overlap=0.2j)
    result = fcc.get(_PID).evaluate(_point(quark, lepton))

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
        "left_uc_coupling",
        "right_uc_coupling",
        "lambda_b",
        "lepton_left_delta_emu",
        "lepton_right_delta_emu",
        "c9_lfv_combination",
        "c10_lfv_combination",
        "c9_lfv_np",
        "c10_lfv_np",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "electron_mass_gev",
        "muon_mass_gev",
        "phase_space_lambda_sqrt",
        "np_shift_branching_fraction",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["wilson_coefficients"]


@pytest.mark.parametrize(
    ("quark", "lepton", "expected_pass"),
    [
        (_uc_couplings(left=0.10), _lfv_proxy(0.10), True),
        (_uc_couplings(left=50.0), _lfv_proxy(5.0), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    quark: QuarkMassBasisCouplings,
    lepton,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    result = constraint.evaluate(_point(quark, lepton))

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    quark = _uc_couplings(left=1.0)
    lepton = _lfv_proxy(left_emu_overlap=0.5)
    default_result = fcc.get(_PID).evaluate(_point(quark, lepton))
    heavy_result = fcc.get(_PID).evaluate(
        _point(quark, lepton, kk_ew_mass_gev=6000.0)
    )

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert heavy_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert heavy_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert abs(heavy_result.diagnostics["c10_lfv_combination"]) == pytest.approx(
        abs(default_result.diagnostics["c10_lfv_combination"]) / 4.0
    )
    assert heavy_result.predicted == pytest.approx(default_result.predicted / 16.0)


def test_evaluate_is_pure_and_deterministic():
    quark = _uc_couplings(left=1.0 + 0.2j, right=0.5j)
    lepton = _lfv_proxy(left_emu_overlap=0.3 + 0.1j, right_emu_overlap=0.2j)
    before_left_up = quark.left_up.copy()
    before_right_up = quark.right_up.copy()
    point = _point(quark, lepton)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(quark.left_up, before_left_up)
    np.testing.assert_array_equal(quark.right_up, before_right_up)
    assert point.get_extra("lepton_mass_basis_couplings") == lepton
