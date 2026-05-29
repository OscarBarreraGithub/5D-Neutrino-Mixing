"""Production tests for C008 (D+ -> pi+ e+- mu-+ LFV)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.rare_charm_lfv_semileptonic import (
    rare_charm_lfv_proxy_input,
)
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_charm_lfv_dilepton import compute_rare_charm_lfv_wilsons
from quarkConstraints.rare_charm_semileptonic import dtopi_fplus, dtopi_fzero

_PID = "C008"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charm" / "C008.yaml"
_METRIC = np.diag([1.0, -1.0, -1.0, -1.0])


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
        source="C008 test lepton proxy",
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


def _manual_kallen(a: float, b: float, c: float) -> float:
    return a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c)


def _manual_dtopi_emu_proxy(
    couplings: QuarkMassBasisCouplings,
    lepton_proxy,
    inputs,
    *,
    m_kk_gev: float | None = None,
    charge_mode: str = "eplus_muminus",
) -> float:
    wilsons = compute_rare_charm_lfv_wilsons(
        couplings,
        lepton_proxy,
        transition="c_u",
        m_kk_gev=m_kk_gev,
        inputs=inputs.rare_charm,
    )
    c9_semileptonic = wilsons.c9_lfv_np + wilsons.c9p_lfv_np
    c10_semileptonic = wilsons.c10_lfv_np + wilsons.c10p_lfv_np

    m_d = inputs.dplus_mass_gev
    m_pi = inputs.piplus_mass_gev
    m_e = inputs.rare_charm.electron.mass_gev
    m_mu = inputs.rare_charm.muon.mass_gev
    if charge_mode == "eplus_muminus":
        m_minus, m_plus = m_mu, m_e
    elif charge_mode == "eminus_muplus":
        m_minus, m_plus = m_e, m_mu
    else:
        raise ValueError(charge_mode)
    q2_min = (m_e + m_mu) ** 2
    q2_max = (m_d - m_pi) ** 2
    tau = inputs.dplus_lifetime_ps * 1.0e-12 / inputs.rare_charm.hbar_gev_s
    gamma0 = (
        inputs.rare_charm.gf_gev_minus2**2
        * inputs.rare_charm.alpha_em_mz**2
        * abs(wilsons.lambda_b) ** 2
        / ((4.0 * math.pi) ** 5 * m_d**3)
    )

    q_nodes, q_weights = np.polynomial.legendre.leggauss(inputs.quadrature_points)
    a_nodes, a_weights = np.polynomial.legendre.leggauss(48)
    mid = 0.5 * (q2_max + q2_min)
    half_width = 0.5 * (q2_max - q2_min)
    total = 0.0
    for q_node, q_weight in zip(q_nodes, q_weights):
        q2 = mid + half_width * float(q_node)
        lam_had = _manual_kallen(m_d * m_d, m_pi * m_pi, q2)
        lam_lep = _manual_kallen(q2, m_minus * m_minus, m_plus * m_plus)
        beta_lfv = math.sqrt(max(0.0, lam_lep)) / q2
        angular = 0.0
        for cos_theta, a_weight in zip(a_nodes, a_weights):
            vector, axial = _manual_tensor_contractions(
                q2,
                float(cos_theta),
                m_minus=m_minus,
                m_plus=m_plus,
                inputs=inputs,
            )
            angular += float(a_weight) * (
                abs(c9_semileptonic) ** 2 * vector
                + abs(c10_semileptonic) ** 2 * axial
            )
        total += (
            float(q_weight)
            * tau
            * gamma0
            * math.sqrt(max(0.0, lam_had))
            * beta_lfv
            * angular
            / 4.0
        )
    return float(half_width * total)


def _manual_tensor_contractions(
    q2: float,
    cos_theta: float,
    *,
    m_minus: float,
    m_plus: float,
    inputs,
) -> tuple[float, float]:
    root_q2 = math.sqrt(q2)
    m_d = inputs.dplus_mass_gev
    m_pi = inputs.piplus_mass_gev
    m_d2 = m_d * m_d
    m_pi2 = m_pi * m_pi
    delta = m_d2 - m_pi2
    lam_had = _manual_kallen(m_d2, m_pi2, q2)
    lam_lep = _manual_kallen(q2, m_minus * m_minus, m_plus * m_plus)
    p_had = math.sqrt(max(0.0, lam_had)) / (2.0 * root_q2)
    p_lep = math.sqrt(max(0.0, lam_lep)) / (2.0 * root_q2)
    e_d = (m_d2 - m_pi2 + q2) / (2.0 * root_q2)
    e_pi = (m_d2 - m_pi2 - q2) / (2.0 * root_q2)
    e_minus = (q2 + m_minus * m_minus - m_plus * m_plus) / (2.0 * root_q2)
    e_plus = (q2 + m_plus * m_plus - m_minus * m_minus) / (2.0 * root_q2)
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))

    p_d = np.array([e_d, 0.0, 0.0, p_had], dtype=np.complex128)
    p_pi = np.array([e_pi, 0.0, 0.0, p_had], dtype=np.complex128)
    q_vec = np.array([root_q2, 0.0, 0.0, 0.0], dtype=np.complex128)
    p_minus = np.array(
        [e_minus, p_lep * sin_theta, 0.0, p_lep * cos_theta],
        dtype=np.complex128,
    )
    p_plus = np.array(
        [e_plus, -p_lep * sin_theta, 0.0, -p_lep * cos_theta],
        dtype=np.complex128,
    )
    fplus = dtopi_fplus(q2, inputs)
    fzero = dtopi_fzero(q2, inputs)
    hadronic = fplus * (p_d + p_pi - delta / q2 * q_vec) + fzero * delta / q2 * q_vec
    hadronic_cov = _METRIC @ hadronic
    p_dot = float(np.real(p_minus @ _METRIC @ p_plus))
    vector_tensor = _manual_lepton_tensor(p_minus, p_plus, p_dot + m_minus * m_plus)
    axial_tensor = _manual_lepton_tensor(p_minus, p_plus, p_dot - m_minus * m_plus)
    return (
        _manual_contract(hadronic_cov, vector_tensor),
        _manual_contract(hadronic_cov, axial_tensor),
    )


def _manual_lepton_tensor(
    p_minus: np.ndarray,
    p_plus: np.ndarray,
    metric_term: float,
) -> np.ndarray:
    tensor = np.zeros((4, 4), dtype=np.complex128)
    for mu in range(4):
        for nu in range(4):
            tensor[mu, nu] = 4.0 * (
                p_minus[mu] * p_plus[nu]
                + p_minus[nu] * p_plus[mu]
                - _METRIC[mu, nu] * metric_term
            )
    return tensor


def _manual_contract(hadronic_cov: np.ndarray, tensor: np.ndarray) -> float:
    value = 0.0j
    for mu in range(4):
        for nu in range(4):
            value += hadronic_cov[mu] * np.conjugate(hadronic_cov[nu]) * tensor[mu, nu]
    return float(np.real(value))


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charm"
    assert constraint.observable == "BR(D+ -> pi+ e+- mu-+) LFV full-q2 proxy"


def test_anchor_matches_yaml_and_budget():
    constraint = fcc.get(_PID)
    data = _yaml()
    pdg = data["pdg_or_equivalent"]
    eplus = pdg["dplus_piplus_eplus_muminus"]
    eminus = pdg["dplus_piplus_eminus_muplus"]
    scope = pdg["lhcb_2021_search_scope"]
    babar = pdg["babar_2011_predecessor"]
    rs = data["paper_era_reference"]["rs_baseline"]

    assert constraint.anchor.eplus_muminus_limit.value == pytest.approx(eplus["value"])
    assert constraint.anchor.eplus_muminus_limit.confidence_level == pytest.approx(
        eplus["confidence_level"]
    )
    assert constraint.anchor.eplus_muminus_limit.source_url == eplus["source_url"]
    assert constraint.anchor.eplus_muminus_limit.companion_95cl_limit == pytest.approx(
        eplus["companion_95cl_limit"]
    )
    assert constraint.anchor.eminus_muplus_limit.value == pytest.approx(eminus["value"])
    assert constraint.anchor.eminus_muplus_limit.companion_95cl_limit == pytest.approx(
        eminus["companion_95cl_limit"]
    )
    assert constraint.anchor.search_scope.decay_modes_investigated == pytest.approx(
        scope["decay_modes_investigated"]
    )
    assert constraint.anchor.search_scope.integrated_luminosity_fb_inv == pytest.approx(
        scope["integrated_luminosity_fb_inv"]
    )
    assert constraint.anchor.babar_predecessor.eplus_muminus_limit == pytest.approx(
        babar["dplus_piplus_eplus_muminus_limit"]
    )
    assert constraint.anchor.babar_predecessor.eminus_muplus_limit == pytest.approx(
        babar["dplus_piplus_eminus_muplus_limit"]
    )
    assert constraint.anchor.rs_baseline.source == rs["source"]
    assert constraint.anchor.rs_baseline.use == rs["use"]
    assert constraint.anchor.budget == pytest.approx(min(eplus["value"], eminus["value"]))
    assert constraint.anchor.budget == pytest.approx(2.1e-7)


def test_anchor_loading_fails_loudly_for_missing_candidate_or_value():
    with pytest.raises(fcc.AnchorError, match="none of the expected anchor keys"):
        fcc.load_anchor(_PID, family="charm", candidates=("missing_c008_anchor",))
    with pytest.raises(fcc.AnchorError, match="has no 'missing_limit' field"):
        fcc.load_anchor(
            _PID,
            family="charm",
            candidates=("dplus_piplus_eplus_muminus",),
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
    assert result.diagnostics["sm_lfv_policy"].startswith("D+ -> pi+ e mu")
    assert result.diagnostics["prediction_per_charge_mode_branching_fraction"] == (
        pytest.approx(0.0)
    )
    assert result.passes is True


def test_lfv_dtopi_proxy_matches_independent_manual_integration():
    constraint = fcc.get(_PID)
    quark = _uc_couplings(left=2.0e-2 + 0.3e-2j, right=0.4e-2j)
    lepton = _lfv_proxy(left_emu_overlap=0.40 + 0.10j, right_emu_overlap=0.05j)
    result = constraint.evaluate(_point(quark, lepton))
    manual = _manual_dtopi_emu_proxy(quark, lepton, constraint.sd_inputs)

    assert result.predicted == pytest.approx(manual, rel=1e-12, abs=0.0)
    assert result.predicted == pytest.approx(
        1.2409725765054486e-11,
        rel=1e-12,
        abs=0.0,
    )
    assert result.diagnostics["eplus_muminus_ratio"] == pytest.approx(
        result.predicted / constraint.anchor.eplus_muminus_limit.value
    )
    assert result.diagnostics["eminus_muplus_ratio"] == pytest.approx(
        result.predicted / constraint.anchor.eminus_muplus_limit.value
    )
    assert result.ratio == pytest.approx(
        max(
            result.diagnostics["eplus_muminus_ratio"],
            result.diagnostics["eminus_muplus_ratio"],
        )
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["rs_matching_assumption"]


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    quark = _uc_couplings(left=1.0e-2 + 0.2e-2j, right=0.5e-2j)
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
        "c9_lfv_semileptonic_np",
        "c10_lfv_semileptonic_np",
        "c9_lfv_np",
        "c10_lfv_np",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "q2_min_gev2",
        "q2_max_gev2",
        "prediction_per_charge_mode_branching_fraction",
        "np_shift_branching_fraction",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["wilson_coefficients"]
    assert result.diagnostics["semileptonic_primed_combination_is_plus"] is True
    assert result.diagnostics["scalar_tensor_operators_included"] is False
    assert result.diagnostics["lhcb_window_acceptance_applied"] is False
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("quark", "lepton", "expected_pass"),
    [
        (_uc_couplings(left=1.0), _lfv_proxy(0.5), True),
        (_uc_couplings(left=5.0), _lfv_proxy(1.0), False),
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
        assert result.ratio > 1.0


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
    assert abs(heavy_result.diagnostics["c10_lfv_semileptonic_np"]) == pytest.approx(
        abs(default_result.diagnostics["c10_lfv_semileptonic_np"]) / 4.0
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
