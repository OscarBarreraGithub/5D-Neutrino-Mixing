"""Production tests for K020 (K+ -> pi+ mu+ e- LFV)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import (
    ConstraintLevel,
    ConstraintProtocol,
    Severity,
)
from flavor_catalog_constraints.secondary.kaon import K020 as k020_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_kaon_dilepton import (
    compute_rare_kaon_dilepton_wilsons,
    g_sm_squared,
)

_PID = "K020"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = (
    _REPO_ROOT
    / "flavor_catalog"
    / "processes"
    / "secondary"
    / "kaon"
    / "K020.yaml"
)
_METRIC = np.diag([1.0, -1.0, -1.0, -1.0])


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _value_entry(value_id: str):
    for entry in _yaml()["pdg_or_equivalent"]["values"]:
        if entry["value_id"] == value_id:
            return entry
    raise AssertionError(f"no K020 value entry for {value_id!r}")


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


def _lepton_proxy(
    left: complex,
    right: complex = 0.0j,
    M_KK: float = 3000.0,
) -> dict[str, complex | float | str]:
    return {
        "left_emu_overlap": left,
        "right_emu_overlap": right,
        "m_kk_gev": M_KK,
        "source": "test e-mu LFV proxy",
    }


def _point(
    quark_couplings: QuarkMassBasisCouplings,
    lepton_proxy: dict[str, complex | float | str],
    *,
    kk_ew_mass_gev: float | None = None,
):
    extras = {
        "quark_mass_basis_couplings": quark_couplings,
        "lepton_mass_basis_couplings": lepton_proxy,
    }
    if kk_ew_mass_gev is not None:
        extras["kk_ew_mass_gev"] = kk_ew_mass_gev
    return point_builder.make_point(**extras)


def _weak_z_coupling(inputs) -> float:
    g_weak = math.sqrt(4.0 * math.pi * inputs.alpha_em_mz / inputs.sin2_theta_w)
    cos_theta_w = math.sqrt(1.0 - inputs.sin2_theta_w)
    return float(g_weak / cos_theta_w)


def _independent_kplus_piplus_emu_branching(
    quark_couplings: QuarkMassBasisCouplings,
    lepton_proxy: dict[str, complex | float | str],
    inputs,
    *,
    m_kk_gev: float | None = None,
    q2_points: int = 160,
    angular_points: int = 48,
) -> float:
    rare_inputs = inputs.rare_kaon
    base = compute_rare_kaon_dilepton_wilsons(
        quark_couplings,
        m_kk_gev=m_kk_gev,
        inputs=rare_inputs,
    )
    g_z = _weak_z_coupling(rare_inputs)
    lepton_left = g_z * complex(lepton_proxy["left_emu_overlap"])
    lepton_right = g_z * complex(lepton_proxy["right_emu_overlap"])
    lepton_vector = lepton_left + lepton_right
    lepton_axial = lepton_right - lepton_left
    quark_vector = base.left_quark_delta + base.right_quark_delta
    normalizer = g_sm_squared(rare_inputs) * base.M_KK**2
    y_vector = quark_vector * lepton_vector / normalizer
    y_axial = quark_vector * lepton_axial / normalizer

    q2_min = (inputs.electron_mass_gev + inputs.muon_mass_gev) ** 2
    q2_max = (inputs.kplus_mass_gev - inputs.piplus_mass_gev) ** 2
    nodes, weights = np.polynomial.legendre.leggauss(q2_points)
    mid = 0.5 * (q2_max + q2_min)
    half_width = 0.5 * (q2_max - q2_min)
    total = 0.0
    for node, weight in zip(nodes, weights):
        q2 = mid + half_width * float(node)
        total += float(weight) * _independent_differential_branching(
            q2,
            y_vector=y_vector,
            y_axial=y_axial,
            inputs=inputs,
            angular_points=angular_points,
        )
    return float(half_width * total)


def _independent_differential_branching(
    q2: float,
    *,
    y_vector: complex,
    y_axial: complex,
    inputs,
    angular_points: int,
) -> float:
    m_minus = inputs.electron_mass_gev
    m_plus = inputs.muon_mass_gev
    m_k = inputs.kplus_mass_gev
    m_pi = inputs.piplus_mass_gev
    lam_had = max(0.0, _kallen(m_k * m_k, m_pi * m_pi, q2))
    lam_lep = max(0.0, _kallen(q2, m_minus * m_minus, m_plus * m_plus))
    if lam_had <= 0.0 or lam_lep <= 0.0:
        return 0.0

    angular = _independent_angular_tensor_integral(
        q2,
        m_minus=m_minus,
        m_plus=m_plus,
        y_vector=y_vector,
        y_axial=y_axial,
        inputs=inputs,
        angular_points=angular_points,
    )
    rare_inputs = inputs.rare_kaon
    c0 = (
        rare_inputs.gf_gev_minus2
        / math.sqrt(2.0)
        * rare_inputs.alpha_em_mz
        / (2.0 * math.pi * rare_inputs.sin2_theta_w)
    )
    tau = inputs.kplus_lifetime_s / inputs.hbar_gev_s
    phase_space = math.sqrt(lam_had) * math.sqrt(lam_lep) / (
        512.0 * math.pi**3 * m_k**3 * q2
    )
    return float(tau * c0 * c0 * phase_space * angular)


def _independent_angular_tensor_integral(
    q2: float,
    *,
    m_minus: float,
    m_plus: float,
    y_vector: complex,
    y_axial: complex,
    inputs,
    angular_points: int,
) -> float:
    nodes, weights = np.polynomial.legendre.leggauss(angular_points)
    total = 0.0
    for cos_theta, weight in zip(nodes, weights):
        vector, axial = _independent_tensor_contractions(
            q2,
            float(cos_theta),
            m_minus=m_minus,
            m_plus=m_plus,
            inputs=inputs,
        )
        total += float(weight) * (
            abs(y_vector) ** 2 * vector + abs(y_axial) ** 2 * axial
        )
    return float(total)


def _independent_tensor_contractions(
    q2: float,
    cos_theta: float,
    *,
    m_minus: float,
    m_plus: float,
    inputs,
) -> tuple[float, float]:
    root_q2 = math.sqrt(q2)
    m_k = inputs.kplus_mass_gev
    m_pi = inputs.piplus_mass_gev
    m_k2 = m_k * m_k
    m_pi2 = m_pi * m_pi
    delta = m_k2 - m_pi2
    lam_had = max(0.0, _kallen(m_k2, m_pi2, q2))
    lam_lep = max(0.0, _kallen(q2, m_minus * m_minus, m_plus * m_plus))
    p_had = math.sqrt(lam_had) / (2.0 * root_q2)
    p_lep = math.sqrt(lam_lep) / (2.0 * root_q2)
    e_k = (m_k2 - m_pi2 + q2) / (2.0 * root_q2)
    e_pi = (m_k2 - m_pi2 - q2) / (2.0 * root_q2)
    e_minus = (q2 + m_minus * m_minus - m_plus * m_plus) / (2.0 * root_q2)
    e_plus = (q2 + m_plus * m_plus - m_minus * m_minus) / (2.0 * root_q2)
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))

    p_k = np.array([e_k, 0.0, 0.0, p_had], dtype=np.complex128)
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
    ff = inputs.form_factor
    fplus = ff.fplus_0 / (1.0 - q2 / ff.vector_pole_mass_gev**2)
    fzero = ff.fzero_0 / (1.0 - q2 / ff.scalar_pole_mass_gev**2)
    hadronic = fplus * (p_k + p_pi - delta / q2 * q_vec) + fzero * delta / q2 * q_vec
    hadronic_cov = _METRIC @ hadronic
    p_dot = float(np.real(p_minus @ _METRIC @ p_plus))
    vector_tensor = _lepton_tensor(p_minus, p_plus, p_dot + m_minus * m_plus)
    axial_tensor = _lepton_tensor(p_minus, p_plus, p_dot - m_minus * m_plus)
    return (
        _contract_hadronic_tensor(hadronic_cov, vector_tensor),
        _contract_hadronic_tensor(hadronic_cov, axial_tensor),
    )


def _lepton_tensor(p_minus: np.ndarray, p_plus: np.ndarray, metric_term: float):
    tensor = np.zeros((4, 4), dtype=np.complex128)
    for mu in range(4):
        for nu in range(4):
            tensor[mu, nu] = 4.0 * (
                p_minus[mu] * p_plus[nu]
                + p_minus[nu] * p_plus[mu]
                - _METRIC[mu, nu] * metric_term
            )
    return tensor


def _contract_hadronic_tensor(hadronic_cov: np.ndarray, tensor: np.ndarray) -> float:
    value = 0.0j
    for mu in range(4):
        for nu in range(4):
            value += hadronic_cov[mu] * np.conjugate(hadronic_cov[nu]) * tensor[mu, nu]
    return float(np.real(value))


def _kallen(a: float, b: float, c: float) -> float:
    return float(a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c))


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.level is ConstraintLevel.SECONDARY
    assert constraint.family == "kaon"
    assert constraint.observable == "BR(K+ -> pi+ mu+ e-) LFV"


def test_anchor_matches_yaml_and_routes_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = k020_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(k020_module, "load_anchor", spy_load_anchor)
    anchor = k020_module._load_k020_anchor(_PID)
    current = _value_entry(k020_module._CURRENT_LIMIT_VALUE_ID)
    sher = _value_entry(k020_module._SHER_E865_VALUE_ID)
    opposite = _value_entry(k020_module._OPPOSITE_CHARGE_VALUE_ID)

    assert calls == [
        ("pdg_or_equivalent.values[0]",),
        ("pdg_or_equivalent.values[1]",),
        ("pdg_or_equivalent.values[2]",),
    ]
    assert anchor.current_limit.value == pytest.approx(float(current["value"]))
    assert anchor.current_limit.observable == current["observable"]
    assert anchor.current_limit.limit_type == current["limit_type"]
    assert anchor.current_limit.confidence_level == current["cl"]
    assert anchor.current_limit.units == current["units"]
    assert anchor.current_limit.source_url == current["source_url"]
    assert anchor.sher_e865_limit.value == pytest.approx(float(sher["value"]))
    assert anchor.sher_e865_limit.source_url == sher["source_url"]
    assert anchor.opposite_charge_limit.value == pytest.approx(float(opposite["value"]))
    assert anchor.opposite_charge_limit.source_url == opposite["source_url"]
    assert anchor.value == pytest.approx(1.3e-11)
    assert anchor.budget == pytest.approx(anchor.value)


def test_k020_anchor_loud_fail_probe(monkeypatch):
    broken = _yaml()
    first = dict(broken["pdg_or_equivalent"]["values"][0])
    first["cl"] = "95% CL"
    broken["pdg_or_equivalent"]["values"][0] = first

    def fake_load_full_yaml(*_args, **_kwargs):
        return broken

    monkeypatch.setattr(k020_module, "load_full_yaml", fake_load_full_yaml)
    with pytest.raises(k020_module.AnchorError):
        k020_module._load_k020_anchor(_PID)


@pytest.mark.parametrize(
    ("point", "missing"),
    [
        (
            point_builder.empty_point(),
            ("quark_mass_basis_couplings", "lepton_mass_basis_couplings"),
        ),
        (
            point_builder.build_from_quark_couplings(_sd_couplings(left=1.0e-5)),
            ("lepton_mass_basis_couplings",),
        ),
    ],
)
def test_evaluate_without_required_inputs_degrades_gracefully(point, missing):
    result = fcc.get(_PID).evaluate(point)

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extras"] == missing
    assert result.diagnostics["evaluated"] is False
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_lfv_branching_fraction_matches_independent_recomputation():
    constraint = fcc.get(_PID)
    quark_couplings = _sd_couplings(left=1.0e-3 + 0.2e-3j, right=0.3e-3j)
    lepton_proxy = _lepton_proxy(left=1.0 + 0.2j, right=0.1j)
    result = constraint.evaluate(_point(quark_couplings, lepton_proxy))
    independent = _independent_kplus_piplus_emu_branching(
        quark_couplings,
        lepton_proxy,
        constraint.sm_inputs,
    )

    assert result.predicted == pytest.approx(independent, rel=2.0e-10, abs=1.0e-24)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.ratio == pytest.approx(result.predicted / constraint.anchor.budget)
    assert result.diagnostics["sm_branching_fraction"] == pytest.approx(0.0)
    assert result.diagnostics["charge_mode"] == "muplus_eminus"
    assert result.diagnostics["charge_conjugate_modes_included"] is False
    assert result.diagnostics["quark_vector_current_uses_left_plus_right"] is True
    assert result.diagnostics["evaluated"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_evaluate_runs_end_to_end_with_real_fields_and_complex_diagnostics():
    quark_couplings = _sd_couplings(left=1.0e-3 + 0.2e-3j, right=0.3e-3j)
    lepton_proxy = _lepton_proxy(left=1.0 + 0.2j, right=0.1j)
    result = fcc.get(_PID).evaluate(_point(quark_couplings, lepton_proxy))

    for value in (
        result.predicted,
        result.sm_prediction,
        result.experimental,
        result.ratio,
        result.budget,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    for key in (
        "left_sd_coupling",
        "right_sd_coupling",
        "quark_vector_delta",
        "lepton_vector_delta_emu",
        "lepton_axial_delta_emu",
        "y_vector_lfv",
        "y_axial_lfv",
        "lambda_c",
        "lambda_t",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "q2_min_gev2",
        "q2_max_gev2",
        "effective_hamiltonian_prefactor_gev_minus2",
        "fplus_0",
        "fplus_q2_min",
        "fplus_q2_mid",
        "fplus_q2_max",
        "np_shift_branching_fraction",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert isinstance(result.diagnostics["wilson_coefficients"], dict)
    assert isinstance(
        result.diagnostics["wilson_coefficients"]["YV_LFV_semileptonic"],
        complex,
    )


@pytest.mark.parametrize(
    ("lepton_proxy", "expected_pass"),
    [
        (_lepton_proxy(left=1.0), True),
        (_lepton_proxy(left=5.0), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(lepton_proxy, expected_pass):
    constraint = fcc.get(_PID)
    quark_couplings = _sd_couplings(left=1.0e-3)
    result = constraint.evaluate(_point(quark_couplings, lepton_proxy))

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    constraint = fcc.get(_PID)
    quark_couplings = _sd_couplings(left=1.0e-3)
    lepton_proxy = _lepton_proxy(left=1.0)
    default_result = constraint.evaluate(_point(quark_couplings, lepton_proxy))
    ew_result = constraint.evaluate(
        _point(quark_couplings, lepton_proxy, kk_ew_mass_gev=6000.0)
    )

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert ew_result.predicted == pytest.approx(default_result.predicted / 16.0)


def test_evaluate_is_pure_and_deterministic():
    quark_couplings = _sd_couplings(left=1.0e-3 + 0.2e-3j, right=0.3e-3j)
    lepton_proxy = _lepton_proxy(left=1.0 + 0.2j, right=0.1j)
    before_left_down = quark_couplings.left_down.copy()
    before_right_down = quark_couplings.right_down.copy()
    before_lepton_proxy = dict(lepton_proxy)
    point = _point(quark_couplings, lepton_proxy)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(quark_couplings.left_down, before_left_down)
    np.testing.assert_array_equal(quark_couplings.right_down, before_right_down)
    assert lepton_proxy == before_lepton_proxy
