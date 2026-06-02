"""Production tests for K019 (K_L -> e+- mu-+ LFV)."""

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
from flavor_catalog_constraints.secondary.kaon import K019 as k019_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.model import RotationParameters, ckm_like_unitary
from quarkConstraints.rare_kaon_dilepton import (
    compute_rare_kaon_dilepton_wilsons,
)

_PID = "K019"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = (
    _REPO_ROOT
    / "flavor_catalog"
    / "processes"
    / "secondary"
    / "kaon"
    / "K019.yaml"
)


def _yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


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


def _independent_lambda_wolfenstein(inputs) -> float:
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=inputs.theta12,
            theta13=inputs.theta13,
            theta23=inputs.theta23,
            delta=inputs.delta,
        )
    )
    return float(abs(matrix[0, 1]))


def _independent_klong_emu_branching(
    quark_couplings: QuarkMassBasisCouplings,
    lepton_proxy: dict[str, complex | float | str],
    inputs,
    *,
    m_kk_gev: float | None = None,
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
    y_base_total = base.y_np_left - base.y_np_right
    y_vector = y_base_total * lepton_vector / base.muon_axial_delta
    y_axial = y_base_total * lepton_axial / base.muon_axial_delta

    lam = _independent_lambda_wolfenstein(rare_inputs)
    kappa_mu = rare_inputs.kappa_mu_ref * (
        lam / rare_inputs.kappa_lambda_ref
    ) ** 8

    r_e = inputs.electron_mass_gev / inputs.kaon_mass_gev
    r_mu = inputs.muon_mass_gev / inputs.kaon_mass_gev
    x = r_e * r_e
    y = r_mu * r_mu
    phase_lambda = 1.0 + x * x + y * y - 2.0 * (x + y + x * y)
    beta_mumu = math.sqrt(1.0 - 4.0 * r_mu * r_mu)
    mass_sum = r_e + r_mu
    mass_diff = r_mu - r_e
    vector_term = (1.0 - mass_sum * mass_sum) * abs(mass_diff * y_vector) ** 2
    axial_term = (1.0 - mass_diff * mass_diff) * abs(mass_sum * y_axial) ** 2
    return float(
        inputs.charge_state_factor
        * kappa_mu
        * (math.sqrt(phase_lambda) / beta_mumu)
        * (vector_term + axial_term)
        / (4.0 * r_mu * r_mu * lam**10)
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.level is ConstraintLevel.SECONDARY
    assert constraint.family == "kaon"
    assert constraint.observable == "BR(K_L -> e+- mu-+) LFV"


def test_anchor_matches_yaml_and_routes_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = k019_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(k019_module, "load_anchor", spy_load_anchor)
    anchor = k019_module._load_k019_anchor(_PID)
    pdg = _yaml()["pdg_or_equivalent"]
    e871 = pdg["experiment_input"]

    assert calls == [("pdg_or_equivalent",)]
    assert anchor.current_limit.value == pytest.approx(float(pdg["value"]))
    assert anchor.current_limit.limit_type == pdg["limit_type"]
    assert anchor.current_limit.confidence_level == pdg["cl"]
    assert anchor.current_limit.units == pdg["units"]
    assert anchor.current_limit.source_url == pdg["source_url"]
    assert anchor.bnl_e871_limit.value == pytest.approx(float(e871["value"]))
    assert anchor.bnl_e871_limit.experiment == e871["experiment"]
    assert anchor.bnl_e871_limit.source_url == e871["source_url"]
    assert anchor.value == pytest.approx(4.7e-12)
    assert anchor.budget == pytest.approx(anchor.value)


def test_k019_anchor_loud_fail_probe(monkeypatch):
    broken = _yaml()
    broken["pdg_or_equivalent"] = dict(
        broken["pdg_or_equivalent"],
        cl="95% CL",
    )

    def fake_load_full_yaml(*_args, **_kwargs):
        return broken

    monkeypatch.setattr(k019_module, "load_full_yaml", fake_load_full_yaml)
    with pytest.raises(k019_module.AnchorError):
        k019_module._load_k019_anchor(_PID)


@pytest.mark.parametrize(
    ("point", "missing"),
    [
        (point_builder.empty_point(), ("quark_mass_basis_couplings", "lepton_mass_basis_couplings")),
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
    lepton_proxy = _lepton_proxy(left=0.1 + 0.02j, right=0.01j)
    result = constraint.evaluate(_point(quark_couplings, lepton_proxy))
    independent = _independent_klong_emu_branching(
        quark_couplings,
        lepton_proxy,
        constraint.sm_inputs,
    )

    assert result.predicted == pytest.approx(independent)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.ratio == pytest.approx(result.predicted / constraint.anchor.budget)
    assert result.diagnostics["sm_branching_fraction"] == pytest.approx(0.0)
    assert result.diagnostics["charge_conjugate_modes_included"] is True
    assert result.diagnostics["evaluated"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_evaluate_runs_end_to_end_with_real_fields_and_complex_diagnostics():
    quark_couplings = _sd_couplings(left=1.0e-3 + 0.2e-3j, right=0.3e-3j)
    lepton_proxy = _lepton_proxy(left=0.1 + 0.02j, right=0.01j)
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
        "y_vector_lfv",
        "y_axial_lfv",
        "lambda_c",
        "lambda_t",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "kappa_mu",
        "lambda_wolfenstein",
        "phase_space_lambda_sqrt",
        "np_shift_branching_fraction",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert isinstance(result.diagnostics["wilson_coefficients"], dict)
    assert isinstance(
        result.diagnostics["wilson_coefficients"]["YV_LFV_total"],
        complex,
    )


@pytest.mark.parametrize(
    ("lepton_proxy", "expected_pass"),
    [
        (_lepton_proxy(left=0.1), True),
        (_lepton_proxy(left=1.0), False),
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
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    constraint = fcc.get(_PID)
    quark_couplings = _sd_couplings(left=1.0e-3)
    lepton_proxy = _lepton_proxy(left=0.1)
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
    lepton_proxy = _lepton_proxy(left=0.1 + 0.02j, right=0.01j)
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
