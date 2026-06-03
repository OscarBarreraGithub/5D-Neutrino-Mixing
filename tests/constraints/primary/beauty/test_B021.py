"""Production tests for B021 (baryonic Lambda_b -> Lambda mu mu)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.beauty import B021 as b021_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.model import RotationParameters, ckm_like_unitary
from quarkConstraints.rare_b_baryon_dilepton import evaluate_lambdab_to_lambda_mumu
from tests.rare_b_phase3d_helpers import (
    core_wilsons_from_rs_coeff,
    rs_coeff,
    sample_rare_b_point,
    scaled_rare_b_point,
    sm_limit_rare_b_point,
)

_PID = "B021"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B021.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _find_entry(entries, key: str, expected: str):
    for entry in entries:
        if str(entry.get(key)) == expected:
            return entry
    raise AssertionError(f"missing {key}={expected!r}")


def _sb_couplings(
    left: complex,
    right: complex = 0.0j,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with only the s-b slot populated."""

    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[1, 2] = left
    left_down[2, 1] = np.conj(left)
    right_down[1, 2] = right
    right_down[2, 1] = np.conj(right)
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


def _manual_fplus_fminus(q2: float, mode) -> tuple[float, float]:
    e_lambda = (
        mode.parent_mass_gev**2 + mode.daughter_mass_gev**2 - q2
    ) / (2.0 * mode.parent_mass_gev)
    delta_e = e_lambda - mode.daughter_mass_gev
    f_plus = mode.n_plus_gev2 / (mode.x_plus_gev + delta_e) ** 2
    f_minus = mode.n_minus_gev2 / (mode.x_minus_gev + delta_e) ** 2
    return float(f_plus), float(f_minus)


def _manual_lambdab_to_lambda_mumu(
    inputs,
    *,
    q2_min_gev2: float,
    q2_max_gev2: float,
    c9_total: complex,
    c10_total: complex,
    grid_points: int = 12001,
) -> float:
    sd = inputs.short_distance_inputs
    mode = inputs.form_factor
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=sd.theta12,
            theta13=sd.theta13,
            theta23=sd.theta23,
            delta=sd.delta,
        )
    )
    lambda_t = complex(matrix[2, 2] * np.conjugate(matrix[2, 1]))
    tau = mode.lifetime_ps * 1.0e-12 / sd.hbar_gev_s
    threshold = 4.0 * sd.muon_mass_gev**2
    m_parent = mode.parent_mass_gev
    m_daughter = mode.daughter_mass_gev
    c_gamma = mode.hqet_c_gamma
    c_v = mode.hqet_c_v

    def dbr(q2: float) -> float:
        if q2 <= threshold or q2 >= mode.q2_max_gev2:
            return 0.0
        beta2 = 1.0 - threshold / q2
        phase_left = (m_parent - m_daughter) ** 2 - q2
        phase_right = (m_parent + m_daughter) ** 2 - q2
        if beta2 <= 0.0 or phase_left <= 0.0 or phase_right <= 0.0:
            return 0.0
        f_plus, f_minus = _manual_fplus_fminus(q2, mode)
        f_combined = phase_left * f_minus**2 + phase_right * f_plus**2
        g_combined = (
            m_parent**6
            - m_parent**4 * (3.0 * m_daughter**2 + q2)
            - m_parent**2
            * (q2 - m_daughter**2)
            * (3.0 * m_daughter**2 + q2)
            + (q2 - m_daughter**2) ** 3
        )
        a_10_10 = (
            (
                (2.0 * c_gamma**2 + 2.0 * c_gamma * c_v + c_v**2)
                * (2.0 * sd.muon_mass_gev**2 + q2)
                * (
                    m_parent**4
                    - 2.0 * m_parent**2 * m_daughter**2
                    + (q2 - m_daughter**2) ** 2
                )
                + 2.0
                * m_parent**2
                * q2
                * (
                    4.0 * c_gamma**2 * (q2 - 4.0 * sd.muon_mass_gev**2)
                    - (2.0 * c_gamma * c_v + c_v**2)
                    * (q2 - 10.0 * sd.muon_mass_gev**2)
                )
            )
            * f_combined
            + 4.0
            * c_gamma
            * (c_gamma + c_v)
            * (2.0 * sd.muon_mass_gev**2 + q2)
            * g_combined
            * f_plus
            * f_minus
        )
        a_9_9 = (
            (
                (2.0 * c_gamma**2 + 2.0 * c_gamma * c_v + c_v**2)
                * (m_parent**4 + (q2 - m_daughter**2) ** 2)
                - 2.0
                * m_parent**2
                * (
                    2.0 * c_gamma**2 * (m_daughter**2 - 2.0 * q2)
                    + (2.0 * c_gamma * c_v + c_v**2)
                    * (m_daughter**2 + q2)
                )
            )
            * f_combined
            + 4.0 * c_gamma * (c_gamma + c_v) * g_combined * f_plus * f_minus
        )
        prefactor = (
            sd.alpha_em_mz**2
            * sd.gf_gev_minus2**2
            * abs(lambda_t) ** 2
            / (6144.0 * math.pi**5 * q2**2 * m_parent**5)
        )
        phase = math.sqrt(beta2) * math.sqrt(phase_left * phase_right)
        return float(
            tau
            * prefactor
            * phase
            * (
                q2 * abs(complex(c10_total)) ** 2 * a_10_10
                + q2
                * (q2 + 2.0 * sd.muon_mass_gev**2)
                * abs(complex(c9_total)) ** 2
                * a_9_9
            )
        )

    xs = np.linspace(q2_min_gev2, q2_max_gev2, grid_points)
    ys = np.array([dbr(float(x)) for x in xs])
    return float(np.trapezoid(ys, xs))


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "beauty"
    assert constraint.observable == "BR(Lambda_b0 -> Lambda mu+ mu-) high-q2"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    total = _find_entry(
        pdg["observables"],
        "name",
        "BR(Lambda_b0 -> Lambda mu+ mu-)",
    )
    high = _find_entry(
        pdg["observables"],
        "name",
        "dBR/dq2(Lambda_b0 -> Lambda mu+ mu-)",
    )
    total_value = float(total["value"]) * 1.0e-6
    total_sigma = 0.28e-6
    high_avg = float(high["value"]) * 1.0e-7
    q2_width = 5.0
    high_value = high_avg * q2_width
    high_sigma_up = math.sqrt(0.09**2 + 0.03**2 + 0.27**2) * 1.0e-7 * q2_width
    high_sigma_down = math.sqrt(0.08**2 + 0.03**2 + 0.27**2) * 1.0e-7 * q2_width
    central = abs(high_value - constraint.sm_result.branching_fraction)
    proxy = (
        b021_module.RARE_B_BARYONIC_DILEPTON_PROXY_THEORY_UNCERTAINTY_FRACTION
        * constraint.sm_result.branching_fraction
    )

    assert constraint.anchor.total_branching_fraction.value == pytest.approx(
        total_value
    )
    assert constraint.anchor.total_branching_fraction.uncertainty == pytest.approx(
        total_sigma
    )
    assert constraint.anchor.high_q2_bin.value == pytest.approx(high_value)
    assert constraint.anchor.high_q2_bin.average_differential_value == pytest.approx(
        high_avg
    )
    assert constraint.anchor.high_q2_bin.uncertainty_upper == pytest.approx(
        high_sigma_up
    )
    assert constraint.anchor.high_q2_bin.uncertainty_lower == pytest.approx(
        high_sigma_down
    )
    assert constraint.anchor.high_q2_bin.source_url == high["source_url"]
    assert constraint.anchor.high_q2_bin.snapshot_path == high["snapshot_path"]
    assert constraint.anchor.high_q2_bin.year == high["year"]
    assert constraint.anchor.q2_min_gev2 == pytest.approx(15.0)
    assert constraint.anchor.q2_max_gev2 == pytest.approx(20.0)
    assert constraint.anchor.budget_band.central_residual == pytest.approx(central)
    assert constraint.anchor.budget_band.proxy_theory_sigma == pytest.approx(proxy)
    assert constraint.anchor.budget == pytest.approx(
        central + high_sigma_up + proxy
    )
    assert constraint.anchor.budget == pytest.approx(3.458207347033031e-07)
    assert "form-factor" in constraint.anchor.budget_band.proxy_theory_rationale
    assert "NEEDS-HUMAN-PHYSICS" in constraint.anchor.budget_band.construction

    with pytest.raises(AnchorError):
        b021_module._load_branching_observable(
            process_id=_PID,
            observable_name="no such B021 observable",
        )


def test_observable_anchors_route_through_scaffold_load_anchor_and_fail_loudly(
    monkeypatch,
):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b021_module.anchor_scaffold.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b021_module.anchor_scaffold, "load_anchor", spy_load_anchor)
    anchor = b021_module._load_b021_anchor(
        _PID,
        formula_sm=fcc.get(_PID).sm_result.branching_fraction,
    )

    assert calls == [
        ("pdg_or_equivalent.observables[0]",),
        ("pdg_or_equivalent.observables[1]",),
    ]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent.observables[99]",
            value=1.0,
            uncertainty=0.1,
        )

    monkeypatch.setattr(
        b021_module.anchor_scaffold,
        "load_anchor",
        mismatched_load_anchor,
    )
    with pytest.raises(AnchorError, match="load_anchor selected"):
        b021_module._load_b021_anchor(
            _PID,
            formula_sm=fcc.get(_PID).sm_result.branching_fraction,
        )


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())
    old_style_result = constraint.evaluate(
        point_builder.build_from_quark_couplings(_sb_couplings(left=1.0e-2))
    )

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.sm_prediction == pytest.approx(
        constraint.sm_result.branching_fraction
    )
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert old_style_result.passes is True
    assert old_style_result.predicted is None
    assert old_style_result.diagnostics["evaluated"] is False
    assert old_style_result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert old_style_result.diagnostics["legacy_quark_mass_basis_couplings_present"] is True


def test_sm_limit_branching_fraction_matches_independent_integral_and_catalog_bin():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(sm_limit_rare_b_point())
    manual = _manual_lambdab_to_lambda_mumu(
        constraint.sm_inputs,
        q2_min_gev2=constraint.anchor.q2_min_gev2,
        q2_max_gev2=constraint.anchor.q2_max_gev2,
        c9_total=constraint.sm_inputs.c9_sm,
        c10_total=constraint.sm_inputs.short_distance_inputs.c10_sm,
    )

    assert result.predicted == pytest.approx(manual, rel=3.0e-4)
    assert result.predicted == pytest.approx(5.532430650131983e-07)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.predicted / 5.0 == pytest.approx(1.1064861300263965e-07)
    assert abs(
        result.predicted / 5.0
        - constraint.anchor.high_q2_bin.average_differential_value
    ) < constraint.anchor.high_q2_bin.average_differential_uncertainty_upper
    assert result.diagnostics[
        "sm_average_differential_branching_fraction"
    ] == pytest.approx(result.sm_prediction / 5.0)
    assert result.diagnostics["form_factor_bundle"] == (
        "lambdab_lambda_fplus_fminus_detmold2013_hqet_dipole_v1"
    )
    assert "form-factor" in result.diagnostics["form_factor_uncertainty"]
    assert result.ratio == pytest.approx(
        abs(result.predicted - constraint.anchor.value) / constraint.anchor.budget
    )
    assert result.passes is True


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    point = sample_rare_b_point()
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
        "left_qb_coupling",
        "right_qb_coupling",
        "lambda_t",
        "c9_total",
        "c10_total",
        "c9_vector_np",
        "c10_axial_np",
        "c9_np",
        "c10_np",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "q2_min_gev2",
        "q2_max_gev2",
        "average_differential_branching_fraction",
        "sm_average_differential_branching_fraction",
        "fplus_q2_mid",
        "fminus_q2_mid",
        "budget_central_residual",
        "budget_experimental_sigma_used",
        "budget_proxy_theory_sigma",
        "budget_proxy_theory_fraction",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["pdg_total_scope"] == (
        "diagnostic_only_not_the_active_high_q2_bin"
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics[
        "budget_proxy_theory_rationale"
    ]
    assert result.diagnostics["c7_nonlocal_charm_omitted"] is True
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert result.diagnostics["rs_semileptonic_matching_status"].endswith(
        "no_second_1_over_M_KK_squared"
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("scale", "expected_pass"),
    [
        (1.0e-3, True),
        (1.0e5, False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    scale: float,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = scaled_rare_b_point(scale=scale)
    result = constraint.evaluate(point)
    direct = evaluate_lambdab_to_lambda_mumu(
        core_wilsons_from_rs_coeff(
            rs_coeff(point, transition="b_s", lepton="mu"),
            matching_scale_gev=point.extras["kk_ew_mass_gev"],
        ),
        q2_min_gev2=constraint.anchor.q2_min_gev2,
        q2_max_gev2=constraint.anchor.q2_max_gev2,
        inputs=constraint.sm_inputs,
    )

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.sm_prediction == pytest.approx(direct.sm_branching_fraction)
    assert result.ratio == pytest.approx(
        abs(direct.branching_fraction - constraint.anchor.value)
        / constraint.anchor.budget
    )
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    default_point = scaled_rare_b_point(scale=20.0)
    ew_point = scaled_rare_b_point(scale=20.0, kk_ew_mass_gev=6000.0)
    default_result = fcc.get(_PID).evaluate(default_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert ew_result.predicted == pytest.approx(default_result.predicted)
    assert ew_result.diagnostics["c9_vector_np"] == pytest.approx(
        default_result.diagnostics["c9_vector_np"]
    )
    assert ew_result.diagnostics["c10_axial_np"] == pytest.approx(
        default_result.diagnostics["c10_axial_np"]
    )
    assert ew_result.diagnostics["second_mkk_suppression_applied"] is False


def test_evaluate_is_pure_and_deterministic():
    point = sample_rare_b_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
