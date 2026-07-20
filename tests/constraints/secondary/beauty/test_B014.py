"""Production tests for B014 (exclusive b -> d gamma modes)."""

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintLevel, ConstraintProtocol, Severity
from flavor_catalog_constraints.secondary.beauty import B014 as b014_module
from qcd.constants import M_TOP_MS
from quarkConstraints import bsgamma as bsgamma_core
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.model import RotationParameters, ckm_like_unitary
from quarkConstraints.modern.inputs import ModernDefaultCKMTarget
from quarkConstraints.qcd_running import run_alpha_s

_PID = "B014"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = (
    _REPO_ROOT
    / "flavor_catalog"
    / "processes"
    / "secondary"
    / "beauty"
    / "B014.yaml"
)


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _bd_couplings(
    left: complex,
    right: complex,
    *,
    bs_left: complex = 0.0j,
    bs_right: complex = 0.0j,
    M_KK: float = 3000.0,
    g_s: float = 1.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings with the d-b slot populated."""

    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_down[0, 2] = left
    left_down[2, 0] = np.conj(left)
    right_down[0, 2] = right
    right_down[2, 0] = np.conj(right)
    left_down[1, 2] = bs_left
    left_down[2, 1] = np.conj(bs_left)
    right_down[1, 2] = bs_right
    right_down[2, 1] = np.conj(bs_right)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=g_s,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=zeros,
        left_down=left_down,
        right_up=zeros,
        right_down=right_down,
    )


def _bs_core_couplings_from_bd(
    couplings: QuarkMassBasisCouplings,
) -> QuarkMassBasisCouplings:
    """Manual d-b to s-b reshuffle for direct core recomputation."""

    left = complex(couplings.left_down[0, 2])
    right = complex(couplings.right_down[0, 2])
    return _bd_couplings(
        left=0.0j,
        right=0.0j,
        bs_left=left,
        bs_right=right,
        M_KK=float(couplings.M_KK),
        g_s=float(couplings.g_s),
    )


def _expected_ckm_power_suppression() -> float:
    target = ModernDefaultCKMTarget()
    matrix = ckm_like_unitary(
        RotationParameters(
            theta12=target.theta12,
            theta13=target.theta13,
            theta23=target.theta23,
            delta=target.delta,
        )
    )
    lambda_t_bd = complex(np.conjugate(matrix[2, 0]) * matrix[2, 2])
    lambda_t_bs = complex(np.conjugate(matrix[2, 1]) * matrix[2, 2])
    return float(abs(lambda_t_bd / lambda_t_bs) ** 2)


def _manual_c7_prediction(
    *,
    exclusive_normalization: float,
    c7_sm: complex,
    c7p_sm: complex,
    reference_scale_gev: float,
    proxy_normalization: float,
    c8_proxy_normalization: float,
    low_scale_gev: float,
    left_bd: complex,
    right_bd: complex,
    m_kk_gev: float = 3000.0,
    g_s: float = 1.0,
) -> dict[str, complex | float]:
    proxy_scale = proxy_normalization * (reference_scale_gev / m_kk_gev) ** 2
    c7_matching = proxy_scale * left_bd / g_s
    c7p_matching = proxy_scale * right_bd / g_s
    c8_matching = c8_proxy_normalization * proxy_scale * left_bd / g_s
    c8p_matching = c8_proxy_normalization * proxy_scale * right_bd / g_s
    u77, u78, u88 = _manual_ll_running_coefficients(
        matching_scale_gev=m_kk_gev,
        low_scale_gev=low_scale_gev,
    )
    c7_np = u77 * c7_matching + u78 * c8_matching
    c7p_np = u77 * c7p_matching + u78 * c8p_matching
    c8_np = u88 * c8_matching
    c8p_np = u88 * c8p_matching
    sm_power = abs(c7_sm) ** 2 + abs(c7p_sm) ** 2
    total_power = abs(c7_sm + c7_np) ** 2 + abs(c7p_sm + c7p_np) ** 2
    ratio_to_sm = total_power / sm_power
    branching_fraction = exclusive_normalization * ratio_to_sm
    return {
        "branching_fraction": float(branching_fraction),
        "ratio_to_sm": float(ratio_to_sm),
        "c7_np": complex(c7_np),
        "c7p_np": complex(c7p_np),
        "c8_np": complex(c8_np),
        "c8p_np": complex(c8p_np),
        "c7_np_matching": complex(c7_matching),
        "c7p_np_matching": complex(c7p_matching),
        "c7_running_from_c7": float(u77),
        "c7_running_from_c8": float(u78),
        "c8_running_from_c8": float(u88),
    }


def _manual_ll_running_coefficients(
    *,
    matching_scale_gev: float,
    low_scale_gev: float,
) -> tuple[float, float, float]:
    boundaries = [
        matching_scale_gev,
        *[
            scale
            for scale in (M_TOP_MS, 4.18, 1.27)
            if low_scale_gev < scale < matching_scale_gev
        ],
        low_scale_gev,
    ]
    u77_total = 1.0
    u78_total = 0.0
    u88_total = 1.0
    for high, low in zip(boundaries, boundaries[1:]):
        nf = 6 if high > M_TOP_MS else 5 if high > 4.18 else 4 if high > 1.27 else 3
        beta0 = (33.0 - 2.0 * nf) / 3.0
        eta = run_alpha_s(high) / run_alpha_s(low)
        u77 = eta ** ((32.0 / 3.0) / (2.0 * beta0))
        u88 = eta ** ((28.0 / 3.0) / (2.0 * beta0))
        u78 = (8.0 / 3.0) * (u88 - u77)
        u78_total = u77 * u78_total + u78 * u88_total
        u77_total *= u77
        u88_total *= u88
    return float(u77_total), float(u78_total), float(u88_total)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.level is ConstraintLevel.SECONDARY
    assert constraint.family == "beauty"
    assert constraint.observable == "BR(B+ -> rho+ gamma), BR(B0 -> rho0/omega gamma)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    charged = pdg["pdg2025_bp_to_rhop_gamma"]
    neutral = pdg["pdg2025_b0_to_rho0_gamma"]
    omega = pdg["pdg2025_b0_to_omega_gamma"]

    assert not any("sm" in key.lower() and "prediction" in key.lower() for key in pdg)
    assert constraint.anchor.charged_rho.value == pytest.approx(charged["value"])
    assert constraint.anchor.charged_rho.uncertainty == pytest.approx(
        charged["uncertainty"]
    )
    assert constraint.anchor.neutral_rho.value == pytest.approx(neutral["value"])
    assert constraint.anchor.neutral_rho.uncertainty == pytest.approx(
        neutral["uncertainty"]
    )
    assert constraint.anchor.omega.value == pytest.approx(omega["value"])
    assert constraint.anchor.omega.uncertainty == pytest.approx(
        max(omega["uncertainty_plus"], omega["uncertainty_minus"])
    )
    assert constraint.anchor.hflav_rho.value == pytest.approx(
        pdg["hflav2024_b_to_rho_gamma"]["value"]
    )
    assert constraint.anchor.hflav_rho_omega.value == pytest.approx(
        pdg["hflav2024_b_to_rho_omega_gamma"]["value"]
    )
    assert constraint.anchor.isospin_rho.value == pytest.approx(
        pdg["hflav2024_isospin_rho_gamma"]["value"]
    )
    assert constraint.anchor.xd_gamma.value == pytest.approx(
        pdg["hflav2024_b_to_xd_gamma"]["value"]
    )
    assert constraint.anchor.lhcb_ratio_rho0_kstar0.value == pytest.approx(
        pdg["lhcb2025_ratio_b0rho0_to_b0kstar0_gamma"]["value"]
    )
    assert constraint.anchor.budget == pytest.approx(neutral["uncertainty"])
    assert constraint.anchor.neutral_rho.budget_band.lower_total_edge == pytest.approx(
        neutral["value"] - neutral["uncertainty"]
    )
    assert constraint.anchor.omega.budget_band.upper_total_edge == pytest.approx(
        omega["value"] + max(omega["uncertainty_plus"], omega["uncertainty_minus"])
    )

    with pytest.raises(anchors.AnchorError):
        anchors.load_anchor(
            _PID,
            family="beauty",
            tier=ConstraintLevel.SECONDARY,
            candidates=("no_such_block",),
        )


def test_anchor_loader_rejects_mismatched_block_key(monkeypatch):
    original_load_anchor = b014_module.load_anchor

    def mismatched_load_anchor(*args, **kwargs):
        anchor = original_load_anchor(*args, **kwargs)
        if kwargs["candidates"] == ("pdg2025_b0_to_rho0_gamma",):
            return replace(anchor, block_key="wrong_block")
        return anchor

    monkeypatch.setattr(b014_module, "load_anchor", mismatched_load_anchor)

    with pytest.raises(anchors.AnchorError, match="load_anchor selected"):
        b014_module._load_b014_anchor(_PID)


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.sm_prediction == pytest.approx(constraint.anchor.reference_value)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert result.diagnostics["budget_sm_theory_sigma"] is None
    assert result.diagnostics["sm_theory_prediction_available"] is False
    assert "measurement-consistency band" in result.diagnostics["budget_policy"]
    assert "theory-only" in result.diagnostics["sm_theory_prediction_gap"]
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics[
        "ckm_power_suppression_vtd_over_vts_squared"
    ] == pytest.approx(_expected_ckm_power_suppression())


def test_sm_limit_branching_fraction_validates_exclusive_scale_and_ckm():
    constraint = fcc.get(_PID)
    couplings = _bd_couplings(left=0.0j, right=0.0j)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))
    manual = _manual_c7_prediction(
        exclusive_normalization=constraint.anchor.neutral_rho.reference_value,
        c7_sm=constraint.sm_inputs.c7_sm_eff,
        c7p_sm=constraint.sm_inputs.c7p_sm_eff,
        reference_scale_gev=constraint.sm_inputs.reference_scale_gev,
        proxy_normalization=constraint.sm_inputs.c7_proxy_normalization,
        c8_proxy_normalization=constraint.sm_inputs.c8_proxy_normalization,
        low_scale_gev=constraint.sm_inputs.low_scale_gev,
        left_bd=0.0j,
        right_bd=0.0j,
    )

    assert result.predicted == pytest.approx(manual["branching_fraction"])
    assert result.predicted == pytest.approx(constraint.anchor.neutral_rho.value)
    assert result.predicted == pytest.approx(8.6e-7)
    assert result.sm_prediction == pytest.approx(constraint.anchor.neutral_rho.value)
    assert result.diagnostics["dominant_mode"] == "b0_to_rho0_gamma"
    assert result.diagnostics[
        "normalization_formula_minus_reference_measurement"
    ] == pytest.approx(0.0)
    assert result.diagnostics["ratio_to_sm_c7_power"] == pytest.approx(1.0)
    assert result.diagnostics[
        "ckm_power_suppression_vtd_over_vts_squared"
    ] == pytest.approx(_expected_ckm_power_suppression())
    assert 0.035 < result.diagnostics["ckm_power_suppression_vtd_over_vts_squared"] < 0.05
    assert result.ratio == pytest.approx(0.0)
    assert result.passes is True


def test_bd_c7_proxy_numerics_match_independent_recomputation():
    constraint = fcc.get(_PID)
    left = 0.03e-1 + 0.02e-1j
    right = 1.0e-2 - 0.3e-2j
    couplings = _bd_couplings(left=left, right=right)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))
    expected = _manual_c7_prediction(
        exclusive_normalization=constraint.anchor.neutral_rho.reference_value,
        c7_sm=constraint.sm_inputs.c7_sm_eff,
        c7p_sm=constraint.sm_inputs.c7p_sm_eff,
        reference_scale_gev=constraint.sm_inputs.reference_scale_gev,
        proxy_normalization=constraint.sm_inputs.c7_proxy_normalization,
        c8_proxy_normalization=constraint.sm_inputs.c8_proxy_normalization,
        low_scale_gev=constraint.sm_inputs.low_scale_gev,
        left_bd=left,
        right_bd=right,
    )

    assert result.predicted == pytest.approx(expected["branching_fraction"])
    assert result.diagnostics["ratio_to_sm_c7_power"] == pytest.approx(
        expected["ratio_to_sm"]
    )
    assert result.diagnostics["c7_np"] == pytest.approx(expected["c7_np"])
    assert result.diagnostics["c7p_np"] == pytest.approx(expected["c7p_np"])
    assert result.diagnostics["c8_np"] == pytest.approx(expected["c8_np"])
    assert result.diagnostics["c8p_np"] == pytest.approx(expected["c8p_np"])
    assert result.diagnostics["c7_np_matching"] == pytest.approx(
        expected["c7_np_matching"]
    )
    assert result.diagnostics["c7p_np_matching"] == pytest.approx(
        expected["c7p_np_matching"]
    )
    assert result.diagnostics["left_bd_coupling"] == pytest.approx(left)
    assert result.diagnostics["right_bd_coupling"] == pytest.approx(right)
    assert result.ratio == pytest.approx(
        abs(result.predicted - constraint.anchor.neutral_rho.reference_value)
        / constraint.anchor.neutral_rho.budget
    )
    assert "BR(B0 -> rho0 gamma)" in result.diagnostics["exclusive_branching_formula"]
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_constraint_matches_shared_bsgamma_core_after_manual_bd_slot_mapping():
    constraint = fcc.get(_PID)
    left = 0.03e-1 + 0.02e-1j
    right = 1.0e-2 - 0.3e-2j
    couplings = _bd_couplings(left=left, right=right)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))

    core_couplings = _bs_core_couplings_from_bd(couplings)
    wilsons = bsgamma_core.compute_bsgamma_wilsons(
        core_couplings,
        inputs=constraint.sm_inputs,
    )
    direct = bsgamma_core.branching_fraction_from_c7(
        c7_np=wilsons.c7_np,
        c7p_np=wilsons.c7p_np,
        sm_branching_fraction=constraint.anchor.neutral_rho.reference_value,
        inputs=constraint.sm_inputs,
    )

    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.diagnostics["c7_np"] == pytest.approx(wilsons.c7_np)
    assert result.diagnostics["c7p_np"] == pytest.approx(wilsons.c7p_np)
    assert wilsons.c7_np_matching == pytest.approx(left)
    assert wilsons.c7p_np_matching == pytest.approx(right)


def test_bsgamma_bs_slot_is_ignored_for_b014():
    constraint = fcc.get(_PID)
    couplings = _bd_couplings(
        left=0.0j,
        right=0.0j,
        bs_left=1.0 + 0.2j,
        bs_right=0.8 - 0.1j,
    )
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))

    assert result.predicted == pytest.approx(constraint.anchor.neutral_rho.value)
    assert result.ratio == pytest.approx(0.0)
    assert result.diagnostics["left_bd_coupling"] == pytest.approx(0.0j)
    assert result.diagnostics["right_bd_coupling"] == pytest.approx(0.0j)
    assert result.diagnostics["down_sector_indices"] == (0, 2)


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    couplings = _bd_couplings(left=2.0e-3 + 0.5e-3j, right=1.0e-2 - 0.3e-2j)
    result = fcc.get(_PID).evaluate(point_builder.build_from_quark_couplings(couplings))

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
        "left_bd_coupling",
        "right_bd_coupling",
        "lambda_t_bd",
        "lambda_t_bs",
        "c7_sm_eff",
        "c7_total",
        "c7_np",
        "c7p_np",
        "c7_np_matching",
        "c7p_np_matching",
        "c8_np",
        "c8p_np",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "low_scale_gev",
        "proxy_scale_factor",
        "dipole_power_sm",
        "dipole_power_total",
        "np_shift_branching_fraction",
        "max_mode_ratio",
        "c7_running_from_c7",
        "c7_running_from_c8",
        "c8_running_from_c8",
        "alpha_s_matching_scale",
        "alpha_s_low_scale",
        "exclusive_normalization_branching_fraction",
        "abs_vtd_over_vts",
        "ckm_power_suppression_vtd_over_vts_squared",
        "b0_to_rho0_gamma_hard_veto_np_shift_budget",
        "bp_to_rhop_gamma_hard_veto_np_shift_budget",
        "b0_to_omega_gamma_hard_veto_np_shift_budget",
        "hflav_xd_gamma_branching_fraction",
        "lhcb_ratio_b0rho0_to_b0kstar0_gamma",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["wilson_coefficients"]["C7_NP"] == pytest.approx(
        result.diagnostics["c7_np"]
    )
    assert result.diagnostics["bsgamma_rg_running_applied"] is True
    assert result.diagnostics["kk_ew_mass_extra_used"] is False
    assert result.diagnostics["down_sector_indices"] == (0, 2)
    assert result.diagnostics["hflav_isospin_rho_gamma_evaluated"] is False
    assert result.diagnostics["budget_sm_theory_sigma"] is None
    assert result.diagnostics["sm_theory_prediction_available"] is False


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_bd_couplings(left=0.0j, right=2.0e-1), True),
        (_bd_couplings(left=0.0j, right=3.5e-1), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    result = fcc.get(_PID).evaluate(point_builder.build_from_quark_couplings(couplings))

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    couplings = _bd_couplings(left=0.0j, right=3.5e-1)
    default_point = point_builder.build_from_quark_couplings(couplings)
    heavy_point = point_builder.make_point(
        quark_mass_basis_couplings=couplings,
        kk_ew_mass_gev=6000.0,
    )
    default_result = fcc.get(_PID).evaluate(default_point)
    heavy_result = fcc.get(_PID).evaluate(heavy_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert heavy_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert heavy_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert abs(heavy_result.diagnostics["c7p_np_matching"]) == pytest.approx(
        abs(default_result.diagnostics["c7p_np_matching"]) / 4.0
    )
    assert heavy_result.diagnostics["c7_running_from_c7"] < (
        default_result.diagnostics["c7_running_from_c7"]
    )


def test_evaluate_is_pure_and_deterministic():
    couplings = _bd_couplings(left=2.0e-3 + 0.5e-3j, right=1.0e-2 - 0.3e-2j)
    before_left_down = couplings.left_down.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)
