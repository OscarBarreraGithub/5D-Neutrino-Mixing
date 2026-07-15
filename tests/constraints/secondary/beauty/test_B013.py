"""Production tests for B013 (exclusive B_s0 -> phi(1020) gamma)."""

from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintLevel, ConstraintProtocol, Severity
from qcd.constants import M_TOP_MS
from flavor_catalog_constraints.secondary.beauty import B013 as b013_module
from quarkConstraints import bsgamma as bsgamma_core
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.qcd_running import run_alpha_s

_PID = "B013"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = (
    _REPO_ROOT
    / "flavor_catalog"
    / "processes"
    / "secondary"
    / "beauty"
    / "B013.yaml"
)


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _bs_couplings(
    left: complex,
    right: complex,
    M_KK: float = 3000.0,
    g_s: float = 1.0,
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
        g_s=g_s,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=zeros,
        left_down=left_down,
        right_up=zeros,
        right_down=right_down,
    )


def _manual_c7_prediction(
    *,
    exclusive_normalization: float,
    c7_sm: complex,
    c7p_sm: complex,
    reference_scale_gev: float,
    proxy_normalization: float,
    c8_proxy_normalization: float,
    low_scale_gev: float,
    left: complex,
    right: complex,
    m_kk_gev: float = 3000.0,
    g_s: float = 1.0,
) -> dict[str, complex | float]:
    proxy_scale = proxy_normalization * (reference_scale_gev / m_kk_gev) ** 2
    c7_matching = proxy_scale * left / g_s
    c7p_matching = proxy_scale * right / g_s
    c8_matching = c8_proxy_normalization * proxy_scale * left / g_s
    c8p_matching = c8_proxy_normalization * proxy_scale * right / g_s
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
        "c8_np_matching": complex(c8_matching),
        "c8p_np_matching": complex(c8p_matching),
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
    assert constraint.observable == "BR(B_s0 -> phi(1020) gamma)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["branching_fraction_hflav_2024"]
    ref = pdg["branching_fraction_pdg_2025"]
    s_phi = pdg["s_phigamma_lhcb_2019"]
    a_delta = pdg["a_delta_phigamma_lhcb_2019"]
    branching_rows = {key for key in pdg if "branching_fraction" in key}
    combined_sigma = math.sqrt(exp["uncertainty"] ** 2 + ref["uncertainty"] ** 2)
    central_residual = abs(exp["value"] - ref["value"])

    assert branching_rows == {
        "branching_fraction_hflav_2024",
        "branching_fraction_pdg_2025",
        "belle_2014_branching_fraction_input",
    }

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.uncertainty == pytest.approx(
        exp["uncertainty"]
    )
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.reference_measurement.value == pytest.approx(ref["value"])
    assert constraint.anchor.reference_measurement.uncertainty == pytest.approx(
        ref["uncertainty"]
    )
    assert constraint.anchor.reference_measurement.source_url == ref["source_url"]
    assert constraint.anchor.s_phigamma.value == pytest.approx(s_phi["value"])
    assert constraint.anchor.s_phigamma.uncertainty == pytest.approx(
        s_phi["uncertainty_stat"]
    )
    assert constraint.anchor.a_delta_phigamma.value == pytest.approx(
        a_delta["value"]
    )
    assert constraint.anchor.a_delta_phigamma.uncertainty == pytest.approx(
        a_delta["uncertainty_stat_plus"]
    )
    assert constraint.anchor.budget_band.central_residual == pytest.approx(
        central_residual
    )
    assert constraint.anchor.budget_band.combined_sigma == pytest.approx(
        combined_sigma
    )
    assert constraint.anchor.budget_band.reference_measurement_sigma == pytest.approx(
        ref["uncertainty"]
    )
    assert constraint.anchor.budget == pytest.approx(
        central_residual + combined_sigma
    )

    with pytest.raises(anchors.AnchorError):
        anchors.load_anchor(
            _PID,
            family="beauty",
            tier=ConstraintLevel.SECONDARY,
            candidates=("no_such_block",),
        )


def test_anchor_loader_rejects_mismatched_block_key(monkeypatch):
    original_load_anchor = b013_module.load_anchor

    def mismatched_load_anchor(*args, **kwargs):
        anchor = original_load_anchor(*args, **kwargs)
        if kwargs["candidates"] == ("branching_fraction_hflav_2024",):
            return replace(anchor, block_key="wrong_block")
        return anchor

    monkeypatch.setattr(b013_module, "load_anchor", mismatched_load_anchor)

    with pytest.raises(anchors.AnchorError, match="load_anchor selected"):
        b013_module._load_b013_anchor(_PID)


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
    assert result.diagnostics["helicity_observables_evaluated"] is False
    assert result.diagnostics["budget_sm_theory_sigma"] is None
    assert result.diagnostics["sm_theory_prediction_available"] is False
    assert "measurement-consistency band" in result.diagnostics["budget_policy"]
    assert "theory-only" in result.diagnostics["sm_theory_prediction_gap"]
    assert result.diagnostics["a_delta_phigamma_uncertainty_stat_plus"] == pytest.approx(
        0.37
    )
    assert result.diagnostics["a_delta_phigamma_uncertainty_stat_minus"] == pytest.approx(
        0.41
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_sm_limit_branching_fraction_validates_exclusive_scale():
    constraint = fcc.get(_PID)
    couplings = _bs_couplings(left=0.0j, right=0.0j)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))
    manual = _manual_c7_prediction(
        exclusive_normalization=constraint.anchor.reference_value,
        c7_sm=constraint.sm_inputs.c7_sm_eff,
        c7p_sm=constraint.sm_inputs.c7p_sm_eff,
        reference_scale_gev=constraint.sm_inputs.reference_scale_gev,
        proxy_normalization=constraint.sm_inputs.c7_proxy_normalization,
        c8_proxy_normalization=constraint.sm_inputs.c8_proxy_normalization,
        low_scale_gev=constraint.sm_inputs.low_scale_gev,
        left=0.0j,
        right=0.0j,
    )

    assert result.predicted == pytest.approx(manual["branching_fraction"])
    assert result.predicted == pytest.approx(constraint.anchor.reference_value)
    assert result.predicted == pytest.approx(3.4e-5)
    assert result.sm_prediction == pytest.approx(constraint.anchor.reference_value)
    assert result.diagnostics[
        "normalization_formula_minus_reference_measurement"
    ] == pytest.approx(0.0)
    assert result.diagnostics["ratio_to_sm_c7_power"] == pytest.approx(1.0)
    assert result.ratio == pytest.approx(0.0)
    assert result.passes is True


def test_c7_proxy_numerics_match_independent_recomputation():
    constraint = fcc.get(_PID)
    left = 0.03e-1 + 0.02e-1j
    right = 1.0e-2 - 0.3e-2j
    couplings = _bs_couplings(left=left, right=right)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))
    expected = _manual_c7_prediction(
        exclusive_normalization=constraint.anchor.reference_value,
        c7_sm=constraint.sm_inputs.c7_sm_eff,
        c7p_sm=constraint.sm_inputs.c7p_sm_eff,
        reference_scale_gev=constraint.sm_inputs.reference_scale_gev,
        proxy_normalization=constraint.sm_inputs.c7_proxy_normalization,
        c8_proxy_normalization=constraint.sm_inputs.c8_proxy_normalization,
        low_scale_gev=constraint.sm_inputs.low_scale_gev,
        left=left,
        right=right,
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
    assert result.diagnostics["c7_running_from_c7"] == pytest.approx(
        expected["c7_running_from_c7"]
    )
    assert result.diagnostics["c7_running_from_c8"] == pytest.approx(
        expected["c7_running_from_c8"]
    )
    assert result.ratio == pytest.approx(
        abs(result.predicted - constraint.anchor.reference_value)
        / constraint.anchor.budget
    )
    assert "BR(B_s -> phi gamma)" in result.diagnostics["exclusive_branching_formula"]
    assert result.diagnostics["helicity_observables_evaluated"] is False
    assert "A_Delta" in result.diagnostics["helicity_observable_note"]
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_constraint_matches_shared_bsgamma_core_recomputation():
    constraint = fcc.get(_PID)
    left = 0.03e-1 + 0.02e-1j
    right = 1.0e-2 - 0.3e-2j
    couplings = _bs_couplings(left=left, right=right)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))

    wilsons = bsgamma_core.compute_bsgamma_wilsons(
        couplings,
        inputs=constraint.sm_inputs,
    )
    direct = bsgamma_core.branching_fraction_from_c7(
        c7_np=wilsons.c7_np,
        c7p_np=wilsons.c7p_np,
        sm_branching_fraction=constraint.anchor.reference_value,
        inputs=constraint.sm_inputs,
    )

    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.diagnostics["c7_np"] == pytest.approx(wilsons.c7_np)
    assert result.diagnostics["c7p_np"] == pytest.approx(wilsons.c7p_np)
    assert wilsons.c7_np_matching == pytest.approx(left)
    assert wilsons.c7p_np_matching == pytest.approx(right)


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    couplings = _bs_couplings(left=2.0e-3 + 0.5e-3j, right=1.0e-2 - 0.3e-2j)
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
        "left_bs_coupling",
        "right_bs_coupling",
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
        "budget_combined_sigma",
        "hard_veto_np_shift_budget",
        "c7_running_from_c7",
        "c7_running_from_c8",
        "c8_running_from_c8",
        "alpha_s_matching_scale",
        "alpha_s_low_scale",
        "exclusive_normalization_branching_fraction",
        "reference_measurement_branching_fraction",
        "normalization_formula_branching_fraction",
        "budget_reference_measurement_sigma",
        "s_phigamma_uncertainty_combined",
        "a_delta_phigamma_uncertainty_combined_plus",
        "a_delta_phigamma_uncertainty_combined_minus",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["wilson_coefficients"]["C7_NP"] == pytest.approx(
        result.diagnostics["c7_np"]
    )
    assert result.diagnostics["bsgamma_rg_running_applied"] is True
    assert result.diagnostics["kk_ew_mass_extra_used"] is False
    assert result.diagnostics["down_sector_indices"] == (1, 2)
    assert result.diagnostics["s_phigamma_experimental"] == pytest.approx(0.43)
    assert result.diagnostics["a_delta_phigamma_experimental"] == pytest.approx(-0.67)
    assert result.diagnostics["a_delta_phigamma_uncertainty_stat_plus"] == pytest.approx(
        0.37
    )
    assert result.diagnostics["a_delta_phigamma_uncertainty_stat_minus"] == pytest.approx(
        0.41
    )
    assert result.diagnostics["a_delta_phigamma_uncertainty_combined_plus"] == pytest.approx(
        math.sqrt(0.37**2 + 0.17**2)
    )
    assert result.diagnostics["a_delta_phigamma_uncertainty_combined_minus"] == pytest.approx(
        math.sqrt(0.41**2 + 0.17**2)
    )
    assert result.diagnostics["budget_sm_theory_sigma"] is None
    assert result.diagnostics["sm_theory_prediction_available"] is False


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_bs_couplings(left=0.0j, right=1.0e-1), True),
        (_bs_couplings(left=0.0j, right=2.0e-1), False),
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
    couplings = _bs_couplings(left=0.0j, right=2.0e-1)
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
    couplings = _bs_couplings(left=2.0e-3 + 0.5e-3j, right=1.0e-2 - 0.3e-2j)
    before_left_down = couplings.left_down.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)
