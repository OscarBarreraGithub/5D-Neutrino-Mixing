"""Production tests for B015 (inclusive B -> X_s ell ell)."""

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
from flavor_catalog_constraints.primary.beauty import B015 as b015_module
from quarkConstraints.bsgamma import compute_bsgamma_wilsons
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_b_dilepton import evaluate_inclusive_b_to_xs_mumu
from tests.rare_b_phase3d_helpers import (
    core_wilsons_from_rs_coeff,
    rs_coeff,
    sample_rare_b_point,
    scaled_rare_b_point,
    sm_limit_rare_b_point,
)

_PID = "B015"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B015.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _find_entry(entries, key: str, expected: str):
    for entry in entries:
        if str(entry.get(key)) == expected:
            return entry
    raise AssertionError(f"missing {key}={expected!r}")


def _bs_couplings(
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


def _manual_kernel(
    q2: float,
    *,
    mb: float,
    c7: complex,
    c7p: complex,
    c9: complex,
    c9p: complex,
    c10: complex,
    c10p: complex,
) -> float:
    shat = q2 / mb**2
    phase = (1.0 - shat) ** 2
    vector_axial = (1.0 + 2.0 * shat) * (
        abs(c9) ** 2 + abs(c10) ** 2 + abs(c9p) ** 2 + abs(c10p) ** 2
    )
    dipole = 4.0 * (1.0 + 2.0 / shat) * (abs(c7) ** 2 + abs(c7p) ** 2)
    interference = 12.0 * (
        (c7 * c9.conjugate()).real + (c7p * c9p.conjugate()).real
    )
    return float(max(0.0, phase * (vector_axial + dipole + interference)))


def _manual_integral(
    inputs,
    *,
    q2_min: float,
    q2_max: float,
    c7: complex,
    c7p: complex,
    c9: complex,
    c9p: complex,
    c10: complex,
    c10p: complex,
    grid_points: int = 12001,
) -> float:
    xs = np.linspace(q2_min, q2_max, grid_points)
    ys = np.array(
        [
            _manual_kernel(
                float(x),
                mb=inputs.partonic_b_mass_gev,
                c7=c7,
                c7p=c7p,
                c9=c9,
                c9p=c9p,
                c10=c10,
                c10p=c10p,
            )
            for x in xs
        ]
    )
    return float(np.trapezoid(ys, xs))


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "beauty"
    assert constraint.observable == "BR(B -> X_s ell+ ell-) low-q2 inclusive"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    obs = pdg["observables"]
    total = _find_entry(obs, "name", "BR(B -> X_s ell+ ell-)")
    low = _find_entry(obs, "name", "Experimental low-q2 weighted average")
    high = _find_entry(obs, "name", "Experimental high-q2 weighted average")
    sm_mumu = _find_entry(obs, "name", "SM BR[1,6]_mumu")
    sm_ee = _find_entry(obs, "name", "SM BR[1,6]_ee")
    low_value = float(low["value"]) * 1.0e-6
    low_sigma = 0.37e-6
    sm_value = float(sm_mumu["value"]) * 1.0e-6
    sm_sigma = 0.09e-6
    central = abs(low_value - sm_value)
    combined = math.sqrt(low_sigma**2 + sm_sigma**2)
    proxy = (
        b015_module.RARE_B_DILEPTON_INCLUSIVE_XS_PROXY_THEORY_UNCERTAINTY_FRACTION
        * sm_value
    )

    assert constraint.anchor.total_hflav.value == pytest.approx(
        float(total["value"]) * 1.0e-6
    )
    assert constraint.anchor.low_q2_experimental.value == pytest.approx(low_value)
    assert constraint.anchor.low_q2_experimental.uncertainty == pytest.approx(
        low_sigma
    )
    assert constraint.anchor.low_q2_experimental.source_url == low["source_url"]
    assert constraint.anchor.low_q2_experimental.snapshot_path == low["snapshot_path"]
    assert constraint.anchor.high_q2_experimental.value == pytest.approx(
        float(high["value"]) * 1.0e-6
    )
    assert constraint.anchor.low_q2_sm_mumu.value == pytest.approx(sm_value)
    assert constraint.anchor.low_q2_sm_mumu.uncertainty == pytest.approx(sm_sigma)
    assert constraint.anchor.low_q2_sm_ee.value == pytest.approx(
        float(sm_ee["value"]) * 1.0e-6
    )
    assert constraint.anchor.q2_min_gev2 == pytest.approx(1.0)
    assert constraint.anchor.q2_max_gev2 == pytest.approx(6.0)
    assert constraint.anchor.budget_band.central_residual == pytest.approx(central)
    assert constraint.anchor.budget_band.combined_exp_sm_sigma == pytest.approx(
        combined
    )
    assert constraint.anchor.budget_band.proxy_theory_sigma == pytest.approx(proxy)
    assert constraint.anchor.budget == pytest.approx(central + combined + proxy)
    assert constraint.anchor.budget == pytest.approx(9.067886552931955e-07)
    assert "NEEDS-HUMAN-PHYSICS" in constraint.anchor.budget_band.construction

    with pytest.raises(AnchorError):
        b015_module._load_branching_observable(
            process_id=_PID,
            observable_name="no such B015 observable",
        )


def test_observable_anchors_route_through_scaffold_load_anchor_and_fail_loudly(
    monkeypatch,
):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b015_module.anchor_scaffold.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b015_module.anchor_scaffold, "load_anchor", spy_load_anchor)
    anchor = b015_module._load_b015_anchor(_PID)

    assert calls == [
        ("pdg_or_equivalent.observables[0]",),
        ("pdg_or_equivalent.observables[1]",),
        ("pdg_or_equivalent.observables[4]",),
        ("pdg_or_equivalent.observables[5]",),
        ("pdg_or_equivalent.observables[8]",),
        ("pdg_or_equivalent.observables[7]",),
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
        b015_module.anchor_scaffold,
        "load_anchor",
        mismatched_load_anchor,
    )
    with pytest.raises(AnchorError, match="load_anchor selected"):
        b015_module._load_b015_anchor(_PID)


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())
    old_style_result = constraint.evaluate(
        point_builder.build_from_quark_couplings(_bs_couplings(left=1.0e-2))
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


def test_sm_limit_branching_fraction_matches_independent_integral_and_yaml_sm():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(sm_limit_rare_b_point())
    inputs = constraint.sm_inputs
    q2_min = constraint.anchor.q2_min_gev2
    q2_max = constraint.anchor.q2_max_gev2
    sm_integral = _manual_integral(
        inputs,
        q2_min=q2_min,
        q2_max=q2_max,
        c7=inputs.dipole_inputs.c7_sm_eff,
        c7p=inputs.dipole_inputs.c7p_sm_eff,
        c9=inputs.c9_sm,
        c9p=0.0j,
        c10=inputs.short_distance_inputs.c10_sm,
        c10p=0.0j,
    )
    manual = constraint.anchor.sm_value * sm_integral / sm_integral

    assert result.predicted == pytest.approx(manual)
    assert result.predicted == pytest.approx(1.62e-06)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.sm_prediction == pytest.approx(constraint.anchor.sm_value)
    assert result.diagnostics["sm_formula_minus_anchor"] == pytest.approx(0.0)
    assert result.diagnostics["sm_shape_integral"] == pytest.approx(
        sm_integral,
        rel=3.0e-4,
    )
    assert result.diagnostics["c7_included"] is True
    assert result.diagnostics["c7_np"] == pytest.approx(0.0j)
    assert result.diagnostics["c9_np"] == pytest.approx(0.0j)
    assert result.diagnostics["c10_np"] == pytest.approx(0.0j)
    assert result.diagnostics["c9_c10_proxy_reused"] is False
    assert result.diagnostics["c9_c10_rs_semileptonic_rewired"] is True
    assert result.diagnostics["wilson_prefactor_reused"] is False
    assert result.diagnostics["second_mkk_suppression_applied"] is False
    assert result.passes is True


def test_np_prediction_matches_core_recomputation():
    constraint = fcc.get(_PID)
    point = scaled_rare_b_point(scale=100.0)
    result = constraint.evaluate(point)
    direct = evaluate_inclusive_b_to_xs_mumu(
        core_wilsons_from_rs_coeff(
            rs_coeff(point, transition="b_s", lepton="mu"),
            matching_scale_gev=point.extras["kk_ew_mass_gev"],
        ),
        sm_branching_fraction=constraint.anchor.sm_value,
        q2_min_gev2=constraint.anchor.q2_min_gev2,
        q2_max_gev2=constraint.anchor.q2_max_gev2,
        inputs=constraint.sm_inputs,
    )

    assert result.predicted == pytest.approx(direct.branching_fraction)
    assert result.sm_prediction == pytest.approx(direct.sm_branching_fraction)
    assert result.ratio == pytest.approx(
        abs(direct.branching_fraction - constraint.anchor.value)
        / constraint.anchor.budget
    )
    assert result.diagnostics["c7_np"] == pytest.approx(direct.c7_np)
    assert result.diagnostics["c9_np"] == pytest.approx(direct.c9_np)
    assert result.diagnostics["c10_np"] == pytest.approx(direct.c10_np)


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
        "c7_total",
        "c9_total",
        "c10_total",
        "c7_np",
        "c9_np",
        "c10_np",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "q2_min_gev2",
        "q2_max_gev2",
        "partonic_b_mass_gev",
        "sm_shape_integral",
        "total_shape_integral",
        "budget_central_residual",
        "budget_experimental_sigma",
        "budget_sm_theory_sigma",
        "budget_proxy_theory_sigma",
        "budget_proxy_theory_fraction",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["high_q2_scope"] == (
        "diagnostic_only_no_matching_SM_high_q2_anchor_in_B015_yaml"
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics[
        "budget_proxy_theory_rationale"
    ]
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert result.diagnostics["rs_semileptonic_matching_status"].endswith(
        "no_second_1_over_M_KK_squared"
    )
    assert result.diagnostics["c9_c10_rs_semileptonic_rewired"] is True
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

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_c7_dipole_path_is_untouched_when_optional_dipole_source_is_present():
    constraint = fcc.get(_PID)
    dipole_couplings = _bs_couplings(left=1.0e-1, right=2.0e-2j)
    no_dipole_point = scaled_rare_b_point(scale=20.0)
    dipole_point = scaled_rare_b_point(
        scale=20.0,
        quark_mass_basis_couplings=dipole_couplings,
    )
    no_dipole_result = constraint.evaluate(no_dipole_point)
    dipole_result = constraint.evaluate(dipole_point)
    direct_dipole = compute_bsgamma_wilsons(
        dipole_couplings,
        m_kk_gev=dipole_point.extras["kk_ew_mass_gev"],
        inputs=constraint.sm_inputs.dipole_inputs,
    )

    assert no_dipole_result.diagnostics["c7_np"] == pytest.approx(0.0j)
    assert "dipole_wilson_coefficients" not in no_dipole_result.diagnostics
    assert dipole_result.diagnostics["c7_np"] == pytest.approx(direct_dipole.c7_np)
    assert dipole_result.diagnostics["c7_np_matching"] == pytest.approx(
        direct_dipole.c7_np_matching
    )
    assert dipole_result.diagnostics["dipole_wilson_coefficients"]["C7_NP"] == (
        pytest.approx(direct_dipole.c7_np)
    )
    assert dipole_result.diagnostics["c9_np"] == pytest.approx(
        no_dipole_result.diagnostics["c9_np"]
    )
    assert dipole_result.diagnostics["c10_np"] == pytest.approx(
        no_dipole_result.diagnostics["c10_np"]
    )
    assert dipole_result.diagnostics["c9_c10_rs_semileptonic_rewired"] is True


def test_evaluate_is_pure_and_deterministic():
    point = sample_rare_b_point()
    before_bundle = point.extras["rs_semileptonic_wilsons"]
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras["rs_semileptonic_wilsons"] is before_bundle
