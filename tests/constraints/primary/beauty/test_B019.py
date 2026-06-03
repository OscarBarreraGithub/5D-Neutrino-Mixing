"""Production tests for B019 (R_K* in B0 -> K*(892)0 ell ell)."""

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
from flavor_catalog_constraints.primary.beauty import B019 as b019_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.model import RotationParameters, ckm_like_unitary
from quarkConstraints.rare_b_kstar_dilepton import evaluate_b_to_kstar_mumu
from tests.rare_b_phase3d_helpers import (
    core_wilsons_from_rs_coeff,
    rs_coeff,
    sample_rare_b_point,
    scaled_rare_b_point,
    sm_limit_rare_b_point,
)

_PID = "B019"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B019.yaml"


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


def _manual_form_factors(q2: float, mode) -> tuple[float, float, float]:
    def pole(f0: float, pole_mass: float, slope: float) -> float:
        return float(
            f0
            * (1.0 + slope * q2 / mode.parent_mass_gev**2)
            / (1.0 - q2 / pole_mass**2)
        )

    return (
        pole(mode.v_0, mode.v_pole_mass_gev, mode.v_slope),
        pole(mode.a1_0, mode.a1_pole_mass_gev, mode.a1_slope),
        pole(mode.a2_0, mode.a2_pole_mass_gev, mode.a2_slope),
    )


def _manual_kstar_helicity_weight(q2: float, mode) -> float:
    m_b = mode.parent_mass_gev
    m_v = mode.daughter_mass_gev
    kallen = (
        m_b**4
        + m_v**4
        + q2**2
        - 2.0 * (m_b**2 * m_v**2 + m_b**2 * q2 + m_v**2 * q2)
    )
    if q2 <= 0.0 or kallen <= 0.0:
        return 0.0
    v_ff, a1_ff, a2_ff = _manual_form_factors(q2, mode)
    sqrt_kallen = math.sqrt(kallen)
    h_perp = sqrt_kallen * v_ff / (m_b + m_v)
    h_parallel = (m_b + m_v) * a1_ff
    h_long = (
        (m_b * m_b - m_v * m_v - q2) * (m_b + m_v) * a1_ff
        - kallen * a2_ff / (m_b + m_v)
    ) / (2.0 * m_v * math.sqrt(q2))
    return float(h_perp * h_perp + h_parallel * h_parallel + h_long * h_long)


def _manual_b_to_kstar_mumu(
    inputs,
    *,
    q2_min_gev2: float,
    q2_max_gev2: float,
    c9_total: complex,
    c10_total: complex,
    grid_points: int = 12001,
) -> float:
    sd = inputs.short_distance_inputs
    mode = inputs.bzero_kstarzero
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

    def dbr(q2: float) -> float:
        if q2 <= threshold:
            return 0.0
        kallen = (
            mode.parent_mass_gev**4
            + mode.daughter_mass_gev**4
            + q2**2
            - 2.0
            * (
                mode.parent_mass_gev**2 * mode.daughter_mass_gev**2
                + mode.parent_mass_gev**2 * q2
                + mode.daughter_mass_gev**2 * q2
            )
        )
        if kallen <= 0.0:
            return 0.0
        beta = math.sqrt(max(0.0, 1.0 - threshold / q2))
        lepton_mass_factor = beta * (1.0 + 2.0 * sd.muon_mass_gev**2 / q2)
        prefactor = (
            tau
            * sd.gf_gev_minus2**2
            * sd.alpha_em_mz**2
            * abs(lambda_t) ** 2
            / (3072.0 * math.pi**5 * mode.parent_mass_gev**3)
        )
        return float(
            prefactor
            * math.sqrt(kallen)
            * lepton_mass_factor
            * _manual_kstar_helicity_weight(q2, mode)
            * (abs(complex(c9_total)) ** 2 + abs(complex(c10_total)) ** 2)
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
    assert constraint.observable == "R_K* central-q2"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    obs = pdg["observables"]
    theory = pdg["theory_inputs"]
    low = _find_entry(obs, "name", "R_K* low-q2")
    active = _find_entry(obs, "name", "R_K* central-q2")
    old_low = _find_entry(obs, "name", "LHCb 2017 low-q2 R_K*")
    old_central = _find_entry(obs, "name", "LHCb 2017 central-q2 R_K*")
    sm = _find_entry(theory, "name", "BIP 2016 central-q2 SM R_K*")
    sm_low = _find_entry(theory, "name", "BIP 2016 low-q2 SM R_K* for LHCb 2023 bin")

    active_value = float(active["value"])
    active_sigma = 0.074
    sm_ratio = float(sm["value"])
    sm_sigma = 0.01

    assert constraint.anchor.rkstar_low.value == pytest.approx(float(low["value"]))
    assert constraint.anchor.rkstar_low.uncertainty == pytest.approx(0.097)
    assert constraint.anchor.rkstar_central.value == pytest.approx(active_value)
    assert constraint.anchor.rkstar_central.uncertainty == pytest.approx(active_sigma)
    assert constraint.anchor.rkstar_central.source_url == active["source_url"]
    assert constraint.anchor.rkstar_central.snapshot_path == active["snapshot_path"]
    assert constraint.anchor.rkstar_central.q2_region == active["q2_region"]
    assert constraint.anchor.lhcb_2017_low.value == pytest.approx(float(old_low["value"]))
    assert constraint.anchor.lhcb_2017_low.uncertainty == pytest.approx(
        0.5 * (math.sqrt(0.11**2 + 0.03**2) + math.sqrt(0.07**2 + 0.03**2))
    )
    assert constraint.anchor.lhcb_2017_central.value == pytest.approx(
        float(old_central["value"])
    )
    assert constraint.anchor.lhcb_2017_central.significance == old_central[
        "significance"
    ]
    assert constraint.anchor.sm_central.value == pytest.approx(sm_ratio)
    assert constraint.anchor.sm_central.uncertainty == pytest.approx(sm_sigma)
    assert constraint.anchor.sm_central.snapshot_path == sm["snapshot_path"]
    assert constraint.anchor.sm_low.value == pytest.approx(float(sm_low["value"]))
    assert constraint.anchor.sm_low.uncertainty == pytest.approx(0.014)
    assert constraint.anchor.q2_min_gev2 == pytest.approx(1.1)
    assert constraint.anchor.q2_max_gev2 == pytest.approx(6.0)
    assert constraint.sm_lfu_ratio == pytest.approx(1.0)
    assert constraint.anchor.budget_band.central_residual == pytest.approx(
        abs(active_value - sm_ratio)
    )
    assert constraint.anchor.budget_band.experimental_sigma == pytest.approx(
        active_sigma
    )
    assert constraint.anchor.budget_band.sm_theory_sigma == pytest.approx(sm_sigma)
    assert constraint.anchor.budget == pytest.approx(
        abs(active_value - sm_ratio) + active_sigma + sm_sigma
    )
    assert constraint.anchor.budget == pytest.approx(0.11200000000000002)
    assert "B019.yaml" in constraint.anchor.budget_band.source
    assert "sigma_SM_QED" in constraint.anchor.budget_band.construction

    with pytest.raises(AnchorError):
        b019_module._load_observable(
            process_id=_PID,
            observable_name="no such B019 observable",
        )
    with pytest.raises(AnchorError):
        b019_module._load_theory(
            process_id=_PID,
            theory_name="no such B019 theory row",
        )


def test_observable_and_theory_anchors_route_through_scaffold_load_anchor_and_fail_loudly(
    monkeypatch,
):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b019_module.anchor_scaffold.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b019_module.anchor_scaffold, "load_anchor", spy_load_anchor)
    observables = b019_module._load_b019_observables(_PID)

    assert calls == [
        ("pdg_or_equivalent.observables[0]",),
        ("pdg_or_equivalent.observables[1]",),
        ("pdg_or_equivalent.observables[2]",),
        ("pdg_or_equivalent.observables[3]",),
        ("pdg_or_equivalent.theory_inputs[0]",),
        ("pdg_or_equivalent.theory_inputs[1]",),
    ]
    assert observables.rkstar_central.value == pytest.approx(fcc.get(_PID).anchor.value)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent.observables[99]",
            value=1.0,
            uncertainty=0.1,
        )

    monkeypatch.setattr(
        b019_module.anchor_scaffold,
        "load_anchor",
        mismatched_load_anchor,
    )
    with pytest.raises(AnchorError, match="load_anchor selected"):
        b019_module._load_b019_observables(_PID)


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
    assert result.sm_prediction == pytest.approx(1.0)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert result.diagnostics["sm_lfu_ratio"] == pytest.approx(1.0)
    assert result.diagnostics["sm_theory_lfu_ratio"] == pytest.approx(1.0)
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert old_style_result.passes is True
    assert old_style_result.predicted is None
    assert old_style_result.diagnostics["evaluated"] is False
    assert old_style_result.diagnostics["missing_extra"] == "rs_semileptonic_wilsons"
    assert old_style_result.diagnostics["legacy_quark_mass_basis_couplings_present"] is True


def test_sm_limit_rkstar_matches_independent_integral_and_yaml_budget():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(sm_limit_rare_b_point())
    sd = constraint.sm_inputs.short_distance_inputs
    manual = _manual_b_to_kstar_mumu(
        constraint.sm_inputs,
        q2_min_gev2=constraint.anchor.q2_min_gev2,
        q2_max_gev2=constraint.anchor.q2_max_gev2,
        c9_total=constraint.sm_inputs.c9_sm,
        c10_total=sd.c10_sm,
    )

    assert result.diagnostics["rkstar_proxy_numerator_branching_fraction"] == (
        pytest.approx(manual, rel=3.0e-4)
    )
    assert result.diagnostics["rkstar_proxy_denominator_branching_fraction"] == (
        pytest.approx(manual, rel=3.0e-4)
    )
    assert result.diagnostics["rkstar_proxy_denominator_branching_fraction"] == (
        pytest.approx(3.676601736716339e-08)
    )
    assert result.predicted == pytest.approx(1.0)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.diagnostics["sm_theory_lfu_ratio"] == pytest.approx(1.0)
    assert abs(result.predicted - constraint.anchor.sm_value) <= (
        constraint.anchor.sm_central.uncertainty
    )
    assert result.ratio == pytest.approx(
        abs(result.predicted - constraint.anchor.value) / constraint.anchor.budget
    )
    assert result.ratio == pytest.approx(0.25)
    assert result.passes is True
    assert result.diagnostics["active_experimental_raw_value"] == pytest.approx(1.028)


def test_np_prediction_matches_underlying_core_recomputation():
    constraint = fcc.get(_PID)
    point = scaled_rare_b_point(scale=100.0, lepton_scales={"e": 0.5})
    result = constraint.evaluate(point)
    direct_mu = evaluate_b_to_kstar_mumu(
        core_wilsons_from_rs_coeff(
            rs_coeff(point, transition="b_s", lepton="mu"),
            matching_scale_gev=point.extras["kk_ew_mass_gev"],
        ),
        q2_min_gev2=constraint.anchor.q2_min_gev2,
        q2_max_gev2=constraint.anchor.q2_max_gev2,
        inputs=constraint.sm_inputs,
    )
    direct_e = evaluate_b_to_kstar_mumu(
        core_wilsons_from_rs_coeff(
            rs_coeff(point, transition="b_s", lepton="e"),
            matching_scale_gev=point.extras["kk_ew_mass_gev"],
        ),
        q2_min_gev2=constraint.anchor.q2_min_gev2,
        q2_max_gev2=constraint.anchor.q2_max_gev2,
        inputs=constraint.sm_inputs,
    )
    expected = direct_mu.branching_fraction / direct_e.branching_fraction

    assert result.predicted == pytest.approx(expected)
    assert result.sm_prediction == pytest.approx(1.0)
    assert result.diagnostics["rkstar_proxy_numerator_branching_fraction"] == (
        pytest.approx(direct_mu.branching_fraction)
    )
    assert result.diagnostics["rkstar_proxy_denominator_branching_fraction"] == (
        pytest.approx(direct_e.branching_fraction)
    )
    assert result.diagnostics["c9_vector_np"] == pytest.approx(direct_mu.c9_vector_np)
    assert result.diagnostics["c10_axial_np"] == pytest.approx(direct_mu.c10_axial_np)
    assert result.diagnostics["denominator_wilson_coefficients"]["C9_NP"] == (
        pytest.approx(direct_e.diagnostics["c9_np"])
    )


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
        "helicity_weight_at_bin_center",
        "v_0",
        "a1_0",
        "a2_0",
        "budget_central_residual",
        "budget_experimental_sigma",
        "budget_sm_lfu_ratio",
        "budget_sm_theory_sigma",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["lepton_universal_rkstar_proxy"] == pytest.approx(1.0)
    assert "lepton-universal" in result.diagnostics[
        "lepton_universal_cancellation_note"
    ]
    assert result.diagnostics["rs_semileptonic_wilsons_present"] is True
    assert result.diagnostics["rs_semileptonic_matching_status"].endswith(
        "no_second_1_over_M_KK_squared"
    )
    assert result.diagnostics["denominator_wilson_coefficients"]
    assert result.diagnostics["c7_nonlocal_charm_omitted"] is True
    assert result.diagnostics["qed_lepton_mass_omitted"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("point", "expected_pass"),
    [
        (scaled_rare_b_point(scale=1.0e5), True),
        (scaled_rare_b_point(scale=1.0e5, lepton_scales={"e": 0.0}), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    point,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point)

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
        assert result.predicted == pytest.approx(1.0)
        assert result.diagnostics["rkstar_proxy_muon_ratio_to_sm"] == pytest.approx(
            result.diagnostics["rkstar_proxy_electron_ratio_to_sm"]
        )
    else:
        assert result.ratio > 1.0
        assert result.predicted != pytest.approx(1.0)


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    default_point = scaled_rare_b_point(scale=100.0, lepton_scales={"e": 0.5})
    ew_point = scaled_rare_b_point(
        scale=100.0,
        lepton_scales={"e": 0.5},
        kk_ew_mass_gev=6000.0,
    )
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
