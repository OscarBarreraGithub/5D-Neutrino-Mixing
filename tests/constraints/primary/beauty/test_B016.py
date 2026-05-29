"""Production tests for B016 (exclusive B -> K ell ell)."""

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
from flavor_catalog_constraints.primary.beauty import B016 as b016_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.model import RotationParameters, ckm_like_unitary
from quarkConstraints.rare_b_dilepton import evaluate_b_to_k_mumu

_PID = "B016"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B016.yaml"


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


def _manual_fplus(q2: float, mode) -> float:
    t_plus = (mode.parent_mass_gev + mode.daughter_mass_gev) ** 2
    t_minus = (mode.parent_mass_gev - mode.daughter_mass_gev) ** 2
    t_0 = t_plus * (1.0 - math.sqrt(1.0 - t_minus / t_plus))

    def z(value: float) -> float:
        return (
            math.sqrt(t_plus - value) - math.sqrt(t_plus - t_0)
        ) / (
            math.sqrt(t_plus - value) + math.sqrt(t_plus - t_0)
        )

    z_q = z(q2)
    z_0 = z(0.0)
    shape = 1.0 + mode.bcl_a1 * (z_q - z_0) + mode.bcl_a2 * (
        z_q * z_q - z_0 * z_0
    )
    return float(mode.fplus_0 * shape / (1.0 - q2 / mode.pole_mass_gev**2))


def _manual_b_to_k_mumu(
    inputs,
    *,
    q2_min_gev2: float,
    q2_max_gev2: float,
    c9_total: complex,
    c10_total: complex,
    grid_points: int = 12001,
) -> float:
    sd = inputs.short_distance_inputs
    mode = inputs.bplus_kplus
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
        fplus = _manual_fplus(q2, mode)
        prefactor = (
            tau
            * sd.gf_gev_minus2**2
            * sd.alpha_em_mz**2
            * abs(lambda_t) ** 2
            / (1536.0 * math.pi**5 * mode.parent_mass_gev**3)
        )
        lepton_mass_factor = beta * (1.0 + 2.0 * sd.muon_mass_gev**2 / q2)
        return float(
            prefactor
            * kallen ** 1.5
            * lepton_mass_factor
            * (
                abs(complex(c9_total) * fplus) ** 2
                + abs(complex(c10_total) * fplus) ** 2
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
    assert constraint.observable == "BR(B+ -> K+ ell+ ell-)"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    charged = _find_entry(
        pdg["observables"],
        "name",
        "BR(B+ -> K+ ell+ ell-)",
    )
    charged_value = float(charged["value"]) * 1.0e-7
    charged_sigma = 0.40e-7
    central = abs(charged_value - constraint.sm_result.branching_fraction)
    proxy_sigma = (
        b016_module._PROXY_THEORY_UNCERTAINTY_FRACTION
        * constraint.sm_result.branching_fraction
    )

    assert constraint.anchor.charged.value == pytest.approx(charged_value)
    assert constraint.anchor.charged.uncertainty == pytest.approx(charged_sigma)
    assert constraint.anchor.charged.source_url == charged["source_url"]
    assert constraint.anchor.charged.snapshot_path == charged["snapshot_path"]
    assert constraint.anchor.charged.year == charged["year"]
    assert not hasattr(constraint.anchor, "neutral")
    assert constraint.anchor.low_q2_formula_benchmark.label == (
        "internal_c9_c10_formula_benchmark"
    )
    assert constraint.anchor.low_q2_formula_benchmark.q2_min_gev2 == pytest.approx(
        1.1
    )
    assert constraint.anchor.low_q2_formula_benchmark.q2_max_gev2 == pytest.approx(
        6.0
    )
    assert "not a data validation" in constraint.anchor.low_q2_formula_benchmark.note
    assert constraint.anchor.budget_band.central_residual == pytest.approx(central)
    assert constraint.anchor.budget_band.experimental_sigma == pytest.approx(
        charged_sigma
    )
    assert constraint.anchor.budget_band.proxy_theory_sigma == pytest.approx(
        proxy_sigma
    )
    assert constraint.anchor.budget == pytest.approx(
        central + charged_sigma + proxy_sigma
    )
    assert constraint.anchor.budget == pytest.approx(2.1348703212043188e-07)
    assert "NEEDS-HUMAN-PHYSICS" in constraint.anchor.budget_band.construction

    with pytest.raises(AnchorError):
        b016_module._load_branching_observable(
            process_id=_PID,
            observable_name="no such B016 observable",
        )


def test_observable_anchors_route_through_scaffold_load_anchor_and_fail_loudly(
    monkeypatch,
):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b016_module.anchor_scaffold.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b016_module.anchor_scaffold, "load_anchor", spy_load_anchor)
    anchor = b016_module._load_b016_anchor(
        _PID,
        formula_sm=fcc.get(_PID).sm_result.branching_fraction,
    )

    assert calls == [
        ("pdg_or_equivalent.observables[0]",),
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
        b016_module.anchor_scaffold,
        "load_anchor",
        mismatched_load_anchor,
    )
    with pytest.raises(AnchorError, match="load_anchor selected"):
        b016_module._load_b016_anchor(
            _PID,
            formula_sm=fcc.get(_PID).sm_result.branching_fraction,
        )


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

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
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_sm_limit_branching_fraction_matches_independent_integral_and_formula_benchmark():
    constraint = fcc.get(_PID)
    couplings = _sb_couplings(left=0.0j, right=0.0j)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))
    sd = constraint.sm_inputs.short_distance_inputs
    mode = constraint.sm_inputs.bplus_kplus
    threshold = 4.0 * sd.muon_mass_gev**2
    endpoint = mode.q2_max_gev2
    manual_full = _manual_b_to_k_mumu(
        constraint.sm_inputs,
        q2_min_gev2=threshold,
        q2_max_gev2=endpoint,
        c9_total=constraint.sm_inputs.c9_sm,
        c10_total=sd.c10_sm,
    )
    manual_low = _manual_b_to_k_mumu(
        constraint.sm_inputs,
        q2_min_gev2=1.1,
        q2_max_gev2=6.0,
        c9_total=constraint.sm_inputs.c9_sm,
        c10_total=sd.c10_sm,
    )

    assert result.predicted == pytest.approx(manual_full, rel=3.0e-4)
    assert result.predicted == pytest.approx(5.750185255422402e-07)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert abs(result.predicted - constraint.anchor.value) < constraint.anchor.uncertainty
    assert result.diagnostics[
        "low_q2_formula_benchmark_branching_fraction"
    ] == pytest.approx(manual_low, rel=3.0e-4)
    assert result.diagnostics[
        "low_q2_formula_benchmark_branching_fraction"
    ] == pytest.approx(1.851950637006692e-07)
    assert result.diagnostics["low_q2_formula_benchmark_label"] == (
        "internal_c9_c10_formula_benchmark"
    )
    assert "not a data validation" in result.diagnostics[
        "low_q2_formula_benchmark_note"
    ]
    assert "low_q2_validation_source_key" not in result.diagnostics
    assert result.ratio == pytest.approx(
        abs(result.predicted - constraint.anchor.value) / constraint.anchor.budget
    )
    assert result.passes is True


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    couplings = _sb_couplings(left=1.0e-2 + 0.2e-2j, right=0.5e-2j)
    point = point_builder.build_from_quark_couplings(couplings)
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
        "fplus_0",
        "fplus_q2_mid",
        "budget_central_residual",
        "budget_experimental_sigma",
        "budget_proxy_theory_sigma",
        "budget_proxy_theory_fraction",
        "low_q2_formula_benchmark_branching_fraction",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["neutral_mode_scope"] == (
        "out_of_scope_separate_entry_future"
    )
    assert "neutral_hflav_branching_fraction" not in result.diagnostics
    assert "neutral_sm_formula_branching_fraction" not in result.diagnostics
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["budget_proxy_theory_rationale"]
    assert result.diagnostics["c7_nonlocal_charm_omitted"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_sb_couplings(left=1.0e-2), True),
        (_sb_couplings(left=5.0e-1), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct = evaluate_b_to_k_mumu(couplings, inputs=constraint.sm_inputs)

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
    couplings = _sb_couplings(left=1.0e-1)
    default_point = point_builder.build_from_quark_couplings(couplings)
    ew_point = point_builder.make_point(
        quark_mass_basis_couplings=couplings,
        kk_ew_mass_gev=6000.0,
    )
    default_result = fcc.get(_PID).evaluate(default_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert abs(ew_result.diagnostics["c9_vector_np"]) == pytest.approx(
        abs(default_result.diagnostics["c9_vector_np"]) / 4.0
    )
    assert abs(ew_result.diagnostics["c10_axial_np"]) == pytest.approx(
        abs(default_result.diagnostics["c10_axial_np"]) / 4.0
    )


def test_evaluate_is_pure_and_deterministic():
    couplings = _sb_couplings(left=1.0e-2 + 0.2e-2j, right=0.5e-2j)
    before_left_down = couplings.left_down.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_down, before_left_down)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)
