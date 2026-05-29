"""Production tests for B018 (R_K in B+ -> K+ ell ell)."""

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
from flavor_catalog_constraints.primary.beauty import B018 as b018_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.model import RotationParameters, ckm_like_unitary
from quarkConstraints.rare_b_dilepton import evaluate_b_to_k_mumu

_PID = "B018"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B018.yaml"


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
    assert constraint.observable == "R_K central-q2"


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    obs = _yaml_pdg_block()["observables"]
    low = _find_entry(obs, "name", "R_K low-q2")
    active = _find_entry(obs, "name", "R_K central-q2")
    full = _find_entry(obs, "name", "R_K full-q2 Belle-only")
    superseded = _find_entry(obs, "name", "superseded LHCb 2021 central-q2 R_K")

    active_value = float(active["value"])
    active_sigma = 0.047
    sm_ratio = constraint.sm_lfu_ratio

    assert constraint.anchor.rk_low.value == pytest.approx(float(low["value"]))
    assert constraint.anchor.rk_low.uncertainty == pytest.approx(0.090)
    assert constraint.anchor.rk_central.value == pytest.approx(active_value)
    assert constraint.anchor.rk_central.uncertainty == pytest.approx(active_sigma)
    assert constraint.anchor.rk_central.source_url == active["source_url"]
    assert constraint.anchor.rk_central.snapshot_path == active["snapshot_path"]
    assert constraint.anchor.rk_central.q2_region == active["q2_region"]
    assert constraint.anchor.rk_full_belle.value == pytest.approx(float(full["value"]))
    assert constraint.anchor.rk_full_belle.uncertainty == pytest.approx(0.16)
    assert constraint.anchor.superseded_rk_2021.value == pytest.approx(
        float(superseded["value"])
    )
    assert constraint.anchor.superseded_rk_2021.significance == superseded[
        "significance"
    ]
    assert constraint.anchor.q2_min_gev2 == pytest.approx(1.1)
    assert constraint.anchor.q2_max_gev2 == pytest.approx(6.0)
    assert sm_ratio == pytest.approx(1.0)
    assert constraint.anchor.budget_band.central_residual == pytest.approx(
        abs(active_value - sm_ratio)
    )
    assert constraint.anchor.budget_band.experimental_sigma == pytest.approx(
        active_sigma
    )
    assert constraint.anchor.budget == pytest.approx(
        abs(active_value - sm_ratio) + active_sigma
    )
    assert constraint.anchor.budget == pytest.approx(0.10000000000000005)
    assert "B018.yaml" in constraint.anchor.budget_band.source
    assert "no B016 proxy-theory envelope" in constraint.anchor.budget_band.construction

    with pytest.raises(AnchorError):
        b018_module._load_observable(
            process_id=_PID,
            observable_name="no such B018 observable",
        )


def test_observable_anchors_route_through_scaffold_load_anchor_and_fail_loudly(
    monkeypatch,
):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b018_module.anchor_scaffold.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b018_module.anchor_scaffold, "load_anchor", spy_load_anchor)
    observables = b018_module._load_b018_observables(_PID)

    assert calls == [
        ("pdg_or_equivalent.observables[0]",),
        ("pdg_or_equivalent.observables[1]",),
        ("pdg_or_equivalent.observables[2]",),
        ("pdg_or_equivalent.observables[3]",),
    ]
    assert observables.rk_central.value == pytest.approx(fcc.get(_PID).anchor.value)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent.observables[99]",
            value=1.0,
            uncertainty=0.1,
        )

    monkeypatch.setattr(
        b018_module.anchor_scaffold,
        "load_anchor",
        mismatched_load_anchor,
    )
    with pytest.raises(AnchorError, match="load_anchor selected"):
        b018_module._load_b018_observables(_PID)


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.sm_prediction == pytest.approx(1.0)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert result.diagnostics["sm_lfu_ratio"] == pytest.approx(1.0)
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_sm_limit_rk_matches_independent_integral_and_yaml_budget():
    constraint = fcc.get(_PID)
    couplings = _sb_couplings(left=0.0j, right=0.0j)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))
    sd = constraint.sm_inputs.short_distance_inputs
    manual = _manual_b_to_k_mumu(
        constraint.sm_inputs,
        q2_min_gev2=constraint.anchor.q2_min_gev2,
        q2_max_gev2=constraint.anchor.q2_max_gev2,
        c9_total=constraint.sm_inputs.c9_sm,
        c10_total=sd.c10_sm,
    )

    assert result.diagnostics["rk_proxy_numerator_branching_fraction"] == (
        pytest.approx(manual, rel=3.0e-4)
    )
    assert result.diagnostics["rk_proxy_denominator_branching_fraction"] == (
        pytest.approx(manual, rel=3.0e-4)
    )
    assert result.diagnostics["rk_proxy_denominator_branching_fraction"] == (
        pytest.approx(1.851950637006692e-07)
    )
    assert result.predicted == pytest.approx(1.0)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.ratio == pytest.approx(
        abs(result.predicted - constraint.anchor.value) / constraint.anchor.budget
    )
    assert result.ratio < 1.0
    assert result.diagnostics["active_experimental_raw_value"] == pytest.approx(0.947)


def test_np_prediction_matches_underlying_core_recomputation():
    constraint = fcc.get(_PID)
    couplings = _sb_couplings(left=0.03e-1 + 0.02e-1j, right=1.0e-2 - 0.3e-2j)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))
    direct = evaluate_b_to_k_mumu(
        couplings,
        q2_min_gev2=constraint.anchor.q2_min_gev2,
        q2_max_gev2=constraint.anchor.q2_max_gev2,
        inputs=constraint.sm_inputs,
    )
    expected = (
        direct.branching_fraction / constraint.sm_electron_proxy_branching_fraction
    )

    assert result.predicted == pytest.approx(expected)
    assert result.sm_prediction == pytest.approx(1.0)
    assert result.diagnostics["rk_proxy_numerator_branching_fraction"] == (
        pytest.approx(direct.branching_fraction)
    )
    assert result.diagnostics["c9_vector_np"] == pytest.approx(direct.c9_vector_np)
    assert result.diagnostics["c10_axial_np"] == pytest.approx(direct.c10_axial_np)


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
        "q2_min_gev2",
        "q2_max_gev2",
        "rk_proxy_numerator_branching_fraction",
        "rk_proxy_denominator_branching_fraction",
        "budget_central_residual",
        "budget_experimental_sigma",
        "budget_sm_lfu_ratio",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["lepton_universal_rk_proxy"] == pytest.approx(1.0)
    assert "lepton-universal" in result.diagnostics[
        "lepton_universal_cancellation_note"
    ]
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["c7_nonlocal_charm_omitted"] is True


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_sb_couplings(left=1.0e-2), True),
        (_sb_couplings(left=1.0e-1), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct = evaluate_b_to_k_mumu(
        couplings,
        q2_min_gev2=constraint.anchor.q2_min_gev2,
        q2_max_gev2=constraint.anchor.q2_max_gev2,
        inputs=constraint.sm_inputs,
    )

    assert result.predicted == pytest.approx(
        direct.branching_fraction / constraint.sm_electron_proxy_branching_fraction
    )
    assert result.passes is expected_pass
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
