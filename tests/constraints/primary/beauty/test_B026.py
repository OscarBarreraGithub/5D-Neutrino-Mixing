"""Production tests for B026 (R_D* in B -> D* tau nu)."""

from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.beauty import B026 as b026_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.semileptonic_lfu import evaluate_rd_lfu_ratio

_PID = "B026"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B026.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _bc_proxy_couplings(
    *,
    charm_left: complex = 1.0 + 0.0j,
    bottom_right: complex = 0.0j,
    M_KK: float = 3000.0,
    g_s: float = 1.0,
) -> QuarkMassBasisCouplings:
    """Minimal valid mass-basis couplings for the c_L-b_R proxy slots."""

    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_up = zeros.copy()
    right_down = zeros.copy()
    left_up[1, 1] = charm_left
    right_down[2, 2] = bottom_right
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=g_s,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=left_up,
        left_down=zeros.copy(),
        right_up=zeros.copy(),
        right_down=right_down,
    )


def _manual_rdstar_proxy(
    couplings: QuarkMassBasisCouplings,
    *,
    sm_rdstar: float,
    m_kk_gev: float | None = None,
    bottom_mass_gev: float = 4.183,
    tau_mass_gev: float = 1.77686,
    gf_gev_minus2: float = 1.1663787e-5,
) -> tuple[float, complex]:
    resolved_mkk = float(couplings.M_KK if m_kk_gev is None else m_kk_gev)
    charm_left_overlap = complex(couplings.left_up[1, 1]) / float(couplings.g_s)
    bottom_right_overlap = complex(couplings.right_down[2, 2]) / float(couplings.g_s)
    scalar_shift = (
        charm_left_overlap
        * bottom_right_overlap
        * bottom_mass_gev
        * tau_mass_gev
        / (2.0 * math.sqrt(2.0) * gf_gev_minus2 * resolved_mkk**2)
    )
    return float(sm_rdstar * abs(1.0 + scalar_shift) ** 2), complex(scalar_shift)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "beauty"
    assert (
        constraint.observable
        == "R_D* = Gamma(B -> D* tau nu) / Gamma(B -> D* l nu)"
    )


def test_anchor_matches_yaml_and_budget_band():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["canonical_average"]
    sm = pdg["sm_reference_used_by_hflav"]
    joint = pdg["joint_fit_context"]
    central = abs(float(exp["value"]) - float(sm["rdstar_value"]))
    combined = math.sqrt(
        float(exp["uncertainty"]) ** 2 + float(sm["rdstar_uncertainty"]) ** 2
    )

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.uncertainty == pytest.approx(exp["uncertainty"])
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.standard_model.value == pytest.approx(sm["rdstar_value"])
    assert constraint.anchor.standard_model.uncertainty == pytest.approx(
        sm["rdstar_uncertainty"]
    )
    assert constraint.anchor.standard_model.source_url == sm["source_url"]
    assert constraint.anchor.joint_fit_context.rd_value == pytest.approx(
        joint["rd_value"]
    )
    assert constraint.anchor.joint_fit_context.correlation_rd_rdstar == pytest.approx(
        joint["correlation_rd_rdstar"]
    )
    assert constraint.anchor.budget_band.central_residual == pytest.approx(central)
    assert constraint.anchor.budget_band.combined_sigma == pytest.approx(combined)
    assert constraint.anchor.budget == pytest.approx(central + combined)
    assert constraint.anchor.budget == pytest.approx(0.0390830459735946)
    assert "B026.yaml" in constraint.anchor.budget_band.source


def test_anchor_loads_through_scaffold_and_fails_loudly(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b026_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b026_module, "load_anchor", spy_load_anchor)
    anchor = b026_module._load_b026_anchor(_PID)

    assert calls == [("canonical_average",), ("sm_reference_used_by_hflav",)]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)
    assert anchor.sm_value == pytest.approx(0.254)

    def missing_load_anchor(*args, **kwargs):
        raise AnchorError("forced missing B026 anchor")

    monkeypatch.setattr(b026_module, "load_anchor", missing_load_anchor)
    with pytest.raises(AnchorError, match="forced missing B026 anchor"):
        b026_module._load_b026_anchor(_PID)

    with pytest.raises(AnchorError):
        original_load_anchor(
            _PID,
            family="beauty",
            candidates=("not_a_b026_block",),
        )


@pytest.mark.parametrize(
    "target_candidates",
    [
        ("canonical_average",),
        ("sm_reference_used_by_hflav",),
    ],
)
def test_anchor_rejects_mismatched_load_anchor_block_key(monkeypatch, target_candidates):
    original_load_anchor = b026_module.load_anchor

    def mismatched_load_anchor(*args, **kwargs):
        anchor = original_load_anchor(*args, **kwargs)
        if tuple(kwargs["candidates"]) == target_candidates:
            return replace(anchor, block_key="wrong_block")
        return anchor

    monkeypatch.setattr(b026_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected 'wrong_block'"):
        b026_module._load_b026_anchor(_PID)


def test_evaluate_without_input_degrades_gracefully():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.sm_prediction == pytest.approx(constraint.anchor.sm_value)
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_sm_limit_rdstar_matches_yaml_sm_reference_and_manual_budget():
    constraint = fcc.get(_PID)
    couplings = _bc_proxy_couplings(bottom_right=0.0j)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))
    expected_ratio = abs(constraint.anchor.sm_value - constraint.anchor.value) / (
        abs(constraint.anchor.value - constraint.anchor.sm_value)
        + math.sqrt(
            constraint.anchor.experimental.uncertainty**2
            + constraint.anchor.standard_model.uncertainty**2
        )
    )

    assert constraint.sm_inputs.mode == "B->D*"
    assert result.predicted == pytest.approx(0.254)
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.predicted == pytest.approx(constraint.anchor.sm_value)
    assert result.ratio == pytest.approx(expected_ratio)
    assert result.ratio == pytest.approx(0.6908366358712635)
    assert result.passes is True
    assert result.diagnostics["sm_formula_minus_anchor"] == pytest.approx(0.0)


def test_np_prediction_matches_underlying_core_and_independent_formula():
    constraint = fcc.get(_PID)
    couplings = _bc_proxy_couplings(charm_left=1.0 + 0.2j, bottom_right=3.8 - 0.1j)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    direct = evaluate_rd_lfu_ratio(couplings, inputs=constraint.sm_inputs)
    manual, scalar_shift = _manual_rdstar_proxy(
        couplings,
        sm_rdstar=constraint.anchor.sm_value,
    )

    assert direct.mode == "B->D*"
    assert result.predicted == pytest.approx(direct.ratio)
    assert result.predicted == pytest.approx(manual)
    assert result.diagnostics["scalar_amplitude_shift"] == pytest.approx(scalar_shift)
    assert result.diagnostics["response_factor"] == pytest.approx(
        abs(1.0 + scalar_shift) ** 2
    )
    assert result.diagnostics["np_shift_rdstar"] == pytest.approx(
        result.predicted - constraint.anchor.sm_value
    )


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    couplings = _bc_proxy_couplings(charm_left=1.0 + 0.2j, bottom_right=3.8 - 0.1j)
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
        "charm_left_coupling",
        "bottom_right_coupling",
        "charm_left_overlap",
        "bottom_right_overlap",
        "scalar_overlap_proxy",
        "scalar_amplitude_shift",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "m_kk_gev",
        "matching_scale_gev",
        "bottom_mass_gev",
        "tau_mass_gev",
        "gf_gev_minus2",
        "response_factor",
        "np_shift_rdstar",
        "budget_combined_sigma",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["up_down_sector_indices"] == (1, 2)
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_bc_proxy_couplings(bottom_right=2.0), True),
        (_bc_proxy_couplings(bottom_right=8.0), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.build_from_quark_couplings(couplings))
    manual, _ = _manual_rdstar_proxy(couplings, sm_rdstar=constraint.anchor.sm_value)

    assert result.passes is expected_pass
    assert result.predicted == pytest.approx(manual)
    assert result.ratio == pytest.approx(
        abs(result.predicted - constraint.anchor.value) / constraint.anchor.budget
    )
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_optional_kk_ew_mass_extra_changes_matching_scale():
    couplings = _bc_proxy_couplings(bottom_right=8.0)
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
    assert abs(ew_result.diagnostics["scalar_amplitude_shift"]) == pytest.approx(
        abs(default_result.diagnostics["scalar_amplitude_shift"]) / 4.0
    )


def test_evaluate_is_pure_and_deterministic():
    couplings = _bc_proxy_couplings(charm_left=1.0 + 0.2j, bottom_right=3.8 - 0.1j)
    before_left_up = couplings.left_up.copy()
    before_right_down = couplings.right_down.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_up, before_left_up)
    np.testing.assert_array_equal(couplings.right_down, before_right_down)
