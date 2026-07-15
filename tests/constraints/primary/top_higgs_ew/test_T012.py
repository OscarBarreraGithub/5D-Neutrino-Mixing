"""Production tests for T012 (Z -> c cbar pole observables)."""

from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.top_higgs_ew import T012 as t012_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.zpole import (
    evaluate_quark_pseudo_observables,
    shifted_couplings,
    sm_couplings,
)
from tests.rs_ew_phase3b_helpers import sample_rs_ew_point, sm_limit_rs_ew_point

_PID = "T012"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "top_higgs_ew" / "T012.yaml"


def _yaml_entries_by_observable():
    with open(_SIDECAR) as handle:
        entries = yaml.safe_load(handle)["pdg_or_equivalent"]
    return {entry["observable"]: entry for entry in entries}


def _zcc_couplings(
    *,
    left_up_diag: tuple[float, float, float] = (0.1, 0.1, 0.1),
    right_up_diag: tuple[float, float, float] = (0.1, 0.1, 0.1),
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    left_up_overlap = np.diag(np.asarray(left_up_diag, dtype=float)).astype(
        np.complex128
    )
    right_up_overlap = np.diag(np.asarray(right_up_diag, dtype=float)).astype(
        np.complex128
    )
    zeros = np.zeros((3, 3), dtype=np.complex128)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=left_up_overlap.copy(),
        right_up_overlap=right_up_overlap,
        right_down_overlap=zeros,
        left_up=left_up_overlap.copy(),
        left_down=left_up_overlap.copy(),
        right_up=right_up_overlap.copy(),
        right_down=zeros,
    )


def _malformed_zcc_rs_ew_couplings(kind: str) -> SimpleNamespace:
    right = np.zeros((3, 3), dtype=np.complex128)
    if kind == "nan":
        left = np.zeros((3, 3), dtype=np.complex128)
        left[1, 1] = complex(np.nan, 0.0)
    elif kind == "shape":
        left = np.zeros((1, 1), dtype=np.complex128)
    else:
        raise AssertionError(f"unknown malformed coupling kind {kind!r}")
    return SimpleNamespace(
        z_delta_g_L_u=left,
        z_delta_g_R_u=right,
        matching_assumption="malformed z_delta_g_* regression fixture",
        model_label="malformed-test",
        kk_ew_mass_gev=3000.0,
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "top_higgs_ew"
    assert constraint.observable == "max Zcc pull from R_c^0 and A_c"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    entries = _yaml_entries_by_observable()
    rc = entries["R_c^0"]
    afb = entries["A_FB^{0,c}"]
    ac = entries["A_c"]

    assert constraint.anchor.r_c.value == pytest.approx(rc["value"])
    assert constraint.anchor.r_c.uncertainty == pytest.approx(rc["uncertainty"])
    assert constraint.anchor.r_c.source_url == rc["source_url"]
    assert constraint.anchor.a_fb.value == pytest.approx(afb["value"])
    assert constraint.anchor.a_fb.uncertainty == pytest.approx(afb["uncertainty"])
    assert constraint.anchor.a_c.value == pytest.approx(ac["value"])
    assert constraint.anchor.a_c.uncertainty == pytest.approx(ac["uncertainty"])

    assert constraint.anchor.budgets["R_c^0"].combined_sigma == pytest.approx(
        rc["uncertainty"]
    )
    assert constraint.anchor.budgets["A_c"].combined_sigma == pytest.approx(
        ac["uncertainty"]
    )

    with pytest.raises(AnchorError):
        t012_module._load_scaffold_list_anchor("not a T012 observable", process_id=_PID)


def test_t012_list_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = t012_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(t012_module, "load_anchor", spy_load_anchor)
    anchor = t012_module._load_t012_anchor(_PID)

    assert calls == [
        ("pdg_or_equivalent[0]",),
        ("pdg_or_equivalent[1]",),
        ("pdg_or_equivalent[2]",),
    ]
    assert anchor.r_c.value == pytest.approx(fcc.get(_PID).anchor.r_c.value)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent[99]",
            value=1.0,
            uncertainty=0.1,
        )

    monkeypatch.setattr(t012_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        t012_module._load_scaffold_list_anchor("R_c^0", process_id=_PID)


def _manual_zcc_observables(constraint, point):
    inputs = constraint.sm_inputs
    couplings = point.extras["rs_ew_couplings"]
    delta_g_left = complex(couplings.z_delta_g_L_u[1, 1])
    delta_g_right = complex(couplings.z_delta_g_R_u[1, 1])
    shifted_charm = shifted_couplings(
        sm_couplings("c", inputs),
        delta_g_left=delta_g_left,
        delta_g_right=delta_g_right,
    )
    observables = evaluate_quark_pseudo_observables(
        "c",
        {"c": shifted_charm},
        inputs=inputs,
    )
    return delta_g_left, delta_g_right, shifted_charm, observables


def _manual_t012_ratio(constraint, observables):
    r_c_pull = (
        observables.r_q - constraint.anchor.r_c.value
    ) / constraint.anchor.budgets["R_c^0"].combined_sigma
    a_c_pull = (
        observables.a_q - constraint.anchor.a_c.value
    ) / constraint.anchor.budgets["A_c"].combined_sigma
    return max(abs(r_c_pull), abs(a_c_pull))


def test_sm_zpole_numbers_match_independent_effective_coupling_recomputation():
    constraint = fcc.get(_PID)
    inputs = constraint.sm_inputs
    s2 = inputs.sin2_theta_eff
    gl_u = 0.5 - (2.0 / 3.0) * s2
    gr_u = -(2.0 / 3.0) * s2
    gl_d = -0.5 + (1.0 / 3.0) * s2
    gr_d = (1.0 / 3.0) * s2
    up_weight = 3.0 * (gl_u * gl_u + gr_u * gr_u)
    down_weight = 3.0 * (gl_d * gl_d + gr_d * gr_d)
    charm_weight = up_weight * inputs.radiator_for("c")
    manual_r_c = charm_weight / (up_weight + charm_weight + 3.0 * down_weight)
    manual_a_c = (gl_u * gl_u - gr_u * gr_u) / (gl_u * gl_u + gr_u * gr_u)

    result = constraint.evaluate(sm_limit_rs_ew_point())

    assert manual_r_c == pytest.approx(constraint.anchor.r_c.value)
    assert manual_r_c == pytest.approx(0.1721)
    assert manual_a_c == pytest.approx(constraint.sm_observables.a_q)
    assert manual_a_c == pytest.approx(0.667, abs=1.0e-3)
    assert result.diagnostics["sm_validation"]["r_c_formula"] == pytest.approx(
        manual_r_c
    )
    assert result.diagnostics["sm_validation"]["a_c_formula"] == pytest.approx(
        manual_a_c
    )
    assert result.diagnostics["a_fb_context"]["experimental"] == pytest.approx(
        constraint.anchor.a_fb.value
    )
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["delta_g_left_c"] == pytest.approx(0.0j)
    assert result.diagnostics["delta_g_right_c"] == pytest.approx(0.0j)
    assert "needs_human_physics" not in result.diagnostics


def test_rigorous_rs_ew_zcc_point_matches_independent_core_recomputation():
    constraint = fcc.get(_PID)
    point = sample_rs_ew_point()
    delta_left, delta_right, shifted_charm, manual_observables = (
        _manual_zcc_observables(constraint, point)
    )
    result = constraint.evaluate(point)

    assert delta_left != pytest.approx(0.0j)
    assert delta_right != pytest.approx(0.0j)
    assert manual_observables.r_q != pytest.approx(constraint.sm_observables.r_q)
    assert manual_observables.a_q != pytest.approx(constraint.sm_observables.a_q)
    assert result.diagnostics["evaluated"] is True
    assert "zcc_proxy" not in result.diagnostics
    assert "needs_human_physics" not in result.diagnostics
    assert result.diagnostics["delta_g_left_c"] == pytest.approx(delta_left)
    assert result.diagnostics["delta_g_right_c"] == pytest.approx(delta_right)
    assert result.diagnostics["shifted_zcc_couplings"]["g_left"] == pytest.approx(
        shifted_charm.g_left
    )
    assert result.diagnostics["shifted_zcc_couplings"]["g_right"] == pytest.approx(
        shifted_charm.g_right
    )
    assert result.diagnostics["observables"]["R_c^0"]["predicted"] == pytest.approx(
        manual_observables.r_q
    )
    assert result.diagnostics["observables"]["A_c"]["predicted"] == pytest.approx(
        manual_observables.a_q
    )
    assert result.ratio == pytest.approx(_manual_t012_ratio(constraint, manual_observables))


def test_old_style_point_reports_unevaluated_rs_ew_missing_extra():
    old_style_point = point_builder.build_from_quark_couplings(
        _zcc_couplings(left_up_diag=(0.0, 20.0, 0.0), right_up_diag=(0.0, 0.0, 0.0))
    )
    result = fcc.get(_PID).evaluate(old_style_point)

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_ew_couplings"
    assert result.diagnostics["legacy_quark_mass_basis_couplings_present"] is True
    assert "rs_ew_couplings not provided" in result.notes


@pytest.mark.parametrize(
    ("kind", "message"),
    [
        ("nan", "z_delta_g_L_u[1,1] must be finite"),
        ("shape", "z_delta_g_L_u[1,1] is not available"),
    ],
)
def test_present_malformed_z_delta_g_coupling_fails_closed(kind, message):
    point = point_builder.make_point(
        rs_ew_couplings=_malformed_zcc_rs_ew_couplings(kind)
    )

    result = fcc.get(_PID).evaluate(point)

    assert result.process_id == _PID
    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["invalid_extra"] == "rs_ew_couplings"
    assert result.diagnostics["invalid_input"] is True
    assert result.diagnostics["invalid_input_policy"] == "fail_closed_m18"
    assert result.diagnostics["exception_type"] == "ValueError"
    assert result.diagnostics["exception"] == message


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_rs_ew_diagnostics():
    point = sample_rs_ew_point()
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
    assert result.diagnostics["evaluated"] is True
    assert isinstance(result.diagnostics["delta_g_left_c"], complex)
    assert isinstance(result.diagnostics["delta_g_right_c"], complex)
    assert result.diagnostics["rs_matching_assumption"]


def test_sm_limit_rs_ew_point_recovers_committed_sm_only_output():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(sm_limit_rs_ew_point())

    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["delta_g_left_c"] == pytest.approx(0.0j)
    assert result.diagnostics["delta_g_right_c"] == pytest.approx(0.0j)
    assert result.diagnostics["observables"]["R_c^0"]["predicted"] == pytest.approx(
        constraint.sm_observables.r_q
    )
    assert result.diagnostics["observables"]["A_c"]["predicted"] == pytest.approx(
        constraint.sm_observables.a_q
    )
    assert result.predicted == pytest.approx(result.sm_prediction)


def test_evaluate_is_pure_and_deterministic():
    point = sample_rs_ew_point()
    couplings = point.extras["rs_ew_couplings"]
    before_left = couplings.z_delta_g_L_u.copy()
    before_right = couplings.z_delta_g_R_u.copy()
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.z_delta_g_L_u, before_left)
    np.testing.assert_array_equal(couplings.z_delta_g_R_u, before_right)
