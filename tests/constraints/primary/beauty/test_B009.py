"""Production tests for B009 (B+ -> tau+ nu_tau)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.beauty import B009 as b009_module
from tests.constraints.charged_current_phase5b_helpers import (
    charged_with_epsilon,
    sample_charged_point,
    universal_charged_point,
)

_PID = "B009"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "beauty" / "B009.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _manual_tree_br(inputs) -> float:
    tau_gev_inv = inputs.lifetime_ps * 1.0e-12 / inputs.hbar_gev_s
    phase = (1.0 - inputs.lepton_mass_gev**2 / inputs.meson_mass_gev**2) ** 2
    return float(
        inputs.gf_gev_minus2**2
        / (8.0 * math.pi)
        * inputs.meson_mass_gev
        * inputs.lepton_mass_gev**2
        * phase
        * inputs.decay_constant_gev**2
        * inputs.ckm_abs**2
        * tau_gev_inv
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "beauty"
    assert constraint.observable == "BR(B+ -> tau+ nu_tau)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["canonical_average"]
    sm = pdg["sm_prediction"]
    f_b = pdg["theory_flag_f_B"]
    v_ub = pdg["theory_ckm_abs_vub"]
    scale = 1.0e-4
    combined = math.sqrt(
        (exp["uncertainty"] * scale) ** 2 + (sm["uncertainty"] * scale) ** 2
    )
    central = abs(exp["value"] * scale - sm["value"] * scale)

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"] * scale)
    assert constraint.anchor.experimental.uncertainty == pytest.approx(
        exp["uncertainty"] * scale
    )
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.standard_model.value == pytest.approx(sm["value"] * scale)
    assert constraint.anchor.standard_model.uncertainty == pytest.approx(
        sm["uncertainty"] * scale
    )
    assert constraint.anchor.standard_model.source_url == sm["source_url"]
    assert constraint.anchor.f_b.value == pytest.approx(f_b["value"])
    assert constraint.anchor.f_b.units == "GeV"
    assert constraint.anchor.f_b.source_url == f_b["source_url"]
    assert constraint.anchor.v_ub.value == pytest.approx(v_ub["value"])
    assert constraint.anchor.v_ub.units == "dimensionless"
    assert constraint.anchor.v_ub.source_url == v_ub["source_url"]
    assert constraint.sm_inputs.decay_constant_gev == pytest.approx(f_b["value"])
    assert constraint.sm_inputs.ckm_abs == pytest.approx(v_ub["value"])
    assert "yaml_anchors" in constraint.sm_inputs.input_bundle
    assert constraint.anchor.budget_band.central_residual == pytest.approx(central)
    assert constraint.anchor.budget_band.combined_sigma == pytest.approx(combined)
    assert constraint.anchor.budget == pytest.approx(central + combined)
    assert constraint.anchor.budget == pytest.approx(4.493733520830467e-05)
    assert len(constraint.anchor.experimental_inputs) == 5

    with pytest.raises(anchors.AnchorError):
        anchors.load_anchor(
            _PID,
            family="beauty",
            candidates=("no_such_block",),
        )


def test_theory_inputs_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = b009_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(b009_module, "load_anchor", spy_load_anchor)
    anchor = b009_module._load_b009_anchor(_PID)

    assert ("theory_flag_f_B",) in calls
    assert ("theory_ckm_abs_vub",) in calls
    assert anchor.f_b.value == pytest.approx(0.1900)
    assert anchor.v_ub.value == pytest.approx(0.00368)


def test_absent_charged_current_degrades_non_vetoing():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.sm_prediction == pytest.approx(constraint.sm_result.branching_fraction)
    assert result.experimental == pytest.approx(constraint.anchor.value)
    assert result.ratio is None
    assert result.budget == pytest.approx(constraint.anchor.budget)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_charged_current"
    assert result.diagnostics["theory_inputs_yaml_backed"] is True
    assert result.diagnostics["f_b_anchor_block"] == "theory_flag_f_B"
    assert result.diagnostics["v_ub_anchor_block"] == "theory_ckm_abs_vub"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_sm_limit_branching_fraction_matches_independent_formula_and_yaml_anchor():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(universal_charged_point())
    manual = _manual_tree_br(constraint.sm_inputs)

    assert constraint.sm_result.branching_fraction == pytest.approx(manual)
    assert constraint.sm_result.branching_fraction == pytest.approx(
        8.630796690372665e-05
    )
    assert abs(constraint.sm_result.branching_fraction - constraint.anchor.sm_value) < (
        constraint.anchor.standard_model.uncertainty
    )
    assert result.sm_prediction == pytest.approx(manual)
    assert result.diagnostics["sm_anchor_branching_fraction"] == pytest.approx(
        constraint.anchor.sm_value
    )
    assert result.diagnostics["epsilon_ub_tau"] == 0.0j
    assert result.diagnostics["np_shift_branching_fraction"] == pytest.approx(0.0)
    assert result.passes is True


def test_charged_current_epsilon_matches_independent_recomputation():
    constraint = fcc.get(_PID)
    charged = charged_with_epsilon(
        universal_charged_point().extras["rs_charged_current"],
        {(0, 2, 2): 0.5 + 0.0j},
    )
    result = constraint.evaluate(point_builder.make_point(rs_charged_current=charged))
    expected = _manual_tree_br(constraint.sm_inputs) * abs(1.0 + 0.5) ** 2

    assert result.predicted == pytest.approx(expected)
    assert result.sm_prediction == pytest.approx(_manual_tree_br(constraint.sm_inputs))
    assert result.diagnostics["epsilon_ub_tau"] == pytest.approx(0.5 + 0.0j)
    assert result.diagnostics["amplitude_multiplier"] == pytest.approx(2.25)
    assert result.ratio == pytest.approx(
        abs(result.predicted - result.sm_prediction) / constraint.anchor.budget
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["theory_inputs_yaml_backed"] is True
    assert result.diagnostics["f_b_anchor_block"] == "theory_flag_f_B"
    assert result.diagnostics["v_ub_anchor_block"] == "theory_ckm_abs_vub"
    assert "mass proxy" not in result.notes


def test_evaluate_runs_end_to_end_with_real_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(sample_charged_point())

    for value in (
        result.predicted,
        result.ratio,
        result.budget,
        result.sm_prediction,
        result.experimental,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    assert isinstance(result.diagnostics["epsilon_ub_tau"], complex)
    for key in (
        "kk_ew_mass_gev",
        "m_w_gev",
        "m_wprime_gev",
        "amplitude_multiplier",
        "np_shift_branching_fraction",
        "budget_combined_sigma",
        "hard_veto_np_shift_budget",
        "meson_mass_gev",
        "decay_constant_gev",
        "ckm_abs",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])


@pytest.mark.parametrize(
    ("epsilon_ub_tau", "expected_pass"),
    [
        (0.1, True),
        (0.5, False),
    ],
)
def test_safe_point_passes_and_large_epsilon_shift_fails(
    epsilon_ub_tau: float,
    expected_pass: bool,
):
    charged = charged_with_epsilon(
        universal_charged_point().extras["rs_charged_current"],
        {(0, 2, 2): complex(epsilon_ub_tau)},
    )
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(rs_charged_current=charged)
    )

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_evaluate_is_pure_and_deterministic():
    point = sample_charged_point()
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.get_extra("rs_charged_current") is not None
