"""Production tests for T010 (Z -> b bbar pole observables)."""

from __future__ import annotations

import math
from pathlib import Path
import re

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.top_higgs_ew import T010 as t010_module
from quarkConstraints.couplings import QuarkMassBasisCouplings

_PID = "T010"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "top_higgs_ew" / "T010.yaml"
_SM_SNAPSHOT = _REPO_ROOT / "flavor_catalog" / "references" / "T010" / "lepslc_2006_z_resonance.txt"
_NUMBER_RE = r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"


def _yaml_entries_by_observable():
    with open(_SIDECAR) as handle:
        entries = yaml.safe_load(handle)["pdg_or_equivalent"]
    return {entry["observable"]: entry for entry in entries}


def _snapshot_sm_fit(label: str) -> tuple[float, float, float]:
    text = _SM_SNAPSHOT.read_text()
    pattern = re.compile(
        rf"{re.escape(label)}\s*=\s*"
        rf"(?P<exp>{_NUMBER_RE})\s*\+/-\s*(?P<exp_unc>{_NUMBER_RE})"
        rf";\s*SM fit value shown there:\s*"
        rf"(?P<sm>{_NUMBER_RE})\s*\+/-\s*(?P<sm_unc>{_NUMBER_RE})"
        rf";\s*pull\s*(?P<pull>{_NUMBER_RE})"
    )
    match = pattern.search(text)
    if match is None:
        raise AssertionError(f"could not parse {label} from {_SM_SNAPSHOT}")
    return (
        float(match.group("sm")),
        float(match.group("sm_unc")),
        float(match.group("pull")),
    )


def _zbb_couplings(
    *,
    left_diag: tuple[float, float, float] = (0.1, 0.1, 0.1),
    right_diag: tuple[float, float, float] = (0.1, 0.1, 0.1),
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    left_overlap = np.diag(np.asarray(left_diag, dtype=float)).astype(np.complex128)
    right_overlap = np.diag(np.asarray(right_diag, dtype=float)).astype(np.complex128)
    zeros = np.zeros((3, 3), dtype=np.complex128)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=left_overlap,
        right_up_overlap=zeros,
        right_down_overlap=right_overlap,
        left_up=zeros,
        left_down=left_overlap.copy(),
        right_up=zeros,
        right_down=right_overlap.copy(),
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "top_higgs_ew"
    assert constraint.observable == "max Zbb pull from R_b^0 and A_b"


def test_anchor_matches_yaml_sm_snapshot_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    entries = _yaml_entries_by_observable()
    rb = entries["R_b^0"]
    afb = entries["A_FB^{0,b}"]
    ab = entries["A_b"]
    afb_pull = entries["LEP/SLC final-combination pull for A_FB^{0,b}"]
    projection = entries["FCC-ee projected relative uncertainty for R_b and A_FB^b"]
    rb_sm, rb_sm_unc, rb_pull = _snapshot_sm_fit("R_b^0")
    afb_sm, afb_sm_unc, afb_sm_pull = _snapshot_sm_fit("A_FB^(0,b)")
    ab_sm, ab_sm_unc, ab_pull = _snapshot_sm_fit("A_b")

    assert constraint.anchor.r_b.value == pytest.approx(rb["value"])
    assert constraint.anchor.r_b.uncertainty == pytest.approx(rb["uncertainty"])
    assert constraint.anchor.r_b.source_url == rb["source_url"]
    assert constraint.anchor.a_fb.value == pytest.approx(afb["value"])
    assert constraint.anchor.a_fb.uncertainty == pytest.approx(afb["uncertainty"])
    assert constraint.anchor.a_b.value == pytest.approx(ab["value"])
    assert constraint.anchor.a_b.uncertainty == pytest.approx(ab["uncertainty"])
    assert constraint.anchor.a_fb_pull.value == pytest.approx(afb_pull["value"])
    assert constraint.anchor.fcc_projection.value == pytest.approx(projection["value"])

    assert constraint.anchor.sm_fit_values["R_b^0"].value == pytest.approx(rb_sm)
    assert constraint.anchor.sm_fit_values["R_b^0"].uncertainty == pytest.approx(
        rb_sm_unc
    )
    assert constraint.anchor.sm_fit_values["R_b^0"].pull == pytest.approx(rb_pull)
    assert constraint.anchor.sm_fit_values["A_FB^{0,b}"].value == pytest.approx(
        afb_sm
    )
    assert constraint.anchor.sm_fit_values["A_FB^{0,b}"].uncertainty == pytest.approx(
        afb_sm_unc
    )
    assert constraint.anchor.sm_fit_values["A_FB^{0,b}"].pull == pytest.approx(
        afb_sm_pull
    )
    assert constraint.anchor.sm_fit_values["A_b"].value == pytest.approx(ab_sm)
    assert constraint.anchor.sm_fit_values["A_b"].uncertainty == pytest.approx(
        ab_sm_unc
    )
    assert constraint.anchor.sm_fit_values["A_b"].pull == pytest.approx(ab_pull)

    rb_budget = math.sqrt(rb["uncertainty"] ** 2 + rb_sm_unc**2)
    ab_budget = math.sqrt(ab["uncertainty"] ** 2 + ab_sm_unc**2)
    assert constraint.anchor.budgets["R_b^0"].combined_sigma == pytest.approx(
        rb_budget
    )
    assert constraint.anchor.budgets["A_b"].combined_sigma == pytest.approx(ab_budget)

    with pytest.raises(AnchorError):
        t010_module._load_scaffold_list_anchor("not a T010 observable", process_id=_PID)


def test_t010_list_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = t010_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(t010_module, "load_anchor", spy_load_anchor)
    anchor = t010_module._load_t010_anchor(_PID)

    assert calls == [
        ("pdg_or_equivalent[0]",),
        ("pdg_or_equivalent[1]",),
        ("pdg_or_equivalent[2]",),
        ("pdg_or_equivalent[3]",),
        ("pdg_or_equivalent[4]",),
    ]
    assert anchor.r_b.value == pytest.approx(fcc.get(_PID).anchor.r_b.value)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent[99]",
            value=1.0,
            uncertainty=0.1,
        )

    monkeypatch.setattr(t010_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        t010_module._load_scaffold_list_anchor("R_b^0", process_id=_PID)


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
    bottom_weight = down_weight * inputs.radiator_for("b")
    manual_r_b = bottom_weight / (2.0 * up_weight + 2.0 * down_weight + bottom_weight)
    manual_a_b = (gl_d * gl_d - gr_d * gr_d) / (gl_d * gl_d + gr_d * gr_d)

    result = constraint.evaluate(point_builder.empty_point())

    assert manual_r_b == pytest.approx(0.21562)
    assert manual_r_b == pytest.approx(constraint.sm_observables.r_q)
    assert manual_a_b == pytest.approx(constraint.sm_observables.a_q)
    assert manual_a_b == pytest.approx(0.935, abs=1.0e-3)
    assert result.diagnostics["sm_validation"]["r_b_formula"] == pytest.approx(
        manual_r_b
    )
    assert result.diagnostics["sm_validation"]["a_b_formula"] == pytest.approx(
        manual_a_b
    )
    assert result.diagnostics["a_fb_legacy_context"]["lep_slc_pull_sigma"] == (
        pytest.approx(2.8)
    )


def test_evaluate_without_input_reports_sm_limit_and_missing_extra():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.r_b.value)
    assert result.ratio < 1.0
    assert result.budget == pytest.approx(
        fcc.get(_PID).anchor.budgets["R_b^0"].combined_sigma
    )
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_proxy_diagnostics():
    couplings = _zbb_couplings(left_diag=(0.0, 0.0, 1.0), right_diag=(0.0, 0.0, 0.0))
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
    assert result.passes is False
    assert result.ratio > 1.0
    assert isinstance(result.diagnostics["delta_g_left_b"], float)
    assert isinstance(result.diagnostics["delta_g_right_b"], float)
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["rs_matching_assumption"]
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["zbb_proxy"]["matching_assumption"]


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_zbb_couplings(left_diag=(0.2, 0.2, 0.2), right_diag=(0.1, 0.1, 0.1)), True),
        (_zbb_couplings(left_diag=(0.0, 0.0, 5.0), right_diag=(0.0, 0.0, 0.0)), False),
    ],
)
def test_safe_point_passes_and_large_proxy_shift_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    point = point_builder.build_from_quark_couplings(couplings)
    result = fcc.get(_PID).evaluate(point)

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 5.0


def test_optional_kk_ew_mass_extra_changes_proxy_scale():
    couplings = _zbb_couplings(left_diag=(0.0, 0.0, 5.0), right_diag=(0.0, 0.0, 0.0))
    default_point = point_builder.build_from_quark_couplings(couplings)
    ew_point = point_builder.make_point(
        quark_mass_basis_couplings=couplings,
        kk_ew_mass_gev=6000.0,
    )
    default_result = fcc.get(_PID).evaluate(default_point)
    ew_result = fcc.get(_PID).evaluate(ew_point)

    assert default_result.diagnostics["zbb_proxy"]["m_kk_gev"] == pytest.approx(3000.0)
    assert ew_result.diagnostics["zbb_proxy"]["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert ew_result.diagnostics["delta_g_left_b"] == pytest.approx(
        default_result.diagnostics["delta_g_left_b"] / 4.0
    )


def test_evaluate_is_pure_and_deterministic():
    couplings = _zbb_couplings(left_diag=(0.0, 0.0, 1.0), right_diag=(0.0, 0.0, 0.0))
    before_left = couplings.left_overlap.copy()
    before_right = couplings.right_down_overlap.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_overlap, before_left)
    np.testing.assert_array_equal(couplings.right_down_overlap, before_right)
