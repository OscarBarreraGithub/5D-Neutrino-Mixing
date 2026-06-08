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
from quarkConstraints.zpole import (
    evaluate_quark_pseudo_observables,
    shifted_couplings,
    sm_couplings,
)
from quarkConstraints.couplings import QuarkMassBasisCouplings
from tests.rs_ew_phase3b_helpers import sample_rs_ew_point, sm_limit_rs_ew_point

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


def _manual_zbb_observables(constraint, point):
    inputs = constraint.sm_inputs
    couplings = point.extras["rs_ew_couplings"]
    delta_g_left = complex(couplings.z_delta_g_L_d[2, 2])
    delta_g_right = complex(couplings.z_delta_g_R_d[2, 2])
    shifted_bottom = shifted_couplings(
        sm_couplings("b", inputs),
        delta_g_left=delta_g_left,
        delta_g_right=delta_g_right,
    )
    observables = evaluate_quark_pseudo_observables(
        "b",
        {"b": shifted_bottom},
        inputs=inputs,
    )
    return delta_g_left, delta_g_right, shifted_bottom, observables


def _manual_t010_ratio(constraint, observables):
    r_b_pull = (
        observables.r_q - constraint.anchor.r_b.value
    ) / constraint.anchor.budgets["R_b^0"].combined_sigma
    a_b_pull = (
        observables.a_q - constraint.anchor.a_b.value
    ) / constraint.anchor.budgets["A_b"].combined_sigma
    return max(abs(r_b_pull), abs(a_b_pull))


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

    result = constraint.evaluate(sm_limit_rs_ew_point())

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
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["delta_g_left_b"] == pytest.approx(0.0j)
    assert result.diagnostics["delta_g_right_b"] == pytest.approx(0.0j)
    assert result.predicted == pytest.approx(result.sm_prediction)


def test_old_style_point_reports_unevaluated_rs_ew_missing_extra():
    old_style_point = point_builder.build_from_quark_couplings(
        _zbb_couplings(left_diag=(0.0, 0.0, 5.0), right_diag=(0.0, 0.0, 0.0))
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


def test_rigorous_rs_ew_zbb_path_matches_independent_core_recomputation():
    point = sample_rs_ew_point()
    constraint = fcc.get(_PID)
    delta_left, delta_right, shifted_bottom, manual_observables = (
        _manual_zbb_observables(constraint, point)
    )
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
    assert "zbb_proxy" not in result.diagnostics
    assert result.diagnostics["delta_g_left_b"] == pytest.approx(delta_left)
    assert result.diagnostics["delta_g_right_b"] == pytest.approx(delta_right)
    assert result.diagnostics["shifted_zbb_couplings"]["g_left"] == pytest.approx(
        shifted_bottom.g_left
    )
    assert result.diagnostics["shifted_zbb_couplings"]["g_right"] == pytest.approx(
        shifted_bottom.g_right
    )
    assert result.diagnostics["observables"]["R_b^0"]["predicted"] == pytest.approx(
        manual_observables.r_q
    )
    assert result.diagnostics["observables"]["A_b"]["predicted"] == pytest.approx(
        manual_observables.a_q
    )
    assert result.ratio == pytest.approx(_manual_t010_ratio(constraint, manual_observables))
    assert result.diagnostics["minimal_rs_tree_complete"] is False
    assert result.diagnostics["minimal_rs_tree_veto_ready"] is False
    assert result.diagnostics["fermion_kk_mixing_included"] is False
    assert result.diagnostics["custodial_variant_deferred"] is True
    assert result.diagnostics["custodial_toppartner_zbL_deferred"] is True
    assert result.diagnostics["brane_kinetic_terms_deferred"] is True
    assert "minimal non-custodial Zbb is incomplete" in result.diagnostics["needs_human_physics"]
    assert "PARTIAL/NEEDS-HUMAN-PHYSICS" not in result.diagnostics["needs_human_physics"]


def test_sm_limit_rs_ew_point_recovers_committed_sm_only_output():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(sm_limit_rs_ew_point())

    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["delta_g_left_b"] == pytest.approx(0.0j)
    assert result.diagnostics["delta_g_right_b"] == pytest.approx(0.0j)
    assert result.diagnostics["observables"]["R_b^0"]["predicted"] == pytest.approx(
        constraint.sm_observables.r_q
    )
    assert result.diagnostics["observables"]["A_b"]["predicted"] == pytest.approx(
        constraint.sm_observables.a_q
    )
    assert result.predicted == pytest.approx(result.sm_prediction)


def test_evaluate_is_pure_and_deterministic():
    point = sample_rs_ew_point()
    couplings = point.extras["rs_ew_couplings"]
    before_left = couplings.z_delta_g_L_d.copy()
    before_right = couplings.z_delta_g_R_d.copy()
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.z_delta_g_L_d, before_left)
    np.testing.assert_array_equal(couplings.z_delta_g_R_d, before_right)
