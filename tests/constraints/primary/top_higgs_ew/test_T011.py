"""Production tests for T011 (Z -> b bbar pole asymmetries)."""

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
from flavor_catalog_constraints.primary.top_higgs_ew import T011 as t011_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.zpole import (
    evaluate_quark_pseudo_observables,
    shifted_couplings,
    sm_couplings,
)
from tests.rs_ew_phase3b_helpers import sample_rs_ew_point, sm_limit_rs_ew_point

_PID = "T011"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "top_higgs_ew" / "T011.yaml"
_NUMBER_RE = r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"


def _sidecar_yaml():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)


def _yaml_entries_by_observable():
    entries = _sidecar_yaml()["pdg_or_equivalent"]
    if isinstance(entries, dict):
        return {
            entry.get("observable", key): entry
            for key, entry in entries.items()
            if isinstance(entry, dict)
        }
    return {entry["observable"]: entry for entry in entries}


def _snapshot_sm_fit(label: str) -> tuple[float, float, float]:
    path = (
        _REPO_ROOT
        / "flavor_catalog"
        / "references"
        / "T011"
        / "lepslc_2006_z_resonance_bottom.txt"
    )
    text = path.read_text()
    pattern = re.compile(
        rf"{re.escape(label)}\s*=\s*"
        rf"(?P<exp>{_NUMBER_RE})\s*\+/-\s*(?P<exp_unc>{_NUMBER_RE})"
        rf";\s*SM fit value shown there:\s*"
        rf"(?P<sm>{_NUMBER_RE})\s*\+/-\s*(?P<sm_unc>{_NUMBER_RE})"
        rf";\s*pull\s*(?P<pull>{_NUMBER_RE})"
    )
    match = pattern.search(text)
    if match is None:
        raise AssertionError(f"could not parse {label} from {path}")
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
    assert constraint.observable == "max Zbb asymmetry NP shift from A_FB^0,b and A_b"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    data = _sidecar_yaml()
    entries = _yaml_entries_by_observable()
    afb = entries["A_FB^{0,b}"]
    ab = entries["A_b"]
    afb_sm, afb_sm_unc, afb_pull = _snapshot_sm_fit("A_FB^(0,b)")
    ab_sm, ab_sm_unc, ab_pull = _snapshot_sm_fit("A_b")

    assert "canonical_home" not in data
    assert "merged_into" not in data
    assert constraint.anchor.source_process_id == _PID
    assert constraint.anchor.canonical_home_fallback is False
    assert constraint.anchor.a_fb.process_id == _PID
    assert constraint.anchor.a_b.process_id == _PID
    assert constraint.anchor.a_fb.value == pytest.approx(afb["value"])
    assert constraint.anchor.a_fb.uncertainty == pytest.approx(afb["uncertainty"])
    assert constraint.anchor.a_fb.source_url == afb["source_url"]
    assert constraint.anchor.a_b.value == pytest.approx(ab["value"])
    assert constraint.anchor.a_b.uncertainty == pytest.approx(ab["uncertainty"])
    assert constraint.anchor.a_b.source_url == ab["source_url"]

    assert constraint.anchor.sm_fit_values["A_FB^{0,b}"].value == pytest.approx(
        afb_sm
    )
    assert constraint.anchor.sm_fit_values["A_FB^{0,b}"].uncertainty == pytest.approx(
        afb_sm_unc
    )
    assert constraint.anchor.sm_fit_values["A_FB^{0,b}"].pull == pytest.approx(
        afb_pull
    )
    assert constraint.anchor.sm_fit_values["A_b"].value == pytest.approx(ab_sm)
    assert constraint.anchor.sm_fit_values["A_b"].uncertainty == pytest.approx(
        ab_sm_unc
    )
    assert constraint.anchor.sm_fit_values["A_b"].pull == pytest.approx(ab_pull)

    afb_budget = abs(afb["value"] - afb_sm) + math.sqrt(
        afb["uncertainty"] ** 2 + afb_sm_unc**2
    )
    ab_budget = abs(ab["value"] - ab_sm) + math.sqrt(
        ab["uncertainty"] ** 2 + ab_sm_unc**2
    )
    assert constraint.anchor.budgets["A_FB^{0,b}"].hard_veto_budget == pytest.approx(
        afb_budget
    )
    assert constraint.anchor.budgets["A_b"].hard_veto_budget == pytest.approx(
        ab_budget
    )

    with pytest.raises(AnchorError):
        t011_module._load_scaffold_entry_anchor(
            "not a T011 observable",
            process_id=_PID,
        )


def test_t010_canonical_home_fallback_is_rejected(monkeypatch):
    legacy_stub = {
        "process_id": _PID,
        "family": "top_higgs_ew",
        "canonical_home": {"process_id": "T010"},
    }

    monkeypatch.setattr(
        t011_module,
        "load_full_yaml",
        lambda *args, **kwargs: legacy_stub,
    )

    with pytest.raises(AnchorError, match="T010 fallback"):
        t011_module._load_t011_anchor(_PID)


def test_t011_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, str, tuple[str, ...]]] = []
    original_load_anchor = t011_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append((args[0], kwargs["family"], tuple(kwargs["candidates"])))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(t011_module, "load_anchor", spy_load_anchor)
    anchor = t011_module._load_t011_anchor(_PID)

    assert calls == [
        (_PID, "top_higgs_ew", ("A_FB^{0,b}",)),
        (_PID, "top_higgs_ew", ("A_b",)),
    ]
    assert anchor.a_b.value == pytest.approx(fcc.get(_PID).anchor.a_b.value)
    assert anchor.source_process_id == _PID
    assert anchor.canonical_home_fallback is False

    def t010_load_anchor(*args, **kwargs):
        return Anchor(
            process_id="T010",
            block_key=kwargs["candidates"][0],
            value=1.0,
            uncertainty=0.1,
        )

    monkeypatch.setattr(t011_module, "load_anchor", t010_load_anchor)
    with pytest.raises(AnchorError, match="standalone T011 process"):
        t011_module._load_scaffold_entry_anchor("A_b", process_id=_PID)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent[99]",
            value=1.0,
            uncertainty=0.1,
        )

    monkeypatch.setattr(t011_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        t011_module._load_scaffold_entry_anchor("A_b", process_id=_PID)


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


def _manual_t011_ratio(constraint, observables):
    manual_a_fb_shift = observables.a_fb - constraint.sm_observables.a_fb
    manual_a_b_shift = observables.a_q - constraint.sm_observables.a_q
    return max(
        abs(manual_a_fb_shift)
        / constraint.anchor.budgets["A_FB^{0,b}"].hard_veto_budget,
        abs(manual_a_b_shift) / constraint.anchor.budgets["A_b"].hard_veto_budget,
    )


def test_sm_zpole_asymmetries_match_independent_effective_coupling_recomputation():
    constraint = fcc.get(_PID)
    inputs = constraint.sm_inputs
    s2 = inputs.sin2_theta_eff
    gl_e = -0.5 + s2
    gr_e = s2
    gl_b = -0.5 + (1.0 / 3.0) * s2
    gr_b = (1.0 / 3.0) * s2
    manual_a_e = (gl_e * gl_e - gr_e * gr_e) / (gl_e * gl_e + gr_e * gr_e)
    manual_a_b = (gl_b * gl_b - gr_b * gr_b) / (gl_b * gl_b + gr_b * gr_b)
    manual_a_fb = 0.75 * manual_a_e * manual_a_b

    result = constraint.evaluate(sm_limit_rs_ew_point())

    assert manual_a_b == pytest.approx(constraint.sm_observables.a_q)
    assert manual_a_fb == pytest.approx(constraint.sm_observables.a_fb)
    assert manual_a_b == pytest.approx(0.935, abs=1.0e-3)
    assert manual_a_fb == pytest.approx(0.103, abs=1.0e-3)
    assert result.diagnostics["sm_validation"]["a_b_formula"] == pytest.approx(
        manual_a_b
    )
    assert result.diagnostics["sm_validation"]["a_fb_formula"] == pytest.approx(
        manual_a_fb
    )
    assert result.diagnostics["legacy_a_fb_context"]["lep_slc_pull_sigma"] == (
        pytest.approx(2.8)
    )
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["delta_g_left_b"] == pytest.approx(0.0j)
    assert result.diagnostics["delta_g_right_b"] == pytest.approx(0.0j)


def test_rigorous_rs_ew_zbb_point_matches_independent_core_recomputation():
    constraint = fcc.get(_PID)
    point = sample_rs_ew_point()
    delta_left, delta_right, shifted_bottom, manual_observables = (
        _manual_zbb_observables(constraint, point)
    )
    result = constraint.evaluate(point)

    assert delta_left != pytest.approx(0.0j)
    assert delta_right != pytest.approx(0.0j)
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
    assert result.diagnostics["observables"]["A_FB^{0,b}"]["predicted"] == pytest.approx(
        manual_observables.a_fb
    )
    assert result.diagnostics["observables"]["A_b"]["predicted"] == pytest.approx(
        manual_observables.a_q
    )
    assert result.ratio == pytest.approx(_manual_t011_ratio(constraint, manual_observables))
    assert result.diagnostics["minimal_rs_tree_complete"] is False
    assert result.diagnostics["fermion_kk_mixing_included"] is False
    assert result.diagnostics["custodial_variant_needs_human"] is True
    assert result.diagnostics["custodial_toppartner_zbL_needs_human"] is True
    assert "top-partner" in result.diagnostics["needs_human_physics"]
    assert "PARTIAL/NEEDS-HUMAN-PHYSICS" not in result.diagnostics["needs_human_physics"]


def test_old_style_point_reports_unevaluated_rs_ew_missing_extra():
    old_style_point = point_builder.build_from_quark_couplings(
        _zbb_couplings(left_diag=(0.0, 0.0, 0.0), right_diag=(0.0, 0.0, 30.0))
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
    assert isinstance(result.diagnostics["delta_g_left_b"], complex)
    assert isinstance(result.diagnostics["delta_g_right_b"], complex)
    assert result.diagnostics["rs_matching_assumption"]


def test_sm_limit_rs_ew_point_recovers_committed_sm_only_output():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(sm_limit_rs_ew_point())

    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["delta_g_left_b"] == pytest.approx(0.0j)
    assert result.diagnostics["delta_g_right_b"] == pytest.approx(0.0j)
    assert result.diagnostics["observables"]["A_FB^{0,b}"]["predicted"] == pytest.approx(
        constraint.sm_observables.a_fb
    )
    assert result.diagnostics["observables"]["A_b"]["predicted"] == pytest.approx(
        constraint.sm_observables.a_q
    )
    assert result.predicted == pytest.approx(result.sm_prediction)
    assert result.ratio == pytest.approx(0.0)


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
