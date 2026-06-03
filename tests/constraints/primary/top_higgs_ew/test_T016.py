"""Production tests for T016 (LFV Z -> e tau)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.top_higgs_ew import T016 as t016_module
from tests.constraints.primary.top_higgs_ew.z_lfv_rewire_helpers import (
    LFV_CHANNELS,
    amplified_tree_point,
    channel_deltas,
    diagonal_rs_ew_point,
    lfv_live_rs_ew_point,
    manual_lfv_br,
    manual_sm_width_weights,
    old_style_lepton_only_point,
)

_PID = "T016"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "top_higgs_ew" / "T016.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _entry(value_id: str):
    values = _yaml_pdg_block()["values"]
    return next(entry for entry in values if entry["value_id"] == value_id)


def _limit(value) -> float:
    return float(str(value).replace("<", "").strip())


def _entry_limit(entry) -> float:
    return _limit(entry["normalized_value"])


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "top_higgs_ew"
    assert constraint.observable == LFV_CHANNELS[_PID]["observable"]


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _entry("PDG2025:T016:zetau_limit")
    atlas = _entry("ATLAS2021:T016:zetau_combined_limit")
    atlas_dataset = _entry("ATLAS2021:T016:dataset")
    cms = _entry("CMS2025:T016:zetau_limit")
    cms_dataset = _entry("CMS2025:T016:dataset")
    limit = _entry_limit(pdg)
    sm_total = sum(manual_sm_width_weights().values())
    coupling_limit = math.sqrt(limit * sm_total / (2.0 * (1.0 - limit)))

    assert constraint.anchor.experimental.value == pytest.approx(limit)
    assert constraint.anchor.experimental.source_url == pdg["source_url"]
    assert constraint.anchor.experimental.units == "branching fraction"
    assert constraint.anchor.pdg_limit.limit == pytest.approx(limit)
    assert constraint.anchor.atlas_combined_limit.limit == pytest.approx(_entry_limit(atlas))
    assert constraint.anchor.cms_limit.limit == pytest.approx(_entry_limit(cms))
    assert constraint.anchor.cms_limit.expected_limit == pytest.approx(
        _limit(cms["expected_value"])
    )
    assert constraint.anchor.atlas_dataset["value"] == atlas_dataset["value"]
    assert constraint.anchor.cms_dataset["value"] == cms_dataset["value"]
    assert constraint.anchor.budget == pytest.approx(limit)
    assert constraint.anchor.sm_total_width_weight == pytest.approx(sm_total)
    assert constraint.anchor.effective_coupling_limit == pytest.approx(coupling_limit)

    with pytest.raises(AnchorError):
        t016_module._load_scaffold_value_anchor(
            "not a T016 value_id",
            process_id=_PID,
        )


def test_t016_value_anchor_routes_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = t016_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(t016_module, "load_anchor", spy_load_anchor)
    anchor = t016_module._load_t016_anchor(_PID)

    assert calls == [("pdg_or_equivalent.values[0]",)]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent.values[99]",
            value=1.0,
            uncertainty=None,
        )

    monkeypatch.setattr(t016_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        t016_module._load_scaffold_value_anchor(
            "PDG2025:T016:zetau_limit",
            process_id=_PID,
        )


def test_diagonal_rs_ew_builder_gives_zero_tree_lfv_and_loop_deferred():
    point = diagonal_rs_ew_point()
    delta_left, delta_right = channel_deltas(point.extras["rs_ew_couplings"], _PID)
    result = fcc.get(_PID).evaluate(point)

    assert delta_left == pytest.approx(0.0j, abs=1.0e-20)
    assert delta_right == pytest.approx(0.0j, abs=1.0e-20)
    assert result.predicted == pytest.approx(0.0)
    assert result.ratio == pytest.approx(0.0)
    assert result.passes is True
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["initial_flavor"] == "e"
    assert result.diagnostics["final_flavor"] == "tau"
    assert result.diagnostics["tree_level_matching_status"] == (
        "rigorous_tree_light_z_from_rs_ew_couplings"
    )
    assert result.diagnostics["loop_lfv_status"] == (
        "deferred_to_phase_7_lfv_dipole_spurion"
    )
    assert "lfv_dipole_spurion" in result.diagnostics["needs_human_physics"]
    assert "proxy" not in result.diagnostics["needs_human_physics"].lower()


def test_lfv_live_toy_nonzero_matches_independent_recompute_and_mkk_scaling():
    constraint = fcc.get(_PID)
    cfg = LFV_CHANNELS[_PID]
    low = constraint.evaluate(lfv_live_rs_ew_point(3000.0))
    high = constraint.evaluate(lfv_live_rs_ew_point(6000.0))
    delta_left, delta_right = channel_deltas(
        lfv_live_rs_ew_point(3000.0).extras["rs_ew_couplings"],
        _PID,
    )
    expected_br, norm, sm_total = manual_lfv_br(delta_left, delta_right)

    assert abs(delta_left) > 0.0 or abs(delta_right) > 0.0
    assert low.predicted == pytest.approx(expected_br)
    assert low.ratio == pytest.approx(expected_br / constraint.anchor.budget)
    assert low.diagnostics[cfg["left_key"]] == pytest.approx(delta_left)
    assert low.diagnostics[cfg["right_key"]] == pytest.approx(delta_right)
    assert low.diagnostics["effective_coupling_norm"] == pytest.approx(norm)
    assert low.diagnostics["sm_total_width_weight"] == pytest.approx(sm_total)
    assert low.diagnostics["initial_flavor"] == cfg["initial_flavor"]
    assert low.diagnostics["final_flavor"] == cfg["final_flavor"]
    assert low.diagnostics["z_delta_g_indices"] == {
        "left": "z_delta_g_L_e[0,2]",
        "right": "z_delta_g_R_e[0,2]",
    }

    assert high.diagnostics["effective_coupling_norm"] / norm == pytest.approx(
        (3000.0 / 6000.0) ** 4,
        rel=1.0e-12,
    )
    assert high.predicted / low.predicted == pytest.approx(
        (3000.0 / 6000.0) ** 4,
        rel=1.0e-8,
    )


def test_absent_rs_ew_couplings_degrades_gracefully_for_empty_and_old_style_points():
    for point, legacy_present in (
        (point_builder.empty_point(), False),
        (old_style_lepton_only_point(), True),
    ):
        result = fcc.get(_PID).evaluate(point)

        assert result.process_id == _PID
        assert result.passes is True
        assert result.predicted is None
        assert result.ratio is None
        assert result.sm_prediction == pytest.approx(0.0)
        assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
        assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
        assert result.notes.startswith("NOT EVALUATED -")
        assert result.diagnostics["evaluated"] is False
        assert result.diagnostics["missing_extra"] == "rs_ew_couplings"
        assert result.diagnostics["legacy_lepton_mass_basis_couplings_present"] is legacy_present
        assert "non-vetoing only" in result.diagnostics["passes_semantics"]


def test_invalid_rs_ew_couplings_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(point_builder.make_point(rs_ew_couplings=object()))

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["invalid_extra"] == "rs_ew_couplings"
    assert result.diagnostics["exception_type"] == "ValueError"


def test_amplified_tree_coupling_fails_limit_in_new_path():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(amplified_tree_point(_PID))
    expected_br, norm, _ = manual_lfv_br(0.01, 0.0j)

    assert result.predicted == pytest.approx(expected_br)
    assert result.diagnostics["effective_coupling_norm"] == pytest.approx(norm)
    assert result.passes is False
    assert result.ratio > 1.0


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(lfv_live_rs_ew_point())

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
    for key in ("delta_g_left_etau", "delta_g_right_etau"):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "effective_coupling_norm",
        "effective_coupling_limit",
        "lfv_width_weight",
        "sm_total_width_weight",
        "total_width_weight",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["evaluated"] is True


def test_evaluate_is_pure_and_deterministic():
    point = lfv_live_rs_ew_point()
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.get_extra("rs_ew_couplings") is not None
