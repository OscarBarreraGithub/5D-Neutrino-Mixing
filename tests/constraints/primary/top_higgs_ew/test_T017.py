"""Production tests for T017 (LFV Z -> mu tau)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.zpole_lfv_mutau import (
    zpole_lfv_mutau_proxy_input,
)
from flavor_catalog_constraints.primary.top_higgs_ew import T017 as t017_module
from quarkConstraints.zpole import default_sm_inputs, partial_width_weight, sm_couplings

_PID = "T017"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "top_higgs_ew" / "T017.yaml"


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


def _proxy_point(
    left_mutau_overlap: complex,
    right_mutau_overlap: complex = 0.0j,
    *,
    m_kk_gev: float = 3000.0,
):
    proxy = zpole_lfv_mutau_proxy_input(
        left_mutau_overlap,
        right_mutau_overlap,
        m_kk_gev,
        source="T017 test proxy",
    )
    return point_builder.make_point(lepton_mass_basis_couplings=proxy)


def _manual_sm_width_weights():
    inputs = default_sm_inputs()
    weights = {}
    for flavor in ("u", "d", "s", "c", "b"):
        weights[flavor] = partial_width_weight(
            sm_couplings(flavor, inputs),
            radiator=inputs.radiator_for(flavor),
        )
    for flavor in ("e", "mu", "tau"):
        weights[flavor] = partial_width_weight(sm_couplings(flavor, inputs))
    for flavor in ("nu_e", "nu_mu", "nu_tau"):
        weights[flavor] = partial_width_weight(sm_couplings("nu", inputs))
    return weights


def _manual_lfv_br(
    left_mutau_overlap: complex,
    right_mutau_overlap: complex = 0.0j,
    *,
    m_kk_gev: float = 3000.0,
) -> tuple[float, complex, complex, float, float]:
    inputs = default_sm_inputs()
    scale = (inputs.m_z_gev / m_kk_gev) ** 2
    delta_left = complex(inputs.proxy_strength * scale * left_mutau_overlap)
    delta_right = complex(inputs.proxy_strength * scale * right_mutau_overlap)
    norm = abs(delta_left) ** 2 + abs(delta_right) ** 2
    sm_total = sum(_manual_sm_width_weights().values())
    lfv_weight = 2.0 * norm
    return (
        float(lfv_weight / (sm_total + lfv_weight)),
        delta_left,
        delta_right,
        float(norm),
        float(sm_total),
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "top_higgs_ew"
    assert constraint.observable == "BR(Z -> mu tau)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _entry("PDG2025:T017:zmutau_limit")
    atlas = _entry("ATLAS2021:T017:zmutau_combined_limit")
    atlas_leptonic = _entry("ATLAS2021:T017:zmutau_leptonic_tau_limit")
    atlas_dataset = _entry("ATLAS2021:T017:dataset")
    cms = _entry("CMS2025:T017:zmutau_limit")
    cms_dataset = _entry("CMS2025:T017:dataset")
    limit = _entry_limit(pdg)
    sm_total = sum(_manual_sm_width_weights().values())
    coupling_limit = math.sqrt(limit * sm_total / (2.0 * (1.0 - limit)))

    assert constraint.anchor.experimental.value == pytest.approx(limit)
    assert constraint.anchor.experimental.source_url == pdg["source_url"]
    assert constraint.anchor.experimental.units == "branching fraction"
    assert constraint.anchor.pdg_limit.limit == pytest.approx(limit)
    assert constraint.anchor.atlas_combined_limit.limit == pytest.approx(_entry_limit(atlas))
    assert constraint.anchor.atlas_leptonic_limit.limit == pytest.approx(
        _entry_limit(atlas_leptonic)
    )
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
        t017_module._load_scaffold_value_anchor(
            "not a T017 value_id",
            process_id=_PID,
        )


def test_t017_value_anchor_routes_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = t017_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(t017_module, "load_anchor", spy_load_anchor)
    anchor = t017_module._load_t017_anchor(_PID)

    assert calls == [("pdg_or_equivalent.values[0]",)]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent.values[99]",
            value=1.0,
            uncertainty=None,
        )

    monkeypatch.setattr(t017_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        t017_module._load_scaffold_value_anchor(
            "PDG2025:T017:zmutau_limit",
            process_id=_PID,
        )


def test_sm_zpole_total_weight_and_zero_lfv_rate_match_independent_recompute():
    sm_total = sum(_manual_sm_width_weights().values())
    result = fcc.get(_PID).evaluate(_proxy_point(0.0j, 0.0j))

    assert sm_total == pytest.approx(3.649563333333334)
    assert result.predicted == pytest.approx(0.0)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.diagnostics["sm_total_width_weight"] == pytest.approx(sm_total)
    assert result.diagnostics["initial_flavor"] == "mu"
    assert result.diagnostics["final_flavor"] == "tau"
    assert result.passes is True


def test_proxy_numerics_match_independent_effective_coupling_recomputation():
    constraint = fcc.get(_PID)
    left = 0.25 + 0.10j
    right = 0.05j
    result = constraint.evaluate(_proxy_point(left, right))
    expected_br, delta_left, delta_right, norm, sm_total = _manual_lfv_br(left, right)

    assert result.predicted == pytest.approx(expected_br)
    assert result.ratio == pytest.approx(expected_br / constraint.anchor.budget)
    assert result.diagnostics["delta_g_left_mutau"] == pytest.approx(delta_left)
    assert result.diagnostics["delta_g_right_mutau"] == pytest.approx(delta_right)
    assert result.diagnostics["effective_coupling_norm"] == pytest.approx(norm)
    assert result.diagnostics["sm_total_width_weight"] == pytest.approx(sm_total)
    assert result.diagnostics["charge_state_factor"] == pytest.approx(2.0)
    assert "mu+- tau-+" in result.diagnostics["branching_formula"]
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["z_lfv_mutau_proxy"][
        "matching_assumption"
    ]


def test_matrix_proxy_uses_mu_tau_offdiagonal_entry():
    left_matrix = [[0.0, 8.0, 0.30], [0.0, 0.0, 0.40], [0.0, 0.0, 0.0]]
    right_matrix = [[0.0, 4.0, 0.02j], [0.0, 0.0, 0.03j], [0.0, 0.0, 0.0]]
    point = point_builder.make_point(
        lepton_mass_basis_couplings={
            "left_charged_lepton_overlap": left_matrix,
            "right_charged_lepton_overlap": right_matrix,
            "m_kk_gev": 3000.0,
            "source": "T017 matrix proxy",
        }
    )
    result = fcc.get(_PID).evaluate(point)
    expected_br, delta_left, delta_right, norm, _ = _manual_lfv_br(0.40, 0.03j)

    assert result.predicted == pytest.approx(expected_br)
    assert result.diagnostics["delta_g_left_mutau"] == pytest.approx(delta_left)
    assert result.diagnostics["delta_g_right_mutau"] == pytest.approx(delta_right)
    assert result.diagnostics["effective_coupling_norm"] == pytest.approx(norm)
    assert result.diagnostics["z_lfv_mutau_proxy"]["left_mutau_overlap"] == pytest.approx(
        0.40
    )


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "lepton_mass_basis_couplings"
    assert "non-vetoing only" in result.diagnostics["passes_semantics"]


def test_invalid_lepton_input_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(lepton_mass_basis_couplings={"left_emu": 1.0})
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["invalid_extra"] == "lepton_mass_basis_couplings"
    assert result.diagnostics["exception_type"] == "KeyError"


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(_proxy_point(0.03 + 0.02j, 0.01j))

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
    for key in ("delta_g_left_mutau", "delta_g_right_mutau"):
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


@pytest.mark.parametrize(
    ("point", "expected_pass"),
    [
        (_proxy_point(0.10), True),
        (_proxy_point(20.0), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(point, expected_pass: bool):
    result = fcc.get(_PID).evaluate(point)

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_changes_proxy_scale():
    proxy = zpole_lfv_mutau_proxy_input(
        0.20,
        0.0j,
        3000.0,
        source="T017 test proxy",
    )
    default_point = point_builder.make_point(lepton_mass_basis_couplings=proxy)
    heavy_point = point_builder.make_point(
        lepton_mass_basis_couplings=proxy,
        kk_ew_mass_gev=6000.0,
    )

    default_result = fcc.get(_PID).evaluate(default_point)
    heavy_result = fcc.get(_PID).evaluate(heavy_point)

    assert default_result.diagnostics["z_lfv_mutau_proxy"]["m_kk_gev"] == pytest.approx(
        3000.0
    )
    assert heavy_result.diagnostics["z_lfv_mutau_proxy"]["m_kk_gev"] == pytest.approx(
        6000.0
    )
    assert heavy_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert heavy_result.diagnostics["delta_g_left_mutau"] == pytest.approx(
        default_result.diagnostics["delta_g_left_mutau"] / 4.0
    )
    assert heavy_result.diagnostics["effective_coupling_norm"] == pytest.approx(
        default_result.diagnostics["effective_coupling_norm"] / 16.0
    )


def test_evaluate_is_pure_and_deterministic():
    proxy = zpole_lfv_mutau_proxy_input(
        0.15 + 0.01j,
        0.03j,
        3000.0,
        source="T017 test proxy",
    )
    point = point_builder.make_point(lepton_mass_basis_couplings=proxy)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.get_extra("lepton_mass_basis_couplings") == proxy
