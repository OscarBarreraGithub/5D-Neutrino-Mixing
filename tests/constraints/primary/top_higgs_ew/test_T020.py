"""Production tests for T020 (LFV h -> e mu)."""

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
from flavor_catalog_constraints.physics_adapters.higgs_lfv import (
    higgs_lfv_yukawa_proxy_input,
)
from flavor_catalog_constraints.primary.top_higgs_ew import T020 as t020_module
from quarkConstraints.higgs_lfv import (
    h_lfv_branching_fraction_from_yukawas as core_h_lfv_branching_fraction_from_yukawas,
)

_PID = "T020"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "top_higgs_ew" / "T020.yaml"
_NUMBER_RE = re.compile(
    r"^\s*<?\s*(?P<base>[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+))"
    r"(?:(?:[eE](?P<e_power>[+-]?[0-9]+))|"
    r"(?:\s*x\s*10\s*\^\s*(?P<x_power>[+-]?[0-9]+)))?"
    r"\s*(?P<percent>%?)\s*$"
)


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _entry(value_id: str):
    values = _yaml_pdg_block()["values"]
    return next(entry for entry in values if entry["value_id"] == value_id)


def _limit(entry, key: str = "normalized_value") -> float:
    raw = entry.get(key, entry["value"])
    match = _NUMBER_RE.match(str(raw))
    if match is None:
        raise AssertionError(f"not a numeric limit string: {raw!r}")
    value = float(match.group("base"))
    power = match.group("e_power") or match.group("x_power")
    if power is not None:
        value *= 10.0 ** int(power)
    if match.group("percent"):
        value /= 100.0
    return value


def _proxy_point(
    y_e_mu: complex,
    y_mu_e: complex = 0.0j,
):
    proxy = higgs_lfv_yukawa_proxy_input(
        "e",
        "mu",
        y_e_mu,
        y_mu_e,
        source="T020 test proxy",
    )
    return point_builder.make_point(lepton_mass_basis_couplings=proxy)


def _core_higgs_lfv_result(
    y_e_mu: complex,
    y_mu_e: complex = 0.0j,
    *,
    constraint,
):
    return core_h_lfv_branching_fraction_from_yukawas(
        initial_flavor="e",
        final_flavor="mu",
        yukawa_ij=y_e_mu,
        yukawa_ji=y_mu_e,
        br_limit=constraint.anchor.budget,
        inputs=constraint.sm_inputs,
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "top_higgs_ew"
    assert constraint.observable == "BR(h -> e mu)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    cms = _entry("CMS2023:T020:emu_limit")
    pdg_cms = _entry("PDG2025:T020:cms_run2")
    pdg_atlas = _entry("PDG2025:T020:atlas_run2")
    atlas = _entry("ATLAS2020:T020:emu_limit")
    higgs_mass = _entry("PDG2025:T020:higgs_mass_hypothesis")
    cms_dataset = _entry("CMS2023:T020:dataset")
    cms_mass_scan = _entry("CMS2023:T020:mass_scan")
    atlas_dataset = _entry("ATLAS2020:T020:dataset")
    inputs = constraint.sm_inputs
    limit = _limit(cms)
    yukawa_limit = math.sqrt(
        limit
        * inputs.total_higgs_width_gev
        * 8.0
        * math.pi
        / inputs.higgs_mass_gev
    )

    assert constraint.anchor.experimental.value == pytest.approx(limit)
    assert constraint.anchor.experimental.source_url == cms["source_url"]
    assert constraint.anchor.experimental.units == "branching fraction"
    assert constraint.anchor.cms_limit.limit == pytest.approx(limit)
    assert constraint.anchor.cms_limit.expected_limit == pytest.approx(
        _limit(cms, "expected_value")
    )
    assert constraint.anchor.pdg_cms_limit.limit == pytest.approx(_limit(pdg_cms))
    assert constraint.anchor.pdg_atlas_limit.limit == pytest.approx(_limit(pdg_atlas))
    assert constraint.anchor.atlas_limit.limit == pytest.approx(_limit(atlas))
    assert constraint.anchor.higgs_mass_hypothesis["value"] == higgs_mass["value"]
    assert constraint.anchor.cms_dataset["value"] == cms_dataset["value"]
    assert constraint.anchor.cms_mass_scan["value"] == cms_mass_scan["value"]
    assert constraint.anchor.atlas_dataset["value"] == atlas_dataset["value"]
    assert constraint.anchor.budget == pytest.approx(limit)
    assert constraint.anchor.effective_yukawa_limit == pytest.approx(yukawa_limit)
    assert constraint.anchor.effective_yukawa_limit == pytest.approx(1.8956345e-4)

    with pytest.raises(AnchorError):
        t020_module._load_scaffold_value_anchor(
            "not a T020 value_id",
            process_id=_PID,
        )


def test_t020_value_anchor_routes_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = t020_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(t020_module, "load_anchor", spy_load_anchor)
    anchor = t020_module._load_t020_anchor(_PID)

    assert calls == [("pdg_or_equivalent.values[3]",)]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent.values[99]",
            value=1.0,
            uncertainty=None,
        )

    monkeypatch.setattr(t020_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        t020_module._load_scaffold_value_anchor(
            "CMS2023:T020:emu_limit",
            process_id=_PID,
        )


def test_sm_higgs_lfv_rate_is_zero_and_formula_inputs_are_documented():
    constraint = fcc.get(_PID)
    result = constraint.evaluate(_proxy_point(0.0j, 0.0j))

    assert constraint.sm_result.branching_fraction == pytest.approx(0.0)
    assert constraint.sm_result.sm_branching_fraction == pytest.approx(0.0)
    assert result.predicted == pytest.approx(0.0)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.diagnostics["m_h_gev"] == pytest.approx(
        constraint.sm_inputs.higgs_mass_gev
    )
    assert result.diagnostics["total_higgs_width_gev"] == pytest.approx(
        constraint.sm_inputs.total_higgs_width_gev
    )
    assert result.passes is True


def test_proxy_numerics_match_core_higgs_lfv_evaluator_and_manual_formula():
    constraint = fcc.get(_PID)
    y_e_mu = 1.0e-4 + 0.4e-4j
    y_mu_e = 0.5e-4j
    result = constraint.evaluate(_proxy_point(y_e_mu, y_mu_e))
    core_result = _core_higgs_lfv_result(y_e_mu, y_mu_e, constraint=constraint)
    manual_width = (
        constraint.sm_inputs.higgs_mass_gev
        * (abs(y_e_mu) ** 2 + abs(y_mu_e) ** 2)
        / (8.0 * math.pi)
    )
    manual_yukawa_norm_squared = abs(y_e_mu) ** 2 + abs(y_mu_e) ** 2
    manual_yukawa_norm = math.sqrt(manual_yukawa_norm_squared)
    manual_br = manual_width / constraint.sm_inputs.total_higgs_width_gev

    assert result.predicted == pytest.approx(core_result.branching_fraction)
    assert result.predicted == pytest.approx(manual_br)
    assert result.ratio == pytest.approx(core_result.ratio_to_limit)
    assert result.passes is core_result.passes
    assert result.diagnostics["partial_width_gev"] == pytest.approx(manual_width)
    assert result.diagnostics["effective_yukawa_norm"] == pytest.approx(
        core_result.yukawa_norm
    )
    assert result.diagnostics["effective_yukawa_norm"] == pytest.approx(
        manual_yukawa_norm
    )
    assert result.diagnostics["effective_yukawa_norm"] < result.diagnostics[
        "effective_yukawa_limit"
    ]
    assert result.diagnostics["effective_yukawa_norm_squared"] == pytest.approx(
        core_result.yukawa_norm_squared
    )
    assert result.diagnostics["effective_yukawa_norm_squared"] == pytest.approx(
        manual_yukawa_norm_squared
    )
    assert result.diagnostics["y_e_mu"] == pytest.approx(y_e_mu)
    assert result.diagnostics["y_mu_e"] == pytest.approx(y_mu_e)
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["higgs_lfv_proxy"][
        "matching_assumption"
    ]


def test_matrix_proxy_extracts_e_mu_entries():
    matrix = np.zeros((3, 3), dtype=np.complex128)
    matrix[0, 1] = 1.1e-4 + 0.2e-4j
    matrix[1, 0] = -0.3e-4j
    point = point_builder.make_point(
        lepton_mass_basis_couplings={
            "higgs_yukawa_matrix": matrix,
            "source": "T020 matrix proxy",
        }
    )
    result = fcc.get(_PID).evaluate(point)
    constraint = fcc.get(_PID)
    core_result = _core_higgs_lfv_result(
        matrix[0, 1], matrix[1, 0], constraint=constraint
    )

    assert result.predicted == pytest.approx(core_result.branching_fraction)
    assert result.diagnostics["partial_width_gev"] == pytest.approx(
        core_result.partial_width_gev
    )
    assert result.diagnostics["effective_yukawa_norm"] == pytest.approx(
        core_result.yukawa_norm
    )
    assert result.diagnostics["effective_yukawa_norm_squared"] == pytest.approx(
        core_result.yukawa_norm_squared
    )
    assert result.diagnostics["y_e_mu"] == pytest.approx(matrix[0, 1])
    assert result.diagnostics["y_mu_e"] == pytest.approx(matrix[1, 0])


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
    result = fcc.get(_PID).evaluate(_proxy_point(0.6e-4 + 0.1e-4j, 0.1e-4j))

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
    for key in ("y_e_mu", "y_mu_e"):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "effective_yukawa_norm",
        "effective_yukawa_norm_squared",
        "effective_yukawa_limit",
        "partial_width_gev",
        "total_higgs_width_gev",
        "m_h_gev",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["evaluated"] is True


@pytest.mark.parametrize(
    ("point", "expected_pass"),
    [
        (_proxy_point(1.0e-4), True),
        (_proxy_point(2.0e-4), False),
        (_proxy_point(1.0e-3), False),
    ],
)
def test_safe_point_passes_and_excluded_points_fail(point, expected_pass: bool):
    result = fcc.get(_PID).evaluate(point)

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_e_mu_limit_is_tighter_than_existing_tau_higgs_lfv_limits():
    t020 = fcc.get(_PID)
    t018 = fcc.get("T018")
    t019 = fcc.get("T019")

    assert t020.anchor.budget < t018.anchor.budget
    assert t020.anchor.budget < t019.anchor.budget
    assert t020.anchor.effective_yukawa_limit < t018.anchor.effective_yukawa_limit
    assert t020.anchor.effective_yukawa_limit < t019.anchor.effective_yukawa_limit


def test_evaluate_is_pure_and_deterministic():
    proxy = higgs_lfv_yukawa_proxy_input(
        "e",
        "mu",
        0.10e-3 + 0.01e-3j,
        0.03e-3j,
        source="T020 test proxy",
    )
    point = point_builder.make_point(lepton_mass_basis_couplings=proxy)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.get_extra("lepton_mass_basis_couplings") == proxy
