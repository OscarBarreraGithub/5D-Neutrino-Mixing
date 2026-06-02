"""Production tests for T014 (down-sector FCNC Z decays)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import (
    ConstraintLevel,
    ConstraintProtocol,
    Severity,
)
from flavor_catalog_constraints.secondary.top_higgs_ew import T014 as t014_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.zpole import (
    default_sm_inputs,
    down_fcnc_branching_fraction_from_couplings,
    down_fcnc_effective_coupling_limit,
    down_fcnc_sm_hadronic_width_weight,
    down_fcnc_sm_total_width_weight,
    partial_width_weight,
    sm_couplings,
)

_PID = "T014"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = (
    _REPO_ROOT
    / "flavor_catalog"
    / "processes"
    / "secondary"
    / "top_higgs_ew"
    / "T014.yaml"
)
_INDEX = {"d": 0, "s": 1, "b": 2}


def _yaml_values():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _entry(value_id: str):
    return next(entry for entry in _yaml_values() if entry["value_id"] == value_id)


def _limit(entry) -> float:
    return float(str(entry["normalized_value"]).replace("<", "").strip())


def _set_pair(matrix, flavor_i: str, flavor_j: str, value: complex) -> None:
    i = _INDEX[flavor_i]
    j = _INDEX[flavor_j]
    matrix[i, j] = value
    matrix[j, i] = np.conj(value)


def _fcnc_couplings(
    *,
    bs_left: complex = 0.0j,
    bs_right: complex = 0.0j,
    bd_left: complex = 0.0j,
    bd_right: complex = 0.0j,
    sd_left: complex = 0.0j,
    sd_right: complex = 0.0j,
    M_KK: float = 3000.0,
) -> QuarkMassBasisCouplings:
    left_overlap = np.zeros((3, 3), dtype=np.complex128)
    right_overlap = np.zeros((3, 3), dtype=np.complex128)
    _set_pair(left_overlap, "b", "s", bs_left)
    _set_pair(right_overlap, "b", "s", bs_right)
    _set_pair(left_overlap, "b", "d", bd_left)
    _set_pair(right_overlap, "b", "d", bd_right)
    _set_pair(left_overlap, "s", "d", sd_left)
    _set_pair(right_overlap, "s", "d", sd_right)
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


def _manual_sm_hadronic_width_weights():
    inputs = default_sm_inputs()
    return {
        flavor: partial_width_weight(
            sm_couplings(flavor, inputs),
            radiator=inputs.radiator_for(flavor),
        )
        for flavor in ("u", "d", "s", "c", "b")
    }


def _manual_sm_total_width_weights():
    inputs = default_sm_inputs()
    weights = dict(_manual_sm_hadronic_width_weights())
    for flavor in ("e", "mu", "tau"):
        weights[flavor] = partial_width_weight(sm_couplings(flavor, inputs))
    for flavor in ("nu_e", "nu_mu", "nu_tau"):
        weights[flavor] = partial_width_weight(sm_couplings("nu", inputs))
    return weights


def _manual_fcnc_br(
    flavor_i: str,
    flavor_j: str,
    left_overlap: complex,
    right_overlap: complex = 0.0j,
    *,
    m_kk_gev: float = 3000.0,
) -> tuple[float, complex, complex, float, float, float, float, float]:
    inputs = default_sm_inputs()
    scale = (inputs.m_z_gev / m_kk_gev) ** 2
    delta_left = complex(inputs.proxy_strength * scale * left_overlap)
    delta_right = complex(inputs.proxy_strength * scale * right_overlap)
    norm = abs(delta_left) ** 2 + abs(delta_right) ** 2
    sm_hadronic = sum(_manual_sm_hadronic_width_weights().values())
    sm_total = sum(_manual_sm_total_width_weights().values())
    radiator = 0.5 * (
        inputs.radiator_for(flavor_i) + inputs.radiator_for(flavor_j)
    )
    fcnc_weight = 2.0 * 3.0 * radiator * norm
    return (
        float(fcnc_weight / (sm_total + fcnc_weight)),
        delta_left,
        delta_right,
        float(norm),
        float(sm_hadronic),
        float(sm_total),
        float(fcnc_weight),
        float(sm_hadronic / sm_total),
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.level is ConstraintLevel.SECONDARY
    assert constraint.family == "top_higgs_ew"
    assert constraint.observable == "max BR(Z -> bs, bd, sd) down-sector FCNC"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    value_ids = {
        "bs": "ECFA2025:T014:zbs_direct_hadronic_width_limit",
        "bd": "ECFA2025:T014:zbd_direct_hadronic_width_limit",
        "sd": "AbuAjamieh2026:T014:zsd_direct_hadronic_width_limit",
    }
    sm_hadronic = sum(_manual_sm_hadronic_width_weights().values())
    sm_total = sum(_manual_sm_total_width_weights().values())

    assert set(constraint.anchor.channels) == {"bs", "bd", "sd"}
    for channel, value_id in value_ids.items():
        entry = _entry(value_id)
        anchor = constraint.anchor.channels[channel]
        limit = _limit(entry)
        manual_coupling_limit = math.sqrt(
            limit * sm_total / (2.0 * 3.0 * (1.0 - limit))
        )
        core_coupling_limit = down_fcnc_effective_coupling_limit(
            limit,
            flavor_i=anchor.flavor_i,
            flavor_j=anchor.flavor_j,
        )

        assert anchor.value_id == value_id
        assert anchor.value == pytest.approx(limit)
        assert anchor.confidence_level == "95% CL"
        assert anchor.units == "dimensionless branching fraction"
        assert anchor.limit_type == entry["limit_type"]
        assert anchor.source_url == entry["source_url"]
        assert anchor.effective_coupling_limit == pytest.approx(
            manual_coupling_limit
        )
        assert anchor.effective_coupling_limit == pytest.approx(core_coupling_limit)

    assert constraint.anchor.value == pytest.approx(
        _limit(_entry(value_ids["bs"]))
    )
    assert constraint.anchor.sm_hadronic_width_weight == pytest.approx(sm_hadronic)
    assert constraint.anchor.sm_total_width_weight == pytest.approx(sm_total)
    assert constraint.anchor.sm_hadronic_to_total_width_ratio == pytest.approx(
        sm_hadronic / sm_total
    )
    assert constraint.anchor.sm_hadronic_width_weight == pytest.approx(
        sum(down_fcnc_sm_hadronic_width_weight().values())
    )
    assert constraint.anchor.sm_total_width_weight == pytest.approx(
        sum(down_fcnc_sm_total_width_weight().values())
    )

    with pytest.raises(AnchorError):
        t014_module._load_scaffold_limit_anchor("not a T014 value_id", process_id=_PID)


def test_t014_value_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = t014_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(t014_module, "load_anchor", spy_load_anchor)
    anchor = t014_module._load_t014_anchor(_PID)

    assert calls == [
        ("pdg_or_equivalent[0]",),
        ("pdg_or_equivalent[1]",),
        ("pdg_or_equivalent[2]",),
    ]
    assert anchor.value == pytest.approx(fcc.get(_PID).anchor.value)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent[99]",
            value=1.0,
            uncertainty=None,
        )

    monkeypatch.setattr(t014_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        t014_module._load_scaffold_limit_anchor(
            "ECFA2025:T014:zbs_direct_hadronic_width_limit",
            process_id=_PID,
        )


def test_zero_fcnc_coupling_and_width_weights_match_independent_recompute():
    point = point_builder.make_point(quark_mass_basis_couplings=_fcnc_couplings())
    result = fcc.get(_PID).evaluate(point)
    sm_hadronic = sum(_manual_sm_hadronic_width_weights().values())
    sm_total = sum(_manual_sm_total_width_weights().values())

    assert sm_hadronic == pytest.approx(2.5225098333333333)
    assert sm_total == pytest.approx(3.649563333333334)
    assert sm_hadronic / sm_total == pytest.approx(0.6911812737414247)
    assert result.predicted == pytest.approx(0.0)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.ratio == pytest.approx(0.0)
    assert result.passes is True
    assert result.diagnostics["sm_validation"][
        "sm_hadronic_width_weight"
    ] == pytest.approx(sm_hadronic)
    assert result.diagnostics["sm_validation"]["sm_total_width_weight"] == pytest.approx(
        sm_total
    )
    assert result.diagnostics["sm_validation"][
        "sm_hadronic_to_total_width_ratio"
    ] == pytest.approx(sm_hadronic / sm_total)
    assert all(
        value == pytest.approx(0.0)
        for value in result.diagnostics["sm_validation"][
            "zero_fcnc_coupling_branching_fractions"
        ].values()
    )


def test_proxy_numerics_match_independent_effective_coupling_recomputation():
    constraint = fcc.get(_PID)
    left = 0.25 + 0.10j
    right = 0.05j
    point = point_builder.make_point(
        quark_mass_basis_couplings=_fcnc_couplings(bs_left=left, bs_right=right)
    )
    result = constraint.evaluate(point)
    (
        expected_br,
        delta_left,
        delta_right,
        norm,
        sm_hadronic,
        sm_total,
        fcnc_weight,
        hadronic_to_total,
    ) = (
        _manual_fcnc_br("b", "s", left, right)
    )
    core_result = down_fcnc_branching_fraction_from_couplings(
        flavor_i="b",
        flavor_j="s",
        delta_g_left=delta_left,
        delta_g_right=delta_right,
        br_limit=constraint.anchor.channels["bs"].value,
    )

    assert result.diagnostics["selected_channel"] == "bs"
    assert result.predicted == pytest.approx(expected_br)
    assert result.predicted == pytest.approx(core_result.branching_fraction)
    assert result.ratio == pytest.approx(expected_br / constraint.anchor.channels["bs"].value)
    assert result.ratio == pytest.approx(core_result.ratio_to_limit)
    assert result.diagnostics["channels"]["bs"]["delta_g_left"] == pytest.approx(
        delta_left
    )
    assert result.diagnostics["channels"]["bs"]["delta_g_right"] == pytest.approx(
        delta_right
    )
    assert result.diagnostics["channels"]["bs"]["effective_coupling_norm"] == pytest.approx(
        norm
    )
    assert result.diagnostics["channels"]["bs"]["sm_hadronic_width_weight"] == pytest.approx(
        sm_hadronic
    )
    assert result.diagnostics["channels"]["bs"]["sm_total_width_weight"] == pytest.approx(
        sm_total
    )
    assert result.diagnostics["channels"]["bs"]["fcnc_width_weight"] == pytest.approx(
        fcnc_weight
    )
    assert result.diagnostics["channels"]["bs"]["total_width_weight"] == pytest.approx(
        sm_total + fcnc_weight
    )
    assert result.diagnostics["channels"]["bs"][
        "sm_hadronic_to_total_width_ratio"
    ] == pytest.approx(hadronic_to_total)
    assert "SM total Z width weight" in result.diagnostics["branching_formula"]
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["channels"]["bs"]["proxy"][
        "matching_assumption"
    ]


def test_evaluate_without_or_with_invalid_input_degrades_gracefully():
    missing = fcc.get(_PID).evaluate(point_builder.empty_point())
    invalid = fcc.get(_PID).evaluate(
        point_builder.make_point(quark_mass_basis_couplings={"left_bs": 1.0})
    )

    assert missing.process_id == _PID
    assert missing.passes is True
    assert missing.predicted is None
    assert missing.ratio is None
    assert missing.sm_prediction == pytest.approx(0.0)
    assert missing.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert missing.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert missing.notes.startswith("NOT EVALUATED -")
    assert missing.diagnostics["evaluated"] is False
    assert missing.diagnostics["missing_extra"] == "quark_mass_basis_couplings"

    assert invalid.passes is True
    assert invalid.predicted is None
    assert invalid.ratio is None
    assert invalid.diagnostics["evaluated"] is False
    assert invalid.diagnostics["invalid_extra"] == "quark_mass_basis_couplings"


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(
            quark_mass_basis_couplings=_fcnc_couplings(
                bs_left=100.0,
                bd_left=0.10 + 0.02j,
                sd_right=0.05j,
            )
        )
    )

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
    assert isinstance(result.diagnostics["channels"]["bs"]["delta_g_left"], complex)
    assert isinstance(result.diagnostics["channels"]["bs"]["delta_g_right"], complex)
    assert isinstance(
        result.diagnostics["channels"]["bs"]["effective_coupling_norm"],
        float,
    )
    assert result.diagnostics["evaluated"] is True


@pytest.mark.parametrize(
    ("couplings", "expected_pass"),
    [
        (_fcnc_couplings(bs_left=1.0, bd_left=0.5, sd_left=0.25), True),
        (_fcnc_couplings(bs_left=100.0), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(
    couplings: QuarkMassBasisCouplings,
    expected_pass: bool,
):
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(quark_mass_basis_couplings=couplings)
    )

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_optional_kk_ew_mass_extra_changes_proxy_scale():
    couplings = _fcnc_couplings(bs_left=20.0)
    default_point = point_builder.make_point(quark_mass_basis_couplings=couplings)
    heavy_point = point_builder.make_point(
        quark_mass_basis_couplings=couplings,
        kk_ew_mass_gev=6000.0,
    )

    default_result = fcc.get(_PID).evaluate(default_point)
    heavy_result = fcc.get(_PID).evaluate(heavy_point)

    assert default_result.diagnostics["channels"]["bs"]["proxy"]["m_kk_gev"] == (
        pytest.approx(3000.0)
    )
    assert heavy_result.diagnostics["channels"]["bs"]["proxy"]["m_kk_gev"] == (
        pytest.approx(6000.0)
    )
    assert heavy_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert heavy_result.diagnostics["channels"]["bs"]["delta_g_left"] == pytest.approx(
        default_result.diagnostics["channels"]["bs"]["delta_g_left"] / 4.0
    )
    assert heavy_result.diagnostics["channels"]["bs"]["effective_coupling_norm"] == (
        pytest.approx(
            default_result.diagnostics["channels"]["bs"]["effective_coupling_norm"]
            / 16.0
        )
    )


def test_evaluate_is_pure_and_deterministic():
    couplings = _fcnc_couplings(
        bs_left=0.15 + 0.01j,
        bd_right=0.03j,
        sd_left=0.02,
    )
    before_left = couplings.left_overlap.copy()
    before_right = couplings.right_down_overlap.copy()
    point = point_builder.make_point(quark_mass_basis_couplings=couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_overlap, before_left)
    np.testing.assert_array_equal(couplings.right_down_overlap, before_right)
