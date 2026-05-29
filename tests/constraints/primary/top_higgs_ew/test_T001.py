"""Production tests for T001 (top FCNC t -> c Z)."""

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
from flavor_catalog_constraints.primary.top_higgs_ew import T001 as t001_module
from quarkConstraints.couplings import QuarkMassBasisCouplings

_PID = "T001"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "top_higgs_ew" / "T001.yaml"
_NUMBER_RE = re.compile(
    r"(?P<number>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(?P<pct>%)?"
)


def _yaml_entries_by_id():
    with open(_SIDECAR) as handle:
        values = yaml.safe_load(handle)["pdg_or_equivalent"]["values"]
    return {entry["value_id"]: entry for entry in values if "value_id" in entry}


def _parse_yaml_number(value: object) -> float:
    match = _NUMBER_RE.search(str(value))
    if match is None:
        raise AssertionError(f"could not parse numeric value from {value!r}")
    number = float(match.group("number"))
    if match.group("pct"):
        number *= 1.0e-2
    return number


def _ct_couplings(
    left: complex,
    right: complex = 0.0j,
    M_KK: float = 3000.0,
    g_s: float = 1.0,
) -> QuarkMassBasisCouplings:
    """Minimal mass-basis couplings with only the c-t up-sector slot populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_up = zeros.copy()
    right_up = zeros.copy()
    left_up[1, 2] = left
    left_up[2, 1] = np.conj(left)
    right_up[1, 2] = right
    right_up[2, 1] = np.conj(right)
    return QuarkMassBasisCouplings(
        M_KK=M_KK,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=g_s,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=left_up,
        left_down=zeros,
        right_up=right_up,
        right_down=zeros,
    )


def _manual_t_to_c_z_proxy_br(
    *,
    left: complex,
    right: complex,
    M_KK: float,
    g_s: float,
    inputs,
) -> float:
    g_weak = math.sqrt(4.0 * math.pi * inputs.alpha_em_mz / inputs.sin2_theta_w)
    cos_theta_w = math.sqrt(1.0 - inputs.sin2_theta_w)
    g_z = g_weak / (2.0 * cos_theta_w)
    suppression = (inputs.m_z_gev / M_KK) ** 2
    vector_left = left / g_s * suppression
    vector_right = right / g_s * suppression
    rho = (inputs.m_z_gev / inputs.m_top_gev) ** 2
    phase = (1.0 - rho) ** 2 * (1.0 + 2.0 * rho)
    partial_width = (
        g_z**2
        * inputs.m_top_gev**3
        / (32.0 * math.pi * inputs.m_z_gev**2)
        * phase
        * (abs(vector_left) ** 2 + abs(vector_right) ** 2)
    )
    return float(partial_width / inputs.total_top_width_gev)


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "top_higgs_ew"
    assert constraint.observable == "BR(t -> c Z)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    entries = _yaml_entries_by_id()
    pdg_left = entries["PDG2025:T001:tZc_left"]
    pdg_right = entries["PDG2025:T001:tZc_right"]
    atlas_left = entries["ATLAS2023:T001:tZc_left"]
    cms = entries["CMS2017:T001:tZc"]
    sm = entries["AguilarSaavedra2004:T001:SM"]

    assert constraint.anchor.pdg_tzc_left.value == pytest.approx(
        _parse_yaml_number(pdg_left["value"])
    )
    assert constraint.anchor.pdg_tzc_left.source_url == pdg_left["source_url"]
    assert constraint.anchor.pdg_tzc_right.value == pytest.approx(
        _parse_yaml_number(pdg_right["value"])
    )
    assert constraint.anchor.atlas_tzc_left.value == pytest.approx(
        _parse_yaml_number(atlas_left["value"])
    )
    assert constraint.anchor.cms_tzc.value == pytest.approx(
        _parse_yaml_number(cms["normalized_value"])
    )
    assert constraint.anchor.sm_value == pytest.approx(_parse_yaml_number(sm["value"]))
    assert constraint.anchor.budget == pytest.approx(1.2e-4)
    assert constraint.anchor.active_limit.value_id == "PDG2025:T001:tZc_right"

    with pytest.raises(AnchorError):
        t001_module._load_limit_anchor("PDG2025:T001:not_present", process_id=_PID)


def test_value_list_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = t001_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(t001_module, "load_anchor", spy_load_anchor)
    anchor = t001_module._load_t001_anchor(_PID)

    assert calls == [
        ("pdg_or_equivalent.values[0]",),
        ("pdg_or_equivalent.values[1]",),
        ("pdg_or_equivalent.values[3]",),
        ("pdg_or_equivalent.values[4]",),
        ("pdg_or_equivalent.values[6]",),
        ("pdg_or_equivalent.values[7]",),
        ("pdg_or_equivalent.values[8]",),
    ]
    assert anchor.budget == pytest.approx(fcc.get(_PID).anchor.budget)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent.values[99]",
            value=1.0,
            uncertainty=None,
        )

    monkeypatch.setattr(t001_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        t001_module._load_limit_anchor("PDG2025:T001:tZc_left", process_id=_PID)


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.sm_prediction == pytest.approx(1.0e-14)
    assert result.budget == pytest.approx(1.2e-4)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_sm_limit_and_manual_width_recomputation():
    constraint = fcc.get(_PID)
    couplings = _ct_couplings(left=0.0j, right=0.0j)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)

    assert result.predicted == pytest.approx(0.0)
    assert result.sm_prediction == pytest.approx(1.0e-14)
    assert result.diagnostics["total_branching_fraction_including_sm"] == pytest.approx(
        1.0e-14
    )
    assert result.ratio == pytest.approx(0.0)
    assert result.passes is True

    nonzero = 1.0 + 0.25j
    couplings = _ct_couplings(left=nonzero, right=0.3j, g_s=1.0)
    point = point_builder.build_from_quark_couplings(couplings)
    result = constraint.evaluate(point)
    manual = _manual_t_to_c_z_proxy_br(
        left=nonzero,
        right=0.3j,
        M_KK=3000.0,
        g_s=1.0,
        inputs=constraint.sm_inputs,
    )

    assert manual == pytest.approx(4.817352363171973e-7)
    assert result.predicted == pytest.approx(manual)
    assert result.ratio == pytest.approx(manual / constraint.anchor.budget)
    assert result.diagnostics["z_mixing_suppression"] == pytest.approx(
        (constraint.sm_inputs.m_z_gev / 3000.0) ** 2
    )


def test_safe_point_passes_and_large_np_point_fails():
    constraint = fcc.get(_PID)
    safe = constraint.evaluate(
        point_builder.build_from_quark_couplings(_ct_couplings(left=1.0))
    )
    excluded = constraint.evaluate(
        point_builder.build_from_quark_couplings(_ct_couplings(left=25.0))
    )

    assert safe.passes is True
    assert safe.ratio < 1.0
    assert excluded.passes is False
    assert excluded.ratio > 1.0
    assert excluded.predicted > constraint.anchor.budget


def test_optional_kk_ew_mass_extra_changes_proxy_scale():
    couplings = _ct_couplings(left=2.0)
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
    assert ew_result.predicted == pytest.approx(default_result.predicted / 16.0)


def test_evaluate_has_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(
        point_builder.build_from_quark_couplings(
            _ct_couplings(left=1.0e-1 + 2.0e-2j, right=5.0e-2j)
        )
    )

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
        "left_qt_coupling",
        "right_qt_coupling",
        "left_qt_overlap",
        "right_qt_overlap",
        "vector_left",
        "vector_right",
    ):
        assert isinstance(result.diagnostics[key], complex)
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["sm_is_negligible_for_limit"] is True


def test_evaluate_is_pure_and_deterministic():
    couplings = _ct_couplings(left=1.0e-1 + 2.0e-2j, right=5.0e-2j)
    before_left_up = couplings.left_up.copy()
    before_right_up = couplings.right_up.copy()
    point = point_builder.build_from_quark_couplings(couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_up, before_left_up)
    np.testing.assert_array_equal(couplings.right_up, before_right_up)
