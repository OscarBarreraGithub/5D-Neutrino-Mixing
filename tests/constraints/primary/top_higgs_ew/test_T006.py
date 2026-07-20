"""Production tests for T006 (top FCNC t -> u g)."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.top_higgs_ew import T006 as t006_module
from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.top_fcnc import gluon_dipole_branching_fraction

_PID = "T006"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "top_higgs_ew" / "T006.yaml"
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


def _ut_couplings(
    left: complex,
    right: complex = 0.0j,
    M_KK: float = 3000.0,
    g_s: float = 1.0,
) -> QuarkMassBasisCouplings:
    """Minimal mass-basis couplings with only the u-t up-sector slot populated."""
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_up = zeros.copy()
    right_up = zeros.copy()
    left_up[0, 2] = left
    left_up[2, 0] = np.conj(left)
    right_up[0, 2] = right
    right_up[2, 0] = np.conj(right)
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


def _core_t_to_u_gluon_recomputation(
    *,
    couplings: QuarkMassBasisCouplings,
    inputs,
    m_kk_gev: float | None = None,
):
    left_qt = complex(np.asarray(couplings.left_up, dtype=np.complex128)[0, 2])
    right_qt = complex(np.asarray(couplings.right_up, dtype=np.complex128)[0, 2])
    resolved_m_kk = float(couplings.M_KK if m_kk_gev is None else m_kk_gev)
    g_s = float(couplings.g_s)
    suppression = (inputs.m_top_gev / resolved_m_kk) ** 2
    return gluon_dipole_branching_fraction(
        light_quark="u",
        dipole_left=left_qt / g_s * suppression,
        dipole_right=right_qt / g_s * suppression,
        inputs=inputs,
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "top_higgs_ew"
    assert constraint.observable == "BR(t -> u g)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    entries = _yaml_entries_by_id()
    atlas_tug = entries["PDG2025:T006:ATLAS2022:t_ug"]
    atlas_cugut = entries["PDG2025:T006:ATLAS2022:CuGut"]
    atlas_xsec = entries["PDG2025:T006:ATLAS2022:ug_to_t_cross_section"]
    cms_tug = entries["CMS2017:T006:t_ug"]
    cms_kappa = entries["CMS2017:T006:kappa_tug"]
    sm = entries["AguilarSaavedra2004:T006:t_ug_SM"]

    assert constraint.anchor.atlas_tug.value == pytest.approx(
        _parse_yaml_number(atlas_tug["value"])
    )
    assert constraint.anchor.atlas_tug.source_url == atlas_tug["source_url"]
    assert constraint.anchor.atlas_cugut.value == pytest.approx(
        _parse_yaml_number(atlas_cugut["value"])
    )
    assert constraint.anchor.atlas_ug_to_t_cross_section.value == pytest.approx(
        _parse_yaml_number(atlas_xsec["value"])
    )
    assert constraint.anchor.cms_tug.value == pytest.approx(
        _parse_yaml_number(cms_tug["value"])
    )
    assert constraint.anchor.cms_kappa_tug.value == pytest.approx(
        _parse_yaml_number(cms_kappa["value"])
    )
    assert constraint.anchor.sm_value == pytest.approx(_parse_yaml_number(sm["value"]))
    assert constraint.anchor.budget == pytest.approx(
        min(_parse_yaml_number(atlas_tug["value"]), _parse_yaml_number(cms_tug["value"]))
    )
    assert constraint.anchor.active_limit.value_id == "CMS2017:T006:t_ug"

    with pytest.raises(AnchorError):
        t006_module._load_limit_anchor("PDG2025:T006:not_present", process_id=_PID)


def test_value_list_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = t006_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(t006_module, "load_anchor", spy_load_anchor)
    anchor = t006_module._load_t006_anchor(_PID)

    assert calls == [
        ("pdg_or_equivalent.values[4]",),
        ("pdg_or_equivalent.values[3]",),
        ("pdg_or_equivalent.values[1]",),
        ("pdg_or_equivalent.values[7]",),
        ("pdg_or_equivalent.values[6]",),
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

    monkeypatch.setattr(t006_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        t006_module._load_limit_anchor(
            "PDG2025:T006:ATLAS2022:t_ug",
            process_id=_PID,
        )


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())
    sm = _parse_yaml_number(_yaml_entries_by_id()["AguilarSaavedra2004:T006:t_ug_SM"]["value"])

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(fcc.get(_PID).anchor.value)
    assert result.sm_prediction == pytest.approx(sm)
    assert result.budget == pytest.approx(fcc.get(_PID).anchor.budget)
    assert result.diagnostics["missing_extra"] == "quark_mass_basis_couplings"
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_sm_limit_and_core_width_recomputation():
    constraint = fcc.get(_PID)
    sm = _parse_yaml_number(_yaml_entries_by_id()["AguilarSaavedra2004:T006:t_ug_SM"]["value"])
    couplings = _ut_couplings(left=0.0j, right=0.0j)
    point = point_builder.make_point(quark_mass_basis_couplings=couplings)
    result = constraint.evaluate(point)

    assert result.predicted == pytest.approx(0.0)
    assert result.sm_prediction == pytest.approx(sm)
    assert result.diagnostics["total_branching_fraction_including_sm"] == pytest.approx(
        sm
    )
    assert result.ratio == pytest.approx(0.0)
    assert result.passes is True

    nonzero = 1.0 + 0.25j
    couplings = _ut_couplings(left=nonzero, right=0.3j, g_s=1.0)
    point = point_builder.make_point(quark_mass_basis_couplings=couplings)
    result = constraint.evaluate(point)
    core = _core_t_to_u_gluon_recomputation(
        couplings=couplings,
        inputs=constraint.sm_inputs,
    )

    assert core.branching_fraction == pytest.approx(5.590559579388595e-5)
    assert result.predicted == pytest.approx(core.branching_fraction)
    assert result.ratio == pytest.approx(core.branching_fraction / constraint.anchor.budget)
    assert result.diagnostics["dipole_scale_suppression"] == pytest.approx(
        (constraint.sm_inputs.m_top_gev / 3000.0) ** 2
    )
    assert result.diagnostics["alpha_s_mt"] == pytest.approx(
        constraint.sm_inputs.alpha_s_mt
    )
    assert result.diagnostics["color_factor"] == pytest.approx(4.0 / 3.0)


def test_safe_point_passes_and_large_np_point_fails():
    constraint = fcc.get(_PID)
    safe = constraint.evaluate(
        point_builder.make_point(quark_mass_basis_couplings=_ut_couplings(left=0.4))
    )
    excluded = constraint.evaluate(
        point_builder.make_point(quark_mass_basis_couplings=_ut_couplings(left=0.7))
    )

    assert safe.passes is True
    assert safe.ratio < 1.0
    assert excluded.passes is False
    assert excluded.ratio > 1.0
    assert excluded.predicted > constraint.anchor.budget


def test_optional_kk_gluon_mass_extra_changes_proxy_scale():
    couplings = _ut_couplings(left=0.5)
    default_point = point_builder.make_point(quark_mass_basis_couplings=couplings)
    gluon_point = point_builder.make_point(
        quark_mass_basis_couplings=couplings,
        kk_gluon_mass_gev=6000.0,
    )
    default_result = fcc.get(_PID).evaluate(default_point)
    gluon_result = fcc.get(_PID).evaluate(gluon_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert default_result.diagnostics["kk_gluon_mass_extra_used"] is False
    assert gluon_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert gluon_result.diagnostics["kk_gluon_mass_extra_used"] is True
    assert gluon_result.predicted == pytest.approx(default_result.predicted / 16.0)


def test_evaluate_has_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(
            quark_mass_basis_couplings=_ut_couplings(
                left=1.0e-1 + 2.0e-2j,
                right=5.0e-2j,
            )
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
        "dipole_left",
        "dipole_right",
    ):
        assert isinstance(result.diagnostics[key], complex)
    assert result.diagnostics["light_quark"] == "u"
    assert result.diagnostics["light_up_index"] == 0
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]
    assert result.diagnostics["sm_is_negligible_for_limit"] is True


def test_evaluate_is_pure_and_deterministic():
    couplings = _ut_couplings(left=1.0e-1 + 2.0e-2j, right=5.0e-2j)
    before_left_up = couplings.left_up.copy()
    before_right_up = couplings.right_up.copy()
    point = point_builder.make_point(quark_mass_basis_couplings=couplings)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    np.testing.assert_array_equal(couplings.left_up, before_left_up)
    np.testing.assert_array_equal(couplings.right_up, before_right_up)
