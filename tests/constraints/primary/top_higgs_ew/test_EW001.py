"""Production tests for EW001 (S,T,U oblique parameters)."""

from __future__ import annotations

import math
import re
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.anchors import Anchor, AnchorError
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.primary.top_higgs_ew import EW001 as ew001_module
from quarkConstraints import oblique_stu as core

_PID = "EW001"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "top_higgs_ew" / "EW001.yaml"


def _yaml_entries():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _entry(observable: str):
    matches = [entry for entry in _yaml_entries() if entry.get("observable") == observable]
    assert len(matches) == 1
    return matches[0]


def _manual_chi2(*, s: float, t: float, fit: core.ObliqueSTFit) -> float:
    ds = float(s - fit.s_central)
    dt = float(t - fit.t_central)
    var_s = float(fit.sigma_s**2)
    var_t = float(fit.sigma_t**2)
    cov_st = float(fit.rho_st * fit.sigma_s * fit.sigma_t)
    det = var_s * var_t - cov_st * cov_st
    return float((var_t * ds * ds - 2.0 * cov_st * ds * dt + var_s * dt * dt) / det)


def _solved_floor_tev(*, ew_model: str = "minimal_rs") -> float:
    constraint = fcc.get(_PID)
    lo = 100.0
    hi = 100000.0

    def chi2(m_kk_gev: float) -> float:
        return core.evaluate_rs_oblique_proxy(
            m_kk_gev=m_kk_gev,
            fit=constraint.anchor.fit,
            s_coefficient=constraint.anchor.warped_s_coefficient.value,
            ew_model=ew_model,
        ).chi2

    assert chi2(lo) > core.CHI2_2DOF_95
    assert chi2(hi) < core.CHI2_2DOF_95
    for _ in range(160):
        mid = 0.5 * (lo + hi)
        if chi2(mid) > core.CHI2_2DOF_95:
            lo = mid
        else:
            hi = mid
    return hi / 1000.0


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "top_higgs_ew"
    assert constraint.observable == "S,T,U oblique electroweak parameters"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    s_entry = _entry("PDG 2025 global oblique fit, S with U fixed to zero")
    t_entry = _entry("PDG 2025 global oblique fit, T with U fixed to zero")
    u_entry = _entry("PDG 2025 global oblique fit, U fixed value")
    rho_entry = _entry("PDG 2025 correlation rho(S,T), U fixed to zero")
    coeff_entry = _entry("PDG 2025 warped-extra-dimension S coefficient")
    upper_entry = _entry("PDG 2025 one-sided S upper value used for warped M_KK context")
    mkk_entry = _entry("PDG 2025 warped M_KK lower-bound context")

    assert constraint.anchor.fit.s_central == pytest.approx(s_entry["value"])
    assert constraint.anchor.fit.t_central == pytest.approx(t_entry["value"])
    assert constraint.anchor.fit.u_fixed == pytest.approx(u_entry["value"])
    assert constraint.anchor.fit.sigma_s == pytest.approx(s_entry["uncertainty"])
    assert constraint.anchor.fit.sigma_t == pytest.approx(t_entry["uncertainty"])
    assert constraint.anchor.fit.rho_st == pytest.approx(rho_entry["value"])
    assert constraint.anchor.fit.source_url == s_entry["source_url"]
    assert constraint.anchor.warped_s_coefficient.value == pytest.approx(
        coeff_entry["value"]
    )
    assert constraint.anchor.warped_s_coefficient.units == coeff_entry["units"]
    assert constraint.anchor.warped_s_upper_95.value == pytest.approx(
        upper_entry["value"]
    )
    assert constraint.anchor.warped_mkk_context.value == pytest.approx(
        mkk_entry["value"]
    )
    assert constraint.anchor.budget == pytest.approx(core.CHI2_2DOF_95)
    gfitter_s = _entry("Gfitter public oblique fit, S")
    gfitter_t = _entry("Gfitter public oblique fit, T")
    gfitter_u = _entry("Gfitter public oblique fit, U")
    assert gfitter_s["year"] == 2018
    assert gfitter_t["year"] == 2018
    assert gfitter_u["year"] == 2018

    with pytest.raises(AnchorError):
        ew001_module._load_value_anchor("not a real EW001 observable", process_id=_PID)


def test_list_anchors_route_through_scaffold_load_anchor(monkeypatch):
    calls: list[tuple[str, ...]] = []
    original_load_anchor = ew001_module.load_anchor

    def spy_load_anchor(*args, **kwargs):
        calls.append(tuple(kwargs["candidates"]))
        return original_load_anchor(*args, **kwargs)

    monkeypatch.setattr(ew001_module, "load_anchor", spy_load_anchor)
    anchor = ew001_module._load_ew001_anchor(_PID)

    assert calls == [
        ("pdg_or_equivalent[0]",),
        ("pdg_or_equivalent[1]",),
        ("pdg_or_equivalent[2]",),
        ("pdg_or_equivalent[3]",),
        ("pdg_or_equivalent[7]",),
        ("pdg_or_equivalent[8]",),
        ("pdg_or_equivalent[9]",),
    ]
    assert anchor.fit.s_central == pytest.approx(fcc.get(_PID).anchor.fit.s_central)

    def mismatched_load_anchor(*args, **kwargs):
        return Anchor(
            process_id=_PID,
            block_key="pdg_or_equivalent[99]",
            value=1.0,
            uncertainty=0.1,
        )

    monkeypatch.setattr(ew001_module, "load_anchor", mismatched_load_anchor)
    with pytest.raises(AnchorError, match="load_anchor selected"):
        ew001_module._load_value_anchor(
            "PDG 2025 global oblique fit, S with U fixed to zero",
            process_id=_PID,
        )


def test_sm_reference_and_numerical_chi2_match_independent_recomputation():
    constraint = fcc.get(_PID)
    fit = constraint.anchor.fit
    sm_manual = _manual_chi2(s=0.0, t=0.0, fit=fit)
    result = constraint.evaluate(point_builder.make_point(kk_ew_mass_gev=6000.0))

    scale = (core.DEFAULT_HIGGS_VEV_GEV / 6000.0) ** 2
    s_proxy = constraint.anchor.warped_s_coefficient.value * scale
    t_proxy = core.minimal_rs_t_coefficient() * scale
    expected = _manual_chi2(s=s_proxy, t=t_proxy, fit=fit)

    assert result.sm_prediction == pytest.approx(sm_manual)
    # sm_prediction (S=T=0 reference) is M_KK-independent -> unchanged by M2.
    assert result.sm_prediction == pytest.approx(0.9627934850901363)
    # ``expected`` is recomputed in-test from core.minimal_rs_t_coefficient(),
    # which now carries the x1^2 physical-convention factor (M2), so this
    # self-consistency check still holds.  The chi2 literal is re-pinned: the
    # corrected Delta T is ~x1^2 ~ 6 larger, and chi2 is quadratic in (S, T),
    # so the 6 TeV chi2 moved 3.6313 -> 519.5502 (PLAN §5/M2).
    assert result.predicted == pytest.approx(expected)
    assert result.predicted == pytest.approx(519.5501818374169)
    assert result.ratio == pytest.approx(expected / core.CHI2_2DOF_95)
    assert result.diagnostics["s_prediction"] == pytest.approx(s_proxy)
    assert result.diagnostics["t_prediction"] == pytest.approx(t_proxy)
    assert result.diagnostics["u_prediction"] == pytest.approx(0.0)


def test_evaluate_without_input_degrades_gracefully():
    result = fcc.get(_PID).evaluate(point_builder.empty_point())

    assert result.process_id == _PID
    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.experimental == pytest.approx(0.0)
    assert result.budget == pytest.approx(core.CHI2_2DOF_95)
    assert result.diagnostics["missing_extras"] == (
        "kk_ew_mass_gev",
        "kk_gluon_mass_gev",
        "quark_mass_basis_couplings",
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


@pytest.mark.parametrize(
    ("m_kk_gev", "expected_pass"),
    [
        # Re-pinned after M2 (Delta T x1^2 physical-convention factor): the
        # minimal-RS oblique floor moved up to ~16 TeV (the proxy is excluded
        # from the STRICT lane, but this gate is now convention-honest).  6 TeV
        # now FAILS (it passed only because Delta T was ~6 too small).
        (18000.0, True),
        (6000.0, False),
        (3000.0, False),
    ],
)
def test_safe_point_passes_and_excluded_point_fails(
    m_kk_gev: float,
    expected_pass: bool,
):
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(kk_ew_mass_gev=m_kk_gev)
    )

    assert result.passes is expected_pass
    assert result.predicted is not None
    assert result.budget == pytest.approx(core.CHI2_2DOF_95)
    assert result.experimental == pytest.approx(0.0)
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_c3_documented_floor_matches_code_solved_floor():
    floor_tev = _solved_floor_tev()
    custodial_floor_tev = _solved_floor_tev(ew_model="custodial_rs_plr")

    assert floor_tev == pytest.approx(15.955484201013507, rel=1e-12)
    assert custodial_floor_tev == pytest.approx(5.796965450093533, rel=1e-12)

    for path in (_REPO_ROOT / "docs" / "FLOOR_SUMMARY.md", _REPO_ROOT / "CLAUDE.md"):
        text = path.read_text()
        match = re.search(r"code-verified\s+\*\*([0-9]+\.[0-9]+)\s+TeV\*\*", text)
        assert match, f"{path} is missing the C-3 code-verified EW001 floor"
        assert float(match.group(1)) == pytest.approx(round(floor_tev, 2))
        assert "18-20 TeV" not in text


def test_mass_resolution_prefers_ew_mass_and_falls_back_to_couplings():
    ew_point = point_builder.make_point(
        kk_ew_mass_gev=6000.0,
        kk_gluon_mass_gev=3000.0,
        quark_mass_basis_couplings=SimpleNamespace(M_KK=2500.0),
    )
    fallback_point = point_builder.make_point(
        quark_mass_basis_couplings=SimpleNamespace(M_KK=6000.0)
    )

    ew_result = fcc.get(_PID).evaluate(ew_point)
    fallback_result = fcc.get(_PID).evaluate(fallback_point)

    assert ew_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert ew_result.diagnostics["mass_source"] == "kk_ew_mass_gev"
    assert fallback_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert fallback_result.diagnostics["mass_source"] == (
        "quark_mass_basis_couplings.M_KK"
    )
    assert fallback_result.predicted == pytest.approx(ew_result.predicted)


def test_invalid_mass_returns_non_crashing_failure():
    result = fcc.get(_PID).evaluate(point_builder.make_point(kk_ew_mass_gev=-1.0))

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["invalid_m_kk_gev"] == pytest.approx(-1.0)


def test_evaluate_is_pure_and_deterministic_with_real_finite_fields():
    point = point_builder.make_point(kk_ew_mass_gev=6000.0)
    before = dict(point.extras)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.extras == before
    for value in (
        first.predicted,
        first.ratio,
        first.budget,
        first.sm_prediction,
        first.experimental,
    ):
        assert isinstance(value, float)
        assert math.isfinite(value)
    assert "NEEDS-HUMAN-PHYSICS" in first.diagnostics["needs_human_physics"]
