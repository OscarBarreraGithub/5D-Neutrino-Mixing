"""Production tests for E001 (electron EDM)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.edm import lepton_edm_proxy_input
from flavor_catalog_constraints.primary.edm_neutrino import E001 as e001_module
from quarkConstraints.edm import edm_e_cm_from_cp_odd_dipole, evaluate_charged_lepton_edm

_PID = "E001"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "edm_neutrino" / "E001.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _proxy_point(
    cp_odd_dipole_coefficient_gev_inv: float,
    *,
    m_kk_gev: float | None = 3000.0,
):
    proxy = lepton_edm_proxy_input(
        cp_odd_dipole_coefficient_gev_inv,
        lepton="e",
        m_kk_gev=m_kk_gev,
        source="E001 test proxy",
    )
    return point_builder.make_point(lepton_mass_basis_couplings=proxy)


def _core_edm_result(
    cp_odd_dipole_coefficient_gev_inv: float,
    *,
    experimental_limit_e_cm: float,
):
    return evaluate_charged_lepton_edm(
        cp_odd_dipole_coefficient_gev_inv,
        experimental_limit_e_cm=experimental_limit_e_cm,
        lepton="e",
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "edm_neutrino"
    assert constraint.observable == "|d_e|"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    limit = pdg["canonical_limit"]

    assert constraint.anchor.experimental.value == pytest.approx(limit["value"])
    assert constraint.anchor.experimental.block_key == "canonical_limit"
    assert constraint.anchor.experimental.source_url == limit["source_url"]
    assert constraint.anchor.experimental.units == "e cm"
    assert constraint.anchor.budget == pytest.approx(4.1e-30)
    assert constraint.anchor.limit_operator == limit["limit_operator"]
    assert constraint.anchor.confidence_level == limit["confidence_level"]
    assert constraint.anchor.used_measurement == limit["used_measurement"]
    assert constraint.anchor.table_value == pytest.approx(limit["table_value"])
    assert constraint.anchor.table_units == limit["table_units"]

    with pytest.raises(anchors.AnchorError):
        anchors.load_anchor(
            _PID,
            family="edm_neutrino",
            candidates=("no_such_block",),
        )


def test_anchor_rejects_superseded_limit_block(monkeypatch):
    superseded = _yaml_pdg_block()["superseded_acme_2018"]

    def superseded_load_anchor(*args, **kwargs):
        return anchors.Anchor(
            process_id=_PID,
            block_key="superseded_acme_2018",
            value=superseded["value"],
            uncertainty=None,
            observable=superseded["observable"],
            units=superseded["units"],
            source=superseded["source"],
            source_url=superseded["source_url"],
            year=superseded["year"],
            snapshot_path=superseded["snapshot_path"],
        )

    monkeypatch.setattr(e001_module, "load_anchor", superseded_load_anchor)
    with pytest.raises(anchors.AnchorError, match="canonical_limit"):
        e001_module._load_e001_anchor(_PID)


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
        point_builder.make_point(lepton_mass_basis_couplings={"y_n_bar": [1.0, 2.0]})
    )

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["invalid_extra"] == "lepton_mass_basis_couplings"
    assert result.diagnostics["exception_type"] in {"KeyError", "ValueError", "TypeError"}


def test_proxy_numerics_match_core_edm_conversion():
    constraint = fcc.get(_PID)
    coefficient = 1.0e-16
    result = constraint.evaluate(_proxy_point(coefficient))
    expected = _core_edm_result(
        coefficient,
        experimental_limit_e_cm=constraint.anchor.budget,
    )

    assert result.predicted == pytest.approx(expected.abs_edm_e_cm)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.ratio == pytest.approx(expected.ratio_to_limit)
    assert result.diagnostics["edm_e_cm_signed"] == pytest.approx(
        edm_e_cm_from_cp_odd_dipole(coefficient)
    )
    assert result.diagnostics["cp_odd_dipole_coefficient_gev_inv"] == pytest.approx(
        coefficient
    )
    assert result.diagnostics["hbarc_gev_cm"] == pytest.approx(expected.hbarc_gev_cm)
    assert result.diagnostics["evaluated"] is True
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_complex_chiral_proxy_uses_imaginary_part():
    constraint = fcc.get(_PID)
    coefficient = 3.0e-16 + 1.0e-16j
    point = point_builder.make_point(
        lepton_mass_basis_couplings={
            "electron_chiral_dipole_coefficient_gev_inv": coefficient,
            "m_kk_gev": 3000.0,
            "source": "E001 chiral test proxy",
        }
    )
    result = constraint.evaluate(point)
    expected = _core_edm_result(
        coefficient.imag,
        experimental_limit_e_cm=constraint.anchor.budget,
    )

    assert result.predicted == pytest.approx(expected.abs_edm_e_cm)
    assert result.ratio == pytest.approx(expected.ratio_to_limit)
    assert result.diagnostics["cp_odd_projection"] == "imaginary_part"
    assert result.diagnostics["chiral_dipole_coefficient_gev_inv"] == coefficient
    assert isinstance(result.diagnostics["chiral_dipole_coefficient_gev_inv"], complex)


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(
            lepton_mass_basis_couplings={
                "electron_chiral_dipole_coefficient_gev_inv": 2.0e-17 + 4.0e-17j,
                "source": "E001 finite-fields proxy",
            }
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
    assert isinstance(result.diagnostics["chiral_dipole_coefficient_gev_inv"], complex)
    for key in (
        "edm_e_cm_signed",
        "abs_edm_e_cm",
        "cp_odd_dipole_coefficient_gev_inv",
        "hbarc_gev_cm",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["used_proxy"] is True


@pytest.mark.parametrize(
    ("point", "expected_pass"),
    [
        (_proxy_point(1.0e-16), True),
        (_proxy_point(1.0e-15), False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(point, expected_pass: bool):
    result = fcc.get(_PID).evaluate(point)

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 1.0


def test_optional_kk_ew_mass_extra_is_diagnostic_only_for_direct_proxy():
    proxy = lepton_edm_proxy_input(
        1.0e-16,
        lepton="e",
        m_kk_gev=3000.0,
        source="E001 test proxy",
    )
    default_point = point_builder.make_point(lepton_mass_basis_couplings=proxy)
    heavy_point = point_builder.make_point(
        lepton_mass_basis_couplings=proxy,
        kk_ew_mass_gev=6000.0,
    )

    default_result = fcc.get(_PID).evaluate(default_point)
    heavy_result = fcc.get(_PID).evaluate(heavy_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert heavy_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert heavy_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert heavy_result.diagnostics["m_kk_diagnostic_only"] is True
    assert heavy_result.predicted == pytest.approx(default_result.predicted)


def test_evaluate_is_pure_and_deterministic():
    proxy = lepton_edm_proxy_input(
        1.0e-16,
        lepton="e",
        m_kk_gev=3000.0,
        source="E001 test proxy",
    )
    point = point_builder.make_point(lepton_mass_basis_couplings=proxy)
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.get_extra("lepton_mass_basis_couplings") == proxy
