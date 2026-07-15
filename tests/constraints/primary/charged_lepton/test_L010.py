"""Production tests for L010 (tau -> 3e)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.lfv_three_body_taue import (
    tau_to_3e_from_lepton_input,
)
from quarkConstraints.lfv_three_body import (
    LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION,
    TAU_TO_E_NUNU_BRANCHING_FRACTION,
    lfv_three_body_contact_amplitudes,
    lfv_three_body_from_components,
)

_PID = "L010"
_BELLEII_2025_LIMIT = 2.5e-8
_BELLEII_2025_ARXIV = "2507.18236"
_BELLEII_2025_JOURNAL = "JHEP 12 (2025) 169"
_BELLEII_2025_URL = "https://arxiv.org/abs/2507.18236"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charged_lepton" / "L010.yaml"


def _yaml_pdg_block():
    with open(_SIDECAR) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _point_from_mapping(mapping):
    return point_builder.make_point(lepton_mass_basis_couplings=dict(mapping))


def _core_tau3e_prediction(
    constraint,
    lepton: dict,
    *,
    m_kk_gev: float,
):
    contact = lfv_three_body_contact_amplitudes(
        lepton,
        initial_flavor="tau",
        final_flavor="e",
        m_kk_gev=m_kk_gev,
        inputs=constraint.sm_inputs,
    )
    return lfv_three_body_from_components(
        dipole_parent_branching_fraction=float(
            lepton.get("dipole_parent_branching_fraction", 0.0)
        ),
        contact_amplitudes=contact,
        br_limit=constraint.anchor.budget,
        inputs=constraint.sm_inputs,
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charged_lepton"
    assert constraint.observable == "BR(tau -> 3e)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["primary_current_limit"]
    primary = pdg["primary_experiment"]
    babar = pdg["supporting_experiment"]

    assert exp["value"] == pytest.approx(_BELLEII_2025_LIMIT)
    assert exp["source"] == "Belle II Collaboration"
    assert exp["year"] == 2025
    assert exp["arxiv_id"] == _BELLEII_2025_ARXIV
    assert exp["journal"] == _BELLEII_2025_JOURNAL
    assert exp["source_url"] == _BELLEII_2025_URL
    assert exp["journal_url"] == "https://doi.org/10.1007/JHEP12(2025)169"
    assert exp["data_sample"] == "428 fb^-1"
    assert exp["observed_limit_90cl"] == pytest.approx(_BELLEII_2025_LIMIT)
    assert exp["expected_limit_90cl"] == pytest.approx(2.7e-8)
    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.primary_experiment.value == pytest.approx(primary["value"])
    assert constraint.anchor.primary_experiment.source_url == primary["source_url"]
    assert constraint.anchor.supporting_experiment.value == pytest.approx(babar["value"])
    assert constraint.anchor.supporting_experiment.source_url == babar["source_url"]
    assert constraint.anchor.budget == pytest.approx(_BELLEII_2025_LIMIT)
    assert constraint.anchor.prospects.source_url == pdg["prospects"]["source_url"]
    assert constraint.anchor.lhcb_context.source_url == pdg["lhcb_context"]["source_url"]

    with pytest.raises(anchors.AnchorError):
        anchors.load_anchor(
            _PID,
            family="charged_lepton",
            candidates=("no_such_block",),
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


def test_invalid_lepton_input_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(lepton_mass_basis_couplings={"source": "empty"})
    )

    assert result.passes is False
    assert result.predicted is None
    assert result.ratio is None
    assert result.notes.startswith("NOT EVALUATED -")
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["invalid_extra"] == "lepton_mass_basis_couplings"
    assert result.diagnostics["exception_type"] == "TypeError"


def test_tau_to_3e_rejects_caller_flavor_override():
    constraint = fcc.get(_PID)
    bad = {
        "initial_flavor": "mu",
        "final_flavor": "e",
        "left_lfv_overlap": 2.0e-4,
        "m_kk_gev": 3000.0,
        "source": "bad mu->e probe",
    }

    with pytest.raises(ValueError, match="pinned"):
        tau_to_3e_from_lepton_input(
            bad,
            br_limit=constraint.anchor.budget,
            inputs=constraint.sm_inputs,
        )

    result = constraint.evaluate(_point_from_mapping(bad))
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["exception_type"] == "ValueError"
    assert "pinned to initial_flavor='tau'" in result.diagnostics["exception"]


@pytest.mark.parametrize("bad_key", ["left_emu_overlap", "left_taumu_overlap"])
def test_tau_to_3e_rejects_mismatched_spurion_aliases(bad_key):
    constraint = fcc.get(_PID)
    bad = {
        bad_key: 1.0e-2,
        "m_kk_gev": 3000.0,
        "source": "bad flavor-alias probe",
    }

    with pytest.raises(ValueError, match="does not accept"):
        tau_to_3e_from_lepton_input(
            bad,
            br_limit=constraint.anchor.budget,
            inputs=constraint.sm_inputs,
        )

    result = constraint.evaluate(_point_from_mapping(bad))
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["exception_type"] == "ValueError"
    assert bad_key in result.diagnostics["exception"]


def test_tau_e_dipole_only_factor_matches_manual_formula():
    constraint = fcc.get(_PID)
    dipole_parent_br = 1.0e-8
    lepton = {
        "dipole_parent_branching_fraction": dipole_parent_br,
        "source": "L010 dipole-only test proxy",
    }
    result = constraint.evaluate(_point_from_mapping(lepton))
    m_tau = constraint.sm_inputs.charged_lepton_mass("tau")
    m_e = constraint.sm_inputs.charged_lepton_mass("e")
    expected_factor = (
        constraint.sm_inputs.alpha_em
        / (3.0 * math.pi)
        * (math.log((m_tau / m_e) ** 2) - 11.0 / 4.0)
    )

    assert expected_factor == pytest.approx(0.0104975811678193)
    assert result.predicted == pytest.approx(dipole_parent_br * expected_factor)
    assert result.diagnostics["dipole_conversion_factor"] == pytest.approx(
        expected_factor
    )
    assert result.diagnostics["leptonic_normalization_branching_fraction"] == pytest.approx(
        TAU_TO_E_NUNU_BRANCHING_FRACTION
    )
    assert result.diagnostics["initial_flavor"] == "tau"
    assert result.diagnostics["final_flavor"] == "e"
    assert result.diagnostics["contact_input_present"] is False
    assert result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert result.passes is True


def test_matrix_contact_is_pinned_to_tau_e_not_core_mu_e_default():
    constraint = fcc.get(_PID)
    left_overlap = [[0.0j for _ in range(3)] for _ in range(3)]
    left_overlap[0][1] = 0.09 + 0.01j
    left_overlap[0][2] = 0.003 + 0.004j
    lepton = {
        "left_charged_lepton_overlap": left_overlap,
        "m_kk_gev": 3000.0,
        "source": "L010 matrix pinning regression",
    }

    result = constraint.evaluate(_point_from_mapping(lepton))
    expected = _core_tau3e_prediction(
        constraint,
        {
            "left_lfv_overlap": 0.003 + 0.004j,
            "m_kk_gev": 3000.0,
        },
        m_kk_gev=3000.0,
    )
    wrong_default = lfv_three_body_from_components(
        dipole_parent_branching_fraction=0.0,
        contact_amplitudes=lfv_three_body_contact_amplitudes(
            lepton,
            initial_flavor="mu",
            final_flavor="e",
            m_kk_gev=3000.0,
            inputs=constraint.sm_inputs,
        ),
        br_limit=constraint.anchor.budget,
        inputs=constraint.sm_inputs,
    )

    assert result.predicted == pytest.approx(expected.branching_fraction)
    assert result.predicted != pytest.approx(wrong_default.branching_fraction)
    assert result.diagnostics["left_lfv_overlap"] == pytest.approx(0.003 + 0.004j)
    assert result.diagnostics["right_lfv_overlap"] == pytest.approx(0.0j)
    assert result.diagnostics["initial_flavor"] == "tau"
    assert result.diagnostics["final_flavor"] == "e"


def test_tau_e_named_alias_is_accepted_and_normalized():
    result = fcc.get(_PID).evaluate(
        _point_from_mapping(
            {
                "left_taue_overlap": 2.0e-3 + 1.0e-3j,
                "right_e_tau_overlap": 0.5e-3j,
                "m_kk_gev": 3000.0,
                "source": "L010 tau-e alias probe",
            }
        )
    )

    assert result.predicted is not None
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["left_lfv_overlap"] == pytest.approx(2.0e-3 + 1.0e-3j)
    assert result.diagnostics["right_lfv_overlap"] == pytest.approx(0.5e-3j)
    assert result.diagnostics["contact_input_normalized_to_tau_e"] is True


def test_proxy_numerics_match_independent_recomputation():
    constraint = fcc.get(_PID)
    lepton = {
        "initial_flavor": "tau",
        "final_flavor": "e",
        "m_kk_gev": 3000.0,
        "left_lfv_overlap": 3.0e-3 + 1.0e-3j,
        "right_lfv_overlap": 1.0e-3j,
        "box_ll": 1.0e-6 + 2.0e-7j,
        "box_lr": -0.5e-6j,
        "box_rl": 0.3e-6,
        "box_rr": -0.2e-6 + 0.1e-6j,
        "dipole_parent_branching_fraction": 2.0e-8,
        "source": "L010 test proxy",
    }
    result = constraint.evaluate(_point_from_mapping(lepton))
    expected = _core_tau3e_prediction(
        constraint,
        lepton,
        m_kk_gev=3000.0,
    )

    assert result.predicted == pytest.approx(expected.branching_fraction)
    assert result.ratio == pytest.approx(result.predicted / constraint.anchor.budget)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.diagnostics["dipole_parent_branching_fraction"] == pytest.approx(
        expected.dipole_parent_branching_fraction
    )
    assert result.diagnostics["leptonic_normalization_branching_fraction"] == pytest.approx(
        TAU_TO_E_NUNU_BRANCHING_FRACTION
    )
    assert result.diagnostics["dipole_component"] == pytest.approx(
        expected.dipole_component
    )
    assert result.diagnostics["z_penguin_component"] == pytest.approx(
        expected.z_penguin_component
    )
    assert result.diagnostics["box_component"] == pytest.approx(expected.box_component)
    assert result.diagnostics["contact_component"] == pytest.approx(
        expected.contact_component
    )
    assert result.diagnostics["dipole_contact_interference_component"] == pytest.approx(
        expected.dipole_contact_interference_component
    )
    assert result.diagnostics["dipole_contact_interference_component"] > 0.0
    assert (
        result.diagnostics["dipole_contact_interference_treatment"]
        == "constructive_sign_envelope_NEEDS-HUMAN-PHYSICS"
    )
    assert result.diagnostics["dipole_contact_interference_convention"] == (
        LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION
    )
    assert result.diagnostics["total_ll"] == pytest.approx(
        expected.contact_amplitudes.total_ll
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(
        _point_from_mapping(
            {
                "left_lfv_overlap": 1.0e-2 + 1.0e-3j,
                "right_lfv_overlap": 0.5e-2j,
                "box_ll": 0.2e-6,
                "m_kk_gev": 3000.0,
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
    for key in (
        "left_lfv_overlap",
        "right_lfv_overlap",
        "z_ll",
        "z_lr",
        "z_rl",
        "z_rr",
        "box_ll",
        "total_ll",
    ):
        assert isinstance(result.diagnostics[key], complex)
    for key in (
        "dipole_component",
        "dipole_contact_interference_component",
        "dipole_contact_interference_lower",
        "dipole_contact_interference_upper",
        "z_penguin_component",
        "box_component",
        "contact_component",
        "m_kk_gev",
        "matching_scale_gev",
        "dipole_conversion_factor",
    ):
        assert isinstance(result.diagnostics[key], float)
        assert math.isfinite(result.diagnostics[key])
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["contact_input_present"] is True
    assert result.diagnostics["dipole_input_present"] is False


@pytest.mark.parametrize(
    ("lepton", "expected_pass"),
    [
        ({"box_ll": 1.0e-4, "m_kk_gev": 3000.0}, True),
        ({"box_ll": 1.0e-3, "m_kk_gev": 3000.0}, False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(lepton, expected_pass: bool):
    result = fcc.get(_PID).evaluate(_point_from_mapping(lepton))

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_optional_kk_ew_mass_extra_changes_z_scaling():
    lepton = {
        "left_lfv_overlap": 0.1,
        "m_kk_gev": 3000.0,
        "source": "L010 scaling test proxy",
    }
    default_point = point_builder.make_point(lepton_mass_basis_couplings=lepton)
    heavy_point = point_builder.make_point(
        lepton_mass_basis_couplings=lepton,
        kk_ew_mass_gev=6000.0,
    )

    default_result = fcc.get(_PID).evaluate(default_point)
    heavy_result = fcc.get(_PID).evaluate(heavy_point)

    assert default_result.diagnostics["m_kk_gev"] == pytest.approx(3000.0)
    assert heavy_result.diagnostics["m_kk_gev"] == pytest.approx(6000.0)
    assert heavy_result.diagnostics["kk_ew_mass_extra_used"] is True
    assert heavy_result.diagnostics["z_penguin_component"] == pytest.approx(
        default_result.diagnostics["z_penguin_component"] / 16.0
    )
    assert heavy_result.diagnostics["contact_component"] == pytest.approx(
        default_result.diagnostics["contact_component"] / 16.0
    )


def test_evaluate_is_pure_and_deterministic():
    lepton = {
        "left_lfv_overlap": 1.0e-2 + 2.0e-3j,
        "right_lfv_overlap": 0.5e-2,
        "box_ll": 1.0e-6,
        "m_kk_gev": 3000.0,
    }
    point = point_builder.make_point(lepton_mass_basis_couplings=dict(lepton))
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.get_extra("lepton_mass_basis_couplings") == lepton
