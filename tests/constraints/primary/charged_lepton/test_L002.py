"""Production tests for L002 (mu -> 3e)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml

import flavor_catalog_constraints as fcc
from flavorConstraints.muToEGamma import (
    check_mu_to_e_gamma_raw,
    coefficient_from_br_limit,
)
from flavor_catalog_constraints import anchors, point_builder
from flavor_catalog_constraints.base import ConstraintProtocol, Severity
from flavor_catalog_constraints.physics_adapters.lfv_three_body import (
    lfv_three_body_contact_amplitudes as adapter_lfv_three_body_contact_amplitudes,
    lfv_three_body_proxy_input,
    mu_to_3e_from_lepton_input,
)
from quarkConstraints.lfv_three_body import (
    LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION,
    lfv_three_body_contact_amplitudes as core_lfv_three_body_contact_amplitudes,
    lfv_three_body_from_components,
)

_PID = "L002"
_REPO_ROOT = Path(__file__).resolve().parents[4]
_SIDECAR = _REPO_ROOT / "flavor_catalog" / "processes" / "charged_lepton" / "L002.yaml"
_L001_SIDECAR = (
    _REPO_ROOT / "flavor_catalog" / "processes" / "charged_lepton" / "L001.yaml"
)


def _yaml_pdg_block(path: Path = _SIDECAR):
    with open(path) as handle:
        return yaml.safe_load(handle)["pdg_or_equivalent"]


def _rotation_pmns() -> np.ndarray:
    theta = 0.5
    c = math.cos(theta)
    s = math.sin(theta)
    return np.asarray(
        [
            [c, s, 0.0],
            [-s, c, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=complex,
    )


def _point_from_mapping(mapping):
    return point_builder.make_point(lepton_mass_basis_couplings=dict(mapping))


def _l001_dipole_parent_br_from_core(
    *,
    y_n_bar: tuple[float, float, float] | None,
    prefactor_br: float,
    dipole_br_limit: float,
    m_kk_gev: float,
) -> tuple[float, complex]:
    if y_n_bar is not None:
        c_lfv = coefficient_from_br_limit(dipole_br_limit, prefactor=prefactor_br)
        core = check_mu_to_e_gamma_raw(
            np.asarray(y_n_bar, dtype=complex),
            _rotation_pmns(),
            m_kk_gev,
            C=c_lfv,
            reference_scale=3000.0,
        )
        return (
            float(
                prefactor_br * float(core["lhs"]) ** 2 * (3000.0 / m_kk_gev) ** 4
            ),
            complex(core["off_diagonal_12"]),
        )
    return 0.0, 0.0j


def _core_mu3e_prediction(
    constraint,
    lepton: dict,
    *,
    y_n_bar: tuple[float, float, float] | None,
    m_kk_gev: float,
):
    dipole_parent_br, off_diagonal = _l001_dipole_parent_br_from_core(
        y_n_bar=y_n_bar,
        prefactor_br=constraint.anchor.dipole_prefactor_br.value,
        dipole_br_limit=constraint.anchor.dipole_br_limit.value,
        m_kk_gev=m_kk_gev,
    )
    contact = core_lfv_three_body_contact_amplitudes(
        lepton,
        initial_flavor="mu",
        final_flavor="e",
        m_kk_gev=m_kk_gev,
        inputs=constraint.sm_inputs,
    )
    return (
        lfv_three_body_from_components(
            dipole_parent_branching_fraction=dipole_parent_br,
            contact_amplitudes=contact,
            br_limit=constraint.anchor.budget,
            inputs=constraint.sm_inputs,
        ),
        off_diagonal,
    )


def test_registration_metadata():
    constraint = fcc.get(_PID)

    assert isinstance(constraint, ConstraintProtocol)
    assert constraint.process_id == _PID
    assert constraint.severity is Severity.HARD
    assert constraint.family == "charged_lepton"
    assert constraint.observable == "BR(mu -> 3e)"


def test_anchor_matches_yaml_and_loud_fail_probe():
    constraint = fcc.get(_PID)
    pdg = _yaml_pdg_block()
    exp = pdg["primary_current_limit"]
    original = pdg["original_experiment"]
    l001 = _yaml_pdg_block(_L001_SIDECAR)

    assert constraint.anchor.experimental.value == pytest.approx(exp["value"])
    assert constraint.anchor.experimental.source_url == exp["source_url"]
    assert constraint.anchor.original_experiment.value == pytest.approx(original["value"])
    assert constraint.anchor.original_experiment.source_url == original["source_url"]
    assert constraint.anchor.budget == pytest.approx(1.0e-12)
    assert constraint.anchor.dipole_br_limit.value == pytest.approx(
        l001["primary_current_limit"]["value"]
    )
    assert constraint.anchor.dipole_prefactor_br.value == pytest.approx(
        l001["repo_default"]["prefac_br"]["value"]
    )

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
    assert result.notes.startswith("NOT EVALUATED --")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "lepton_mass_basis_couplings"
    assert "needs_human_physics" not in result.diagnostics


def test_invalid_lepton_input_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(lepton_mass_basis_couplings={"source": "empty"})
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.notes.startswith("NOT EVALUATED --")
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["invalid_extra"] == "lepton_mass_basis_couplings"
    assert result.diagnostics["exception_type"] == "TypeError"


def test_mu_to_3e_rejects_caller_flavor_override():
    constraint = fcc.get(_PID)
    bad = {
        "initial_flavor": "tau",
        "final_flavor": "mu",
        "left_lfv_overlap": 2.0e-4,
        "m_kk_gev": 3000.0,
        "source": "bad tau->mu probe",
    }

    with pytest.raises(ValueError, match="pinned"):
        mu_to_3e_from_lepton_input(
            bad,
            br_limit=constraint.anchor.budget,
            dipole_br_limit=constraint.anchor.dipole_br_limit.value,
            dipole_prefactor_br=constraint.anchor.dipole_prefactor_br.value,
            inputs=constraint.sm_inputs,
        )

    result = constraint.evaluate(_point_from_mapping(bad))
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["exception_type"] == "ValueError"
    assert "pinned to initial_flavor='mu'" in result.diagnostics["exception"]


def test_absent_tree_ignores_legacy_overlap_but_preserves_dipole_and_box():
    constraint = fcc.get(_PID)
    y_n_bar = (0.03, 0.07, 0.11)
    lepton = {
        "y_n_bar": y_n_bar,
        "pmns": _rotation_pmns(),
        "m_kk_gev": 3000.0,
        "left_emu_overlap": 2.0e-4 + 1.0e-4j,
        "right_emu_overlap": 1.0e-4j,
        "box_ll": 1.0e-7 + 2.0e-8j,
        "box_lr": -0.5e-7j,
        "box_rl": 0.3e-7,
        "box_rr": -0.2e-7 + 0.1e-7j,
        "source": "L002 test proxy",
    }
    result = constraint.evaluate(_point_from_mapping(lepton))
    expected_dipole_br, expected_off_diagonal = _l001_dipole_parent_br_from_core(
        y_n_bar=y_n_bar,
        prefactor_br=constraint.anchor.dipole_prefactor_br.value,
        dipole_br_limit=constraint.anchor.dipole_br_limit.value,
        m_kk_gev=3000.0,
    )
    expected_contact = adapter_lfv_three_body_contact_amplitudes(
        lepton,
        initial_flavor="mu",
        final_flavor="e",
        m_kk_gev=3000.0,
        inputs=constraint.sm_inputs,
    )
    expected = lfv_three_body_from_components(
        dipole_parent_branching_fraction=expected_dipole_br,
        contact_amplitudes=expected_contact,
        br_limit=constraint.anchor.budget,
        inputs=constraint.sm_inputs,
    )

    assert result.predicted == pytest.approx(expected.branching_fraction)
    assert result.ratio == pytest.approx(result.predicted / constraint.anchor.budget)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.diagnostics["dipole_parent_branching_fraction"] == pytest.approx(
        expected_dipole_br
    )
    assert result.diagnostics["dipole_component"] == pytest.approx(
        expected.dipole_component
    )
    assert result.diagnostics["dipole_conversion_factor"] == pytest.approx(
        expected.dipole_conversion_factor
    )
    assert result.diagnostics["z_penguin_component"] == pytest.approx(
        0.0
    )
    assert result.diagnostics["tree_contact_missing_extra"] == "rs_ew_couplings"
    assert result.diagnostics["legacy_overlap_tree_proxy_ignored"] is True
    assert result.diagnostics["box_component"] == pytest.approx(
        expected.box_component
    )
    assert result.diagnostics["contact_component"] == pytest.approx(
        expected.contact_component
    )
    assert result.diagnostics["z_box_interference_component"] == pytest.approx(
        expected.z_box_interference_component
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
    assert result.diagnostics["dipole_off_diagonal_12"] == pytest.approx(
        expected_off_diagonal
    )
    assert result.diagnostics["total_ll"] == pytest.approx(
        expected.contact_amplitudes.total_ll
    )
    assert "NEEDS-HUMAN-PHYSICS" in result.diagnostics["needs_human_physics"]


def test_core_dipole_contact_interference_uses_chiral_structure():
    inputs = fcc.get(_PID).sm_inputs
    contact = core_lfv_three_body_contact_amplitudes(
        {
            "box_ll": 3.0e-7 + 1.0e-7j,
            "box_lr": -2.0e-7j,
            "box_rl": 0.5e-7 + 0.25e-7j,
            "box_rr": -1.5e-7,
            "m_kk_gev": 3000.0,
            "source": "explicit chiral interference test",
        },
        initial_flavor="mu",
        final_flavor="e",
        inputs=inputs,
    )
    a_left = 1.2e-8 + 0.4e-8j
    a_right = -0.7e-8 + 0.2e-8j
    dipole_parent_br = 384.0 * math.pi**2 * (abs(a_left) ** 2 + abs(a_right) ** 2)

    result = lfv_three_body_from_components(
        dipole_parent_branching_fraction=dipole_parent_br,
        contact_amplitudes=contact,
        inputs=inputs,
        dipole_amplitude_left=a_left,
        dipole_amplitude_right=a_right,
    )

    electric_charge = math.sqrt(4.0 * math.pi * inputs.alpha_em)
    left_combo = 2.0 * contact.total_rr + contact.total_rl
    right_combo = 2.0 * contact.total_ll + contact.total_lr
    # Kuno-Okada Eq. (2.14) fixes the coefficient/sign to -8e in this convention.
    expected_interference = float(
        -8.0
        * electric_charge
        * (a_right * right_combo.conjugate() + a_left * left_combo.conjugate()).real
    )

    assert result.dipole_contact_interference_component == pytest.approx(
        expected_interference,
        abs=0.0,
    )
    assert result.dipole_contact_interference_treatment == (
        "explicit_chiral_dipole_amplitudes"
    )
    assert result.diagnostics["dipole_contact_chiral_combination_left"] == pytest.approx(
        left_combo
    )
    assert result.diagnostics["dipole_contact_chiral_combination_right"] == pytest.approx(
        right_combo
    )


def test_evaluate_runs_end_to_end_with_real_finite_fields_and_complex_diagnostics():
    result = fcc.get(_PID).evaluate(
        _point_from_mapping(
            {
                "left_emu_overlap": 1.0e-4 + 1.0e-4j,
                "right_emu_overlap": 0.5e-4j,
                "box_ll": 0.2e-7,
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
    assert result.diagnostics["dipole_input_present"] is False
    assert result.diagnostics["tree_contact_missing_extra"] == "rs_ew_couplings"
    assert result.diagnostics["legacy_overlap_tree_proxy_ignored"] is True
    assert result.diagnostics["z_penguin_component"] == pytest.approx(0.0)


@pytest.mark.parametrize(
    ("lepton", "expected_pass"),
    [
        ({"box_ll": 1.0e-7, "m_kk_gev": 3000.0}, True),
        ({"box_ll": 1.0e-5, "m_kk_gev": 3000.0}, False),
    ],
)
def test_safe_point_passes_and_large_np_point_fails(lepton, expected_pass: bool):
    result = fcc.get(_PID).evaluate(_point_from_mapping(lepton))

    assert result.passes is expected_pass
    if expected_pass:
        assert result.ratio < 1.0
    else:
        assert result.ratio > 10.0


def test_legacy_overlap_only_proxy_is_unevaluated_not_real_pass():
    result = fcc.get(_PID).evaluate(
        _point_from_mapping({"left_emu_overlap": 1.0e-4, "m_kk_gev": 3000.0})
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["exception_type"] == "TypeError"


def test_typed_legacy_overlap_only_proxy_is_unevaluated_not_real_pass():
    proxy = lfv_three_body_proxy_input(
        1.0e-4,
        0.0j,
        3000.0,
        source="typed stale L002 overlap proxy",
    )
    result = fcc.get(_PID).evaluate(
        point_builder.make_point(lepton_mass_basis_couplings=proxy)
    )

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["exception_type"] == "TypeError"


def test_optional_kk_ew_mass_extra_changes_dipole_while_absent_tree_stays_zero():
    lepton = {
        "y_n_bar": (0.03, 0.07, 0.11),
        "pmns": _rotation_pmns(),
        "m_kk_gev": 3000.0,
        "left_emu_overlap": 2.0e-4,
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
    assert heavy_result.diagnostics["dipole_component"] == pytest.approx(
        default_result.diagnostics["dipole_component"] / 16.0
    )
    assert heavy_result.diagnostics["z_penguin_component"] == pytest.approx(
        0.0
    )
    assert default_result.diagnostics["z_penguin_component"] == pytest.approx(0.0)
    assert default_result.diagnostics["legacy_overlap_tree_proxy_ignored"] is True
    assert heavy_result.diagnostics["legacy_overlap_tree_proxy_ignored"] is True
    assert heavy_result.diagnostics["dipole_contact_interference_component"] == pytest.approx(
        default_result.diagnostics["dipole_contact_interference_component"] / 16.0
    )


def test_evaluate_is_pure_and_deterministic():
    lepton = {
        "left_emu_overlap": 1.0e-4 + 2.0e-5j,
        "right_emu_overlap": 0.5e-4,
        "box_ll": 1.0e-7,
        "m_kk_gev": 3000.0,
    }
    point = point_builder.make_point(lepton_mass_basis_couplings=dict(lepton))
    constraint = fcc.get(_PID)

    first = constraint.evaluate(point)
    second = constraint.evaluate(point)

    assert first == second
    assert point.get_extra("lepton_mass_basis_couplings") == lepton
