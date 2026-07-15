"""Phase-4c LFV-leptonic rewire tests for L002/L009/L003-L005."""

from __future__ import annotations

import math

import pytest

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.physics_adapters.lfv_three_body import (
    LFV_THREE_BODY_DEFERRED_PIECES_V1,
    LFV_THREE_BODY_TREE_CONTACT_RIGOROUS_V1,
)
from flavor_catalog_constraints.physics_adapters.mu_e_conversion import (
    MU_E_CONVERSION_DEFERRED_PIECES_V1,
    MU_E_CONVERSION_VECTOR_TREE_RIGOROUS_V1,
)
from quarkConstraints.lfv_three_body import TAU_TO_MU_NUNU_BRANCHING_FRACTION
from quarkConstraints.mu_e_conversion import G_F_GEV_MINUS2, MUON_MASS_GEV
from tests.constraints.primary.charged_lepton.test_L002 import _rotation_pmns
from tests.constraints.primary.top_higgs_ew.z_lfv_rewire_helpers import (
    diagonal_rs_ew_point,
    lfv_live_rs_ew_point,
)


_LEPTON_INDEX = {"e": 0, "mu": 1, "tau": 2}
_THREE_BODY_LEPTONIC_NORMALIZATION = {
    ("mu", "e"): 1.0,
    ("tau", "mu"): TAU_TO_MU_NUNU_BRANCHING_FRACTION,
}


def _point_with_rs_ew(base_point):
    return point_builder.make_point(
        lepton_mass_basis_couplings=base_point.extras["lepton_mass_basis_couplings"],
        rs_ew_couplings=base_point.extras["rs_ew_couplings"],
        kk_ew_mass_gev=base_point.extras["kk_ew_mass_gev"],
    )


def _lfv_live_point(mkk_gev: float = 3000.0):
    diagonal = diagonal_rs_ew_point(mkk_gev)
    live = lfv_live_rs_ew_point(mkk_gev)
    return point_builder.make_point(
        lepton_mass_basis_couplings=diagonal.extras["lepton_mass_basis_couplings"],
        rs_ew_couplings=live.extras["rs_ew_couplings"],
        kk_ew_mass_gev=diagonal.extras["kk_ew_mass_gev"],
    )


def _tree_only_lfv_live_point(mkk_gev: float = 3000.0):
    live = lfv_live_rs_ew_point(mkk_gev)
    return point_builder.make_point(
        rs_ew_couplings=live.extras["rs_ew_couplings"],
        kk_ew_mass_gev=mkk_gev,
    )


def _manual_three_body_tree_component(couplings, *, initial: str, final: str) -> float:
    i = _LEPTON_INDEX[initial]
    f = _LEPTON_INDEX[final]
    z_ll = 2.0 * couplings.z_delta_g_L_e[f, i] * couplings.z_total_g_L_e[f, f]
    z_lr = 2.0 * couplings.z_delta_g_L_e[f, i] * couplings.z_total_g_R_e[f, f]
    z_rl = 2.0 * couplings.z_delta_g_R_e[f, i] * couplings.z_total_g_L_e[f, f]
    z_rr = 2.0 * couplings.z_delta_g_R_e[f, i] * couplings.z_total_g_R_e[f, f]
    width_component = float(
        2.0 * (abs(z_ll) ** 2 + abs(z_rr) ** 2) + abs(z_lr) ** 2 + abs(z_rl) ** 2
    )
    return float(width_component * _THREE_BODY_LEPTONIC_NORMALIZATION[(initial, final)])


def _manual_mu_e_vector_component(couplings, nuclear) -> float:
    norm = math.sqrt(2.0) * G_F_GEV_MINUS2

    def contact(sector: str, q_ch: str, l_ch: str) -> complex:
        return complex(couplings.contact(sector, q_ch, l_ch, 0, 0, 0, 1))

    g_lv_u = (contact("u", "L", "L") + contact("u", "R", "L")) / norm
    g_rv_u = (contact("u", "L", "R") + contact("u", "R", "R")) / norm
    g_lv_d = (contact("d", "L", "L") + contact("d", "R", "L")) / norm
    g_rv_d = (contact("d", "L", "R") + contact("d", "R", "R")) / norm
    g_lv_p = 2.0 * g_lv_u + g_lv_d
    g_lv_n = g_lv_u + 2.0 * g_lv_d
    g_rv_p = 2.0 * g_rv_u + g_rv_d
    g_rv_n = g_rv_u + 2.0 * g_rv_d
    left = g_lv_p * nuclear.V_p + g_lv_n * nuclear.V_n
    right = g_rv_p * nuclear.V_p + g_rv_n * nuclear.V_n
    prefactor = 2.0 * G_F_GEV_MINUS2**2 * MUON_MASS_GEV**5
    return float(prefactor * (abs(left) ** 2 + abs(right) ** 2) / nuclear.capture_rate_gev)


@pytest.mark.parametrize(
    ("pid", "initial", "final"),
    [("L002", "mu", "e"), ("L009", "tau", "mu")],
)
def test_three_body_tree_contact_v1_zero_and_lfv_live_nonzero_scaling(pid, initial, final):
    constraint = fcc.get(pid)
    diagonal = _point_with_rs_ew(diagonal_rs_ew_point())
    live = _lfv_live_point()
    high = _lfv_live_point(6000.0)

    diagonal_result = constraint.evaluate(diagonal)
    live_result = constraint.evaluate(live)
    high_result = constraint.evaluate(high)
    manual_live = _manual_three_body_tree_component(
        live.extras["rs_ew_couplings"],
        initial=initial,
        final=final,
    )
    manual_high = _manual_three_body_tree_component(
        high.extras["rs_ew_couplings"],
        initial=initial,
        final=final,
    )

    assert diagonal_result.diagnostics["tree_contact_rigorous"] is True
    assert diagonal_result.diagnostics["tree_contact_zero_for_diagonal_fit"] is True
    assert diagonal_result.diagnostics["z_penguin_component"] == pytest.approx(0.0)
    assert (
        diagonal_result.diagnostics["matching_assumption"]
        == LFV_THREE_BODY_TREE_CONTACT_RIGOROUS_V1
    )
    assert "overlap proxy" not in diagonal_result.diagnostics["matching_assumption"]
    assert live_result.diagnostics["z_penguin_component"] == pytest.approx(manual_live)
    assert live_result.diagnostics["z_penguin_component"] > 0.0
    assert (
        high_result.diagnostics["z_penguin_component"]
        / live_result.diagnostics["z_penguin_component"]
        == pytest.approx(manual_high / manual_live)
    )
    assert live_result.diagnostics["needs_human_physics"] == LFV_THREE_BODY_DEFERRED_PIECES_V1


@pytest.mark.parametrize(
    ("pid", "component_key", "rigorous_key"),
    [
        ("L002", "z_penguin_component", "tree_contact_rigorous"),
        ("L009", "z_penguin_component", "tree_contact_rigorous"),
        ("L003", "vector_component", "vector_tree_rigorous"),
        ("L004", "vector_component", "vector_tree_rigorous"),
        ("L005", "vector_component", "vector_tree_rigorous"),
    ],
)
def test_tree_only_rs_ew_extra_evaluates_without_lepton_extra(
    pid,
    component_key,
    rigorous_key,
):
    result = fcc.get(pid).evaluate(_tree_only_lfv_live_point())

    assert result.predicted is not None
    assert result.ratio is not None
    assert math.isfinite(result.predicted)
    assert math.isfinite(result.ratio)
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["rs_ew_couplings_extra_present"] is True
    assert result.diagnostics[rigorous_key] is True
    assert result.diagnostics[component_key] > 0.0
    assert result.diagnostics["dipole_input_present"] is False


@pytest.mark.parametrize("pid", ["L003", "L004", "L005"])
def test_mu_e_conversion_vector_v1_zero_and_lfv_live_nonzero_scaling(pid):
    constraint = fcc.get(pid)
    diagonal = _point_with_rs_ew(diagonal_rs_ew_point())
    live = _lfv_live_point()
    high = _lfv_live_point(6000.0)

    diagonal_result = constraint.evaluate(diagonal)
    live_result = constraint.evaluate(live)
    high_result = constraint.evaluate(high)
    manual_live = _manual_mu_e_vector_component(
        live.extras["rs_ew_couplings"],
        constraint.nuclear_inputs,
    )
    manual_high = _manual_mu_e_vector_component(
        high.extras["rs_ew_couplings"],
        constraint.nuclear_inputs,
    )

    assert diagonal_result.diagnostics["vector_tree_rigorous"] is True
    assert diagonal_result.diagnostics["vector_tree_zero_for_diagonal_fit"] is True
    assert diagonal_result.diagnostics["vector_component"] == pytest.approx(0.0)
    assert (
        diagonal_result.diagnostics["vector_tree_matching"]
        == MU_E_CONVERSION_VECTOR_TREE_RIGOROUS_V1
    )
    assert live_result.diagnostics["vector_component"] == pytest.approx(manual_live)
    assert live_result.diagnostics["vector_component"] > 0.0
    assert (
        high_result.diagnostics["vector_component"]
        / live_result.diagnostics["vector_component"]
        == pytest.approx(manual_high / manual_live)
    )
    assert live_result.diagnostics["needs_human_physics"] == MU_E_CONVERSION_DEFERRED_PIECES_V1


def test_absent_rs_ew_tree_contacts_degrade_to_zero_with_partial_dipoles_preserved():
    lepton = diagonal_rs_ew_point().extras["lepton_mass_basis_couplings"]
    point = point_builder.make_point(lepton_mass_basis_couplings=lepton)
    tau_dipole_point = point_builder.make_point(
        lepton_mass_basis_couplings={
            "dipole_parent_branching_fraction": 1.0e-8,
            "source": "tau dipole-only absent tree regression",
        }
    )

    l002 = fcc.get("L002").evaluate(point)
    l009 = fcc.get("L009").evaluate(tau_dipole_point)
    l003 = fcc.get("L003").evaluate(point)

    assert l002.diagnostics["tree_contact_missing_extra"] == "rs_ew_couplings"
    assert "not evaluated" in l002.diagnostics["tree_contact_matching"]
    assert l002.diagnostics["z_penguin_component"] == pytest.approx(0.0)
    assert l002.diagnostics["dipole_component"] > 0.0
    assert l009.diagnostics["tree_contact_missing_extra"] == "rs_ew_couplings"
    assert "not evaluated" in l009.diagnostics["tree_contact_matching"]
    assert l009.diagnostics["z_penguin_component"] == pytest.approx(0.0)
    assert l009.diagnostics["box_component"] == pytest.approx(0.0)
    assert l009.diagnostics["dipole_component"] > 0.0
    assert l003.diagnostics["vector_tree_missing_extra"] == "rs_ew_couplings"
    assert "not evaluated" in l003.diagnostics["vector_tree_matching"]
    assert l003.diagnostics["vector_component"] == pytest.approx(0.0)
    assert l003.diagnostics["dipole_component"] > 0.0


@pytest.mark.parametrize("pid", ["L003", "L004", "L005"])
def test_mu_e_conversion_scalar_dipole_partial_flags_survive_rigorous_vector(pid):
    live = lfv_live_rs_ew_point()
    lepton = {
        "y_n_bar": (1.0e-4, 2.0e-4, 3.0e-4),
        "pmns": _rotation_pmns(),
        "m_kk_gev": 3000.0,
        "g_ls_p": 1.5e-12,
        "g_rs_n": -0.4e-12j,
        "source": "scalar/dipole preservation with rigorous vector",
    }
    point = point_builder.make_point(
        lepton_mass_basis_couplings=lepton,
        rs_ew_couplings=live.extras["rs_ew_couplings"],
        kk_ew_mass_gev=3000.0,
    )
    result = fcc.get(pid).evaluate(point)

    assert result.diagnostics["vector_tree_rigorous"] is True
    assert result.diagnostics["vector_component"] > 0.0
    assert result.diagnostics["scalar_component"] > 0.0
    assert result.diagnostics["dipole_input_present"] is True
    assert result.diagnostics["dipole_contact_relative_phase_status"] == "NEEDS-HUMAN-PHYSICS"
    assert result.diagnostics["needs_human_physics"] == MU_E_CONVERSION_DEFERRED_PIECES_V1


def test_l002_box_and_dipole_contact_phase_partial_flags_survive_rigorous_tree():
    live = lfv_live_rs_ew_point()
    lepton = {
        "y_n_bar": (0.03, 0.07, 0.11),
        "pmns": _rotation_pmns(),
        "m_kk_gev": 3000.0,
        "box_ll": 1.0e-7 + 2.0e-8j,
        "box_lr": -0.5e-7j,
        "box_rl": 0.3e-7,
        "box_rr": -0.2e-7 + 0.1e-7j,
        "source": "box preservation with rigorous tree",
    }
    point = point_builder.make_point(
        lepton_mass_basis_couplings=lepton,
        rs_ew_couplings=live.extras["rs_ew_couplings"],
        kk_ew_mass_gev=3000.0,
    )
    result = fcc.get("L002").evaluate(point)

    assert result.diagnostics["tree_contact_rigorous"] is True
    assert result.diagnostics["z_penguin_component"] > 0.0
    assert result.diagnostics["box_component"] > 0.0
    assert result.diagnostics["dipole_input_present"] is True
    assert result.diagnostics["dipole_contact_interference_treatment"] == (
        "constructive_sign_envelope_NEEDS-HUMAN-PHYSICS"
    )
    assert result.diagnostics["needs_human_physics"] == LFV_THREE_BODY_DEFERRED_PIECES_V1
