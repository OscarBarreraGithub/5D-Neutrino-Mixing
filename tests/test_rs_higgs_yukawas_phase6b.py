import math
from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

import flavor_catalog_constraints as fcc
from flavor_catalog_constraints import point_builder
from quarkConstraints.higgs_lfv import (
    h_lfv_branching_fraction_from_yukawas,
    h_lfv_branching_fraction_with_proxy,
)
from quarkConstraints.rs_higgs_yukawas import (
    RS_HIGGS_YUKAWA_MATCHING_ASSUMPTION_V1,
    RSHiggsYukawaCouplings,
    build_rs_higgs_yukawas,
)
from tests.constraints.primary.top_higgs_ew.higgs_lfv_phase6b_helpers import (
    diagonal_higgs_point,
    live_lepton_couplings,
    live_higgs_yukawas,
)


def _offdiag(matrix):
    arr = np.array(matrix, dtype=np.complex128, copy=True)
    np.fill_diagonal(arr, 0.0)
    return arr


def test_builder_emits_declared_immutable_api_and_diagonal_v1_exact_zero():
    point = diagonal_higgs_point()
    higgs = point.extras["rs_higgs_yukawas"]

    assert "rs_higgs_yukawas" in point_builder.KNOWN_EXTRA_KEYS
    assert isinstance(higgs, RSHiggsYukawaCouplings)
    assert higgs.units == "dimensionless"
    assert higgs.matching_assumption == RS_HIGGS_YUKAWA_MATCHING_ASSUMPTION_V1
    assert higgs.includes_fermion_kk_mixing is True
    assert higgs.Y_h_mass is higgs.higgs_yukawa_matrix
    assert np.array_equal(_offdiag(higgs.delta_L), np.zeros((3, 3), dtype=np.complex128))
    assert np.array_equal(_offdiag(higgs.delta_E), np.zeros((3, 3), dtype=np.complex128))
    assert np.array_equal(
        _offdiag(higgs.higgs_yukawa_matrix),
        np.zeros((3, 3), dtype=np.complex128),
    )
    assert higgs.diagnostics["diagonal_v1_tree_lfv_zero"] is True
    assert higgs.diagnostics["higgs_lfv_offdiag_max_abs"] == pytest.approx(0.0)

    for arr in (
        higgs.charged_lepton_masses_gev,
        higgs.c_L,
        higgs.c_E,
        higgs.F_L,
        higgs.F_E,
        higgs.profile_B_L,
        higgs.profile_B_E,
        higgs.x_e,
        higgs.delta_L,
        higgs.delta_E,
        higgs.higgs_yukawa_matrix,
        higgs.Y_E_bar_matrix,
        higgs.U_e_L,
        higgs.U_e_R,
    ):
        assert np.all(np.isfinite(arr))
        assert arr.flags.writeable is False
    with pytest.raises(ValueError):
        higgs.higgs_yukawa_matrix[0, 1] = 1.0
    with pytest.raises(TypeError):
        higgs.diagnostics["new"] = "blocked"
    with pytest.raises(FrozenInstanceError):
        higgs.units = "GeV"


@pytest.mark.parametrize("pid", ["T018", "T019", "T020"])
def test_diagonal_v1_t018_t019_t020_evaluate_rigorous_zero(pid):
    result = fcc.get(pid).evaluate(
        point_builder.make_point(
            rs_higgs_yukawas=diagonal_higgs_point().extras["rs_higgs_yukawas"]
        )
    )

    assert result.passes is True
    assert result.predicted == pytest.approx(0.0)
    assert result.ratio == pytest.approx(0.0)
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.diagnostics["evaluated"] is True
    assert result.diagnostics["rs_higgs_yukawas_units"] == "dimensionless"
    assert result.diagnostics["includes_fermion_kk_mixing"] is True
    assert result.diagnostics["rs_higgs_yukawas_diagnostics"][
        "diagonal_v1_tree_lfv_zero"
    ] is True
    assert "NEEDS-HUMAN-PHYSICS" not in result.diagnostics["rs_matching_assumption"]


def test_lfv_live_toy_nonzero_matches_core_and_scales_as_inverse_lambda_ir_squared():
    higgs_3tev = live_higgs_yukawas(3000.0)
    higgs_6tev = live_higgs_yukawas(6000.0)
    adapter_result, proxy = h_lfv_branching_fraction_with_proxy(
        higgs_3tev,
        initial_flavor="mu",
        final_flavor="tau",
        br_limit=1.0e-3,
    )
    core_result = h_lfv_branching_fraction_from_yukawas(
        initial_flavor="mu",
        final_flavor="tau",
        yukawa_ij=higgs_3tev.higgs_yukawa_matrix[1, 2],
        yukawa_ji=higgs_3tev.higgs_yukawa_matrix[2, 1],
        br_limit=1.0e-3,
    )

    assert higgs_3tev.diagnostics["charged_lepton_yukawa_is_diagonal"] is False
    assert higgs_3tev.diagnostics["higgs_lfv_offdiag_exact_zero"] is False
    assert abs(higgs_3tev.higgs_yukawa_matrix[1, 2]) > 1.0e-5
    assert adapter_result.branching_fraction == pytest.approx(
        core_result.branching_fraction
    )
    assert adapter_result.partial_width_gev == pytest.approx(core_result.partial_width_gev)
    assert proxy.matching_assumption == higgs_3tev.matching_assumption

    y_ratio = higgs_6tev.higgs_yukawa_matrix[1, 2] / higgs_3tev.higgs_yukawa_matrix[1, 2]
    norm_3tev = abs(higgs_3tev.higgs_yukawa_matrix[1, 2]) ** 2 + abs(
        higgs_3tev.higgs_yukawa_matrix[2, 1]
    ) ** 2
    norm_6tev = abs(higgs_6tev.higgs_yukawa_matrix[1, 2]) ** 2 + abs(
        higgs_6tev.higgs_yukawa_matrix[2, 1]
    ) ** 2
    assert y_ratio == pytest.approx(0.25 + 0.0j, rel=1.0e-14, abs=1.0e-18)
    assert norm_6tev / norm_3tev == pytest.approx(0.0625, rel=1.0e-14)


def test_alignment_and_mass_to_zero_sanity():
    diagonal = diagonal_higgs_point()
    lepton = diagonal.extras["lepton_mass_basis_couplings"]
    spectrum = diagonal.extras["rs_ew_spectrum"]
    aligned = replace(
        lepton,
        c_L=np.array([0.51, 0.65, 0.74], dtype=float),
        c_E=np.array([0.68, 0.57, 0.48], dtype=float),
    )
    aligned_higgs = build_rs_higgs_yukawas(aligned, spectrum=spectrum)
    zero_mu_higgs = build_rs_higgs_yukawas(
        live_lepton_couplings(),
        spectrum=spectrum,
        charged_lepton_masses_gev=(0.51099895e-3, 0.0, 1.77686),
    )

    assert np.array_equal(_offdiag(aligned_higgs.higgs_yukawa_matrix), np.zeros((3, 3)))
    assert np.array_equal(zero_mu_higgs.higgs_yukawa_matrix[1, :], np.zeros(3))
    assert np.array_equal(zero_mu_higgs.higgs_yukawa_matrix[:, 1], np.zeros(3))


@pytest.mark.parametrize("pid", ["T018", "T019", "T020"])
def test_absent_rs_higgs_yukawas_path_is_non_vetoing(pid):
    result = fcc.get(pid).evaluate(point_builder.empty_point())

    assert result.passes is True
    assert result.predicted is None
    assert result.ratio is None
    assert result.sm_prediction == pytest.approx(0.0)
    assert result.diagnostics["evaluated"] is False
    assert result.diagnostics["missing_extra"] == "rs_higgs_yukawas"


def test_higgs_yukawa_matching_fails_loudly_on_missing_point_v():
    point = diagonal_higgs_point()
    lepton = point.extras["lepton_mass_basis_couplings"]
    bad = replace(lepton, params={k: v for k, v in lepton.params.items() if k != "v"})

    with pytest.raises(ValueError, match="params\\['v'\\]"):
        build_rs_higgs_yukawas(bad, spectrum=point.extras["rs_ew_spectrum"])
