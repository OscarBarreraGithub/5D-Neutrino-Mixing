"""Regression tests for the rare-dilepton Z-prime proxy matching."""

from __future__ import annotations

import math

import numpy as np
import pytest

from quarkConstraints import rare_b_dilepton, rare_charm_dilepton
from quarkConstraints.couplings import QuarkMassBasisCouplings


def _couplings(
    *,
    left_down_12: complex = 0.0j,
    right_down_12: complex = 0.0j,
    left_up_01: complex = 0.0j,
    right_up_01: complex = 0.0j,
    m_kk_gev: float = 3000.0,
) -> QuarkMassBasisCouplings:
    zeros = np.zeros((3, 3), dtype=np.complex128)
    left_down = zeros.copy()
    right_down = zeros.copy()
    left_up = zeros.copy()
    right_up = zeros.copy()
    left_down[1, 2] = left_down_12
    right_down[1, 2] = right_down_12
    left_up[0, 1] = left_up_01
    right_up[0, 1] = right_up_01
    return QuarkMassBasisCouplings(
        M_KK=m_kk_gev,
        xi_KK=1.0,
        alpha_s=0.09,
        g_s=1.0,
        left_overlap=zeros,
        right_up_overlap=zeros,
        right_down_overlap=zeros,
        left_up=left_up,
        left_down=left_down,
        right_up=right_up,
        right_down=right_down,
    )


def _z_couplings(alpha_em_mz: float, sin2_theta_w: float) -> tuple[float, float, float]:
    g_weak = math.sqrt(4.0 * math.pi * alpha_em_mz / sin2_theta_w)
    g_z = g_weak / math.sqrt(1.0 - sin2_theta_w)
    lepton_left = g_z * (-0.5 + sin2_theta_w)
    lepton_right = g_z * sin2_theta_w
    return g_z / 2.0, lepton_left + lepton_right, lepton_right - lepton_left


def test_rare_b_zprime_proxy_matching_uses_wet_sign_and_half_factor():
    inputs = rare_b_dilepton.default_sm_inputs()
    source = _couplings(left_down_12=2.0e-3 + 1.0e-4j, right_down_12=-8.0e-4j)

    wilsons = rare_b_dilepton.compute_rare_b_dilepton_wilsons(
        source,
        transition="b_s",
        inputs=inputs,
    )

    lambda_t = rare_b_dilepton.ckm_factors("b_s", inputs).lambda_t
    neutral_delta, mu_vector, mu_axial = _z_couplings(
        inputs.alpha_em_mz,
        inputs.sin2_theta_w,
    )
    prefactor = math.pi / (
        2.0
        * math.sqrt(2.0)
        * inputs.gf_gev_minus2
        * inputs.alpha_em_mz
        * lambda_t
        * source.M_KK**2
    )
    old_prefactor = -math.pi / (
        math.sqrt(2.0)
        * inputs.gf_gev_minus2
        * inputs.alpha_em_mz
        * lambda_t
        * source.M_KK**2
    )

    expected_c9 = prefactor * neutral_delta * source.left_down[1, 2] * mu_vector
    expected_c10 = prefactor * neutral_delta * source.left_down[1, 2] * mu_axial
    expected_c9p = prefactor * neutral_delta * source.right_down[1, 2] * mu_vector
    assert wilsons.c9_np == pytest.approx(expected_c9)
    assert wilsons.c10_np == pytest.approx(expected_c10)
    assert wilsons.c9p_np == pytest.approx(expected_c9p)
    assert wilsons.c9_np == pytest.approx(
        -0.5 * old_prefactor * neutral_delta * source.left_down[1, 2] * mu_vector
    )


def test_rare_charm_zprime_proxy_matching_uses_wet_sign_and_half_factor():
    inputs = rare_charm_dilepton.default_sm_inputs()
    source = _couplings(left_up_01=1.5e-3 - 2.0e-4j, right_up_01=7.0e-4j)

    wilsons = rare_charm_dilepton.compute_rare_charm_dilepton_wilsons(
        source,
        transition="c_u",
        inputs=inputs,
    )

    lambda_b = rare_charm_dilepton.ckm_factors("c_u", inputs).lambda_b
    neutral_delta, lep_vector, lep_axial = _z_couplings(
        inputs.alpha_em_mz,
        inputs.sin2_theta_w,
    )
    prefactor = math.pi / (
        2.0
        * math.sqrt(2.0)
        * inputs.gf_gev_minus2
        * inputs.alpha_em_mz
        * lambda_b
        * source.M_KK**2
    )
    old_prefactor = -math.pi / (
        math.sqrt(2.0)
        * inputs.gf_gev_minus2
        * inputs.alpha_em_mz
        * lambda_b
        * source.M_KK**2
    )

    expected_c9 = prefactor * neutral_delta * source.left_up[0, 1] * lep_vector
    expected_c10 = prefactor * neutral_delta * source.left_up[0, 1] * lep_axial
    expected_c10p = prefactor * neutral_delta * source.right_up[0, 1] * lep_axial
    assert wilsons.c9_np == pytest.approx(expected_c9)
    assert wilsons.c10_np == pytest.approx(expected_c10)
    assert wilsons.c10p_np == pytest.approx(expected_c10p)
    assert wilsons.c10_np == pytest.approx(
        -0.5 * old_prefactor * neutral_delta * source.left_up[0, 1] * lep_axial
    )
