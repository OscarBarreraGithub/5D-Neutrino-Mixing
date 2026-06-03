"""Shared Phase-4c LFV rare-K/charm rewire test helpers."""

from __future__ import annotations

from dataclasses import replace
import math
from functools import lru_cache

import numpy as np

from flavor_catalog_constraints import point_builder
from quarkConstraints.rare_kaon_dilepton import ckm_factors as kaon_ckm_factors
from quarkConstraints.rare_kaon_dilepton import g_sm_squared as kaon_g_sm_squared
from quarkConstraints.rs_ew_couplings import build_rs_ew_couplings
from quarkConstraints.rs_semileptonic_wilsons import build_rs_semileptonic_wilsons
from quarkConstraints.rare_charm_semileptonic import dtopi_fplus, dtopi_fzero
from tests.constraints.primary.top_higgs_ew.z_lfv_rewire_helpers import (
    _rot12,
    _rot13,
    _rot23,
    diagonal_rs_ew_point,
)
from tests.rs_ew_phase3b_helpers import (
    MAX_OVERLAP_MODES,
    MIN_OVERLAP_MODES,
    OVERLAP_REL_TOL,
    _sample_fit,
)

_METRIC = np.diag([1.0, -1.0, -1.0, -1.0])


@lru_cache(maxsize=None)
def diagonal_lfv_rare_point(mkk_gev: float = 3000.0):
    return diagonal_rs_ew_point(mkk_gev)


@lru_cache(maxsize=None)
def lfv_live_rare_point(mkk_gev: float = 3000.0):
    base = diagonal_lfv_rare_point(mkk_gev)
    lepton = replace(
        base.extras["lepton_mass_basis_couplings"],
        c_L=np.array([0.51, 0.65, 0.74], dtype=float),
        c_E=np.array([0.68, 0.57, 0.48], dtype=float),
        U_e_L=_rot12(0.28) @ _rot23(-0.22) @ _rot13(0.18, 0.35),
        U_e_R=_rot12(-0.19) @ _rot23(0.24) @ _rot13(-0.16, -0.45),
    )
    couplings = build_rs_ew_couplings(
        _sample_fit(),
        spectrum=base.extras["rs_ew_spectrum"],
        lepton_mass_basis_couplings=lepton,
        overlap_rel_tol=OVERLAP_REL_TOL,
        min_overlap_modes=MIN_OVERLAP_MODES,
        max_overlap_modes=MAX_OVERLAP_MODES,
    )
    return point_builder.make_point(
        kk_ew_mass_gev=base.extras["kk_ew_mass_gev"],
        rs_ew_couplings=couplings,
        rs_semileptonic_wilsons=build_rs_semileptonic_wilsons(couplings),
    )


def lfv_coeff(point, transition: str, lepton_pair: str = "e_mu"):
    return point.extras["rs_semileptonic_wilsons"].lfv_llqq[transition][lepton_pair]


def scaled_lfv_rare_point(point, transition: str, scale: float, lepton_pair: str = "e_mu"):
    bundle = point.extras["rs_semileptonic_wilsons"]
    blocks = {name: dict(block) for name, block in bundle.lfv_llqq.items()}
    coeff = blocks[transition][lepton_pair]
    blocks[transition][lepton_pair] = replace(
        coeff,
        contact_LL=coeff.contact_LL * scale,
        contact_LR=coeff.contact_LR * scale,
        contact_RL=coeff.contact_RL * scale,
        contact_RR=coeff.contact_RR * scale,
        c9_lfv_np=coeff.c9_lfv_np * scale,
        c10_lfv_np=coeff.c10_lfv_np * scale,
        c9p_lfv_np=coeff.c9p_lfv_np * scale,
        c10p_lfv_np=coeff.c10p_lfv_np * scale,
    )
    extras = dict(point.extras)
    extras["rs_semileptonic_wilsons"] = replace(bundle, lfv_llqq=blocks)
    return point_builder.make_point(**extras)


def kaon_lfv_y_two_body(coeff, inputs):
    norm = _kaon_y_norm(inputs.rare_kaon)
    lam = complex(coeff.lambda_ckm)
    return (
        -lam * norm * complex(coeff.c9_lfv_np - coeff.c9p_lfv_np),
        -lam * norm * complex(coeff.c10_lfv_np - coeff.c10p_lfv_np),
    )


def kaon_lfv_y_semileptonic(coeff, inputs):
    norm = _kaon_y_norm(inputs.rare_kaon)
    lam = complex(coeff.lambda_ckm)
    return (
        -lam * norm * complex(coeff.c9_lfv_np + coeff.c9p_lfv_np),
        -lam * norm * complex(coeff.c10_lfv_np + coeff.c10p_lfv_np),
    )


def manual_klong_emu_rate(coeff, inputs):
    y_vector, y_axial = kaon_lfv_y_two_body(coeff, inputs)
    rare_inputs = inputs.rare_kaon
    lam_wolf = kaon_ckm_factors(rare_inputs).lambda_wolfenstein
    kappa_mu = rare_inputs.kappa_mu_ref * (
        lam_wolf / rare_inputs.kappa_lambda_ref
    ) ** 8
    r_e = inputs.electron_mass_gev / inputs.kaon_mass_gev
    r_mu = inputs.muon_mass_gev / inputs.kaon_mass_gev
    phase_lambda = _kallen_dimensionless(r_e * r_e, r_mu * r_mu)
    beta_mumu = math.sqrt(1.0 - 4.0 * r_mu * r_mu)
    mass_sum = r_e + r_mu
    mass_diff = r_mu - r_e
    vector_term = (1.0 - mass_sum * mass_sum) * abs(mass_diff * y_vector) ** 2
    axial_term = (1.0 - mass_diff * mass_diff) * abs(mass_sum * y_axial) ** 2
    return float(
        inputs.charge_state_factor
        * kappa_mu
        * (math.sqrt(phase_lambda) / beta_mumu)
        * (vector_term + axial_term)
        / (4.0 * r_mu * r_mu * lam_wolf**10)
    )


def manual_kaon_ktopi_rate(coeff, inputs, *, charge_mode: str):
    y_vector, y_axial = kaon_lfv_y_semileptonic(coeff, inputs)
    q2_min = (inputs.electron_mass_gev + inputs.muon_mass_gev) ** 2
    q2_max = (inputs.kplus_mass_gev - inputs.piplus_mass_gev) ** 2
    nodes, weights = np.polynomial.legendre.leggauss(inputs.quadrature_points)
    mid = 0.5 * (q2_max + q2_min)
    half_width = 0.5 * (q2_max - q2_min)
    total = 0.0
    for node, weight in zip(nodes, weights):
        total += float(weight) * _manual_kaon_differential(
            mid + half_width * float(node),
            y_vector=y_vector,
            y_axial=y_axial,
            inputs=inputs,
            charge_mode=charge_mode,
        )
    return float(half_width * total)


def manual_d0_emu_rate(coeff, inputs, *, charge_state_factor: float = 2.0):
    c9_combo = complex(coeff.c9_lfv_np - coeff.c9p_lfv_np)
    c10_combo = complex(coeff.c10_lfv_np - coeff.c10p_lfv_np)
    meson = inputs.d0
    electron = inputs.electron
    muon = inputs.muon
    r_e = electron.mass_gev / meson.meson_mass_gev
    r_mu = muon.mass_gev / meson.meson_mass_gev
    phase_lambda = _kallen_dimensionless(r_e * r_e, r_mu * r_mu)
    mass_sum = r_e + r_mu
    mass_diff = r_mu - r_e
    vector_term = (1.0 - mass_sum * mass_sum) * abs(mass_diff * c9_combo) ** 2
    axial_term = (1.0 - mass_diff * mass_diff) * abs(mass_sum * c10_combo) ** 2
    tau = meson.lifetime_ps * 1.0e-12 / inputs.hbar_gev_s
    return float(
        charge_state_factor
        * tau
        * inputs.gf_gev_minus2**2
        * inputs.alpha_em_mz**2
        / (64.0 * math.pi**3)
        * meson.decay_constant_gev**2
        * meson.meson_mass_gev**3
        * abs(coeff.lambda_ckm) ** 2
        * math.sqrt(phase_lambda)
        * (vector_term + axial_term)
    )


def manual_dtopi_emu_rate(coeff, inputs, *, charge_mode: str = "eplus_muminus"):
    c9_semileptonic = complex(coeff.c9_lfv_np + coeff.c9p_lfv_np)
    c10_semileptonic = complex(coeff.c10_lfv_np + coeff.c10p_lfv_np)
    m_d = inputs.dplus_mass_gev
    m_pi = inputs.piplus_mass_gev
    m_e = inputs.rare_charm.electron.mass_gev
    m_mu = inputs.rare_charm.muon.mass_gev
    if charge_mode == "eplus_muminus":
        m_minus, m_plus = m_mu, m_e
    elif charge_mode == "eminus_muplus":
        m_minus, m_plus = m_e, m_mu
    else:
        raise ValueError(charge_mode)
    q2_min = (m_e + m_mu) ** 2
    q2_max = (m_d - m_pi) ** 2
    tau = inputs.dplus_lifetime_ps * 1.0e-12 / inputs.rare_charm.hbar_gev_s
    gamma0 = (
        inputs.rare_charm.gf_gev_minus2**2
        * inputs.rare_charm.alpha_em_mz**2
        * abs(coeff.lambda_ckm) ** 2
        / ((4.0 * math.pi) ** 5 * m_d**3)
    )
    nodes, weights = np.polynomial.legendre.leggauss(inputs.quadrature_points)
    mid = 0.5 * (q2_max + q2_min)
    half_width = 0.5 * (q2_max - q2_min)
    total = 0.0
    for node, weight in zip(nodes, weights):
        q2 = mid + half_width * float(node)
        lam_had = _kallen(m_d * m_d, m_pi * m_pi, q2)
        lam_lep = _kallen(q2, m_minus * m_minus, m_plus * m_plus)
        beta_lfv = math.sqrt(max(0.0, lam_lep)) / q2
        angular = _manual_dtopi_angular(
            q2,
            m_minus=m_minus,
            m_plus=m_plus,
            c9_lfv_semileptonic=c9_semileptonic,
            c10_lfv_semileptonic=c10_semileptonic,
            inputs=inputs,
        )
        total += (
            float(weight)
            * tau
            * gamma0
            * math.sqrt(max(0.0, lam_had))
            * beta_lfv
            * angular
            / 4.0
        )
    return float(half_width * total)


def _manual_kaon_differential(q2, *, y_vector, y_axial, inputs, charge_mode):
    if charge_mode == "muplus_eminus":
        m_minus = inputs.electron_mass_gev
        m_plus = inputs.muon_mass_gev
    elif charge_mode == "muminus_eplus":
        m_minus = inputs.muon_mass_gev
        m_plus = inputs.electron_mass_gev
    else:
        raise ValueError(charge_mode)
    m_k = inputs.kplus_mass_gev
    m_pi = inputs.piplus_mass_gev
    lam_had = max(0.0, _kallen(m_k * m_k, m_pi * m_pi, q2))
    lam_lep = max(0.0, _kallen(q2, m_minus * m_minus, m_plus * m_plus))
    if lam_had <= 0.0 or lam_lep <= 0.0:
        return 0.0
    angular = _manual_kaon_angular(
        q2,
        m_minus=m_minus,
        m_plus=m_plus,
        y_vector=y_vector,
        y_axial=y_axial,
        inputs=inputs,
    )
    rare = inputs.rare_kaon
    c0 = rare.gf_gev_minus2 / math.sqrt(2.0) * rare.alpha_em_mz / (
        2.0 * math.pi * rare.sin2_theta_w
    )
    tau = inputs.kplus_lifetime_s / inputs.hbar_gev_s
    phase_space = math.sqrt(lam_had) * math.sqrt(lam_lep) / (
        512.0 * math.pi**3 * m_k**3 * q2
    )
    return float(inputs.charge_state_factor * tau * c0 * c0 * phase_space * angular)


def _manual_kaon_angular(q2, *, m_minus, m_plus, y_vector, y_axial, inputs):
    nodes, weights = np.polynomial.legendre.leggauss(48)
    total = 0.0
    for cos_theta, weight in zip(nodes, weights):
        vector, axial = _kaon_tensor_contractions(
            q2,
            float(cos_theta),
            m_minus=m_minus,
            m_plus=m_plus,
            inputs=inputs,
        )
        total += float(weight) * (
            abs(y_vector) ** 2 * vector + abs(y_axial) ** 2 * axial
        )
    return float(total)


def _manual_dtopi_angular(
    q2,
    *,
    m_minus,
    m_plus,
    c9_lfv_semileptonic,
    c10_lfv_semileptonic,
    inputs,
):
    nodes, weights = np.polynomial.legendre.leggauss(48)
    total = 0.0
    for cos_theta, weight in zip(nodes, weights):
        vector, axial = _dtopi_tensor_contractions(
            q2,
            float(cos_theta),
            m_minus=m_minus,
            m_plus=m_plus,
            inputs=inputs,
        )
        total += float(weight) * (
            abs(c9_lfv_semileptonic) ** 2 * vector
            + abs(c10_lfv_semileptonic) ** 2 * axial
        )
    return float(total)


def _kaon_tensor_contractions(q2, cos_theta, *, m_minus, m_plus, inputs):
    root_q2 = math.sqrt(q2)
    m_k = inputs.kplus_mass_gev
    m_pi = inputs.piplus_mass_gev
    m_k2 = m_k * m_k
    m_pi2 = m_pi * m_pi
    delta = m_k2 - m_pi2
    p_had = math.sqrt(max(0.0, _kallen(m_k2, m_pi2, q2))) / (2.0 * root_q2)
    p_lep = math.sqrt(max(0.0, _kallen(q2, m_minus * m_minus, m_plus * m_plus))) / (
        2.0 * root_q2
    )
    e_k = (m_k2 - m_pi2 + q2) / (2.0 * root_q2)
    e_pi = (m_k2 - m_pi2 - q2) / (2.0 * root_q2)
    e_minus = (q2 + m_minus * m_minus - m_plus * m_plus) / (2.0 * root_q2)
    e_plus = (q2 + m_plus * m_plus - m_minus * m_minus) / (2.0 * root_q2)
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
    p_k = np.array([e_k, 0.0, 0.0, p_had], dtype=np.complex128)
    p_pi = np.array([e_pi, 0.0, 0.0, p_had], dtype=np.complex128)
    q_vec = np.array([root_q2, 0.0, 0.0, 0.0], dtype=np.complex128)
    p_minus, p_plus = _lepton_momenta(e_minus, e_plus, p_lep, cos_theta, sin_theta)
    ff = inputs.form_factor
    fplus = ff.fplus_0 / (1.0 - q2 / ff.vector_pole_mass_gev**2)
    fzero = ff.fzero_0 / (1.0 - q2 / ff.scalar_pole_mass_gev**2)
    hadronic = fplus * (p_k + p_pi - delta / q2 * q_vec) + fzero * delta / q2 * q_vec
    return _contract_vector_axial(hadronic, p_minus, p_plus, m_minus, m_plus)


def _dtopi_tensor_contractions(q2, cos_theta, *, m_minus, m_plus, inputs):
    root_q2 = math.sqrt(q2)
    m_d = inputs.dplus_mass_gev
    m_pi = inputs.piplus_mass_gev
    m_d2 = m_d * m_d
    m_pi2 = m_pi * m_pi
    delta = m_d2 - m_pi2
    p_had = math.sqrt(max(0.0, _kallen(m_d2, m_pi2, q2))) / (2.0 * root_q2)
    p_lep = math.sqrt(max(0.0, _kallen(q2, m_minus * m_minus, m_plus * m_plus))) / (
        2.0 * root_q2
    )
    e_d = (m_d2 - m_pi2 + q2) / (2.0 * root_q2)
    e_pi = (m_d2 - m_pi2 - q2) / (2.0 * root_q2)
    e_minus = (q2 + m_minus * m_minus - m_plus * m_plus) / (2.0 * root_q2)
    e_plus = (q2 + m_plus * m_plus - m_minus * m_minus) / (2.0 * root_q2)
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
    p_d = np.array([e_d, 0.0, 0.0, p_had], dtype=np.complex128)
    p_pi = np.array([e_pi, 0.0, 0.0, p_had], dtype=np.complex128)
    q_vec = np.array([root_q2, 0.0, 0.0, 0.0], dtype=np.complex128)
    p_minus, p_plus = _lepton_momenta(e_minus, e_plus, p_lep, cos_theta, sin_theta)
    fplus = dtopi_fplus(q2, inputs)
    fzero = dtopi_fzero(q2, inputs)
    hadronic = fplus * (p_d + p_pi - delta / q2 * q_vec) + fzero * delta / q2 * q_vec
    return _contract_vector_axial(hadronic, p_minus, p_plus, m_minus, m_plus)


def _lepton_momenta(e_minus, e_plus, p_lep, cos_theta, sin_theta):
    p_minus = np.array(
        [e_minus, p_lep * sin_theta, 0.0, p_lep * cos_theta],
        dtype=np.complex128,
    )
    p_plus = np.array(
        [e_plus, -p_lep * sin_theta, 0.0, -p_lep * cos_theta],
        dtype=np.complex128,
    )
    return p_minus, p_plus


def _contract_vector_axial(hadronic, p_minus, p_plus, m_minus, m_plus):
    hadronic_cov = _METRIC @ hadronic
    p_dot = float(np.real(p_minus @ _METRIC @ p_plus))
    vector_tensor = _lepton_tensor(p_minus, p_plus, p_dot + m_minus * m_plus)
    axial_tensor = _lepton_tensor(p_minus, p_plus, p_dot - m_minus * m_plus)
    return (
        _contract_hadronic_tensor(hadronic_cov, vector_tensor),
        _contract_hadronic_tensor(hadronic_cov, axial_tensor),
    )


def _lepton_tensor(p_minus, p_plus, metric_term):
    tensor = np.zeros((4, 4), dtype=np.complex128)
    for mu in range(4):
        for nu in range(4):
            tensor[mu, nu] = 4.0 * (
                p_minus[mu] * p_plus[nu]
                + p_minus[nu] * p_plus[mu]
                - _METRIC[mu, nu] * metric_term
            )
    return tensor


def _contract_hadronic_tensor(hadronic_cov, tensor):
    value = 0.0j
    for mu in range(4):
        for nu in range(4):
            value += hadronic_cov[mu] * np.conjugate(hadronic_cov[nu]) * tensor[mu, nu]
    return float(np.real(value))


def _kaon_y_norm(inputs):
    return float(
        math.sqrt(2.0)
        * inputs.gf_gev_minus2
        * inputs.alpha_em_mz
        / (math.pi * kaon_g_sm_squared(inputs))
    )


def _kallen(a, b, c):
    return float(a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c))


def _kallen_dimensionless(x, y):
    return float(1.0 + x * x + y * y - 2.0 * (x + y + x * y))
