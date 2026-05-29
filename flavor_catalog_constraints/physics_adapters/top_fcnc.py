"""Adapter over :mod:`quarkConstraints.top_fcnc`.

This is the catalog boundary for top-FCNC two-body branching-ratio machinery.
Constraint modules import this adapter only; the effective-coupling formulas
and the documented RS top-Z proxy remain isolated in ``quarkConstraints``.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.top_fcnc import (
    TOP_FCNC_GLUON_DIPOLE_CONVENTION,
    TOP_FCNC_HIGGS_SCALAR_CONVENTION,
    TOP_FCNC_INPUT_BUNDLE_V1,
    TOP_FCNC_MODEL_V1,
    TOP_FCNC_PHOTON_DIPOLE_CONVENTION,
    TOP_FCNC_RS_Z_MATCHING_ASSUMPTION_V1,
    TOP_FCNC_Z_VECTOR_CONVENTION,
    TopFCNCBranchingResult,
    TopFCNCSMInputs,
    TopZFCNCProxyCouplings,
    compute_top_z_fcnc_proxy as _compute_top_z_fcnc_proxy,
    default_sm_inputs as _default_sm_inputs,
    evaluate_t_to_q_z as _evaluate_t_to_q_z,
    gluon_dipole_branching_fraction as _gluon_dipole_branching_fraction,
    higgs_scalar_branching_fraction as _higgs_scalar_branching_fraction,
    photon_dipole_branching_fraction as _photon_dipole_branching_fraction,
    top_gluon_dipole_partial_width as _top_gluon_dipole_partial_width,
    top_higgs_scalar_partial_width as _top_higgs_scalar_partial_width,
    top_photon_dipole_partial_width as _top_photon_dipole_partial_width,
    top_z_vector_partial_width as _top_z_vector_partial_width,
    weak_neutral_current_coupling as _weak_neutral_current_coupling,
    z_vector_branching_fraction as _z_vector_branching_fraction,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "TOP_FCNC_MODEL_V1",
    "TOP_FCNC_INPUT_BUNDLE_V1",
    "TOP_FCNC_Z_VECTOR_CONVENTION",
    "TOP_FCNC_PHOTON_DIPOLE_CONVENTION",
    "TOP_FCNC_GLUON_DIPOLE_CONVENTION",
    "TOP_FCNC_HIGGS_SCALAR_CONVENTION",
    "TOP_FCNC_RS_Z_MATCHING_ASSUMPTION_V1",
    "TopFCNCSMInputs",
    "TopZFCNCProxyCouplings",
    "TopFCNCBranchingResult",
    "top_fcnc_default_sm_inputs",
    "top_fcnc_weak_neutral_current_coupling",
    "top_fcnc_z_vector_partial_width",
    "top_fcnc_photon_dipole_partial_width",
    "top_fcnc_gluon_dipole_partial_width",
    "top_fcnc_higgs_scalar_partial_width",
    "top_fcnc_z_vector_branching_fraction",
    "top_fcnc_photon_dipole_branching_fraction",
    "top_fcnc_gluon_dipole_branching_fraction",
    "top_fcnc_higgs_scalar_branching_fraction",
    "top_z_fcnc_proxy_from_couplings",
    "t_to_q_z_from_couplings",
]


def top_fcnc_default_sm_inputs() -> TopFCNCSMInputs:
    """Return the default top-FCNC input bundle."""
    return _default_sm_inputs()


def top_fcnc_weak_neutral_current_coupling(
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``g/(2 c_W)`` from the top-FCNC input bundle."""
    return _weak_neutral_current_coupling(inputs)


def top_fcnc_z_vector_partial_width(
    vector_left: complex = 0.0j,
    vector_right: complex = 0.0j,
    *,
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``Gamma(t -> q Z)`` for vector effective couplings."""
    return _top_z_vector_partial_width(vector_left, vector_right, inputs=inputs)


def top_fcnc_photon_dipole_partial_width(
    dipole_left: complex = 0.0j,
    dipole_right: complex = 0.0j,
    *,
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``Gamma(t -> q gamma)`` for photon dipole couplings."""
    return _top_photon_dipole_partial_width(dipole_left, dipole_right, inputs=inputs)


def top_fcnc_gluon_dipole_partial_width(
    dipole_left: complex = 0.0j,
    dipole_right: complex = 0.0j,
    *,
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``Gamma(t -> q g)`` for chromomagnetic dipole couplings."""
    return _top_gluon_dipole_partial_width(dipole_left, dipole_right, inputs=inputs)


def top_fcnc_higgs_scalar_partial_width(
    scalar_left: complex = 0.0j,
    scalar_right: complex = 0.0j,
    *,
    inputs: TopFCNCSMInputs | None = None,
) -> float:
    """Return ``Gamma(t -> q h)`` for scalar FCNC couplings."""
    return _top_higgs_scalar_partial_width(scalar_left, scalar_right, inputs=inputs)


def top_fcnc_z_vector_branching_fraction(
    *,
    light_quark: str,
    vector_left: complex = 0.0j,
    vector_right: complex = 0.0j,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q Z)`` from vector effective couplings."""
    return _z_vector_branching_fraction(
        light_quark=light_quark,
        vector_left=vector_left,
        vector_right=vector_right,
        inputs=inputs,
    )


def top_fcnc_photon_dipole_branching_fraction(
    *,
    light_quark: str,
    dipole_left: complex = 0.0j,
    dipole_right: complex = 0.0j,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q gamma)`` from photon dipole couplings."""
    return _photon_dipole_branching_fraction(
        light_quark=light_quark,
        dipole_left=dipole_left,
        dipole_right=dipole_right,
        inputs=inputs,
    )


def top_fcnc_gluon_dipole_branching_fraction(
    *,
    light_quark: str,
    dipole_left: complex = 0.0j,
    dipole_right: complex = 0.0j,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q g)`` from chromomagnetic dipole couplings."""
    return _gluon_dipole_branching_fraction(
        light_quark=light_quark,
        dipole_left=dipole_left,
        dipole_right=dipole_right,
        inputs=inputs,
    )


def top_fcnc_higgs_scalar_branching_fraction(
    *,
    light_quark: str,
    scalar_left: complex = 0.0j,
    scalar_right: complex = 0.0j,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q h)`` from scalar effective couplings."""
    return _higgs_scalar_branching_fraction(
        light_quark=light_quark,
        scalar_left=scalar_left,
        scalar_right=scalar_right,
        inputs=inputs,
    )


def top_z_fcnc_proxy_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    light_quark: str,
    light_up_index: int,
    m_kk_gev: float | None = None,
    inputs: TopFCNCSMInputs | None = None,
) -> TopZFCNCProxyCouplings:
    """Return the v1 documented ``t-q-Z`` proxy from mass-basis couplings."""
    return _compute_top_z_fcnc_proxy(
        couplings,
        light_quark=light_quark,
        light_up_index=light_up_index,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def t_to_q_z_from_couplings(
    couplings: QuarkMassBasisCouplings | TopZFCNCProxyCouplings | None = None,
    *,
    light_quark: str,
    light_up_index: int,
    m_kk_gev: float | None = None,
    inputs: TopFCNCSMInputs | None = None,
) -> TopFCNCBranchingResult:
    """Evaluate ``BR(t -> q Z)`` from mass-basis couplings or a proxy object."""
    return _evaluate_t_to_q_z(
        couplings,
        light_quark=light_quark,
        light_up_index=light_up_index,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
