"""Adapter over :mod:`quarkConstraints.rare_charm_semileptonic`.

This is the catalog boundary for semileptonic rare-charm modes.  It reuses
the shared ``c -> u l+ l-`` Wilson machinery from
``quarkConstraints.rare_charm_dilepton`` through the dedicated
``D+ -> pi+ mu+ mu-`` short-distance rate module.
"""

from __future__ import annotations

from quarkConstraints.couplings import QuarkMassBasisCouplings
from quarkConstraints.rare_charm_semileptonic import (
    RARE_CHARM_DTOPI_MUMU_FORM_FACTOR_MODEL_V1,
    RARE_CHARM_DTOPI_MUMU_INPUT_BUNDLE_V1,
    RARE_CHARM_DTOPI_MUMU_MODEL_V1,
    RARE_CHARM_DTOPI_MUMU_OPERATOR_CONVENTION,
    RARE_CHARM_DTOPI_MUMU_PARAMETRIZATION_CITATION,
    RARE_CHARM_DTOPI_MUMU_RESONANCE_LIMITATION_V1,
    RareCharmDToPiFormFactorInputs,
    RareCharmDToPiMuMuBranchingResult,
    RareCharmDToPiMuMuInputs,
    default_dtopi_mumu_inputs as _default_dtopi_mumu_inputs,
    dtopi_fplus as _dtopi_fplus,
    dtopi_fzero as _dtopi_fzero,
    dtopi_mumu_differential_branching_fraction as _differential_branching,
    dtopi_mumu_q2_range as _q2_range,
    dtopi_mumu_sm as _dtopi_mumu_sm,
    evaluate_dplus_to_piplus_mumu as _evaluate_dplus_to_piplus_mumu,
)

__all__ = [
    "QuarkMassBasisCouplings",
    "RARE_CHARM_DTOPI_MUMU_MODEL_V1",
    "RARE_CHARM_DTOPI_MUMU_INPUT_BUNDLE_V1",
    "RARE_CHARM_DTOPI_MUMU_FORM_FACTOR_MODEL_V1",
    "RARE_CHARM_DTOPI_MUMU_OPERATOR_CONVENTION",
    "RARE_CHARM_DTOPI_MUMU_PARAMETRIZATION_CITATION",
    "RARE_CHARM_DTOPI_MUMU_RESONANCE_LIMITATION_V1",
    "RareCharmDToPiFormFactorInputs",
    "RareCharmDToPiMuMuInputs",
    "RareCharmDToPiMuMuBranchingResult",
    "rare_charm_dtopi_mumu_default_inputs",
    "rare_charm_dtopi_fplus",
    "rare_charm_dtopi_fzero",
    "rare_charm_dtopi_mumu_q2_range",
    "rare_charm_dtopi_mumu_differential_branching_fraction",
    "dplus_piplus_mumu_sm",
    "dplus_piplus_mumu_from_couplings",
]


def rare_charm_dtopi_mumu_default_inputs() -> RareCharmDToPiMuMuInputs:
    """Return the default ``D+ -> pi+ mu+ mu-`` short-distance input bundle."""
    return _default_dtopi_mumu_inputs()


def rare_charm_dtopi_fplus(
    q2_gev2: float,
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> float:
    """Return the exclusive ``D -> pi`` vector form factor ``f_+(q2)``."""
    return _dtopi_fplus(q2_gev2, inputs)


def rare_charm_dtopi_fzero(
    q2_gev2: float,
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> float:
    """Return the exclusive ``D -> pi`` scalar form factor ``f_0(q2)``."""
    return _dtopi_fzero(q2_gev2, inputs)


def rare_charm_dtopi_mumu_q2_range(
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> tuple[float, float]:
    """Return the kinematic dimuon ``q2`` range for ``D+ -> pi+ mu+ mu-``."""
    return _q2_range(inputs)


def rare_charm_dtopi_mumu_differential_branching_fraction(
    q2_gev2: float,
    *,
    c9_semileptonic: complex,
    c10_semileptonic: complex,
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> float:
    """Evaluate ``dBR_SD(D+ -> pi+ mu+ mu-) / dq2``."""
    return _differential_branching(
        q2_gev2,
        c9_semileptonic=c9_semileptonic,
        c10_semileptonic=c10_semileptonic,
        inputs=inputs,
    )


def dplus_piplus_mumu_sm(
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> RareCharmDToPiMuMuBranchingResult:
    """Evaluate the zero-Wilson SM-limit of the smooth SD proxy."""
    return _dtopi_mumu_sm(inputs)


def dplus_piplus_mumu_from_couplings(
    couplings: QuarkMassBasisCouplings,
    *,
    m_kk_gev: float | None = None,
    inputs: RareCharmDToPiMuMuInputs | None = None,
) -> RareCharmDToPiMuMuBranchingResult:
    """Evaluate smooth ``BR_SD(D+ -> pi+ mu+ mu-)`` from mass-basis couplings."""
    return _evaluate_dplus_to_piplus_mumu(
        couplings,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )
