"""Semileptonic Wilson coefficients from RS-EW neutral contacts.

The coefficients in this module are matched directly from physical contact
terms in ``GeV^-2``.  They intentionally do not route through the older proxy
helpers that insert a single ``1/M_KK^2`` factor.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

import numpy as np

from . import rare_b_dilepton as _rare_b
from . import rare_charm_dilepton as _rare_charm
from . import rare_kaon_dilepton as _rare_kaon
from .rs_ew_couplings import LEPTON_FLAVORS, RSEWMassBasisCouplings


RS_SEMILEPTONIC_WILSONS_MODEL_V1 = "RS_EW_SEMILEPTONIC_WILSONS_PHASE3A_V1"
RS_SEMILEPTONIC_OPERATOR_CONVENTION = (
    "H_eff=-4 G_F/sqrt(2) lambda alpha/(4 pi) "
    "[C9 O9 + C10 O10 + C9p O9p + C10p O10p]"
)
RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1 = (
    "C9/C10 matched from Phase-3a light-Z neutral contacts in GeV^-2; "
    "no _wilson_prefactor helper and no second 1/M_KK^2"
)


@dataclass(frozen=True)
class RSSemileptonicWilsonCoefficients:
    """One same-flavor charged-lepton semileptonic Wilson block."""

    model_label: str
    operator_convention: str
    matching_assumption: str
    transition_key: str
    quark_sector: str
    final_quark_index: int
    initial_quark_index: int
    lepton_key: str
    lepton_index: int
    lambda_ckm_name: str
    lambda_ckm: complex
    gf_gev_minus2: float
    alpha_em_mz: float
    contact_units: str
    contact_LL: complex
    contact_LR: complex
    contact_RL: complex
    contact_RR: complex
    c9_np: complex
    c10_np: complex
    c9p_np: complex
    c10p_np: complex

    def __post_init__(self) -> None:
        if self.contact_units != "GeV^-2":
            raise ValueError("contact_units must be 'GeV^-2'")
        if self.quark_sector not in {"u", "d"}:
            raise ValueError("quark_sector must be 'u' or 'd'")
        if self.lepton_key not in LEPTON_FLAVORS:
            raise ValueError(f"unsupported lepton_key {self.lepton_key!r}")
        for name in ("gf_gev_minus2", "alpha_em_mz"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        for name in (
            "lambda_ckm",
            "contact_LL",
            "contact_LR",
            "contact_RL",
            "contact_RR",
            "c9_np",
            "c10_np",
            "c9p_np",
            "c10p_np",
        ):
            _require_finite_complex(getattr(self, name), name)

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {
            "C9_NP": complex(self.c9_np),
            "C10_NP": complex(self.c10_np),
            "C9p_NP": complex(self.c9p_np),
            "C10p_NP": complex(self.c10p_np),
        }

    @property
    def contacts(self) -> Mapping[str, complex]:
        return {
            "C_LL": complex(self.contact_LL),
            "C_LR": complex(self.contact_LR),
            "C_RL": complex(self.contact_RL),
            "C_RR": complex(self.contact_RR),
        }


@dataclass(frozen=True)
class RSSemileptonicWilsonBundle:
    """Phase-3a Wilson coefficients for charged-lepton rare FCNC transitions."""

    model_label: str
    operator_convention: str
    matching_assumption: str
    contact_units: str
    includes_heavy_neutral_exchange: bool
    includes_heavy_neutral_lepton: bool
    b_to_s_ll: Mapping[str, RSSemileptonicWilsonCoefficients]
    b_to_d_ll: Mapping[str, RSSemileptonicWilsonCoefficients]
    s_to_d_ll: Mapping[str, RSSemileptonicWilsonCoefficients]
    c_to_u_ll: Mapping[str, RSSemileptonicWilsonCoefficients]
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.contact_units != "GeV^-2":
            raise ValueError("contact_units must be 'GeV^-2'")
        for name in ("b_to_s_ll", "b_to_d_ll", "s_to_d_ll", "c_to_u_ll"):
            block = dict(getattr(self, name))
            if set(block) != set(LEPTON_FLAVORS):
                raise ValueError(f"{name} must contain e, mu, and tau entries")
            object.__setattr__(self, name, MappingProxyType(block))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    def transition(self, name: str) -> Mapping[str, RSSemileptonicWilsonCoefficients]:
        """Return a transition block by public bundle field name."""

        if name not in {"b_to_s_ll", "b_to_d_ll", "s_to_d_ll", "c_to_u_ll"}:
            raise ValueError(f"unsupported transition {name!r}")
        return getattr(self, name)


def build_rs_semileptonic_wilsons(
    couplings: RSEWMassBasisCouplings,
    *,
    rare_b_inputs: _rare_b.RareBDileptonSMInputs | None = None,
    rare_kaon_inputs: _rare_kaon.RareKaonDileptonSMInputs | None = None,
    rare_charm_inputs: _rare_charm.RareCharmDileptonSMInputs | None = None,
) -> RSSemileptonicWilsonBundle:
    """Build the Phase-3a charged-lepton semileptonic Wilson bundle."""

    b_inputs = _rare_b.default_sm_inputs() if rare_b_inputs is None else rare_b_inputs
    k_inputs = (
        _rare_kaon.default_sm_inputs() if rare_kaon_inputs is None else rare_kaon_inputs
    )
    c_inputs = (
        _rare_charm.default_sm_inputs()
        if rare_charm_inputs is None
        else rare_charm_inputs
    )

    lambda_bs = _rare_b.ckm_factors("b_s", b_inputs).lambda_t
    lambda_bd = _rare_b.ckm_factors("b_d", b_inputs).lambda_t
    lambda_sd = _rare_kaon.ckm_factors(k_inputs).lambda_t
    lambda_cu = _rare_charm.ckm_factors("c_u", c_inputs).lambda_b

    return RSSemileptonicWilsonBundle(
        model_label=RS_SEMILEPTONIC_WILSONS_MODEL_V1,
        operator_convention=RS_SEMILEPTONIC_OPERATOR_CONVENTION,
        matching_assumption=RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
        contact_units=couplings.contact_units,
        includes_heavy_neutral_exchange=bool(couplings.includes_heavy_neutral_exchange),
        includes_heavy_neutral_lepton=bool(couplings.includes_heavy_neutral_lepton),
        b_to_s_ll=_transition_block(
            couplings,
            transition_key="b_s",
            quark_sector="d",
            final_quark_index=1,
            initial_quark_index=2,
            lambda_ckm_name="lambda_t",
            lambda_ckm=lambda_bs,
            gf_gev_minus2=b_inputs.gf_gev_minus2,
            alpha_em_mz=b_inputs.alpha_em_mz,
        ),
        b_to_d_ll=_transition_block(
            couplings,
            transition_key="b_d",
            quark_sector="d",
            final_quark_index=0,
            initial_quark_index=2,
            lambda_ckm_name="lambda_t",
            lambda_ckm=lambda_bd,
            gf_gev_minus2=b_inputs.gf_gev_minus2,
            alpha_em_mz=b_inputs.alpha_em_mz,
        ),
        s_to_d_ll=_transition_block(
            couplings,
            transition_key="s_d",
            quark_sector="d",
            final_quark_index=0,
            initial_quark_index=1,
            lambda_ckm_name="lambda_t",
            lambda_ckm=lambda_sd,
            gf_gev_minus2=k_inputs.gf_gev_minus2,
            alpha_em_mz=k_inputs.alpha_em_mz,
        ),
        c_to_u_ll=_transition_block(
            couplings,
            transition_key="c_u",
            quark_sector="u",
            final_quark_index=0,
            initial_quark_index=1,
            lambda_ckm_name="lambda_b",
            lambda_ckm=lambda_cu,
            gf_gev_minus2=c_inputs.gf_gev_minus2,
            alpha_em_mz=c_inputs.alpha_em_mz,
        ),
        metadata={
            "uses_physical_contact_units": True,
            "wilson_prefactor_reused": False,
            "second_mkk_suppression_applied": False,
        },
    )


def _transition_block(
    couplings: RSEWMassBasisCouplings,
    *,
    transition_key: str,
    quark_sector: str,
    final_quark_index: int,
    initial_quark_index: int,
    lambda_ckm_name: str,
    lambda_ckm: complex,
    gf_gev_minus2: float,
    alpha_em_mz: float,
) -> dict[str, RSSemileptonicWilsonCoefficients]:
    return {
        lepton_key: _coefficient_block(
            couplings,
            transition_key=transition_key,
            quark_sector=quark_sector,
            final_quark_index=final_quark_index,
            initial_quark_index=initial_quark_index,
            lepton_key=lepton_key,
            lepton_index=lepton_index,
            lambda_ckm_name=lambda_ckm_name,
            lambda_ckm=lambda_ckm,
            gf_gev_minus2=gf_gev_minus2,
            alpha_em_mz=alpha_em_mz,
        )
        for lepton_index, lepton_key in enumerate(LEPTON_FLAVORS)
    }


def _coefficient_block(
    couplings: RSEWMassBasisCouplings,
    *,
    transition_key: str,
    quark_sector: str,
    final_quark_index: int,
    initial_quark_index: int,
    lepton_key: str,
    lepton_index: int,
    lambda_ckm_name: str,
    lambda_ckm: complex,
    gf_gev_minus2: float,
    alpha_em_mz: float,
) -> RSSemileptonicWilsonCoefficients:
    i = int(final_quark_index)
    j = int(initial_quark_index)
    a = int(lepton_index)
    c_ll = couplings.contact(quark_sector, "L", "L", i, j, a, a)
    c_lr = couplings.contact(quark_sector, "L", "R", i, j, a, a)
    c_rl = couplings.contact(quark_sector, "R", "L", i, j, a, a)
    c_rr = couplings.contact(quark_sector, "R", "R", i, j, a, a)
    prefactor = _dimensionless_wilson_prefactor(
        lambda_ckm=lambda_ckm,
        gf_gev_minus2=gf_gev_minus2,
        alpha_em_mz=alpha_em_mz,
    )
    c9_np = prefactor * (c_ll + c_lr)
    c10_np = prefactor * (c_lr - c_ll)
    c9p_np = prefactor * (c_rl + c_rr)
    c10p_np = prefactor * (c_rr - c_rl)
    return RSSemileptonicWilsonCoefficients(
        model_label=RS_SEMILEPTONIC_WILSONS_MODEL_V1,
        operator_convention=RS_SEMILEPTONIC_OPERATOR_CONVENTION,
        matching_assumption=RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
        transition_key=transition_key,
        quark_sector=quark_sector,
        final_quark_index=i,
        initial_quark_index=j,
        lepton_key=lepton_key,
        lepton_index=a,
        lambda_ckm_name=lambda_ckm_name,
        lambda_ckm=complex(lambda_ckm),
        gf_gev_minus2=float(gf_gev_minus2),
        alpha_em_mz=float(alpha_em_mz),
        contact_units=couplings.contact_units,
        contact_LL=complex(c_ll),
        contact_LR=complex(c_lr),
        contact_RL=complex(c_rl),
        contact_RR=complex(c_rr),
        c9_np=complex(c9_np),
        c10_np=complex(c10_np),
        c9p_np=complex(c9p_np),
        c10p_np=complex(c10p_np),
    )


def _dimensionless_wilson_prefactor(
    *,
    lambda_ckm: complex,
    gf_gev_minus2: float,
    alpha_em_mz: float,
) -> complex:
    _require_finite_complex(lambda_ckm, "lambda_ckm")
    if abs(lambda_ckm) <= 0.0:
        raise ValueError("lambda_ckm must be non-zero")
    return complex(
        -math.pi
        / (
            math.sqrt(2.0)
            * float(gf_gev_minus2)
            * float(alpha_em_mz)
            * complex(lambda_ckm)
        )
    )


def _require_finite_complex(value: object, name: str) -> None:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")


__all__ = [
    "RSSemileptonicWilsonBundle",
    "RSSemileptonicWilsonCoefficients",
    "RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1",
    "RS_SEMILEPTONIC_OPERATOR_CONVENTION",
    "RS_SEMILEPTONIC_WILSONS_MODEL_V1",
    "build_rs_semileptonic_wilsons",
]
