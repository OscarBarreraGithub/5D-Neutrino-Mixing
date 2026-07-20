"""Semileptonic Wilson coefficients from RS-EW neutral contacts.

The coefficients in this module are matched directly from physical contact
terms in ``GeV^-2``.  They intentionally do not route through the older proxy
helpers that insert a single ``1/M_KK^2`` factor.  The charged-lepton
``C9/C10`` blocks use the same WET convention as ``rare_b_dilepton`` and
``rare_charm_dilepton``; chiral lepton contacts contribute to vector/axial
operators with the explicit ``1/2`` from ``P_{L,R}``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

import numpy as np

from . import rare_b_dilepton as _rare_b
from . import rare_b_nunu as _rare_b_nunu
from . import rare_charm_dilepton as _rare_charm
from . import rare_kaon_dilepton as _rare_kaon
from . import rare_kaon_snd as _rare_kaon_nunu
from .rs_ew_couplings import LEPTON_FLAVORS, RSEWMassBasisCouplings

RS_SEMILEPTONIC_WILSONS_MODEL_V1 = "RS_EW_SEMILEPTONIC_WILSONS_PHASE4A_V1"
RS_SEMILEPTONIC_OPERATOR_CONVENTION = (
    "H_eff=-4 G_F/sqrt(2) lambda alpha/(4 pi) "
    "[C9 O9 + C10 O10 + C9p O9p + C10p O10p]"
)
RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1 = (
    "C9/C10 matched from Phase-4a full light-Z neutral contacts in GeV^-2; "
    "charged-lepton blocks use +pi/(2 sqrt(2) G_F alpha lambda) with the "
    "P_L/P_R -> (V,A)/2 factor; active-nu blocks use X_NP=C/g_SM^2 directly "
    "with no second 1/M_KK^2"
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
class RSLFVSemileptonicWilsonCoefficients:
    """One charged-LFV semileptonic Wilson block."""

    model_label: str
    operator_convention: str
    matching_assumption: str
    transition_key: str
    quark_sector: str
    final_quark_index: int
    initial_quark_index: int
    lepton_pair_key: str
    final_lepton_index: int
    initial_lepton_index: int
    lambda_ckm_name: str
    lambda_ckm: complex
    gf_gev_minus2: float
    alpha_em_mz: float
    contact_units: str
    contact_LL: complex
    contact_LR: complex
    contact_RL: complex
    contact_RR: complex
    c9_lfv_np: complex
    c10_lfv_np: complex
    c9p_lfv_np: complex
    c10p_lfv_np: complex

    def __post_init__(self) -> None:
        if self.contact_units != "GeV^-2":
            raise ValueError("contact_units must be 'GeV^-2'")
        if self.quark_sector not in {"u", "d"}:
            raise ValueError("quark_sector must be 'u' or 'd'")
        if self.final_lepton_index == self.initial_lepton_index:
            raise ValueError("LFV lepton indices must be off diagonal")
        for index_name in ("final_lepton_index", "initial_lepton_index"):
            index = int(getattr(self, index_name))
            if index < 0 or index >= len(LEPTON_FLAVORS):
                raise ValueError(f"{index_name} is out of range")
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
            "c9_lfv_np",
            "c10_lfv_np",
            "c9p_lfv_np",
            "c10p_lfv_np",
        ):
            _require_finite_complex(getattr(self, name), name)

    @property
    def wilsons(self) -> Mapping[str, complex]:
        return {
            "C9_LFV_NP": complex(self.c9_lfv_np),
            "C10_LFV_NP": complex(self.c10_lfv_np),
            "C9p_LFV_NP": complex(self.c9p_lfv_np),
            "C10p_LFV_NP": complex(self.c10p_lfv_np),
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
class RSNuNuWilsonCoefficients:
    """One active-neutrino semileptonic Wilson matrix block."""

    model_label: str
    operator_convention: str
    matching_assumption: str
    transition_key: str
    quark_sector: str
    final_quark_index: int
    initial_quark_index: int
    g_sm_squared_gev_minus2: float
    contact_units: str
    contact_LL: np.ndarray
    contact_RL: np.ndarray
    x_np_left: np.ndarray
    x_np_right: np.ndarray

    def __post_init__(self) -> None:
        if self.contact_units != "GeV^-2":
            raise ValueError("contact_units must be 'GeV^-2'")
        if self.quark_sector not in {"u", "d"}:
            raise ValueError("quark_sector must be 'u' or 'd'")
        g_sm_squared = float(self.g_sm_squared_gev_minus2)
        if not math.isfinite(g_sm_squared) or g_sm_squared <= 0.0:
            raise ValueError("g_sm_squared_gev_minus2 must be positive and finite")
        for name in ("contact_LL", "contact_RL", "x_np_left", "x_np_right"):
            object.__setattr__(self, name, _readonly_complex_matrix(getattr(self, name), name))

    @property
    def wilsons(self) -> Mapping[str, np.ndarray]:
        return {
            "X_NP_L": self.x_np_left,
            "X_NP_R": self.x_np_right,
            "X_NP_total": self.x_np_left + self.x_np_right,
        }


@dataclass(frozen=True)
class RSSemileptonicWilsonBundle:
    """Phase-4a Wilson coefficients for semileptonic rare FCNC transitions."""

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
    lfv_llqq: Mapping[str, Mapping[str, RSLFVSemileptonicWilsonCoefficients]] = field(
        default_factory=dict
    )
    b_to_s_nunu: RSNuNuWilsonCoefficients | None = None
    s_to_d_nunu: RSNuNuWilsonCoefficients | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.contact_units != "GeV^-2":
            raise ValueError("contact_units must be 'GeV^-2'")
        for name in ("b_to_s_ll", "b_to_d_ll", "s_to_d_ll", "c_to_u_ll"):
            block = dict(getattr(self, name))
            if set(block) != set(LEPTON_FLAVORS):
                raise ValueError(f"{name} must contain e, mu, and tau entries")
            object.__setattr__(self, name, MappingProxyType(block))
        lfv_blocks: dict[str, Mapping[str, RSLFVSemileptonicWilsonCoefficients]] = {}
        for name, block in self.lfv_llqq.items():
            lfv_blocks[str(name)] = MappingProxyType(dict(block))
        object.__setattr__(self, "lfv_llqq", MappingProxyType(lfv_blocks))
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
        lfv_llqq={
            "b_to_s": _lfv_transition_block(
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
            "b_to_d": _lfv_transition_block(
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
            "s_to_d": _lfv_transition_block(
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
            "c_to_u": _lfv_transition_block(
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
        },
        b_to_s_nunu=_nunu_block(
            couplings,
            transition_key="b_s",
            quark_sector="d",
            final_quark_index=1,
            initial_quark_index=2,
            g_sm_squared_gev_minus2=_rare_b_nunu.g_sm_squared(),
        ),
        s_to_d_nunu=_nunu_block(
            couplings,
            transition_key="s_d",
            quark_sector="d",
            final_quark_index=0,
            initial_quark_index=1,
            g_sm_squared_gev_minus2=_rare_kaon_nunu.g_sm_squared(),
        ),
        metadata={
            "uses_physical_contact_units": True,
            "wilson_prefactor_reused": False,
            "second_mkk_suppression_applied": False,
            "lfv_llqq_blocks": "off-diagonal charged-lepton pairs",
            "nunu_mapping": "X_NP=C/g_SM^2",
            "nunu_wilson_prefactor_reused": False,
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


def _lfv_transition_block(
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
) -> dict[str, RSLFVSemileptonicWilsonCoefficients]:
    block: dict[str, RSLFVSemileptonicWilsonCoefficients] = {}
    for final_lepton_index, final_key in enumerate(LEPTON_FLAVORS):
        for initial_lepton_index, initial_key in enumerate(LEPTON_FLAVORS):
            if final_lepton_index == initial_lepton_index:
                continue
            pair_key = f"{final_key}_{initial_key}"
            block[pair_key] = _lfv_coefficient_block(
                couplings,
                transition_key=transition_key,
                quark_sector=quark_sector,
                final_quark_index=final_quark_index,
                initial_quark_index=initial_quark_index,
                lepton_pair_key=pair_key,
                final_lepton_index=final_lepton_index,
                initial_lepton_index=initial_lepton_index,
                lambda_ckm_name=lambda_ckm_name,
                lambda_ckm=lambda_ckm,
                gf_gev_minus2=gf_gev_minus2,
                alpha_em_mz=alpha_em_mz,
            )
    return block


def _lfv_coefficient_block(
    couplings: RSEWMassBasisCouplings,
    *,
    transition_key: str,
    quark_sector: str,
    final_quark_index: int,
    initial_quark_index: int,
    lepton_pair_key: str,
    final_lepton_index: int,
    initial_lepton_index: int,
    lambda_ckm_name: str,
    lambda_ckm: complex,
    gf_gev_minus2: float,
    alpha_em_mz: float,
) -> RSLFVSemileptonicWilsonCoefficients:
    i = int(final_quark_index)
    j = int(initial_quark_index)
    a = int(final_lepton_index)
    b = int(initial_lepton_index)
    c_ll = couplings.contact(quark_sector, "L", "L", i, j, a, b)
    c_lr = couplings.contact(quark_sector, "L", "R", i, j, a, b)
    c_rl = couplings.contact(quark_sector, "R", "L", i, j, a, b)
    c_rr = couplings.contact(quark_sector, "R", "R", i, j, a, b)
    prefactor = _dimensionless_wilson_prefactor(
        lambda_ckm=lambda_ckm,
        gf_gev_minus2=gf_gev_minus2,
        alpha_em_mz=alpha_em_mz,
    )
    c9_lfv_np = prefactor * (c_ll + c_lr)
    c10_lfv_np = prefactor * (c_lr - c_ll)
    c9p_lfv_np = prefactor * (c_rl + c_rr)
    c10p_lfv_np = prefactor * (c_rr - c_rl)
    return RSLFVSemileptonicWilsonCoefficients(
        model_label=RS_SEMILEPTONIC_WILSONS_MODEL_V1,
        operator_convention=RS_SEMILEPTONIC_OPERATOR_CONVENTION,
        matching_assumption=RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
        transition_key=transition_key,
        quark_sector=quark_sector,
        final_quark_index=i,
        initial_quark_index=j,
        lepton_pair_key=str(lepton_pair_key),
        final_lepton_index=a,
        initial_lepton_index=b,
        lambda_ckm_name=lambda_ckm_name,
        lambda_ckm=complex(lambda_ckm),
        gf_gev_minus2=float(gf_gev_minus2),
        alpha_em_mz=float(alpha_em_mz),
        contact_units=couplings.contact_units,
        contact_LL=complex(c_ll),
        contact_LR=complex(c_lr),
        contact_RL=complex(c_rl),
        contact_RR=complex(c_rr),
        c9_lfv_np=complex(c9_lfv_np),
        c10_lfv_np=complex(c10_lfv_np),
        c9p_lfv_np=complex(c9p_lfv_np),
        c10p_lfv_np=complex(c10p_lfv_np),
    )


def _nunu_block(
    couplings: RSEWMassBasisCouplings,
    *,
    transition_key: str,
    quark_sector: str,
    final_quark_index: int,
    initial_quark_index: int,
    g_sm_squared_gev_minus2: float,
) -> RSNuNuWilsonCoefficients:
    i = int(final_quark_index)
    j = int(initial_quark_index)
    g_sm_sq = float(g_sm_squared_gev_minus2)
    contact_ll = np.array(
        [
            [
                couplings.nunu_contact(quark_sector, "L", i, j, a, b)
                for b in range(3)
            ]
            for a in range(3)
        ],
        dtype=np.complex128,
    )
    contact_rl = np.array(
        [
            [
                couplings.nunu_contact(quark_sector, "R", i, j, a, b)
                for b in range(3)
            ]
            for a in range(3)
        ],
        dtype=np.complex128,
    )
    return RSNuNuWilsonCoefficients(
        model_label=RS_SEMILEPTONIC_WILSONS_MODEL_V1,
        operator_convention="O_LR=(qbar_i gamma_mu P_{L,R} q_j)(nubar_a gamma^mu P_L nu_b)",
        matching_assumption=RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1,
        transition_key=transition_key,
        quark_sector=quark_sector,
        final_quark_index=i,
        initial_quark_index=j,
        g_sm_squared_gev_minus2=g_sm_sq,
        contact_units=couplings.contact_units,
        contact_LL=contact_ll,
        contact_RL=contact_rl,
        x_np_left=contact_ll / g_sm_sq,
        x_np_right=contact_rl / g_sm_sq,
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
        math.pi
        / (
            2.0
            * math.sqrt(2.0)
            * float(gf_gev_minus2)
            * float(alpha_em_mz)
            * complex(lambda_ckm)
        )
    )


def _require_finite_complex(value: object, name: str) -> None:
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")


def _readonly_complex_matrix(values: object, name: str) -> np.ndarray:
    arr = np.array(values, dtype=np.complex128, copy=True)
    if arr.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    arr.setflags(write=False)
    return arr


__all__ = [
    "RSSemileptonicWilsonBundle",
    "RSSemileptonicWilsonCoefficients",
    "RSLFVSemileptonicWilsonCoefficients",
    "RSNuNuWilsonCoefficients",
    "RS_SEMILEPTONIC_MATCHING_ASSUMPTION_V1",
    "RS_SEMILEPTONIC_OPERATOR_CONVENTION",
    "RS_SEMILEPTONIC_WILSONS_MODEL_V1",
    "build_rs_semileptonic_wilsons",
]
