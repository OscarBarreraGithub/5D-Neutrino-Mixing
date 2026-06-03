"""Shared rare-kaon Phase-3d RS-Wilson test helpers."""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Mapping

from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ParameterPoint
from quarkConstraints.rare_kaon_dilepton import (
    KLongPi0EEWilsonCoefficients,
    RareKaonDileptonSMInputs,
    RareKaonDileptonWilsonCoefficients,
    g_sm_squared,
)
from quarkConstraints.rs_semileptonic_wilsons import (
    RSSemileptonicWilsonBundle,
    RSSemileptonicWilsonCoefficients,
)
from tests.rs_ew_phase3b_helpers import sample_rs_ew_point, sm_limit_rs_ew_point


def sample_rare_kaon_point() -> ParameterPoint:
    """Return a Phase-3a point with non-zero rare-kaon FCNC Wilsons."""

    return sample_rs_ew_point()


def sm_limit_rare_kaon_point() -> ParameterPoint:
    """Return a Phase-3a point whose rare-kaon FCNC Wilsons vanish."""

    return sm_limit_rs_ew_point()


def scaled_rare_kaon_point(
    *,
    scale: float = 1.0,
    lepton_scales: Mapping[str, float] | None = None,
    base: ParameterPoint | None = None,
    **extra_overrides,
) -> ParameterPoint:
    """Return ``base`` with rare-kaon C9/C10/C9p/C10p multiplied."""

    source = sample_rare_kaon_point() if base is None else base
    extras = dict(source.extras)
    extras["rs_semileptonic_wilsons"] = scaled_rare_kaon_wilson_bundle(
        extras["rs_semileptonic_wilsons"],
        scale=scale,
        lepton_scales=lepton_scales,
    )
    extras.update(extra_overrides)
    return point_builder.make_point(raw=source.raw, **extras)


def rare_kaon_point_with_muon_y_np_total(
    y_np_total: complex,
    *,
    base: ParameterPoint | None = None,
    inputs: RareKaonDileptonSMInputs | None = None,
) -> ParameterPoint:
    """Return a point whose muon C10 block maps to the requested ``Y_NP``."""

    source = sample_rare_kaon_point() if base is None else base
    p = RareKaonDileptonSMInputs() if inputs is None else inputs
    bundle = source.extras["rs_semileptonic_wilsons"]
    coeff = bundle.s_to_d_ll["mu"]
    c10_np = -complex(y_np_total) / (
        complex(coeff.lambda_ckm) * c9_c10_to_y_norm(p)
    )
    new_mu = replace(
        coeff,
        c9_np=0.0j,
        c10_np=complex(c10_np),
        c9p_np=0.0j,
        c10p_np=0.0j,
    )
    extras = dict(source.extras)
    extras["rs_semileptonic_wilsons"] = replace(
        bundle,
        s_to_d_ll={**bundle.s_to_d_ll, "mu": new_mu},
    )
    return point_builder.make_point(raw=source.raw, **extras)


def scaled_rare_kaon_wilson_bundle(
    bundle: RSSemileptonicWilsonBundle,
    *,
    scale: float = 1.0,
    lepton_scales: Mapping[str, float] | None = None,
) -> RSSemileptonicWilsonBundle:
    """Scale rare-kaon Wilson coefficients, preserving contacts for diagnostics."""

    return replace(
        bundle,
        s_to_d_ll={
            lepton: _scale_coeff(
                coeff,
                scale=float(scale) * float((lepton_scales or {}).get(lepton, 1.0)),
            )
            for lepton, coeff in bundle.s_to_d_ll.items()
        },
    )


def _scale_coeff(
    coeff: RSSemileptonicWilsonCoefficients,
    *,
    scale: float,
) -> RSSemileptonicWilsonCoefficients:
    return replace(
        coeff,
        c9_np=complex(coeff.c9_np) * scale,
        c10_np=complex(coeff.c10_np) * scale,
        c9p_np=complex(coeff.c9p_np) * scale,
        c10p_np=complex(coeff.c10p_np) * scale,
    )


def rs_coeff(
    point: ParameterPoint,
    *,
    lepton: str = "mu",
) -> RSSemileptonicWilsonCoefficients:
    """Return one rare-kaon 3a Wilson coefficient from a test point."""

    return point.extras["rs_semileptonic_wilsons"].s_to_d_ll[lepton]


def c9_c10_to_y_norm(inputs: RareKaonDileptonSMInputs) -> float:
    """Return the normalization mapping dimensionless C9/C10 to kaon Y/y7."""

    return float(
        math.sqrt(2.0)
        * inputs.gf_gev_minus2
        * inputs.alpha_em_mz
        / (math.pi * g_sm_squared(inputs))
    )


def core_y_wilsons_from_rs_coeff(
    coeff: RSSemileptonicWilsonCoefficients,
    *,
    inputs: RareKaonDileptonSMInputs,
    matching_scale_gev: float = 3000.0,
) -> RareKaonDileptonWilsonCoefficients:
    """Build the rare-kaon Y core dataclass directly from a 3a coefficient."""

    norm = c9_c10_to_y_norm(inputs)
    lambda_ckm = complex(coeff.lambda_ckm)
    left_axial_contact = complex(coeff.contact_LR - coeff.contact_LL)
    right_axial_contact = complex(coeff.contact_RR - coeff.contact_RL)
    return RareKaonDileptonWilsonCoefficients(
        model_label=coeff.model_label,
        operator_convention=coeff.operator_convention,
        matching_assumption=coeff.matching_assumption,
        M_KK=float(matching_scale_gev),
        matching_scale=float(matching_scale_gev),
        left_sd_coupling=left_axial_contact,
        right_sd_coupling=right_axial_contact,
        left_sd_overlap=0.0j,
        right_sd_overlap=0.0j,
        left_quark_delta=left_axial_contact,
        right_quark_delta=right_axial_contact,
        muon_axial_delta=0.0,
        y_np_left=-lambda_ckm * norm * complex(coeff.c10_np),
        y_np_right=-lambda_ckm * norm * complex(coeff.c10p_np),
    )


def core_y7_wilsons_from_rs_coeff(
    coeff: RSSemileptonicWilsonCoefficients,
    *,
    inputs: RareKaonDileptonSMInputs,
    matching_scale_gev: float = 3000.0,
) -> KLongPi0EEWilsonCoefficients:
    """Build the rare-kaon y7 core dataclass directly from a 3a coefficient."""

    norm = c9_c10_to_y_norm(inputs)
    lambda_ckm = complex(coeff.lambda_ckm)
    left_vector_contact = complex(coeff.contact_LL + coeff.contact_LR)
    right_vector_contact = complex(coeff.contact_RL + coeff.contact_RR)
    return KLongPi0EEWilsonCoefficients(
        model_label=coeff.model_label,
        operator_convention=coeff.operator_convention,
        matching_assumption=coeff.matching_assumption,
        M_KK=float(matching_scale_gev),
        matching_scale=float(matching_scale_gev),
        low_scale_gev=1.0,
        left_sd_coupling=left_vector_contact,
        right_sd_coupling=right_vector_contact,
        left_sd_overlap=0.0j,
        right_sd_overlap=0.0j,
        left_quark_delta=left_vector_contact,
        right_quark_delta=right_vector_contact,
        quark_vector_delta=left_vector_contact + right_vector_contact,
        electron_vector_delta=0.0,
        electron_axial_delta=0.0,
        lambda_y7v_np_proxy=-lambda_ckm
        * norm
        * complex(coeff.c9_np + coeff.c9p_np),
        lambda_y7a_np_proxy=-lambda_ckm
        * norm
        * complex(coeff.c10_np + coeff.c10p_np),
    )
