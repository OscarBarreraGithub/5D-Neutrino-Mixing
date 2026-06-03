"""Shared rare-B Phase-3d RS-Wilson test helpers."""

from __future__ import annotations

from dataclasses import replace
from typing import Mapping

from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ParameterPoint
from quarkConstraints.rare_b_dilepton import RareBDileptonWilsonCoefficients
from quarkConstraints.rs_semileptonic_wilsons import (
    RSSemileptonicWilsonBundle,
    RSSemileptonicWilsonCoefficients,
)
from tests.rs_ew_phase3b_helpers import sample_rs_ew_point, sm_limit_rs_ew_point

_RARE_B_BLOCKS = ("b_to_s_ll", "b_to_d_ll")


def sm_limit_rare_b_point() -> ParameterPoint:
    """Return a Phase-3a point whose rare-B FCNC Wilsons vanish."""

    return sm_limit_rs_ew_point()


def sample_rare_b_point() -> ParameterPoint:
    """Return a Phase-3a point with non-zero rare-B FCNC Wilsons."""

    return sample_rs_ew_point()


def scaled_rare_b_point(
    *,
    scale: float = 1.0,
    lepton_scales: Mapping[str, float] | None = None,
    base: ParameterPoint | None = None,
    **extra_overrides,
) -> ParameterPoint:
    """Return ``base`` with rare-B C9/C10/C9p/C10p multiplied."""

    source = sample_rare_b_point() if base is None else base
    extras = dict(source.extras)
    extras["rs_semileptonic_wilsons"] = scaled_rare_b_wilson_bundle(
        extras["rs_semileptonic_wilsons"],
        scale=scale,
        lepton_scales=lepton_scales,
    )
    extras.update(extra_overrides)
    return point_builder.make_point(raw=source.raw, **extras)


def scaled_rare_b_wilson_bundle(
    bundle: RSSemileptonicWilsonBundle,
    *,
    scale: float = 1.0,
    lepton_scales: Mapping[str, float] | None = None,
) -> RSSemileptonicWilsonBundle:
    """Scale rare-B Wilson coefficients, preserving contacts for diagnostics."""

    replacements = {}
    for block_name in _RARE_B_BLOCKS:
        block = getattr(bundle, block_name)
        replacements[block_name] = {
            lepton: _scale_coeff(
                coeff,
                scale=float(scale) * float((lepton_scales or {}).get(lepton, 1.0)),
            )
            for lepton, coeff in block.items()
        }
    return replace(bundle, **replacements)


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


def core_wilsons_from_rs_coeff(
    coeff: RSSemileptonicWilsonCoefficients,
    *,
    matching_scale_gev: float = 3000.0,
) -> RareBDileptonWilsonCoefficients:
    """Build the rare-B core Wilson dataclass directly from a 3a coefficient."""

    left_contact = complex(coeff.contact_LL + coeff.contact_LR)
    right_contact = complex(coeff.contact_RL + coeff.contact_RR)
    return RareBDileptonWilsonCoefficients(
        model_label=coeff.model_label,
        operator_convention=coeff.operator_convention,
        matching_assumption=coeff.matching_assumption,
        transition_key=coeff.transition_key,
        M_KK=float(matching_scale_gev),
        matching_scale=float(matching_scale_gev),
        lambda_t=complex(coeff.lambda_ckm),
        left_qb_coupling=left_contact,
        right_qb_coupling=right_contact,
        left_qb_overlap=0.0j,
        right_qb_overlap=0.0j,
        left_quark_delta=left_contact,
        right_quark_delta=right_contact,
        muon_left_delta=0.0,
        muon_right_delta=0.0,
        muon_vector_delta=0.0,
        muon_axial_delta=0.0,
        c9_np=complex(coeff.c9_np),
        c10_np=complex(coeff.c10_np),
        c9p_np=complex(coeff.c9p_np),
        c10p_np=complex(coeff.c10p_np),
    )


def rs_coeff(
    point: ParameterPoint,
    *,
    transition: str = "b_s",
    lepton: str = "mu",
) -> RSSemileptonicWilsonCoefficients:
    """Return one rare-B 3a Wilson coefficient from a test point."""

    bundle = point.extras["rs_semileptonic_wilsons"]
    block = bundle.b_to_s_ll if transition == "b_s" else bundle.b_to_d_ll
    return block[lepton]
