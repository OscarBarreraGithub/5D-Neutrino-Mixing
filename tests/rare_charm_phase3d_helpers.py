"""Shared rare-charm Phase-3d RS-Wilson test helpers."""

from __future__ import annotations

from dataclasses import replace
from typing import Mapping

from flavor_catalog_constraints import point_builder
from flavor_catalog_constraints.base import ParameterPoint
from quarkConstraints.rare_charm_dilepton import RareCharmDileptonWilsonCoefficients
from quarkConstraints.rs_semileptonic_wilsons import (
    RSSemileptonicWilsonBundle,
    RSSemileptonicWilsonCoefficients,
)
from tests.rs_ew_phase3b_helpers import sample_rs_ew_point, sm_limit_rs_ew_point


def sample_rare_charm_point() -> ParameterPoint:
    """Return a Phase-3a point with non-zero rare-charm FCNC Wilsons."""

    return sample_rs_ew_point()


def sm_limit_rare_charm_point() -> ParameterPoint:
    """Return a Phase-3a point whose rare-charm FCNC Wilsons vanish."""

    return sm_limit_rs_ew_point()


def scaled_rare_charm_point(
    *,
    scale: float = 1.0,
    lepton_scales: Mapping[str, float] | None = None,
    base: ParameterPoint | None = None,
    **extra_overrides,
) -> ParameterPoint:
    """Return ``base`` with rare-charm C9/C10/C9p/C10p multiplied."""

    source = sample_rare_charm_point() if base is None else base
    extras = dict(source.extras)
    extras["rs_semileptonic_wilsons"] = scaled_rare_charm_wilson_bundle(
        extras["rs_semileptonic_wilsons"],
        scale=scale,
        lepton_scales=lepton_scales,
    )
    extras.update(extra_overrides)
    return point_builder.make_point(raw=source.raw, **extras)


def rare_charm_point_with_wilsons(
    *,
    lepton: str,
    c9_np: complex = 0.0j,
    c10_np: complex = 0.0j,
    c9p_np: complex = 0.0j,
    c10p_np: complex = 0.0j,
    base: ParameterPoint | None = None,
) -> ParameterPoint:
    """Return a point whose one lepton entry has explicit rare-charm Wilsons."""

    source = sample_rare_charm_point() if base is None else base
    bundle = source.extras["rs_semileptonic_wilsons"]
    coeff = bundle.c_to_u_ll[lepton]
    replacement = replace(
        coeff,
        c9_np=complex(c9_np),
        c10_np=complex(c10_np),
        c9p_np=complex(c9p_np),
        c10p_np=complex(c10p_np),
    )
    extras = dict(source.extras)
    extras["rs_semileptonic_wilsons"] = replace(
        bundle,
        c_to_u_ll={**bundle.c_to_u_ll, lepton: replacement},
    )
    return point_builder.make_point(raw=source.raw, **extras)


def scaled_rare_charm_wilson_bundle(
    bundle: RSSemileptonicWilsonBundle,
    *,
    scale: float = 1.0,
    lepton_scales: Mapping[str, float] | None = None,
) -> RSSemileptonicWilsonBundle:
    """Scale rare-charm Wilson coefficients, preserving contacts for diagnostics."""

    return replace(
        bundle,
        c_to_u_ll={
            lepton: _scale_coeff(
                coeff,
                scale=float(scale) * float((lepton_scales or {}).get(lepton, 1.0)),
            )
            for lepton, coeff in bundle.c_to_u_ll.items()
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
    """Return one rare-charm 3a Wilson coefficient from a test point."""

    return point.extras["rs_semileptonic_wilsons"].c_to_u_ll[lepton]


def core_dilepton_wilsons_from_rs_coeff(
    coeff: RSSemileptonicWilsonCoefficients,
    *,
    matching_scale_gev: float = 3000.0,
) -> RareCharmDileptonWilsonCoefficients:
    """Build the rare-charm core dataclass directly from a 3a coefficient."""

    left_vector_contact = complex(coeff.contact_LL + coeff.contact_LR)
    right_vector_contact = complex(coeff.contact_RL + coeff.contact_RR)
    return RareCharmDileptonWilsonCoefficients(
        model_label=coeff.model_label,
        operator_convention=coeff.operator_convention,
        matching_assumption=coeff.matching_assumption,
        transition_key="c_u",
        M_KK=float(matching_scale_gev),
        matching_scale=float(matching_scale_gev),
        lambda_b=complex(coeff.lambda_ckm),
        left_uc_coupling=left_vector_contact,
        right_uc_coupling=right_vector_contact,
        left_uc_overlap=0.0j,
        right_uc_overlap=0.0j,
        left_quark_delta=left_vector_contact,
        right_quark_delta=right_vector_contact,
        lepton_left_delta=0.0,
        lepton_right_delta=0.0,
        lepton_vector_delta=0.0,
        lepton_axial_delta=0.0,
        c9_np=complex(coeff.c9_np),
        c10_np=complex(coeff.c10_np),
        c9p_np=complex(coeff.c9p_np),
        c10p_np=complex(coeff.c10p_np),
    )
