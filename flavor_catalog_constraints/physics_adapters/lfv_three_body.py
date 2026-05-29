"""Adapter for charged-LFV three-body lepton decays.

Constraint modules import this adapter only.  It composes the existing L001
dipole adapter with the new :mod:`quarkConstraints.lfv_three_body` contact
formula and proxy matching.

NEEDS-HUMAN-PHYSICS: ``ParameterPoint`` does not carry the full lepton-sector
RS neutral-current and box-matching data needed for rigorous ``l_i -> 3 l_j``.
The contact terms are therefore caller-supplied/documented proxies.
"""

from __future__ import annotations

from typing import Any, Mapping

from quarkConstraints.lfv_three_body import (
    LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION,
    LFV_THREE_BODY_INPUT_BUNDLE_V1,
    LFV_THREE_BODY_MODEL_V1,
    LFV_THREE_BODY_OPERATOR_CONVENTION,
    LFV_THREE_BODY_PROXY_V1,
    LFVThreeBodyBranchingResult,
    LFVThreeBodyContactAmplitudes,
    LFVThreeBodyContactProxyInput,
    LFVThreeBodySMInputs,
    default_sm_inputs as _default_sm_inputs,
    lfv_three_body_contact_amplitudes as _contact_amplitudes,
    lfv_three_body_from_components as _from_components,
    lfv_three_body_has_contact_proxy as _has_contact_proxy,
    lfv_three_body_proxy_input as _proxy_input,
)

from .lepton import mu_to_e_gamma_from_lepton_input

__all__ = [
    "LFV_THREE_BODY_INPUT_BUNDLE_V1",
    "LFV_THREE_BODY_MODEL_V1",
    "LFV_THREE_BODY_OPERATOR_CONVENTION",
    "LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION",
    "LFV_THREE_BODY_PROXY_V1",
    "LFVThreeBodySMInputs",
    "LFVThreeBodyContactProxyInput",
    "LFVThreeBodyContactAmplitudes",
    "LFVThreeBodyBranchingResult",
    "lfv_three_body_default_sm_inputs",
    "lfv_three_body_proxy_input",
    "lfv_three_body_contact_amplitudes",
    "mu_to_3e_from_lepton_input",
]


def lfv_three_body_default_sm_inputs() -> LFVThreeBodySMInputs:
    """Return the default low-energy input bundle."""

    return _default_sm_inputs()


def lfv_three_body_proxy_input(
    left_lfv_overlap: complex,
    right_lfv_overlap: complex,
    m_kk_gev: float,
    *,
    initial_flavor: str = "mu",
    final_flavor: str = "e",
    box_ll: complex = 0.0j,
    box_lr: complex = 0.0j,
    box_rl: complex = 0.0j,
    box_rr: complex = 0.0j,
    source: str = "caller-supplied lfv three-body contact proxy",
) -> LFVThreeBodyContactProxyInput:
    """Build a contact proxy input accepted by this adapter."""

    return _proxy_input(
        left_lfv_overlap,
        right_lfv_overlap,
        m_kk_gev,
        initial_flavor=initial_flavor,
        final_flavor=final_flavor,
        box_ll=box_ll,
        box_lr=box_lr,
        box_rl=box_rl,
        box_rr=box_rr,
        source=source,
    )


def lfv_three_body_contact_amplitudes(
    lepton_input: Any,
    *,
    initial_flavor: str = "mu",
    final_flavor: str = "e",
    m_kk_gev: float | None = None,
    inputs: LFVThreeBodySMInputs | None = None,
) -> LFVThreeBodyContactAmplitudes:
    """Return documented Z-penguin plus box contact amplitudes."""

    return _contact_amplitudes(
        lepton_input,
        initial_flavor=initial_flavor,
        final_flavor=final_flavor,
        m_kk_gev=m_kk_gev,
        inputs=inputs,
    )


def mu_to_3e_from_lepton_input(
    lepton_input: Any,
    *,
    br_limit: float,
    dipole_br_limit: float,
    dipole_prefactor_br: float,
    reference_scale_gev: float = 3000.0,
    m_kk_gev: float | None = None,
    inputs: LFVThreeBodySMInputs | None = None,
) -> LFVThreeBodyBranchingResult:
    """Evaluate ``BR(mu -> 3e)`` from dipole and contact proxy inputs."""

    p = _default_sm_inputs() if inputs is None else inputs
    has_dipole = _has_dipole_input(lepton_input)
    has_contact = _has_contact_proxy(lepton_input)
    if not has_dipole and not has_contact:
        raise TypeError(
            "lepton input must provide a dipole source and/or lfv three-body "
            "contact proxy inputs"
        )
    _assert_mu_to_3e_flavors(lepton_input)

    dipole_parent_br = 0.0
    dipole_diagnostics: dict[str, Any] = {
        "dipole_input_present": bool(has_dipole),
        "dipole_component_reuses": "flavor_catalog_constraints.physics_adapters.lepton",
    }
    resolved_m_kk = m_kk_gev
    if has_dipole:
        dipole_source = _dipole_source(lepton_input)
        try:
            dipole = mu_to_e_gamma_from_lepton_input(
                dipole_source,
                br_limit=dipole_br_limit,
                prefactor_br=dipole_prefactor_br,
                reference_scale_gev=reference_scale_gev,
                m_kk_gev=m_kk_gev,
            )
        except (KeyError, TypeError, ValueError) as exc:
            if not has_contact:
                raise
            dipole_diagnostics.update(
                {
                    "dipole_evaluation_failed": True,
                    "dipole_exception_type": type(exc).__name__,
                    "dipole_exception": str(exc),
                }
            )
        else:
            dipole_parent_br = float(dipole.branching_fraction)
            resolved_m_kk = float(dipole.m_kk_gev)
            dipole_diagnostics.update(
                {
                    "dipole_evaluation_failed": False,
                    "dipole_parent_branching_fraction": float(
                        dipole.branching_fraction
                    ),
                    "dipole_parent_ratio_to_limit": float(dipole.ratio_to_limit),
                    "dipole_parent_limit": float(dipole.br_limit),
                    "dipole_lhs": float(dipole.dipole_lhs),
                    "dipole_rhs": float(dipole.dipole_rhs),
                    "dipole_off_diagonal_12": complex(dipole.off_diagonal_12),
                    "dipole_product_matrix": dipole.product_matrix,
                    "dipole_input_kind": dipole.input_kind,
                    "dipole_used_proxy": bool(dipole.used_proxy),
                }
            )

    contact = _contact_amplitudes(
        lepton_input,
        initial_flavor="mu",
        final_flavor="e",
        m_kk_gev=resolved_m_kk,
        inputs=p,
    )
    result = _from_components(
        dipole_parent_branching_fraction=dipole_parent_br,
        contact_amplitudes=contact,
        br_limit=br_limit,
        inputs=p,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics.update(dipole_diagnostics)
    return LFVThreeBodyBranchingResult(
        model_label=result.model_label,
        input_bundle=result.input_bundle,
        initial_flavor=result.initial_flavor,
        final_flavor=result.final_flavor,
        branching_fraction=result.branching_fraction,
        sm_branching_fraction=result.sm_branching_fraction,
        np_shift_branching_fraction=result.np_shift_branching_fraction,
        dipole_component=result.dipole_component,
        z_penguin_component=result.z_penguin_component,
        box_component=result.box_component,
        z_box_interference_component=result.z_box_interference_component,
        contact_component=result.contact_component,
        dipole_conversion_factor=result.dipole_conversion_factor,
        dipole_parent_branching_fraction=result.dipole_parent_branching_fraction,
        dipole_contact_interference_component=(
            result.dipole_contact_interference_component
        ),
        dipole_contact_interference_lower=result.dipole_contact_interference_lower,
        dipole_contact_interference_upper=result.dipole_contact_interference_upper,
        dipole_contact_interference_treatment=(
            result.dipole_contact_interference_treatment
        ),
        dipole_amplitude_left=result.dipole_amplitude_left,
        dipole_amplitude_right=result.dipole_amplitude_right,
        ratio_to_limit=result.ratio_to_limit,
        br_limit=result.br_limit,
        passes=result.passes,
        contact_amplitudes=result.contact_amplitudes,
        diagnostics=diagnostics,
    )


def _has_dipole_input(value: Any) -> bool:
    if isinstance(value, Mapping):
        if "dipole" in value or "yukawa_result" in value:
            return True
        has_y = any(key in value for key in ("y_n_bar", "Y_N_bar"))
        has_u = any(key in value for key in ("pmns", "PMNS", "pmns_matrix"))
        return bool(has_y and has_u)
    return hasattr(value, "Y_N_matrix") or (
        hasattr(value, "y_n_bar") and hasattr(value, "pmns")
    )


def _dipole_source(value: Any) -> Any:
    if isinstance(value, Mapping) and "dipole" in value:
        return value["dipole"]
    return value


_MU_TO_3E_INITIAL = "mu"
_MU_TO_3E_FINAL = "e"
_INITIAL_FLAVOR_KEYS = ("initial_flavor", "parent_flavor")
_FINAL_FLAVOR_KEYS = ("final_flavor", "daughter_flavor")
_CHARGED_LEPTON_ALIASES = {
    "electron": "e",
    "e-": "e",
    "e+": "e",
    "muon": "mu",
    "mu-": "mu",
    "mu+": "mu",
    "tau-": "tau",
    "tau+": "tau",
}
_CHARGED_LEPTONS = frozenset({"e", "mu", "tau"})


def _assert_mu_to_3e_flavors(value: Any) -> None:
    initial = _first_explicit_flavor(value, _INITIAL_FLAVOR_KEYS)
    final = _first_explicit_flavor(value, _FINAL_FLAVOR_KEYS)
    if initial is not None and _canonical_flavor(initial) != _MU_TO_3E_INITIAL:
        raise ValueError(
            "mu_to_3e_from_lepton_input is pinned to initial_flavor='mu'; "
            f"got {initial!r}"
        )
    if final is not None and _canonical_flavor(final) != _MU_TO_3E_FINAL:
        raise ValueError(
            "mu_to_3e_from_lepton_input is pinned to final_flavor='e'; "
            f"got {final!r}"
        )


def _first_explicit_flavor(value: Any, keys: tuple[str, ...]) -> Any:
    if isinstance(value, LFVThreeBodyContactProxyInput):
        if keys == _INITIAL_FLAVOR_KEYS:
            return value.initial_flavor
        return value.final_flavor
    if isinstance(value, Mapping):
        for key in keys:
            if key in value:
                return value[key]
        return None
    for key in keys:
        if hasattr(value, key):
            return getattr(value, key)
    return None


def _canonical_flavor(value: Any) -> str:
    key = _CHARGED_LEPTON_ALIASES.get(str(value), str(value))
    if key not in _CHARGED_LEPTONS:
        raise ValueError(f"unsupported charged-lepton flavor {value!r}")
    return key
