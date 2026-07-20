"""Adapter for charged-LFV three-body lepton decays.

Constraint modules import this adapter only.  It composes the existing L001
dipole adapter with the new :mod:`quarkConstraints.lfv_three_body` contact
formula and proxy matching.

Tree-level light-Z contacts are read from the Phase-4a ``rs_ew_couplings``
charged-lepton Z matrices when available.  Dipole matching, the
dipole-contact phase, and four-lepton boxes remain deferred loop-level inputs.
"""

from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np

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
)
from quarkConstraints.lfv_three_body import (
    default_sm_inputs as _default_sm_inputs,
)
from quarkConstraints.lfv_three_body import (
    lfv_three_body_from_components as _from_components,
)
from quarkConstraints.lfv_three_body import (
    lfv_three_body_proxy_input as _proxy_input,
)

from .lepton import mu_to_e_gamma_from_lepton_input

__all__ = [
    "LFV_THREE_BODY_INPUT_BUNDLE_V1",
    "LFV_THREE_BODY_MODEL_V1",
    "LFV_THREE_BODY_OPERATOR_CONVENTION",
    "LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION",
    "LFV_THREE_BODY_PROXY_V1",
    "LFV_THREE_BODY_TREE_CONTACT_RIGOROUS_V1",
    "LFV_THREE_BODY_DEFERRED_PIECES_V1",
    "LFVThreeBodySMInputs",
    "LFVThreeBodyContactProxyInput",
    "LFVThreeBodyContactAmplitudes",
    "LFVThreeBodyBranchingResult",
    "lfv_three_body_default_sm_inputs",
    "lfv_three_body_has_tree_contact_input",
    "lfv_three_body_proxy_input",
    "lfv_three_body_contact_amplitudes",
    "mu_to_3e_from_lepton_input",
]

LFV_THREE_BODY_TREE_CONTACT_RIGOROUS_V1 = (
    "tree-level light-Z LFV four-lepton contact from "
    "rs_ew_couplings.z_delta_g_{L,R}_e and z_total_g_{L,R}_e; zero for the "
    "Phase-4a diagonal charged-lepton fit with U_e=I and universal c_L"
)
LFV_THREE_BODY_DEFERRED_PIECES_V1 = (
    "NEEDS-HUMAN-PHYSICS: loop dipole matching, dipole-contact relative phase "
    "when chiral dipoles are not supplied, heavy neutral exchange, and "
    "four-lepton box matching remain deferred"
)
_TREE_CONTACT_MISSING_ASSUMPTION = (
    "tree-level light-Z LFV four-lepton contact not evaluated because "
    "rs_ew_couplings is absent; zero contact placeholder used while any "
    "available partial dipole/box pieces are preserved"
)
_ZERO_CONTACT_REFERENCE_SCALE_GEV = 3000.0
_BOX_LL_KEYS = ("box_ll", "box_left_left", "box_L_L")
_BOX_LR_KEYS = ("box_lr", "box_left_right", "box_L_R")
_BOX_RL_KEYS = ("box_rl", "box_right_left", "box_R_L")
_BOX_RR_KEYS = ("box_rr", "box_right_right", "box_R_R")


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
    """Return tree light-Z plus deferred box contact amplitudes."""

    if _has_rs_ew_couplings(lepton_input):
        return _rigorous_tree_contact_amplitudes(
            lepton_input,
            initial_flavor=initial_flavor,
            final_flavor=final_flavor,
            m_kk_gev=m_kk_gev,
            inputs=inputs,
        )
    return _zero_missing_tree_contact_amplitudes(
        lepton_input,
        initial_flavor=initial_flavor,
        final_flavor=final_flavor,
        m_kk_gev=m_kk_gev,
    )


def lfv_three_body_has_tree_contact_input(value: Any) -> bool:
    """Return true for rigorous tree contacts or deferred contact pieces."""

    return _has_rs_ew_couplings(value) or _has_box_amplitude_input(value)


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
    _assert_mu_to_3e_flavors(lepton_input)
    has_dipole = _has_dipole_input(lepton_input)
    has_contact = lfv_three_body_has_tree_contact_input(lepton_input)
    if not has_dipole and not has_contact:
        raise TypeError(
            "lepton input must provide a dipole source and/or lfv three-body "
            "tree contact inputs"
        )

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

    contact = lfv_three_body_contact_amplitudes(
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
        if "yukawa_result" in value:
            return True
        if "dipole" in value and _has_dipole_input(value["dipole"]):
            return True
        if "lepton_mass_basis_couplings" in value and _has_dipole_input(
            value["lepton_mass_basis_couplings"]
        ):
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
    if isinstance(value, Mapping) and "lepton_mass_basis_couplings" in value:
        return value["lepton_mass_basis_couplings"]
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


def _flavor_index(flavor: str) -> int:
    return {"e": 0, "mu": 1, "tau": 2}[_canonical_flavor(flavor)]


def _has_rs_ew_couplings(value: Any) -> bool:
    return _rs_ew_couplings(value) is not None


def _rs_ew_couplings(value: Any) -> Any | None:
    if _looks_like_rs_ew_couplings(value):
        return value
    if isinstance(value, Mapping):
        candidate = value.get("rs_ew_couplings")
        return candidate if _looks_like_rs_ew_couplings(candidate) else None
    candidate = getattr(value, "rs_ew_couplings", None)
    return candidate if _looks_like_rs_ew_couplings(candidate) else None


def _looks_like_rs_ew_couplings(value: Any) -> bool:
    return all(
        hasattr(value, name)
        for name in (
            "z_delta_g_L_e",
            "z_delta_g_R_e",
            "z_total_g_L_e",
            "z_total_g_R_e",
        )
    )


def _rigorous_tree_contact_amplitudes(
    source: Any,
    *,
    initial_flavor: str,
    final_flavor: str,
    m_kk_gev: float | None,
    inputs: LFVThreeBodySMInputs | None,
) -> LFVThreeBodyContactAmplitudes:
    del inputs
    couplings = _rs_ew_couplings(source)
    if couplings is None:
        raise TypeError("rs_ew_couplings tree-contact source is absent")

    initial = _canonical_flavor(_first_present(source, _INITIAL_FLAVOR_KEYS) or initial_flavor)
    final = _canonical_flavor(_first_present(source, _FINAL_FLAVOR_KEYS) or final_flavor)
    if initial == final:
        raise ValueError("LFV three-body decay requires distinct lepton flavors")
    i = _flavor_index(initial)
    f = _flavor_index(final)

    z_delta_l = _finite_matrix(couplings.z_delta_g_L_e, "z_delta_g_L_e")
    z_delta_r = _finite_matrix(couplings.z_delta_g_R_e, "z_delta_g_R_e")
    z_total_l = _finite_matrix(couplings.z_total_g_L_e, "z_total_g_L_e")
    z_total_r = _finite_matrix(couplings.z_total_g_R_e, "z_total_g_R_e")

    delta_left = complex(z_delta_l[f, i])
    delta_right = complex(z_delta_r[f, i])
    final_left = complex(z_total_l[f, f])
    final_right = complex(z_total_r[f, f])
    z_ll = complex(2.0 * delta_left * final_left)
    z_lr = complex(2.0 * delta_left * final_right)
    z_rl = complex(2.0 * delta_right * final_left)
    z_rr = complex(2.0 * delta_right * final_right)
    boxes = _box_amplitudes(source)
    m_kk = _resolved_tree_m_kk(couplings, source, m_kk_gev)
    return _contact_amplitude_bundle(
        initial=initial,
        final=final,
        m_kk_gev=m_kk,
        matching_assumption=LFV_THREE_BODY_TREE_CONTACT_RIGOROUS_V1,
        left_lfv_overlap=delta_left,
        right_lfv_overlap=delta_right,
        delta_left=delta_left,
        delta_right=delta_right,
        z_ll=z_ll,
        z_lr=z_lr,
        z_rl=z_rl,
        z_rr=z_rr,
        box_ll=boxes["box_ll"],
        box_lr=boxes["box_lr"],
        box_rl=boxes["box_rl"],
        box_rr=boxes["box_rr"],
        source=_source_label(source) or "rs_ew_couplings charged-lepton Z matrices",
        diagnostics={
            "tree_contact_rigorous": True,
            "tree_contact_source": "rs_ew_couplings",
            "tree_contact_input_path": (
                "z_delta_g_{L,R}_e[final,initial] x "
                "z_total_g_{L,R}_e[final,final]"
            ),
            "tree_contact_zero": _all_zero(z_ll, z_lr, z_rl, z_rr),
            "tree_contact_zero_for_diagonal_fit": _all_zero(z_ll, z_lr, z_rl, z_rr),
            "left_lfv_z_delta_g": delta_left,
            "right_lfv_z_delta_g": delta_right,
            "z_final_lepton_total_g_left": final_left,
            "z_final_lepton_total_g_right": final_right,
            "box_matching_status": (
                "explicit_box_amplitudes_or_zero_NEEDS-HUMAN-PHYSICS"
            ),
            "deferred_matching_needs_human_physics": (
                LFV_THREE_BODY_DEFERRED_PIECES_V1
            ),
        },
    )


def _zero_missing_tree_contact_amplitudes(
    source: Any,
    *,
    initial_flavor: str,
    final_flavor: str,
    m_kk_gev: float | None,
) -> LFVThreeBodyContactAmplitudes:
    initial = _canonical_flavor(_first_present(source, _INITIAL_FLAVOR_KEYS) or initial_flavor)
    final = _canonical_flavor(_first_present(source, _FINAL_FLAVOR_KEYS) or final_flavor)
    m_kk = _positive_float(
        _first_not_none(
            m_kk_gev,
            _first_present(source, ("m_kk_gev", "M_KK_gev", "M_KK")),
            _ZERO_CONTACT_REFERENCE_SCALE_GEV,
        ),
        "m_kk_gev",
    )
    boxes = _box_amplitudes(source)
    return _contact_amplitude_bundle(
        initial=initial,
        final=final,
        m_kk_gev=m_kk,
        matching_assumption=_TREE_CONTACT_MISSING_ASSUMPTION,
        left_lfv_overlap=0.0j,
        right_lfv_overlap=0.0j,
        delta_left=0.0j,
        delta_right=0.0j,
        z_ll=0.0j,
        z_lr=0.0j,
        z_rl=0.0j,
        z_rr=0.0j,
        box_ll=boxes["box_ll"],
        box_lr=boxes["box_lr"],
        box_rl=boxes["box_rl"],
        box_rr=boxes["box_rr"],
        source="zero tree-contact placeholder; rs_ew_couplings absent",
        diagnostics={
            "tree_contact_rigorous": False,
            "tree_contact_source": "absent",
            "tree_contact_missing_extra": "rs_ew_couplings",
            "tree_contact_zero": True,
            "legacy_overlap_tree_proxy_ignored": _has_legacy_tree_overlap(source),
            "box_matching_status": (
                "explicit_box_amplitudes_or_zero_NEEDS-HUMAN-PHYSICS"
            ),
            "missing_extra": "rs_ew_couplings",
            "deferred_matching_needs_human_physics": (
                LFV_THREE_BODY_DEFERRED_PIECES_V1
            ),
        },
    )


def _contact_amplitude_bundle(
    *,
    initial: str,
    final: str,
    m_kk_gev: float,
    matching_assumption: str,
    left_lfv_overlap: complex,
    right_lfv_overlap: complex,
    delta_left: complex,
    delta_right: complex,
    z_ll: complex,
    z_lr: complex,
    z_rl: complex,
    z_rr: complex,
    box_ll: complex,
    box_lr: complex,
    box_rl: complex,
    box_rr: complex,
    source: str,
    diagnostics: Mapping[str, Any],
) -> LFVThreeBodyContactAmplitudes:
    total_ll = complex(z_ll + box_ll)
    total_lr = complex(z_lr + box_lr)
    total_rl = complex(z_rl + box_rl)
    total_rr = complex(z_rr + box_rr)
    full_diagnostics = {
        "m_kk_gev": float(m_kk_gev),
        "matching_scale_gev": float(m_kk_gev),
        "scale_factor": 1.0,
        "left_lfv_overlap": complex(left_lfv_overlap),
        "right_lfv_overlap": complex(right_lfv_overlap),
        "delta_g_left_lfv": complex(delta_left),
        "delta_g_right_lfv": complex(delta_right),
        "z_ll": complex(z_ll),
        "z_lr": complex(z_lr),
        "z_rl": complex(z_rl),
        "z_rr": complex(z_rr),
        "box_ll": complex(box_ll),
        "box_lr": complex(box_lr),
        "box_rl": complex(box_rl),
        "box_rr": complex(box_rr),
        "total_ll": total_ll,
        "total_lr": total_lr,
        "total_rl": total_rl,
        "total_rr": total_rr,
        "proxy_source": str(source),
        "matching_assumption": str(matching_assumption),
        **dict(diagnostics),
    }
    return LFVThreeBodyContactAmplitudes(
        model_label=LFV_THREE_BODY_MODEL_V1,
        matching_assumption=str(matching_assumption),
        operator_convention=LFV_THREE_BODY_OPERATOR_CONVENTION,
        initial_flavor=initial,
        final_flavor=final,
        M_KK=float(m_kk_gev),
        matching_scale=float(m_kk_gev),
        scale_factor=1.0,
        left_lfv_overlap=complex(left_lfv_overlap),
        right_lfv_overlap=complex(right_lfv_overlap),
        delta_g_left_lfv=complex(delta_left),
        delta_g_right_lfv=complex(delta_right),
        z_ll=complex(z_ll),
        z_lr=complex(z_lr),
        z_rl=complex(z_rl),
        z_rr=complex(z_rr),
        box_ll=complex(box_ll),
        box_lr=complex(box_lr),
        box_rl=complex(box_rl),
        box_rr=complex(box_rr),
        total_ll=total_ll,
        total_lr=total_lr,
        total_rl=total_rl,
        total_rr=total_rr,
        source=str(source),
        diagnostics=full_diagnostics,
    )


def _box_amplitudes(source: Any) -> dict[str, complex]:
    return {
        "box_ll": _optional_complex(_first_present(source, _BOX_LL_KEYS), "box_ll") or 0.0j,
        "box_lr": _optional_complex(_first_present(source, _BOX_LR_KEYS), "box_lr") or 0.0j,
        "box_rl": _optional_complex(_first_present(source, _BOX_RL_KEYS), "box_rl") or 0.0j,
        "box_rr": _optional_complex(_first_present(source, _BOX_RR_KEYS), "box_rr") or 0.0j,
    }


def _has_box_amplitude_input(source: Any) -> bool:
    keys = (*_BOX_LL_KEYS, *_BOX_LR_KEYS, *_BOX_RL_KEYS, *_BOX_RR_KEYS)
    if isinstance(source, Mapping):
        if any(_is_live_complex(source[key]) for key in keys if key in source):
            return True
        nested = source.get("lepton_mass_basis_couplings")
        return bool(
            nested is not None
            and nested is not source
            and _has_box_amplitude_input(nested)
        )
    if any(
        _is_live_complex(getattr(source, key))
        for key in keys
        if hasattr(source, key)
    ):
        return True
    nested = getattr(source, "lepton_mass_basis_couplings", None)
    return bool(
        nested is not None
        and nested is not source
        and _has_box_amplitude_input(nested)
    )


def _is_live_complex(value: Any) -> bool:
    if value is None:
        return False
    try:
        number = complex(value)
    except (TypeError, ValueError):
        return True
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        return True
    return abs(number) > 0.0


def _has_legacy_tree_overlap(source: Any) -> bool:
    keys = (
        "left_lfv_overlap",
        "right_lfv_overlap",
        "left_emu_overlap",
        "right_emu_overlap",
        "left_emu",
        "right_emu",
        "left_charged_lepton_overlap",
        "right_charged_lepton_overlap",
        "left_lepton_overlap",
        "right_lepton_overlap",
        "left_overlap",
        "right_overlap",
    )
    if isinstance(source, Mapping):
        return any(key in source for key in keys)
    return any(hasattr(source, key) for key in keys)


def _resolved_tree_m_kk(couplings: Any, source: Any, override: float | None) -> float:
    return _positive_float(
        _first_not_none(
            override,
            getattr(couplings, "kk_ew_mass_gev", None),
            _first_present(source, ("m_kk_gev", "M_KK_gev", "M_KK")),
            _ZERO_CONTACT_REFERENCE_SCALE_GEV,
        ),
        "m_kk_gev",
    )


def _first_present(value: Any, keys: tuple[str, ...]) -> Any:
    if isinstance(value, Mapping):
        for key in keys:
            if key in value:
                return value[key]
        nested = value.get("lepton_mass_basis_couplings")
        if nested is not None and nested is not value:
            return _first_present(nested, keys)
        return None
    for key in keys:
        if hasattr(value, key):
            return getattr(value, key)
    return None


def _first_not_none(*values: Any) -> Any:
    for value in values:
        if value is not None:
            return value
    return None


def _source_label(value: Any) -> str | None:
    raw = _first_present(value, ("source",))
    return None if raw is None else str(raw)


def _finite_matrix(value: Any, name: str) -> np.ndarray:
    matrix = np.asarray(value, dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(matrix.real)) or not np.all(np.isfinite(matrix.imag)):
        raise ValueError(f"{name} entries must be finite")
    return matrix


def _positive_float(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _optional_complex(value: Any, name: str) -> complex | None:
    if value is None:
        return None
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    return number


def _all_zero(*values: complex) -> bool:
    return all(abs(complex(value)) <= 1.0e-30 for value in values)
