"""Tau-electron adapter for charged-LFV three-body decays.

This adapter pins the shared :mod:`quarkConstraints.lfv_three_body` machinery
to ``tau -> 3e`` for L010.  Constraint modules import this wrapper only; the
lower-level physics module remains behind the adapter boundary.

NEEDS-HUMAN-PHYSICS: a rigorous RS prediction for ``tau -> 3e`` requires
charged-lepton neutral-current, EW KK/Z/Z', loop-level dipole, and four-lepton
box matching inputs that are not present on ``ParameterPoint``.  The contact
terms reuse the documented L002 proxy convention, and any dipole component is
accepted only as an explicit caller-supplied ``tau -> e gamma`` branching
fraction or chiral dipole amplitude in the low-energy convention.
"""

from __future__ import annotations

import math
from dataclasses import replace
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
    lfv_three_body_contact_amplitudes as _contact_amplitudes,
)
from quarkConstraints.lfv_three_body import (
    lfv_three_body_from_components as _from_components,
)
from quarkConstraints.lfv_three_body import (
    lfv_three_body_has_contact_proxy as _has_contact_proxy,
)
from quarkConstraints.lfv_three_body import (
    lfv_three_body_proxy_input as _proxy_input,
)

__all__ = [
    "LFV_THREE_BODY_INPUT_BUNDLE_V1",
    "LFV_THREE_BODY_MODEL_V1",
    "LFV_THREE_BODY_OPERATOR_CONVENTION",
    "LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION",
    "LFV_THREE_BODY_PROXY_V1",
    "TAU_TO_3E_PROXY_V1",
    "LFVThreeBodySMInputs",
    "LFVThreeBodyContactProxyInput",
    "LFVThreeBodyContactAmplitudes",
    "LFVThreeBodyBranchingResult",
    "tau_to_3e_default_sm_inputs",
    "tau_to_3e_proxy_input",
    "tau_to_3e_from_lepton_input",
]

TAU_TO_3E_PROXY_V1 = (
    "NEEDS-HUMAN-PHYSICS: tau->3e reuses the shared l_i->3l_j "
    "Z-penguin/box contact proxy with initial_flavor='tau' and "
    "final_flavor='e'.  The tau->e gamma dipole contribution is included "
    "only when the caller supplies an explicit dipole parent branching "
    "fraction or chiral dipole amplitudes; missing chiral phase information "
    "uses the constructive sign-envelope convention."
)

_INITIAL_FLAVOR = "tau"
_FINAL_FLAVOR = "e"
_ZERO_CONTACT_REFERENCE_SCALE_GEV = 3000.0
_INITIAL_FLAVOR_KEYS = ("initial_flavor", "parent_flavor")
_FINAL_FLAVOR_KEYS = ("final_flavor", "daughter_flavor")
_LEFT_TAUE_OVERLAP_KEYS = (
    "left_lfv_overlap",
    "left_taue_overlap",
    "left_tau_e_overlap",
    "left_etau_overlap",
    "left_e_tau_overlap",
)
_RIGHT_TAUE_OVERLAP_KEYS = (
    "right_lfv_overlap",
    "right_taue_overlap",
    "right_tau_e_overlap",
    "right_etau_overlap",
    "right_e_tau_overlap",
)
_MISMATCHED_SPURION_KEYS = (
    "left_emu_overlap",
    "right_emu_overlap",
    "left_emu",
    "right_emu",
    "left_mue_overlap",
    "right_mue_overlap",
    "left_mu_e_overlap",
    "right_mu_e_overlap",
    "left_e_mu_overlap",
    "right_e_mu_overlap",
    "left_taumu_overlap",
    "right_taumu_overlap",
    "left_mutau_overlap",
    "right_mutau_overlap",
    "left_tau_mu_overlap",
    "right_tau_mu_overlap",
    "left_mu_tau_overlap",
    "right_mu_tau_overlap",
)
_LEFT_MATRIX_KEYS = (
    "left_charged_lepton_overlap",
    "left_lepton_overlap",
    "left_overlap",
)
_RIGHT_MATRIX_KEYS = (
    "right_charged_lepton_overlap",
    "right_lepton_overlap",
    "right_overlap",
)
_BOX_LL_KEYS = ("box_ll", "box_left_left", "box_L_L")
_BOX_LR_KEYS = ("box_lr", "box_left_right", "box_L_R")
_BOX_RL_KEYS = ("box_rl", "box_right_left", "box_R_L")
_BOX_RR_KEYS = ("box_rr", "box_right_right", "box_R_R")
_DIPOLE_BR_KEYS = (
    "dipole_parent_branching_fraction",
    "tau_to_e_gamma_branching_fraction",
    "br_tau_to_e_gamma",
    "br_l_to_lgamma",
)
_DIPOLE_LEFT_KEYS = ("dipole_amplitude_left", "a_left", "A_L")
_DIPOLE_RIGHT_KEYS = ("dipole_amplitude_right", "a_right", "A_R")
_M_KK_KEYS = ("m_kk_gev", "M_KK_gev", "M_KK")
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


def tau_to_3e_default_sm_inputs() -> LFVThreeBodySMInputs:
    """Return the default low-energy input bundle."""

    return _default_sm_inputs()


def tau_to_3e_proxy_input(
    left_taue_overlap: complex,
    right_taue_overlap: complex,
    m_kk_gev: float,
    *,
    box_ll: complex = 0.0j,
    box_lr: complex = 0.0j,
    box_rl: complex = 0.0j,
    box_rr: complex = 0.0j,
    source: str = "caller-supplied tau->3e contact proxy",
) -> LFVThreeBodyContactProxyInput:
    """Build a contact proxy pinned to ``tau -> 3e``."""

    return _proxy_input(
        left_taue_overlap,
        right_taue_overlap,
        m_kk_gev,
        initial_flavor=_INITIAL_FLAVOR,
        final_flavor=_FINAL_FLAVOR,
        box_ll=box_ll,
        box_lr=box_lr,
        box_rl=box_rl,
        box_rr=box_rr,
        source=source,
    )


def tau_to_3e_from_lepton_input(
    lepton_input: Any,
    *,
    br_limit: float,
    m_kk_gev: float | None = None,
    inputs: LFVThreeBodySMInputs | None = None,
) -> LFVThreeBodyBranchingResult:
    """Evaluate ``BR(tau -> 3e)`` from explicit dipole/contact proxies."""

    p = _default_sm_inputs() if inputs is None else inputs
    _assert_tau_to_3e_flavors(lepton_input)
    _assert_no_mismatched_spurion_aliases(lepton_input)
    has_contact = _has_tau_to_3e_contact_proxy(lepton_input)
    dipole = _coerce_dipole_proxy(lepton_input)
    if not has_contact and not dipole["present"]:
        raise TypeError(
            "lepton input must provide tau->3e contact proxy inputs and/or "
            "an explicit tau->e gamma dipole branching fraction"
        )

    if has_contact:
        contact_proxy = _tau_e_contact_proxy(lepton_input, m_kk_gev=m_kk_gev)
        contact = _contact_amplitudes(
            contact_proxy,
            initial_flavor=_INITIAL_FLAVOR,
            final_flavor=_FINAL_FLAVOR,
            inputs=p,
        )
    else:
        contact = _zero_contact_amplitudes(lepton_input, m_kk_gev=m_kk_gev, inputs=p)

    result = _from_components(
        dipole_parent_branching_fraction=float(dipole["branching_fraction"]),
        contact_amplitudes=contact,
        br_limit=br_limit,
        inputs=p,
        dipole_amplitude_left=dipole["amplitude_left"],
        dipole_amplitude_right=dipole["amplitude_right"],
    )
    diagnostics = dict(result.diagnostics)
    diagnostics.update(
        {
            "tau_to_3e_proxy": TAU_TO_3E_PROXY_V1,
            "pinned_initial_flavor": _INITIAL_FLAVOR,
            "pinned_final_flavor": _FINAL_FLAVOR,
            "contact_input_present": bool(has_contact),
            "contact_input_normalized_to_tau_e": bool(has_contact),
            "dipole_input_present": bool(dipole["present"]),
            "explicit_tau_to_e_gamma_branching_fraction": float(
                dipole["branching_fraction"]
            ),
            "dipole_input_kind": dipole["input_kind"],
            "dipole_source": dipole["source"],
            "dipole_amplitudes_supplied": bool(dipole["amplitudes_supplied"]),
            "dipole_component_needs_human_physics": bool(dipole["present"]),
        }
    )
    return replace(result, diagnostics=diagnostics)


def _has_tau_to_3e_contact_proxy(source: Any) -> bool:
    if _has_contact_proxy(source):
        return True
    return any(
        _has_key_or_attr(source, key)
        for key in (*_LEFT_TAUE_OVERLAP_KEYS, *_RIGHT_TAUE_OVERLAP_KEYS)
    )


def _zero_contact_amplitudes(
    lepton_input: Any,
    *,
    m_kk_gev: float | None,
    inputs: LFVThreeBodySMInputs,
) -> LFVThreeBodyContactAmplitudes:
    m_kk = _resolve_m_kk(
        lepton_input,
        m_kk_gev,
        default=_ZERO_CONTACT_REFERENCE_SCALE_GEV,
    )
    return _contact_amplitudes(
        _proxy_input(
            0.0j,
            0.0j,
            m_kk,
            initial_flavor=_INITIAL_FLAVOR,
            final_flavor=_FINAL_FLAVOR,
            source="zero contact placeholder for explicit tau->e gamma dipole",
        ),
        initial_flavor=_INITIAL_FLAVOR,
        final_flavor=_FINAL_FLAVOR,
        inputs=inputs,
    )


def _tau_e_contact_proxy(
    source: Any,
    *,
    m_kk_gev: float | None,
) -> LFVThreeBodyContactProxyInput:
    """Normalize contact inputs before the shared core can apply mu->e defaults."""

    if isinstance(source, LFVThreeBodyContactProxyInput):
        return _proxy_input(
            source.left_lfv_overlap,
            source.right_lfv_overlap,
            source.m_kk_gev if m_kk_gev is None else m_kk_gev,
            initial_flavor=_INITIAL_FLAVOR,
            final_flavor=_FINAL_FLAVOR,
            box_ll=source.box_ll,
            box_lr=source.box_lr,
            box_rl=source.box_rl,
            box_rr=source.box_rr,
            source=source.source,
        )

    left, right = _tau_e_overlaps(source)
    return _proxy_input(
        0.0j if left is None else left,
        0.0j if right is None else right,
        _resolve_m_kk(source, m_kk_gev),
        initial_flavor=_INITIAL_FLAVOR,
        final_flavor=_FINAL_FLAVOR,
        box_ll=_optional_complex(_first_present(source, _BOX_LL_KEYS), "box_ll") or 0.0j,
        box_lr=_optional_complex(_first_present(source, _BOX_LR_KEYS), "box_lr") or 0.0j,
        box_rl=_optional_complex(_first_present(source, _BOX_RL_KEYS), "box_rl") or 0.0j,
        box_rr=_optional_complex(_first_present(source, _BOX_RR_KEYS), "box_rr") or 0.0j,
        source=_source_label(source) or "tau->3e contact proxy normalized by adapter",
    )


def _assert_no_mismatched_spurion_aliases(value: Any) -> None:
    present = [key for key in _MISMATCHED_SPURION_KEYS if _has_key_or_attr(value, key)]
    if present:
        keys = ", ".join(present)
        raise ValueError(
            "tau->3e does not accept mu-e or tau-mu overlap aliases "
            f"({keys}); use tau-e generic LFV overlaps or a flavor matrix"
        )


def _tau_e_overlaps(value: Any) -> tuple[complex | None, complex | None]:
    left = _first_present(value, _LEFT_TAUE_OVERLAP_KEYS)
    right = _first_present(value, _RIGHT_TAUE_OVERLAP_KEYS)
    if left is not None or right is not None:
        return (
            _optional_complex(left, "left_lfv_overlap"),
            _optional_complex(right, "right_lfv_overlap"),
        )

    left_matrix = _first_present(value, _LEFT_MATRIX_KEYS)
    right_matrix = _first_present(value, _RIGHT_MATRIX_KEYS)
    if left_matrix is not None:
        left = _matrix_offdiag(left_matrix, "left")
    if right_matrix is not None:
        right = _matrix_offdiag(right_matrix, "right")
    return (
        None if left is None else complex(left),
        None if right is None else complex(right),
    )


def _matrix_offdiag(value: Any, name: str) -> complex:
    matrix = np.asarray(value, dtype=np.complex128)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} charged-lepton overlap matrix must have shape (3, 3)")
    if not np.all(np.isfinite(matrix.real)) or not np.all(np.isfinite(matrix.imag)):
        raise ValueError(f"{name} charged-lepton overlap matrix entries must be finite")
    return complex(matrix[_flavor_index(_FINAL_FLAVOR), _flavor_index(_INITIAL_FLAVOR)])


def _coerce_dipole_proxy(source: Any) -> dict[str, Any]:
    root = source
    nested = _nested_dipole_source(source)
    if nested is not None:
        root = nested
        _assert_tau_to_3e_flavors(root)

    left = _optional_complex(_first_present(root, _DIPOLE_LEFT_KEYS), "dipole_amplitude_left")
    right = _optional_complex(
        _first_present(root, _DIPOLE_RIGHT_KEYS),
        "dipole_amplitude_right",
    )
    amplitudes_supplied = left is not None or right is not None
    br_value = _first_present(root, _DIPOLE_BR_KEYS)

    if br_value is None and _is_numeric(root):
        br_value = root

    if br_value is None and amplitudes_supplied:
        left_for_norm = 0.0j if left is None else left
        right_for_norm = 0.0j if right is None else right
        br_value = 384.0 * math.pi**2 * (
            abs(left_for_norm) ** 2 + abs(right_for_norm) ** 2
        )

    present = br_value is not None or amplitudes_supplied
    branching_fraction = 0.0 if br_value is None else _bounded_nonnegative(
        br_value,
        "tau->e gamma dipole branching fraction",
    )

    return {
        "present": bool(present),
        "branching_fraction": float(branching_fraction),
        "amplitude_left": left,
        "amplitude_right": right,
        "amplitudes_supplied": bool(amplitudes_supplied),
        "input_kind": _dipole_input_kind(source, nested),
        "source": _source_label(root),
    }


def _nested_dipole_source(source: Any) -> Any:
    if isinstance(source, Mapping):
        return source.get("dipole")
    return getattr(source, "dipole", None)


def _source_label(value: Any) -> str | None:
    if isinstance(value, Mapping):
        raw = value.get("source")
        return None if raw is None else str(raw)
    raw = getattr(value, "source", None)
    return None if raw is None else str(raw)


def _dipole_input_kind(source: Any, nested: Any) -> str:
    if nested is not None:
        return "nested_dipole_proxy"
    if _is_numeric(source):
        return "numeric_branching_fraction"
    if isinstance(source, Mapping):
        return "mapping"
    return type(source).__name__


def _resolve_m_kk(
    source: Any,
    override: float | None,
    *,
    default: float | None = None,
) -> float:
    if override is not None:
        return _positive_float(override, "m_kk_gev")
    value = _first_present(source, _M_KK_KEYS)
    if value is None:
        if default is not None:
            return _positive_float(default, "m_kk_gev")
        raise KeyError("tau->3e proxy must provide m_kk_gev, M_KK_gev, or M_KK")
    return _positive_float(value, "m_kk_gev")


def _assert_tau_to_3e_flavors(value: Any) -> None:
    initial = _first_explicit_flavor(value, _INITIAL_FLAVOR_KEYS)
    final = _first_explicit_flavor(value, _FINAL_FLAVOR_KEYS)
    if initial is not None and _canonical_flavor(initial) != _INITIAL_FLAVOR:
        raise ValueError(
            "tau_to_3e_from_lepton_input is pinned to initial_flavor='tau'; "
            f"got {initial!r}"
        )
    if final is not None and _canonical_flavor(final) != _FINAL_FLAVOR:
        raise ValueError(
            "tau_to_3e_from_lepton_input is pinned to final_flavor='e'; "
            f"got {final!r}"
        )


def _first_explicit_flavor(value: Any, keys: tuple[str, ...]) -> Any:
    if isinstance(value, LFVThreeBodyContactProxyInput):
        if keys == _INITIAL_FLAVOR_KEYS:
            return value.initial_flavor
        return value.final_flavor
    return _first_present(value, keys)


def _first_present(value: Any, keys: tuple[str, ...]) -> Any:
    if isinstance(value, Mapping):
        for key in keys:
            if key in value:
                return value[key]
        return None
    for key in keys:
        if hasattr(value, key):
            return getattr(value, key)
    return None


def _has_key_or_attr(value: Any, key: str) -> bool:
    if isinstance(value, Mapping):
        return key in value
    return hasattr(value, key)


def _canonical_flavor(value: Any) -> str:
    key = _CHARGED_LEPTON_ALIASES.get(str(value), str(value))
    if key not in _CHARGED_LEPTONS:
        raise ValueError(f"unsupported charged-lepton flavor {value!r}")
    return key


def _flavor_index(flavor: str) -> int:
    return {"e": 0, "mu": 1, "tau": 2}[_canonical_flavor(flavor)]


def _positive_float(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _bounded_nonnegative(value: Any, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0 or number >= 1.0:
        raise ValueError(f"{name} must be finite and satisfy 0 <= value < 1")
    return number


def _optional_complex(value: Any, name: str) -> complex | None:
    if value is None:
        return None
    number = complex(value)
    if not math.isfinite(number.real) or not math.isfinite(number.imag):
        raise ValueError(f"{name} must be finite")
    return number


def _is_numeric(value: Any) -> bool:
    return isinstance(value, (int, float, complex)) and not isinstance(value, bool)
