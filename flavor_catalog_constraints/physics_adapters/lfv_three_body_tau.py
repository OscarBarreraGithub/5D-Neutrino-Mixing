"""Tau-specific adapter for charged-LFV three-body decays.

This adapter pins the shared :mod:`quarkConstraints.lfv_three_body` machinery
to ``tau -> 3mu`` for L009.  Constraint modules import this wrapper only; the
lower-level physics module remains behind the adapter boundary.

Tree-level light-Z contacts are read from the Phase-4a ``rs_ew_couplings``
charged-lepton Z matrices when available.  The tau dipole, dipole-contact
phase when only a parent dipole rate is supplied, heavy neutral exchange, and
four-lepton boxes remain deferred.
"""

from __future__ import annotations

import math
from dataclasses import replace
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

from .lfv_three_body import (
    LFV_THREE_BODY_DEFERRED_PIECES_V1,
    LFV_THREE_BODY_TREE_CONTACT_RIGOROUS_V1,
)
from .lfv_three_body import (
    lfv_three_body_contact_amplitudes as _adapter_contact_amplitudes,
)
from .lfv_three_body import (
    lfv_three_body_has_tree_contact_input as _has_tree_contact_input,
)

__all__ = [
    "LFV_THREE_BODY_INPUT_BUNDLE_V1",
    "LFV_THREE_BODY_MODEL_V1",
    "LFV_THREE_BODY_OPERATOR_CONVENTION",
    "LFV_THREE_BODY_DIPOLE_CONTACT_INTERFERENCE_CONVENTION",
    "LFV_THREE_BODY_PROXY_V1",
    "LFV_THREE_BODY_TREE_CONTACT_RIGOROUS_V1",
    "LFV_THREE_BODY_DEFERRED_PIECES_V1",
    "TAU_TO_3MU_PROXY_V1",
    "LFVThreeBodySMInputs",
    "LFVThreeBodyContactProxyInput",
    "LFVThreeBodyContactAmplitudes",
    "LFVThreeBodyBranchingResult",
    "tau_to_3mu_default_sm_inputs",
    "tau_to_3mu_proxy_input",
    "tau_to_3mu_from_lepton_input",
]

TAU_TO_3MU_PROXY_V1 = (
    "tree-level light-Z tau->3mu contact is rigorous when rs_ew_couplings is "
    "present and is zero for the Phase-4a diagonal charged-lepton fit; "
    "tau->mu gamma dipole, dipole-contact phase, heavy neutral exchange, and "
    "box matching remain NEEDS-HUMAN-PHYSICS"
)

_INITIAL_FLAVOR = "tau"
_FINAL_FLAVOR = "mu"
_ZERO_CONTACT_REFERENCE_SCALE_GEV = 3000.0
_INITIAL_FLAVOR_KEYS = ("initial_flavor", "parent_flavor")
_FINAL_FLAVOR_KEYS = ("final_flavor", "daughter_flavor")
_EMU_ALIAS_SPURION_KEYS = (
    "left_emu_overlap",
    "right_emu_overlap",
    "left_emu",
    "right_emu",
)
_DIPOLE_BR_KEYS = (
    "dipole_parent_branching_fraction",
    "tau_to_mu_gamma_branching_fraction",
    "br_tau_to_mu_gamma",
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


def tau_to_3mu_default_sm_inputs() -> LFVThreeBodySMInputs:
    """Return the default low-energy input bundle."""

    return _default_sm_inputs()


def tau_to_3mu_proxy_input(
    left_taumu_overlap: complex,
    right_taumu_overlap: complex,
    m_kk_gev: float,
    *,
    box_ll: complex = 0.0j,
    box_lr: complex = 0.0j,
    box_rl: complex = 0.0j,
    box_rr: complex = 0.0j,
    source: str = "caller-supplied tau->3mu contact proxy",
) -> LFVThreeBodyContactProxyInput:
    """Build a contact proxy pinned to ``tau -> 3mu``."""

    return _proxy_input(
        left_taumu_overlap,
        right_taumu_overlap,
        m_kk_gev,
        initial_flavor=_INITIAL_FLAVOR,
        final_flavor=_FINAL_FLAVOR,
        box_ll=box_ll,
        box_lr=box_lr,
        box_rl=box_rl,
        box_rr=box_rr,
        source=source,
    )


def tau_to_3mu_from_lepton_input(
    lepton_input: Any,
    *,
    br_limit: float,
    m_kk_gev: float | None = None,
    inputs: LFVThreeBodySMInputs | None = None,
) -> LFVThreeBodyBranchingResult:
    """Evaluate ``BR(tau -> 3mu)`` from explicit dipole/contact proxies."""

    p = _default_sm_inputs() if inputs is None else inputs
    _assert_tau_to_3mu_flavors(lepton_input)
    _assert_no_emu_alias_spurions(lepton_input)
    has_contact = _has_tree_contact_input(lepton_input)
    dipole = _coerce_dipole_proxy(lepton_input)
    if not has_contact and not dipole["present"]:
        raise TypeError(
            "lepton input must provide tau->3mu contact proxy inputs and/or "
            "an explicit tau->mu gamma dipole branching fraction"
        )

    if has_contact:
        contact = _adapter_contact_amplitudes(
            lepton_input,
            initial_flavor=_INITIAL_FLAVOR,
            final_flavor=_FINAL_FLAVOR,
            m_kk_gev=m_kk_gev,
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
            "tau_to_3mu_proxy": TAU_TO_3MU_PROXY_V1,
            "pinned_initial_flavor": _INITIAL_FLAVOR,
            "pinned_final_flavor": _FINAL_FLAVOR,
            "contact_input_present": bool(has_contact),
            "contact_input_normalized_to_tau_mu": bool(has_contact),
            "tree_contact_matching": contact.matching_assumption,
            "dipole_input_present": bool(dipole["present"]),
            "explicit_tau_to_mu_gamma_branching_fraction": float(
                dipole["branching_fraction"]
            ),
            "dipole_input_kind": dipole["input_kind"],
            "dipole_source": dipole["source"],
            "dipole_amplitudes_supplied": bool(dipole["amplitudes_supplied"]),
            "dipole_component_needs_human_physics": bool(dipole["present"]),
            "deferred_matching_needs_human_physics": (
                LFV_THREE_BODY_DEFERRED_PIECES_V1
            ),
        }
    )
    return replace(result, diagnostics=diagnostics)


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
    return _adapter_contact_amplitudes(
        {
            "m_kk_gev": m_kk,
            "source": "zero tree-contact placeholder; rs_ew_couplings absent",
        },
        initial_flavor=_INITIAL_FLAVOR,
        final_flavor=_FINAL_FLAVOR,
        m_kk_gev=m_kk,
        inputs=inputs,
    )


def _assert_no_emu_alias_spurions(value: Any) -> None:
    present = [key for key in _EMU_ALIAS_SPURION_KEYS if _has_key_or_attr(value, key)]
    if present:
        keys = ", ".join(present)
        raise ValueError(
            "tau->3mu does not accept mu->e/e-mu overlap aliases "
            f"({keys}); use tau-mu generic LFV overlaps or a flavor matrix"
        )


def _coerce_dipole_proxy(source: Any) -> dict[str, Any]:
    root = source
    nested = _nested_dipole_source(source)
    if nested is not None:
        root = nested

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
        "tau->mu gamma dipole branching fraction",
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
        raise KeyError("tau->3mu proxy must provide m_kk_gev, M_KK_gev, or M_KK")
    return _positive_float(value, "m_kk_gev")


def _assert_tau_to_3mu_flavors(value: Any) -> None:
    initial = _first_explicit_flavor(value, _INITIAL_FLAVOR_KEYS)
    final = _first_explicit_flavor(value, _FINAL_FLAVOR_KEYS)
    if initial is not None and _canonical_flavor(initial) != _INITIAL_FLAVOR:
        raise ValueError(
            "tau_to_3mu_from_lepton_input is pinned to initial_flavor='tau'; "
            f"got {initial!r}"
        )
    if final is not None and _canonical_flavor(final) != _FINAL_FLAVOR:
        raise ValueError(
            "tau_to_3mu_from_lepton_input is pinned to final_flavor='mu'; "
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
