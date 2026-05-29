"""Adapter for coherent mu->e conversion constraints.

Constraint modules import this adapter only.  It composes the L001 lepton
dipole adapter with the reusable :mod:`quarkConstraints.mu_e_conversion`
overlap-integral core.

NEEDS-HUMAN-PHYSICS: the current ``ParameterPoint`` does not carry full
charged-lepton scalar/vector lepton-quark RS matching data.  Scalar/vector
terms are therefore explicit low-energy coefficient proxies, flagged in the
result diagnostics.
"""

from __future__ import annotations

from typing import Any, Mapping

from quarkConstraints.mu_e_conversion import (
    MU_E_CONVERSION_DIPOLE_CONVENTION,
    MU_E_CONVERSION_INPUT_BUNDLE_V1,
    MU_E_CONVERSION_MODEL_V1,
    MU_E_CONVERSION_OPERATOR_CONVENTION,
    MU_E_CONVERSION_PROXY_V1,
    MuEConversionCoefficientProxyInput,
    MuEConversionNuclearInputs,
    MuEConversionResult,
    aluminum_nuclear_inputs as _aluminum_nuclear_inputs,
    gold_nuclear_inputs as _gold_nuclear_inputs,
    mu_e_conversion_coefficients as _coefficients,
    mu_e_conversion_from_components as _from_components,
    mu_e_conversion_has_coefficient_proxy as _has_coefficient_proxy,
    mu_e_conversion_proxy_input as _proxy_input,
    nuclear_inputs_for_target as _nuclear_inputs_for_target,
    titanium_nuclear_inputs as _titanium_nuclear_inputs,
    zero_mu_e_conversion_coefficients as _zero_coefficients,
)

from .lepton import mu_to_e_gamma_from_lepton_input

__all__ = [
    "MU_E_CONVERSION_MODEL_V1",
    "MU_E_CONVERSION_INPUT_BUNDLE_V1",
    "MU_E_CONVERSION_OPERATOR_CONVENTION",
    "MU_E_CONVERSION_DIPOLE_CONVENTION",
    "MU_E_CONVERSION_PROXY_V1",
    "MuEConversionCoefficientProxyInput",
    "MuEConversionNuclearInputs",
    "MuEConversionResult",
    "mu_e_conversion_aluminum_nuclear_inputs",
    "mu_e_conversion_titanium_nuclear_inputs",
    "mu_e_conversion_gold_nuclear_inputs",
    "mu_e_conversion_nuclear_inputs_for_target",
    "mu_e_conversion_proxy_input",
    "mu_e_conversion_from_lepton_input",
]


def mu_e_conversion_aluminum_nuclear_inputs() -> MuEConversionNuclearInputs:
    """Return bundled aluminum nuclear inputs."""

    return _aluminum_nuclear_inputs()


def mu_e_conversion_titanium_nuclear_inputs() -> MuEConversionNuclearInputs:
    """Return bundled titanium nuclear inputs."""

    return _titanium_nuclear_inputs()


def mu_e_conversion_gold_nuclear_inputs() -> MuEConversionNuclearInputs:
    """Return bundled gold nuclear inputs."""

    return _gold_nuclear_inputs()


def mu_e_conversion_nuclear_inputs_for_target(
    target: str,
) -> MuEConversionNuclearInputs:
    """Return bundled nuclear inputs for a supported target."""

    return _nuclear_inputs_for_target(target)


def mu_e_conversion_proxy_input(
    *,
    g_lv_p: complex = 0.0j,
    g_lv_n: complex = 0.0j,
    g_rv_p: complex = 0.0j,
    g_rv_n: complex = 0.0j,
    g_ls_p: complex = 0.0j,
    g_ls_n: complex = 0.0j,
    g_rs_p: complex = 0.0j,
    g_rs_n: complex = 0.0j,
    dipole_amplitude_left: complex | None = None,
    dipole_amplitude_right: complex | None = None,
    m_kk_gev: float | None = None,
    source: str = "caller-supplied mu-e conversion coefficient proxy",
) -> MuEConversionCoefficientProxyInput:
    """Build a coefficient proxy accepted by this adapter."""

    return _proxy_input(
        g_lv_p=g_lv_p,
        g_lv_n=g_lv_n,
        g_rv_p=g_rv_p,
        g_rv_n=g_rv_n,
        g_ls_p=g_ls_p,
        g_ls_n=g_ls_n,
        g_rs_p=g_rs_p,
        g_rs_n=g_rs_n,
        dipole_amplitude_left=dipole_amplitude_left,
        dipole_amplitude_right=dipole_amplitude_right,
        m_kk_gev=m_kk_gev,
        source=source,
    )


def mu_e_conversion_from_lepton_input(
    lepton_input: Any,
    *,
    conversion_rate_limit: float,
    dipole_br_limit: float,
    dipole_prefactor_br: float,
    reference_scale_gev: float = 3000.0,
    m_kk_gev: float | None = None,
    nuclear_inputs: MuEConversionNuclearInputs | None = None,
) -> MuEConversionResult:
    """Evaluate ``CR(mu -> e, target)`` from dipole and coefficient proxies."""

    has_dipole = _has_dipole_input(lepton_input)
    has_coefficients = _has_coefficient_proxy(lepton_input)
    if not has_dipole and not has_coefficients:
        raise TypeError(
            "lepton input must provide a dipole source and/or mu-e conversion "
            "coefficient proxy inputs"
        )

    dipole_parent_br = 0.0
    resolved_m_kk = m_kk_gev
    dipole_diagnostics: dict[str, Any] = {
        "dipole_input_present": bool(has_dipole),
        "dipole_component_reuses": "flavor_catalog_constraints.physics_adapters.lepton",
    }
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
            if not has_coefficients:
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

    if has_coefficients:
        coefficients = _coefficients(lepton_input)
    else:
        coefficients = _zero_coefficients(source="dipole-only mu-e conversion proxy")
    if coefficients.m_kk_gev is None and resolved_m_kk is not None:
        coefficients = _proxy_input(
            g_lv_p=coefficients.g_lv_p,
            g_lv_n=coefficients.g_lv_n,
            g_rv_p=coefficients.g_rv_p,
            g_rv_n=coefficients.g_rv_n,
            g_ls_p=coefficients.g_ls_p,
            g_ls_n=coefficients.g_ls_n,
            g_rs_p=coefficients.g_rs_p,
            g_rs_n=coefficients.g_rs_n,
            dipole_amplitude_left=coefficients.dipole_amplitude_left,
            dipole_amplitude_right=coefficients.dipole_amplitude_right,
            m_kk_gev=resolved_m_kk,
            source=coefficients.source,
        )

    result = _from_components(
        dipole_parent_branching_fraction=dipole_parent_br,
        coefficients=coefficients,
        conversion_rate_limit=conversion_rate_limit,
        nuclear_inputs=nuclear_inputs,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics.update(dipole_diagnostics)
    return MuEConversionResult(
        model_label=result.model_label,
        input_bundle=result.input_bundle,
        target=result.target,
        conversion_rate=result.conversion_rate,
        conversion_rate_lower=result.conversion_rate_lower,
        conversion_rate_upper=result.conversion_rate_upper,
        sm_conversion_rate=result.sm_conversion_rate,
        np_shift_conversion_rate=result.np_shift_conversion_rate,
        conversion_width_gev=result.conversion_width_gev,
        conversion_width_lower_gev=result.conversion_width_lower_gev,
        conversion_width_upper_gev=result.conversion_width_upper_gev,
        capture_rate_gev=result.capture_rate_gev,
        dipole_component=result.dipole_component,
        scalar_component=result.scalar_component,
        vector_component=result.vector_component,
        vector_scalar_interference_component=(
            result.vector_scalar_interference_component
        ),
        contact_component=result.contact_component,
        dipole_contact_interference_component=(
            result.dipole_contact_interference_component
        ),
        dipole_contact_interference_lower=(
            result.dipole_contact_interference_lower
        ),
        dipole_contact_interference_upper=(
            result.dipole_contact_interference_upper
        ),
        dipole_contact_interference_treatment=(
            result.dipole_contact_interference_treatment
        ),
        dipole_parent_branching_fraction=result.dipole_parent_branching_fraction,
        dipole_amplitude_left=result.dipole_amplitude_left,
        dipole_amplitude_right=result.dipole_amplitude_right,
        left_nuclear_amplitude=result.left_nuclear_amplitude,
        right_nuclear_amplitude=result.right_nuclear_amplitude,
        ratio_to_limit=result.ratio_to_limit,
        ratio_to_limit_lower=result.ratio_to_limit_lower,
        ratio_to_limit_upper=result.ratio_to_limit_upper,
        conversion_rate_limit=result.conversion_rate_limit,
        passes=result.passes,
        coefficients=result.coefficients,
        nuclear_inputs=result.nuclear_inputs,
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
