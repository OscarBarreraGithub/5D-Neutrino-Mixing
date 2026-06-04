"""RS electroweak ``ParameterPoint`` extras builder."""

from __future__ import annotations

from typing import Any

from quarkConstraints.rs_ew_couplings import (
    RSEWNeutralCurrentInputs,
    build_rs_ew_couplings,
    build_rs_lepton_mass_basis_couplings,
)
from quarkConstraints.rs_charged_current import (
    RSChargedCurrentInputs,
    build_rs_charged_current,
)
from quarkConstraints.rs_ew_spectrum import (
    DEFAULT_N_GAUGE_MODES,
    DEFAULT_OVERLAP_RTOL,
    DEFAULT_QUADRATURE_ORDER,
    RSEWSpectrum,
)
from quarkConstraints.rs_semileptonic_wilsons import (
    build_rs_semileptonic_wilsons,
)
from warpConfig.baseParams import MPL


def build_rs_ew_extras(
    quark_fit_result: Any,
    *,
    Lambda_IR: float,
    k: float = MPL,
    n_gauge_modes: int = DEFAULT_N_GAUGE_MODES,
    ew_model: str = "minimal_rs",
    include_fermion_kk_mixing: bool = False,
    include_loop_dipoles: bool = False,
    lepton_yukawa_result: Any | None = None,
    lepton_sweep_inputs: dict[str, Any] | None = None,
    neutral_current_inputs: RSEWNeutralCurrentInputs | None = None,
    include_charged_current: bool = False,
    charged_current_inputs: RSChargedCurrentInputs | None = None,
    overlap_rel_tol: float = DEFAULT_OVERLAP_RTOL,
    min_overlap_modes: int = 16,
    max_overlap_modes: int | None = None,
    quadrature_order: int = DEFAULT_QUADRATURE_ORDER,
) -> dict[str, Any]:
    """Build RS-EW neutral-current extras without rewiring constraints."""

    if ew_model != "minimal_rs":
        raise ValueError("RS-EW extras currently support ew_model='minimal_rs' only")
    if include_loop_dipoles:
        raise ValueError("loop dipoles are deferred beyond Phase 3a")

    spectrum = RSEWSpectrum.build(
        lambda_ir_gev=float(Lambda_IR),
        k_gev=float(k),
        n_gauge_modes=int(n_gauge_modes),
        quadrature_order=int(quadrature_order),
        model_label=str(ew_model),
    )
    resolved_lepton_yukawas = _resolve_lepton_yukawa_result(
        Lambda_IR=float(Lambda_IR),
        k=float(k),
        lepton_yukawa_result=lepton_yukawa_result,
        lepton_sweep_inputs=lepton_sweep_inputs,
    )
    lepton_couplings = (
        None
        if resolved_lepton_yukawas is None
        else build_rs_lepton_mass_basis_couplings(
            resolved_lepton_yukawas,
            kk_ew_mass_gev=float(spectrum.kk_ew_mass_gev),
        )
    )
    couplings = build_rs_ew_couplings(
        quark_fit_result,
        spectrum=spectrum,
        lepton_mass_basis_couplings=lepton_couplings,
        inputs=neutral_current_inputs,
        include_fermion_kk_mixing=bool(include_fermion_kk_mixing),
        overlap_rel_tol=float(overlap_rel_tol),
        min_overlap_modes=int(min_overlap_modes),
        max_overlap_modes=max_overlap_modes,
        model_label=str(ew_model),
    )
    charged_current = None
    if include_charged_current:
        if lepton_couplings is None:
            raise ValueError(
                "include_charged_current=True requires lepton_yukawa_result "
                "or lepton_sweep_inputs with c_L"
            )
        resolved_charged_inputs = charged_current_inputs
        if resolved_charged_inputs is None:
            resolved_charged_inputs = RSChargedCurrentInputs(
                a_ref_c=float(couplings.a_ref_c)
            )
        if float(resolved_charged_inputs.a_ref_c) != float(couplings.a_ref_c):
            raise ValueError("charged_current_inputs.a_ref_c must match rs_ew_couplings.a_ref_c")
        charged_current = build_rs_charged_current(
            quark_fit_result,
            spectrum=spectrum,
            lepton_mass_basis_couplings=lepton_couplings,
            inputs=resolved_charged_inputs,
            shared_a_ref=float(couplings.a_ref),
            overlap_rel_tol=float(overlap_rel_tol),
            min_overlap_modes=int(min_overlap_modes),
            max_overlap_modes=max_overlap_modes,
            model_label=str(ew_model),
        )
    wilsons = build_rs_semileptonic_wilsons(couplings)
    return {
        "kk_ew_mass_gev": float(spectrum.kk_ew_mass_gev),
        "rs_ew_spectrum": spectrum,
        "rs_ew_couplings": couplings,
        "rs_semileptonic_wilsons": wilsons,
        **(
            {}
            if lepton_couplings is None
            else {"lepton_mass_basis_couplings": lepton_couplings}
        ),
        **({} if charged_current is None else {"rs_charged_current": charged_current}),
    }


def _resolve_lepton_yukawa_result(
    *,
    Lambda_IR: float,
    k: float,
    lepton_yukawa_result: Any | None,
    lepton_sweep_inputs: dict[str, Any] | None,
) -> Any | None:
    if lepton_yukawa_result is not None and lepton_sweep_inputs is not None:
        raise ValueError("provide only one of lepton_yukawa_result or lepton_sweep_inputs")
    if lepton_yukawa_result is not None:
        return lepton_yukawa_result
    if lepton_sweep_inputs is None:
        return None

    from yukawa.compute_yukawas import compute_all_yukawas

    inputs = dict(lepton_sweep_inputs)
    if "M_N" not in inputs:
        if "MN_over_k" not in inputs:
            raise ValueError("lepton_sweep_inputs must include M_N or MN_over_k")
        inputs["M_N"] = float(inputs.pop("MN_over_k")) * float(k)
    inputs.setdefault("Lambda_IR", float(Lambda_IR))
    inputs.setdefault("k", float(k))
    return compute_all_yukawas(
        Lambda_IR=float(inputs["Lambda_IR"]),
        c_L=inputs["c_L"],
        c_E=inputs["c_E"],
        c_N=inputs["c_N"],
        M_N=float(inputs["M_N"]),
        lightest_nu_mass=float(inputs.get("lightest_nu_mass", 0.0)),
        ordering=str(inputs.get("ordering", "normal")),
        majorana_alpha=float(inputs.get("majorana_alpha", 0.0)),
        majorana_beta=float(inputs.get("majorana_beta", 0.0)),
        k=float(inputs["k"]),
        v=float(inputs.get("v", 174.0)),
    )


__all__ = ["build_rs_ew_extras"]
