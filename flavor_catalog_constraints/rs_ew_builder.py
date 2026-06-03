"""Phase-3a RS electroweak ``ParameterPoint`` extras builder."""

from __future__ import annotations

from typing import Any

from quarkConstraints.rs_ew_couplings import (
    RSEWNeutralCurrentInputs,
    build_rs_ew_couplings,
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
    neutral_current_inputs: RSEWNeutralCurrentInputs | None = None,
    overlap_rel_tol: float = DEFAULT_OVERLAP_RTOL,
    min_overlap_modes: int = 16,
    max_overlap_modes: int | None = None,
    quadrature_order: int = DEFAULT_QUADRATURE_ORDER,
) -> dict[str, Any]:
    """Build the 3a neutral-current extras without rewiring constraints."""

    if ew_model != "minimal_rs":
        raise ValueError("Phase 3a supports ew_model='minimal_rs' only")
    if include_fermion_kk_mixing:
        raise ValueError("fermion-KK mixing is deferred beyond Phase 3a")
    if include_loop_dipoles:
        raise ValueError("loop dipoles are deferred beyond Phase 3a")

    spectrum = RSEWSpectrum.build(
        lambda_ir_gev=float(Lambda_IR),
        k_gev=float(k),
        n_gauge_modes=int(n_gauge_modes),
        quadrature_order=int(quadrature_order),
        model_label=str(ew_model),
    )
    couplings = build_rs_ew_couplings(
        quark_fit_result,
        spectrum=spectrum,
        inputs=neutral_current_inputs,
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
    }


__all__ = ["build_rs_ew_extras"]
