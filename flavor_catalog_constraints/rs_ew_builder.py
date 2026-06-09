"""RS electroweak ``ParameterPoint`` extras builder."""

from __future__ import annotations

import math
from typing import Any, Callable

import numpy as np

from flavor_catalog_constraints.physics_adapters.lepton import (
    lmfv_lepton_parameters_from_yukawa_result,
)
from quarkConstraints.rs_ew_couplings import (
    DEFAULT_CUSTODIAL_B_R_REP,
    DEFAULT_CUSTODIAL_B_R_STRATEGY,
    DEFAULT_CUSTODIAL_PROTECT_SCOPE,
    DEFAULT_CUSTODIAL_Q_L_REP,
    DEFAULT_CUSTODIAL_T_R_REP,
    RSEWNeutralCurrentInputs,
    SUPPORTED_RS_EW_MODELS,
    build_rs_ew_couplings,
    build_rs_lepton_mass_basis_couplings,
)
from quarkConstraints.rs_charged_current import (
    RSChargedCurrentInputs,
    build_rs_charged_current,
)
from quarkConstraints.rs_ew_spectrum import (
    DEFAULT_MAX_TRUNCATION_MODES,
    DEFAULT_N_GAUGE_MODES,
    DEFAULT_OVERLAP_RTOL,
    DEFAULT_QUADRATURE_ORDER,
    RSEWSpectrum,
)
from quarkConstraints.rs_higgs_yukawas import build_rs_higgs_yukawas
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
    qL_rep: str = DEFAULT_CUSTODIAL_Q_L_REP,
    tR_rep: str = DEFAULT_CUSTODIAL_T_R_REP,
    bR_rep: str = DEFAULT_CUSTODIAL_B_R_REP,
    protect_scope: Any = DEFAULT_CUSTODIAL_PROTECT_SCOPE,
    bR_strategy: str = DEFAULT_CUSTODIAL_B_R_STRATEGY,
    kappa_b: float = 0.0,
    custodial_PLR_breaking_residual: bool = False,
    include_top_partner_loops: bool = False,
    include_fermion_kk_mixing: bool = False,
    include_higgs_yukawas: bool = True,
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
    spectrum: RSEWSpectrum | None = None,
    a_of_c: Callable[[float], float] | None = None,
    omega_of_c: Callable[[float], Any] | None = None,
    rs_ew_cache: Any | None = None,
) -> dict[str, Any]:
    """Build RS-EW neutral-current extras without rewiring constraints."""

    if str(ew_model) not in SUPPORTED_RS_EW_MODELS:
        raise ValueError(
            f"unsupported ew_model {ew_model!r}; supported models are "
            f"{SUPPORTED_RS_EW_MODELS}"
        )
    if include_top_partner_loops:
        raise ValueError("include_top_partner_loops=True is deferred to PR2")
    if include_loop_dipoles:
        raise ValueError("loop dipoles are deferred beyond Phase 3a")

    base_spectrum = _resolve_spectrum(
        Lambda_IR=float(Lambda_IR),
        k=float(k),
        n_gauge_modes=int(n_gauge_modes),
        quadrature_order=int(quadrature_order),
        ew_model=str(ew_model),
        spectrum=spectrum,
        rs_ew_cache=rs_ew_cache,
    )
    overlap_provider, resolved_max_overlap_modes = _resolve_overlap_provider(
        base_spectrum,
        min_overlap_modes=int(min_overlap_modes),
        max_overlap_modes=max_overlap_modes,
        overlap_rel_tol=float(overlap_rel_tol),
        include_charged_current=bool(include_charged_current),
        a_of_c=a_of_c,
        omega_of_c=omega_of_c,
        rs_ew_cache=rs_ew_cache,
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
            kk_ew_mass_gev=float(base_spectrum.kk_ew_mass_gev),
        )
    )
    lmfv_leptons = (
        None
        if resolved_lepton_yukawas is None
        else lmfv_lepton_parameters_from_yukawa_result(
            resolved_lepton_yukawas,
            m_kk_gev=float(base_spectrum.kk_ew_mass_gev),
        )
    )
    couplings = build_rs_ew_couplings(
        quark_fit_result,
        spectrum=overlap_provider,
        lepton_mass_basis_couplings=lepton_couplings,
        inputs=neutral_current_inputs,
        include_fermion_kk_mixing=bool(include_fermion_kk_mixing),
        overlap_rel_tol=float(overlap_rel_tol),
        min_overlap_modes=int(min_overlap_modes),
        max_overlap_modes=resolved_max_overlap_modes,
        model_label=str(ew_model),
        qL_rep=qL_rep,
        tR_rep=tR_rep,
        bR_rep=bR_rep,
        protect_scope=protect_scope,
        bR_strategy=bR_strategy,
        kappa_b=float(kappa_b),
        custodial_PLR_breaking_residual=bool(custodial_PLR_breaking_residual),
        include_top_partner_loops=bool(include_top_partner_loops),
    )
    charged_current = None
    higgs_yukawas = None
    if include_higgs_yukawas and lepton_couplings is not None:
        higgs_yukawas = build_rs_higgs_yukawas(
            lepton_couplings,
            spectrum=base_spectrum,
        )
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
            spectrum=overlap_provider,
            lepton_mass_basis_couplings=lepton_couplings,
            inputs=resolved_charged_inputs,
            shared_a_ref=float(couplings.a_ref),
            overlap_rel_tol=float(overlap_rel_tol),
            min_overlap_modes=int(min_overlap_modes),
            max_overlap_modes=resolved_max_overlap_modes,
            model_label=str(ew_model),
        )
    wilsons = build_rs_semileptonic_wilsons(couplings)
    return {
        "kk_ew_mass_gev": float(base_spectrum.kk_ew_mass_gev),
        "rs_ew_spectrum": base_spectrum,
        "rs_ew_couplings": couplings,
        "rs_semileptonic_wilsons": wilsons,
        **(
            {}
            if lepton_couplings is None
            else {"lepton_mass_basis_couplings": lepton_couplings}
        ),
        **(
            {}
            if lmfv_leptons is None
            else {"lepton_lmfv_parameters": lmfv_leptons}
        ),
        **({} if charged_current is None else {"rs_charged_current": charged_current}),
        **({} if higgs_yukawas is None else {"rs_higgs_yukawas": higgs_yukawas}),
    }


class _CallableOverlapProvider:
    """Delegate spectrum data while replacing raw overlap evaluations."""

    def __init__(
        self,
        *,
        spectrum: RSEWSpectrum,
        a_of_c: Callable[[float], float],
        omega_of_c: Callable[[float], Any] | None,
        rel_tol: float,
        min_modes: int,
        max_modes: int,
    ) -> None:
        self.spectrum = spectrum
        self._a_of_c = a_of_c
        self._omega_of_c = omega_of_c
        self.rel_tol = float(rel_tol)
        self.min_modes = int(min_modes)
        self.max_modes = int(max_modes)

    def __getattr__(self, name: str) -> Any:
        return getattr(self.spectrum, name)

    def a(
        self,
        c: float,
        *,
        a_ref: float = 0.0,
        rel_tol: float = DEFAULT_OVERLAP_RTOL,
        min_modes: int = 16,
        max_modes: int = DEFAULT_N_GAUGE_MODES,
    ) -> float:
        _validate_overlap_request(
            rel_tol=float(rel_tol),
            min_modes=int(min_modes),
            max_modes=int(max_modes),
            expected_rel_tol=self.rel_tol,
            expected_min_modes=self.min_modes,
            expected_max_modes=self.max_modes,
            source="a_of_c",
        )
        value = float(self._a_of_c(float(c))) - float(a_ref)
        if not math.isfinite(value):
            raise ValueError("a_of_c returned a non-finite value")
        return value

    def omega(self, c: float, max_modes: int | None = None) -> np.ndarray:
        if self._omega_of_c is None:
            return self.spectrum.omega(float(c), max_modes=max_modes)
        n = self.spectrum.n_gauge_modes if max_modes is None else int(max_modes)
        if n < 1 or n > int(self.spectrum.n_gauge_modes):
            raise ValueError("max_modes must be between 1 and spectrum.n_gauge_modes")
        values = np.array(self._omega_of_c(float(c)), dtype=float, copy=True)
        if values.ndim != 1 or values.size < n:
            raise ValueError("omega_of_c must return a one-dimensional vector covering max_modes")
        if not np.all(np.isfinite(values)):
            raise ValueError("omega_of_c returned non-finite values")
        return values[:n]


def _resolve_spectrum(
    *,
    Lambda_IR: float,
    k: float,
    n_gauge_modes: int,
    quadrature_order: int,
    ew_model: str,
    spectrum: RSEWSpectrum | None,
    rs_ew_cache: Any | None,
) -> RSEWSpectrum:
    cache_spectrum = None if rs_ew_cache is None else getattr(rs_ew_cache, "spectrum", None)
    if rs_ew_cache is not None and cache_spectrum is None:
        raise TypeError("rs_ew_cache must expose the RSEWSpectrum as .spectrum")
    if cache_spectrum is not None and not isinstance(cache_spectrum, RSEWSpectrum):
        raise TypeError("rs_ew_cache.spectrum must be an RSEWSpectrum")
    if spectrum is not None and not isinstance(spectrum, RSEWSpectrum):
        raise TypeError("spectrum must be an RSEWSpectrum")
    if spectrum is not None and cache_spectrum is not None and spectrum is not cache_spectrum:
        raise ValueError("spectrum and rs_ew_cache.spectrum must be the same object")

    active = spectrum if spectrum is not None else cache_spectrum
    if active is None:
        active = RSEWSpectrum.build(
            lambda_ir_gev=float(Lambda_IR),
            k_gev=float(k),
            n_gauge_modes=int(n_gauge_modes),
            quadrature_order=int(quadrature_order),
            model_label=str(ew_model),
        )
    else:
        _validate_injected_spectrum(
            active,
            Lambda_IR=Lambda_IR,
            k=k,
            n_gauge_modes=n_gauge_modes,
            quadrature_order=quadrature_order,
            ew_model=ew_model,
        )
    return active


def _validate_injected_spectrum(
    spectrum: RSEWSpectrum,
    *,
    Lambda_IR: float,
    k: float,
    n_gauge_modes: int,
    quadrature_order: int,
    ew_model: str,
) -> None:
    checks = (
        ("lambda_ir_gev", float(Lambda_IR)),
        ("k_gev", float(k)),
    )
    for attr, expected in checks:
        actual = float(getattr(spectrum, attr))
        if not math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-9):
            raise ValueError(f"injected spectrum {attr} does not match builder inputs")
    if int(spectrum.n_gauge_modes) != int(n_gauge_modes):
        raise ValueError("injected spectrum n_gauge_modes does not match builder inputs")
    if int(spectrum.quadrature_order) != int(quadrature_order):
        raise ValueError("injected spectrum quadrature_order does not match builder inputs")
    if str(spectrum.model_label) != str(ew_model):
        raise ValueError("injected spectrum model_label does not match ew_model")


def _resolve_overlap_provider(
    spectrum: RSEWSpectrum,
    *,
    min_overlap_modes: int,
    max_overlap_modes: int | None,
    overlap_rel_tol: float,
    include_charged_current: bool,
    a_of_c: Callable[[float], float] | None,
    omega_of_c: Callable[[float], Any] | None,
    rs_ew_cache: Any | None,
) -> tuple[Any, int | None]:
    if rs_ew_cache is not None and (a_of_c is not None or omega_of_c is not None):
        raise ValueError("provide either rs_ew_cache or a_of_c/omega_of_c, not both")

    resolved_max = (
        _default_max_overlap_modes(spectrum)
        if max_overlap_modes is None and (rs_ew_cache is not None or a_of_c is not None)
        else max_overlap_modes
    )
    if rs_ew_cache is not None:
        if not hasattr(rs_ew_cache, "a"):
            raise TypeError("rs_ew_cache must provide an a(c, ...) method")
        if include_charged_current and not hasattr(rs_ew_cache, "omega"):
            raise TypeError("include_charged_current=True requires rs_ew_cache.omega(c)")
        if getattr(rs_ew_cache, "spectrum", None) is not spectrum:
            raise ValueError("rs_ew_cache.spectrum must be the active injected spectrum")
        _validate_overlap_metadata(
            rs_ew_cache,
            source="rs_ew_cache",
            min_modes=int(min_overlap_modes),
            max_modes=int(resolved_max),
            rel_tol=float(overlap_rel_tol),
            require_domain=True,
        )
        if include_charged_current:
            omega_grid = getattr(rs_ew_cache, "omega_grid", None)
            if omega_grid is None:
                raise ValueError(
                    "include_charged_current=True requires rs_ew_cache with Omega_n(c) precomputed"
                )
        cache_max = getattr(rs_ew_cache, "max_modes", None)
        if cache_max is not None:
            cache_max_int = int(cache_max)
            if max_overlap_modes is None:
                resolved_max = cache_max_int
            elif cache_max_int != int(max_overlap_modes):
                raise ValueError("rs_ew_cache.max_modes must match max_overlap_modes")
        return rs_ew_cache, resolved_max

    if a_of_c is not None:
        if include_charged_current and omega_of_c is None:
            raise ValueError(
                "include_charged_current=True with a_of_c requires omega_of_c"
            )
        _validate_overlap_metadata(
            a_of_c,
            source="a_of_c",
            min_modes=int(min_overlap_modes),
            max_modes=int(resolved_max),
            rel_tol=float(overlap_rel_tol),
            require_domain=True,
        )
        a_source = _metadata_source(a_of_c)
        _validate_callable_spectrum(a_source, spectrum=spectrum, source="a_of_c")
        if omega_of_c is not None:
            _validate_callable_max_modes(
                omega_of_c,
                source="omega_of_c",
                max_modes=int(resolved_max),
                require_domain=True,
            )
            omega_source = _metadata_source(omega_of_c)
            if omega_source is not a_source:
                raise ValueError("a_of_c and omega_of_c must come from the same cache object")
            _validate_callable_spectrum(
                omega_source,
                spectrum=spectrum,
                source="omega_of_c",
            )
        return (
            _CallableOverlapProvider(
                spectrum=spectrum,
                a_of_c=a_of_c,
                omega_of_c=omega_of_c,
                rel_tol=float(overlap_rel_tol),
                min_modes=int(min_overlap_modes),
                max_modes=int(resolved_max),
            ),
            resolved_max,
        )
    if omega_of_c is not None:
        raise ValueError("omega_of_c requires a_of_c so the injected overlap path is complete")
    return spectrum, resolved_max


def _default_max_overlap_modes(spectrum: RSEWSpectrum) -> int:
    max_modes = min(int(spectrum.n_gauge_modes), DEFAULT_MAX_TRUNCATION_MODES)
    max_modes = 1 << (max_modes.bit_length() - 1)
    if max_modes < 32:
        raise ValueError("max_overlap_modes must be at least 32")
    return max_modes


def _metadata_source(candidate: Any) -> Any:
    return getattr(candidate, "__self__", candidate)


def _validate_overlap_metadata(
    candidate: Any,
    *,
    source: str,
    min_modes: int,
    max_modes: int,
    rel_tol: float,
    require_domain: bool,
) -> None:
    metadata_source = _metadata_source(candidate)
    for attr, expected in (
        ("min_modes", int(min_modes)),
        ("max_modes", int(max_modes)),
    ):
        actual = getattr(metadata_source, attr, None)
        if actual is None or int(actual) != expected:
            raise ValueError(f"{source}.{attr} must match the requested overlap setup")
    actual_rel_tol = getattr(metadata_source, "rel_tol", None)
    if actual_rel_tol is None or not math.isclose(
        float(actual_rel_tol),
        float(rel_tol),
        rel_tol=0.0,
        abs_tol=0.0,
    ):
        raise ValueError(f"{source}.rel_tol must match the requested overlap setup")
    if require_domain:
        _validate_callable_domain(metadata_source, source=source)


def _validate_callable_max_modes(
    candidate: Any,
    *,
    source: str,
    max_modes: int,
    require_domain: bool,
) -> None:
    metadata_source = _metadata_source(candidate)
    actual = getattr(metadata_source, "max_modes", None)
    if actual is None or int(actual) != int(max_modes):
        raise ValueError(f"{source}.max_modes must match the requested overlap setup")
    if require_domain:
        _validate_callable_domain(metadata_source, source=source)


def _validate_callable_domain(candidate: Any, *, source: str) -> None:
    c_min = getattr(candidate, "c_min", None)
    c_max = getattr(candidate, "c_max", None)
    if c_min is None or c_max is None:
        raise ValueError(f"{source} must expose c_min/c_max metadata")
    if float(c_min) > 0.3 or float(c_max) < 0.9:
        raise ValueError(f"{source} must cover the W6 c-domain [0.3, 0.9]")


def _validate_callable_spectrum(
    candidate: Any,
    *,
    spectrum: RSEWSpectrum,
    source: str,
) -> None:
    actual = getattr(candidate, "spectrum", None)
    if actual is None or actual is not spectrum:
        raise ValueError(f"{source}.spectrum must be the active injected spectrum")


def _validate_overlap_request(
    *,
    rel_tol: float,
    min_modes: int,
    max_modes: int,
    expected_rel_tol: float,
    expected_min_modes: int,
    expected_max_modes: int,
    source: str,
) -> None:
    if int(min_modes) != int(expected_min_modes) or int(max_modes) != int(expected_max_modes):
        raise ValueError(f"{source} request does not match the injected overlap modes")
    if not math.isclose(float(rel_tol), float(expected_rel_tol), rel_tol=0.0, abs_tol=0.0):
        raise ValueError(f"{source} request does not match the injected overlap rel_tol")


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
