import math
import time
from collections.abc import Mapping
from dataclasses import dataclass, fields, is_dataclass

import numpy as np
import pytest

from flavor_catalog_constraints import point_builder
from quarkConstraints.rs_ew_spectrum import RSEWOverlapSplineCache, RSEWSpectrum
from warpConfig.wavefuncs import f_IR
from yukawa.compute_yukawas import compute_all_yukawas


GAUGE_ROOT_EPS_1E_MINUS_15 = 2.450509663813736
EPSILON_RS = 1.0e-15
N_GAUGE_MODES = 64
QUADRATURE_ORDER = 512
MIN_OVERLAP_MODES = 16
MAX_OVERLAP_MODES = 64
OVERLAP_REL_TOL = 1.0e-3


@dataclass(frozen=True)
class _BulkState:
    c_Q: np.ndarray
    c_u: np.ndarray
    c_d: np.ndarray
    F_Q: np.ndarray
    F_d: np.ndarray
    Y_d_bulk_basis: np.ndarray


@dataclass(frozen=True)
class _QuarkFit:
    bulk_state: _BulkState
    U_L_u: np.ndarray
    U_L_d: np.ndarray
    U_R_u: np.ndarray
    U_R_d: np.ndarray
    masses_down: np.ndarray


def _scales_for_mkk(mkk_gev: float = 3000.0) -> tuple[float, float]:
    lambda_ir = float(mkk_gev) / GAUGE_ROOT_EPS_1E_MINUS_15
    return lambda_ir, lambda_ir / EPSILON_RS


def _rot12(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.complex128)


def _rot23(theta: float) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]], dtype=np.complex128)


def _fit(offset: float = 0.0) -> _QuarkFit:
    c_q = np.array([0.641 + offset, 0.563 + offset, 0.431 + offset], dtype=float)
    c_u = np.array([0.623 - offset, 0.347 + offset, 0.319 + offset], dtype=float)
    c_d = np.array([0.667 - offset, 0.571 + offset, 0.353 + offset], dtype=float)
    y_d = np.array(
        [
            [1.0, 0.02, 1.0e-5],
            [0.01, 1.0, 2.0e-5],
            [3.0e-5, 4.0e-5, 1.0],
        ],
        dtype=np.complex128,
    )
    return _QuarkFit(
        bulk_state=_BulkState(
            c_Q=c_q,
            c_u=c_u,
            c_d=c_d,
            F_Q=np.asarray(f_IR(c_q, EPSILON_RS), dtype=float),
            F_d=np.asarray(f_IR(c_d, EPSILON_RS), dtype=float),
            Y_d_bulk_basis=y_d,
        ),
        U_L_u=_rot12(0.11 + offset) @ _rot23(-0.07),
        U_L_d=_rot12(-0.19 + offset) @ _rot23(0.08),
        U_R_u=_rot12(-0.09 + offset) @ _rot23(0.045),
        U_R_d=_rot12(0.05 + offset) @ _rot23(-0.12),
        masses_down=np.array([0.004, 0.09, 4.18], dtype=float),
    )


@pytest.fixture(scope="module")
def spectrum() -> RSEWSpectrum:
    lambda_ir, k = _scales_for_mkk()
    return RSEWSpectrum.build(
        lambda_ir_gev=lambda_ir,
        k_gev=k,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        model_label="minimal_rs",
    )


@pytest.fixture(scope="module")
def overlap_cache(spectrum: RSEWSpectrum) -> RSEWOverlapSplineCache:
    return RSEWOverlapSplineCache.build(
        spectrum,
        c_min=0.3,
        c_max=0.9,
        grid_size=241,
        verify_points=41,
        rel_tol=OVERLAP_REL_TOL,
        min_modes=MIN_OVERLAP_MODES,
        max_modes=MAX_OVERLAP_MODES,
        include_omega=True,
    )


@pytest.fixture(scope="module")
def lepton_yukawa_result():
    lambda_ir, k = _scales_for_mkk()
    return compute_all_yukawas(
        Lambda_IR=lambda_ir,
        c_L=0.587,
        c_E=[0.587, 0.587, 0.587],
        c_N=0.27,
        M_N=1.22e18,
        lightest_nu_mass=0.002,
        ordering="normal",
        majorana_alpha=0.0,
        majorana_beta=0.0,
        k=k,
        v=174.0,
    )


def _point_kwargs(lepton_yukawa_result) -> dict[str, object]:
    lambda_ir, k = _scales_for_mkk()
    return {
        "Lambda_IR": lambda_ir,
        "k": k,
        "n_gauge_modes": N_GAUGE_MODES,
        "quadrature_order": QUADRATURE_ORDER,
        "min_overlap_modes": MIN_OVERLAP_MODES,
        "max_overlap_modes": MAX_OVERLAP_MODES,
        "overlap_rel_tol": OVERLAP_REL_TOL,
        "include_fermion_kk_mixing": True,
        "include_charged_current": True,
        "include_higgs_yukawas": True,
        "lepton_yukawa_result": lepton_yukawa_result,
    }


def test_a_spline_accuracy_against_direct_overlap(
    spectrum: RSEWSpectrum,
    overlap_cache: RSEWOverlapSplineCache,
):
    c_values = np.linspace(0.307, 0.893, 47, dtype=float)
    direct = np.array(
        [
            spectrum.a(
                float(c),
                rel_tol=OVERLAP_REL_TOL,
                min_modes=MIN_OVERLAP_MODES,
                max_modes=MAX_OVERLAP_MODES,
            )
            for c in c_values
        ],
        dtype=float,
    )
    interpolated = np.array([overlap_cache.a_spline(float(c)) for c in c_values], dtype=float)
    max_rel = _max_array_rel_diff(direct, interpolated)
    omega_c_values = np.array([0.321, 0.487, 0.587, 0.731, 0.887], dtype=float)
    omega_direct = np.stack(
        [spectrum.omega(float(c), max_modes=MAX_OVERLAP_MODES) for c in omega_c_values],
        axis=0,
    )
    omega_interpolated = np.stack(
        [overlap_cache.omega(float(c), max_modes=MAX_OVERLAP_MODES) for c in omega_c_values],
        axis=0,
    )
    omega_max_rel = _max_array_rel_diff(omega_direct, omega_interpolated)

    print(f"W6A spline max_rel_err={max_rel:.6e}")
    print(f"W6A omega_spline max_rel_err={omega_max_rel:.6e}")
    assert overlap_cache.max_a_rel_err < 1.0e-3
    assert max_rel < 1.0e-3
    assert omega_max_rel < 1.0e-3
    assert omega_interpolated.shape == (omega_c_values.size, MAX_OVERLAP_MODES)
    assert np.all(np.isfinite(omega_interpolated))


def test_injected_spectrum_spline_path_matches_rebuild_path(
    spectrum: RSEWSpectrum,
    overlap_cache: RSEWOverlapSplineCache,
    lepton_yukawa_result,
):
    kwargs = _point_kwargs(lepton_yukawa_result)
    rebuild = point_builder.build_from_rs_ew_inputs(_fit(), **kwargs)
    injected = point_builder.build_from_rs_ew_inputs(
        _fit(),
        **kwargs,
        spectrum=spectrum,
        rs_ew_cache=overlap_cache,
    )
    callable_injected = point_builder.build_from_rs_ew_inputs(
        _fit(),
        **kwargs,
        spectrum=spectrum,
        a_of_c=overlap_cache.a_spline,
        omega_of_c=overlap_cache.omega_spline,
    )

    max_rel = 0.0
    for key in (
        "rs_ew_spectrum",
        "lepton_mass_basis_couplings",
        "rs_ew_couplings",
        "rs_semileptonic_wilsons",
        "rs_charged_current",
        "rs_higgs_yukawas",
    ):
        max_rel = max(max_rel, _max_numeric_rel_diff(rebuild.extras[key], injected.extras[key]))
        max_rel = max(
            max_rel,
            _max_numeric_rel_diff(rebuild.extras[key], callable_injected.extras[key]),
        )

    print(f"W6A rebuild_vs_injected max_rel_diff={max_rel:.6e}")
    assert injected.extras["rs_ew_spectrum"] is spectrum
    assert callable_injected.extras["rs_ew_spectrum"] is spectrum
    assert max_rel < 1.0e-3


def test_injected_cache_and_callable_metadata_fail_closed(
    spectrum: RSEWSpectrum,
    overlap_cache: RSEWOverlapSplineCache,
    lepton_yukawa_result,
):
    kwargs = _point_kwargs(lepton_yukawa_result)
    with pytest.raises(TypeError, match="rs_ew_cache must expose"):
        point_builder.build_from_rs_ew_inputs(
            _fit(),
            **kwargs,
            spectrum=spectrum,
            rs_ew_cache=object(),
        )
    with pytest.raises(ValueError, match="rel_tol"):
        point_builder.build_from_rs_ew_inputs(
            _fit(),
            **{**kwargs, "overlap_rel_tol": 5.0e-4},
            spectrum=spectrum,
            rs_ew_cache=overlap_cache,
        )
    with pytest.raises(ValueError, match="requires omega_of_c"):
        point_builder.build_from_rs_ew_inputs(
            _fit(),
            **kwargs,
            spectrum=spectrum,
            a_of_c=overlap_cache.a_spline,
        )
    with pytest.raises(ValueError, match="min_modes"):
        point_builder.build_from_rs_ew_inputs(
            _fit(),
            **kwargs,
            spectrum=spectrum,
            a_of_c=lambda c: overlap_cache.a_spline(c),
            omega_of_c=overlap_cache.omega_spline,
        )
    def a_without_domain(c):
        return overlap_cache.a_spline(c)

    a_without_domain.min_modes = MIN_OVERLAP_MODES
    a_without_domain.max_modes = MAX_OVERLAP_MODES
    a_without_domain.rel_tol = OVERLAP_REL_TOL
    with pytest.raises(ValueError, match="c_min/c_max"):
        point_builder.build_from_rs_ew_inputs(
            _fit(),
            **kwargs,
            spectrum=spectrum,
            a_of_c=a_without_domain,
            omega_of_c=overlap_cache.omega_spline,
        )
    lambda_ir, k = _scales_for_mkk()
    duplicate_spectrum = RSEWSpectrum.build(
        lambda_ir_gev=lambda_ir,
        k_gev=k,
        n_gauge_modes=N_GAUGE_MODES,
        quadrature_order=QUADRATURE_ORDER,
        model_label="minimal_rs",
    )
    with pytest.raises(ValueError, match="active injected spectrum"):
        point_builder.build_from_rs_ew_inputs(
            _fit(),
            **kwargs,
            spectrum=duplicate_spectrum,
            a_of_c=overlap_cache.a_spline,
            omega_of_c=overlap_cache.omega_spline,
        )


def test_injected_cache_reduces_amortized_per_point_build_time(
    spectrum: RSEWSpectrum,
    overlap_cache: RSEWOverlapSplineCache,
    lepton_yukawa_result,
):
    kwargs = _point_kwargs(lepton_yukawa_result)
    fits = [_fit(0.0005 * idx) for idx in range(4)]

    start = time.perf_counter()
    for fit in fits:
        point_builder.build_from_rs_ew_inputs(fit, **kwargs)
    rebuild_seconds = time.perf_counter() - start

    start = time.perf_counter()
    for fit in fits:
        point_builder.build_from_rs_ew_inputs(
            fit,
            **kwargs,
            spectrum=spectrum,
            rs_ew_cache=overlap_cache,
        )
    injected_seconds = time.perf_counter() - start
    speedup = (rebuild_seconds / len(fits)) / (injected_seconds / len(fits))

    print(f"W6A amortized_per_point_speedup={speedup:.3f}x")
    assert speedup > 2.0


def _max_numeric_rel_diff(left, right) -> float:
    if isinstance(left, (bool, np.bool_)) or isinstance(right, (bool, np.bool_)):
        return 0.0
    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        return _max_array_rel_diff(np.asarray(left), np.asarray(right))
    if isinstance(left, (float, int, complex, np.number)) and isinstance(
        right, (float, int, complex, np.number)
    ):
        return _max_array_rel_diff(np.asarray([left]), np.asarray([right]))
    if is_dataclass(left) and is_dataclass(right):
        max_rel = 0.0
        for field in fields(left):
            if field.name.startswith("_"):
                continue
            max_rel = max(
                max_rel,
                _max_numeric_rel_diff(getattr(left, field.name), getattr(right, field.name)),
            )
        return max_rel
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        assert set(left) == set(right)
        max_rel = 0.0
        for key in left:
            max_rel = max(max_rel, _max_numeric_rel_diff(left[key], right[key]))
        return max_rel
    if isinstance(left, (tuple, list)) and isinstance(right, (tuple, list)):
        assert len(left) == len(right)
        max_rel = 0.0
        for left_item, right_item in zip(left, right, strict=True):
            max_rel = max(max_rel, _max_numeric_rel_diff(left_item, right_item))
        return max_rel
    return 0.0


def _max_array_rel_diff(left, right) -> float:
    left_arr = np.asarray(left, dtype=np.complex128)
    right_arr = np.asarray(right, dtype=np.complex128)
    assert left_arr.shape == right_arr.shape
    assert np.all(np.isfinite(left_arr))
    assert np.all(np.isfinite(right_arr))
    denom = np.maximum(np.maximum(np.abs(left_arr), np.abs(right_arr)), 1.0e-12)
    return float(np.max(np.abs(left_arr - right_arr) / denom))
