"""RS electroweak gauge-KK spectrum and overlap kernel.

This module implements Phase 2 of the RS electroweak sector design.  It is
standalone: it builds the exact gauge-NN KK tower from ``solvers.bessel``,
normalizes the relative gauge profiles ``chi_n(t)``, and computes the
zero-mode overlap form factor

    a(c) = sum_n (M_KK**2 / m_n**2) chi_n(1) Omega_n(c)

with the doubled truncation criterion pinned in
``derivations/rs_ew_gauge_kk_coupling.tex``.  Coupling-level use of
``a(c) - a_ref`` is intentionally left to later phases.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from functools import lru_cache
from typing import Any, Mapping

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.special import j0, j1, y0, y1

from solvers.bessel import solve_kk
from warpConfig.baseParams import DEFAULT_LAMBDA_IR, MPL, get_warp_params


DEFAULT_N_GAUGE_MODES = 512
DEFAULT_OVERLAP_RTOL = 1.0e-3
DEFAULT_MIN_TRUNCATION_MODES = 16
DEFAULT_MAX_TRUNCATION_MODES = 512
DEFAULT_QUADRATURE_ORDER = 4096
_C_HALF_ATOL = 1.0e-12
_DENOM_FLOOR = 1.0e-12
_ROOT_MIN_X = 1.0e-10
_ROOT_SCAN_STEP = math.pi / 16.0
_ROOT_SEPARATION_ATOL = 1.0e-8


class RSEWOverlapConvergenceError(RuntimeError):
    """Raised when the approved KK truncation test does not converge."""


@dataclass(frozen=True)
class RSEWOverlapResult:
    """Converged value and diagnostics for the numerical overlap form factor."""

    value: float
    raw_value: float
    a_ref: float
    modes: int
    previous_modes: int
    relative_delta: float
    partial_values: Mapping[int, float]


def g0_squared(c: float, t: float | np.ndarray, epsilon: float) -> float | np.ndarray:
    """Return the normalized fermion zero-mode profile squared, g0(c,t)^2."""

    c_float = float(c)
    eps = float(epsilon)
    t_arr = np.asarray(t, dtype=float)
    warp_log = -math.log(eps)
    power = 1.0 - 2.0 * c_float

    if abs(power) < _C_HALF_ATOL:
        out = np.full_like(t_arr, 1.0 / warp_log, dtype=float)
    else:
        denom = -math.expm1(power * math.log(eps))
        prefactor = power / denom
        out = prefactor * np.power(t_arr, power)

    if np.isscalar(t):
        return float(out)
    return out


def g0(c: float, t: float | np.ndarray, epsilon: float) -> float | np.ndarray:
    """Return the normalized fermion zero-mode profile g0(c,t)."""

    out = np.sqrt(np.maximum(g0_squared(c, t, epsilon), 0.0))
    if np.isscalar(t):
        return float(out)
    return out


def w0(c: float, t: float | np.ndarray, epsilon: float) -> float | np.ndarray:
    """Return the normalized zero-mode quadrature weight w0(c,t)=g0(c,t)^2/t."""

    out = np.asarray(g0_squared(c, t, epsilon), dtype=float) / np.asarray(t, dtype=float)
    if np.isscalar(t):
        return float(out)
    return out


def _validate_power_of_two(value: int, *, name: str) -> None:
    if value <= 0 or value & (value - 1):
        raise ValueError(f"{name} must be a positive power of two")


@lru_cache(maxsize=16)
def _log_quadrature(epsilon: float, order: int) -> tuple[np.ndarray, np.ndarray, float]:
    if order < 64:
        raise ValueError("quadrature order must be at least 64")
    eps = float(epsilon)
    if not (0.0 < eps < 1.0):
        raise ValueError("epsilon must satisfy 0 < epsilon < 1")

    nodes, weights = leggauss(int(order))
    y = 0.5 * (nodes + 1.0)
    w = 0.5 * weights
    log_eps = math.log(eps)
    t = np.exp(log_eps * (1.0 - y))
    return t, w, -log_eps


def _gauge_profile_norms(
    roots_x: np.ndarray,
    bessel_b: np.ndarray,
    *,
    epsilon: float,
    quadrature_order: int,
) -> np.ndarray:
    t, weights, warp_log = _log_quadrature(float(epsilon), int(quadrature_order))
    u = np.asarray(roots_x, dtype=float)[:, None] * t[None, :]
    raw = t[None, :] * (j1(u) + np.asarray(bessel_b, dtype=float)[:, None] * y1(u))
    norm_sq = warp_log * np.sum(weights[None, :] * raw * raw, axis=1)
    if not np.all(np.isfinite(norm_sq)) or np.any(norm_sq <= 0.0):
        raise ValueError("non-finite or non-positive gauge profile normalization")
    return 1.0 / np.sqrt(norm_sq)


def _gauge_nn_cross_product(x: float, epsilon: float) -> float:
    return float(j0(x) * y0(epsilon * x) - j0(epsilon * x) * y0(x))


def _solve_ordered_gauge_nn_roots(
    *,
    epsilon: float,
    n_roots: int,
    x_max: float,
    tol: float = 1.0e-12,
) -> np.ndarray:
    """Return the first gauge-NN roots in strict ascending order."""

    eps = float(epsilon)
    n = int(n_roots)
    if n < 1:
        raise ValueError("n_roots must be positive")
    if not (0.0 < eps < 1.0):
        raise ValueError("epsilon must satisfy 0 < epsilon < 1")

    roots: list[float] = []
    x_limit = max(float(x_max), (n + 2.0) * math.pi)
    x_prev = _ROOT_MIN_X
    f_prev = _gauge_nn_cross_product(x_prev, eps)

    while len(roots) < n and x_prev < x_limit:
        x_next = min(x_prev + _ROOT_SCAN_STEP, x_limit)
        f_next = _gauge_nn_cross_product(x_next, eps)

        if np.sign(f_prev) == 0.0:
            root = x_prev
        elif np.sign(f_prev) != np.sign(f_next):
            root = brentq(
                _gauge_nn_cross_product,
                x_prev,
                x_next,
                args=(eps,),
                xtol=tol,
                rtol=tol,
                maxiter=200,
            )
        else:
            root = None

        if root is not None:
            if roots and root <= roots[-1] + _ROOT_SEPARATION_ATOL:
                raise ValueError("gauge-NN root scan produced a duplicate root")
            roots.append(float(root))

        x_prev, f_prev = x_next, f_next

    if len(roots) != n:
        raise ValueError(f"requested {n} gauge roots, found {len(roots)} up to x={x_limit}")

    return np.asarray(roots, dtype=float)


def _exact_gauge_nn_tower(
    geometry: Mapping[str, float],
    *,
    n_roots: int,
    x_max: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    roots = _solve_ordered_gauge_nn_roots(
        epsilon=float(geometry["epsilon"]),
        n_roots=int(n_roots),
        x_max=float(x_max),
    )
    bvals = -j0(roots) / y0(roots)
    masses = roots * float(geometry["Lambda_IR"])
    return masses, roots, np.asarray(bvals, dtype=float)


@dataclass(frozen=True)
class RSEWSpectrum:
    """Exact gauge-NN KK tower and normalized relative profiles."""

    model_label: str
    k_gev: float
    lambda_ir_gev: float
    epsilon: float
    warp_log: float
    n_gauge_modes: int
    gauge_roots_x: np.ndarray
    gauge_masses_gev: np.ndarray
    bessel_b: np.ndarray
    gauge_profile_norms: np.ndarray
    exact_bessel: bool
    profile_normalization: str = "int_epsilon_1_dt_over_t_hm_hn_delta"
    quadrature_order: int = DEFAULT_QUADRATURE_ORDER
    metadata: Mapping[str, Any] = field(default_factory=dict)
    _omega_cache: dict[float, np.ndarray] = field(default_factory=dict, init=False, repr=False)
    _a_cache: dict[tuple[float, float, int, int], RSEWOverlapResult] = field(
        default_factory=dict, init=False, repr=False
    )

    def __post_init__(self) -> None:
        roots = np.asarray(self.gauge_roots_x, dtype=float)
        masses = np.asarray(self.gauge_masses_gev, dtype=float)
        bvals = np.asarray(self.bessel_b, dtype=float)
        norms = np.asarray(self.gauge_profile_norms, dtype=float)
        n_modes = int(self.n_gauge_modes)

        for name, values in (
            ("gauge_roots_x", roots),
            ("gauge_masses_gev", masses),
            ("bessel_b", bvals),
            ("gauge_profile_norms", norms),
        ):
            if values.shape != (n_modes,):
                raise ValueError(f"{name} must have shape ({n_modes},)")
            if not np.all(np.isfinite(values)):
                raise ValueError(f"{name} contains non-finite values")

        if n_modes < 1:
            raise ValueError("n_gauge_modes must be positive")
        if not (0.0 < float(self.epsilon) < 1.0):
            raise ValueError("epsilon must satisfy 0 < epsilon < 1")
        if not np.all(roots > 0.0):
            raise ValueError("gauge roots must be positive")
        if roots.size > 1 and not np.all(np.diff(roots) > _ROOT_SEPARATION_ATOL):
            raise ValueError("gauge roots must be strictly increasing and unique")
        if not np.all(masses > 0.0):
            raise ValueError("gauge masses must be positive")
        if not np.all(norms > 0.0):
            raise ValueError("gauge profile normalizations must be positive")

        object.__setattr__(self, "gauge_roots_x", roots)
        object.__setattr__(self, "gauge_masses_gev", masses)
        object.__setattr__(self, "bessel_b", bvals)
        object.__setattr__(self, "gauge_profile_norms", norms)
        object.__setattr__(self, "k_gev", float(self.k_gev))
        object.__setattr__(self, "lambda_ir_gev", float(self.lambda_ir_gev))
        object.__setattr__(self, "epsilon", float(self.epsilon))
        object.__setattr__(self, "warp_log", float(self.warp_log))
        object.__setattr__(self, "n_gauge_modes", n_modes)
        object.__setattr__(self, "quadrature_order", int(self.quadrature_order))

    @classmethod
    def build(
        cls,
        *,
        lambda_ir_gev: float = DEFAULT_LAMBDA_IR,
        k_gev: float = MPL,
        n_gauge_modes: int = DEFAULT_N_GAUGE_MODES,
        exact_bessel: bool = True,
        quadrature_order: int = DEFAULT_QUADRATURE_ORDER,
        model_label: str = "minimal_rs",
    ) -> "RSEWSpectrum":
        """Build a spectrum from the repo geometry and exact gauge-NN solver."""

        geometry = get_warp_params(k=float(k_gev), Lambda_IR=float(lambda_ir_gev))
        x_max = max(200.0, (int(n_gauge_modes) + 2.0) * math.pi)
        if exact_bessel:
            masses, roots, bvals = _exact_gauge_nn_tower(
                geometry,
                n_roots=int(n_gauge_modes),
                x_max=x_max,
            )
            solver_labels = {
                "species": "gauge",
                "bc": "NN",
                "exact": True,
                "root_ordering": "ordered_cross_product_scan",
            }
        else:
            masses, extras = solve_kk(
                "gauge",
                "NN",
                geometry,
                n_roots=int(n_gauge_modes),
                exact=False,
                x_max=x_max,
            )
            roots = np.asarray(extras["x"], dtype=float)
            bvals = np.asarray(extras["b"], dtype=float)
            solver_labels = extras.get("labels", {})
        if roots.shape != (int(n_gauge_modes),):
            raise ValueError(
                f"requested {n_gauge_modes} gauge roots, got {roots.shape[0]}"
            )
        norms = _gauge_profile_norms(
            roots,
            bvals,
            epsilon=float(geometry["epsilon"]),
            quadrature_order=int(quadrature_order),
        )
        return cls(
            model_label=model_label,
            k_gev=float(geometry["k"]),
            lambda_ir_gev=float(geometry["Lambda_IR"]),
            epsilon=float(geometry["epsilon"]),
            warp_log=float(geometry["warp_log"]),
            n_gauge_modes=int(n_gauge_modes),
            gauge_roots_x=roots,
            gauge_masses_gev=np.asarray(masses, dtype=float),
            bessel_b=bvals,
            gauge_profile_norms=norms,
            exact_bessel=bool(exact_bessel),
            quadrature_order=int(quadrature_order),
            metadata={"solver_labels": solver_labels},
        )

    @property
    def kk_ew_mass_gev(self) -> float:
        """First physical electroweak gauge KK mass, x_1 * Lambda_IR."""

        return float(self.gauge_masses_gev[0])

    @property
    def kk_ew_mass(self) -> float:
        """Alias for ``kk_ew_mass_gev``."""

        return self.kk_ew_mass_gev

    def _profile_values(
        self, mode_slice: slice | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        t, _, warp_log = _log_quadrature(self.epsilon, self.quadrature_order)
        roots = self.gauge_roots_x[mode_slice]
        bvals = self.bessel_b[mode_slice]
        norms = self.gauge_profile_norms[mode_slice]
        u = roots[:, None] * t[None, :]
        raw = t[None, :] * (j1(u) + bvals[:, None] * y1(u))
        chi = math.sqrt(warp_log) * norms[:, None] * raw
        return t, chi

    def chi(self, mode_index: int, t: float | np.ndarray) -> float | np.ndarray:
        """Return chi for a zero-based massive-mode index."""

        idx = int(mode_index)
        if idx < 0 or idx >= self.n_gauge_modes:
            raise IndexError("mode_index out of range")
        t_arr = np.asarray(t, dtype=float)
        u = self.gauge_roots_x[idx] * t_arr
        raw = t_arr * (j1(u) + self.bessel_b[idx] * y1(u))
        out = math.sqrt(self.warp_log) * self.gauge_profile_norms[idx] * raw
        if np.isscalar(t):
            return float(out)
        return out

    def chi_n(self, mode_number: int, t: float | np.ndarray) -> float | np.ndarray:
        """Return chi_n for a one-based massive-mode number n >= 1."""

        return self.chi(int(mode_number) - 1, t)

    def chi_ir(self, max_modes: int | None = None) -> np.ndarray:
        """Return chi_n(1) for the first ``max_modes`` massive modes."""

        n = self.n_gauge_modes if max_modes is None else int(max_modes)
        if n < 1 or n > self.n_gauge_modes:
            raise ValueError("max_modes must be between 1 and n_gauge_modes")
        return np.asarray([self.chi(i, 1.0) for i in range(n)], dtype=float)

    def omega(self, c: float, max_modes: int | None = None) -> np.ndarray:
        """Compute Omega_n(c) for the first ``max_modes`` massive modes."""

        n = self.n_gauge_modes if max_modes is None else int(max_modes)
        if n < 1 or n > self.n_gauge_modes:
            raise ValueError("max_modes must be between 1 and n_gauge_modes")

        c_key = float(c)
        cached = self._omega_cache.get(c_key)
        if cached is None:
            t, weights, warp_log = _log_quadrature(self.epsilon, self.quadrature_order)
            _, chi = self._profile_values(slice(None))
            density_y = warp_log * np.asarray(g0_squared(c_key, t, self.epsilon), dtype=float)
            omega = np.sum(weights[None, :] * density_y[None, :] * chi, axis=1)
            if not np.all(np.isfinite(omega)):
                raise ValueError("Omega_n contains non-finite values")
            self._omega_cache[c_key] = omega
            cached = omega
        return np.array(cached[:n], dtype=float, copy=True)

    def omega_n(self, mode_number: int, c: float) -> float:
        """Return Omega_n(c) for a one-based massive-mode number n >= 1."""

        n = int(mode_number)
        if n < 1 or n > self.n_gauge_modes:
            raise IndexError("mode_number out of range")
        return float(self.omega(c, max_modes=n)[n - 1])

    def a_terms(self, c: float, max_modes: int | None = None) -> np.ndarray:
        """Return the individual summands entering the numerical a(c)."""

        n = self.n_gauge_modes if max_modes is None else int(max_modes)
        omega = self.omega(c, max_modes=n)
        chi_ir = self.chi_ir(max_modes=n)
        weights = (self.kk_ew_mass_gev * self.kk_ew_mass_gev) / (
            self.gauge_masses_gev[:n] * self.gauge_masses_gev[:n]
        )
        terms = weights * chi_ir * omega
        if not np.all(np.isfinite(terms)):
            raise ValueError("a(c) terms contain non-finite values")
        return terms

    def a_with_diagnostics(
        self,
        c: float,
        *,
        a_ref: float = 0.0,
        rel_tol: float = DEFAULT_OVERLAP_RTOL,
        min_modes: int = DEFAULT_MIN_TRUNCATION_MODES,
        max_modes: int = DEFAULT_MAX_TRUNCATION_MODES,
    ) -> RSEWOverlapResult:
        """Return a(c)-a_ref with the approved doubled-N convergence test."""

        min_n = int(min_modes)
        max_n = int(max_modes)
        _validate_power_of_two(min_n, name="min_modes")
        _validate_power_of_two(max_n, name="max_modes")
        if max_n <= min_n:
            raise ValueError("max_modes must be greater than min_modes")
        if max_n > self.n_gauge_modes:
            raise ValueError("spectrum does not contain enough modes for max_modes")
        if rel_tol < 0.0:
            raise ValueError("rel_tol must be non-negative")

        cache_key = (float(c), float(rel_tol), min_n, max_n)
        raw_result = self._a_cache.get(cache_key)
        if raw_result is None:
            terms = self.a_terms(float(c), max_modes=max_n)
            partial_sums = np.cumsum(terms)
            partial_values: dict[int, float] = {}
            previous = min_n
            while previous * 2 <= max_n:
                doubled = previous * 2
                a_previous = float(partial_sums[previous - 1])
                a_doubled = float(partial_sums[doubled - 1])
                partial_values[previous] = a_previous
                partial_values[doubled] = a_doubled
                rel_delta = abs(a_doubled - a_previous) / max(
                    abs(a_doubled), _DENOM_FLOOR
                )
                if rel_delta < rel_tol:
                    raw_result = RSEWOverlapResult(
                        value=a_doubled,
                        raw_value=a_doubled,
                        a_ref=0.0,
                        modes=doubled,
                        previous_modes=previous,
                        relative_delta=float(rel_delta),
                        partial_values=partial_values,
                    )
                    self._a_cache[cache_key] = raw_result
                    break
                previous = doubled

            if raw_result is None:
                last = float(partial_sums[max_n - 1])
                raise RSEWOverlapConvergenceError(
                    "RS-EW overlap a(c) did not converge by "
                    f"N={max_n}; last a_N={last:.16e}"
                )

        if a_ref == 0.0:
            return raw_result
        return replace(
            raw_result,
            value=float(raw_result.raw_value - float(a_ref)),
            a_ref=float(a_ref),
        )

    def a(
        self,
        c: float,
        *,
        a_ref: float = 0.0,
        rel_tol: float = DEFAULT_OVERLAP_RTOL,
        min_modes: int = DEFAULT_MIN_TRUNCATION_MODES,
        max_modes: int = DEFAULT_MAX_TRUNCATION_MODES,
    ) -> float:
        """Return the converged numerical overlap a(c)-a_ref."""

        return float(
            self.a_with_diagnostics(
                c,
                a_ref=a_ref,
                rel_tol=rel_tol,
                min_modes=min_modes,
                max_modes=max_modes,
            ).value
        )


@lru_cache(maxsize=8)
def _cached_spectrum(
    lambda_ir_gev: float,
    k_gev: float,
    n_gauge_modes: int,
    exact_bessel: bool,
    quadrature_order: int,
) -> RSEWSpectrum:
    return RSEWSpectrum.build(
        lambda_ir_gev=float(lambda_ir_gev),
        k_gev=float(k_gev),
        n_gauge_modes=int(n_gauge_modes),
        exact_bessel=bool(exact_bessel),
        quadrature_order=int(quadrature_order),
    )


def kk_ew_mass(
    *,
    lambda_ir_gev: float = DEFAULT_LAMBDA_IR,
    k_gev: float = MPL,
    exact_bessel: bool = True,
) -> float:
    """Return x_1 * Lambda_IR from the exact gauge-NN Bessel equation."""

    spectrum = _cached_spectrum(
        float(lambda_ir_gev),
        float(k_gev),
        1,
        bool(exact_bessel),
        DEFAULT_QUADRATURE_ORDER,
    )
    return spectrum.kk_ew_mass_gev


def a(
    c: float,
    *,
    spectrum: RSEWSpectrum | None = None,
    lambda_ir_gev: float = DEFAULT_LAMBDA_IR,
    k_gev: float = MPL,
    a_ref: float = 0.0,
    rel_tol: float = DEFAULT_OVERLAP_RTOL,
    min_modes: int = DEFAULT_MIN_TRUNCATION_MODES,
    max_modes: int = DEFAULT_MAX_TRUNCATION_MODES,
    quadrature_order: int = DEFAULT_QUADRATURE_ORDER,
) -> float:
    """Return the converged numerical overlap a(c)-a_ref."""

    active_spectrum = spectrum
    if active_spectrum is None:
        active_spectrum = _cached_spectrum(
            float(lambda_ir_gev),
            float(k_gev),
            int(max_modes),
            True,
            int(quadrature_order),
        )
    return active_spectrum.a(
        c,
        a_ref=a_ref,
        rel_tol=rel_tol,
        min_modes=min_modes,
        max_modes=max_modes,
    )


__all__ = [
    "DEFAULT_MAX_TRUNCATION_MODES",
    "DEFAULT_MIN_TRUNCATION_MODES",
    "DEFAULT_N_GAUGE_MODES",
    "DEFAULT_OVERLAP_RTOL",
    "RSEWOverlapConvergenceError",
    "RSEWOverlapResult",
    "RSEWSpectrum",
    "a",
    "g0",
    "g0_squared",
    "kk_ew_mass",
    "w0",
]
