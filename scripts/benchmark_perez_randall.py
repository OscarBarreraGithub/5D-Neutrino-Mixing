#!/usr/bin/env python
"""Benchmark check for Perez–Randall point.

Run:
  python scripts/benchmark_perez_randall.py
"""

import sys
import numpy as np

try:
    from yukawa import compute_all_yukawas
    from flavorConstraints import check_mu_to_e_gamma
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "Cannot import packages. Install the repo first: `pip install -e .`."
    ) from exc


def main() -> int:
    result = compute_all_yukawas(
        Lambda_IR=3000,
        c_L=0.58,
        c_E=[0.75, 0.60, 0.50],
        c_N=0.27,
        M_N=1.22e18,
        lightest_nu_mass=0.002,
        ordering='normal',
    )

    expected = {
        "f_L": 0.01598,
        "f_N": 0.4796,
        "f_N_UV": 1.232e-4,
        "Y_E_bar": np.array([2.94, 4.37, 5.42]),
        "Y_N_bar": np.array([0.204, 0.431, 1.024]),
    }

    tol = 0.02  # 2% relative tolerance

    checks = {
        "f_L": np.isclose(result.f_L, expected["f_L"], rtol=tol),
        "f_N": np.isclose(result.f_N, expected["f_N"], rtol=tol),
        "f_N_UV": np.isclose(result.f_N_UV, expected["f_N_UV"], rtol=tol),
        "Y_E_bar": np.allclose(result.Y_E_bar, expected["Y_E_bar"], rtol=tol),
        "Y_N_bar": np.allclose(result.Y_N_bar, expected["Y_N_bar"], rtol=tol),
    }

    print("Perez–Randall benchmark:")
    print(f"  f_L    = {result.f_L:.6f}")
    print(f"  f_N    = {result.f_N:.6f}")
    print(f"  f_N_UV = {result.f_N_UV:.6e}")
    print(f"  Y_E_bar = {result.Y_E_bar}")
    print(f"  Y_N_bar = {result.Y_N_bar}")
    print()

    failed = [k for k, ok in checks.items() if not ok]
    if failed:
        print("FAILED checks:", ", ".join(failed))
        return 1

    print("All checks passed within 2% relative tolerance.")
    print()

    # --- μ→eγ constraint check (informational) ---
    lfv = check_mu_to_e_gamma(result)
    print("μ→eγ constraint (at M_KK = 3 TeV):")
    print(f"  |(Ȳ_N Ȳ_N†)₁₂| = {lfv['lhs']:.6f}")
    print(f"  Bound            = {lfv['rhs']:.6f}")
    print(f"  Ratio            = {lfv['ratio']:.4f}")
    print(f"  Passes           = {lfv['passes']}")
    if not lfv['passes']:
        # The benchmark point requires M_KK > 3 TeV to satisfy μ→eγ.
        min_mkk = 3000.0 * np.sqrt(lfv['ratio'])
        print(f"  → Minimum M_KK for this point ≈ {min_mkk:.0f} GeV")

    return 0


if __name__ == "__main__":
    sys.exit(main())
