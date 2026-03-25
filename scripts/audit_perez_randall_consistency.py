#!/usr/bin/env python
"""Audit the internal consistency of the displayed Perez–Randall benchmark.

This script checks whether the numbers shown in arXiv:0805.4652 Eq. (10) can
reproduce the example neutrino spectrum in Eq. (7) when combined with Eq. (6)
and the quoted overlap data in Eq. (11) / Table I.
"""

from __future__ import annotations

import numpy as np

from warpConfig.baseParams import MBARPL, get_warp_params
from warpConfig.wavefuncs import f_IR, f_UV


def paper_like_geometry():
    """Return the paper-like geometry inferred from Eq. (11)/Table I."""
    # The paper text uses epsilon = exp(-k pi rc) with k pi rc = log(Mbar_Pl / TeV),
    # so Lambda_IR = k * epsilon ~ 1 TeV when k = Mbar_Pl.
    k = MBARPL
    Lambda_IR = 1000.0
    geometry = get_warp_params(k=k, Lambda_IR=Lambda_IR)
    epsilon = geometry["epsilon"]
    return {
        "k": k,
        "Lambda_IR": Lambda_IR,
        "epsilon": epsilon,
        "f_L": float(f_IR(0.58, epsilon)),
        "f_N": float(f_IR(0.27, epsilon)),
        "f_N_UV": float(f_UV(0.27, epsilon)),
        "M_N": k / 10.0,
    }


def masses_from_kyn(
    kY_N: np.ndarray,
    f_L: float,
    f_N: float,
    f_N_UV: float,
    M_N: float,
) -> np.ndarray:
    """Compute light neutrino masses implied by Eq. (6), returned in eV."""
    return 2.0 * (174.0**2) * (f_L**2) * (f_N**2) * (kY_N**2) / ((f_N_UV**2) * M_N) * 1.0e9


def required_kyn(
    masses_eV: np.ndarray,
    f_L: float,
    f_N: float,
    f_N_UV: float,
    M_N: float,
) -> np.ndarray:
    """Invert Eq. (6) in the universal limit to obtain kY_N."""
    return np.sqrt(
        masses_eV * 1.0e-9 * (f_N_UV**2) * M_N / (2.0 * (174.0**2) * (f_L**2) * (f_N**2))
    )


def main() -> int:
    state = paper_like_geometry()

    eq7_masses = np.array([0.002, 0.009, 0.05])
    eq10_kyn = np.array([0.02, 0.03, 0.07])

    implied_masses = masses_from_kyn(
        eq10_kyn,
        f_L=state["f_L"],
        f_N=state["f_N"],
        f_N_UV=state["f_N_UV"],
        M_N=state["M_N"],
    )
    needed_kyn = required_kyn(
        eq7_masses,
        f_L=state["f_L"],
        f_N=state["f_N"],
        f_N_UV=state["f_N_UV"],
        M_N=state["M_N"],
    )

    print("Perez–Randall displayed-benchmark audit")
    print(f"  k = Mbar_Pl = {state['k']:.6e} GeV")
    print(f"  Lambda_IR = {state['Lambda_IR']:.0f} GeV")
    print(f"  f_L = {state['f_L']:.6f}  (paper: 0.016)")
    print(f"  f_N = {state['f_N']:.6f}  (paper: 0.48)")
    print(f"  f_N_UV = {state['f_N_UV']:.6e}  (paper: 1.6e-4)")
    print(f"  M_N = {state['M_N']:.6e} GeV  (= Mbar_Pl/10)")
    print()
    print(f"  Eq. (7) masses [eV]:        {eq7_masses}")
    print(f"  Eq. (10) kY_N:             {eq10_kyn}")
    print(f"  Eq. (10) -> Eq. (6) masses: {implied_masses}")
    print(f"  Eq. (7) -> Eq. (6) kY_N:    {needed_kyn}")
    print()
    print("Conclusion:")
    print("  Under paper-like geometry choices, the displayed Eq. (10) neutrino Yukawas")
    print("  do not reproduce the example spectrum in Eq. (7).")
    print("  The discrepancy is O(3-4) in kY_N, not an O(1) convention drift.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
