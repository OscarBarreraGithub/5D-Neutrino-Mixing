"""SU(2)_R vs P_LR decomposition for the custodial report (Section 2).

Custodial = SU(2)_L x SU(2)_R x U(1)_X x P_LR bundles two protections:
  - SU(2)_R protects the oblique T parameter (CGHNP Eq. 153),
  - the discrete P_LR protects the Z b_L vertex (ACDP).
This script runs the three EW models and prints the decomposition table used in
reports/collaborator_2026-06/custodial_2x2_comparison.tex:

  EW model            oblique S,T,U floor      Z b_L shift delta_g_L^b
  minimal_rs          ~16 TeV (T problem)      unprotected
  custodial_rs_su2r   ~5.8 TeV (T fixed, S)    = minimal (still unprotected)
  custodial_rs_plr    ~5.8 TeV (unchanged)     0 (protected)

The oblique floor is the largest M_KK excluded at 95% against the PDG-2025 U-fixed
ellipse (EW001 anchor, s_coefficient = 30.0).  delta_g_L^b is z_delta_g_L_d[2,2] of
the RS-EW couplings at M_KK = 3 TeV, built from a representative quark fit.
"""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tests"))

from quarkConstraints.oblique_stu import (  # noqa: E402
    evaluate_rs_oblique_proxy,
    ObliqueSTFit,
    CHI2_2DOF_95,
)

MODELS = ("minimal_rs", "custodial_rs_su2r", "custodial_rs_plr")
S_COEFFICIENT = 30.0  # EW001.yaml warped S coefficient
# PDG-2025 U-fixed ellipse (EW001 anchor).
FIT = ObliqueSTFit(
    s_central=0.026, t_central=0.047, sigma_s=0.075, sigma_t=0.066,
    rho_st=0.90, u_fixed=0.0, chi2_budget=CHI2_2DOF_95, confidence_level=0.95,
)
M_GRID_TEV = np.linspace(1.0, 40.0, 39001)  # 1 GeV steps


def oblique_floor_tev(ew_model: str) -> float | None:
    excluded = np.array([
        not evaluate_rs_oblique_proxy(
            m_kk_gev=float(m * 1000.0), fit=FIT,
            s_coefficient=S_COEFFICIENT, ew_model=ew_model,
        ).passes
        for m in M_GRID_TEV
    ])
    idx = np.where(excluded)[0]
    return float(M_GRID_TEV[idx[-1]]) if len(idx) else None


def zbb_delta_gL_b(ew_model: str) -> float:
    # Reuse the custodial-test build harness for a representative quark fit.
    from test_rs_ew_custodial_pr1 import _build_point, _sample_fit  # noqa: E402
    pt = _build_point(_sample_fit(), mkk_gev=3000.0, ew_model=ew_model)
    return float(pt.extras["rs_ew_couplings"].z_delta_g_L_d[2, 2].real)


def main() -> None:
    print("=== SU(2)_R vs P_LR decomposition ===\n")
    print(f"{'EW model':22s}{'oblique S,T,U floor':>22s}{'delta_g_L^b (Z b_L)':>24s}")
    for m in MODELS:
        f = oblique_floor_tev(m)
        bL = zbb_delta_gL_b(m)
        floor = f"{f:.2f} TeV" if f is not None else "<1 TeV"
        print(f"{m:22s}{floor:>22s}{bL:>24.3e}")
    print(
        "\nSU(2)_R alone pulls the oblique floor down (T cure); P_LR alone zeroes "
        "delta_g_L^b (Z b_L cure).  The two are independent."
    )


if __name__ == "__main__":
    main()
