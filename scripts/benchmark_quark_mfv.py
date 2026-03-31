#!/usr/bin/env python
"""Benchmark the quark-sector MFV fit and observable diagnostics."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from quarkConstraints.benchmarks import default_quark_targets
from quarkConstraints.couplings import compute_quark_kk_gluon_couplings
from quarkConstraints.deltaf2 import evaluate_delta_f2_constraints
from quarkConstraints.fit import fit_quark_sector
from quarkConstraints.proxies import summarize_flavor_diagnostics
from quarkConstraints.scales import (
    DEFAULT_QUARK_BENCHMARK_XI_KK,
    DEFAULT_QUARK_BENCHMARK_H_RS_MAX,
    DEFAULT_QUARK_PAPER_H_RS_MAX,
    default_quark_m_kk_from_lambda_ir,
)
from quarkConstraints.validation import benchmark_fit_summary


def main() -> int:
    targets = default_quark_targets()
    solution = fit_quark_sector(targets, r=0.25, overall_scale=3.0, max_nfev=150)
    result = solution.result
    xi_kk = DEFAULT_QUARK_BENCHMARK_XI_KK
    m_kk = default_quark_m_kk_from_lambda_ir(result.point.Lambda_IR, xi_KK=xi_kk)
    diagnostics = summarize_flavor_diagnostics(result, m_kk=m_kk)
    couplings = compute_quark_kk_gluon_couplings(result, M_KK=m_kk, xi_KK=xi_kk)
    deltaf2 = evaluate_delta_f2_constraints(couplings, M_KK=m_kk, xi_KK=xi_kk)
    fit_summary = benchmark_fit_summary(solution, xi_KK=xi_kk)

    print("Quark-sector MFV benchmark")
    print(f"  success       = {solution.success}")
    print(f"  message       = {solution.message}")
    print(f"  nfev          = {solution.nfev}")
    print(f"  score         = {result.score:.6e}")
    print(f"  residual_norm = {result.residual_norm:.6e}")
    print(f"  target label  = {targets.label}")
    print(f"  Lambda_IR     = {result.point.Lambda_IR:.2f} GeV")
    print(f"  xi_KK         = {xi_kk:.6f} (explicit repo bookkeeping convention)")
    print(f"  M_KK          = {m_kk:.2f} GeV")
    print()
    print("Target vs fitted up-sector masses [GeV]")
    for target, fitted in zip(targets.up_masses, result.masses_up):
        print(f"  target={target:10.6f}  fitted={fitted:10.6f}")
    print("Target vs fitted down-sector masses [GeV]")
    for target, fitted in zip(targets.down_masses, result.masses_down):
        print(f"  target={target:10.6f}  fitted={fitted:10.6f}")
    print()
    print("CKM observables [|Vus|, |Vcb|, |Vub|, J]")
    print(f"  target = {targets.ckm_observables}")
    print(f"  fit    = {result.ckm_observables}")
    print()
    print(f"  c_Q = {np.array2string(result.state.c_Q, precision=5)}")
    print(f"  c_u = {np.array2string(result.state.c_u, precision=5)}")
    print(f"  c_d = {np.array2string(result.state.c_d, precision=5)}")
    print(f"  F_Q = {np.array2string(result.state.F_Q, precision=5)}")
    print(f"  F_u = {np.array2string(result.state.F_u, precision=5)}")
    print(f"  F_d = {np.array2string(result.state.F_d, precision=5)}")
    print()
    print("Flavor diagnostics")
    print(f"  h_RS proxy      = {diagnostics.h_rs_proxy:.6e}")
    print(
        "  down misalign.  = "
        f"{diagnostics.diagnostics.down_misalignment_in_q_basis:.6e}"
    )
    print(f"  up misalign.    = {diagnostics.diagnostics.up_misalignment_in_q_basis:.6e}")
    print(
        "  down/up stress  = "
        f"{diagnostics.diagnostics.down_to_up_misalignment_ratio:.6e}"
    )
    print()
    print("Mass-basis KK-gluon couplings |g_ij|")
    for label, matrix in (
        ("left_down", couplings.left_down),
        ("right_down", couplings.right_down),
        ("left_up", couplings.left_up),
        ("right_up", couplings.right_up),
    ):
        print(
            f"  {label:10s} = "
            f"(|12|={abs(matrix[0, 1]):.6e}, |13|={abs(matrix[0, 2]):.6e}, |23|={abs(matrix[1, 2]):.6e})"
        )
    print()
    print("Delta F = 2 summary")
    for system in ("K", "B_d", "B_s", "D"):
        item = deltaf2.by_system[system]
        print(
            f"  {system:3s} ratio={item.ratio_to_bound:.6e} "
            f"pass={item.passes} dominant={item.dominant_operator}"
        )
    print(f"  overall pass    = {deltaf2.passes_all}")
    print()
    print("Validation gates")
    print(f"  benchmark h_RS < {DEFAULT_QUARK_BENCHMARK_H_RS_MAX:.2f} = {fit_summary.passes_proxy}")
    print(f"  paper h_RS < {DEFAULT_QUARK_PAPER_H_RS_MAX:.2f}     = {fit_summary.passes_paper_proxy}")
    print(f"  misalignment gate = {fit_summary.passes_misalignment}")
    print(f"  Delta F = 2 gate  = {fit_summary.passes_deltaf2}")

    if not solution.success or not fit_summary.passes_all:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
