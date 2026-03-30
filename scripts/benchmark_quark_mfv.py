#!/usr/bin/env python
"""Benchmark the quark-sector MFV fit and lightweight diagnostics."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from quarkConstraints.benchmarks import default_quark_targets
from quarkConstraints.fit import fit_quark_sector
from quarkConstraints.proxies import summarize_flavor_diagnostics


def main() -> int:
    targets = default_quark_targets()
    solution = fit_quark_sector(targets, r=0.25, overall_scale=3.0, max_nfev=150)
    result = solution.result
    diagnostics = summarize_flavor_diagnostics(result)

    print("Quark-sector MFV benchmark")
    print(f"  success       = {solution.success}")
    print(f"  message       = {solution.message}")
    print(f"  nfev          = {solution.nfev}")
    print(f"  score         = {result.score:.6e}")
    print(f"  residual_norm = {result.residual_norm:.6e}")
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
        "  down alignment  = "
        f"{diagnostics.diagnostics.down_offdiag_ratio_in_q_basis:.6e}"
    )
    print(f"  up alignment    = {diagnostics.diagnostics.up_offdiag_ratio_in_q_basis:.6e}")

    max_mass_residual = np.max(
        np.abs(np.concatenate([result.mass_residuals_up, result.mass_residuals_down]))
    )
    max_ckm_residual = np.max(np.abs(result.ckm_residuals))
    if not solution.success:
        return 1
    if max_mass_residual > 1.0e-2 or max_ckm_residual > 1.0e-2:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
