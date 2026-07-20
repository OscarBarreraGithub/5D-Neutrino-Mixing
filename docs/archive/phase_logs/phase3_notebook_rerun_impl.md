# Phase 3 Notebook Re-execution — Implementation Report
**Date**: 2026-05-16
**Branch**: paper/quark-scan-2026q2

## Per-notebook summary

- `notebooks/dense_scan_2sigma_vs_1sigma_comparison.ipynb` — CLEAN. Commit `e2b83492df97579fa85393dbd4c36ae745bfa8ec`. Re-executed end-to-end with the requested nbconvert command. Key output cells: `2σ accepted points: 83,961`; `1σ accepted        : 83,958 / 83,961  (100.0%)`; `dropped 2σ → 1σ    : 0`; `saved fig2_mkk_bound_2sigma.{pdf,png}` and `saved fig2_mkk_bound_1sigma.{pdf,png}`. No notebook bug fix was needed.
- `notebooks/dense_scan_mkk_constraints_pdg2024.ipynb` — CLEAN. Commit `5471ccd9257629250b89abaceb8e37f3cf75078c`. Re-executed end-to-end with the requested nbconvert command. Key output cells: `records (total scanned)            :  100000`; `accepted (publication ξ_KK rescaled):   83961`; `data["M_KK"] range : 1224 … 48974 GeV`; `epsilon_K     34096 ( 91.4 %)`. No notebook bug fix was needed.
- `notebooks/pdg_quark_target_fix_verification.ipynb` — CLEAN. Commit `a4ee3bea34fcfe91c94c7cbddd77673a78f63195`. Re-executed end-to-end with the requested nbconvert command. Key output cells: `PDG edition  : PDG 2024`; `NEW GATE on synthetic point (m_s=0.065 GeV)`; `-> new verdict: REJECT`; `REJECTED by new strange gate alone :   940  (94.0%)`. No notebook bug fix was needed.
- `notebooks/rs_anarchy_analysis.ipynb` — FIXED-AND-CLEAN. Commit `2121e4e545827f45998df4d9f0a91ce4440378db`. Fixed the stale pre-audit run directory and histogram logic, then re-executed end-to-end with the requested nbconvert command. Key output cells: `Run dir: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scan_outputs/rs_anarchy_runA_20260515T085316`; `PDG-passing draws pooled across RUNA tiles: n = 1,532,640`; `M_KK_min percentiles (TeV, g_s*=3): p5=12.38, p50=47.26, p95=127.13`; `HEADLINE: M_KK at 50% acceptance among PDG-passing anarchic NMFV draws is 47.26 TeV at g_s*=3 (16.54 TeV at perturbative g_s≈1.05).`

## Bug fixes

- `notebooks/rs_anarchy_analysis.ipynb:69`: replaced stale historical pre-audit run path `scan_outputs/rs_anarchy_20260507T030811` with canonical post-audit RUNA path `scan_outputs/rs_anarchy_runA_20260515T085316` from the rc1 artifact manifest.
- `notebooks/rs_anarchy_analysis.ipynb:74`: redirected notebook-generated figures to `results/figures/quark/notebook_reruns` instead of a `scan_outputs/.../figures` subdirectory so the rerun does not write new scan-output artifacts.
- `notebooks/rs_anarchy_analysis.ipynb:308-335`: replaced the old single 10 TeV tile histogram with pooled streaming over all post-audit RUNA PDG-passing draws.
- `notebooks/rs_anarchy_analysis.ipynb:459-468`: added the `g_s*=3` rescaling and literal headline p50/p95 output cell, matching the rc1-locked numbers.

## Verification

- `dense_scan_2sigma_vs_1sigma_comparison.ipynb`: `2σ accepted points: 83,961`; `1σ accepted        : 83,958 / 83,961  (100.0%)`; `dropped 2σ → 1σ    : 0`.
- `dense_scan_mkk_constraints_pdg2024.ipynb`: `accepted (publication ξ_KK rescaled):   83961`; `data["M_KK"] range : 1224 … 48974 GeV`; `epsilon_K     34096 ( 91.4 %)`.
- `pdg_quark_target_fix_verification.ipynb`: `PDG edition  : PDG 2024`; `Strange exceeds its 2sigma window by ~14x  --  this is the failure mode.`; `REJECTED by new strange gate alone :   940  (94.0%)`.
- `rs_anarchy_analysis.ipynb`: `M_KK_min percentiles (TeV, perturbative g_s≈1.05): p5=4.33, p50=16.54, p95=44.50`; `M_KK_min percentiles (TeV, g_s*=3): p5=12.38, p50=47.26, p95=127.13`; `HEADLINE: M_KK at 50% acceptance among PDG-passing anarchic NMFV draws is 47.26 TeV at g_s*=3 (16.54 TeV at perturbative g_s≈1.05).`

## Final state

- 4 notebook commits pushed to origin.
- SHA range: `e2b83492df97579fa85393dbd4c36ae745bfa8ec..2121e4e545827f45998df4d9f0a91ce4440378db`.
- Notebook commit SHAs: `e2b83492df97579fa85393dbd4c36ae745bfa8ec`, `5471ccd9257629250b89abaceb8e37f3cf75078c`, `a4ee3bea34fcfe91c94c7cbddd77673a78f63195`, `2121e4e545827f45998df4d9f0a91ce4440378db`.
- All notebooks re-execute without error.
- Note: the shared paper branch received two intervening review commits between the third and fourth notebook commits; the four notebook commits above are the Phase 3 rerun deliverable.

===PHASE_3_NOTEBOOK_RERUN_IMPL_END===
