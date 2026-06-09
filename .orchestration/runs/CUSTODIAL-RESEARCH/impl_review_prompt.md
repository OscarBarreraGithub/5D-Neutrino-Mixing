# Custodial RS-EW PR1 — INDEPENDENT REVIEW (codex, gpt-5.x xhigh). Verify the full uncommitted implementation.
Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Review the uncommitted custodial PR1 (`git diff` + new `tests/test_rs_ew_custodial_pr1.py`) against the DUAL-APPROVED plan `.orchestration/runs/CUSTODIAL-RESEARCH/plan_rev_codex_out.md` and spec `.orchestration/runs/CUSTODIAL-RESEARCH/CORRECTED_PRESCRIPTION.md`. Do NOT trust the author. Full suite is currently green (1723 passed, 1 skipped). Files: quarkConstraints/rs_ew_couplings.py, quarkConstraints/oblique_stu.py, flavor_catalog_constraints/physics_adapters/oblique_stu.py, flavor_catalog_constraints/primary/top_higgs_ew/{EW001,T010,T011}.py, rs_ew_builder.py, point_builder.py, scripts/run_full_catalog_scan.py (EW001 manifest), tests/test_rs_ew_custodial_pr1.py.

VERIFY (recompute / read code):
1. **Diagonal-only Zb_L zeroing**: custodial mode zeros ONLY protected diagonal down-left z_delta_l_d[0,0],[1,1],[2,2]; OFF-diagonal (i!=j) entries LEFT at minimal. Independently confirm with a custodial-vs-minimal build that off-diagonal z_delta_g_L_d (bs/bd/sd) and T014 results are UNCHANGED, only diagonals differ. (This was the prior blocking error - confirm it's fixed.)
2. **Residual**: z_delta_l_d[2,2] = kappa_b*(1/L)*minimal[2,2] (minimal = gauge+admixture SUM captured before zeroing); default kappa_b=0 => exact zero. Hermiticity preserved.
3. **Oblique**: ew_model plumbed rs_minimal_oblique_proxy->evaluate_rs_oblique_proxy->EW001 _resolve_ew_model (reads rs_ew_couplings.metadata, minimal fallback); custodial T = -pi/(4 cW^2 L); c_S UNCHANGED; U=0.
4. **Data-driven deferred flag**: T010/T011 custodial_variant_deferred=False when custodial active (tagged rigorous/active, not partial); top-partner loop stays deferred-flagged.
5. **Minimal byte-identity**: ew_model default => arrays + config hash + all minimal results identical (the existing quark-only/Zbb numbers unchanged).
6. **EW001 manifest**: rs_ew_couplings added to QUARK_ONLY_ALLOWLIST_EXTRAS["EW001"], confirmed quark-side/not-forbidden.
7. **Benchmark test**: a ~15 TeV point that minimal Z->bb vetoes now SURVIVES custodial (T010/T011 pass) - confirm the test asserts this and it's physically the intended relaxation.
8. Run `python -m pytest tests/ -q` (env: source ~/.bashrc && conda activate ising_bootstrap && export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"). Confirm green.

OUTPUT (<=14 lines): per-item PASS/CONCERN with code/numeric evidence; the custodial-vs-minimal off-diagonal/T014 unchanged check (numbers); minimal byte-identity confirm; pytest counts; any defect. END with EXACTLY ONE line: `CUSTODIAL-IMPL-REVIEW: APPROVE` or `CUSTODIAL-IMPL-REVIEW: NEEDS-FIXES`.
