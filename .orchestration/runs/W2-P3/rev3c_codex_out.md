1. BLOCKER: none. Rewire matches plan at `T014.py:70` and `T014.py:462`: bs/bd/sd use `[1,2]`, `[0,2]`, `[0,1]` via the zpole adapter.
2. SHOULD-FIX: none. Proxy/NEEDS-HUMAN removed from T014 diagnostics; physics path reaches `zpole_down_fcnc_branching_fraction_from_couplings` only.
3. NIT: worktree has untracked `.orchestration/runs/W2-P3/impl3c_*` notes, but `git diff --name-status` code isolation is only `T014.py` + `test_T014.py`.
4. Recompute bs: `z_delta_g_L_d[1,2]=-2.038738022510e-04-2.203656917467e-08j`, `z_delta_g_R_d[1,2]=-1.968364441525e-04+0j`, `BR=1.320307590878e-07`.
5. SM-limit: `max|z_delta_g_L/R_d|=0`, `predicted=0`, `ratio=0`, `passes=True`, `evaluated=True`.
6. Degradation: old-style point gives `evaluated=False`, `missing_extra=rs_ew_couplings`, `predicted=None`, `ratio=None`, `passes=True`.
7. Limits/policy intact: pure-NP SM-zero/direct-width policy remains; all three YAML limits stay `< 2.9e-3`.
8. Test-count OK: T014 `11 -> 11`; proxy tests replaced by SM-limit, direct bs recompute, old-style absent/invalid, no-rescale, direct finite-field/purity coverage.
9. Pytest: `1646 passed, 1 skipped in 816.91s`.

P3C-OK