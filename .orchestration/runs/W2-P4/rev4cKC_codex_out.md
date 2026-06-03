1. BLOCKER: none found; no outside tests broke or were silently dropped.
2. NIT: author’s scoped-count claim is wrong. Actual scoped tests are `57 -> 50` (-7), not `49 -> 50` (+1).
3. FULL-SUITE COUNT: `1645 -> 1638 passed` is fully explained by those five scoped files; `pytest tests/ -q`: `1638 passed, 1 skipped`.
4. Removed/added exactly: C006,C008 each `12 -> 10`; removed old missing/invalid proxy, zero_np, proxy-manual, optional_kk_ew_mass, safe/large param(2); added absent_rs, invalid_rs, diagonal-zero, live-scaling, amplified.
5. K019,K020,K021 each `11 -> 10`; removed old missing param(2), proxy-manual, optional_kk_ew_mass, safe/large param(2); added absent_rs, invalid_rs, diagonal-zero, live-scaling, amplified.
6. Isolation confirmed: only 4 LFV adapters + K019/K020/K021/C006/C008 + five tests changed; new helper has no tests; no YAML/anchor changes.
7. Rewire OK: constraints consume `rs_semileptonic_wilsons.lfv_llqq`; no `_wilson_prefactor`, no second `1/M_KK^2`; diagnostics mark both false.
8. Old `NEEDS-HUMAN-PHYSICS` proxy path is removed from production notes/diagnostics and replaced by tree-zero/deferred-loop status.
9. v1-zero OK: diagonal point gives all five `predicted=0`, `ratio=0`, `passes=True`, with LFV C9/C9p zero.
10. Live recompute OK: K019 `6.871e-17`, K020 `4.704e-19`, K021 `4.333e-18`, C006 `3.383e-18`, C008 `1.172e-16`.
11. Scaling OK: `M_KK 3000 -> 6000` gives rate ratio `3.90625e-3 = (3000/6000)^8`; amplified toys veto for all five.
12. Degradation OK: old-style/proxy point is non-vetoing with `evaluated=False`, `missing_extra=rs_semileptonic_wilsons`, `predicted=None`, `ratio=None`.
13. Public numeric fields checked as real floats for all five live evaluations.

P4CKC-OK