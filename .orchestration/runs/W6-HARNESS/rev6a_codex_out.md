1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: local builder speedup was `3.98-4.02x`, not author `11.34x`; overlap microbench was `42.3x`, so the perf hook is real but timing is environment-dependent.
4. Equivalence: recursive all-extra compare max-rel-diff `4.038244873e-06` (`rs_charged_current`); all listed extras `<1e-3`.
5. Spline accuracy: targeted test `a=6.056183e-07`, `Omega=3.766492e-08`; denser off-knot probe `a=8.792121e-07`, `Omega=1.513259e-06`, all `<<1e-3`.
6. Default unchanged: no-injection rebuild path still uses direct spectrum/overlap behavior; full suite passed.
7. Fail-closed: rel_tol, domain, max_modes, and spectrum-identity mismatches all raised.
8. Perf scope: splined expensive numerical `a(c)`/`Omega_n(c)`; closed-form `B(c)` remains direct.
9. Isolation: code diff limited to `rs_ew_spectrum.py`, `rs_ew_builder.py`, `point_builder.py`; new W6a test present; untracked W6 harness prompt/out files also present.
10. Pytest: `1700 passed, 1 skipped in 790.69s`.

W6A-OK