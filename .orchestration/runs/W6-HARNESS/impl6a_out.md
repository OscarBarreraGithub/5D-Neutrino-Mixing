Plan executed: cache API, opt-in builder injection, equivalence/perf tests, full suite, review.

API: `RSEWOverlapSplineCache.build(...).a_spline(c)/omega_spline(c)` plus `spectrum=`, `a_of_c=`, `omega_of_c=`, `rs_ew_cache=` on both RS-EW builder paths.

Validation is fail-closed: injected cache/callables must match spectrum identity, modes, `rel_tol`, and `[0.3,0.9]` domain. Default rebuild path unchanged.

Spline max-rel-err: `a(c)=6.056183e-07`; `Omega_n(c)=3.766492e-08`.

Rebuild-vs-injected equivalence max-rel-diff: `4.038245e-06`.

Isolated per-point perf speedup: `11.338x`.

Splined: expensive numerical `a(c)` and `Omega_n(c)`. Left direct: closed-form Zbb/Higgs `B(c)` profiles.

Tests: focused RS-EW `40 passed`; full suite `1700 passed, 1 skipped`.

Review: Codex review approved; two independent agent reviews approved. Opus-specific selector was not exposed in this environment.

W6A-DONE.