1. BLOCKER: `RSEWSpectrum` has exact-float `_omega_cache/_a_cache` and `_cached_spectrum`, but `build_rs_ew_extras` calls `RSEWSpectrum.build(...)` directly; the plan must specify per-geometry spectrum reuse plus precomputed interpolation.
2. BLOCKER: continuous `c` means per-c memoization will not hit across points; production needs spline/lookup for `a(c)` and likely `Omega_n(c)` because charged-current contacts use full `omega(c)` vectors.
3. PERFORMANCE: local split is ~9.0s spectrum build + ~1.87s for 16 fresh `a(c)` calls; exact repeats are microseconds. Spectrum-cache-only still looks ~4-5s/point incl. 2.0s eval.
4. REALISTIC FIXED COST: with spectrum cache + `a/Omega` lookup + evaluator de-dup, expect ~1.2-1.8s/point; `1e8` is ~33k-50k core-hours. Feasible only with large cluster allocation and a hard smoke benchmark gate.
5. SHOULD-FIX: `evaluate_all` does not rebuild RS-EW matrices; semileptonic Wilsons are shared, but Delta-F=2 adapters recompute Wilson/evolution work across related constraints, so the 2.0s eval is partly reducible.
6. SHOULD-FIX: schema is mostly sufficient, but add explicit `hard_not_evaluated`/`coverage_complete`; `survives_all_HARD` must mean evaluated, active, vetoable HARD only, with missing/exception constraints separate from `excluded_by`.
7. BLOCKER: `NEEDS-HUMAN/PARTIAL` is not encoded in `Severity`; many such constraints are still `HARD`, so a machine policy must quarantine/advisory them or they will veto.
8. OK: advisory lists match live registry exactly: SOFT `EW002,EW003,K018,L003`; INFO `B032,B033,B034,C003,CR011,E004,E006,E007,E008,E009,K003,K013`; 103 constraints, zero import failures.
9. OK with caveat: all-include RS-EW builder plus quark KK-gluon merge covers declared extras, but `run_rs_anarchy.py` rows are not enough for `include_fermion_kk_mixing`; use real `QuarkFitResult`.
10. SHOULD-FIX: “reuse `run_rs_anarchy.py`” is only style reuse; current script lacks full catalog, leptons, fitted quarks, per-tile atomic resume, and optimized cache hooks.
11. OK: per-tile seed, atomic rename, resume hashes, and registry exception isolation are sound; universal-c SM sanity is the right physics guard.
12. OK: smoke plan `1e4 -> 1e5/1e6` is adequate if it verifies zero HARD coverage gaps, optimized sec/point, advisory policy, and exclusion rankings before `1e8`.
PLAN-NEEDS-FIXES.