W6a — RS-EW builder perf hook: a(c)/Omega spline + spectrum injection (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement the perf-foundation of the DUAL-APPROVED W6 plan (.orchestration/runs/W6-HARNESS/plan.md): make per-point RS-EW building fast by reusing a per-tile prebuilt spectrum + a precomputed a(c)/Omega spline. First a SHORT plan, then implement code+tests. Codex + Opus dual-review (both must APPROVE). This touches the RS-EW core — the injected path MUST produce IDENTICAL couplings to the rebuild path (physics unchanged).

BUILD:
1. a(c)/Omega SPLINE: add a helper (on RSEWSpectrum or a sibling) that, given a built RSEWSpectrum, precomputes a CubicSpline (scipy.interpolate.CubicSpline) of a(c) over c in [0.3,0.9] (and Omega_n(c) vectors if charged-current/contacts need them) on a dense grid, with a verified max-rel-err < 1e-3 vs direct a(c) (target ~1e-6). Expose a fast `a_spline(c)` callable. One-time build cost ~tens of seconds; per-eval ~microseconds.
2. INJECTION HOOK: add OPTIONAL `spectrum=` and `a_of_c=` (and/or a combined `rs_ew_cache=`) kwargs to `flavor_catalog_constraints.rs_ew_builder.build_rs_ew_extras` (and the public `build_from_rs_ew_inputs` path) so that when provided, it REUSES the prebuilt spectrum + spline instead of calling RSEWSpectrum.build()/direct a(c). When NOT provided, behavior is UNCHANGED (existing tests + the rebuild path keep working).
3. CRITICAL EQUIVALENCE: the injected (spline) path must yield rs_ew_spectrum/rs_ew_couplings/rs_charged_current/rs_higgs_yukawas/lepton_mass_basis_couplings IDENTICAL to the rebuild path within the spline tolerance (assert all z_delta_g matrices, contacts, Wilsons agree to <1e-3 rel, ideally <1e-5). Physics MUST NOT change.
4. The Zbb fermion-mixing B(c) and Higgs B(c) also use the (1/(1-2c))(1/F^2-1+F^2/(3+2c)) form — if they're cheap (closed form, no integral) leave them; only a(c)/Omega (the expensive numerical overlaps) need splining. State which pieces are splined vs left.

TESTS:
- spline accuracy: max rel err of a_spline(c) vs direct a(c) over c-grid < 1e-3 (report actual).
- equivalence: build a point via rebuild vs via injected spectrum+spline; assert all extras' matrices/Wilsons agree to tolerance (report max rel diff).
- perf: time N points via rebuild vs via 1 prebuilt spectrum+spline reused over N; report the per-point speedup (expect >>10x once amortized).
- existing suite unchanged: `python -m pytest tests/ -q` stays green (the rebuild path is default).

CONSTRAINTS: additive/opt-in; do NOT change default behavior or physics; numeric outputs real finite floats. Touch rs_ew_spectrum.py + rs_ew_builder.py (+ point_builder if needed) + new tests only.

OUTPUT (<=14 lines): short plan; the spline API + injection kwargs; spline max-rel-err; rebuild-vs-injected equivalence max-rel-diff; per-point speedup; what's splined vs left; pytest counts. End with: W6A-DONE.
