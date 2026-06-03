# W2 PHASE 2 — spectrum + overlap kernel (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Implement PHASE 2 of the APPROVED design `.orchestration/rs_ew_sector_design_CONSENSUS.md` (§5 item 2; plus the relevant schema/point_builder sections). First a SHORT plan, then implement code + tests. A codex reviewer and an Opus reviewer will independently verify (both must APPROVE). This is the FIRST production-code phase — get the kernel right; later phases build on it.

SCOPE (exactly Phase 2 — do NOT re-wire any constraints yet, that is Phase 3+):
Implement the RS-EW spectrum + overlap kernel as a NEW module (follow the design's file placement + repo conventions, e.g. a new `quarkConstraints/rs_ew_*.py` or the package the design names):
- `RSEWSpectrum`: gauge-KK masses `m_n` and profiles `chi_n` via `solvers/bessel.py` (`solve_kk("gauge","NN", exact=True)`); `kk_ew_mass_gev = x_1 * Lambda_IR` with x_1 the EXACT gauge NN root (~2.45 at physical eps), sourced numerically (NOT the literal 2.4048).
- Profile normalization; the zero-mode weight `w0(c,t)=g0(c,t)^2/t` consistent with `warpConfig/wavefuncs.py` (f_IR/f_UV) — must satisfy the endpoint identities `g0(c,1)^2=2 f_IR^2`, `g0(c,eps)^2=2 f_UV^2`.
- `Omega_n(c)` overlap integrals and the cached numerical overlap `a(c)=Sum_n (M_KK^2/m_n^2) chi_n(1) Omega_n(c)`, using the Phase-1 truncation criterion (double N=16,32,64,…; relative tolerance <1e-3; RAISE if not converged by N=512).
- `a_ref` handling (the universal piece to be subtracted at the coupling level in Phase 3).
Ground every formula in `derivations/rs_ew_gauge_kk_coupling.tex` (the dual-approved Phase-1 derivation) and §6 of the design.

TESTS (must independently verify, not just call-through):
- KK root: exact gauge NN root ~2.449–2.451 (NOT 2.4048); `kk_ew_mass_gev = x_1*Lambda_IR`.
- `a(c)` convergence: the truncation criterion actually converges to rel<1e-3 and raises past N=512 if forced.
- Endpoint identities g0(c,1)^2=2 f_IR^2 and g0(c,eps)^2=2 f_UV^2 to ~machine precision across several c.
- Universal-c precursor: for universal c across a species set, `a(c)-a_ref = 0` for all (the building block of delta_g=0).
- Sign precursor: IR-localized b_R (small c_R) gives `a(c_bR) > a(c_light)`.
- Determinism + finite (no NaN/Inf); caching does not change results.

CONSTRAINTS: do NOT modify the 103 constraints, the scaffold, or existing physics cores; ADD a new module + new tests only. `python -m pytest tests/ -q` (or the repo's runner) must stay green plus your new tests pass. Numeric outputs real finite floats.

OUTPUT (<=16 lines): short plan; the new module path + key public API (RSEWSpectrum, a(c), kk_ew_mass); the ACTUAL numbers your tests assert (exact root, a(c) sample + convergence N, endpoint-identity residuals, the universal-c=0 and b_R-sign checks); pytest counts. End with: PHASE2-DONE.
