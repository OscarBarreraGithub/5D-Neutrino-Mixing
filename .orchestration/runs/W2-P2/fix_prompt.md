# W2 PHASE 2 FIX (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. The dual review of `quarkConstraints/rs_ew_spectrum.py` returned a BLOCKER: the gauge-KK multi-root TOWER is wrong (the first root x_1=2.4505 is correct, but higher modes are mis-extracted). Fix it. A codex reviewer and an Opus reviewer will re-verify (both must APPROVE), and they will independently recompute the tower + a(c) — so the fix must be genuinely correct, not just matched to the old values.

THE BUG (from review):
- The KK tower roots beyond n=1 are wrong/duplicated/skipped. Evidence: module gives n6 root = 21.261, but the correct 6th gauge-NN root is ≈18.12; and n10 (27.544) REGRESSES below/duplicates after n9 (30.685) — a clear root-finding ordering/duplication failure.
- Physics check: successive roots of the gauge NN cross-product `J_0(x)Y_0(εx)−J_0(εx)Y_0(x)=0` are spaced ≈π apart, so x_n ≈ x_1 + (n−1)π ≈ 2.45, 5.6, 8.7, 11.8, 15.0, 18.1, 21.3, … (asymptotically (n−1/4)π-ish). The module's tower must reproduce this: strictly increasing, unique, no skips/dupes.
- This contaminates a(c): correct-root Q=1024 gives a(0.2)=21.8605, a(0.65)=−1.4571 (vs the current wrong module 21.9165, −1.4761). After the fix, a(c) should match the CORRECT values (recompute and report).

FIX:
1. Correct the KK-tower root extraction in `rs_ew_spectrum.py` so `m_n = x_n·Λ_IR` are the correctly-ordered, unique successive roots of the gauge NN cross-product (use `solvers/bessel.py`'s tower API correctly, or bracket successive roots at ≈π spacing; ensure no skipped/duplicated roots and strict monotonic increase). Confirm x_1=2.4505 unchanged.
2. Recompute the cached `a(c)` with the corrected tower; update any pinned test values to the CORRECT numbers.
3. **Test (codex SHOULD-FIX):** `tests/test_rs_ew_spectrum.py` must add explicit assertions that the tower is strictly increasing + unique (no dupes/regressions) and that the first ~6 roots match the ≈(n−1)π-spaced expected values (e.g. n6≈18.1, NOT 21.3). The overlap test must recompute the tower independently (bracketed roots), not call through the module's `omega()/chi_n()`, so it cannot lock in a bad tower again.

CONSTRAINTS: ADD/fix only this module + its test; do not touch the 103 constraints/scaffold/other cores. `python -m pytest tests/ -q` must stay green (≈1644+). Numeric outputs real finite floats.

OUTPUT (<=12 lines): the root-finding fix (file:line), the corrected first ~6 roots + a(0.2)/a(0.65), confirmation the tower is now strictly-increasing/unique with the new test, pytest counts. End with: PHASE2-FIX-DONE.
