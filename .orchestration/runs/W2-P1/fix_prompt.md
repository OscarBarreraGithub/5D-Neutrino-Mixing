# W2 PHASE 1 derivation FIX (codex, gpt-5.x xhigh). Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. Revise `derivations/rs_ew_gauge_kk_coupling.tex` to fix the two issues BOTH reviewers (codex + Opus) raised. No production code. Re-review (codex+opus) follows.

FIX 1 (both reviewers — the root pin is wrong): the note boxes `x_1 ≈ 2.4048`, which is the ε→0 J_0 zero, NOT the actual lightest gauge NN (++) KK root. At physical ε≈1e-15 the exact root from `solvers/bessel.py` (`solve_kk(...,exact=True)`, the cross-product `J_0(x)Y_0(εx)−J_0(εx)Y_0(x)=0`) is ≈2.45 (reviewers got 2.4487–2.4505 depending on ε). 
- Replace the boxed literal: state `x_1 = j_{0,1} + O(1/L)` with `j_{0,1}=2.4048` only the ε→0 limit, and the PRODUCTION value `M_KK = x_1 Λ_IR` must take `x_1` NUMERICALLY from `solve_kk("gauge","NN", exact=True)` (≈2.45 at physical ε), NOT the hardcoded 2.4048. Make clear Phase 2 sources the root from `solve_kk`, not a literal.

FIX 2 (codex SHOULD-FIX): the note defines `a_N` (KK-truncated overlap) but no truncation criterion. Add an explicit convergence/truncation rule for Phase 2: state how many KK modes N (or what relative tolerance on `a_N`) must be summed, given the `Σ 1/m_n²` (∝ 1/n²) convergence, so Phase 2 has a definite stopping rule.

Keep everything else intact (the overlap formula, w0↔f_IR/f_UV endpoint identities, the zpole additive-no-g_Z normalization, a_ref subtraction, the b_R sign/magnitude check, and the deferred-items section were all APPROVED by both reviewers — do NOT change them).

OUTPUT (<=10 lines): what you changed (the new x_1 statement + numerical sourcing, the truncation criterion), confirm pdflatex still compiles, confirm no other section changed. End with: PHASE1-FIX-DONE.
