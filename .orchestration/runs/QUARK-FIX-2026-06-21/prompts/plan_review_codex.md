You are an INDEPENDENT adversarial reviewer (Codex/gpt-5.4, xhigh) in a dual-signoff gate.
An Opus agent authored a fix plan for three BLOCKER-level quark-sector bugs. Your job: try to
BREAK the plan. Do not trust it. Re-derive the disputed physics yourself from primary sources
and verify every claimed code location against the actual files. You have read-only intent —
do NOT modify production code; only write your review file.

Repo root: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing

## Read
- THE PLAN under review: `.orchestration/runs/QUARK-FIX-2026-06-21/PLAN.md`
- The audit it implements: `.orchestration/runs/QUARK-AUDIT-2026-06-10/SYNTHESIS.md` and
  slices `slice-1-deltaf2.md`, `slice-2-model-couplings.md`, `slice-3-rsew.md`, `slice-4-adapters.md`, `slice-5-harness.md`
- The downstream notes: `review_local/z_to_bb_review.tex`, `epsilon_k_review.tex`,
  `deltaF2_framework_review.tex`, `constraint_formulas.tex`
  NOTE: the plan asserts these notes DOCUMENT THE BUGGY CODE and are NOT independent authority.
  Independently decide whether that is correct.
- The code to be changed: `quarkConstraints/fit.py` (B2), `quarkConstraints/deltaf2.py` +
  `quarkConstraints/modern/phenomenology.py` (B3), `quarkConstraints/rs_ew_couplings.py` (B1),
  and the tests/adapters they feed.

## Stress-test these specifically (the plan author flagged them)
1. **B3 direction (highest stakes).** The whole ε_K ×1.72 correction rides on this. Independently
   re-derive the ΔS=2 LR operator matrix elements (GGMS / Ciuchini–Franco–Masiero–Silvestrini
   "CFW" basis). Is the colour-SINGLET O4 really the LARGE coefficient (R/4 + 1/24) and the
   crossed O5 the small (R/12 + 1/8)? Is the code genuinely (a) basis-swapped AND (b) missing
   the 1/(2m_M) state normalization (factor 2)? Does the review_local note's "R5 ratio = 1.0000
   vs CFW" actually validate the O4/O5 split, or only the VLL sector (as the plan claims)?
   Confirm ALL SIX matrix-element sites are located (2 in deltaf2.py + 4 in modern/phenomenology.py).
2. **B1 Zbb.** Verify the CGHNP 0807.4937 retranslation: convention dictionary (F²=2 f_IR²,
   c_CGHNP = −c_repo), the wrong c-sign in 1/(1−2c), the missing factor 2 in F², and the
   flavour-sum index (1/F²(c_bR) for ALL light-gen terms, not 1/F²(c_d_i)). Is "correct it"
   (vs "drop it") the right call? Check the sign and ~190× magnitude claim.
3. **B2 ε_K rephasing.** Is rephasing the SVD columns to a PDG CKM convention correct and, crucially,
   DETERMINISTIC across a 1M paired scan (no np.angle(0) ambiguity, no sign-flip nondeterminism that
   would break minimal-vs-custodial draw pairing)? Would a rephasing-invariant ε_K construction be safer?
4. **Fix coupling/order.** The audit says B3's two sub-bugs partially cancel and must be fixed
   together; and B1/B2/B3 move floors in opposite directions. Is the plan's order safe? Are the
   literature-anchored ABSOLUTE test pins genuinely anchored (paper numbers), not snapshots of
   current (buggy) output?
5. **Scope calls.** Is it acceptable to defer M6 (scale pairing), the re-scan, and the
   STATE_OF_PROJECT restatement? Is treating M1 (T010 gate policy) as a separate sign-off right?

## Output
Write your full review to `.orchestration/runs/QUARK-FIX-2026-06-21/review_plan_codex_v1.md`.
End your STDOUT reply with a single line:
`VERDICT: APPROVE` or `VERDICT: REVISE` followed by a numbered list of BLOCKING issues
(each: location, what's wrong, what the plan must change) and a separate list of non-blocking suggestions.
Be terse in STDOUT; put detail in the file.
