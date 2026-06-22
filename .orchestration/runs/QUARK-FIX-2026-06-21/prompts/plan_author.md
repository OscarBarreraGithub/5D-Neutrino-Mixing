# PLAN AUTHOR — quark-sector audit fixes (B1/B2/B3 + entangled items)

You are the **plan author** in a dual-signoff gate. Your job is to produce a single,
complete, implementation-ready fix plan. **Do NOT modify any production code.** You may
read anything and you may run read-only commands. Think very hard (ultrathink); this plan
will be adversarially reviewed by an independent Codex (gpt-5.4 xhigh) agent and you will
iterate with it until both of you APPROVE.

## Background

A read-only adversarial audit of the quark sector (2026-06-10) found three independent
BLOCKER-level implementation errors plus several Majors. The headline RS floors (minimal
25–30 TeV, custodial strict 2–3 TeV, inclusive 7 TeV) are NOT defensible as computed. Your
plan must fix the code so the floors become trustworthy.

## REQUIRED READING (read all before planning)

Audit (authoritative finding list):
- `.orchestration/runs/QUARK-AUDIT-2026-06-10/SYNTHESIS.md`  ← start here
- `.orchestration/runs/QUARK-AUDIT-2026-06-10/STATUS.md`
- `.orchestration/runs/QUARK-AUDIT-2026-06-10/slice-1-deltaf2.md`   (B3, M6)
- `.orchestration/runs/QUARK-AUDIT-2026-06-10/slice-2-model-couplings.md` (B2)
- `.orchestration/runs/QUARK-AUDIT-2026-06-10/slice-3-rsew.md`     (B1, M1, M2)
- `.orchestration/runs/QUARK-AUDIT-2026-06-10/slice-4-adapters.md` (M6, anchors)
- `.orchestration/runs/QUARK-AUDIT-2026-06-10/slice-5-harness.md`  (M4, M5, floor prose)

Independent derivations / source translations (from a separate review effort):
- `review_local/z_to_bb_review.tex`        (B1 — Zbb / CGHNP 0807.4937)
- `review_local/epsilon_k_review.tex`      (B2 — ε_K rephasing)
- `review_local/deltaF2_framework_review.tex` (B3 — O4/O5 coefficients, normalization)
- `review_local/constraint_formulas.tex`   (cross-check of formulas as implemented)

Code to be changed:
- `quarkConstraints/fit.py`            (B2 — SVD column phases ~lines 245-258)
- `quarkConstraints/deltaf2.py`        (B3 — O4/O5 coeffs + 1/(2m_M) ~lines 719-724, 908-913)
- `quarkConstraints/modern/phenomenology.py` (B3 — vendored copies; audit says "3 copies")
- `quarkConstraints/rs_ew_couplings.py` (B1 — `_casagrande_zbb_B_profile` ~1835-1847,
  `build_rs_zbb_fermion_kk_mixing` ~940-971)
- relevant tests under `tests/` and the constraint adapters they feed (T010/T011, K001,
  B002/B004, C002, D0)

Context / governance:
- `docs/STATE_OF_PROJECT.md`, `.orchestration/PHASE2_PROGRAM_LEDGER.md`
- `CLAUDE.md` (conventions: c = M_5/k, f_IR/f_UV, v=174)

## What the plan must cover

For EACH fix item, the plan must give:
1. **Exact location(s)** — file + function + line range, and ALL duplicate/vendored copies
   (the audit flags multiple copies for B3; locate every one and list them explicitly).
2. **The wrong expression as currently coded** and **the corrected expression**, written in
   repo variables, with the derivation or source equation it comes from (cite the review_local
   .tex and/or CGHNP/CFW equation numbers). State the convention dictionary used
   (e.g. CGHNP F²=2 f_IR², c_CGHNP=−c_repo for B1).
3. **Numerical sanity check** to confirm at a representative point (give the expected
   before/after numbers from the audit, e.g. B1 δg_L^b sign flip and ~190×; B3 ε_K ×1.72).
4. **Interaction with other fixes** — the audit warns B3's two sub-bugs partially cancel and
   must be fixed TOGETHER; and that fixing B1 + B2 + B3 move floors in OPPOSITE directions.
   Make the ordering and any "fix-together" couplings explicit.

The plan must also decide and justify (these are YOUR design calls — the orchestrator makes none):
- **Scope & order.** The audit's recommended order is B2 → (B3+M6 together) → B1 → M1/M2 →
  M5 + tests. Confirm or revise, with reasoning. Decide whether M6 (bag/Wilson scale pairing),
  M1 (T010 1σ vs combined-fit gate policy), M2 (EW001 ΔT Λ_IR-vs-M_KK convention), and M5
  (tag_result substring bug) are IN SCOPE for this fix pass or explicitly deferred — justify each.
- **Test strategy.** The audit's M7 says the test suite pins code-against-itself and would
  never catch these bugs. The plan MUST add literature-anchored ABSOLUTE pins (a fixed known
  input → a number from the paper/derivation) for B1, B2, B3 — not snapshot-of-current-output
  pins. Specify the anchor values and their sources.
- **Re-scan / downstream.** The two 1M scans, the minimal-vs-custodial comparison, the Scan
  Explorer page, and the STATE_OF_PROJECT/methodology-note floor prose all depend on these
  numbers. State precisely what must be re-run/restated AFTER the code fix, and whether that is
  in scope for this pass or a documented follow-up. (Re-running 1M-point scans may be compute-
  heavy — flag the cost; the orchestrator will surface the go/no-go to the user.)
- **Commit structure.** Per the ledger, one commit per logical item with a message recording
  plan/impl approvers. Propose the commit breakdown.

## Output

Write the full plan to `.orchestration/runs/QUARK-FIX-2026-06-21/PLAN.md` (create it).
Structure it so a reviewer can check each item independently. Then RETURN to me (the
orchestrator) a CONCISE summary only: the fix list, your scope/order decision, the key
open questions you want the reviewer to stress, and a one-line "PLAN DRAFT v1 ready" — do
not paste the whole plan into your reply.
