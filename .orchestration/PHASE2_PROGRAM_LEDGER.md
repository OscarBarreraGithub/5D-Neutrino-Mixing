# Phase-2 Program Ledger — RS-EW rigor + full-catalog scan readiness

Durable state for the post-rebuild program. Survives orchestrator context resets.
**On resume/compaction, read this FIRST, then `rs_ew_sector_design.md` and `NEEDS_HUMAN_PHYSICS.md`.**
The constraint rebuild (103 constraints) is COMPLETE (see `REBUILD_LEDGER.md`). This program
makes the catalog's NEW-PHYSICS rigorous (close the G1/G2 proxy gap) so a definitive 100M+
cluster scan is meaningful, plus the items the triage marked "we build, not human input".

## ⛔ GOVERNANCE — the dual-signoff gate (NON-NEGOTIABLE, user mandate 2026-06-02)

The orchestrator (Claude main loop) makes **NO design decisions and writes NO production code**.
It only routes work, records verdicts, and maintains this ledger. EVERY work item — plan AND
implementation — must be independently approved by **BOTH a codex agent AND a Claude/Opus agent**.
Nothing is committed without dual APPROVE.

Per work item:
1. **PLAN** — a codex (gpt-5.x xhigh) authors an implementation plan, grounded in the approved design.
2. **PLAN REVIEW (dual)** — a *second* codex independently critiques it **and** an Opus agent independently
   critiques it. Route critiques back to the plan author; iterate until **both a codex and Opus APPROVE the plan**.
3. **IMPLEMENT** — codex implements per the approved plan (code + tests).
4. **IMPLEMENT REVIEW (dual)** — codex code/physics review **and** Opus independent review (re-derive numbers,
   check scaffold contract / isolation / honesty / determinism). Route fixes; iterate until **both APPROVE**.
5. **COMMIT** — one commit per item. Message records: plan approvers (codex+opus SHAs/verdicts), impl approvers,
   what changed. Update this ledger + push.
6. Orchestrator never overrides an agent verdict with its own opinion. Disagreements → another review round.

**Retroactive review (R-items):** anything committed since the rebuild began that did NOT get dual (codex+Claude)
review gets it now — codex review + Opus review of the existing code/doc; both APPROVE → mark "retro-OK"; any
finding → fix gate. Scope: from the rebuild start; includes a "similar implementation plan from last week" to be located.

Codex via `~/bin/codex_worker.sh` (CODEX_MAX_CONCURRENCY=6, own background task, never chained-heredoc).
Build prompts under `.orchestration/runs/<ITEM>/`. Keep orchestrator context lean: delegate all reading; read terse verdicts only.

## Scope split (from NEEDS_HUMAN_PHYSICS triage)
- **IN SCOPE (we build, dual-gated):** G1 RS-EW couplings, G2 lepton couplings, G4 CKM phase (B002/B004),
  exclusive form factors (B013/B014 from literature), re-wire affected constraints proxy→rigorous, full-catalog harness + smoke scan.
- **OUT (genuine human input — stays advisory/non-vetoing, surfaced for the user):** EDM rigor (E001/E004/E006–E009),
  ε′/ε (K003) + charm/nonleptonic CPV (C003/B032–B034) SM adoption, collider σ×BR recast scope (CR*). DO NOT fake these.

## Work items
| ID | Item | Plan (codex+opus) | Impl review (codex+opus) | Committed |
|----|------|-------------------|--------------------------|-----------|
| W1 | G1 design — opus∥codex drafts → consensus → dual-review (codex caught contact-formula BLOCKER: dropped light-Z g_q^SM·δg_l^LFV) → revise → **RE-REVIEW DUAL APPROVE** (codex DESIGN-OK + Opus DESIGN-OK). `rs_ew_sector_design_CONSENSUS.md` is the implementation-ready spec; 7 phases; 24→26 FULL/17 PARTIAL/7 HUMAN. | dual ✅ | n/a (design) | `4e5eff0` |
| W2-P1 | Derivation pins (`derivations/rs_ew_gauge_kk_coupling.tex`) — dual gate caught x_1=2.4048(ε→0) vs true gauge NN root ~2.45 + missing KK-sum truncation → fixed → **DUAL APPROVE** (codex+Opus PHASE1-OK). | dual ✅ | n/a (deriv) | `ffb4ad2` |
| W2-P2 | Spectrum+overlap kernel (RSEWSpectrum, kk_ew_mass=x_1·Λ_IR, Ω_n(c), cached a(c), a_ref) — codex author IN FLIGHT → dual review | author in flight | — | — |
| W2 | G1 implement, sub-phased per approved design (S1 KK gauge mass+a(c); S2 quark Z matrices; S3 lepton sector; S4 charged-current+oblique; S5 numeric oracle) | — | — | — |
| W3 | Re-wire constraints proxy→rigorous (fan-out by family; update NEEDS_HUMAN flags+tests) | — | — | — |
| W4 | G4 CKM phase (B002 sin2β, B004 φ_s) in-core | — | — | — |
| W5 | Exclusive form factors (B013, B014) from cited literature | — | — | — |
| W6 | Full-catalog cluster harness (sweep→point_builder→evaluate_all→serialized) + smoke scan | — | — | — |
| R1 | Scaffold hardening `02e2424`→`f82036a` — retro-review found 3 framework gaps (NaN/Inf accepted; load_anchor couldn't validate value_id/block_key/units/CL; mutably-shared extras). Fixed: finite/bool/Severity guards, optional anchor validators, immutable extras, reset_for_tests, TEMPLATE. **RETRO-OK + HARDENED** ✅ (dual: codex SCAFFOLD-FIX-OK + Opus SCAFFOLD-FIX-OK; 1054→1061 passed, backward-compat verified) | dual ✅ | dual ✅ | `f82036a` |
| R2 | ΔF=2 adapter running-wrappers `fd2f46a` — Opus-only | — | — | — |
| R3 | Complex-M12 phase helpers (B002/C002) `e08977d` — Opus-only | — | — | — |
| R4 | mu_e_conversion m_μ⁵ core fix `c6c949c` — **RETRO-OK** ✅ (codex MUE-OK + Opus MUE-OK; recomputed rate/BR matches to 15 sig figs; KKO m_μ⁵ dimensionally correct; L003/L004/L005 consistent; 1054 passed) | dual ✅ | n/a (review-only) | retro-OK |
| R5 | Pre-rebuild cores never re-reviewed in a gate: `quarkConstraints/deltaf2.py` (5206fc8), `scales.py` (c540830) | — | — | — |
| R6 | G1 design doc `rs_ew_sector_design.md` — folded into W1 cross-review | — | — | — |
| R7 | NEEDS_HUMAN_PHYSICS triage + PHASE3_SCAFFOLDING_PLAN.md (last-week plan, no codex sign-off) — doc-review | — | — | — |
| R8 | Orchestration/status docs (ledgers, wave-done commits) — low stakes, doc-review | — | — | — |
| -- | (tooling ~/bin/codex_worker.sh + codex_usage.sh: NOT in git, operator tooling — out of scope unless user includes) | — | — | — |

## ⏯️ CURRENT STATE / NEXT ACTION (updated 2026-06-02)
- Governance gate active. Discovery audit done → R-worklist populated (R1–R8 above). "Last week's plan" = `PHASE3_SCAFFOLDING_PLAN.md` (3292d68, no codex sign-off; R7).
- **R4 mu_e core: RETRO-OK** (dual ✅). **R1 scaffold: dual gate FAILED** (codex blockers vs Opus OK) → hardening fix IN FLIGHT (`by6va1dim`): NaN/Inf guard on ConstraintResult, optional load_anchor value_id/block_key/units/CL validation, immutable extras, reset_for_tests fix, TEMPLATE update — backward-compatible, suite must stay ≥1054.
- **W1 G1 design: consensus synthesis IN FLIGHT** (`boyrrvgq2`) merging opus+codex drafts → then dual-review (codex+opus) → approved design.
- NEXT when these land: (1) dual-review the R1 scaffold fix → retro-OK + commit; (2) dual-review the G1 consensus design → approved → begin W2 implementation sub-phases (each through the full gate: codex plan+impl → codex review + opus review → both APPROVE → commit). (3) Then remaining R-items R2/R3/R5/R7 (R8 low-stakes), then W3 re-wire, W4 CKM phase, W5 form factors, W6 harness+smoke→100M.
- Concurrency: codex ≤6 via wrapper; launch each as OWN background task. Watch for orphan codex procs after agent-spawned codex (kill wrappers first, then children).
