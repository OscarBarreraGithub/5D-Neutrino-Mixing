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
| W2-P2 | Spectrum+overlap kernel `rs_ew_spectrum.py` — dual gate: codex caught a KK-TOWER root-extraction BLOCKER Opus missed (n6=21.3 vs physical 18.1, n9→n10 regression) → fixed (ordered bracketed scan) → both independently recomputed tower+a(c) → **DUAL APPROVE** (1645 passed) | dual ✅ | dual ✅ | `389de37` |
| W2-P3 | Quark NC **COMPLETE** (dual-approved). 3a builder `7813e6c`; 3b Z-pole `7c0ecff`; 3c FCNC-Z `651389e`; 3d-B rare-B `4834d50`; 3d-K rare-K `f04ae1c`; 3d-C rare-charm (this commit). Z-pole/FCNC-Z/rare-B-K-charm all rigorous; graceful degradation everywhere; 1646 passed. | plan dual ✅ | all sub-steps dual ✅ | DONE |
| W2 | G1 implement, sub-phased per approved design (S1 KK gauge mass+a(c); S2 quark Z matrices; S3 lepton sector; S4 charged-current+oblique; S5 numeric oracle) | — | — | — |
| W2-P4 | Lepton NC **COMPLETE** (dual-approved). 4a lepton builder `7265c8d`; 4b Z-LFV `a585265`; 4c-KC LFV-rare-K/charm `6665ccb`; 4c-L LFV-leptonic `711b56d`; 4d nunu (this commit). Tree LFV=0 for diagonal v1 fit (loop deferred P7); nunu bites; all proxy lepton NEEDS-HUMAN resolved; 1664 passed. | plan dual ✅ | all sub-steps dual ✅ | DONE |
| W2-P5 | Charged-current **COMPLETE** (dual-approved). 5a builder `2b89ae2`; 5b rewire EW002/K018/K017/B009/B025 `3f7762a`; 5c EW003 (this commit). v1 near-SM (universal absorbed); EW002 SOFT, B025 PARTIAL, EW003 data-level; 1677 passed. | plan dual ✅ | all sub-steps dual ✅ | DONE |
| W2-P6 | Fermion-KK/Higgs. **PLAN DUAL-APPROVED** (arbitrations: minimal-RS Casagrande ZMA default [no full tower in repo]; top-partner/custodial Zb_L = NEEDS-HUMAN [dominant but rep-dependent]; T010/T011 FULL only minimal_rs_tree + custodial_variant_needs_human; Higgs-LFV zero in diagonal v1 fit). Sub-steps 6a Zbb fermion-mixing+T010/T011 → 6b Higgs-LFV(rs_higgs_yukawas)+T018/19/20. **6a IMPL IN FLIGHT.** | plan dual ✅ | 6a in flight | — |
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
- Governance gate active and working (it caught a real defect at EVERY stage: 3 scaffold blockers, design LFV cross-term, P1 KK-root limit error, P2 KK-tower bug — several flagged by codex after Opus had OK'd).
- DONE (dual-approved + committed): **W1** `4e5eff0`; **W2-P1** `ffb4ad2`; **W2-P2** kernel `389de37`; **W2-P3 3a** builder `7813e6c`; **3b** Z-pole `7c0ecff`; **3c** FCNC-Z `651389e`; **R1** scaffold `f82036a`; **R4** mu_e retro-OK. Suite 1646 passed.
- IN FLIGHT: **W2-P3 3d-B** rare-B vector rewire (`bpcayesyw`). Remaining P3: 3d-K (rare kaon), 3d-C (rare charm). Then P4 lepton NC, P5 charged-current, P6 fermion-KK/Higgs, P7 loops.
- NEXT after P3: W2-P4 lepton NC → P5 charged-current → P6 fermion-KK/Higgs (24→26) → P7 loops (mostly human-deferred). Then retro R2/R3/R5/R7 (R8 low-stakes). Then W4 CKM phase, W5 form factors, W6 full-catalog harness + smoke scan → 100M.
- HUMAN-INPUT items still open (surfaced, NOT decided): custodial/BKT model choice; dipole-loop normalization; EDM basis/CP/matrix-elements; ε′/ε & charm/nonleptonic CPV SM adoption; collider σ×BR recast scope.
- Concurrency: codex ≤6 via wrapper; launch each as OWN background task (never chained-heredoc). Watch for orphan codex after agent-spawned codex (kill wrappers before children). NOTE row "W2" (generic) is superseded by W2-P1..P7.
