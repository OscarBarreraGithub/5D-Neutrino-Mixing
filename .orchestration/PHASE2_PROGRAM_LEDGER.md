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
| W2-P6 | Fermion-KK/Higgs **COMPLETE** (dual-approved). 6a Zbb fermion-mixing->T010/T011 FULL `85c06ea`; 6b Higgs-LFV T018/19/20 (this commit). Custodial/top-partner Zb_L = NEEDS-HUMAN; Higgs-LFV zero in diagonal v1. 1692 passed. **>> ENTIRE W2 RS-EW BUILD (P1-P6) COMPLETE <<** | plan dual ✅ | all sub-steps dual ✅ | DONE |
| W3 | Re-wire constraints proxy→rigorous (fan-out by family; update NEEDS_HUMAN flags+tests) | — | — | — |
| W4 | G4 CKM phase (B002 sin2β, B004 φ_s) in-core | — | — | — |
| W5 | Exclusive form factors (B013, B014) from cited literature | — | — | — |
| W6 | Full-catalog scan harness. **PLAN DUAL-APPROVED** (perf: a(c)+Ω spline + per-tile spectrum injection hook => ~1.2-3 s/point => 1e8 ~4-8 core-yr feasible; honest veto: tag rigorous|proxy|partial|stub, survives_all_HARD strict+inclusive, excluded_by_rigorous vs _proxy). Sub-steps: W6a builder spectrum/spline-injection hook (code change) → W6b harness driver+tagging+checkpoint → smoke 1e4-1e6 (post-cache gate) → 1e8. **W6a IMPL IN FLIGHT.** | plan dual ✅ | W6a in flight | — |
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

## ⏯️ CURRENT STATE / NEXT ACTION (updated 2026-06-04 — pre-compaction snapshot)
- Governance gate active: codex AND opus must APPROVE every plan+impl; orchestrator only ROUTES — NO design decisions, NO production code. The gate caught a real defect at nearly every stage.
- **DONE (all dual-approved + committed) — HEAD `7956106`, full suite 1700 passed:**
  - **W2 RS-EW BUILD COMPLETE (P1→P6):** P1 derivation `ffb4ad2`; P2 spectrum/overlap kernel `389de37`; P3 quark-NC (3a builder `7813e6c` / 3b Z-pole `7c0ecff` / 3c FCNC-Z `651389e` / 3d-B rare-B `4834d50` / 3d-K rare-K `f04ae1c` / 3d-C rare-charm `0e0dc0c`); P4 lepton-NC (4a `7265c8d` / 4b Z-LFV `a585265` / 4c-KC `6665ccb` / 4c-L `711b56d` / 4d νν `4eac4ce`); P5 charged-current (5a `2b89ae2` / 5b `3f7762a` / 5c `ec071c8`); P6 fermion-KK/Higgs (6a Zbb `85c06ea` / 6b Higgs-LFV `ba1ed58` + test-helper fixup `3c52b3e`).
  - **Retro R1–R5 all dual-OK:** R1 scaffold hardened `f82036a`; R4 mu_e retro-OK; R2/R3 dual-OK; **R5 ΔF=2 factor-2 scare → RESOLVED = CORRECT** (deltaf2 matches published Csaki-Falkowski-Weiler to machine precision; `me_vll=(2/3)f²mB` + Wilson `/6` is the right ½(prop)×⅓(color-Fierz) convention — see the "R5 RESOLVED" section below). **DO NOT "fix" this — a ×2 would BREAK agreement with the literature.**
  - **W4 CKM phase** `c5970a7` (B002/B004 compute 2β/φ_s in-core; NEEDS-HUMAN removed). **W6a** a(c)/Ω spline + per-tile spectrum-injection hook `7956106` (~4-11x builder speedup; injected path == rebuild to 4e-6; default path byte-identical).
- **W6b DONE + dual-approved + committed `df947f3`** (codex APPROVE + Opus APPROVE, both substantive w/ line cites; reviews in `.orchestration/runs/W6-HARNESS/rev6b_*`). Harness `scripts/run_full_catalog_scan.py` + `tests/test_full_catalog_scan_harness.py` (6 pass); full suite 1706 passed/1 skipped. Smoke (10k draws, `.orchestration/runs/W6-HARNESS/smoke_w6b_fresh/`): **533 evaluated / 9,467 skipped (94.7% `nonperturbative_lepton_yukawa`)**; 0.386 s/draw, **7.24 s/evaluated point**; 1e8 ≈ 1.1e4 core-h (draw basis) / 2e5 core-h (evaluated basis); universal-c SM sanity CLEAN (rigorous=[], proxy=[]); top rigorous vetoes B022/K004/L001/B003/B004/B023/K001; top proxy CR009/CR006/B016/CR001/CR005/CR012/CR013/EW001/K010/B015; **5 NEEDS-HUMAN constraints E001/E002/L006/L010/L023 throw TypeError on 100% of points → honestly tagged `stub`/non-vetoing/`hard_not_evaluated` (no proxy input on ParameterPoint — genuine human residual, NOT a bug).**
- **AWAITING USER GO/NO-GO on the 100M run** (presented 2026-06-04). Three decisions are the user's, not mine: (a) is "100M points" 100M *draws* (~5.3M evaluated, ~1.2 core-yr — trivial) or 100M *evaluated* (~1.9e9 draws, ~23 core-yr)? (b) tune priors to cut the 94.7% nonperturbative waste first (physics decision — intended prior)? (c) accept the 5 NEEDS-HUMAN constraints as declared coverage gaps, or wire proxy inputs first? Perf 7.24 s/evaluated pt MISSES the plan §18 1.2–3 s/pt gate → plan says "defer 1e8 + profile"; feasible regardless on draw basis. **Do NOT launch the 100M run without explicit user go.**
- **NEXT after user go:** the 100M scan. Then: **W5** exclusive form factors B013/B014 (DEFERRED — needs user sign-off on source; default Bharucha-Straub-Zwicky LCSR; secondary, stays PARTIAL); **R7/R8** doc dual-review (last-week plan + triage + status docs); **P7** loop residual (dipoles/EDMs) = genuine NEEDS-HUMAN, DOCUMENT only (no build).
- **HOW TO RUN A WORK ITEM (the gate, verbatim):** write prompt under `.orchestration/runs/<ID>/`; codex via `~/bin/codex_worker.sh` (`CODEX_TIMEOUT=<s> CODEX_MAX_CONCURRENCY=6 bash ~/bin/codex_worker.sh --out <file> -C "$PWD" "$(cat prompt)"`, as its OWN `run_in_background` Bash, NEVER a chained-heredoc). For each plan+impl: codex authors → codex review ∥ Opus review (spawn Opus as Agent model=opus, instruct "work SYNCHRONOUSLY, run pytest in foreground, END with the exact verdict line") → route fixes to author → iterate until BOTH a codex AND an Opus return their OK verdict → commit ONE item, update this ledger, push. Reviewers must INDEPENDENTLY recompute, not trust the author's self-report (codex authors often can't self-run Opus — the orchestrator runs the real dual gate).
- **HUMAN-INPUT items open (surfaced, NOT decided — stay advisory/non-vetoing in the scan):** custodial/BKT & top-partner Zb_L; dipole-loop normalization (μ→eγ, b→sγ NP, top dipoles); EDM basis/CP-phases/matrix-elements; ε′/ε & charm/nonleptonic CPV SM adoption; collider σ×BR recast scope; (minor) W5 form-factor source.
- **OPERATIONAL HAZARDS:** codex ≤6 concurrent (flock). Orphan-codex after an agent spawns its own codex: the wrapper watchdog RESPAWNS codex on kill — kill the wrapper(s) FIRST, then children, iterate to zero (`ps -eo cmd|grep '[c]odex exec'` should hit 0). The old generic "W2"/"W3" ledger rows are SUPERSEDED by the W2-P1..P7 / per-item rows.

## ⚠️ CRITICAL OPEN ITEM (R5) — possible factor-2 in ΔF=2 M12^NP normalization (2026-06-03)

Retro-review of `quarkConstraints/deltaf2.py` (the core behind the 5 fully-rigorous anchors K001/K002/B001/B003/C001 + B002/C002) surfaced a possible factor-2 in the NP M12 normalization. UNRESOLVED — agents split, internally inconsistent:
- The constraint path: B003 → `bs_mixing_from_wilsons_with_running` → `compute_m12_np` (deltaf2.py:987,931) → `M12 = c1_vll*me_vll + ...`, `me_vll=(2/3) f² m B` (line 910). Budget = experimental Δm/2 (B003.py:339).
- **HALF (bound 2× too loose) camp** [Opus decisive+tiebreaker, codex SURGICAL]: `me_vll=(2/3)f²mB` is half the standard `⟨O1⟩/(2m)=(4/3)f²mB`; feeding an SM-calibrated C1 through `compute_m12_np` gives Δm=5.4e-12 (half of 1.087e-11). So |M12^NP| is 2× small vs a full-convention budget → NP bound 2× too LOOSE.
- **CORRECT camp** [codex decisive+tiebreaker, Opus surgical]: deltaf2's OWN SM path reproduces experiment-agreeing Δm_s≈1.09–1.14e-11; self-consistent. The "half" tests injected an EXTERNAL standard-convention C1 that doesn't match deltaf2's own C1 matching (line 342, `/6`), so they tested a mismatched C1×me product.
- CRUX: is deltaf2's RS NP Wilson `c1_vll` (line 342) matched in the SAME convention as `me_vll` (line 910), so their PRODUCT M12^NP satisfies Δm=2|M12|? Convention reconstruction is unreliable (agents quoted (1/3),(2/3),(4/3),(8/3) for ⟨O1⟩ across runs). NEEDS the definitive independent cross-check below.
- DEFINITIVE TEST (pending): compute the RS KK-gluon Δm_K (and Δm_Bs) for a benchmark RS point via an EXPLICIT textbook RS-flavor formula (cited paper), compare to deltaf2's 2|M12^NP| for the SAME point. That bypasses the C1-vs-me convention split.
- IMPACT IF REAL: ΔF=2 NP bounds 2× too loose → the 5 anchors admit points they should exclude in the 100M scan. MUST resolve before W6 scan.

### ✅ R5 RESOLVED — ΔF=2 normalization is CORRECT (2026-06-03, dual-confirmed)
Convention-independent cross-check (the coefficient-reconstruction split was settled by physical comparison to PUBLISHED RS formulas):
- **Opus**: deltaf2 `2|M12^NP|` vs Csaki-Falkowski-Weiler 0804.1954 Eq.3.3-3.7 (KK-gluon t-channel + exact SU(3) Fierz ∑T^A⊗T^A=½(δδ-δδ/N_c)→(1/3)Q1) = **ratio 1.0000 to machine precision**.
- **Codex**: deltaf2 vs Blanke et al. 0809.1073 Eq.4.33 = ratio ~1.69 (basis/running diff, "NOT the 0.5× failure mode; ×2 would be the WRONG fix").
- RESOLUTION: deltaf2's `1/6` Wilson factor = ½(propagator)×⅓(color-Fierz); `me_vll=(2/3)f²mB` is the CORRECT matrix element IN THAT convention (color factor in the Wilson, not the ME). The "half" verdicts wrongly imposed the textbook ⟨O1⟩=(8/3)/(4/3) convention (color in the ME). **deltaf2 + the 5 anchors (K001/K002/B001/B003/C001) + B002/C002 are correctly normalized. NO FIX.**
- DOC HAZARD (recorded, not a bug): the repo's `(2/3)`-ME/`/6`-Wilson convention must NEVER be mixed with textbook `(8/3)`-ME-convention Wilsons. The `paper_0710_1869` path uses the SAME `/6` (matching_kkgluon.py:523) and is test/script-only; no live constraint mixes them.
- LESSON: the gate's value here was forcing the convention-independent published-formula cross-check; the orchestrator's refusal to blind-×2 on a split vote prevented BREAKING the bedrock.

**R2/R3/R5 ALL RETRO-OK (dual). Remaining retro: R7/R8 (docs). Then W4/W5/W6.**
