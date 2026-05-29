# Flavor-Constraint Rebuild — Orchestration Ledger

Durable state for the from-scratch constraint implementation. Survives
orchestrator context resets. Updated after every per-constraint gate.

## Contract (decided 2026-05-28)

- **Scope:** Build REAL, computing physics for every catalog constraint
  (~102). Where repo machinery is missing (ΔS=1 RG, EDM loops, collider
  recasts, LFV, …), agent1 builds it. Multi-session project.
- **Per-constraint pipeline (serial, one at a time):**
  1. **agent1** = codex gpt-5.5 xhigh — plan, then implement (code + tests).
  2. **agent2** = codex gpt-5.5 xhigh — physics fact-check vs sources/intuition.
     **agent3** = codex gpt-5.5 xhigh — code review + numerical/unit-test verification.
     agent2 ∥ agent3 run in PARALLEL.
  3. Issues from agent2+agent3 → back to agent1 → agent1 fixes. Iterate until both clear.
  4. **Final review** = Claude Opus 4.8 — confirms physics + code + tests all check out.
  5. Record verdict here, commit, next constraint.
- **Reviewers:** agent2 & agent3 both codex; final review Claude Opus 4.8.
- **Budget rule:** check `~/bin/codex_usage.sh` between constraints. When WEEKLY
  ≥95% used (≤5% left), finish in-flight constraint, PAUSE, report. Resume after reset.
- **Codex invocation:** `bash ~/bin/codex_worker.sh [--out f] [-C repo] "prompt"`
  (hardened: hard timeout, flock semaphore ≤4, loud failures). Never >4 parallel.
- **Honesty gate:** a constraint that can't be gotten right in a few iterations →
  mark `NEEDS-HUMAN-PHYSICS`, do NOT pass a wrong impl through review.
- **Orchestrator hygiene:** agents return terse verdicts only; orchestrator reads
  --out summaries, not full code. Keep context lean.

## Recovery

- Pre-revert snapshot of the prior (distrusted) run: tag
  `phase3-constraints-snapshot-2026-05-28` (local + origin).
- Clean baseline before rebuild: commit `163728c`.

## Phase status

- [x] Scaffold built (Opus #1) — flavor_catalog_constraints/ (Protocol + auto-discovery + schema-flex typed anchors + adapter boundary)
- [x] Scaffold reviewed (Opus #2) — APPROVE-WITH-CHANGES; isolation held under all break-it probes
- [x] Scaffold fixes applied (real-number contract guard, K001 e2e physics test, registry comment) — 16/16 tests pass
- [x] Scaffold committed

## Systemic conventions (from K001 review — apply to ALL ΔF=2 constraints)

- **QCD running MANDATORY**: use the `*_with_running` evaluator (Wilsons → μ=2 GeV). Non-running underpredicts ~10×.
- **Budgets uncertainty-aware bands, not central residuals**: ε_K band per `docs/audits/epsilon_k_sm_decision.md` (loose edge ~3e-4). Long-distance Δm (kaon, charm): `|M12^NP| ≤ Δm^exp/2`. Bd/Bs: SM-vs-exp room w/ uncertainties.
- Next waves: **pre-stage adapter wrappers in one commit before fanning out** (avoid shared-file write races).

## Per-constraint status

| ID | Observable | agent1 | agent2 | agent3 | opus | committed |
|----|-----------|--------|--------|--------|------|-----------|
| K001 | ε_K | ✅ | ✅ | ✅ | ✅ APPROVE | ✅ 2fdfaa0 |
| K002 | Δm_K | ✅ | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 67953fe |
| B001 | Δm_d | ✅ | PHYSICS-OK | fixed | ✅ APPROVE | ✅ 8060cfb |
| B003 | Δm_s | ✅ | PHYSICS-OK | fixed | ✅ APPROVE | ✅ 6e65846 |
| C001 | D⁰ mix | ✅ | PHYSICS-OK | fixed (blocker) | ✅ APPROVE | ✅ 9d72b6b |
| K004 | BR(K⁺→π⁺νν̄) | ✅ built | PHYSICS-OK | fixed (2 blockers) | ✅ APPROVE | ✅ 4ff15a3 |

**Done: 6 / ~102.** K004 = first machinery-building constraint (new ΔS=1 SM core + adapter; RS-NP proxy flagged in NEEDS_HUMAN_PHYSICS.md). Calibration of the hard path: PASS. Adapter wrappers: fd2f46a. Wave-1 reconciliation PASS (no clobber, suite green).
Per-constraint audit trail (plans, reviews, fixes, verdicts) committed under `.orchestration/runs/<ID>/`.

## Next waves (machinery-building — slower, costlier, pre-stage adapters first)

Remaining ΔF=2-adjacent CP phases (B002 S_ψKs, B005-type) then the machinery families:
rare ΔS=1 kaon (K003-K018), EDMs (E*), charged-lepton LFV (L*), charm rare (C*),
top/Higgs/EW (T*/EW*), collider (CR*), secondary. Raise CODEX_MAX_CONCURRENCY→6; pipeline waves.

## Planned order (validation-first)

K001 ε_K → K002 Δm_K → B(Δm_d) → B(Δm_s) → C(D⁰ mixing)  [existing physics; validate scaffold]
→ then machinery-building constraints, family by family.
