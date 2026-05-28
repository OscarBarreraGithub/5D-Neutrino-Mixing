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

## Per-constraint status

| ID | Observable | agent1 | agent2(phys) | agent3(code) | opus | committed | notes |
|----|-----------|--------|--------------|--------------|------|-----------|-------|
| (none yet) | | | | | | | |

## Planned order (validation-first)

K001 ε_K → K002 Δm_K → B(Δm_d) → B(Δm_s) → C(D⁰ mixing)  [existing physics; validate scaffold]
→ then machinery-building constraints, family by family.
