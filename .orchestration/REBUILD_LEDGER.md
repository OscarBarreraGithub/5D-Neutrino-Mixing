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
| B002 | S_ψKs/sin2β | ✅ | PHYSICS-OK | fixed | ✅ APPROVE | ✅ 80d4224 |
| C002 | charm CPV | ✅ | fixed (budget) | fixed | ✅ APPROVE | ✅ b85f482 |
| K005 | K_L→π⁰νν̄ | ✅ built | PHYSICS-OK | fixed | ✅ APPROVE | ✅ ed736c3 |
| B004 | φ_s (Bs) | ✅ | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 6034c96 |
| EW002 | CKM 1st-row | ✅ built | PHYSICS-OK | fixed | ✅ APPROVE | ✅ d640d2e |
| EW003 | \|Vcb\|/\|Vub\| | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 4947894 |
| L001 | μ→eγ | ✅ | fixed (semantics) | CODE-OK | ✅ APPROVE | ✅ 99063a2 |
| B005 | B_s→μμ | ✅ built | fixed (A_ΔΓ) | CODE-OK | ✅ APPROVE | ✅ be82e38 |
| B011 | B→X_sγ | ✅ built | fixed (2 blockers) | fixed | ✅ APPROVE | ✅ f746d56 |
| B022 | B⁺→K⁺νν̄ | ✅ built | fixed (LD) | fixed | ✅ APPROVE | ✅ e4c64d2 |
| B006 | B⁰→μμ | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 30420b6 |
| B012 | B→K*γ | ✅ reuse | PHYSICS-OK | fixed | ✅ APPROVE | ✅ 6731ac0 |
| B023 | B→K*νν̄ | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 93fd563 |

**Done: 19 / ~102.** Wave-5 = first fan-out reusing wave-4 modules (B006/B012/B023), one-per-shared-module-per-wave to avoid contention. NP proxies inherit from the shared modules (already in NEEDS_HUMAN_PHYSICS.md).
Remaining rare_b_dilepton consumers (B015/B016/B017/B018/B019/B021 + secondary B007/B008/B013/B014): pre-extend the module once, then fan out as pure consumers. Wave-4 = B-rare PIONEERS that built 3 shared modules → now fan out:
- `rare_b_dilepton` (B005) → B006, B015, B016, B017, B021 (+B018/B019 ratios, secondary B007/B008/B013/B014)
- `bsgamma` (B011) → B012
- `rare_b_nunu` (B022) → B023
Review caught real physics bugs (B011 chirality + RG running, B005 A_ΔΓ, B022 long-distance) — pipeline working. Wave-3 = reuse/data-driven (B004 Bs phase, EW002/EW003 CKM data, L001 μ→eγ dipole). New adapters: ckm_unitarity, semileptonic_ckm, lepton. EW002/EW003 are SOFT (SM-vs-data tension). K004/K005 = machinery-building (RS-NP proxy flagged). B002/C002 = CP-phase reuse of ΔF=2 M₁₂ (SM-phase input flagged). All NEEDS-HUMAN items tracked in NEEDS_HUMAN_PHYSICS.md.
Reviewer template (_agent3_common.md) hardened: isolation judged by what a constraint ADDS (shared worktree is dirty by design in parallel waves); cross-check must be independent of the adapter under test. Adapter wrappers: fd2f46a. Wave-1 reconciliation PASS (no clobber, suite green).
Per-constraint audit trail (plans, reviews, fixes, verdicts) committed under `.orchestration/runs/<ID>/`.

## Next waves (machinery-building — slower, costlier, pre-stage adapters first)

Remaining ΔF=2-adjacent CP phases (B002 S_ψKs, B005-type) then the machinery families:
rare ΔS=1 kaon (K003-K018), EDMs (E*), charged-lepton LFV (L*), charm rare (C*),
top/Higgs/EW (T*/EW*), collider (CR*), secondary. Raise CODEX_MAX_CONCURRENCY→6; pipeline waves.

## Planned order (validation-first)

K001 ε_K → K002 Δm_K → B(Δm_d) → B(Δm_s) → C(D⁰ mixing)  [existing physics; validate scaffold]
→ then machinery-building constraints, family by family.
