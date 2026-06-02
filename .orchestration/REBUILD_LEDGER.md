# Flavor-Constraint Rebuild — Orchestration Ledger

Durable state for the from-scratch constraint implementation. Survives
orchestrator context resets. Updated after every per-constraint gate.

## ✅ REBUILD COMPLETE (updated 2026-06-02)

- **95 / 95 PRIMARY + 8 SECONDARY = 103 constraints committed + pushed** (main in sync). Registry: 0 import_failures. Full suite **1054 passed**. By severity: 87 HARD / 12 INFO / 4 SOFT. By family: beauty 26, top_higgs_ew 21, kaon 16, collider_rs 14, charged_lepton 11, charm 8, edm_neutrino 7.
- **Every constraint went through the full gate** (codex agent1 → codex agent2 physics ∥ codex agent3 code → route fixes → Claude Opus 4.8 independent final review → commit). The review gate caught real physics/architecture bugs that green tests missed (e.g. T014 hadronic-vs-total-Z denominator + adapter/core boundary; B013 measurement-vs-theory budget; B007 silent anchor block_key accept; earlier waves: ε_K running, K012 Im-sensitivity, L005 m_μ⁵, CR007 graviton root, …).
- **WAVE-20 (5):** L006 `a6f0d74`, CR010 `75afcfc`, E009 `45e9477`, K019 `6e329c5`, B007 `ca580ca`.
- **WAVE-21 (4):** CR012 `69648e2`, B008 `6ece7d8`, K020 `966af3b`, T014 `0183016`.
- **WAVE-22 (3):** CR013 `cc5e457`, K021 `481369c`, B013 `9a6db0b`.
- **WAVE-23 (2, FINAL):** CR014 `6910282` (four-top top-philic vector, 95th PRIMARY, Opus APPROVE), B014 `67b2f46` (B→ργ/ωγ exclusive b→dγ, last constraint, Opus APPROVE). Full suite **1054 passed**.
- **REMAINING WORK:** none in the build. The standing human-input items (RS new-physics matching proxies, incalculable-SM honest stubs) are catalogued in `.orchestration/NEEDS_HUMAN_PHYSICS.md` — the final deliverable. All proxies are non-fabricated, flagged in each constraint's docstring + diagnostics + commit message.
- Build prompts via `cat runs/<ID>/_hdr.md runs/_agent1_common_v2.md > runs/<ID>/agent1_prompt.md` (reviewers: `_agent2_common.md`/`_agent3_common.md`); codex via `~/bin/codex_worker.sh` (CODEX_MAX_CONCURRENCY=6). NOTE: launch each agent1 as its OWN background Bash task (NOT a chained-heredoc eval — that wedges). See ORCHESTRATOR_RUNBOOK.md.
- Build prompts via `cat runs/<ID>/_hdr.md runs/_agent1_common_v2.md > runs/<ID>/agent1_prompt.md` (reviewers: `_agent2_common.md`/`_agent3_common.md`); codex via `~/bin/codex_worker.sh` (CODEX_MAX_CONCURRENCY=6). See ORCHESTRATOR_RUNBOOK.md for the full pattern.

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

| K006 | K_L→μμ | ✅ built | fixed (provenance) | fixed | ✅ APPROVE | ✅ df0eb34 |
| C004 | D⁰→μμ | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ cf65c4b |
| T010 | Z→bb̄ | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ f8ad10d |
| T001 | t→cZ | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 72a22e4 |
| B016 | B⁺→K⁺μμ | ✅ ext | fixed (3) | CODE-OK | ✅ APPROVE | ✅ 3cafb41 |

| K008 | K_L→π⁰ee | ✅ ext | fixed (blocker) | CODE-OK | ✅ APPROVE | ✅ |
| C005 | D⁰→ee | ✅ reuse | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |
| T012 | Z→cc̄ | ✅ reuse | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |
| T002 | t→uZ | ✅ reuse | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |
| B015 | B→X_sℓℓ | ✅ ext | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |

| K009 | K_L→π⁰μμ | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| C007 | D⁺→π⁺μμ | ✅ built | fixed (q²+FF) | CODE-OK | ✅ APPROVE | ✅ |
| T015 | Z→eμ LFV | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| T003 | t→cγ | ✅ reuse | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |
| B017 | B→K(*)ℓℓ | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |

| K010 | K_S→π⁰ee | ✅ reuse | fixed (a_S sign) | CODE-OK | ✅ APPROVE | ✅ |
| C006 | D⁰→eμ LFV | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| T004 | t→uγ | ✅ reuse | fixed (budget) | CODE-OK | ✅ APPROVE | ✅ |
| T016 | Z→eτ LFV | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| B018 | R_K LFU | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |

| K012 | K_S→μμ | ✅ reuse | fixed (Im proj) | CODE-OK | ✅ APPROVE | ✅ |
| C008 | D⁺→πeμ LFV | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| T005 | t→cg | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| T017 | Z→μτ LFV | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| B019 | R_K* LFU | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |

| B009 | B⁺→τν | ✅ built | fixed (anchors) | CODE-OK | ✅ APPROVE | ✅ |
| K018 | K_l3/V_us | ✅ built | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |
| CR001 | KK-g→tt̄ | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| E001 | electron EDM | ✅ built | PHYSICS-OK | fixed (anchor) | ✅ APPROVE | ✅ |
| L002 | μ→3e | ✅ built | fixed (interf) | fixed (flavor) | ✅ APPROVE | ✅ |

| E002 | muon EDM | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| L009 | τ→3μ | ✅ reuse | PHYSICS-OK | fixed (leak) | ✅ APPROVE | ✅ |
| CR002 | VLQ pair | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| B021 | Λ_b→Λμμ | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| T018 | h→μτ LFV | ✅ built | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |

| B025 | R_D | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| L003 | μ→e conv Al | ✅ built | fixed (SOFT) | CODE-OK | ✅ APPROVE | ✅ |
| T019 | h→eτ LFV | ✅ reuse | fixed (norm) | CODE-OK | ✅ APPROVE | ✅ |
| L010 | τ→3e | ✅ reuse | fixed (budget) | CODE-OK | ✅ APPROVE | ✅ |
| CR003 | VLQ 3rd-gen | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |

| B026 | R_D* | ✅ reuse | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |
| T020 | h→eμ LFV | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| L004 | μ→e conv Au | ✅ reuse | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |
| CR004 | VLQ B→Hb | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| T006 | t→ug | ✅ reuse | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |

| L005 | μ→e conv Ti | ✅ reuse | fixed (m_μ⁵ core) | fixed | ✅ APPROVE | ✅ |
| K017 | R_K LFU | ✅ reuse | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |
| T007 | t→cH | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| CR005 | Z' dilepton | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| C003 | ΔA_CP (STUB) | ✅ stub | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |

| T008 | t→uH | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| CR006 | KK W' CC | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| E004 | neutron EDM (STUB) | ✅ stub | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| K013 | K_L→π⁰γγ (STUB) | ✅ stub | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| B032 | B→πK (STUB) | ✅ stub | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |

| CR007 | KK-graviton | ✅ reuse | fixed (spin-2 root) | CODE-OK | ✅ APPROVE | ✅ |
| E006 | Hg EDM (STUB) | ✅ stub | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| E008 | quark cEDM (STUB) | ✅ stub | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |
| B033 | S_φKs (STUB) | ✅ stub | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| CR011 | VBS long. (STUB) | ✅ stub | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |

| EW001 | S,T,U oblique | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| K003 | Re(ε'/ε) (STUB) | ✅ stub | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| L007 | τ→μγ | ✅ reuse | PHYSICS-OK | fixed | ✅ APPROVE | ✅ |
| CR008 | VLQ T singlet | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| T011 | Z→bb̄ asym | ✅ promoted | fixed (standalone) | fixed | ✅ APPROVE | ✅ |

| B034 | φ_s^sss (STUB) | ✅ stub | fixed (budget) | CODE-OK | ✅ APPROVE | ✅ |
| CR009 | DY contact EFT | ✅ reuse | fixed (floor) | CODE-OK | ✅ APPROVE | ✅ |
| E007 | Ra/Xe EDM (STUB) | ✅ stub | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| L008 | τ→eγ | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| L023 | ν trident | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ |
| L006 | muonium M→M̄ | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ a6f0d74 |
| CR010 | VLQ (T,B) pair | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 75afcfc |
| E009 | Weinberg 3g (STUB) | ✅ stub | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 45e9477 |
| K019 | K_L→eμ LFV (SEC) | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 6e329c5 |
| B007 | B0/Bs→ee (SEC) | ✅ reuse | PHYSICS-OK | fixed (1 blocker) | ✅ APPROVE | ✅ ca580ca |
| CR012 | diboson spin-1 | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 69648e2 |
| B008 | Bs/B0→ττ (SEC) | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 6ece7d8 |
| K020 | K+→π+eμ LFV (SEC) | ✅ built | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 966af3b |
| T014 | FCNC Z→bs/bd/sd (SEC) | ✅ built | fixed (blocker) | fixed (blocker) | ✅ APPROVE | ✅ 0183016 |
| CR013 | diphoton spin-0/2 | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ cc5e457 |
| K021 | K_L→π⁰eμ LFV (SEC) | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 481369c |
| B013 | Bs→φγ excl (SEC) | ✅ reuse | fixed (SHOULD-FIX) | CODE-OK | ✅ APPROVE | ✅ 9a6db0b |
| CR014 | four-top vector | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 6910282 |
| B014 | B→ργ/ωγ b→dγ (SEC) | ✅ reuse | PHYSICS-OK | CODE-OK | ✅ APPROVE | ✅ 67b2f46 |

**Done: 89 / 95 PRIMARY.** ⏸️ PAUSED at weekly 94% (spend-to-95% rule).
REMAINING 6 PRIMARY: CR010, CR012, CR013, CR014 (collider reuses), E009 (Weinberg 3-gluon stub), L006 (muonium-antimuonium, niche).
THEN SECONDARY TIER (~8): B007, B008, B013, B014, K019, K020, K021, T014.
RESUME after weekly reset 2026-05-31 16:36 — read ORCHESTRATOR_RUNBOOK.md; same one-per-shared-module wave pattern. Wave-18. T011 promoted from merged-stub to standalone. REMAINING (~11): B034, CR009/CR010/CR012/CR013/CR014, E007/E009, L006/L008/L023. SHAs in git log. Wave-17. SHAs in git log. REMAINING (~16): reuses CR008-CR010/CR012-CR014, T011?/T014, secondaries B007/B008/B013/B014, K019/K020/K021, L006/L023, E007/E009 stubs, B034 stub. Wave-16: 2 reuses + 3 HARD honest-stubs (E004/K013/B032 — INFO/non-vetoing, dual NEEDS-HUMAN, no faked hadronic/ChPT). SHAs in git log. Wave-15. NOTE: L005 review found a shared mu_e_conversion m_μ⁵ normalization bug → fixed core + re-verified committed L003/L004 (Opus). C003 = first HARD honest-stub (INFO, dual NEEDS-HUMAN). SHAs in git log. Wave-14 (all reuses). SHAs in git log. Wave-13. New modules: semileptonic_lfu, mu_e_conversion. SHAs in git log. Wave-12. New modules: rare_b_baryon_dilepton, higgs_lfv. SHAs in git log. Wave-11 = NEW-module pioneers (leptonic_tree, ckm_extraction, collider_resonance, edm, lfv_three_body) → unlock B025/B026?, T018-20 via higgs_lfv next, CR fan-out, E002, L009/L010/L003-005. SHAs in git log. Wave-10 fan-out. SHAs in git log. New module: rare_b_kstar_dilepton. Wave-9 fan-out. SHAs in git log. (Target is ~95 catalog entries, not 102 — some IDs merged/removed.) Wave-8 fan-out. SHAs in git log. New modules: rare_charm_semileptonic, zpole_lfv. Wave-7 = fan-out of wave-6 modules (one per shared module). Review caught K008 direct-CP structure bug (y7V/y7A) + C005 stale-golden. SHAs in git log. Wave-6 = multi-family PIONEERS — 5 new shared modules built in parallel:
`rare_kaon_dilepton`(→K008/K009/K010/K012), `rare_charm_dilepton`(→C005/C007), `zpole`(→T011/T012/T015-T017),
`top_fcnc`(→T002-T008), exclusive `rare_b_dilepton`(→B015/B017-B019/B021). Ready for wide fan-out next.

**[prior wave-5 note below]** Wave-5 = first fan-out reusing wave-4 modules (B006/B012/B023), one-per-shared-module-per-wave to avoid contention. NP proxies inherit from the shared modules (already in NEEDS_HUMAN_PHYSICS.md).
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
