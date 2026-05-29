# ORCHESTRATOR RESUME RUNBOOK — read this FIRST on resume/compaction

You (Claude, the orchestrator) are rebuilding ~95 flavor-physics constraints from scratch
under a strict review pattern. If your context was compacted, do NOT improvise — follow this
file + `REBUILD_LEDGER.md` exactly. Do not invent prompts; reuse the templates below.

## 0. Re-orient (do these reads, nothing else first)
1. `.orchestration/REBUILD_LEDGER.md` — per-constraint status table (what's done, SHAs), conventions.
2. `git log --oneline -15` — confirm last committed constraint matches the ledger.
3. `bash ~/bin/codex_usage.sh` — check WEEKLY %. If ≥95%, PAUSE and tell the user.
4. `python -m pytest tests/constraints/ -q` — confirm green baseline before launching anything.

## 1. The per-constraint pipeline (NON-NEGOTIABLE — user's pattern)
For EACH constraint: agent1 (codex, implement) → agent2 (codex, physics) ∥ agent3 (codex, code)
review IN PARALLEL → route issues to agent1-fix (codex) → Claude **Opus 4.8** final review → commit.
- agent1/agent2/agent3 = codex gpt-5.5 xhigh via the wrapper. Final review = a fresh Opus 4.8 Agent.
- Opus reviewers must INDEPENDENTLY verify (numbers reproduced), not trust the codex verdict.

## 2. Codex invocation (hardened wrapper — never call codex directly)
`CODEX_TIMEOUT=<s> CODEX_MAX_CONCURRENCY=6 bash ~/bin/codex_worker.sh --out <runs/ID/x_out.md> -C <repo> "$(cat <prompt>)"`
- Run as background Bash (`run_in_background: true`); never >6 concurrent (flock limit). Concurrency held at 6 (do NOT raise — error risk, per user).
- Build prompts as: `cat .orchestration/runs/<ID>/_hdr.md .orchestration/runs/_agent1_common_v2.md > .../agent1_prompt.md`
  (reviewers: `_agent2_common.md` / `_agent3_common.md`). Write a per-ID `_hdr.md` with the observable + which module to build/reuse.

## 3. Wave methodology (speed without error risk)
- **Pre-stage shared modules**: build a family's physics module ONCE (a "pioneer"), commit, then fan out
  consumers as separate files. RULE: at most ONE constraint touching a given shared module per wave
  (else concurrent writes clobber). Different modules → safe parallel (5-ish per wave).
- After every wave: run full `tests/constraints/ -q` BEFORE trusting any per-agent "green" (concurrent
  edits can clobber); reconcile if a wrapper is missing.

## 4. Mandatory conventions (from reviews — apply always)
- QCD **running** mandatory for ΔF=2 / dipole (use `*_with_running`, μ=2 GeV).
- Budgets are **uncertainty-aware bands**, never bare central residuals; for proxy-NP constraints
  inflate the budget with a documented proxy-theory uncertainty (~30% used for b→sℓℓ).
- Numeric `ConstraintResult` fields are **real floats** (complex → diagnostics). Anchors via scaffold
  `load_anchor` (loud fail), never `load_pdg_block`/hardcoded. Tests cross-check **independently of the
  adapter** (recompute from the core).
- **Honesty gate**: if rigorous RS-NP matching needs EW KK/Z/Z′/lepton/ν couplings (NOT on
  ParameterPoint), build the rigorous SM side + a documented proxy NP, flag `NEEDS-HUMAN-PHYSICS` in
  docstring + diagnostics, and append a bullet to `.orchestration/NEEDS_HUMAN_PHYSICS.md`. Never fake numbers.
- Reviewer false-alarm: a dirty shared worktree / non-empty `git diff quarkConstraints/` during a parallel
  wave is EXPECTED, not an isolation violation (the templates already say this).

## 5. Commit + bookkeeping per wave
- ONE commit per constraint: stage its `<ID>.py` + test + its module/adapter additions + `.orchestration/runs/<ID>/`.
  Message: observable, evaluator, budget, review findings, fixes, Opus verdict.
- Update `REBUILD_LEDGER.md` table (+SHA) and `NEEDS_HUMAN_PHYSICS.md`; commit those; `git push origin main`.

## 6. What's left (update as you go; cross-check the ledger table)
Built shared modules (reuse, don't rebuild): deltaf2 (ΔF=2), rare_kaon_snd (K→πνν), rare_kaon_dilepton
(s→dℓℓ), rare_charm_dilepton (c→uℓℓ), rare_charm_semileptonic (D→πℓℓ), rare_b_dilepton (b→sℓℓ excl+incl),
rare_b_nunu (b→sνν), bsgamma (C7), zpole (Z→ff), top_fcnc (t→qV/qH/qγ), ckm_unitarity, semileptonic_ckm, lepton (μ→eγ).
Still NEEDS NEW modules: leptonic_tree (B009), bsgamma exclusive done, eps_prime (K003), rare_kaon radiative
(K013), lfv_three_body (L002/L009/L010), mu_e_conversion (L003/L004/L005), edm (E*), higgs_lfv (T018/T019/T020),
collider_resonance + vlq_pair + contact (CR*), ckm_extraction (K018), semileptonic_lfu (B025/B026),
charm direct-CPV (C003), nonleptonic (B032/B033/B034). HARD/likely-NEEDS-HUMAN-on-SM-side too:
K003 ε′/ε, K013, B032-B034, C003, E004/E006-E009, CR011. Ship those as documented stubs (dual NEEDS-HUMAN).
Remaining fan-out consumers of built modules: K010,K012 (rare_kaon_dilepton); C006,C008 (rare_charm; LFV);
B007,B008,B013,B014,B018,B019,B021 (rare_b_dilepton); T004,T005,T006,T007,T008 (top_fcnc); T011?,T014,T016,T017 (zpole);
secondary K019,K020,K021. (Verify each ID against `flavor_catalog/processes/<family>/` + the ledger before starting.)

## 7. Recovery anchors
- Pre-rebuild snapshot of the DISTRUSTED prior run: tag `phase3-constraints-snapshot-2026-05-28`. Do not reuse its code.
- Clean baseline: commit `163728c`. Everything since is the trusted rebuild.
