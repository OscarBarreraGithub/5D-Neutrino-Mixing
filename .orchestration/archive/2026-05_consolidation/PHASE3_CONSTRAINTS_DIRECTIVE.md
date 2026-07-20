# Phase 3 — Implement all ~100 flavor-catalog constraints

User directive received 2026-05-28. **MUST survive compactification.**

## Goal

Implement every constraint in `flavor_catalog/processes/` (94 PRIMARY + 8 SECONDARY = ~102 entries). Each constraint must be:
- **Isolated** to its own file or class
- **Atomic** — modifying or deleting one does not affect others
- **Plug-and-play** — easy to add or remove

## Execution pattern

### Step 1 — Scaffolding (BEFORE any constraint implementation)
1. Spawn Opus agent to design the constraint scaffolding (directory structure, base class, plugin interface, how each constraint plugs in)
2. Spawn second Opus to review the scaffolding plan
3. Iterate until both agree
4. Set up the directory structure + template files

### Step 2 — Per-constraint implementation (SERIAL across constraints)

For EACH constraint (one at a time, NOT parallel across constraints):

- **agent1** (codex gpt-5.5 xhigh): plan + implement the constraint
- **agent2** (codex gpt-5.5 xhigh, parallel with agent3 if safe): physics fact-check against source literature
- **agent3** (codex gpt-5.5 xhigh, parallel with agent2 if safe): code review + numerical verification + unit tests
- agent2 and agent3 deliver findings to agent1
- agent1 fixes the code based on the two reviews
- **opus reviewer**: final sign-off that everything makes sense
- Then move to next constraint

### Codex orchestration caveat (LEARNED THE HARD WAY)
- Codex CLI processes can hang silently with 0 CPU time, never notifying completion
- Foreground sequential execution works fine (~1-3 min per task)
- Parallel >2 codex workers hits file-lock issue ("No locks available, os error 37") on RHEL 8 kernel
- **agent2 + agent3 in parallel within a single constraint is OK** (only 2 simultaneous)
- **MUST actively monitor codex sessions** — task notifications do NOT fire for stuck 0-CPU procs
- Check process CPU time periodically; if stuck >5min with 0 CPU, kill and retry

## Output requirements

- Each constraint gets a self-contained file or class
- Suspicious/problematic constraints flagged in a dedicated document
- Do NOT stop until all constraints are in
- Offload tokens to agents — orchestrator stays lean

## Order
1. Finish Batch B (B3 calibration in progress)
2. Opus pass writing consolidation report appendix
3. Final Codex verification round
4. Push everything to GitHub
5. **THEN** begin Phase 3: scaffolding → per-constraint serial implementation

## Constraint inventory location
- `flavor_catalog/processes/{beauty,charged_lepton,charm,collider_rs,edm_neutrino,kaon,top_higgs_ew}/` — PRIMARY entries
- `flavor_catalog/processes/secondary/{...}/` — SECONDARY entries (8)
- Each yaml/tex pair holds the physics anchor, the source citations, and the constraint spec
