# CLEANUP_REVIEW_R1 — Reviewer Round-1 Critique

**Verdict: APPROVE-WITH-CHANGES.** One more planner round is required before execution.

## Critical issues (must-fix)

### CR-1. C16's v0.2 tag-creation task is bogus
- `flavor-catalog-v0.2` already exists; points to `835cf48` (actual master-compile v0.2 commit), NOT the planner-cited `42ac647` (which is the factcheck consolidation).
- R17-I3 is already satisfied (tag exists). Drop the tag-creation task in C16; reclassify R17-I3 as closed.
- Without this fix, the cleanup agent will hit "tag already exists" and either force-overwrite (destroying provenance) or block.

### CR-2. C01 doesn't account for `tests/test_modern_phenomenology.py` (which exists)
- File is 174+ lines and exercises `build_modern_point_phenomenology_artifact`.
- Doesn't pin specific ε_K numerical ratios (structural fields only) — should not break.
- C01 must (a) acknowledge it exists, (b) include in Check-2 pytest selection explicitly, (c) state it's not expected to require numerical update, (d) also include `tests/test_modern_scan.py` (monkeypatch stub).
- Letting executing agent "search for callers" of a Tier-1 high-stakes function is an unforced error.

### CR-3. C02a is mis-sized for a single Opus session
- Combines a 15-min Python CLI change with 3 SLURM reruns (30-60 min each, queue-dependent) in one atom.
- Either the agent burns hours idle waiting for SLURM, or exits partial.
- Recommended split: **C02a-code** (CLI flag + tests + manifest placeholders) is one atom; **C02a-scans** (SLURM submit + wait + ingest in a separate session) is a different flow. Mirrors planner's own §E.4 contingency but makes it default.
- Combined with G1 below: prefer to defer SLURM reruns entirely.

## Recommended changes (should-fix)

- **RC-1.** R07-I2 / R22-I1: spell out scipy version reference (>=1.10 for `BinomTestResult.proportion_ci`) or provide closed-form for k>0. State vitest snapshot dir (`flavor_catalog/website/src/lib/__tests__/__snapshots__/`).
- **RC-2.** C01 should NOT bundle export-artifact re-run into the bugfix commit. Split into two commits within C01 (C01a bugfix, C01b artifact regen) for clean revert/cherry-pick.
- **RC-3.** (no-op) count off-by-one: planner says 6, actual is 7. Reconcile.
- **RC-4.** C02a's "3 SLURM reruns" vs "2 new + 1 existing" inconsistency — reconcile.
- **RC-5.** C18 must `grep -n` for error patterns rather than absolute line numbers (MERGE_PLAN.md has been edited; line refs may be stale).
- **RC-6.** C20's R01-I2 rationale: state explicitly that `scripts/calibrate_phase0.py` exists but hasn't been re-run since the audit.
- **RC-7.** C13 (R02-I2 spurion seed) should note "no test re-run required because values don't change".

## Planner open questions — reviewer's calls

- **G1 (RUNA rerun scope):** Reject planner's lean (A). Choose hybrid (A'): keep CLI flag + methodology band paragraph using symbolic 1/√|Δε_K| scaling already in `phase2_h5_signoff.md:100-101`; defer actual SLURM reruns to a separately-tracked unit (C02c). Decouples Opus throughput from SLURM queue.
- **G2 (collaborator artifact re-export):** Accept planner's lean (A) BUT split into two commits within the same C01 unit. Re-export, but in separate commit from the bugfix.
- **G3 (factcheck consolidation):** Accept planner's lean (B). Per-master-compile annotation is cheap and explicit; `tools/aggregate_factchecks.py` is nice-to-have but not load-bearing.

## Things the planner missed

- **M-1.** No check whether `quarkConstraints/paper_0710_1869/` also carries duplicate kaon constants. Quick grep needed.
- **M-2.** No mandate that working tree is clean (`git status --porcelain`) before each unit starts.
- **M-3.** No specification of how ISSUES.md `## Closed / Accepted-risk` section is organized (chronological vs. by unit).
- **M-4.** No checkpoint annotated tag after Tier 1 closes (e.g. `cleanup-tier1-complete`). Tier 1 is the only physics-affecting tier; user wants to gate on it.
- **M-5.** `_evaluate_delta_mk_from_bridge` regression directions ("VLL ~−23%, LR ~+15%") need one-line derivation footnote.
- **M-6.** Wall-time compressibility: C18 (MERGE_PLAN.md edits), C15 (pytest_selection backfill), C17 (.gitkeep cleanup), C20 (5 INFO ACCEPTED-RISK closures) can be done by the orchestrator inline via shell scripts. Cuts ~1-2h Opus budget.
- **M-7.** R04-I2 import lift verified safe (no circular with `qcd_running.py`); update C03 wording from "try to lift" to unconditional.
- **M-8.** Final tag name `v2026q2-post-cleanup` is fine; state it sits above `v2026q2-catalog-complete` chronologically. Confirm with user re: signing convention.

## Confidence

One more planner round needed to address CR-1/CR-2/CR-3 and adopt G1/G2/G3 dispositions. Without it, executing agents will hit confused states on C01 (tests), C16 (tag collision), and C02a will partially fail or block on SLURM.
