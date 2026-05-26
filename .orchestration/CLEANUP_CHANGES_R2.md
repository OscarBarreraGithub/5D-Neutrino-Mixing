# CLEANUP_PLAN.md — Round-1 → Round-2 changelog

Author: Claude Opus 4.7 (planner-r2)
Date: 2026-05-25
Inputs: `.orchestration/CLEANUP_PLAN.md` (r1), `.orchestration/CLEANUP_REVIEW_R1.md`.

This changelog enumerates every substantive change from planner-r1 to planner-r2.
Pure copy-edits (rephrasing for clarity) are NOT logged unless they change semantics.

## Critical-issue fixes

### CR-1 — Drop C16's v0.2 tag-creation task
- **r1 §C C16**: planned to create annotated tag `flavor-catalog-v0.2` on
  `42ac647`.
- **r2 §C C16**: task dropped. Pre-flight verification by planner-r2 confirmed
  via `git tag -l 'flavor-catalog-v0*'` that `flavor-catalog-v0.2` ALREADY EXISTS
  at `835cf48` (not the r1-cited `42ac647`, which was wrong).
- **r2 §A triage**: R17-I3 row's Unit column changed from `C16` to `(closed)`
  with annotation "TAG ALREADY EXISTS — verified planner-r2".
- **r2 §C C21**: C21 now records R17-I3 closure with rationale
  "PRE-EXISTING-TAG-VERIFIED at 835cf48".
- **r2 §G G5**: disposition updated from "lean A — tag retroactively" to
  "dropped per CR-1".
- **r2 §H final state**: tag list updated.
- **r2 Appendix**: C16 issue count drops from 3 → 2.

### CR-2 — Bake test surface into C01
- **r1 §C C01**: Check-2 said
  "`pytest tests/test_quark_deltaf2.py tests/test_modern_phenomenology.py` (if
  it exists)". The conditional and the "search for callers" wording were both
  problematic.
- **r2 §C C01**: explicit selection
  `pytest tests/test_quark_deltaf2.py tests/test_modern_phenomenology.py
  tests/test_modern_scan.py` (mandatory, no conditional). Planner-r2 verified
  both files exist (249 + 474 lines).
- **r2 §D D1**: added subsection
  "Tests already on disk that MUST stay green (CR-2)" explicitly listing both
  files and noting neither is expected to require numerical update (structural
  assertions only).
- **r2 §E.2.2**: updated to cite C01's explicit triple as the canonical example.

### CR-3 — Split C02a
- **r1 §C C02a**: combined CLI flag + 3 SLURM reruns + ingest + paper-text in
  one Opus atom — too large.
- **r2 §C**: split into:
  - **C02a-code** (Opus, in-scope): CLI flag + tests + methodology paragraph
    using **symbolic** `1/sqrt(|Δε_K|)` scaling + manifest placeholders.
  - **C02c** (deferred): SLURM RUNA reruns + ingest + paragraph refinement
    with measured values.
  - C02b retained as deferred (paper finalization).
- **r2 §B Tier 1**: serial estimate reduced (~3-5h → ~3-4h Opus interactive;
  SLURM wall-time moved to C02c).
- **r2 §B Tier 6 (new)**: "DEFERRED-BUT-TRACKED" tier introduced containing
  C02b and C02c.
- **r2 §F.1**: CLEANUP_QUEUE.md rows for C02b and C02c added with status
  `deferred`.
- **r2 §F.2**: `cleanup_progress.json` schema gains entries for C02b/C02c with
  `lane: "deferred"`.
- **r2 §G G1**: disposition updated to "hybrid (A')".
- **r2 Appendix**: C02 row split into three rows; C02a-code = partial,
  C02b/C02c = 0 in cleanup.

## Recommended-change fixes

### RC-1 — scipy floor + vitest snapshot dir
- **r1 §C C04**: ambiguous re: scipy version + snapshot dir.
- **r2 §C C04**: explicit requirement "scipy>=1.10 for
  `BinomTestResult.proportion_ci`, OR a closed-form Wilson UL fallback embedded
  in the test" (closed-form provided inline). Vitest snapshot dir stated as
  `flavor_catalog/website/src/lib/__tests__/__snapshots__/`.

### RC-2 — Split C01 into 2 commits
- **r1 §E.4**: noted C02a may need 2-commit split; nothing on C01.
- **r2 §C C01**: explicit two-commit split — **C01a** (bugfix + tests) and
  **C01b** (collaborator artifact regen). Both inside the same Opus session.
- **r2 §E.4**: C01 listed as the only 2-commit unit (no longer C02a).
- **r2 §F.1 / §F.2**: tracking files updated.

### RC-3 — (no-op) off-by-one
- **r1**: claimed 6 no-ops; the actual count was 7 (R08-I1, R14-I2, R15-I3,
  R16-I1, R16-I2, R18-I2, R18-I4).
- **r2 §A bookkeeping checkpoint**: count corrected to 7.
- **r2 §G G4**: count corrected.
- **r2 Appendix no-op row**: count 7.

### RC-4 — "3 SLURM reruns" vs "2 new + 1 existing"
- **r1**: planner-r1 was inconsistent — some sentences said "3", some said "3 at
  three edges" implying central was new vs existing.
- **r2 §D D3**: clarified — all three budget edges are NEW SLURM scans
  (central, low, high). Existing scan was at a single budget value with no
  `--epsilon-k-budget` flag, and is NOT re-used as one of the three. C02c
  submits three NEW jobs.

### RC-5 — C18 grep for error patterns, not absolute line numbers
- **r1 §C C18**: implicitly relied on ISSUES.md cite-text line numbers.
- **r2 §C C18**: re-classified as orchestrator-deterministic (per M-6). Action
  is `grep -n` for error pattern + targeted-replace; absolute line numbers from
  ISSUES.md cite-text are **never** used.

### RC-6 — C20 R01-I2 rationale needs calibrate_phase0.py footnote
- **r1 §C C20 (R01-I2)**: rationale "calibration deferred".
- **r2 §C C20 (R01-I2)**: rationale extended — explicitly states
  `scripts/calibrate_phase0.py` exists in the repo but has not been re-run since
  the post-audit kaon-constants update. Planner-r2 verified file exists.

### RC-7 — C13 note "no test re-run required"
- **r1 §C C13**: Check-2 was "pytest -k benchmarks passes".
- **r2 §C C13**: Check-2 reads "no test re-run required because no values
  change; pure-docstring change; `python -c "import quarkConstraints.benchmarks"`
  imports cleanly".

## Open-question dispositions

- **G1**: r1 lean A → r2 hybrid (A') — see CR-3.
- **G2**: r1 lean A → r2 lean A with two-commit split — see RC-2.
- **G3**: r1 lean B retained → r2 lean B retained.
- **G4**: r1 count 6 → r2 count 7 — see RC-3.
- **G5**: r1 lean A → r2 dropped — see CR-1.
- **G6**: r1 lean retained (drop only redundant).
- **G7**: r1 lean B retained (paper-archive companion).

## Missed items now incorporated

### M-1 — Pre-flight grep for kaon constants in paper_0710_1869
- **r2 §C C01**: added "Pre-flight (M-1)" subsection with the exact grep command
  and the planner-r2 verification result (no duplicate hadronic literals; only
  `default_kaon` artifact IDs).

### M-2 — git status --porcelain clean-tree check
- **r2 §C "Common preflight (M-2)"**: every unit must verify
  `git status --porcelain` is empty before edits. BLOCK on non-empty.

### M-3 — ISSUES.md Closed section organized by C##
- **r2 §E.5**: updated protocol — issues organized into `### Closed by C##`
  subsections rather than chronologically. Within each unit subsection,
  original R##-I# order preserved.
- **r2 §F.3 / §H**: reflects M-3.

### M-4 — Checkpoint tag cleanup-tier1-complete
- **r2 §B Tier 1**: introduced annotated tag created post-C03; pushed to
  origin.
- **r2 §E.6**: new subsection "Tier-1 checkpoint (M-4)" with exact commands.
- **r2 §H final state**: tag listed.

### M-5 — One-line derivation for VLL ~−23%, LR ~+15%
- **r2 §D D1 "M-5 one-line derivation footnote"**: explicit derivation provided
  using post-audit B-factor ratios. To be included in the C01a commit body.

### M-6 — Reclassify C15, C17, C18, C20 as orchestrator-deterministic
- **r2 §B Tier 5**: four units annotated as "orchestrator-deterministic"; Opus
  reviews commit only.
- **r2 §C**: each of C15, C17, C18, C20 has "Execution lane:
  orchestrator-deterministic" plus an explicit shell action.
- **r2 §E.1**: "Two lanes" section added.
- **r2 §F.1 / §F.2**: tracking files include `lane` column / field.
- **r2 §B grand total**: ~1-2 hours of Opus budget moved to shell scripts.

### M-7 — R04-I2 import lift unconditional
- **r1 §C C03**: wording was "try to lift; retain function-local if circular".
- **r2 §C C03**: wording is "**lift unconditionally**, verified non-circular
  with `qcd_running.py` by reviewer-r1 M-7".

### M-8 — v2026q2-post-cleanup tag chronology
- **r2 §C C21**: confirmed `v2026q2-catalog-complete` already exists
  (`git tag -l 'v2026q2*'` returned only it). `v2026q2-post-cleanup` will
  chronologically follow it (descendant in commit DAG). User asked to confirm
  tag-signing convention (annotated, unsigned — same as existing).
- **r2 §H tag list**: chronology stated.

## Other r1 → r2 housekeeping

- **r2 §A**: Unit column updated for R03-I3 (now `C02a-code (+ C02c deferred)`)
  and R17-I3 (now `(closed)`).
- **r2 §B Tier 1**: drops R17-I3 from "Issues" enumeration.
- **r2 §B Tier 5**: drops R17-I3 from C16; C16 now resolves R07-I4 + R14-I3 only.
- **r2 §B Tier 6 (new)**: separate section for deferred-but-tracked units.
- **r2 §F.1**: CLEANUP_QUEUE.md gains a `Lane` column (opus | shell | deferred).
- **r2 §F.2**: cleanup_progress.json gains `lane`, `tier1_checkpoint_tag`,
  `final_tag`, and `dropped_tasks` fields.
- **r2 §H**: end-state §7 (tags) now distinguishes mid-cleanup tag from final
  tag.

## Bookkeeping invariants — both rounds satisfied

- Total issues accounted for: 78 in both rounds.
- Closed-by-cleanup-unit count: 70 in r1, 70 in r2 (R17-I3 moves from "C16
  closure" to "PRE-EXISTING-TAG-VERIFIED no-op closure"; net delta on
  ISSUES.md is zero).
- No-op count: r1 said 6 (off-by-one); r2 says 7 (correct).
