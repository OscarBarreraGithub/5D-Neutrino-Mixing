# PLAN_CHANGES_R2 — Planner round-2 changelog

Author: Claude Opus 4.7 (planner role, round 2)
Date: 2026-05-25
Inputs: `MERGE_PLAN.md` (round 1) + `PLAN_REVIEW_R1.md` (reviewer critique).
Output: `MERGE_PLAN.md` rewritten in place; this file is the per-change audit log.

---

## Critical fixes (reviewer CR-1 / CR-2 / CR-3)

1. **CR-1 fixed**: Removed `82a96f0` from R19's commit list in §B.1. That SHA's subject is `docs(paper): document hadronic input provenance` (verified via `git log -1 --format='%h %s' 82a96f0`) and it now appears in R03 only. R19 count drops from 34 to 33 commits. Added an explicit "intentionally NOT listed here" annotation in the R19 row and a uniqueness invariant in §B.1 with a one-line `awk` enforcement.

2. **CR-2 fixed**: Replaced all references to the non-existent SHA `82b35c0` with the verified canonical SHA `82daa9b` (`flavor-catalog(wave9): CA-v2 ew_tail cycle-2 verdict`). Confirmed `82b35c0` does not exist (`git log -1 82b35c0` returns exit 128). Confirmed `82daa9b` exists and shares its subject with the backup branch tip, justifying the tag-and-delete plan. Added a verification snippet in §A.1 and Appendix A.

3. **CR-3 fixed**: Promoted untracked-file disposition to a hard pre-merge user gate as new **Phase 0** (before any safety tagging or mutation). Added an 8-row decision matrix (`commit` / `gitignore-add` / `rm`) with planner-recommended defaults. Required output: `.orchestration/PRE_MERGE_DECISIONS.md` must exist with all 8 rows resolved before Phase 0.5 runs. Documented the `.gitignore` ordering trap (`!artifacts/**` re-include) and the required explicit deny-line position.

## High-priority recommended fixes

4. **RC-1 (split R10) done**: Split Wave-1 unit into three family-coherent sub-units after reading actual commit subjects via `git log -1 --format='%h %s'`:
   - **R10a** — kaon family (K001-K006, K013): 11 commits, 7 processes.
   - **R10b** — beauty family (B002, B005, B009, B011, B015): 9 commits, 5 processes.
   - **R10c** — top/EW + EDM + charm + lepton (T001, T002, T010, E001, C001, L001): 10 commits, 6 processes across 4 families.
   Total unit count: 22 -> 24. Wall-time estimate: 16.5 h -> 18 h serial (unchanged guidance: multi-session).

5. **RC-2 (literature manifest) done**: §C.1 input #5 now mandates an orchestrator-generated `.orchestration/REFERENCES.md` with one section per unit listing exact arXiv IDs, eq/section numbers, and local PDF paths under `flavor_catalog/external_research/`. Generation moved into §E step 5 (before first reviewer dispatch).

6. **RC-3 (per-unit pytest selection) done**: Added `.orchestration/pytest_selection/<unit>.txt` to the directory layout (§D). Generation step is §E step 6 (one `git log --name-only RANGE -- tests/` pass per unit, plus a default suite for high-risk units).

7. **RC-4 (verify safety-tag push) done**: §A.4 Phase 2 now includes three required `git ls-remote --tags origin | grep ...` verification commands with explicit ABORT semantics. `progress.json` schema (§D.2) gains `safety_tags_pushed` field; field flips to `true` only after verification. The user's request that the safety tag be pushed to origin BEFORE branch deletion is enforced by phase ordering (Phase 2 < Phase 6).

8. **RC-5 (worktree prune in Phase 0) done**: Added new Phase 0.5 with `git worktree prune` + `git worktree list` verification. Runs before safety tagging (so the snapshot reflects the post-prune state) and well before Phase 6 branch deletion (so `git branch -D wa_w5b_charm_edm_work` will not be blocked by the prunable worktree at `/tmp/wa_w5b_charm_edm_worktree`). Row 17 of §A.1 subsumption matrix updated.

9. **RC-6 (per-unit grep checklists) done**: §C.2.3 now contains explicit grep targets for high-risk units (R04, R06, R07, R19) and a global rule for the frozen MEG II constant (per OQ-7). Targets pre-seeded based on the reviewer's suggestions and verified against the file paths in §B.1.

## Recommended fixes / "things the planner missed"

10. **M-1 (PRE_MERGE_STATE.md) done**: Phase 0.5 writes the snapshot containing `git status`, `git branch -avv`, `git tag -l`, `git worktree list`, and `git log -1 HEAD`.

11. **M-2 (GitHub Issues mirroring) done**: New Phase 9 in §A.4 and step 38 in §E mirror each non-ACCEPTED-RISK non-INFO issue to GitHub via `gh issue create`. Issue URLs appended back into local `ISSUES.md` as `Mirrored: <url>`.

12. **M-3 (orchestrator commit policy) done**: New §D.4 lists the exact commits the orchestrator will create, with branch and authorship convention (`Co-Authored-By: Claude (orchestrator)`). Reviewer agents create zero commits.

13. **M-4 (rollback runbook) done**: New §A.6 with the three canonical revert commands. Phase 0.5 writes a richer copy to `.orchestration/ROLLBACK.md`.

14. **M-5 (test-timeout policy) done**: §C.2.2 explicitly says >10 min pytest = PARTIAL, not FAIL. Report template (§C.3) gains a `Wall time` field under Check 2.

15. **M-6 (reviewer-vs-internal-loop mandate) done**: §C.1 adds "Reviewer mandate clarification" paragraph telling reviewers their job is to cross-check external literature and cross-cutting consistency, NOT to re-run the project's Writer-Checker-Opus loop.

## Minor nits

16. **N-1 done**: §A.4 Phase 6 now starts with `git checkout main` before the local delete loop.
17. **N-2 done**: §A.4 Phase 6 remote-delete loop wraps each `git push origin --delete` in a `|| echo WARN ... continuing` so one rejection (e.g., branch protection) does not abort the loop.
18. **N-3 done**: §A.4 Phase 7 now uses a single `git push --tags origin` (the redundant individual push removed).
19. **N-4 noted**: R19 commit count is 33 after CR-1 fix (was reviewer-counted 34). Noted inline in R19 row.
20. **N-5 done**: `progress.json` schema (§D.2) gains `safety_tags_pushed`, `pre_merge_decisions_recorded`, `worktree_pruned`, `github_issues_mirrored` fields; schema version bumped to 2.
21. **N-6 done**: §C.1 input #9 says reviewer MAY use `gh issue list` / `gh pr list` to cross-reference open GitHub issues.

## Open-question answers — incorporated

22. **OQ-1 incorporated**: Phase 0 decision matrix uses the reviewer's recommended defaults (commit the 6 collaborator files + 2 export scripts; do NOT commit `sample_points_may18.*`).
23. **OQ-2 incorporated**: §A.1 rows 4/5 cite `82daa9b` as the canonical replacement.
24. **OQ-3 incorporated**: R10 split into R10a/R10b/R10c (see #4 above).
25. **OQ-4 incorporated**: §C.2.1 spot-check rule tightened to "at least one PRIMARY-tier process from EACH family touched in the wave". Reviewers must document sample selection in the report (new `Sample selection` field in template).
26. **OQ-5 incorporated**: R21+R22 stay separate; R22 explicitly gains `npm run build && grep -r '\\\\(' dist/ | head` check in §C.2.2.
27. **OQ-6 incorporated**: CLAUDE.md update is Phase 10 / §E step 39 (post-merge). Subject confirmed: replace lines 12 and 15.
28. **OQ-7 incorporated**: §C.2.3 MEG II grep rule covers code AND docs; HIGH-severity flag, no code changes proposed.

## Nothing declined

All reviewer recommendations and critical fixes were incorporated. The only judgment call: the planner did NOT add parallelization defaults for R10a/R10b/R10c (they are independent and could batch) — kept the default-serial guidance per the user's 5-hour budget. The parallel-batch option is documented in §E as a noted-but-not-recommended optimization.

## Verification

- `82daa9b` exists and has the expected subject (`flavor-catalog(wave9): CA-v2 ew_tail cycle-2 verdict`).
- `82a96f0` exists and has subject `docs(paper): document hadronic input provenance` (assigned to R03 only).
- `82b35c0` does not exist (`git log -1` exit 128) — all references removed.
- `git worktree list` confirms the prunable `/tmp/wa_w5b_charm_edm_worktree` entry exists; Phase 0.5 prune is correct.
- R10 split SHAs verified against `git log -1 --format='%h %s'` for all 30 Wave-1 commits; family attribution matches the commit-subject prefix in each row.
- `git status` confirms the 8 untracked files match the §A.4 Phase 0 decision matrix.

Plan is now believed execution-ready.
