# PLAN REVIEW — Round 1 (Reviewer)

Reviewer: Claude Opus 4.7 (reviewer role)
Date: 2026-05-25
Plan reviewed: `.orchestration/MERGE_PLAN.md` (planner round 1)

---

## Verdict

**APPROVE-WITH-CHANGES.**

The merge-consolidation half of the plan (Section A) is sound and I independently verified its critical claims. The atomic-review half (Section B) is mostly well-shaped but has **two outright bugs** and one **questionable seam** that must be fixed before R-units are dispatched. None of the issues require re-architecting the plan — they are localised corrections plus three small policy clarifications. With the changes in §C below applied, the plan is execution-ready.

---

## Critical issues (must-fix before any execution)

### CR-1. Commit `82a96f0` is double-assigned to R03 AND R19

In §B.1 the same SHA appears in two units:
- R03 row lists `82a96f0` as one of the bag-parameter audit commits.
- R19 row lists `82a96f0` again among the collider_rs Wave-9 commits.

I checked with `git log -1 --format='%h %s' 82a96f0`:
> `82a96f0 docs(paper): document hadronic input provenance`

That subject clearly belongs in R03 (hadronic provenance is bag-parameter audit, not collider_rs). The R19 listing is a copy-paste error.

**Action**: drop `82a96f0` from R19. While at it, scan all 22 rows for similar duplicates with `awk -F'`' '{for(i=2;i<=NF;i+=2) print $i, NR}' MERGE_PLAN.md | sort | uniq -c | sort -rn | head` (or equivalent) — when COMMIT_INDEX.csv is generated, enforce uniqueness as a precondition.

### CR-2. Planner references a non-existent SHA `82b35c0`

In the `backup/ca-v2-ew-tail-399f116` row the planner writes: "the canonical version lives at `82a96f0`/`82b35c0` on trunk". `git log --format='%h' -1 82b35c0` returns exit 128 — no such object exists in this repo. The planner likely meant `82daa9b` (`flavor-catalog(wave9): CA-v2 ew_tail cycle-2 verdict` — same subject as the backup tip and confirmed in trunk).

This is not just a typo — it is the load-bearing claim that justifies deleting the backup branches. Until that ref is corrected, the "obsolete intermediate state" argument is unverified.

**Action**: replace `82b35c0` with `82daa9b` (verified — same subject as backup tip, present on trunk). Also explicitly show the diff `git diff 399f116 82daa9b` in a verification snippet so a future reader can re-check.

For the record, I confirmed the backup branches share parent `e1aec33` (which IS on trunk) and add a single "CA-v2 ew_tail cycle-2 verdict" verdict commit — losslessly captured by the proposed annotated tags.

### CR-3. The `.gitignore` will silently swallow `sample_points_may18.*` if Phase 1 is run as written

Root `.gitignore` line 71 reads `!artifacts/**` (an un-ignore that re-includes everything under `artifacts/`). The planner's Phase 1 says "Discard the two `sample_points_may18.*` files unless reviewer confirms they belong". But `git check-ignore -v artifacts/sample_points_may18.csv` returns `.gitignore:71:!artifacts/** ...` — i.e. these files are NOT ignored. If the orchestrator runs Phase 1 without first deciding their fate, a subsequent `git add artifacts/` (or even some IDE auto-stage) will sweep them in. They are a 23 KB scratch CSV with a column header that does NOT match the `collaborator_*` schema (different column set, no `claim_level`/`representative_sample` fields), so they look like an earlier export iteration, not a peer-relevant artifact.

**Action**: make Phase 1 explicitly a **pre-merge gate**:
1. Orchestrator stops and asks the user: "commit, ignore (add `artifacts/sample_points_may18.*` to .gitignore), or `rm`?"
2. Only after a decision is recorded in `.orchestration/PRE_MERGE_DECISIONS.md` does any other phase proceed.
3. If the answer is "ignore", the orchestrator must add an explicit `artifacts/sample_points_may18.*` line under `!artifacts/**` (since `!` rules don't support nested re-ignore without ordering).

This gate must run BEFORE Phase 0 safety tagging — otherwise the tag captures a working state that is about to change anyway.

---

## Recommended changes (should-fix)

### RC-1. R10 is too coarse — split into R10a/R10b/R10c

R10 covers 30 commits / 18 PKA processes / ~17,500 added lines across 199 files. That is roughly 4-5x larger than the next-biggest unit (R09 is 8 commits / scaffold only). The planner asks about this in Open Question 3 and tentatively defers. My recommendation: **split**.

Suggested split (each ~10 commits, ~6 processes, evenly distributed across families):
- **R10a — Wave-1 kaon PKAs + WA + CA + v2** (K001, K002, K003, K004, K005, K006, K013) — 13 commits
- **R10b — Wave-1 beauty PKAs + WA + CA + v2** (B002, B005, B009, B011, B015) — 9 commits
- **R10c — Wave-1 EW/top + edm + charm + lepton PKAs + WA + CA + v2** (T001, T002, T010, E001, C001, L001) — 8 commits

Rationale: the planner's own per-unit physics check protocol (§C.2.1) says "sample 2-3 processes per unit". 18 processes / 2 samples = 9 sampled, but the reviewer agent then must hold the full 17.5k-line diff in context. Splitting into family-coherent units lets each Opus reviewer (a) actually load the relevant literature for that family and (b) cover the family in depth rather than skipping the rest.

The total unit count becomes 24, not 22. Wall-time estimate goes from 16.5 h to ~18 h — still fine for a multi-session run.

### RC-2. Each reviewer agent needs a literature manifest, not a verbal list

§C.2.1 names a handful of papers (Agashe hep-ph/0412089, CFW 0804.1954, BMU hep-ph/0005183, FLAG 2024, PDG 2024) but does not say where they live or how the agent retrieves them. The reviewer is a 5-hour-budget Opus run — it will not have time to scavenge papers.

**Action**: add to §C.1 an "input file" `flavor_catalog/external_research/REFERENCES.md` (the planner already mentions `flavor_catalog/external_research/` exists from Wave-7) with one section per unit (R01–R22) listing the EXACT arXiv ID, the EXACT section/eq number to consult, and — where the project already imported a Deep Research PDF — the local path. Without this, "check against literature" becomes "guess and hope", which defeats the purpose of the review.

### RC-3. Numerical-correctness check needs a concrete `pytest` selection per unit, not "and others touched"

§C.2.2 says "Run `pytest tests/test_qcd_running.py tests/test_quark_deltaf2.py tests/test_wilson_rg_audit.py tests/test_cfw_comparison.py` (and others touched by the unit)." I confirmed those four files exist and pytest collects 37 tests cleanly. But "and others touched" is the agent's job, and the agent may miss tests added in the unit's own commits.

**Action**: the COMMIT_INDEX.csv generator should also emit per-unit a `pytest_selection.txt` consisting of all `tests/*` files touched in the unit's commits (`git log --name-only RANGE -- tests/`) plus any default suite. This is one extra awk pass and removes guesswork.

### RC-4. Pre-merge safety tag should be pushed to origin BEFORE anything else, and orchestrator should refuse to proceed if push fails

Plan §A.4 Phase 0 includes `git push origin safety/pre-merge-2026-05-25 ...` — good. But Phase 0 is presented as one atomic block, so if the push silently fails (network, auth, branch-protection), the orchestrator might still proceed to Phase 1+. Also, the local tag alone is NOT a safety net — it sits on this same filesystem alongside everything being mutated.

**Action**: make the Phase 0 contract explicit: "after pushing, run `git ls-remote --tags origin | grep safety/pre-merge-2026-05-25` and abort if the output is empty. No deletion or FF until origin holds the tag."

### RC-5. Worktree at `/tmp/wa_w5b_charm_edm_worktree` needs removal in Phase 0, not Phase 5

`git worktree list` shows `/tmp/wa_w5b_charm_edm_worktree   4eca813 [wa_w5b_charm_edm_work] prunable`. The "prunable" marker means the path on disk is gone (the `/tmp` is wiped). The branch is referenced as a worktree, so `git branch -D wa_w5b_charm_edm_work` (currently in Phase 5) will fail until `git worktree prune` is run. The plan notes this in passing but does not put it in the command sequence.

**Action**: insert at the very top of Phase 0:
```
git worktree prune
git worktree list      # expect only the main worktree at the repo root
```

### RC-6. Code-consistency check should include explicit grep targets, not only "any symbol or constant"

§C.2.3 says "`git grep` for any symbol or constant added/changed by the unit". For an Opus reviewer with no prior context, this is too open-ended.

**Action**: for each high-risk unit (R04 BMU sign fix, R06 CFW comparison, R07 zero-pass UL), add an inline "grep checklist" of 3-5 specific strings/symbols to search for. For R04 the obvious targets are `evolve_deltaf2`, `C4_LR`, `C5_LR`, `bmu_to_scalar_lr`, `wilson_rg`; for R06 `cfw_matched`, `0.11`, `2.2`, `factor_2p2`, the convention-override flag name.

---

## Minor nits

- N-1. §A.4 Phase 5 deletes `flavor-catalog/2026q2` from BOTH remote and local lists, but the orchestrator currently has it checked out (`*` in `git branch -a`). The plan should explicitly include a `git checkout main` before the branch delete loop.
- N-2. The Phase 5 remote-delete loop will fail for branches that have GitHub branch protection rules; add a try/log/continue around each `git push origin --delete` so one rejection doesn't abort the loop.
- N-3. §A.4 Phase 6 has `git push --tags origin` AFTER an individual `git push origin v2026q2-catalog-complete`. The individual push is redundant.
- N-4. R19 commit count is 34 (the planner says "14 collider_rs PRIMARY processes" but the SHA list runs much longer); after removing `82a96f0` (CR-1) it is 33. Acceptable size, but worth noting.
- N-5. §D.2 `progress.json` schema does not include `safety_tags_pushed`. Add it so a resumed run can re-verify before mutation.
- N-6. Plan never tells the reviewer agent it MAY use `gh` to look at the GitHub Issues/PR view — if the project has any open issues touching a unit, they should be cross-referenced.

---

## Answers to planner's open questions

**OQ-1 (untracked files).** Commit the four `collaborator_*` files + provenance JSONs + the two `export_collaborator_*.py` scripts as the planner suggests — they have a `claim_level: operational_scan_only` provenance block and a structured schema. **Do NOT commit `sample_points_may18.csv|txt`**: column set differs from the collaborator exports, the `.txt` is a free-form Markdown-flavoured description (not a structured schema), the timestamp (May 18) and the absence of any provenance JSON pattern-match earlier scratch outputs. Move them to `snapshots/2026q2/` (gitignored) or `rm`; do not bury them under `artifacts/`. This is a user-blocking decision (see CR-3).

**OQ-2 (backup branches).** Tagging is correct. The unique commits in `backup/ca-v2-ew-tail-{399f116,776eff9}` are obsolete intermediate WA-v2 verdict states; trunk's `82daa9b` (canonical, NOT the `82b35c0` the planner cited) supersedes them. Tags `safety/ca-v2-ew-tail-399f116` and `safety/ca-v2-ew-tail-combined-776eff9` are perfectly adequate preservation since the two-line diff is fully recoverable.

**OQ-3 (R10 wave granularity).** Split. See RC-1 above.

**OQ-4 (catalog physics depth).** Spot-check default is correct, with one tightening: REQUIRE the spot-check to include at least one PRIMARY-tier process from each family touched in the wave. Otherwise a wave that adds 3 kaon + 1 beauty + 1 charm processes could pass review with only kaon sampled, leaving 2 untested families.

**OQ-5 (website R21/R22 collapsing).** Keep separate as the planner recommends. The phases-1–5 vs phases-6–12 split corresponds to a real seam (build pipeline vs presentation polish). I would also add a single concrete check to R22: `npm run build && grep -r '\\\\(' dist/ | head` to confirm LaTeX delimiters are actually being rewritten (not raw-rendered).

**OQ-6 (CLAUDE.md staleness).** The planner has this scheduled as Step 33 (post-merge). I confirmed lines 12 and 15 of CLAUDE.md reference `paper/quark-scan-2026q2` as the active paper branch and call lepton-sector "follow-up scope". The update should be: change line 12 to "active paper work lives on `main`" and keep the lepton-sector caveat. Doing it post-merge is fine; including it pre-merge would generate a no-op diff because the branch still exists.

**OQ-7 (MEG II default `br_limit = 1.5e-13`).** Apply the prohibition to BOTH `scanParams/scan.py` AND any docs that quote the number. Inconsistencies (any place that uses `4.2e-13` or `2.0e-13`) should be FLAGGED in the report (severity HIGH) but the reviewer must not propose a code change — the user explicitly wants this number frozen.

---

## Anything the planner missed

### M-1. Pre-merge `git status` snapshot

After the safety tag is created and pushed (Phase 0) but before any mutation, the orchestrator should write `.orchestration/PRE_MERGE_STATE.md` containing the exact output of `git status`, `git branch -avv`, `git tag -l`, and `git worktree list`. This is the "if everything goes wrong, what did it look like before" record. Costs nothing; pays back hugely if Phase 5 misbehaves.

### M-2. No mention of GitHub Issues mirroring local `ISSUES.md`

The user's request item 4 was "issue tracking". The plan's `ISSUES.md` is a local Markdown file. Item 5 was "GitHub mirrors local at the end". Implication: each non-ACCEPTED-RISK issue should be filed as a GitHub Issue (via `gh issue create`) after all units complete. This belongs as Step 32.5 in §E. Without it the user has no GitHub-side trail.

### M-3. No commit-policy clarification for the orchestrator itself

The plan says the reviewer agents "commit nothing" and the orchestrator does all git mutations. But there is no explicit list of WHAT commits the orchestrator will create, with whose authorship. Suggest: one commit per phase, authored as `Claude (orchestrator) <claude@anthropic.com>` co-authored with the user, message format `merge: phase-N <description>`. The Phase 1 collaborator-artifact commit (and the CLAUDE.md commit in Step 33) need this stated explicitly.

### M-4. No rollback runbook

The Risk Register (§Appendix B) mentions "revert via the safety tag is one command" but does not write that command. Suggest a 5-line `.orchestration/ROLLBACK.md` with the exact commands: `git reset --hard safety/pre-merge-2026-05-25`, `git push --force-with-lease origin main`, undo deletions via `git push origin <sha>:refs/heads/<branch>` from the tags. Pre-writing this saves the orchestrator from inventing it under pressure.

### M-5. Reviewer-agent failure mode for tests that need network or large datasets

I confirmed `pytest --collect-only` collects 37 tests in 101 s — that is suspiciously slow for collection alone (likely import-time work in `quarkConstraints/`). If a reviewer agent's per-test wall time blows past the per-unit budget (~30-60 min total), it will run out of time. The plan should say: "if `pytest` for the unit's selection exceeds 10 min, the reviewer marks numerics as PARTIAL with a note, NOT FAIL." Otherwise long-running tests force false failures.

### M-6. R10 (and likely R12, R15) overlap with the project's internal `Writer+Checker+Opus` workflow

The planner notes in OQ-4 that catalog units are already gated by an internal workflow. There is a real risk of duplicating work: if the project's own Opus round-5 already signed off 14/14 on collider_rs (per R19's commit message), then the new R19 reviewer is doing redundant work unless its mandate is explicitly "verify the verifier". I recommend §C.1 add: "The reviewer's job is NOT to re-run the project's internal Writer-Checker-Opus loop; it is to cross-check selected high-impact outputs against external literature and to verify cross-cutting code/schema consistency that single-process reviewers cannot see."

---

## Summary table for orchestrator

| Issue | Severity | Where to fix | Estimated fix time |
|-------|----------|--------------|--------------------|
| CR-1 (82a96f0 double-assigned) | CRITICAL | §B.1 R19 row | 1 min |
| CR-2 (non-existent 82b35c0 SHA) | CRITICAL | §A.1 row 4 | 1 min |
| CR-3 (untracked-file pre-merge gate) | CRITICAL | §A.4 Phase 1 + §E | 10 min |
| RC-1 (split R10) | HIGH | §B.1 + commit index | 15 min |
| RC-2 (literature manifest) | HIGH | new file under §C.1 | 30 min (one-time write) |
| RC-3 (per-unit pytest_selection) | MEDIUM | COMMIT_INDEX.csv generator | 15 min |
| RC-4 (verify tag push) | MEDIUM | §A.4 Phase 0 | 2 min |
| RC-5 (worktree prune in Phase 0) | MEDIUM | §A.4 Phase 0 | 1 min |
| RC-6 (grep checklists) | MEDIUM | §C.2.3 + per-unit rows | 20 min |
| M-1..M-6 | LOW-MEDIUM | various | 30 min total |

Total estimated planner-rework time: ~2 hours. Recommend one more planner round (Round 2) before execution begins.
