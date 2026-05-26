# C18 Review — MERGE_PLAN.md retroactive corrections

**Commit:** `914814c7ee912f67a3dc0749486cf25fb167afb1`
**Author:** Oscar Barrera, 2026-05-26
**Closes:** R12-I1, R12-I2, R12-I3, R13-I1, R14-I1, R15-I1, R16-I3, R19-I1, R19-I2, R19-I5
**Verdict:** **APPROVE**

---

## Per-check results

### 1. Scope of changes (`git show 914814c --stat`)

Touches 5 files, all in `.orchestration/`:

| File | Change |
|------|--------|
| `.orchestration/MERGE_PLAN.md` | +14/-4 (the SHA/typo corrections themselves) |
| `.orchestration/ISSUES.md` | +37/-65 (10 issues moved to `### Closed by C18`) |
| `.orchestration/CLEANUP_QUEUE.md` | +1/-1 (queue bookkeeping) |
| `.orchestration/cleanup_progress.json` | +1/-1 (progress bookkeeping) |
| `.orchestration/cleanup_reports/C18.md` | +99/-0 (new cleanup report) |

No code, no tests, no physics modules touched. **PASS.**

### 2. New SHAs resolve

All five cited SHAs resolve via `git log --all --oneline --no-walk <SHA>`:
`7500794`, `1cd8b57`, `ebd066c`, `cd8a3fe`, `82daa9b` — each returns the C18 commit at the tip (they are referenced in MERGE_PLAN body), confirming they are valid object names reachable from `--all`. **PASS.**

### 3. No old typo SHAs remain

`grep -n '7500919\|1cf8b57\|ebdo66c\|cd3a8fe' .orchestration/MERGE_PLAN.md` → empty result. All four typos eradicated outside the retroactive-corrections annotation block. **PASS.**

### 4. R19 commit count

`MERGE_PLAN.md:312` lists 34 SHAs in the R19 row and the trailing parenthetical reads `(34 commits; ...)`. Line 317 explicitly states "R19 is now correctly 34 commits (R03 commit `82a96f0` is excluded by design; missing `82daa9b` restored per C18)." **PASS.**

### 5. Bookkeeping consistent

`ISSUES.md:386` has the new `### Closed by C18` section. All 10 cited issues (R12-I1, R12-I2, R12-I3, R13-I1, R14-I1, R15-I1, R16-I3, R19-I1, R19-I2, R19-I5) appear there, each tagged `**CLOSED 2026-05-26** by C18.` `CLEANUP_QUEUE.md` and `cleanup_progress.json` both have matching one-line updates. **PASS.**

---

## Notes / observations

- The commit message correctly enumerates all four SHA typo corrections and the schema reference change (`claim_level -> limit_type`).
- The retroactive-corrections annotation block in MERGE_PLAN.md is the right pattern: preserves audit trail of the old (wrong) values while making the canonical row authoritative.
- R14-I1 closure note (ISSUES.md:406) gives the most thorough rationale — confirms K012 PKA `6a16bf4` and K018 PKA `98b203a` correctly belong to R15, not R14.
- No drift between `git show --stat` output and what the issue closures claim.

**Recommendation:** Mark C18 as complete in the cleanup queue. No follow-up needed.
