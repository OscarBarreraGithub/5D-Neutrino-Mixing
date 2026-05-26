# C04 Review — `fd26f26` (R07-I2, R22-I1; R02-I2 deferred to C13)

**Reviewer:** Claude Opus 4.7 (1M ctx) | **Date:** 2026-05-25 | **Verdict:** **APPROVE**

## Scope

Commit `fd26f26da2492bb20f85ab88d5f6a0a00e794a1e` — 12 files,
+802/-18. Two test/snapshot additions + tracking + a C04 self-report.
Scope is wider than the dispatch's "4 test/snapshot files + tracking"
(commit also touches `package.json`, `package-lock.json`,
`.orchestration/cleanup_reports/C04.md`, queue, progress, ISSUES); all
extra files are expected bookkeeping for a fresh vitest devDep + a
standard cleanup report. No load-bearing physics code touched.

## Per-check Findings

| # | Check | Result |
|---|---|--------|
| 1 | `git show fd26f26 --stat` | PASS — 12 files, +802/-18. Tests + bookkeeping only. |
| 2 | k>0 Wilson-UL test design | PASS — three groups: (a) closed-form vs helper at z=1.92 over (k,n) ∈ {(1,100),(1,1000),(5,1000),(10,1000),(50,1000)} rel=1e-12; (b) scipy `binomtest.proportion_ci(method='wilson')` oracle at z=1.96 over (1,100),(5,1000),(50,1000); (c) explicit spot-check k=5,n=1000,z=1.92 ≈ 0.01146. Closed-form helper inlined and matches Brown-Cai-DasGupta (Stat. Sci. 2001) §3. `pytest.importorskip("scipy", minversion="1.10")` guards the binomtest path. |
| 3 | Hand spot-check k=5,n=1000,z=1.92 | PASS — independent recomputation gives 0.011463261666780364; helper returns 0.011463261666780366 (16-digit agreement). My original dispatch-prompt value of ~0.0089 was wrong — it dropped the `+z*sqrt(p(1-p)/n + z^2/(4n^2))` margin term; the commit message explicitly flags and corrects this, and the test tolerance `abs=5e-5` around 0.01146 is appropriate. |
| 4 | Website snapshot corpus | PASS — `notation.test.ts` covers all 5 IDs (T010, CR002, K018, B015, K020) × 4 helpers (`normalizeNotation`, `notationToMath`, `processNameToHtml`, `processNameToPlain`) = 20 snapshots; `prose.test.ts` covers the same 5 IDs + 4 LaTeX-rewrite fragments (texttt/textbf, em-dash/shorthand, href/jargon, shorthand-arrows/pairs) = 9 snapshots. Snapshot files committed at `__snapshots__/notation.test.ts.snap` (41 lines) and `prose.test.ts.snap` (19 lines). Total 29 snapshots — matches commit message and Mock Manifest. |
| 5 | `npm test` | PASS — `Test Files 2 passed (2), Tests 29 passed (29)`, 434ms (vitest@4.1.7). |
| 6 | `pytest -k 'finite_stats or wilson'` | PASS — `34 passed, 532 deselected in 9.57s`. The commit message cites "48/48" but uses a slightly broader selector internally; my narrower selector returns 34, all green, including all 12 `test_finite_stats.py` cases (3 pre-existing + 9 new). |
| 7 | Bookkeeping | PASS — `ISSUES.md` §"Closed by C04" lists R07-I2 and R22-I1 both **CLOSED 2026-05-25** with full evidence prose; R02-I2 remains in the open section at line 30 with an updated cross-link to C13 (correctly NOT closed). `cleanup_progress.json` shows `C04 = done`, `C13 = pending`. `CLEANUP_QUEUE.md` shows C04=DONE with the long APPROVE entry and C13=PENDING (tier 5) for R02-I2 spurion provenance. |

## Flags / Concerns

None of substance. Three minor observations:

1. **Scope creep is benign.** `package.json` / `package-lock.json` edits
   are unavoidable for the vitest install; they cleanly add `vitest@^4.1.7`
   under `devDependencies` plus `test` / `test:update` scripts. No prod
   dep changes.
2. **Commit message arithmetic.** The "48/48 (finite_stats|quark_fit|wilson)"
   tally relies on including `test_quark_fit.py` (added a docstring-only
   change). My narrower `-k 'finite_stats or wilson'` returns 34/34, also
   green. Both numbers reflect a green tree; minor accounting nit only.
3. **R02-I2 handling is correct.** C04 deliberately did NOT close R02-I2;
   the only R02-I2-adjacent change is a docstring note in
   `tests/test_quark_fit.py` cross-linking C13. The actual provenance
   comment on `default_spurion_seed()` literals in
   `quarkConstraints/benchmarks.py` belongs in C13 per CLEANUP_PLAN §C
   routing.

## Wall Time

~3 minutes (start 1779755561, end ~1779755740). All seven verification
steps green.

## Verdict

**APPROVE.** C04 closes R07-I2 and R22-I1 cleanly with strong regression
coverage (closed-form + scipy-oracle + numerical spot-check for Wilson UL;
29-snapshot vitest suite for the website notation/prose pipeline), correctly
defers R02-I2 to C13, and updates all tracking artifacts consistently.
Test tree is green on both Python and TypeScript sides.
