# Cleanup Queue — 22 atomic units (17 Opus + 4 deterministic + 2 deferred)

Statuses: PENDING / IMPL-IN-PROGRESS / REVIEW-IN-PROGRESS / DONE / DEFERRED / BLOCKED
Lane: opus / det / deferred

| Unit | Title | Tier | Lane | Issues Closed | Status | Verdict |
|------|-------|------|------|---------------|--------|---------|
| C01      | kaon constants SSOT + collaborator re-export | 1 | opus | R03-I1, R03-I2 | DONE | APPROVE-WITH-NOTES: option-(b) retry; commits `3ab1f8f` (C01a, literals+pin test) + `7899205` (C01b, re-export); 179 tests pass; epsilon_K ratio tightens ~7x in collaborator CSVs (planner forecast ~6.24x; realized number folds in operator-mix weighting); firewall preserved; full evidence in `.orchestration/cleanup_reports/C01.md` |
| C02a-code| --epsilon-k-budget CLI flag + band paragraph  | 1 | opus | R03-I3 (tasks 1-4) | DONE | APPROVE: `--epsilon-k-budget {central,low,high}` plumbed end-to-end (script → EnsembleConfig → `_worker_init` → `evaluate_delta_f2_constraints` → `evaluate_epsilon_k`); 12/12 new unit tests + 125/125 regression tests pass; smoke run confirms `central` reproduces pre-flag bit-for-bit and `low`/`high` rescale ε_K ratio by exactly 6.7/0.2233; methodology-note paragraph adds symbolic-scaling band $47.26^{+75.10}_{-25.05}$ TeV with explicit deferral of measured edges to C02c; PDF rebuilt cleanly (20 pp); evidence in `.orchestration/cleanup_reports/C02a-code.md` |
| C03      | Wilson-RG / BMU follow-ups                    | 1 | opus | R04-I1..I4 | DONE | APPROVE: R04-I1 (LOW) defensive upper-tri guard added to `scripts/audit_wilson_rg.py::scalar_lr_segment_matrix`; R04-I2 (LOW) function-local `from .qcd_running import evolve_deltaf2_wilsons` in `quarkConstraints/deltaf2.py::_evolve_wilsons` lifted to module top (unconditional per M-7, non-circular); R04-I3 (INFO) verified `paper_0710_1869` LR-sign convention consistent with post-C01 (no mismatch, file:line evidence in report); R04-I4 (INFO) appended `Status as of 2026-05-25 post-cleanup` block to `docs/audits/wilson_rg_inventory.md` tracking 4 deferred follow-ups (forwarded to C19); focused suite 103 passed, broader sweep 196 passed + 1 skipped; audit script reproduces `1.300e-16` max relative discrepancy; full evidence in `.orchestration/cleanup_reports/C03.md` |
| C04      | k>0 test generalization + website snapshots   | 2 | opus | R07-I2, R22-I1 (R02-I2 deferred to C13 per CLEANUP_PLAN §C) | DONE | APPROVE: R07-I2 closed via three new test groups in `tests/test_finite_stats.py` — (a) `test_wilson_upper_limit_k_gt_zero_matches_closed_form` parametrized over (k,n) ∈ {(1,100),(1,1000),(5,1000),(10,1000),(50,1000)} cross-checking helper vs in-file closed-form at z=1.92 (rel 1e-12); (b) `test_wilson_upper_limit_k_gt_zero_matches_scipy_at_z_eq_1p96` cross-checking against `scipy.stats.binomtest(k,n).proportion_ci(method='wilson', confidence_level=0.95).high` for (1,100),(5,1000),(50,1000) with z=`norm.ppf(0.975)`=1.959964; (c) spot-check `wilson_upper_limit(5,1000,z=1.92) ≈ 0.01146` (dispatch-prompt hand-computation of 0.0089 was off by the `+z*sqrt(...)` margin term — actual closed-form 0.011463); R22-I1 closed via vitest snapshot suite at `flavor_catalog/website/src/lib/__tests__/{notation,prose}.test.ts.snap` (29 snapshots — T010/CR002/K018/B015/K020 raw `process_name` × 4 notation helpers + 5 process_name + 4 LaTeX-rewrite prose fragments × `texProseToHtml`), `vitest@^4.1.7` added to devDependencies, `npm test`/`npm run test:update` scripts added; R02-I2 deferred to C13 (CLEANUP_PLAN §C routes R02-I2 to C13 "Spurion-seed provenance comment"; C04 added a supplementary docstring note in `tests/test_quark_fit.py` cross-referencing C13 but did not close the issue); focused suite 48 passed, `npm test` 29 passed; full evidence in `.orchestration/cleanup_reports/C04.md` |
| C05      | external_research MANIFEST + sha256           | 3 | opus | R17-I1 | PENDING | - |
| C06      | git_sha embedded in tile_summary.json         | 3 | opus | R07-I3 | PENDING | - |
| C07      | per-master-compile annotations + agg script   | 3 | opus | R17-I2, R18-I3, R19-I4 | PENDING | - |
| C08      | catalog open_issues thread closures           | 3 | opus | R10a-I3, R10b-I3, ... | PENDING | - |
| C09      | YAML/tex typos + status_history gaps          | 4 | opus | R10b-I1, R10b-I2, R11-I3 | PENDING | - |
| C10      | T010/T011 + sha256 backfill + worklog names   | 4 | opus | R10c-I1, R10c-I2, R10c-I3 | PENDING | - |
| C11      | DA-N worklog closure addenda                  | 4 | opus | R11-I2, R12-I4, R15-I4 | PENDING | - |
| C12      | v0.X compile-report metadata reconcile        | 4 | opus | R18-I1, R19-I3, R21-I1, R22-I3 | PENDING | - |
| C13      | R02-I2 spurion provenance comment             | 5 | opus | R02-I1 | PENDING | - |
| C14      | codex model version drift fix                 | 5 | opus | R14-I3, R20-I2 | PENDING | - |
| C15      | pytest_selection backfill                     | 5 | det  | R05-I2 | PENDING | - |
| C16      | SESSION_NOTES tally drift                     | 5 | opus | R20-I1, R20-I3 | PENDING | - |
| C17      | redundant .gitkeep cleanup                    | 5 | det  | R09-I2 | PENDING | - |
| C18      | MERGE_PLAN retroactive corrections            | 5 | det  | R12-I1..I3, R13-I1, R14-I1, R15-I1, R16-I3, R19-I1, R19-I2, R19-I5 | PENDING | - |
| C19      | final docs sweep + methodology PDF rebuild    | 5 | opus | R08-I3, residual docs | PENDING | - |
| C20      | move INFO issues to ACCEPTED-RISK             | 5 | det  | R01-I2, R04-I4, R07-I3, R07-I4, others | PENDING | - |
| C02b     | sensitivity-band figure                       | 6 | deferred | R03-I3 (task 5) | DEFERRED | - |
| C02c     | SLURM RUNA reruns at 3 budget edges           | 6 | deferred | R03-I3 (followup) | DEFERRED | - |
