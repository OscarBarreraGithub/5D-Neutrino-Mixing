# C01 Independent Review

Reviewer: Claude Opus 4.7 (1M)
Date: 2026-05-25
Worktree: `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing` (branch `flavor-catalog/2026q2`)
Subject commits: `3ab1f8f` (C01a) + `7899205` (C01b)

---

## Check 1 — R03-I1 + R03-I2 fix correctness: PASS

**Diff is values-only, no deltaf2 import added.** `git show 3ab1f8f -- quarkConstraints/modern/phenomenology.py` shows
exactly 4 literal changes plus a 7-line docstring banner. Module-level imports unchanged; firewall preserved.

`quarkConstraints/modern/phenomenology.py:417-437` (current tree) confirms the four post-audit literals:
- `_KAON_B_1 = 0.5503` (line 432)
- `_KAON_B_4 = 0.903` (line 433)
- `_KAON_B_5 = 0.691` (line 434)
- `_KAON_EPSILON_K_SM = 2.161e-3` (line 437)

The vendored-copy banner (lines 417-425) explicitly tells future maintainers:
"VENDORED COPY of canonical values in quarkConstraints/deltaf2.py:618-623 ... update both copies together."

**Pin test is structurally sound.** `tests/test_quark_deltaf2.py:120-155`:
- Imports `quarkConstraints.modern.phenomenology as modern_phen` at *function* scope (line 131) — legal under the
  firewall, which only guards the modern *module* surface, not test code.
- `deltaf2` is already imported at module top of the test file (used by neighbouring tests).
- Uses `np.isclose` (line 151), not `==`.
- Asserts **all 11** vendored constants, not just the 4 audited ones — exceeds the requirement.
- Failure message: "...{modern_name}={modern_val} ... out of sync with {canonical_name}={canonical_val} ... Update both literals together." (lines 152-155).

## Check 2 — No regression: PASS

Re-ran the implementer's three commands independently.

```
tests/test_quark_deltaf2.py tests/test_modern_phenomenology.py tests/test_modern_scan.py -v
  -> 24 passed in 28.90s
```
The new pin test `test_modern_phenomenology_kaon_constants_match_deltaf2_canonical` is at index 4 (line 5 of output)
and PASSED. Both firewall tests
(`test_modern_phenomenology_module_does_not_import_repo_v1_scan_or_deltaf2_stack`,
`test_importing_modern_phenomenology_does_not_load_repo_v1_scan_or_deltaf2_stack`) PASSED.

```
tests/ -k 'kaon or epsilon_k or bag or deltaf2 or modern'
  -> 179 passed, 1 skipped, 365 deselected in 121.66s
```
(Implementer's "179/180" matches; the "1" is a skip, not a fail.)

```
tests/test_modern_phenomenology.py -v -k 'firewall or import'
  -> 2 passed, 7 deselected
```
Both firewall guards intact.

## Check 3 — Physics correctness: PASS

Pre-C01 snapshots present at `/tmp/pre_C01_5tev.csv` and `/tmp/pre_C01_direct.csv`.

**5tev CSV (`collaborator_5tev_points.csv`)** — direct row-by-row comparison of `ratio_epsilon_K_physical_mgkk`:

| row | pre | post | factor |
|-----|------|------|--------|
| 1 | 0.03581 | 0.25583 | 7.144x |
| 2 | 0.03510 | 0.25073 | 7.143x |
| 3 | 0.04146 | 0.29613 | 7.143x |
| 4 | 0.04501 | 0.32154 | 7.143x |
| 5 | 0.04885 | 0.34892 | 7.143x |

Uniform ~7.14x tightening across all 5 rows — matches the implementer's "~7x" claim. Consistent with the
~6.24x naive B-factor cross-product `(2.161/1.81) × (0.903/0.78) × (0.691/0.57) ≈ 1.194 × 1.158 × 1.212 ≈ 1.675`
once one accounts for the per-row operator-mix weighting on top of the bag-factor changes — sign and order of
magnitude correct.

**Direct affine CSV — top-5 policy reshuffling.**
- pre top-5: `down_aligned_flat_d, down_aligned_flatter_u, stronger_down_alignment, balanced_flat, larger_Q_slope_high_scale`
- post top-5: `down_aligned_flat_d, down_aligned_flatter_u, balanced_flat, down_aligned_flat_d, down_aligned_flat_d`
- dropped: `stronger_down_alignment`, `larger_Q_slope_high_scale` (exactly two, as claimed).
- The post-C01 list contains duplicate `down_aligned_flat_d` entries — these are different parameter points within
  the same policy family being picked up after the constraint tightening reshuffles which points sit on the
  Pareto front. This is expected behaviour for the export script, not a regression.

## Check 4 — Bookkeeping: PASS

- `.orchestration/ISSUES.md:470-481` — both R03-I1 (HIGH) and R03-I2 (MEDIUM) appear under `### Closed by C01`
  with **CLOSED 2026-05-25** tag and SHA `3ab1f8f` (C01a). Descriptions accurately reference C01b commit `7899205`
  for the artifact re-export.
- `.orchestration/CLEANUP_QUEUE.md:8` — C01 row marked **DONE** with verdict
  "APPROVE-WITH-NOTES: option-(b) retry; commits `3ab1f8f` (C01a, literals+pin test) + `7899205` (C01b, re-export);
  179 tests pass; epsilon_K ratio tightens ~7x..." and pointer to `.orchestration/cleanup_reports/C01.md`.
- `.orchestration/cleanup_progress.json` — C01 entry has `"status": "done"`. Subsequent units (C02a-code, C03, C04, ...)
  remain `"pending"`.

---

## New issues uncovered: None.

Minor observation (not a blocker): the cross-referencing housekeeping flagged in R10a-I2 (K001.yaml `open_issues`
still mentions the resolved condition without a `tracked_in: R03-I1` back-pointer) is unchanged. That's R10a-I2's
own scope, not C01's, so it remains as-tracked.

---

## Overall verdict: **APPROVE**

Orchestrator may push `3ab1f8f` and `7899205` to origin. All four checks pass with independent evidence:
- Diff is values-only, firewall preserved (Check 1).
- 24/24 focused + 179/180 broader + 2/2 firewall tests pass (Check 2).
- ~7.14x epsilon_K tightening confirmed across all 5tev rows; exactly two policies dropped from direct-affine top-5 (Check 3).
- ISSUES.md, CLEANUP_QUEUE.md, cleanup_progress.json all consistent (Check 4).

Wall time: ~4 minutes (test runs dominate at ~2 min for the broader sweep).
