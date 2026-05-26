# C08 Review — open_issues thread closures (Wave-1 yamls)

**Commit:** `6e1fd9d` ("cleanup(C08): close open_issues threads across wave-1 yamls (R10a-I2, R10a-I3, R10b-I3)")
**Reviewer:** opus / det
**Date:** 2026-05-25
**Closes:** R10a-I2, R10a-I3, R10b-I3

## Verdict: **APPROVE**

Docs-only catalog cleanup. All 7 Wave-1 yamls touched cleanly; every `open_issues` entry now carries an unambiguous disposition marker with verbatim original text preserved. Bookkeeping (CLEANUP_QUEUE, ISSUES, cleanup_progress.json) consistent.

## Per-check results

### Check 1 — `git show --stat`: 7 yamls + bookkeeping
PASS. Diff touches exactly the expected 7 yamls (K001, K003, B002, B005, B009, B011, B015) plus 4 bookkeeping files (`CLEANUP_QUEUE.md`, `ISSUES.md`, `cleanup_progress.json`, new report `C08.md`). 11 files / +146 / -30. No code, no test, no schema change.

### Check 2 — Status markers in every `open_issues` entry
PASS. Every list item now carries one of four explicit markers, with original text preserved verbatim after `Original note: `:

| File | Item | Marker |
|---|---|---|
| K001 | 1 (PDG review vs pdgLive snapshot path) | `[RESOLVED-by-C08 2026-05-25]` |
| K001 | 2 (legacy/modern epsilon_K mismatch) | `[RESOLVED-by-C08 ..., tracked-in R03-I1, closed-by cleanup-C01 commit 3ab1f8f]` |
| K003 | 1 (eps'/eps listing vs DataBlock headline) | `[RESOLVED-by-C08 2026-05-25]` |
| B002 | 1 (sin(2 beta) all-charmonium vs J/psi K_S) | `[RESOLVED-by-C08 2026-05-25]` |
| B005 | 1 (PDG live vs HFLAV Apr-2023) | `[RESOLVED-by-C08 2026-05-25]` |
| B009 | 1 (HFLAV Dec-2025 vs PDG 2025) | `[RESOLVED-by-C08 2026-05-25]` |
| B009 | 2 (direct f_B|V_ub| construction) | `[DEFERRED-to-post-paper-finalization 2026-05-25 (C08)]` |
| B011 | 1 (Misiak 2020 vs PDG-reviewed 2015) | `[RESOLVED-by-C08 2026-05-25]` |
| B015 | 1 (rejected-anchor provenance note) | `[CLOSED-as-N/A-by-C08 2026-05-25]` |
| B015 | 2 (B016-B019 scope boundary) | `[CLOSED-as-by-design 2026-05-25 (C08)]` |

K001 item 2's `tracked-in R03-I1, closed-by cleanup-C01 commit 3ab1f8f` cross-link is exactly the fix R10a-I2 requested. The `[DEFERRED-...]` and `[CLOSED-as-...]` markers are well-chosen for the B009 item 2 and B015 sub-items that are not strictly "resolved by decision."

### Check 3 — YAML parse
PASS. `python -c "import yaml; [yaml.safe_load(...) for p in [...]]; print('OK')"` → `OK`. Schema stays `flavor_catalog.process.v1`; `open_issues` remains `list[str]`; item counts unchanged.

### Check 4 — Bookkeeping consistency
PASS.
- `CLEANUP_QUEUE.md`: C08 row flips `PENDING → DONE` with full APPROVE notes; tracked-issues column corrected from `R10a-I3, R10b-I3, ...` to the precise `R10a-I2, R10a-I3, R10b-I3`.
- `ISSUES.md`: R10a-I2, R10a-I3, R10b-I3 stanzas moved from "Open" to a new `### Closed by C08` block, each with **CLOSED 2026-05-25** marker, closure note, and evidence pointer to `cleanup_reports/C08.md`.
- `cleanup_progress.json`: C08 entry flipped `"pending" → "done"`.
- Commit message correctly carves out R10b-I1, R10b-I2, R11-I3 (left for C09) and is consistent with the `CLEANUP_QUEUE.md` C09 row.

## Notes
- Forward-looking: the per-item markers (`[RESOLVED-by-C08]`, `[DEFERRED-to-...]`, `[CLOSED-as-...]`) establish a clean convention future cleanup units can reuse for catalog-level disposition records.
- No website / no test impact: no consumer of `open_issues` in `tools/` or `flavor_catalog/website/` — purely human-readable catalog metadata.

**Result:** APPROVE C08.
