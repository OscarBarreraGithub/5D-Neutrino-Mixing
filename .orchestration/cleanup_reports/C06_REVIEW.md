# C06 Review — `git_sha` embedded in `tile_summary.json` (R07-I3)

**Commit:** `df61de2`
**Reviewer:** opus
**Date:** 2026-05-25
**Verdict:** **APPROVE**

---

## Per-check results

1. **Scope (`git show --stat`).** PASS. Six files: `scripts/run_rs_anarchy.py`
   (+29 / -0), `tests/test_run_rs_anarchy.py` (+98 / -0),
   `.orchestration/cleanup_reports/C06.md` (new, 112 lines), plus three
   bookkeeping files (`CLEANUP_QUEUE.md`, `ISSUES.md`,
   `cleanup_progress.json`). Total: +247 / -8. No collateral edits.

2. **`_resolve_git_sha()` helper.** PASS. Module-level function in
   `scripts/run_rs_anarchy.py` (new lines 159-184). Calls
   `subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=_REPO_ROOT,
   stderr=subprocess.DEVNULL)` inside a `try/except
   (CalledProcessError, FileNotFoundError, OSError)` block that returns
   `"unknown"`. Empty-output guard (`or "unknown"`) also present.
   Embedded into `summary_payload` at line 908 as a top-level
   `"git_sha"` key, alongside `config` / `elapsed_seconds` / `tiles` /
   `schema`.

3. **Test exists and exercises the field.** PASS. Four new tests in
   `tests/test_run_rs_anarchy.py`:
   - `test_resolve_git_sha_inside_repo_returns_hex` — in-repo SHA shape
     (40 hex, or `"unknown"` sentinel).
   - `test_resolve_git_sha_outside_git_tree_returns_unknown` —
     monkeypatches `_REPO_ROOT` to `tmp_path`, asserts `"unknown"`.
   - `test_tile_summary_embeds_git_sha` — drives `main()` with
     `--n-draws 1 --m-kk-tev 10 --n-workers 1`, reads the produced
     `tile_summary.json`, asserts `"git_sha"` present and well-formed.
   - `test_legacy_tile_summary_without_git_sha_still_loads` — synthesises
     a pre-C06 payload (no `git_sha`), confirms `.get("git_sha",
     "unknown")` defaults cleanly and all legacy top-level keys remain.

4. **Targeted pytest sweep.** PASS. `pytest -k 'tile or summary or scan
   or dispatch'` → **55 passed, 515 deselected, 120.68s** (no failures,
   no errors). All four C06 tests green inside the
   `tests/test_run_rs_anarchy.py` group.

5. **Backward compatibility.** PASS. Grepped for `tile_summary` readers:
   `scripts/{plot_rs_anarchy_summary,rs_anarchy_cfw_comparison,
   rs_anarchy_gate_sensitivity,rs_anarchy_mkk_min_hist_by_cvals,
   rs_anarchy_mkk_min_hist_by_yprior,compare_quark_bulk_mass_maps}.py`.
   None reference `git_sha`; all do plain `json.load` and extract
   specific keys. Adding a new top-level field is inert to them. The
   legacy-payload test exercises the `.get("git_sha", "unknown")`
   contract documented in the commit message.

6. **Bookkeeping consistency.** PASS.
   `CLEANUP_QUEUE.md` C06 PENDING → DONE with APPROVE evidence row;
   `cleanup_progress.json` C06 `"pending"` → `"done"`;
   `ISSUES.md` R07-I3 moved from open list to a new "Closed by C06"
   section with the full closure narrative. All three views agree.

---

## Wall time

Targeted pytest sweep: **123s wall** (pytest reported 120.68s in-loop).
Review activity (inspection + sweep): under five minutes.

---

**Summary (<100 words):** Forward-looking provenance fix is clean. New
helper has the required `try/except → "unknown"` fallback and an empty-
output guard. Four new tests cover in-repo, out-of-tree, end-to-end, and
legacy-payload paths. No existing reader references the new field, so
backward compatibility is structural. Targeted sweep 55/55 green in
~121s. Bookkeeping in `CLEANUP_QUEUE.md`, `cleanup_progress.json`, and
`ISSUES.md` is mutually consistent and R07-I3 is correctly relocated to
the closed-issues section. **APPROVE.**
