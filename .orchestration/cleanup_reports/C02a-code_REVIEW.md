# C02a-code REVIEW

**Reviewer:** Opus 4.7 (1M context)
**Implementer commit:** `3e1e07c1a4285a40f63a1294aee8328bc5ed2040`
**Date:** 2026-05-25
**Verdict:** **APPROVE-WITH-CONCERN** (concern is a +0.5% band-edge rounding
discrepancy on the low edge, explicitly flagged in the orchestration brief
as "worth flagging but not blocking"; the symbolic-scaling paragraph already
labels itself approximate and defers the direct verification to C02c.)

---

## Check 1 — Diff scope reasonable? PASS

`git show 3e1e07c --stat` reports 5 files / +436 / -5, matching the claimed
surfaces exactly:

- **`scripts/run_rs_anarchy.py`** (+78/-1): Only adds the
  `EPSILON_K_BUDGET_EDGES` constant, the two new `EnsembleConfig` fields
  (`epsilon_k_budget_edge`, `epsilon_k_np_budget_override`), plumbing through
  `_worker_init` and `_evaluate_one_draw`, the `--epsilon-k-budget`
  argparse flag, and the `convention` / `config` JSON keys. **`central`
  resolves to `None`, preserving bit-for-bit pre-flag behavior** (verified
  in diff and in test `test_evaluate_epsilon_k_default_budget_matches_implicit`).
- **`quarkConstraints/deltaf2.py`** (+47/-3): Adds keyword-only
  `epsilon_k_np_budget_override` to `evaluate_epsilon_k`,
  `evaluate_epsilon_k_with_running`, `_hadronic_eval_for_system`, and
  `evaluate_delta_f2_constraints`. `None` → no-op default;
  positive float → replaces `abs(EPSILON_K_EXP - EPSILON_K_SM)`;
  zero/negative → `ValueError`. **No other constants touched**
  (`EPSILON_K_EXP`, `EPSILON_K_SM`, `KAPPA_EPSILON`, `DELTA_M_K` unchanged;
  no value drift). **B_d, B_s, D systems explicitly unaffected** (override
  is gated on `key == "epsilon_k"`).
- **`tests/test_run_rs_anarchy.py`** (new, 240 lines): 12 tests covering
  the edge table contract, argparser default + parsing + rejection of
  unknown labels, `EnsembleConfig` round-trip, deltaf2 override
  propagation (including the nonpositive `ValueError`), and **explicitly
  the `_worker_init` rehydration round-trip** (lines 203-241) — the
  multiprocessing seam the implementer self-reports as the worker-init
  bug they caught and fixed.
- **`docs/quark_scan_methodology_note.tex`** (+64): Single new
  `\paragraph{Symbolic NP-budget band on $\Mkk^{\min}$}` block at
  lines 741-793. Cites C02c by name (line 791) for the deferred direct
  scans. C02b reference comes via the ISSUES.md closure note rather
  than the tex itself, which is acceptable (the tex defers the
  *measurement*; the figure deferral lives in bookkeeping).
- `docs/quark_scan_methodology_note.pdf` rebuilt (620832 → 626853 bytes,
  +6 KB ≈ one new paragraph + page).

## Check 2 — Tests pass? PASS

`python -m pytest tests/ -k 'scan or budget or epsilon or rs_anarchy or
deltaf2 or wilson' --tb=short` → **125 passed, 432 deselected in 72.07s**.
Matches the claimed count exactly. The 12 new `test_run_rs_anarchy.py`
tests are all in the passing set (line `tests/test_run_rs_anarchy.py
............ [ 84%]`).

## Check 3 — Physics paragraph correct? APPROVE-WITH-CONCERN

Verified in `docs/quark_scan_methodology_note.tex:741-793`:

- **Central 47.26 TeV** cited (line 773). MATCH.
- **Symbolic $1/\sqrt{|\Delta\eK|}$ scaling** written explicitly as
  `\sqrt{\Delta\eK^{\mathrm{central}}/\Delta\eK}` in lines 759-763, with
  the $\propto 1/\Mkk^2$ Wilson-coefficient proof sketch citing
  `phase2_h5_signoff.md:90-101`. MATCH.
- **C02c citation:** explicit at line 791 (`"deferred cleanup unit C02c
  in .orchestration/CLEANUP_PLAN.md"`), with the verb "will replace this
  paragraph's symbolic edges with measured band ends in the paper
  finalization pass." MATCH.
- **PDF rebuilt:** `docs/quark_scan_methodology_note.pdf` mtime
  `2026-05-25 19:37:16`, post-tex mtime `19:32:37`, pre-commit
  `19:46:21` — the committed PDF is the up-to-date two-pass rebuild.
  Size matches the +6 KB the commit message anticipates.

**FLAG (non-blocking):** Reproducing the band arithmetic from first
principles:

| Edge   | Ratio        | Computed value | Computed delta | Reported delta |
|--------|--------------|----------------|----------------|----------------|
| High   | 2.5884       | 122.33 TeV     | +75.07         | **+75.10**     |
| Low    | 0.4726       | 22.33 TeV      | −24.93         | **−25.05**     |

The high edge agrees to 0.03 TeV (printing precision). The **low edge
shows a 0.12 TeV / 0.5% discrepancy** (22.21 reported vs 22.33 computed,
i.e. −25.05 vs −24.93). This is exactly the discrepancy the orchestration
brief predicted, and it is well inside the "approximate, symbolic-scaling"
self-label of the paragraph (lines 781-785). Most likely cause: rounding
inside the implementer's two-step product (`47.26 * 0.4725 = 22.3304` if
the loose ratio is truncated to 0.4725 instead of 0.4726). Not worth
blocking — the direct C02c scans will replace these numbers — but the
implementer should be aware before final paper submission.

## Check 4 — Bookkeeping PASS

- **`.orchestration/ISSUES.md:474-482`** — `### Closed by C02a-code`
  section present; R03-I3 marked `**CLOSED 2026-05-25** by C02a-code
  (commit 3e1e07c)`. Closure SHA, date, AND the explicit forwarding of
  task 5 → C02b and task 2 (three RUNA reruns) → C02c are all written
  out (lines 480-481). Evidence pointer to `C02a-code.md` and the
  changed surfaces is included.
- **`.orchestration/CLEANUP_QUEUE.md:9`** — C02a-code row marked DONE
  with full per-check evidence string; C02b (line 28) and C02c (line 29)
  rows present and marked DEFERRED. MATCH.
- **`.orchestration/cleanup_progress.json`** — C02a-code → `"done"`
  (line 15); C02b and C02c → `"deferred"` (lines 126-136). MATCH.
- **C02c-deferred check:** confirmed in both `CLEANUP_QUEUE.md:29` and
  `cleanup_progress.json:131-136`.

---

## Summary

The CLI flag is plumbed cleanly through every layer the implementer
claims (script → EnsembleConfig → `_worker_init` → deltaf2), the
no-op `central` path is structurally guaranteed by the `None` sentinel,
the new tests pin the worker-init rehydration seam that was the
self-reported bug class, the broader 125-test regression suite is
green, and the methodology-note paragraph correctly labels itself
symbolic and defers measured edges to C02c. The only finding is a
0.5% rounding mismatch on the low band edge, already anticipated by
the orchestration brief and shielded by the paragraph's explicit
"approximate, symbolic-scaling" disclaimer.

**Verdict: APPROVE-WITH-CONCERN.** Ship.

Suggested follow-up (non-blocking, for C02c): when the direct RUNA
reruns replace the symbolic numbers, also recompute the low-edge
product at full precision so the paper's final band is internally
self-consistent to all printed digits.
