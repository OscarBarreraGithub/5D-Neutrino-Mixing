# C13 Review — Spurion-seed provenance comment (R02-I2)

**Commit:** `0a545f4` — `cleanup(C13): add spurion-seed provenance comment (R02-I1, R02-I2)`
**Reviewer:** Claude Opus 4.7 (1M)
**Date:** 2026-05-26
**Verdict:** APPROVE

## Scope verified

C13 closes R02-I2 (INFO, tag:code, deferred from C04): "Re-derived spurion seed
values would benefit from an in-source provenance pointer." Pure-comment change
inside `quarkConstraints/benchmarks.py::default_spurion_seed`.

## Checks performed

1. **`git show 0a545f4 --stat`** — 5 files changed, 103 insertions, 9 deletions.
   The only source file touched is `quarkConstraints/benchmarks.py` (+11, -1).
   Other four are bookkeeping: `CLEANUP_QUEUE.md`, `ISSUES.md`,
   `cleanup_progress.json`, `cleanup_reports/C13.md` (new, 84 lines).

2. **Provenance pointer present** — `quarkConstraints/benchmarks.py` lines
   190-203 carry both (a) an extended docstring referencing
   `docs/phase_logs/phase2_h4_impl.md` and noting `tests/test_quark_fit.py`
   pins the values, and (b) a 3-line `#`-comment immediately above the
   `return SpurionSeed(` literal block. Both land inside the function body,
   so grep/AST-jump on the literals hits the provenance.

3. **Literals byte-identical** — diff of the literal-bearing lines between
   `0a545f4^` and `0a545f4` is empty. Runtime confirmation:
   `up_singular_values = [0.1434281, 0.33368792, 1.23945796]`,
   `down_singular_values = [0.22903978, 0.16350921, 0.28800013]`,
   `overall_scale = 2.8`, `up_left.theta12 = -2.4040939384155906` —
   all reproduced to machine precision post-edit.

4. **Referenced artifacts exist** — `docs/phase_logs/phase2_h4_impl.md`
   (4360 bytes) and `tests/test_quark_fit.py` (23973 bytes) both resolve.

5. **Bookkeeping consistent** —
   - `CLEANUP_QUEUE.md`: C13 row flipped PENDING → DONE; Issues-Closed
     corrected from earlier typo `R02-I1` to `R02-I2` (R02-I1 is the
     unrelated `pytest.ini` MERGE_PLAN nit, routed to C18).
   - `ISSUES.md`: R02-I2 entry moved from open list to `### Closed by C13`
     with `CLOSED 2026-05-26` stamp and full closure rationale.
   - `cleanup_progress.json`: C13 status pending → done.
   - `cleanup_reports/C13.md` present (84 lines, evidence captured).

## Risk

None. Pure-comment change, no runtime/test impact, no physics or numerical
fields touched. No rebuild or test re-run required.

## Verdict

**APPROVE** — C13 closes R02-I2 cleanly; literals unchanged; bookkeeping
consistent across all four `.orchestration` files.
