# C03 Review — Wilson-RG / BMU follow-ups (R04-I1..I4)

**Commit:** `9dc72d7` ("cleanup(C03): Wilson-RG / BMU follow-ups (R04-I1..I4)")
**Author:** Oscar Barrera, 2026-05-25 20:05:11 -0400
**Reviewer:** Claude Opus 4.7 (1M context), 2026-05-25
**Verdict:** **APPROVE**
**Wall time:** ~62 s

---

## 1. Scope (`git show 9dc72d7 --stat`)

7 files changed, +497 / -45. All within declared scope:

| File | Δ | Purpose |
|------|---|---------|
| `scripts/audit_wilson_rg.py` | +18 | R04-I1 defensive triangular guard |
| `quarkConstraints/deltaf2.py` | +6/-1 | R04-I2 import lift |
| `docs/audits/wilson_rg_inventory.md` | +57 | R04-I4 status block |
| `.orchestration/cleanup_reports/C03.md` | +212 | unit report |
| `.orchestration/ISSUES.md` | +43/-38 | R04-I1..I4 moved to Closed-by-C03 |
| `.orchestration/CLEANUP_QUEUE.md` | +29 | C03 row → DONE |
| `.orchestration/cleanup_progress.json` | +138 | status bump |

No out-of-scope files touched. **PASS.**

## 2. `scripts/audit_wilson_rg.py` diff (R04-I1)

The audit script previously hard-coded `gamma44, gamma45, gamma55` and used the
closed-form upper-triangular propagator without asserting the precondition. C03
adds a named local `gamma54 = 0.0`, materializes the 2x2 ADM, and raises
`ValueError("audit_wilson_rg.scalar_lr_segment_matrix requires an
upper-triangular scalar LR ADM (gamma_54 == 0); …")` if either
`gamma_lr_local[1,0] != 0.0` or `np.allclose(gamma_lr_local - np.triu(...), 0)`
fails.

The check is correct: the closed-form `f44 = eta**(gamma44/(2 beta0))`,
`f55 = eta**(gamma55/(2 beta0))` form assumes strictly upper-triangular ADM, so
the guard exactly matches the precondition. Numeric output is preserved
(reported `1.300e-16` max relative discrepancy is unchanged because
`gamma54 = 0.0` is the same literal value). **PASS.**

## 3. `quarkConstraints/deltaf2.py` diff (R04-I2)

The function-local `from .qcd_running import evolve_deltaf2_wilsons` inside
`_evolve_wilsons` is lifted to the module-top import block (line 20,
alongside the other intra-package imports). The docstring is updated with an
inline provenance note explaining why the lift is safe.

**Circularity check:** grepped `quarkConstraints/qcd_running.py` for any
`deltaf2` reference at module scope — none exists (`qcd_running.py` module-top
imports are `math` and `numpy` only; the string `deltaf2` appears solely in the
function name `evolve_deltaf2_wilsons` and the `__all__` entry). The lift is
unconditional and non-circular. **PASS.**

Minor note: the commit message line "lift function-local deltaf2 import to
module level" is mildly ambiguous about which direction (the orchestrator
prompt's R04-I2 wording originally flipped the file). C03 explicitly resolves
this in both the diff and the ISSUES.md closure: the offender lived in
`deltaf2.py` (importing from `qcd_running`), which matches `CLEANUP_PLAN.md`.

## 4. `docs/audits/wilson_rg_inventory.md` status block (R04-I4)

A new `## Status as of 2026-05-25 post-cleanup (C03)` section is inserted
immediately before the existing `## Invalidation-gate note`. It enumerates the
four explicitly deferred follow-ups (items 2, 4, 7, 10 of the existing 10-item
checklist), records the R04-I3 verification with file:line evidence
(`matching_kkgluon.py:525-526`, `hadronic.py:2237,2244`,
`observables.py:917-918`, `rg.py:98-101`, `qcd_running.py:50-57`,
`audit_wilson_rg.py:85-117`), and pins the audit-script numeric output
(`C4_LR=1 -> 3.53816397486`, `1.300e-16` max relative discrepancy). Reassigns
substantive resolution to **C19** per `CLEANUP_PLAN.md` §C C19. **PASS.**

## 5. Test re-run

```
$ timeout 600 python -m pytest tests/ -k 'wilson or rg or qcd or audit or bmu' --tb=short
…
tests/test_qcd_running.py ............................                   [ 72%]
tests/test_wilson_rg_audit.py ..                                         [100%]
===================== 103 passed, 454 deselected in 24.25s =====================
```

All 103 selected tests green in 24.25 s, well under the 600 s budget. **PASS.**

## 6. Bookkeeping

- `ISSUES.md` §"Closed by C03" present with explicit closure prose for all
  four issues (lines 460–476). Each closure record names file, function/line,
  and the test evidence. **PASS.**
- `CLEANUP_QUEUE.md` C03 row marked **DONE** with full APPROVE blurb
  (line 10). **PASS.**
- `cleanup_progress.json` queue entry `{"unit": "C03", "lane": "opus",
  "tier": 1, "status": "done"}`. **PASS.**
- `.orchestration/cleanup_reports/C03.md` (212 lines) is the unit report. **PASS.**

## Flags

None. The four R04 deferred follow-ups are now formally tracked under C19 in
`CLEANUP_PLAN.md`, with the C03 doc addendum acting as the bridge — this is
the agreed handoff and is not a regression.

## Verdict

**APPROVE.** All R04-I1..I4 closures match the commit, the diffs are minimal
and defensive, the import lift is non-circular, the audit script's numeric
contract is preserved, and the test suite is green.
