# Invalidation-Gate Re-run Sign-off

**Verdict**: PASS

Independent verification that the corrected ΔF=2 constants (see phase2_h5
and phase2_h6 sign-offs, 22.5× cumulative ε_K-budget shift) have been
propagated end-to-end through code → scans → plots → methodology note on
branch `paper/quark-scan-2026q2`.

## Verification checklist (1-10)

1. **Branch state — PASS.** `git log --oneline c83c16d..HEAD` returns
   exactly 5 commits in the expected order:
   `29803ff` (snapshot pre-audit figures), `db02223` (scan rerun),
   `87c728f` (figure rerun), `38a4586` (methodology-note headline update),
   `217af80` (invalidation-gate rerun clearance doc). Branch is
   `paper/quark-scan-2026q2` and is up-to-date with `origin/`.

2. **All 9 new scan dirs exist — PASS.** All
   `scan_outputs/rs_anarchy_*_20260515T085*/tile_summary.json` present
   (runA, runB×3, runC, run3×4). The 9 old `_20260507T*` dirs are still
   on disk (not deleted), so the pre-audit cache is preserved.

3. **Pre-audit figure snapshot — PASS.** Both
   `results/figures/quark_pre_audit_constants/` and
   `results/figures/quark/` list 49 entries and `diff` shows identical
   filenames. The pre-audit copy faithfully snapshots the canonical dir's
   structure prior to regeneration.

4. **Methodology PDF — PASS.** `pdfinfo` reports `Pages: 18`. PDF rebuilt
   cleanly via `pdflatex`.

5. **Headline numbers updated — PASS.** `grep` for `47\.26|127\.13` in
   `quark_scan_methodology_note.tex` returns 12 hits across the
   executive summary, results table, abstract, and discussion. `grep` for
   the OLD numbers `10\.58` and `28\.5` returns zero hits in the live
   methodology note — they survive only in `invalidation_gate_rerun.md`
   and phase2_h5/h6 logs, as policy requires.

6. **Pre-audit values preserved — PASS.** `followup_crossings_summary.json`
   has top-level keys `convention`, `invalidation_gate`,
   `pre_audit_reference`, `runs`; the pre-audit reference block is
   intact for diff/audit reproduction.

7. **PDG-pass count at 3 TeV unchanged — PASS.** RUNA 3 TeV tile reports
   `n_pdg_pass = 164080`, identical to pre-audit (as expected — the audit
   shifts ΔF=2 magnitudes but not the PDG-mass-and-CKM gate).

8. **ε_K still binds — PASS (semantic check).** At RUNA 3 TeV the per-tile
   `epsilon_K` pass-fraction collapses from 51.22% pre-audit to 4.09%
   post-audit (binding fraction over all draws rises from 48.78% to
   95.91%). When normalised over rejected PDG-passing draws (the report's
   convention), ε_K is the binding system in 99.500% of failures vs the
   pre-audit 93%. ε_K binds *more strongly* than before, as expected.

9. **p50/p95 magnitudes plausible — PASS.** `10.58 × 4.47 = 47.29`,
   `28.50 × 4.46 = 127.11` reproduce the new TeV headlines to <0.05 TeV.
   The empirical p50 multiplier 4.47 vs naive √22.5 = 4.74 is consistent
   with the report's caveat that the 22.5× benchmark is computed at a
   single c-point while the pooled distribution averages a broader set.

10. **Tests still pass — PASS.** `pytest -q tests/test_quark_deltaf2.py
    tests/test_wilson_rg_audit.py tests/test_qcd_running.py` →
    **35 passed in 60.57 s**.

## Physics-decision sanity (a, b, c)

(a) **CFW consistency — PASS.** New band lower edge = 47.26 − 24.98 =
    **22.28 TeV**, which sits just above CFW's 21 TeV "generic anarchic"
    bound and overlaps with the bottom of the CFW band within
    well-understood systematic differences (BGS 2020 vs older SM ε_K;
    factor-of-3 vs factor-of-1.5 PDG gate; LO vs effective-NLO Wilson
    running). The methodology note documents this overlap explicitly
    (lines 848–852, 921–922).

(b) **moreUV / moreIR / Run C finite zero-pass bounds — PASS.** Verified
    via summed `n_pdg_pass` across all tiles: all three have k=0 in their
    finite ensembles.  The Wilson-score 95% upper limits are
    `p_UL=2.304e-6` for Run 3 moreUV (N=1.6M),
    `p_UL=2.304e-6` for Run 3 moreIR (N=1.6M), and
    `p_UL=9.216e-7` for Run C (N=4.0M).  These remain c-pattern ×
    Y-prior finite-ensemble bounds independent of ε_K rescaling, as
    expected.

(c) **BGS + LO conditional caveat — PASS.** The note carries the caveat
    in multiple load-bearing places: executive paragraph
    ("conditional on the BGS 2020 central", line 578), discussion of
    BGS-budget systematic dominance over NLO (lines 706–722), CFW
    overlay disclaimer (lines 849, 876, 922), and a dedicated LO-only
    NDR-evolution disclaimer (lines 1047–1050). The asymmetric budget
    band `47.26^{+69.37}_{-24.98}` is quoted as the canonical headline.

## Outstanding deferred items

The orchestrator MUST address the following before final tag:

- **Endpoint mismatch (B/D systems at wrong scale)** — deferred from
  hole #6; running is still done at K-system endpoint conventions and is
  not yet aligned for B_d, B_s, D⁰.
- **LO-only Wilson running** — full LO is in place, but NLO/NDR
  evolution matrices for Q_LR and Q_SLL towers are not yet implemented;
  this is the dominant non-budget systematic remaining.
- **NLO BMU running (deferred from hole #6)** — implementation deferred
  to a subsequent hole; the methodology note flags this explicitly.
- **Bag-parameter endpoint alignment (deferred from hole #6)** — the
  FLAG bag parameters are quoted at scales that must be brought into a
  single consistent convention with the Wilson running.

## Final-claim status

The headline numbers in `docs/quark_scan_methodology_note.tex` are now
the LIVE post-audit values:
`M_KK^{min}(p50) = 16.54 / 47.26 TeV` (perturbative / g_s* = 3) and
`M_KK^{min}(p95) = 44.50 / 127.13 TeV`, with the asymmetric pooled-budget
band `47.26^{+69.37}_{-24.98} TeV` quoted at g_s* = 3. The BGS 2020 +
LO-only Wilson-RG caveat is documented in the executive summary,
robustness, CFW-comparison, and dedicated-LO sections. The hole-#5 and
hole-#6 deferred items (NLO BMU running, B/D endpoint alignment, bag
parameter conventions) are flagged in-text and feed directly into the
next holes; they do not invalidate the present sign-off.

## Recommendation

**Proceed to hole #7 (CFW comparison).** All 10 verification checks PASS,
all three physics-decision sanity checks PASS, and the band overlaps
CFW's 21 TeV edge — making hole #7 (a literature-comparison hole) the
natural next step. No halt warranted; deferred items are tracked and
will be picked up in subsequent holes per the
`paper_execution_decisions.md` policy.
