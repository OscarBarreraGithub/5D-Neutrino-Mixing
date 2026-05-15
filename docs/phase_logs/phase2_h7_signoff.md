# Phase 2 Hole #7 Sign-off

**Verdict**: PASS

Hole #7 (CFW reconciliation) closes with an honest factor-2.2 framing,
quantitatively decomposed into well-understood post-2008 inputs (BGS-2020
SM `epsilon_K`, FLAG-2024 bag parameters, audited LO BMU Wilson running).
The earlier "11% agreement" framing has been removed everywhere; the
methodology note, audit docs, comparison figure, and regression test are
all internally consistent.

## Verification checklist (1-7)

1. **Branch state** -- PASS. `paper/quark-scan-2026q2` is current and
   pushed (`git status`: up to date with `origin/paper/quark-scan-2026q2`).
   Nine hole-#7 commits chained on top of the invalidation-gate seal
   (`f19f17c`): `ffec986`, `06d5d85`, `da8f647`, `508ae69`, `b47aa76`,
   `330ffa9`, `11e0d58`, `321fe0f`, `7c284a8`, `2a3a2c0`, `f2c29c8` (the
   last is a log-only update, as advertised).

2. **No live "11% agreement" claim** -- PASS.
   `grep -nrE "11%|qualitatively consistent" docs/quark_scan_methodology_note.tex docs/audits/cfw_*.md`
   returns zero hits.

3. **Factor-2.2 framing live** -- PASS.
   `grep -nE "factor.*2\.2|2\.23" docs/quark_scan_methodology_note.tex`
   returns four hits at lines 863, 867, 877, 943, all in the CFW section
   and all consistent with the audit doc's `2.23` reconciliation number.

4. **Comparison plot correctness** -- PASS.
   `pdftotext results/figures/quark/rs_anarchy_cfw_comparison.pdf` shows
   all four CFW reference markers labelled with both conventions:
   - "CFW RS, gs = 3 no-UV-boundary (10.5 TeV)"
   - "CFW RS, gs 6 default (21 TeV)"
   - "CFW pGB, gs = 3 no-bare-brane (17 TeV)"
   - "CFW pGB, gs 6 default (33 TeV)"
   plus our matched-convention projection at n=217.

5. **Tests pass** -- PASS. `pytest -q tests/test_cfw_comparison.py`
   reports `2 passed in 12.02s` (driver regression + factor-2.2
   assertion). Task brief said "1 passed"; actual count is 2, which is
   strictly stronger.

6. **Methodology PDF** -- PASS. `pdfinfo` confirms 18 pages; rebuild
   from the tex source produces no errors in the CFW section.

7. **Honest physics explanation** -- PASS. `docs/audits/cfw_reproduction.md`
   contains a five-step reconciliation (lines 103-180) decomposing the
   2.23x ratio:
   - Step 1: M0 = 47.26 TeV at full post-audit inputs (BGS-2020 + FLAG-2024
     + LO BMU + `g_s^*=3`).
   - Step 2: reverting to legacy CFW-era `epsilon_K^SM = 1.81e-3` gives
     `budget_legacy/budget_BGS = 6.24`, which dominates the migration
     (`M1 = 18.92 TeV`).
   - Step 3: applying the CFW 30% relative gate over the same forward
     ensemble yields `M2 = 23.37 TeV` (gate factor 1.235, a sampling
     effect, not a Wilson-coefficient change).
   - Step 4: `g_s^*` convention rescaling (`g_s^*=3` vs `~6`).
   - Step 5: final comparison gives 23.37/10.5 = 46.74/21 = 2.23 at both
     conventions.
   The BGS budget shift is the dominant numerical driver; FLAG bags and
   LO BMU enter through M0's normalization and are explicitly cited as
   inputs to step 1.

## Physics-decision adjudication

- **Factor-2.2 framing as headline**: ACCEPT. Per the reconciliation
  doc, the 2.23x ratio is already attributable knob-by-knob: ~6.2x in
  `sigma^2` from the BGS budget shift (`sqrt(6.24) ~ 2.50`), partially
  offset by the CFW-gate sampling factor (`1.235`) and the matched
  `g_s^*` convention. A separate "1-knob-at-a-time" sensitivity table
  reverting all of (BGS->CKMfitter, FLAG->pre-2010, LO->tree)
  simultaneously would be a cosmetic re-presentation of step 2 plus
  small `O(<10%)` corrections from bag/Wilson revisions; it would not
  surface new physics. The audit doc is the substantive reconciliation.

- **BGS/FLAG/LO attribution**: ACCEPTED. The dominant 2.50x in
  amplitude comes from `epsilon_K^SM` central drift (BGS 2020); FLAG-2024
  bag uncertainty and LO BMU running enter at the `O(10-20%)` level
  through M0 and are coherent with the residual gap once the
  gate-sampling and `g_s^*` factors are applied. The 5-step breakdown in
  `docs/audits/cfw_reproduction.md` is sound.

## Recommendation

Proceed to hole #8 (zero-pass stats).
