# Phase 3 rc1.1 Text-Fix Opus Re-Review Sign-Off

**Reviewer**: Opus (independent re-review)
**Date**: 2026-05-16
**Range reviewed**: `e3a0f1e..8deb291` on `paper/quark-scan-2026q2`

## Verdict
PASS

## Summary
I independently re-ran every check (A-K) without trusting the Codex peer review.
The diff `git diff e3a0f1e..8deb291 -- docs/quark_scan_methodology_note.tex`
matches the implementation report exactly: the three WARNING-class fixes are
in place at the right line numbers, the two NIT cleanups landed, no
load-bearing numbers were dropped, no physics code or scan outputs were
touched, the PDF rebuilt at 19 pages, the new sha256 matches
`artifacts/checksums.sha256`, and no aux/log/bbl detritus was committed.
The rewordings are cosmetic and physically faithful; no new claim is
introduced that was not in rc1.

## Item-by-item (A-K)

A. **PASS** — Headline paragraph at `docs/quark_scan_methodology_note.tex:896-917`.
   Line 899 reads `$23.37$~TeV`; line 900 reads `CFW 30\%-relative gate, $n=217$
   surviving draws`; line 902 reads `$46.74$~TeV versus $21$~TeV`; the footnote
   at lines 902-905 explicitly separates the matched rescale `46.74` TeV from
   the live RUNA p50 `47.26` TeV (`unmatched post-audit factor-of-three gate
   ensemble`). Wording, numerics, and footnote all match the prior sign-off
   prescription.

B. **PASS** — Recommendation at `docs/quark_scan_methodology_note.tex:983-993`.
   Live RUNA p50 `47.26` TeV is now explicitly labeled `for the unmatched
   default ensemble` (line 985), and the factor-2.2 framing immediately
   re-quotes `23.37` TeV versus `10.5` TeV (lines 986-988) and `46.74` versus
   `21` TeV at `\gs\simeq6` (lines 988-990). A reader can no longer compute
   `47.26 / 10.5 = 4.5` and infer inconsistency.

C. **PASS** — Y-prior labeling. `tex:484-492` now labels the `<30\%` figure
   as a `broad literature-prior heuristic` (line 484-485), the `\lesssim 15\%`
   as the empirical `median-scale $\Mkk^{\min}$ bound` (line 488), and
   foreshadows that `95\%`-acceptance crossings have an even smaller
   narrow-to-wide spread (lines 491-492). `tex:819-824` closes the loop:
   the `\approx 5\%` figure is explicitly the `$95\%$-acceptance version`
   of the Run-B check, sitting below both the `\lesssim 15\%` median-scale
   and the `\le 30\%` literature heuristic. The three metrics are now
   disambiguated.

D. **PASS** — NIT-1 applied at `tex:522`: `f_u \approx (0.0011,\;0.16,\;1.00)`,
   replacing `0.117`. At warp log $L \approx 37$, $1/\sqrt{L} \approx 0.164$,
   so 0.16 is the c=1/2 limit one obtains directly without the
   convention-dependent $\sqrt{2}$ factor — internally consistent with the
   rest of the explanatory chain in this paragraph.

E. **PASS** — NIT-2 applied at `tex:530-532`: replaced
   `(1\,\mathrm{MeV},\,500\,\mathrm{MeV},\,200\,\mathrm{GeV})` with
   `\bigl(\mathcal{O}(1\,\mathrm{MeV}),\,\mathcal{O}(1\,\mathrm{GeV}),\,
   \mathcal{O}(100\,\mathrm{GeV})\bigr)`. The order-of-magnitude framing
   removes the loose factor-of-two/three slack the prior sign-off flagged.

F. **PASS** — All load-bearing numbers still present in current .tex:
   `47.26` (9 occurrences), `127.13` (4), `23.37` (3), `10.5` (7), `22.49` (2,
   at `tex:1203,1208`), `99.50` (1, at `tex:1214`), `0.5503` (1), `0.903` (1),
   `0.691` (1), `2.161` (1). `46.74` is newly introduced (3 occurrences),
   consistent with the WARNING-1/2 fixes.

G. **PASS** — `git diff e3a0f1e..8deb291 -- quarkConstraints/ qcd/
   flavorConstraints/ neutrinos/ yukawa/ warpConfig/ solvers/ scanParams/
   tests/ | wc -l` returns `0`. No physics code touched.

H. **PASS** — `git diff e3a0f1e..8deb291 -- scan_outputs/ results/figures/
   | wc -l` returns `0`. No scan outputs or figures touched.

I. **PASS** — `pdfinfo docs/quark_scan_methodology_note.pdf | grep Pages`
   returns `Pages: 19`. `sha256sum docs/quark_scan_methodology_note.pdf`
   returns `5f544e5d1654f52add06fd8adfe97fed77b968e5e1bea79e119882b1ac898883`,
   which matches `artifacts/checksums.sha256:21`.

J. **PASS** — `git show --stat e0b24e8 8deb291 | grep -E
   '\.(aux|log|out|toc|synctex\.gz|bbl|blg)'` is empty. Commits e0b24e8 and
   8deb291 touch only `.tex`, `.pdf`, `checksums.sha256`, and the
   implementation report `.md`. No LaTeX detritus committed.

K. **PASS** — Read each rewritten paragraph end-to-end. The new claims that
   were not in rc1 are: (i) the footnote distinguishing 46.74 TeV from 47.26
   TeV (well-supported by the existing audit `docs/audits/cfw_reproduction.md`);
   (ii) the explicit `n=217` figure for surviving matched-gate draws (already
   in audit and prior sign-off "Specifically validated"); (iii) the Wilson
   confidence-interval style framing of the three Y-prior metrics. None of
   these introduce unsupported physics. No new load-bearing numerical
   claim appears; the rewordings only re-quote numbers already established
   by the audit chain.

## Findings

No findings. The wider diff range `e3a0f1e..8deb291` also contains the
unrelated commit `ebe8a3c` (`flavor_catalog_plan_v0.md` planner doc), which
is honestly disclosed in the implementation report at lines 53-56 and is
out of scope for this re-review — it does not touch the paper TeX, PDF,
checksum, or any physics surface.

## Recommendation

Ready to tag `quarkscan-paper-rc1.1`.

The three WARNING-class presentational issues that the prior end-to-end
Opus read flagged are resolved with the exact wording prescribed, the two
NIT cleanups are in, the rebuild is clean at 19 pages, the checksum is
registered, no physics or scan surfaces were perturbed, and no LaTeX
build artifacts leaked into the commit. The PI's acceptance criteria
(no physics errors, no missing systematics, internally consistent
numbers, no reader-misleading framing) are now fully met.

===PHASE_3_RC1P1_TEXTFIX_OPUS_SIGNOFF_END===
