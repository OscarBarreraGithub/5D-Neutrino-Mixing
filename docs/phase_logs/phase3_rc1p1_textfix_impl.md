# Phase 3 rc1.1 Text-Fix Implementation Report
**Date**: 2026-05-16
**Branch**: paper/quark-scan-2026q2

## Edits applied
WARNING-1:
- Section / line range before: `docs/quark_scan_methodology_note.tex` lines 889-895.
- Diff summary: relabeled the matched CFW comparison as the CFW 30%-relative gate with `n=217` surviving draws, changed the matched number from rounded `23.4` to `23.37` TeV, changed the rescaled value from `47` to `46.74` TeV, and added a footnote separating this matched projection from the live RUNA p50.
- Before: "at common `\gs=3`, we obtain `23.4` TeV (50th percentile, BGS 2020 + FLAG 2024 + factor-3 PDG gate) versus the CFW no-UV-boundary-term result of `10.5` TeV, or equivalently `47` TeV versus `21` TeV at common `\gs\simeq6`."
- After: "at common `\gs=3`, we obtain `23.37` TeV (50th percentile, CFW 30%-relative gate, `n=217` surviving draws) versus the CFW no-UV-boundary-term result of `10.5` TeV, or equivalently `46.74` TeV versus `21` TeV at common `\gs\simeq6`." Footnote added: "`46.74` TeV is the same CFW-matched projection rescaled by `6/3`, distinct from the live RUNA p50 `47.26` TeV at `\gs=3`.

WARNING-2:
- Section / line range before: `docs/quark_scan_methodology_note.tex` lines 970-976.
- Diff summary: separated the live RUNA p50 (`47.26` TeV at `\gs=3`) from the matched-convention factor-2.2 comparison (`23.37` TeV vs `10.5` TeV), and restated the `46.74` vs `21` TeV equivalent at `\gs\simeq6`.
- Before: "Under the central BGS/FLAG inputs and the scope conventions in Appendix A, the p50 value is `47.26` TeV at `\gs=3`; the CFW no-UV-boundary-term RS marker at the same convention is `10.5` TeV, while its default `\gs\simeq6` marker is `21` TeV."
- After: "Under the central BGS/FLAG inputs and the scope conventions in Appendix A, the live RUNA p50 value is `47.26` TeV at `\gs=3` for the unmatched default ensemble. The factor-2.2 CFW framing instead uses the matched-convention comparison, `23.37` TeV versus the CFW no-UV-boundary-term RS marker `10.5` TeV at the same `\gs=3` convention; equivalently, the same matched projection gives `46.74` TeV versus `21` TeV at the CFW default `\gs\simeq6` convention."

WARNING-3:
- Section / line range before: `docs/quark_scan_methodology_note.tex` lines 481-489 and 813-818.
- Diff summary: explicitly labeled the three Y-prior spread statements as different metrics: the broad literature-prior heuristic (`<30%`), the empirical median-scale Run B spread (`<15%`), and the Run B 95%-acceptance crossing spread (`~5%` narrow-to-wide).
- Before: "Our choice of `1.5` sits in the middle of this band; in practice the `\Mkk` bound moves by less than `30%` between half-range `1.0` and half-range `3.0`. An empirical sensitivity check ... confirms that the `\Mkk^{\min}` bound moves by `<15%` across ... priors."
- After: "Our choice of `1.5` sits in the middle of this band; the broad literature-prior heuristic is that the `\Mkk` bound moves by less than `30%` when the half-range is swept from `1.0` to `3.0`. An empirical sensitivity check ... sharpens this using a different metric: the median-scale `\Mkk^{\min}` bound moves by `<15%` across ... priors, while the `95%`-acceptance crossings below have an even smaller narrow-to-wide spread."
- Before: "The full spread between narrow and wide is therefore `~5%` at `95%` acceptance ... This confirms (and slightly improves on) the `<30%` statement made in setup-b..."
- After: "The full spread between narrow and wide is therefore `~5%` at `95%` acceptance ... This is the `95%`-acceptance version of the Run B sensitivity check advertised in setup-b: it sits below both the `<15%` median-scale spread in the same empirical scan and the broader `<30%` literature-prior heuristic."

NIT-1:
- Section / line range before: `docs/quark_scan_methodology_note.tex` lines 516-520.
- Diff summary: changed only the explanatory IR-overlap chain for `f_{u,2}` from `0.117` to `0.16`.
- Before: "`f_u \approx (0.0011,\;0.117,\;1.00)`."
- After: "`f_u \approx (0.0011,\;0.16,\;1.00)`."

NIT-2:
- Section / line range before: `docs/quark_scan_methodology_note.tex` lines 524-528.
- Diff summary: changed only the explanatory up-quark mass-scaling tuple from specific loose values to order estimates consistent with the adjusted overlap chain.
- Before: "`(m_u, m_c, m_t) ... \sim (1 MeV, 500 MeV, 200 GeV)`."
- After: "`(m_u, m_c, m_t) ... \sim (O(1 MeV), O(1 GeV), O(100 GeV))`."

## PDF rebuild
- Command run: `cd docs && rm -f *.aux *.log *.out *.toc *.synctex.gz && pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex && pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex && rm -f *.aux *.log *.out *.toc *.synctex.gz`
- Final page count: 19.
- Any unresolved LaTeX warnings: no unresolved reference or rerun warnings remained on the second pass. Non-fatal warnings still present on the second pass were the existing float-placement warnings, hyperref PDF-string warnings, overfull boxes at lines 954-963, 995-1001, and 1082-1091, and underfull table boxes around lines 1031-1049.

## Checksum refresh
- Old sha256 for `docs/quark_scan_methodology_note.pdf`: `b095faa716c8939496b5086ed78aecd64cb38616d5472badfd3cc1b99eb8abed`.
- New sha256 for `docs/quark_scan_methodology_note.pdf`: `5f544e5d1654f52add06fd8adfe97fed77b968e5e1bea79e119882b1ac898883`.
- `artifacts/checksums.sha256` had no pre-existing methodology-PDF line; added the PDF entry only and left all scan-output checksum lines unchanged.

## Commit SHAs
- text-fix commit: `e0b24e8`.
- checksum commit: this commit (`docs(artifacts): refresh PDF checksum after rc1.1 text fixes`; exact immutable SHA is printed by the implementation runner after commit creation because embedding a commit's own SHA in a tracked file would change that SHA).
- pushed: yes.

## No-touch verification
- Requested check: `git diff --stat HEAD~2..HEAD` should show only `docs/quark_scan_methodology_note.tex`, `docs/quark_scan_methodology_note.pdf`, `artifacts/checksums.sha256`, and `docs/phase_logs/phase3_rc1p1_textfix_impl.md`.
- Workspace note: after the text-fix commit was created, `origin/paper/quark-scan-2026q2` advanced with unrelated commit `ebe8a3c` (`docs(flavor-catalog): planner v0 -- design for discovery-mode catalog subdir`), so the literal `HEAD~2..HEAD` range in this final branch state also includes `docs/phase_logs/flavor_catalog_plan_v0.md`. That file was not created or modified by this implementation.
- Implementation-owned commits/files: `e0b24e8` changes only `docs/quark_scan_methodology_note.tex` and `docs/quark_scan_methodology_note.pdf`; this checksum/report commit changes only `artifacts/checksums.sha256` and `docs/phase_logs/phase3_rc1p1_textfix_impl.md`.
- Confirmed: no physics code touched, no constants touched, no scan outputs touched.

===PHASE_3_RC1P1_TEXTFIX_IMPL_END===
