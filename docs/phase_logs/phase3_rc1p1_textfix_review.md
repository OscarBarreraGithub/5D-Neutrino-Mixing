# Phase 3 rc1.1 Text-Fix Peer Review

**Reviewer**: independent Codex peer reviewer
**Date**: 2026-05-16
**Range reviewed**: `e3a0f1e..8deb291` on `paper/quark-scan-2026q2`

## Verdict
APPROVE

## Item-by-item (1-12)
1. **PASS** - `git show e0b24e8 -- docs/quark_scan_methodology_note.tex` adds the corrected matched-convention wording, and the current Headline paragraph now says `23.37` TeV with `CFW 30%-relative gate, n=217 surviving draws`, gives `46.74` TeV at `\gs\simeq6`, and footnotes that this is distinct from live RUNA p50 `47.26` TeV (`docs/quark_scan_methodology_note.tex:899`, `:902`, `:903`).
2. **PASS** - the recommendation now separates the live unmatched RUNA p50 `47.26` TeV from the factor-`2.2` matched comparison, explicitly quoting `23.37` TeV versus `10.5` TeV and `46.74` versus `21` TeV (`docs/quark_scan_methodology_note.tex:984`, `:986`, `:987`, `:989`).
3. **PASS** - the Y-prior language now labels the `30%` statement as a broad literature-prior heuristic, the `15%` statement as an empirical median-scale Run B metric, and the `5%` statement as the 95%-acceptance narrow-to-wide crossing spread (`docs/quark_scan_methodology_note.tex:484`, `:488`, `:819`, `:821`).
4. **PASS** - optional NIT-1 was tightened: `f_{u,2}` is now `0.16`, not retained as `0.117` (`docs/quark_scan_methodology_note.tex:522`; implementation report `docs/phase_logs/phase3_rc1p1_textfix_impl.md:26`).
5. **PASS** - optional NIT-2 was tightened by replacing the specific `(1 MeV, 500 MeV, 200 GeV)` scaling tuple with order estimates `O(1 MeV)`, `O(1 GeV)`, and `O(100 GeV)` (`docs/quark_scan_methodology_note.tex:530`).
6. **PASS** - the canonical load-bearing TeX numbers still grep in the new file: `47.26`/`127.13` (`docs/quark_scan_methodology_note.tex:587`), `23.37` (`:899`), `10.5` (`:901`), `22.49` (`:1203`), `99.50` (`:1214`), `0.5503`/`0.903`/`0.691`/`2.161` (`:1033`, `:1035`, `:1037`, `:1040`); requested string `543` is absent from both current and `e3a0f1e` TeX, so no rc1.1 removal occurred.
7. **PASS** - `git diff e3a0f1e..8deb291 -- quarkConstraints/ qcd/ flavorConstraints/ neutrinos/ yukawa/ warpConfig/ solvers/ scanParams/ tests/ | wc -c` returned `0`.
8. **PASS** - `git diff e3a0f1e..8deb291 -- scan_outputs/ results/figures/ | wc -c` returned `0`.
9. **PASS** - `pdfinfo docs/quark_scan_methodology_note.pdf | grep -i Pages` returned `Pages:          19`.
10. **PASS** - `sha256sum docs/quark_scan_methodology_note.pdf` returned `5f544e5d1654f52add06fd8adfe97fed77b968e5e1bea79e119882b1ac898883`, matching `artifacts/checksums.sha256:21` for the PDF.
11. **PASS** - `git show --stat e0b24e8 8deb291 | grep -E '\.(aux|log|out|toc|synctex\.gz|bbl|blg)'` produced no output.
12. **PASS** - the implementation report accurately describes WARNING-1/2/3 and NIT-1/2 (`docs/phase_logs/phase3_rc1p1_textfix_impl.md:6`, `:12`, `:18`, `:26`, `:32`) and honestly notes the unrelated `flavor_catalog_plan_v0.md` commit present in the wider reviewed range (`docs/phase_logs/phase3_rc1p1_textfix_impl.md:53`).

## Findings
No findings. The `543` grep target appears to belong to phase acceptance/test logs rather than `docs/quark_scan_methodology_note.tex`; it was absent before and after this implementation range and is not a text-fix regression. Later local notebook commits are outside the requested `e3a0f1e..8deb291` review range and do not alter the paper TeX, PDF, or checksum evidence reviewed here.

## Recommendation
Ready for Opus re-review and rc1.1 tag

===PHASE_3_RC1P1_TEXTFIX_REVIEW_END===
