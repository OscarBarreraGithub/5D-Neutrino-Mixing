# Phase 2 Hole #9 Sign-off

**Verdict**: PASS

## Verification 1-5 with evidence

1. **Branch pushed with the 4 hole-#9 commits.** `git ls-remote origin
   paper/quark-scan-2026q2` returns
   `1a107c8a71d8f047648dcb0f731751d0177f13b1`, matching local HEAD. The four
   hole-#9 commits are present at the tip of the branch in order:
   - `f911a13` audit(scope): inventory scope-wording statements across the
     methodology note
   - `bd8ded1` docs(paper): add Scope and approximations appendix subsection
   - `e359faa` docs(paper): align inline body wording with scope appendix
   - `1a107c8` docs(repo): update README + CLAUDE.md with paper-branch pointer

2. **Methodology note has a "Scope and approximations" subsection.**
   `grep -n "Scope and approximations" docs/quark_scan_methodology_note.tex`
   returns `1121:\subsection{Scope and approximations}` with label
   `\label{app:scope-and-approximations}` on line 1122.

3. **Subsection covers all 4 topics.** Spot-check around the section header
   (lines 1121-1172) confirms four bold paragraph headings:
   - `\textbf{KK tower truncation.}` (line 1124) - first-mode-only matching,
     quotes the +20-30% systematic relative to tower-summed CFW / CPSW.
   - `\textbf{EWPO floor.}` (line 1140) - states the 3 TeV custodial floor is
     external literature (CPSW), not computed here, and notes the quark
     bound is independent.
   - `\textbf{Gauge-coupling convention.}` (line 1154) - declares
     `gs* = 3` as the paper headline, perturbative `~1.05` as code default,
     and `~6` as a CFW-comparison-only convention.
   - `\textbf{Lepton-sector scope.}` (line 1166) - explicitly excludes
     `neutrinos/`, `flavorConstraints/`, `yukawa/` from the current paper.

4. **README and CLAUDE.md reference the paper branch + quark-only scope +
   methodology note location.** `README.md` lines 46-49 read
   "Paper branch note: the active methodology note is
   `docs/quark_scan_methodology_note.tex` on `paper/quark-scan-2026q2`. The
   current paper scope is quark-sector only; lepton-sector packages remain
   repo tools for follow-up work." `CLAUDE.md` lines 12-15 add the matching
   agent-facing pointer and quark-only scope statement.

5. **Methodology PDF builds, page count noted.** `pdfinfo
   docs/quark_scan_methodology_note.pdf` reports `Pages: 19`, mtime
   2026-05-15 09:12, file size 595742 bytes - PDF rebuilt after the
   subsection was added.

## Recommendation

Proceed to hole #2 (canonical run manifest).
