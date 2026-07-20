# Phase 3 Notebook Re-execution Opus Sign-Off

**Reviewer**: Opus (independent final sign-off)
**Date**: 2026-05-16
**Range reviewed**: `e2b8349..2121e4e` on `paper/quark-scan-2026q2`

## Verdict
PASS

## Summary
Task B was executed cleanly. All four target notebooks are re-executed under
the post-audit ΔF=2 constants (`quarkConstraints/deltaf2.py:618-623`), every
notebook is error-free, and the rs-anarchy headline (47.26 / 127.13 TeV at
g_s*=3) literally appears in output cells. The notebook-commit range touches
nothing outside `notebooks/`. The methodology PDF (sha256
`5f544e5d…898883`, 19 pages) is byte-identical to the checksum listed in
`artifacts/checksums.sha256:1`, and the rc1.1 numbers (47.26, 127.13, 23.37,
10.5, 22.49) are all present in `docs/quark_scan_methodology_note.tex`.

## Item-by-item (A–J)

A. PASS. `git log --oneline -- notebooks/...` lists `2121e4e`, `a4ee3be`,
   `5471ccd`, `e2b8349` in order, and `git branch -r --contains <sha>` returns
   `origin/paper/quark-scan-2026q2` for each.

B. PASS. `jq` extraction of `notebooks/rs_anarchy_analysis.ipynb` output cells
   yields the literal strings `M_KK_min percentiles (TeV, g_s*=3): p5=12.38,
   p50=47.26, p95=127.13` and `HEADLINE: M_KK at 50% acceptance among
   PDG-passing anarchic NMFV draws is 47.26 TeV at g_s*=3 (16.54 TeV at
   perturbative g_s≈1.05).` — these live in `outputs[].text`, not source.

C. PASS. `jq '[.cells[]? | .outputs[]? | select(.output_type=="error")] |
   length'` returns 0 for all four notebooks.

D. PASS. `git diff e2b8349^..2121e4e -- quarkConstraints/ qcd/
   flavorConstraints/ neutrinos/ yukawa/ warpConfig/ solvers/ scanParams/
   tests/ scan_outputs/` produces 0 lines of output.

E. PASS. `git show --stat 2121e4e` reports a single file changed
   (`notebooks/rs_anarchy_analysis.ipynb`, +119/-116). Diff confirms the
   four documented edits: stale path `rs_anarchy_20260507T030811` →
   `rs_anarchy_runA_20260515T085316`, `FIG_DIR` redirected to
   `results/figures/quark/notebook_reruns`, single-tile histogram replaced
   with pooled RUNA streaming, and a g_s*=3 rescaling + headline cell added.
   Notebook now contains 6 `g_s*=3` strings and 1 `rs_anarchy_runA_*` /
   1 `notebook_reruns` references — consistent with a minimal scoped fix,
   not a rewrite.

F. PASS. `pdfinfo docs/quark_scan_methodology_note.pdf | grep Pages` returns
   `Pages: 19`. PDF sha256 is `5f544e5d1654f52add06fd8adfe97fed77b968e5e1bea
   79e119882b1ac898883`, matching the single recorded entry in
   `artifacts/checksums.sha256`. All five rc1.1 numbers (47.26 × 9, 127.13 × 4,
   23.37 × 3, 10.5 × 7, 22.49 × 2) appear in `docs/quark_scan_methodology_note.tex`.

G. PASS.
   - `dense_scan_2sigma_vs_1sigma_comparison.ipynb`: outputs show
     `2σ accepted points: 83,961`, `1σ accepted: 83,958 / 83,961  (100.0%)`,
     `dropped 2σ → 1σ: 0`.
   - `dense_scan_mkk_constraints_pdg2024.ipynb`: outputs show
     `accepted (publication ξ_KK rescaled): 83961`, `data["M_KK"] range :
     1224 … 48974 GeV`, `epsilon_K     34096 ( 91.4 %)`.
   - `pdg_quark_target_fix_verification.ipynb`: outputs show
     `PDG edition  : PDG 2024` and `REJECTED by new strange gate alone :
     940  (94.0%)`.

H. PASS. Implementation report (`docs/phase_logs/phase3_notebook_rerun_impl.md`)
   explicitly labels three notebooks "CLEAN" (no bug fix needed) and the
   fourth "FIXED-AND-CLEAN", with the four substantive edits listed at
   notebook lines 69, 74, 308–335, 459–468. No physics-code change is claimed
   — and none exists (item D).

I. PASS. The Codex peer-reviewer WARNING is about prompt-range hygiene
   (`e3a0f1e..2121e4e` versus the correct `e2b8349^..2121e4e`), not about
   the notebook commits' content; the corrected range I used here is empty
   for physics/artifact paths. The NIT (markdown title edits in
   `rs_anarchy_analysis.ipynb`) is purely cosmetic and non-load-bearing.

J. PASS. Spot-reading the rs-anarchy notebook output cells shows the
   pooled-draws line (`PDG-passing draws pooled across RUNA tiles: n =
   1,532,640`) immediately followed by the perturbative and g_s*=3
   percentile lines. The numbers reproduce the rc1.1 headline used in the
   methodology .tex. No new physics is introduced — the g_s*=3 rescaling
   is the same convention already documented in the methodology note.

## Findings

WARNING (carried over, non-blocking): the prompt's broader aggregate range
`e3a0f1e..2121e4e` does include unrelated rc1.1 doc commits; the audit-correct
range is `e2b8349^..2121e4e` (used above). Fix: future sign-off prompts
should cite the notebook-only range explicitly. No code action required.

NIT (carried over): `rs_anarchy_analysis.ipynb` commit `2121e4e` also makes
small markdown/title text changes adjacent to the histogram-logic change.
Fix: optional one-line addendum to the impl report mentioning this; no
behavior change.

## Recommendation
Task B is sealed. Ready to tag `quarkscan-paper-rc1.1`.

===PHASE_3_NOTEBOOK_RERUN_OPUS_SIGNOFF_END===
