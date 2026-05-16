# Phase 3 Notebook Re-execution Peer Review

**Reviewer**: independent Codex peer reviewer
**Date**: 2026-05-16
**Range reviewed**: `e2b8349..2121e4e` on `paper/quark-scan-2026q2`

## Verdict
APPROVE

The notebook rerun implementation is consistent with the stated Task B scope: the target commits are present on the paper branch and on origin, the notebooks contain fresh output cells with no recorded execution errors, and the expected post-audit rs-anarchy headline numbers are present in output cells. I did not execute notebooks. I inspected notebook JSON, commit stats, source diffs, output text, and the implementation report. The only non-blocking caveat is range hygiene: one prompt-suggested aggregate diff includes unrelated rc1.1 documentation commits, but the four explicit notebook commits themselves do not alter rc1 artifacts, physics code, constants, scan dispatchers, tests, or scan outputs.

## Item-by-item (1-10)

1. PASS: `git log --oneline -- notebooks/... | head -8` shows, in order, `2121e4e`, `a4ee3be`, `5471ccd`, and `e2b8349`, and `git branch -r --contains` reports `origin/paper/quark-scan-2026q2` for each SHA.
2. PASS: non-empty code-output cell counts are reasonable and nonzero: 8 for `dense_scan_2sigma_vs_1sigma_comparison.ipynb`, 5 for `dense_scan_mkk_constraints_pdg2024.ipynb`, 15 for `pdg_quark_target_fix_verification.ipynb`, and 8 for `rs_anarchy_analysis.ipynb`.
3. PASS: `rs_anarchy_analysis.ipynb` output cells, not merely source, contain `M_KK_min percentiles (TeV, g_s*=3): p5=12.38, p50=47.26, p95=127.13` and the headline `47.26 TeV at g_s*=3`.
4. PASS: `git diff --stat e3a0f1e..2121e4e -- quarkConstraints/ qcd/ flavorConstraints/ neutrinos/ yukawa/ warpConfig/ solvers/ scanParams/ tests/ scan_outputs/` is empty, and each target commit has no stat output for those paths.
5. PASS: `git show 2121e4e --name-only` touches only `notebooks/rs_anarchy_analysis.ipynb`; source changes replace the stale `rs_anarchy_20260507T030811` path with `rs_anarchy_runA_20260515T085316`, redirect `FIG_DIR` to `results/figures/quark/notebook_reruns`, switch the documented stale single-tile histogram to pooled RUNA rescaling, and update headline formatting without adding physics-package computation.
6. PASS: `git diff --stat e3a0f1e..2121e4e -- scan_outputs/` is empty, and `git status --short scan_outputs/` reports no working-tree scan-output changes.
7. PASS: all four notebooks have `kernelspec=python3`, `language=python`, `version=3.12.5`, and the error-output count is 0 for each notebook.
8. PASS: the four notebook commits do not modify rc1 artifacts or dispatchers; per-commit artifact checks are empty and `git diff e2b8349^..2121e4e -- docs/quark_scan_methodology_note.tex docs/quark_scan_methodology_note.pdf artifacts/checksums.sha256 artifacts/quarkscan_paper_rc1_manifest.json scripts/run_rs_anarchy*.py scripts/run_rs_anarchy*.sbatch` is empty, though the broader prompt range `e3a0f1e..2121e4e` is non-empty because it also includes unrelated rc1.1 document commits before the notebook rerun.
9. PASS: every implementation-report key-output string appears in notebook output cells, including the dense-scan `fig2_mkk_bound_*.{pdf,png}` saves, PDG 2024 gate strings, and the rs-anarchy pooled RUNA and headline strings.
10. PASS: `rg -n "0\.5500|0\.917|0\.747|2\.16e-3" notebooks/...` returns no matches, so I found no notebook reference to the superseded pre-audit FLAG values or old epsilon string.

## Findings

WARNING: The aggregate check `git diff e3a0f1e..2121e4e -- docs/quark_scan_methodology_note.tex docs/quark_scan_methodology_note.pdf artifacts/checksums.sha256 ...` is not empty, showing changes to the methodology `.tex`, methodology `.pdf`, and `artifacts/checksums.sha256`. Specific fix: for Task B sign-off, use the four explicit notebook SHAs or the corrected notebook-commit range `e2b8349^..2121e4e`; do not describe `e3a0f1e..2121e4e` as containing only notebook reruns.

NIT: `rs_anarchy_analysis.ipynb` also updates aligned explanatory markdown and a histogram title around the documented histogram logic change. Specific fix: if exact source-line containment is audited again, amend the implementation report to mention these non-executable text/title edits; no code change is required.

## Recommendation
Ready for Opus sign-off

===PHASE_3_NOTEBOOK_RERUN_REVIEW_END===
