# Phase 3 Final Acceptance Log

**Date**: 2026-05-15
**Branch**: `paper/quark-scan-2026q2`
**Target tag**: `quarkscan-paper-rc1`

## 10-item checklist

1. **PASS** - `paper/quark-scan-2026q2` contains the tracked code, docs, notebooks, artifacts, and whitelisted figures needed for the quark scan. Evidence: branch `paper/quark-scan-2026q2` at pre-acceptance `e3c0c79e545ded087e08036fb3ce69675997b47d`; `git ls-files | wc -l` returned `378`. Per-extension tracked counts: `py 151`, `tex 4`, `md 68`, `json 24`, `sha256 2`, `tar.gz 1`, `ipynb 21`, `sbatch 10`; the requested combined grep returned `281`.

2. **PASS** - no xfailed tests remain in `tests/test_quark_fit.py`, and the old benchmark surface is preserved in regression fixtures. Evidence: `/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python -m pytest -q tests/test_quark_fit.py` returned `14 passed in 5.12s`. `ls tests/baselines/` returned `pre-fix-xfail-output.txt` and `post-fix-test-output.txt`.

3. **PASS** - bag parameters and Delta F=2 RG conventions are source-traceable in the methodology note. Evidence: `grep -nE "FLAG 2024|BGS 2020|BMU|hep-ph/0102316|2411\.04268|1911\.06822" docs/quark_scan_methodology_note.tex` returned hits for every requested source, including FLAG/BGS/BMU at lines 893, 896, 898; FLAG arXiv:2411.04268 at line 1017; BGS arXiv:1911.06822 at line 1024; and BMU arXiv:hep-ph/0102316 at line 1070.

4. **PASS** - bag/RG audits tripped the invalidation gate, the rerun completed, and post-audit headlines are live. Evidence: `phase2_h5_signoff.md` records `INVALIDATION_GATE_TRIPPED`; `phase2_h6_signoff.md` confirms the cumulative `22.49x` invalidation factor and mandatory reruns; `invalidation_gate_signoff.md` verifies the scan rerun, figure rerun, and methodology headline update. `grep -nE "47\.26|127\.13|22\.49|22\.5" docs/quark_scan_methodology_note.tex` returned live hits at lines 582-583, 649, 698-700, 706, 719, 920-922, 974, 1157, 1188, and 1193.

5. **PASS** - CFW comparison is quantitatively narrowed to the factor-2.2 statement, not the retired 11% agreement framing. Evidence: `grep -nE "factor.*2\.2" docs/quark_scan_methodology_note.tex` returned hits at lines 891, 895, 905, and 971.

6. **PASS** - `g_s^*`, EWPO floor, KK tower truncation, and lepton-sector scope are unambiguous. Evidence: `grep -nE "Scope and approximations|KK tower|EWPO|lepton" docs/quark_scan_methodology_note.tex` returned the cross-reference at lines 93-94 and the appendix subsection at line 1128, with KK tower, EWPO, and lepton-sector entries at lines 1131, 1147, and 1174-1175.

7. **PASS** - prior/gate robustness and zero-pass statements use finite-ensemble and Wilson-score upper-limit wording. Evidence: `grep -nE "p_UL|Wilson.*upper|N=.*1\.6|finite-ensemble" docs/quark_scan_methodology_note.tex` returned the finite-ensemble caveat at line 497, `N=1,600,000` Wilson upper-limit statements at lines 743 and 761, additional finite-ensemble wording at lines 771 and 775, and the Run C Wilson upper limit `p_{\rm pass}\le 9.2\times10^{-7}` at lines 858-859.

8. **PASS** - a tracked final run manifest maps paper claims to code SHA, command, output, checksum, and figure hash. Evidence: `git ls-files | grep artifact_manifest` returned `docs/artifact_manifest.md`. `git ls-files artifacts/` returned `artifacts/README.md`, `artifacts/checksums.sha256`, `artifacts/quarkscan_paper_rc1_manifest.json`, `artifacts/quarkscan_paper_rc1_manifest_bundle.tar.gz`, and `artifacts/quarkscan_paper_rc1_manifest_bundle.tar.gz.sha256`. Spot check: `sha256sum -c artifacts/checksums.sha256 2>&1 | head -5` returned five `OK` lines, through `scan_outputs/rs_anarchy_run3_qtop_shifted_20260515T085324/tile_summary.json: OK`.

9. **PASS** - submission figures are pruned, named, and regenerable from tracked scripts. Evidence: `ls results/figures/quark/*.pdf | wc -l` returned `9`; `ls results/figures/quark/exploratory/ | wc -l` returned `40`. `docs/audits/figure_prune_inventory.md` maps all nine cited PDFs to tracked regeneration scripts or commands at lines 44-52.

10. **PASS** - clean-tree tests, PDF builds, release-tag preconditions, and external scan-output snapshot are verified. Evidence: full suite `/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python -m pytest -q` returned `543 passed, 1 skipped in 960.06s (0:16:00)` with no failures and no xfails. The methodology PDF was rebuilt from `docs/` by removing aux/log/out/toc/synctex files, running `pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex` twice, and removing aux files again; `pdfinfo docs/quark_scan_methodology_note.pdf | grep Pages` returned `Pages:          19`. The compact assumptions PDF was also rebuilt twice per roadmap Phase 3; `pdfinfo docs/quark_scan_assumptions_compact.pdf | grep Pages` returned `Pages:          7`. The artifact tarball and checksums are tracked under `artifacts/`. The release tag is created and pushed after this acceptance commit, with remote confirmation recorded in `phase3_final_summary.md`.

## Headline numbers (locked at rc1)

| Quantity | Value | Convention |
|---|---|---|
| M_KK^min p50 (perturbative g_s) | 16.54 TeV | g_s ≈ 1.05 |
| M_KK^min p50 (g_s*=3 headline) | 47.26 +69.4 −25.0 TeV | g_s*=3, BGS 2020 + LO + factor-3 |
| M_KK^min p95 (g_s*=3) | 127.13 TeV | same |
| CFW reconciliation | factor 2.2 stronger at matched g_s*=3 | -- |
| Zero-pass UL (moreUV/moreIR) | p ≤ 2.3×10⁻⁶ | 95% Wilson |
| Zero-pass UL (Run C factor-1.5) | p ≤ 9.2×10⁻⁷ | 95% Wilson |
| ε_K binding fraction | 99.5% | post-audit |
| Total commits on paper branch | 154 |
| Methodology note pages | 19 |
| Test suite | 543 passed, 1 skipped, 0 xfails |

## Deferred follow-up items (not blockers)

1. Endpoint mismatch (B-system bags at m_b vs Wilson endpoint at 2 GeV) → 10–30% systematic, documented in appendix
2. LO-only Wilson running (NLO not implemented) → 10–30% systematic, documented
3. VLL/VRR NLO sanity check
4. Quantitative BMU/BBL convention table
5. Run C high-statistics rerun (current n=217 at CFW gate)

## Sign-off chain (commit SHAs)

- `docs/phase_logs/phase2_h2_signoff.md` - `80641cf`, `cb27631`, `b86e66d`, `2729659`, `b946abf`, `5e60bc0`
- `docs/phase_logs/phase2_h4_signoff.md` - `16e6b36`, `2871fce`, `559b851`, `fd83d96`, `7e679ea`
- `docs/phase_logs/phase2_h5_signoff.md` - `dc9c498`, `82a96f0`, `695f35e`
- `docs/phase_logs/phase2_h6_signoff.md` - `7f71908`, `561d848`, `c80a8e0`
- `docs/phase_logs/invalidation_gate_signoff.md` - `29803ff`, `db02223`, `87c728f`, `38a4586`, `217af80`
- `docs/phase_logs/phase2_h7_signoff.md` - `ffec986`, `06d5d85`, `da8f647`, `508ae69`, `b47aa76`, `330ffa9`, `11e0d58`, `321fe0f`, `7c284a8`, `2a3a2c0`, `f2c29c8`
- `docs/phase_logs/phase2_h8_signoff.md` - `f34f9df`, `a9d1058`, `212320f`, `1d9e17e`
- `docs/phase_logs/phase2_h9_signoff.md` - `f911a13`, `bd8ded1`, `e359faa`, `1a107c8`
- `docs/phase_logs/phase3_h10_signoff.md` - `e7d824d`, `3951f41`, `bf6186c`, `00b222d`
