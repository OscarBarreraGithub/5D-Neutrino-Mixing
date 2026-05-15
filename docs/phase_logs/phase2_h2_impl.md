# Phase 2 Hole #2 Implementation: Canonical Run-Manifest Blessing

Date: 2026-05-15  
Branch: `docs/canonical-manifest`  
Target branch: `paper/quark-scan-2026q2`  
Run-time code SHA: `29803fffa70ede4ffab736b4bebc705eeceef0f6` (inferred from `20260515T0853*` UTC scan stamps after the Wilson-RG fix; `tile_summary.json` does not embed git metadata)  
Manifest: `docs/artifact_manifest.md`  
Machine-readable manifest: `artifacts/quarkscan_paper_rc1_manifest.json`  
Checksums: `artifacts/checksums.sha256`

## Canonical scans and checksums

| Scan | Canonical dir | Used for | tile_summary SHA256 | draws SHA256 | Size MB |
|---|---|---|---|---|---:|
| Baseline 8M (RUNA) | `scan_outputs/rs_anarchy_runA_20260515T085316` | All headline p50/p95 quotes | `b146fefc0fba7b1fdb0fa488b19c8ccb04066eed860bb02a59a40fabf37bd3ce` | `a2d43049d161e041ad9b790e2c95777e72ff0e1af9480dffea4454905d0f8640` | 8208.01 |
| Run3 baseline | `scan_outputs/rs_anarchy_run3_baseline_20260515T085324` | c-pattern figure | `246cf3dfb073aa12e56a27aa3b77a03917623f6c1703f4af91700c97c6c21e24` | `d9dabdcb9543d211f88713e77bd9a65a49145ed3386a9102e2f8721134f4f864` | 1640.90 |
| Run3 qtop_shifted | `scan_outputs/rs_anarchy_run3_qtop_shifted_20260515T085324` | c-pattern figure | `ef0fbc2b8238fc57f1e697027db5c790070c460a671e612cd0723d6e6bc25f58` | `badf0a538efdcdb1620794fd37393e138ff803e84c2e16aa299feb7250af083d` | 1640.96 |
| Run3 moreUV | `scan_outputs/rs_anarchy_run3_moreUV_20260515T085324` | c-pattern figure + zero-pass UL | `3150aadc477183574d07d0be227edc3d646a08c81ce3e668fd5ef8daae3158cd` | `d199744d5a92cef521051c911b3e750fd48dd9ec6bf431a8d819e0a0d46f79cb` | 1660.10 |
| Run3 moreIR | `scan_outputs/rs_anarchy_run3_moreIR_20260515T085324` | c-pattern figure + zero-pass UL | `e6899f1af63821ab7afcd6389126d6522eff2f626ac154541c3e679890aa30d5` | `66ba3a70c5be0f2779bb95853e2c326a9b754f9ae367804e68fb0babbb91b0de` | 1619.43 |
| RunB narrow_uniform | `scan_outputs/rs_anarchy_runB_narrow_uniform_20260515T085324` | Y-prior figure | `17a0a27ddcc494677720077527dd6a263f2fbfc3b16711e31d1e00aadefb4edb` | `f73c976a3c5aee29458594b423cd7cc65fefdda6fc33c8381cecc093eed36e83` | 1641.79 |
| RunB wide_uniform | `scan_outputs/rs_anarchy_runB_wide_uniform_20260515T085325` | Y-prior figure | `9ecbe08487ba4aedb2e862f058311e9c81a85fc7bedcd42508ffb38bd62d795a` | `90500218a8b6c6541508a43814671a307921ecf51de255fabbb4be75bf9931c1` | 1638.95 |
| RunB gaussian_3sigma | `scan_outputs/rs_anarchy_runB_gaussian_3sigma_20260515T085324` | Y-prior figure | `6d180ca1e5417fe3c095ea2447058b003d26abf18dbe731f4439d4c4a5d45695` | `aeadd9e1cfba4a8115f3582178a4b964041bb09152812ab72234dd4ebb2f614e` | 1641.79 |
| RunC CFW gate | `scan_outputs/rs_anarchy_runC_20260515T085323` | CFW comparison + zero-pass UL | `92f57214f1e8164eac1dfffeec4c8e26154520d8f8543784f7ea5cead460cf3d` | `3f8dee6037c51a344e77cd0992fcb3c51b4c06032c9052bc331679e68384b53f` | 4098.97 |
| Old 800k baseline | `scan_outputs/rs_anarchy_20260507T030811` | Historical cross-check only; NOT used for paper claims | `28d94014acda9663afbfd87578750be94132110cb38b5d1a626a59684d877b11` | `ffe12d1225309d04275119353e2b345e8b31daa4cdf08e8ac73c2e44d8a6d11a` | 803.57 |

## Paper claim mapping

| Claim location | Quoted number / claim | Scan dir | Plot file |
|---|---|---|---|
| `docs/quark_scan_methodology_note.tex:582` | `47.26 TeV` at p50, `g_s*=3` | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.pdf` |
| `docs/quark_scan_methodology_note.tex:583` | `127.13 TeV` at p95, `g_s*=3` | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.pdf` |
| `docs/quark_scan_methodology_note.tex:649` | Headline table `16.54/47.26/44.50/127.13 TeV` for `ratio < 1` | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_gate_sensitivity.pdf` |
| `docs/quark_scan_methodology_note.tex:692` | RUNA `8M` draws, `1,532,640` PDG-passing; p50/p95 Wilson MC intervals | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.pdf` |
| `docs/quark_scan_methodology_note.tex:718` | BGS-budget band `47.26^{+69.37}_{-24.98} TeV` | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.pdf` |
| `docs/quark_scan_methodology_note.tex:736` | Run3 c-pattern figure; qtop shifted `45.73/120.39 TeV`; moreUV/moreIR `p_UL=2.3e-6` | Run3 baseline, qtop_shifted, moreUV, moreIR | `results/figures/quark/rs_anarchy_mkk_min_hist_by_cvals.pdf` |
| `docs/quark_scan_methodology_note.tex:806` | RunB 95%-acceptance crossings `46.32/48.77/48.84 TeV` perturbative and `132.35/139.34/139.55 TeV` at `g_s*=3` | RunB narrow_uniform, wide_uniform, gaussian_3sigma | `results/figures/quark/rs_anarchy_mkk_min_hist_by_yprior_gsstar.pdf` |
| `docs/quark_scan_methodology_note.tex:836` | PDG-gate counts `1,532,640`, `110,519`, `1,164` | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_mkk_min_hist_by_pdg_tightness.pdf` |
| `docs/quark_scan_methodology_note.tex:854` | RunC zero-pass finite-ensemble UL `p_pass <= 9.2e-7` | `rs_anarchy_runC_20260515T085323` | `results/figures/quark/rs_anarchy_cfw_comparison.pdf` |
| `docs/quark_scan_methodology_note.tex:889` | CFW comparison: factor `2.2` stronger; matched projection `23.37 TeV`; `n=217` | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_cfw_comparison.pdf` |

## Verification

- `python -m json.tool artifacts/quarkscan_paper_rc1_manifest.json` passed.
- `sha256sum -c artifacts/checksums.sha256` passed for the canonical scans and the historical 800k cross-check.
- `pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex` was run twice after adding the methodology-note pointer; output is `docs/quark_scan_methodology_note.pdf` with 19 pages.
- Local cluster bundle created at `artifacts/quarkscan_paper_rc1_manifest_bundle.tar.gz`; SHA256 sidecar records `5402d2ac035801b3a64711e22332ced7cf9bc7e926c0d4414f2a996a9b4b27b7`.

Hole #2 ready for peer review.
