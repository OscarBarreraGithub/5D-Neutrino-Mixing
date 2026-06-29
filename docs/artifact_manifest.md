# Canonical Artifact Manifest

> **LEGACY (May 2026 ΔF=2-only / anarchy artifacts, pre-audit).** This manifest
> seals the `quarkscan-paper-rc1` ΔF=2 scan artifacts. It predates the June 2026
> audit and does not reflect the corrected minimal-RS floors (typical ~30 TeV
> from `epsilon_K`, existence ~18-20 TeV from oblique S,T,U, Z→bb ~5 TeV). For
> the current picture see `docs/FLOOR_SUMMARY.md` and
> `reports/collaborator_2026-06/CONTENT.md`.

**Tag**: `quarkscan-paper-rc1` (to be created only at Phase 3 sign-off)  
**Date sealed**: 2026-05-15  
**Branch**: `paper/quark-scan-2026q2`  
**Implementation branch**: `docs/canonical-manifest`  
**Code SHA at run time**: `29803fffa70ede4ffab736b4bebc705eeceef0f6` (inferred post-Wilson-RG-fix SHA; scan dirs were stamped `20260515T0853*` UTC after this commit and before rerun-record commit `db0222366cf53184ed90e7009739005580facec9`; no git SHA is embedded in `tile_summary.json`)  
**Code SHA at manifest seal**: `80641cf653a3a304afb40b65ffa0287222eea945` (docs-only hygiene parent for this manifest; no physics code changes after the run-time SHA)  
**Bundle path**: `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/artifacts/quarkscan_paper_rc1_manifest_bundle.tar.gz` with sidecar SHA256 at the same path plus `.sha256`  
**DOI**: deferred until submission.

## Canonical scan dirs (cluster paths)

| Scan | Path | tile_summary SHA256 | draws SHA256 | Size MB | Run command |
|---|---|---|---|---:|---|
| Baseline 8M (RUNA) | `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scan_outputs/rs_anarchy_runA_20260515T085316` | `b146fefc0fba7b1fdb0fa488b19c8ccb04066eed860bb02a59a40fabf37bd3ce` | `a2d43049d161e041ad9b790e2c95777e72ff0e1af9480dffea4454905d0f8640` | 8208.01 | `sbatch scripts/run_rs_anarchy_runA.sbatch` |
| Run3 baseline | `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scan_outputs/rs_anarchy_run3_baseline_20260515T085324` | `246cf3dfb073aa12e56a27aa3b77a03917623f6c1703f4af91700c97c6c21e24` | `d9dabdcb9543d211f88713e77bd9a65a49145ed3386a9102e2f8721134f4f864` | 1640.90 | `sbatch --export=ALL,PATTERN_TAG=baseline,C_Q="0.63_0.57_0.20",C_U="0.66_0.50_-0.50",C_D="0.66_0.61_0.55" scripts/run_rs_anarchy_run3.sbatch` |
| Run3 qtop_shifted | `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scan_outputs/rs_anarchy_run3_qtop_shifted_20260515T085324` | `ef0fbc2b8238fc57f1e697027db5c790070c460a671e612cd0723d6e6bc25f58` | `badf0a538efdcdb1620794fd37393e138ff803e84c2e16aa299feb7250af083d` | 1640.96 | `sbatch --export=ALL,PATTERN_TAG=qtop_shifted,C_Q="0.63_0.57_0.30",C_U="0.66_0.50_-0.55",C_D="0.66_0.61_0.55" scripts/run_rs_anarchy_run3.sbatch` |
| Run3 moreUV | `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scan_outputs/rs_anarchy_run3_moreUV_20260515T085324` | `3150aadc477183574d07d0be227edc3d646a08c81ce3e668fd5ef8daae3158cd` | `d199744d5a92cef521051c911b3e750fd48dd9ec6bf431a8d819e0a0d46f79cb` | 1660.10 | `sbatch --export=ALL,PATTERN_TAG=moreUV,C_Q="0.68_0.62_0.25",C_U="0.71_0.55_-0.45",C_D="0.71_0.66_0.60" scripts/run_rs_anarchy_run3.sbatch` |
| Run3 moreIR | `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scan_outputs/rs_anarchy_run3_moreIR_20260515T085324` | `e6899f1af63821ab7afcd6389126d6522eff2f626ac154541c3e679890aa30d5` | `66ba3a70c5be0f2779bb95853e2c326a9b754f9ae367804e68fb0babbb91b0de` | 1619.43 | `sbatch --export=ALL,PATTERN_TAG=moreIR,C_Q="0.58_0.52_0.15",C_U="0.61_0.45_-0.55",C_D="0.61_0.56_0.50" scripts/run_rs_anarchy_run3.sbatch` |
| RunB narrow_uniform | `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scan_outputs/rs_anarchy_runB_narrow_uniform_20260515T085324` | `17a0a27ddcc494677720077527dd6a263f2fbfc3b16711e31d1e00aadefb4edb` | `f73c976a3c5aee29458594b423cd7cc65fefdda6fc33c8381cecc093eed36e83` | 1641.79 | `sbatch --export=ALL,PRIOR_TAG=narrow_uniform scripts/run_rs_anarchy_runB.sbatch` |
| RunB wide_uniform | `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scan_outputs/rs_anarchy_runB_wide_uniform_20260515T085325` | `9ecbe08487ba4aedb2e862f058311e9c81a85fc7bedcd42508ffb38bd62d795a` | `90500218a8b6c6541508a43814671a307921ecf51de255fabbb4be75bf9931c1` | 1638.95 | `sbatch --export=ALL,PRIOR_TAG=wide_uniform scripts/run_rs_anarchy_runB.sbatch` |
| RunB gaussian_3sigma | `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scan_outputs/rs_anarchy_runB_gaussian_3sigma_20260515T085324` | `6d180ca1e5417fe3c095ea2447058b003d26abf18dbe731f4439d4c4a5d45695` | `aeadd9e1cfba4a8115f3582178a4b964041bb09152812ab72234dd4ebb2f614e` | 1641.79 | `sbatch --export=ALL,PRIOR_TAG=gaussian_3sigma scripts/run_rs_anarchy_runB.sbatch` |
| RunC CFW gate | `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scan_outputs/rs_anarchy_runC_20260515T085323` | `92f57214f1e8164eac1dfffeec4c8e26154520d8f8543784f7ea5cead460cf3d` | `3f8dee6037c51a344e77cd0992fcb3c51b4c06032c9052bc331679e68384b53f` | 4098.97 | `sbatch scripts/run_rs_anarchy_runC.sbatch` |

Historical cross-check only: `scan_outputs/rs_anarchy_20260507T030811` is the old 800k pre-audit baseline and is **not used for paper claims**. Its hashes are `tile_summary.json` = `28d94014acda9663afbfd87578750be94132110cb38b5d1a626a59684d877b11`; `draws.jsonl` = `ffe12d1225309d04275119353e2b345e8b31daa4cdf08e8ac73c2e44d8a6d11a`; total size = 803.57 MB.

## Mapping: paper claim -> scan dir

| Claim location (file:line) | Quoted number | Scan dir | Plot file |
|---|---|---|---|
| `docs/quark_scan_methodology_note.tex:582` | `47.26 TeV` at p50, `g_s*=3` | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.pdf` (`153594583823fb3f4d20ec46484194464a0691593d520036dcf35171ef61fa4a`) |
| `docs/quark_scan_methodology_note.tex:583` | `127.13 TeV` at p95, `g_s*=3` | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.pdf` (`153594583823fb3f4d20ec46484194464a0691593d520036dcf35171ef61fa4a`) |
| `docs/quark_scan_methodology_note.tex:649` | Table row `16.54/47.26/44.50/127.13 TeV` for `ratio < 1` | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_gate_sensitivity.pdf` (`e2782e2a77f09b7929f9e36d66e8998b3cff8063e297e34b8505ee90f3915d91`) |
| `docs/quark_scan_methodology_note.tex:692` | RUNA `8M` draws, `1,532,640` PDG-passing; p50/p95 and Wilson MC intervals | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.pdf` (`153594583823fb3f4d20ec46484194464a0691593d520036dcf35171ef61fa4a`) |
| `docs/quark_scan_methodology_note.tex:718` | BGS-budget band `47.26^{+69.37}_{-24.98} TeV` | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_mkk_min_hist_gsstar.pdf` (`153594583823fb3f4d20ec46484194464a0691593d520036dcf35171ef61fa4a`) |
| `docs/quark_scan_methodology_note.tex:736` | Run3 c-pattern figure; qtop shifted `45.73/120.39 TeV`; moreUV/moreIR zero-pass UL `2.3e-6` | `rs_anarchy_run3_baseline_20260515T085324`, `rs_anarchy_run3_qtop_shifted_20260515T085324`, `rs_anarchy_run3_moreUV_20260515T085324`, `rs_anarchy_run3_moreIR_20260515T085324` | `results/figures/quark/rs_anarchy_mkk_min_hist_by_cvals.pdf` (`40f3263b61051574375aa5eca41d188d9f6d56852f9966a35f5bd6bc2a97e6da`) |
| `docs/quark_scan_methodology_note.tex:806` | RunB 95%-acceptance crossings `46.32/48.77/48.84 TeV` perturbative and `132.35/139.34/139.55 TeV` at `g_s*=3` | `rs_anarchy_runB_narrow_uniform_20260515T085324`, `rs_anarchy_runB_wide_uniform_20260515T085325`, `rs_anarchy_runB_gaussian_3sigma_20260515T085324` | `results/figures/quark/rs_anarchy_mkk_min_hist_by_yprior_gsstar.pdf` (`613cc5fd69a8e7d9fbe63450a32e964ea7e81ffd3e3f73f409a6fc6b9e880aa0`) |
| `docs/quark_scan_methodology_note.tex:836` | PDG-gate counts `1,532,640`, `110,519`, `1,164` | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_mkk_min_hist_by_pdg_tightness.pdf` (`6f53e6aebbb90230df487da6c5a3a2743322f6903fb647d19148c35b114ae8c3`) |
| `docs/quark_scan_methodology_note.tex:854` | RunC zero-pass finite-ensemble UL `p_pass <= 9.2e-7` | `rs_anarchy_runC_20260515T085323` | `results/figures/quark/rs_anarchy_cfw_comparison.pdf` (`b78b5d3c9d20fe22fc1136124da1057565a04a687b38d11068a0f81e23da6f95`) |
| `docs/quark_scan_methodology_note.tex:889` | CFW comparison: factor `2.2` stronger; matched projection `23.37 TeV`; `n=217` | `rs_anarchy_runA_20260515T085316` | `results/figures/quark/rs_anarchy_cfw_comparison.pdf` (`b78b5d3c9d20fe22fc1136124da1057565a04a687b38d11068a0f81e23da6f95`) |
| `docs/quark_scan_methodology_note.tex:892` | CFW-matched p50 projection quoted as `23.4 TeV` (rounded from `23.37 TeV`) | `rs_anarchy_runA_20260515T085316` (matched comparison computed from RUNA with plot-time CFW-era overrides) | `results/figures/quark/rs_anarchy_cfw_comparison.pdf` (`b78b5d3c9d20fe22fc1136124da1057565a04a687b38d11068a0f81e23da6f95`) |
| `docs/quark_scan_methodology_note.tex:894` | CFW no-UV-boundary-term `g_s^*=3` reference quoted as `10.5 TeV` | external reference, CFW 0804.1954, Section 2 | literature value; no local scan dir |

## Audit dependency chain

- PDG benchmark fit and strict tests (hole #4): re-derived default spurion seed and removed stale xfails (signoff: `docs/phase_logs/phase2_h4_signoff.md`)
- Bag params (hole #5): FLAG 2024 + BGS 2020 (signoff: `docs/phase_logs/phase2_h5_signoff.md`)
- Wilson RG (hole #6): BMU LO + scalar-LR sign (signoff: `docs/phase_logs/phase2_h6_signoff.md`)
- Invalidation gate clearance: `docs/phase_logs/invalidation_gate_signoff.md`
- CFW reconciliation (hole #7): factor-2.2 stronger (signoff: `docs/phase_logs/phase2_h7_signoff.md`)
- Zero-pass treatment (hole #8): Wilson upper limits (signoff: `docs/phase_logs/phase2_h8_signoff.md`)
- Scope wording (hole #9): appendix subsection (signoff: `docs/phase_logs/phase2_h9_signoff.md`)

## Environment

Conda env: `/n/home09/obarrera/.conda/envs/ising_bootstrap`  
Python: `3.11.14`  
Key package versions: `numpy==2.4.2`, `scipy==1.17.0`, `matplotlib==3.10.8`, `pandas==3.0.1`.

## Re-running this artifact

To regenerate any of the canonical scans:

```bash
sbatch <dispatcher.sbatch>  # see scripts/run_followups_README.md
```

The dispatcher pins thread-count env vars and the Python entry point; reproducibility is expected to be byte-equivalent under the same conda env plus the `quarkscan-paper-rc1` SHA. Scan outputs themselves are intentionally not committed; use the cluster paths and checksums above, mirrored in `artifacts/checksums.sha256`.
