# Phase 2 Hole #7 Implementation Log: CFW convention reconciliation

Date: 2026-05-15  
Branch: `audit/cfw-comparison`, fast-forwarded into `paper/quark-scan-2026q2`  
Scope: CFW 0804.1954 extraction, convention mapping, comparison-driver overrides, regenerated CFW comparison figure, methodology-note update, and regression test.

## Summary of CFW extraction findings

Read arXiv:0804.1954v2 from the fetched PDF/text in `/tmp/cfw`. The extracted conventions are recorded in `docs/audits/cfw_convention_extract.md`.

Key findings:

- CFW's KK-gluon operator basis is compatible with the post-audit code basis: conventional scalar-LR `C4/C5`, with `C4=-g_R g_L/M_G^2` and `C5=+g_R g_L/(3 M_G^2)`. This matches the signed-off `deltaf2.py` matching and the BMU-to-scalar map in `docs/audits/wilson_rg_inventory.md`.
- CFW's ordinary-RS scan samples localization/Yukawa choices directly: `c_q3 in [0.4,0.45]`, `c_u3 in [-0.3,-0.05]`, `c_d3=-0.55`, `|Y_*| in [1,3]`. Their pGB scan samples `c_q3 in [0.2,0.48]`, fixes `c_-d3=-0.55`, solves `c_u3` from EWSB, and requires `N_CFT >= 5`.
- CFW's Yukawa normalization maps directly to ours after the usual `246/sqrt(2)=174 GeV` convention; the difference is the prior and inverse localization scan design, not a missing normalization factor.
- CFW's published `21 TeV` RS marker is tied to a boundary-term color-coupling convention around `g_s^*~6`; their no-UV-boundary-term variant lowers `g_s^*` to `~3` and the bound to `10.5 TeV`.
- CFW use a `30%` maximum relative distance on masses, CKM magnitudes, and `J`. In our stored residuals this is implemented as a symmetric log gate factor `1/(1-0.30)=1.428571`.
- CFW do not publish explicit `B_K`, `B_4^K`, `B_5^K`, or `epsilon_K^SM` values; they import UTfit model-independent `Delta F=2` scale bounds. No live `deltaf2.py` constants were changed.

## Convention-mapping table summary

The full table is `docs/audits/cfw_vs_ours.md`.

Summary:

- Operator basis: compatible after the hole #6 BMU/scalar sign audit.
- `g_s` convention: partially compatible; comparison must state whether values are in common `g_s^*=3` units or literal CFW tree-level `g_s^*~6` units.
- Yukawa/c sampling: only partially compatible; CFW uses inverse c/localization sampling, while our RUNA sample fixes c and forward-draws `Y`.
- PDG gate: CFW's 30% relative gate maps to a symmetric factor `1.43` in our stored residual implementation.
- `epsilon_K^SM`: CFW does not quote one; the closest CFW-era proxy is the repo's legacy `1.81e-3` central value.
- Bag inputs: CFW supplies no bag table; our FLAG 2024 constants remain live and are only bypassed by plot-layer ratio rescaling where requested.
- Headline: live post-audit p50 is `47.26 TeV` at `g_s^*=3`; convention-matched p50 is `23.37 TeV` at `g_s^*=3`.

## Reproduction result and reconciliation

Commands run:

```bash
python scripts/rs_anarchy_cfw_comparison.py \
  --run scan_outputs/rs_anarchy_runA_20260515T085316 \
  --out-dir results/figures/quark \
  --summary-out /tmp/cfw_default_summary.json \
  --no-plot

python scripts/rs_anarchy_cfw_comparison.py \
  --run scan_outputs/rs_anarchy_runA_20260515T085316 \
  --out-dir results/figures/quark \
  --eps-k-sm 1.81e-3 \
  --pdg-relative-tolerance 0.30 \
  --summary-out /tmp/cfw_matched_summary.json
```

Results:

- Post-audit default: `n=1,532,640`, p50 `47.25701 TeV`, p95 `127.13273 TeV` at `g_s^*=3`.
- CFW-matched projection: `n=217`, p50 `23.37157 TeV`, p95 `52.80281 TeV` at `g_s^*=3`.

Arithmetic:

1. Start: `47.25701 TeV` (BGS 2020, FLAG 2024, LO RG, factor-3 mass/CKM gate, factor-5 `J`, `g_s^*=3`).
2. Switch to legacy `epsilon_K^SM=1.81e-3`: budget changes from `6.70e-5` to `4.18e-4`, so `47.25701 / sqrt(4.18e-4 / 6.70e-5) = 18.92 TeV`.
3. Apply CFW 30% relative mass/CKM/`J` gate: empirical forward subset gives `23.37157 TeV` (`217` rows). This gate hardens the fixed-c forward subset instead of lowering it.
4. Coupling convention: CFW's RS marker is `10.5 TeV` in the no-UV-boundary-term `g_s^*=3` variant and `21 TeV` in the default boundary-term `g_s^*~6` setup. The matched projection is `23.37 TeV` at `g_s^*=3`, or `46.74 TeV` after rescaling to `g_s^*~6`.
5. Comparison: `23.37157/10.5 = 2.23` and `46.74314/21 = 2.23`. The corrected conclusion is a factor-2.2 stronger post-audit bound at matched conventions, not a percent-level reproduction of CFW's 2008 number.

## Methodology-note update summary

Updated `docs/quark_scan_methodology_note.tex` and rebuilt `docs/quark_scan_methodology_note.pdf`.

Changes:

- Replaced the prior Run C zero-pass/qualitative-only CFW wording with the quantitative convention-matched p50 value `23.37 TeV`.
- Updated Fig. `robust-cfw` caption to describe the blue post-audit curve, red CFW-matched curve, and explicit `10.5/21/17/33 TeV` CFW convention markers.
- Cited `docs/audits/cfw_convention_extract.md` and `docs/audits/cfw_reproduction.md` for provenance.
- Corrected the earlier coupling-convention paragraph so it no longer says CFW simply uses `g_s^*=3` by default.
- Rebuilt PDF: `pdfinfo docs/quark_scan_methodology_note.pdf` reports `Pages: 18`.

## Verification

- `python scripts/rs_anarchy_cfw_comparison.py --run scan_outputs/rs_anarchy_runA_20260515T085316 --out-dir results/figures/quark --summary-out /tmp/cfw_default_summary.json --no-plot` reproduced p50 `47.26 TeV` and p95 `127.13 TeV`.
- `python scripts/rs_anarchy_cfw_comparison.py --run scan_outputs/rs_anarchy_runA_20260515T085316 --out-dir results/figures/quark --eps-k-sm 1.81e-3 --pdg-relative-tolerance 0.30 --summary-out /tmp/cfw_matched_summary.json` regenerated `results/figures/quark/rs_anarchy_cfw_comparison.{pdf,png}`.
- `pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex` from `docs/` exited 0 after fixing the `\gs` double-superscript issue.
- `/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/pytest -q tests/test_cfw_comparison.py` -> `2 passed in 16.23s`.
- `pytest -q tests/test_cfw_comparison.py tests/test_quark_deltaf2.py tests/test_wilson_rg_audit.py tests/test_qcd_running.py` -> `36 passed in 9.45s`.

Hole #7 ready for peer review.

## Peer-review revision

Reviewer finding: the first CFW reconciliation mixed `g_s^*` conventions.
CFW's published RS `21 TeV` bound belongs to the boundary-term choice
associated with `g_s^*~6`; the no-UV-boundary-term/running variant
associated with `g_s^*~3` is `10.5 TeV`. The pGB `33 TeV` marker is
likewise the small-boundary-term `g_s^*~6` setup, with the no-bare-brane
variant around `17 TeV`.

Corrected comparison numbers:

- CFW-matched projection at common `g_s^*=3`: `23.37157 TeV`.
- CFW RS marker at common `g_s^*=3`: `10.5 TeV`.
- Ratio: `23.37157 / 10.5 = 2.23`.
- Same comparison at common `g_s^*~6`: `46.74314 / 21 = 2.23`.
- pGB markers now shown separately at `17 TeV` (`g_s^*=3` no-bare-brane)
  and `33 TeV` (`g_s^*~6` default).

Physical interpretation: the reconciliation is honest but not an
agreement-within-errors reproduction. The post-audit pipeline gives a
factor-2.2 stronger RS-anarchy bound at matched conventions. The residual
is attributable to newer BGS-2020 `epsilon_K^SM` inputs, FLAG-2024 bag
parameters, and the audited LO BMU-corrected sign convention relative to
CFW's UTfit-era constraint inputs. Both analyses still point to
`M_KK > O(10 TeV)` under anarchic flavor.

Revision artifacts:

- `scripts/rs_anarchy_cfw_comparison.py` now draws four CFW vertical
  markers: RS `10.5/21 TeV` and pGB `17/33 TeV`.
- `results/figures/quark/rs_anarchy_cfw_comparison.{pdf,png}` were
  regenerated with `--eps-k-sm 1.81e-3 --pdg-relative-tolerance 0.30`.
- `docs/quark_scan_methodology_note.tex` now leads the CFW subsection with
  the factor-2.2 finding and includes the `n=217` finite-statistics caveat
  (`[21, 26] TeV` Wilson-score p50 interval at `g_s^*=3`).
- `docs/audits/cfw_reproduction.md` and `docs/audits/cfw_vs_ours.md` now
  state the matched-convention factor-2.2 result.
- `tests/test_cfw_comparison.py` now asserts the factor-2.2 matched
  convention result.

Verification:

- `pdflatex -interaction=nonstopmode quark_scan_methodology_note.tex` from
  `docs/` exited 0.
- `/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/pytest -q tests/test_cfw_comparison.py`
  exited 0 with `2 passed`.

Hole #7 ready for re-review.
