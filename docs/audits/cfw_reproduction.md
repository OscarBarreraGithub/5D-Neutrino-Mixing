# CFW Convention-Matched Reproduction

Source extraction: `docs/audits/cfw_convention_extract.md`.  
Comparison table: `docs/audits/cfw_vs_ours.md`.  
Driver: `scripts/rs_anarchy_cfw_comparison.py`.

## CFW headline and assumptions

CFW quote a generic RS anarchic lower bound near `21 TeV` and a
pseudo-Goldstone Higgs lower bound near `33 TeV` (0804.1954, abstract).
Those abstract markers use the boundary-term choice associated with
`g_s^*~6`.  In the no-UV-boundary-term/no-bare-brane variants associated
with `g_s^*~3`, the corresponding CFW markers are `10.5 TeV` for generic
RS and about `17 TeV` for the pGB Higgs case.
For the generic RS number, their load-bearing ingredients are:

- conventional scalar-LR `C4/C5` operators with `C4=-g_R g_L/M_G^2` and
  `C5=+g_R g_L/(3M_G^2)`;
- UTfit `Delta F=2` model-independent bounds, especially
  `Lambda_Im C4K = 1.6e5 TeV`;
- direct localization/Yukawa sampling with
  `c_q3 in [0.4,0.45]`, `c_u3 in [-0.3,-0.05]`, `c_d3=-0.55`,
  `|Y_*| in [1,3]`;
- a relative `30%` maximum-distance gate on masses, CKM magnitudes, and
  `J`;
- a boundary-term color-coupling convention whose published `21 TeV`
  number is tied to `g_s^*~6`, while their no-UV-boundary-term variant
  lowers `g_s^*` to `~3` and the bound to `10.5 TeV`.

CFW do not publish explicit `B_K`, `B_4^K`, `B_5^K`, or
`epsilon_K^SM` values.  Their hadronic inputs are hidden inside the UTfit
scale bounds, so no bag constants were copied into the live code.

## Commands

Default post-audit check:

```bash
python scripts/rs_anarchy_cfw_comparison.py \
  --run scan_outputs/rs_anarchy_runA_20260515T085316 \
  --out-dir results/figures/quark \
  --summary-out /tmp/cfw_default_summary.json \
  --no-plot
```

Convention-matched CFW projection and replacement figure:

```bash
python scripts/rs_anarchy_cfw_comparison.py \
  --run scan_outputs/rs_anarchy_runA_20260515T085316 \
  --out-dir results/figures/quark \
  --eps-k-sm 1.81e-3 \
  --pdg-relative-tolerance 0.30 \
  --summary-out /tmp/cfw_matched_summary.json
```

The `--pdg-relative-tolerance 0.30` option implements the CFW relative
window `0.7--1.3` as a symmetric stored-residual factor
`1/(1-0.30)=1.428571` on masses, CKM entries, and `J`.

## Numerical result

The default post-audit curve exactly reproduces the signed-off RUNA
headline:

| Curve | Accepted rows | p50 perturbative | p50 at `g_s^*=3` | p95 at `g_s^*=3` |
|---|---:|---:|---:|---:|
| Post-audit default | 1,532,640 | 16.53995 TeV | 47.25701 TeV | 127.13273 TeV |
| CFW-matched projection | 217 | 8.18005 TeV | 23.37157 TeV | 52.80281 TeV |

Thus, under the plotted `g_s^*=3` convention and with the CFW-era epsilon
budget plus CFW 30% relative gate, our forward-only pipeline gives

```text
M_KK^min(p50) = 23.37 TeV.
```

The corrected convention-matched comparison is against CFW's
no-UV-boundary-term `10.5 TeV` RS marker, not the default `21 TeV` marker:

```text
23.37157 / 10.5 = 2.2259.
```

Equivalently, rescaling the same projection to `g_s^*~6` gives

```text
46.74314 / 21.0 = 2.2259.
```

The reconciliation is therefore not a percent-level agreement claim.  It is a
factor-`2.2` stronger RS-anarchy bound at matched `g_s^*` conventions.  The
factor-`2.2` enhancement is attributable to newer BGS-2020
`epsilon_K^SM`, FLAG-2024 bag parameters, and the audited LO BMU-corrected
sign convention relative to CFW's UTfit-era constraint inputs.  The
directions agree: both analyses imply an `M_KK > O(10 TeV)` lower bound
under anarchic flavor, but the modern pipeline does not reproduce CFW's
2008 numerical value.  The matched-gate subset has only `n=217` accepted
draws; its approximate 95% Wilson-score interval on the p50 is
`[21, 26] TeV` at `g_s^*=3`, which is much too small to erase the factor-2.2
residual.

## Step-by-step reconciliation

1. Start from the live post-audit central:

   ```text
   M0 = 47.25701 TeV
   ```

   This is BGS 2020 `epsilon_K^SM = 2.161e-3`, FLAG 2024 bags, audited LO
   Wilson running, stored factor-3 mass/CKM gate with factor-5 `J`, and
   `g_s^*=3`.

2. Switch only the central `epsilon_K` budget to the legacy CFW-era proxy:

   ```text
   budget_BGS = |2.228e-3 - 2.161e-3| = 6.70e-5
   budget_legacy = |2.228e-3 - 1.81e-3| = 4.18e-4
   budget_legacy / budget_BGS = 6.2388
   M1 = M0 / sqrt(6.2388) = 18.92 TeV
   ```

   This reproduces the expected order-of-magnitude migration from the BGS
   central budget back to the older CKMfitter-style budget.

3. Apply the CFW 30% relative gate to masses, CKM, and `J`:

   ```text
   relative window = 0.7--1.3
   symmetric stored-residual factor = 1/(1-0.30) = 1.428571
   accepted forward rows = 217
   M2 = 23.37157 TeV
   gate factor = M2 / M1 = 1.235
   ```

   The gate does not lower the bound in this fixed-c forward ensemble.  It
   selects a small, statistically rougher subset with a somewhat harder
   `M_KK^min` distribution.  This is a sampling-design effect, not a
   Wilson-coefficient convention.

4. Gauge-coupling convention:

   The regenerated figure is in the common `g_s^*=3` convention, so the
   plotted comparison uses

   ```text
   M3(g_s^*=3) = M2 = 23.37 TeV.
   ```

   If one instead multiplies the same forward projection by the literal CFW
   tree-level `g_s^*~6` value, the number becomes

   ```text
   M3(g_s^*=6) = 23.37157 * (6/3) = 46.74 TeV.
   ```

   For `g_s^*=4`, it would be

   ```text
   M3(g_s^*=4) = 23.37157 * (4/3) = 31.16 TeV.
   ```

   These are not the plotted common-convention values; they quantify the
   remaining coupling-convention ambiguity that CFW themselves discuss through
   boundary kinetic terms.

5. Compare to CFW:

   ```text
   common g_s^*=3 comparison:
     our projection 23.37 TeV / CFW no-UV-boundary RS 10.5 TeV = 2.23

   common g_s^*~6 comparison:
     our rescaled projection 46.74 TeV / CFW default RS 21 TeV = 2.23
   ```

   The paper figure now displays both CFW conventions explicitly.  The
   conclusion is a factor-2.2 stronger post-audit bound, not a percent-level
   reproduction of the CFW 2008 marker.
