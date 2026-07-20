# Phase 2 Hole #7 Review: CFW reconciliation

### Verdict
REJECT-WITH-REVISIONS

The implementation mechanics are mostly solid: the branch is pushed, the source extraction is sound, the matched-driver rerun reproduces `n=217` and `23.37 TeV`, the plot exists, and the requested tests pass. The blocker is the physics comparison: the paper still calls `23.37 TeV` vs CFW `21 TeV` an 11% matched-convention agreement in `g_s^*=3` units, but CFW Section 2 says the `21 TeV` number is tied to the boundary-term convention with `g_s^* ~ 6`; their `g_s^* ~ 3` no-UV-boundary-term variant is `10.5 TeV`. The apples-to-apples comparison is therefore either `23.37 TeV` vs `10.5 TeV` at `g_s^*=3`, or `46.74 TeV` vs `21 TeV` at `g_s^*~6`. Both are a factor `2.23`, not 11%.

### Item-by-item (1-8)

1. PASS - Branch and push. Current branch is `paper/quark-scan-2026q2`; `HEAD` and `origin/paper/quark-scan-2026q2` both resolve to `330ffa9fbcd85048957a98391453e3a24f1f232d`. `audit/cfw-comparison` and `origin/audit/cfw-comparison` both resolve to `b47aa769990febaa7761e3cad56e618394496f02`, and that tip is an ancestor of the paper branch. The five requested commits appear linearly before the implementation-log commit.

2. PASS - CFW PDF extraction. `/tmp/cfw/0804.1954.pdf` and `/tmp/cfw/0804.1954.txt` exist. I independently pulled `https://arxiv.org/pdf/0804.1954` to `/tmp/cfw_webcheck.pdf` and converted it with `pdftotext`; the arXiv abstract also confirms the `21 TeV` and `33 TeV` headline values. The independent PDF text confirms: Section 2/Fig. 1 gives `c_q3 in [0.4,0.45]`, `c_u3 in [-0.3,-0.05]`, `c_d3=-0.55`, `|Y_*| in [1,3]`; Section 5.3 gives `c_q3 in [0.2,0.48]`, `c_-d3=-0.55`, and the `30%` maximum relative-distance gate. Section 2 around eqs. 2.20-2.21 states `g_s^* ~ 6`, then says no UV boundary kinetic term changes `g_s^*` from about `6` to about `3` and reduces the bound from `21 TeV` to `10.5 TeV`. Section 5.3 similarly ties the `33 TeV` pGB estimate to `g_s^* ~ 6` and quotes a relaxed no-bare-brane case around `17 TeV`.

3. FAIL - Convention-mapping table. The operator-basis row is acceptable: CFW's eq. 2.18 has `C4=-g_R g_L/M_G^2`, `C5=+g_R g_L/(3 M_G^2)`, and the code uses `C4_LR=-LR/M_KK^2`, `C5_LR=LR/(3 M_KK^2)`. The code's `_kaon_matrix_elements()` documents positive conventional `<O5_LR>` with positive `B_5_K`, and CFW does not publish a separate bag/matrix-element table; their use of UTfit conventional `C4/C5` bounds supports the no-sign-flip conclusion. The PDG gate conversion `30% -> 1/(1-0.30)=1.428571` is correct. The `epsilon_K^SM` and bag-input rows are correct: CFW use UTfit `Delta F=2` scale bounds and do not quote explicit `epsilon_K^SM`, `B_K`, `B_4^K`, or `B_5^K`. The failure is the table's headline row: it says the `23.37 TeV` result is within about `11%` of CFW's `21 TeV` after convention matching, while the same table correctly states that CFW's `g_s^*=3` equivalent is `10.5 TeV`.

4. FAIL - Reconciliation arithmetic. Steps 1-3 check out. `scan_outputs/followup_crossings_summary.json` gives RUNA p50 `47.25701174470522 TeV` at `g_s^*=3`. The budget switch factor is `sqrt(4.18e-4 / 6.7e-5)=2.497760`, so `47.2570117/2.497760=18.919755 TeV`. The exact requested rerun,

```bash
python scripts/rs_anarchy_cfw_comparison.py --run scan_outputs/rs_anarchy_runA_20260515T085316 --eps-k-sm 1.81e-3 --pdg-relative-tolerance 0.30 --no-plot
```

returned `cfw_matched: n=217`, p50 `23.37 TeV`, p95 `52.80 TeV` at `g_s^*=3`. The gate hardening is therefore empirically reproduced. The coupling step is wrong: if the figure is in common `g_s^*=3` units, the CFW marker to compare against is `10.5 TeV`, not `21 TeV`. Equivalently, comparing to CFW's literal `21 TeV` convention requires scaling our matched curve to `46.74 TeV`. The report does acknowledge the small subset as `statistically rougher`, but the paper text should make the resulting percentile uncertainty explicit; `n=217` gives a broad median-rank interval, not a precision 11% test.

5. PASS - Comparison plot integrity, with physics caveat. `pdftotext` and image inspection of `results/figures/quark/rs_anarchy_cfw_comparison.pdf/.png` show the blue post-audit default curve at `g_s^*=3`, the red dashed CFW-matched projection at `g_s^*=3`, and legend entries for `CFW generic anarchic (21 TeV)` and `CFW PGB Higgs (33 TeV)`. The 21 TeV marker is grey dotted and the 33 TeV marker is black dashed, so the styling is not exactly two dashed markers, but the labels are present. The physics caveat remains: those vertical markers are CFW abstract/default-convention values on an axis labeled `g_s^*=3 convention`.

6. PASS - Tests. The requested command passed:

```text
pytest -q tests/test_cfw_comparison.py tests/test_quark_deltaf2.py tests/test_wilson_rg_audit.py tests/test_qcd_running.py
....................................                                     [100%]
36 passed in 11.35s
```

7. FAIL - Methodology-note text update. The old qualitative-band placeholder is gone; `23.37 TeV`, the 11% claim, and `docs/audits/cfw_*` citations are present. The `g_s^*=3` vs `g_s^*~6` ambiguity is documented in adjacent prose. However, the text still says `23.37 TeV` is `11%` above CFW's `21 TeV` marker in the common `g_s^*=3` plotting convention, then admits the `21 TeV` marker is tied to `g_s^*~6`. That is not a resolved ambiguity; it is a contradiction. The note should either rescale CFW's marker to `10.5 TeV` for the `g_s^*=3` plot, or rescale our curve to `46.74 TeV` before comparing to CFW's `21 TeV` value, and then narrow/remove the agreement claim.

8. PASS - No drift into live code/constants/scan outputs. `git show --stat ffec986 06d5d85 da8f647 508ae69 b47aa76` touches only `docs/audits/cfw_*`, `scripts/rs_anarchy_cfw_comparison.py`, `results/figures/quark/rs_anarchy_cfw_comparison.{pdf,png}`, `docs/quark_scan_methodology_note.{tex,pdf}`, and `tests/test_cfw_comparison.py`. It does not touch `deltaf2.py`, bag constants, Wilson-RG implementation, or scan-output data.

### Physics-decision flags

- `g_s^*` ambiguity in CFW: CFW's published `21 TeV` RS bound is under the boundary-term choice associated with `g_s^* ~ 6`; the no-UV-boundary-term/running variant associated with `g_s^* ~ 3` is `10.5 TeV`. The pGB `33 TeV` marker is also quoted in the small-boundary-term `g_s^* ~ 6` setup, with a relaxed value around `17 TeV` in the no-bare-brane case. The paper must state which CFW convention is being placed on the plot axis.

- `23.37` vs `21 TeV`: the 11% agreement claim is not honest as a matched-convention comparison. It compares our `g_s^*=3` curve to CFW's `g_s^*~6` marker. At matched `g_s^*=3`, CFW's RS marker is `10.5 TeV`, so `23.37/10.5=2.23`. At matched CFW default `g_s^*~6`, our curve is `46.74 TeV`, so `46.74/21=2.23`. This is a real factor-of-two residual, not a percent-level reproduction.

- Small matched subset: the implementation report and `docs/audits/cfw_reproduction.md` acknowledge that the `n=217` subset is statistically rougher, and the rerun confirms that count. The methodology note should not present the median as a high-precision 11% validation without an explicit small-sample percentile caveat.

### Final
Send back because the `g_s^*` convention-matched comparison uses the wrong CFW marker. The source extraction, script machinery, plot generation, and tests are otherwise ready once the paper/report replace the 11% agreement claim with the correct convention-rescaled comparison.

===PHASE_2_H7_REVIEW_END===
