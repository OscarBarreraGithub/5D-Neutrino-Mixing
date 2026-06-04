# WQ-1M ANALYSIS-OUTPUT INDEPENDENT REVIEW (codex reviewer, gpt-5.x xhigh)

Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. The 1M quark-only r×M_KK scan COMPLETED; data merged to `scan_outputs/wq_quarkonly_1M_full/` (5 r-labels × 10 shards = 500 tile summaries, 1,000,000 rows), analysis written to `scan_outputs/wq_quarkonly_1M_full/analysis/` (`analysis_report.md` + 11 PNGs). INDEPENDENTLY verify the OUTPUTS are faithful to the raw data — do NOT trust the report; recompute from the JSONL.

## VERIFY (recompute from `scan_outputs/wq_quarkonly_1M_full/**/tile-*.jsonl`)
1. **Completeness/integrity:** exactly 1,000,000 rows; all 5 r values {0.05,0.1,0.25,0.5,1} × 10 M_KK {1,2,3,5,7,10,15,20,30,50} TeV present with 20,000 rows each; NO duplicate (seed) within or across groups (spot-check seeds are unique). Confirm the merge lost nothing.
2. **Survival numbers:** independently recompute strict & inclusive HARD survival fractions for a sample of (r, M_KK) groups directly from the rows (`survives_all_HARD_strict/_inclusive`) and confirm they MATCH the report's "Survival By Group" table. In particular verify the headline transitions: r=0.05 strict survives even at 1 TeV; r=0.25 strict dies at 1 TeV, survives ≥2 TeV; r=1 strict dies at ≤2 TeV, survives ≥3 TeV; inclusive survival only reaches 1 at M_KK ≥ 7 TeV (EW001 oblique reach).
3. **No scale-mixing:** confirm the per-r evolution plots and survival/veto curves group by EXACT (r, M_KK) and never pool across r or across M_KK. Confirm rigorous vs proxy vetoes are reported SEPARATELY and honestly (rigorous = ΔF=2 K001/B003/B004…; proxy = EW001 oblique, B011/B012/B013 radiative).
4. **r-evolution physics:** confirm the data really shows the rigorous ΔF=2 bound STRENGTHENING as r grows (more up-Yukawa weight in C_Q = r·YuYu†+YdYd† → tighter eps_K/Delta-m → higher M_KK floor). Note the evaluated-count drop at small r (r=0.05 ~48-50% evaluated vs ~100% at r≥0.25) and CONFIRM it is an honest skip (quark-fit / c=0.5 singularity), not silently miscounted survival.
5. **Plots render & are labeled correctly** (axes, per-r facets, scale labels); the report is labeled as a quark-only run and does NOT present the intentionally-skipped universal-c sanity as a failure.

## OUTPUT (≤14 lines)
Recomputed integrity (row count, group coverage, dup check); your independently-recomputed survival for ≥4 sample groups vs the report; confirm no scale-mixing + honest rigorous/proxy split; confirm the r-evolution trend + the small-r skip is honest; any discrepancy. END with EXACTLY ONE line:
`WQ-1M-OUTPUT-REVIEW: APPROVE`  or  `WQ-1M-OUTPUT-REVIEW: NEEDS-FIXES`
