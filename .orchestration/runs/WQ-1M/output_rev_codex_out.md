Integrity: recomputed 500 JSONL tiles -> 1,000,000 rows, 878,787 evaluated, 121,213 skipped, 0 malformed; all 50 `(r,M_KK)` groups have exactly 20,000 rows.
Seeds: 1,000,000 unique seed values; no duplicate seed within or across groups; every tile had 2,000 rows.
Survival samples vs report: r=0.05,M=1 strict 9958/9958=1 incl 0; r=0.25,M=1 strict 0/19994=0 incl 0; r=0.25,M=2 strict 19995/19995=1 incl 0; r=1,M=2 strict 0/19994=0 incl 0; r=1,M=3 strict 19995/19995=1 incl 0.
All 50 Survival By Group entries match the report; inclusive survival reaches 1 only at `M_KK >= 7 TeV` for every r.
Strict floor rises with r: r=0.05 -> 1 TeV, r=0.1/0.25 -> 2 TeV, r=0.5/1 -> 3 TeV; rigorous veto IDs are Delta-F=2 K001/B003/B004/B001.
Veto split is honest: rigorous and proxy are separate in rows/report/plots; proxy veto IDs are EW001 and B011/B012/B013-family, not mixed into strict survival.
No scale-mixing found: analysis code and rendered plots use exact `(quark_fit_r,M_KK)` keys; survival/veto curves are per-r/fixed-M_KK or exact-group heatmaps.
Small-r skip is honest: r=0.05 evaluated 96,838/200,000=48.419%, r=0.1 91.071%, r>=0.25 ~99.96-99.97%; skipped rows are excluded from survival denominators.
Skip reasons in raw rows: nonperturbative_quark_yukawa 120,962 and quark_fit_failed 251; no c_half_singularity rows appeared in merged JSONL.
Plots: all 11 PNGs exist, render nonblank, have correct axes/facets; report is quark-only and universal-c sanity is logged as skipped for quark-only mode, not a failure.
Discrepancy: none affecting output faithfulness.
WQ-1M-OUTPUT-REVIEW: APPROVE