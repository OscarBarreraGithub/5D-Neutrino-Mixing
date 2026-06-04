# WQ Quark-Only Analysis Report

- input JSONL files: 500
- rows: 1000000
- evaluated rows: 878787
- skipped rows: 121213
- malformed rows skipped: 0
- r values: 0.05, 0.1, 0.25, 0.5, 1
- M_KK values [TeV]: 1, 2, 3, 5, 7, 10, 15, 20, 30, 50
- grouping rule: every fraction and distribution is grouped by exact (quark_fit_r, M_KK); no rows are pooled across r or M_KK.
- singular-value source counts: {'fitted': 878787}

## Plots

- plots/survival_vs_mkk.png
- plots/constraint_veto_fraction_rigorous.png
- plots/constraint_veto_fraction_proxy.png
- plots/constraining_power_ranking_rigorous.png
- plots/constraining_power_ranking_proxy.png
- plots/yukawa_singular_values.png
- plots/bulk_c_distributions.png
- plots/max_abs_quark_yukawa.png
- plots/per_r_survival_by_mkk.png
- plots/per_r_cq_localization_by_mkk.png
- plots/per_r_yukawa_singular_values_by_mkk.png

## Survival By Group

| r | M_KK [TeV] | rows | evaluated | strict survival | inclusive survival |
|---:|---:|---:|---:|---:|---:|
| 0.05 | 1 | 20000 | 9958 | 1 | 0 |
| 0.05 | 2 | 20000 | 9947 | 1 | 0 |
| 0.05 | 3 | 20000 | 9846 | 1 | 0 |
| 0.05 | 5 | 20000 | 9753 | 1 | 0 |
| 0.05 | 7 | 20000 | 9728 | 1 | 1 |
| 0.05 | 10 | 20000 | 9661 | 1 | 1 |
| 0.05 | 15 | 20000 | 9494 | 1 | 1 |
| 0.05 | 20 | 20000 | 9518 | 1 | 1 |
| 0.05 | 30 | 20000 | 9551 | 1 | 1 |
| 0.05 | 50 | 20000 | 9382 | 1 | 1 |
| 0.1 | 1 | 20000 | 18265 | 0 | 0 |
| 0.1 | 2 | 20000 | 18328 | 1 | 0 |
| 0.1 | 3 | 20000 | 18279 | 1 | 0 |
| 0.1 | 5 | 20000 | 18231 | 1 | 0 |
| 0.1 | 7 | 20000 | 18215 | 1 | 1 |
| 0.1 | 10 | 20000 | 18212 | 1 | 1 |
| 0.1 | 15 | 20000 | 18211 | 1 | 1 |
| 0.1 | 20 | 20000 | 18206 | 1 | 1 |
| 0.1 | 30 | 20000 | 18145 | 1 | 1 |
| 0.1 | 50 | 20000 | 18049 | 1 | 1 |
| 0.25 | 1 | 20000 | 19994 | 0 | 0 |
| 0.25 | 2 | 20000 | 19995 | 1 | 0 |
| 0.25 | 3 | 20000 | 19990 | 1 | 0 |
| 0.25 | 5 | 20000 | 19999 | 1 | 0 |
| 0.25 | 7 | 20000 | 19996 | 1 | 1 |
| 0.25 | 10 | 20000 | 19997 | 1 | 1 |
| 0.25 | 15 | 20000 | 19995 | 1 | 1 |
| 0.25 | 20 | 20000 | 19998 | 1 | 1 |
| 0.25 | 30 | 20000 | 19989 | 1 | 1 |
| 0.25 | 50 | 20000 | 19995 | 1 | 1 |
| 0.5 | 1 | 20000 | 19993 | 0 | 0 |
| 0.5 | 2 | 20000 | 19995 | 0 | 0 |
| 0.5 | 3 | 20000 | 19999 | 1 | 0 |
| 0.5 | 5 | 20000 | 19998 | 1 | 0 |
| 0.5 | 7 | 20000 | 19992 | 1 | 1 |
| 0.5 | 10 | 20000 | 19992 | 1 | 1 |
| 0.5 | 15 | 20000 | 19992 | 1 | 1 |
| 0.5 | 20 | 20000 | 19991 | 1 | 1 |
| 0.5 | 30 | 20000 | 19994 | 1 | 1 |
| 0.5 | 50 | 20000 | 19994 | 1 | 1 |
| 1 | 1 | 20000 | 19990 | 0 | 0 |
| 1 | 2 | 20000 | 19994 | 0 | 0 |
| 1 | 3 | 20000 | 19995 | 1 | 0 |
| 1 | 5 | 20000 | 19991 | 1 | 0 |
| 1 | 7 | 20000 | 19993 | 1 | 1 |
| 1 | 10 | 20000 | 19993 | 1 | 1 |
| 1 | 15 | 20000 | 19990 | 1 | 1 |
| 1 | 20 | 20000 | 19992 | 1 | 1 |
| 1 | 30 | 20000 | 19992 | 1 | 1 |
| 1 | 50 | 20000 | 19990 | 1 | 1 |

## Top HARD Vetoes By Group

| r | M_KK [TeV] | rigorous top vetoes | proxy top vetoes |
|---:|---:|---|---|
| 0.05 | 1 | none | B012 1.000, EW001 1.000 |
| 0.05 | 2 | none | EW001 1.000 |
| 0.05 | 3 | none | EW001 1.000 |
| 0.05 | 5 | none | EW001 1.000 |
| 0.05 | 7 | none | none |
| 0.05 | 10 | none | none |
| 0.05 | 15 | none | none |
| 0.05 | 20 | none | none |
| 0.05 | 30 | none | none |
| 0.05 | 50 | none | none |
| 0.1 | 1 | K001 1.000, B003 0.639 | B012 1.000, EW001 1.000, B011 0.331 |
| 0.1 | 2 | none | B012 1.000, EW001 1.000 |
| 0.1 | 3 | none | EW001 1.000 |
| 0.1 | 5 | none | EW001 1.000 |
| 0.1 | 7 | none | none |
| 0.1 | 10 | none | none |
| 0.1 | 15 | none | none |
| 0.1 | 20 | none | none |
| 0.1 | 30 | none | none |
| 0.1 | 50 | none | none |
| 0.25 | 1 | B003 1.000, B004 1.000, K001 1.000 | B011 1.000, B012 1.000, EW001 1.000 |
| 0.25 | 2 | none | B012 1.000, EW001 1.000 |
| 0.25 | 3 | none | EW001 1.000 |
| 0.25 | 5 | none | EW001 1.000 |
| 0.25 | 7 | none | none |
| 0.25 | 10 | none | none |
| 0.25 | 15 | none | none |
| 0.25 | 20 | none | none |
| 0.25 | 30 | none | none |
| 0.25 | 50 | none | none |
| 0.5 | 1 | B003 1.000, B004 1.000, K001 1.000 | B011 1.000, B012 1.000, B013 1.000 |
| 0.5 | 2 | K001 1.000 | B012 1.000, EW001 1.000 |
| 0.5 | 3 | none | EW001 1.000 |
| 0.5 | 5 | none | EW001 1.000 |
| 0.5 | 7 | none | none |
| 0.5 | 10 | none | none |
| 0.5 | 15 | none | none |
| 0.5 | 20 | none | none |
| 0.5 | 30 | none | none |
| 0.5 | 50 | none | none |
| 1 | 1 | B003 1.000, B004 1.000, K001 1.000 | B011 1.000, B012 1.000, B013 1.000 |
| 1 | 2 | K001 1.000, B004 0.573 | B012 1.000, EW001 1.000 |
| 1 | 3 | none | B012 1.000, EW001 1.000 |
| 1 | 5 | none | EW001 1.000 |
| 1 | 7 | none | none |
| 1 | 10 | none | none |
| 1 | 15 | none | none |
| 1 | 20 | none | none |
| 1 | 30 | none | none |
| 1 | 50 | none | none |
