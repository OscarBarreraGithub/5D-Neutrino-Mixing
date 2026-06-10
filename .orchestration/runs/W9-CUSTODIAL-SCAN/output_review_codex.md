Audited from raw JSONL, not the CSVs.

Survival counts are correct. I recomputed all 50 `(r, M_KK)` cells from raw minimal and custodial JSONL: `survival_by_r_mkk.csv` had `0` mismatches. The requested cells matched exactly:

| cell | minimal rows/eval/skip/strict/incl | custodial rows/eval/skip/strict/incl |
|---|---:|---:|
| `(1.0, 5)` | `20000/19992/8/0/0` | `20000/19992/8/19992/0` |
| `(0.25, 3)` | `20000/19989/11/0/0` | `20000/19989/11/19989/0` |
| `(1.0, 30)` | `20000/19994/6/19994/19994` | `20000/19994/6/19994/19994` |

Pairing is correct: raw keys were `1,000,000` minimal, `1,000,000` custodial, common `1,000,000`, `only_min=0`, `only_cust=0`, duplicate keys `0`.

Physics interpretation is correct. In raw rows, custodial branch is active for all custodial rows: top-level and provenance `ew_model='custodial_rs_plr'` for `1,000,000/1,000,000`. `T010` is active/evaluated but does not veto custodial anywhere; `paired_vetoes.parquet` has `T010` vetoes only for `minimal_rs` (`694123`) and none for custodial. Same-seed checks:

- `(1.0, 5)`: minimal `T010` fails `19992/19992`; custodial `T010` fails `0/19992`.
- `(0.25, 3)`: minimal `T010` fails `19989/19989`; custodial `T010` fails `0/19989`.
- `(1.0, 30)`: both pass, as expected at high scale.

The strict-vs-inclusive custodial gap is also correct. Raw custodial strict floors are `r=0.05: 1 TeV`, `r=0.1: 2 TeV`, `r=0.25: 2 TeV`, `r=0.5: 3 TeV`, `r=1.0: 3 TeV`; inclusive floor is `7 TeV` for every r. The gap is driven by proxy HARD vetoes, especially `EW001` plus CR collider proxies: e.g. `(1.0,5)` has custodial strict `19992/19992` but inclusive `0/19992`, with `EW001` and `CR001` vetoing all evaluated rows; `(0.25,3)` is vetoed by `EW001`, `CR001`, `CR012`, `CR013`.

`constraint_veto_by_r_mkk.csv` also matched the raw recomputation for the 66 checked key-veto rows in the requested cells. `paired_draws.parquet` has 1,000,000 rows, duplicate keys `0`, tile/draw ID mismatches `0`; `paired_vetoes.parquet` has 8,824,191 veto records and duplicate records `0`.

VERDICT: APPROVE