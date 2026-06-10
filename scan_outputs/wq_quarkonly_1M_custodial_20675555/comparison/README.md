# WQ Quark-Only Minimal vs Custodial Comparison

This compares baseline non-custodial `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scan_outputs/wq_quarkonly_1M_20128400` against custodial `/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/scan_outputs/wq_quarkonly_1M_custodial_20675555`.

Rows are paired by normalized `(r, mkk_tev, draw_seed)`.
Raw scan rows use `seed`, not `draw_seed`; this builder normalizes `row["seed"]`.
Raw scan rows use `params.M_KK` in GeV; this builder writes `mkk_tev` in TeV.
Raw quark-only rows use top-level `quark_fit_r`; `params.quark_fit_r` is only a fallback.

Survival uses the existing row semantics: `survives_all_HARD_strict` and `survives_all_HARD_inclusive` from each raw row.
Skipped fit rows remain in `paired_draws.parquet` but are not counted as evaluated survival in `survival_by_r_mkk.csv`.

W8 top-partner loop status for this W9 scan: deferred.
