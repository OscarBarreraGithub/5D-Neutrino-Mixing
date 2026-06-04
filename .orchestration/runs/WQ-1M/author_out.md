Plan executed: scoped serialization, sharded 1M launcher, grouped analysis/plots, tests.

r grid: `0.05,0.10,0.25,0.50,1.00`, spanning down-dominated through default `0.25` to up-dominated LH-doublet localization. Total: `10 M_KK x 5 r x 20,000 = 1,000,000` draws.

Sharding: 50 SLURM array tasks; each task runs all 10 locked `M_KK` values for one `(r, draw-shard)` with `2,000` draws/tile.

Seed proof: `base=202606040000 + 20000000*task_id`; `tile_seed=base + 1000003*tile_id`; `draw_seed=tile_seed+draw_idx`; `20000000 > 1000003*10 + 2000 = 10002030`, so shard/tile draw intervals are disjoint.

Changed/added: `scripts/run_full_catalog_scan.py`, `scripts/analyze_wq_quarkonly.py`, `scripts/wq_quarkonly_1m_plan.py`, `scripts/run_wq_quarkonly_1m_array.sbatch`, plus tests.

Serialization: quark-only rows now include `quark_fit_r` in row/params/provenance and fitted `up/down` Yukawa singular values in `fit_diagnostics`; full mode unchanged.

Smoke analysis ran on `.orchestration/runs/WQ-QUARKONLY/smoke-1k/` into `scan_outputs/wq_quarkonly_smoke_analysis_wq1m/`.

Produced 11 PNGs: survival, rigorous/proxy veto curves, rigorous/proxy rankings, Yukawa singulars, bulk-c, max-|Yq|, and per-r evolution plots; report: `analysis_report.md`.

Note: old smoke rows lacked fitted singular values, so analysis used the seed-matrix fallback only for that validation input.

Tests: `python -m pytest tests/ -q` passed: `1713 passed, 1 skipped`.

No commit made and no SLURM array submitted.

WQ-1M-AUTHOR-DONE.