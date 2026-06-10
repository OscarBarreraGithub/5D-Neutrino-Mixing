# Scan Explorer — DATA extraction build (codex)

Read `.orchestration/runs/WEB-EXPLORER/DESIGN.md` (Part A is authoritative). Build the Python
extraction script and PRODUCE the JSON. Do NOT touch the front-end in this task.

## Deliverables
1. `flavor_catalog/website/scripts/build_scan_explorer.py` — reads BOTH scan roots
   (`scan_outputs/wq_quarkonly_1M_20128400` minimal, `scan_outputs/wq_quarkonly_1M_custodial_20675555`
   custodial), streams the raw JSONL (do NOT load all rows; aggregate incrementally — these are 1M
   rows each), and writes `flavor_catalog/website/src/content/scan_explorer.json` per the DESIGN schema.
2. Run it and produce the JSON. Confirm size < ~1 MB.

## Correctness requirements (these get dual-reviewed against the raw scans)
- `veto[ew_model][r][cid]` = veto_fraction over EVALUATED (skipped==False) rows at each grid M_KK,
  aligned to `meta.mkk_grid_tev`. veto_fraction = (# rows where constraint failed & active & HARD) /
  (# evaluated rows) in that (ew_model,r,mkk) cell. Cross-check a few cells against
  `scan_outputs/wq_quarkonly_1M_custodial_20675555/comparison/constraint_veto_by_r_mkk.csv`.
- Only include a constraint in `constraints[]` if veto_fraction > meta.floor_threshold (0.5) for at
  least one (ew_model,r,mkk) — i.e. it actually sets a floor somewhere. Give each a human label
  (e.g. T010="Z → b b̄", EW001="Oblique S,T,U", K001="ε_K", B003/B004 ΔF=2 labels, CR00x collider
  labels). Read the real labels from `flavor_catalog/website/src/content/entries/<ID>.json` (title)
  where possible; tag (rigorous/proxy) from the row's constraint tag.
- `bare_floor_tev[ew_model][r]` = smallest grid M_KK with at least one non-skipped (fit-succeeded)
  evaluated row (the fit/masses/CKM baseline).
- `yukawa[ew_model][r][mkk]` = per the schema: p25/p50/p75 of the 3 ascending up and 3 ascending down
  fitted singular values over EVALUATED rows in that cell, plus n. (Stream a reservoir or exact
  percentiles — exact is fine via per-cell lists if memory allows per-cell, but DO NOT hold all 1M
  at once; process cell-by-cell or accumulate compact per-cell arrays.)
- `rep_matrices`: pick 3 representative SURVIVING (survives_all_HARD_strict==True) points spanning r
  (e.g. r=0.05, 0.25, 1.0), reconstruct the full 3x3 complex Y_u, Y_d from `params.quark_yukawa_seed`
  using the SAME draw/seed code path the scan harness uses (find it in
  `scripts/run_full_catalog_scan.py` / the draw module — do NOT invent a new RNG). Emit |Y_u|,|Y_d|
  (abs, rounded to ~4 sig figs). VERIFY: numpy SVD of the reconstructed Y_u/Y_d must match the row's
  `fitted_up/down_yukawa_singular_values` to < 1e-6 — assert this in the script and report it.

## After building
Run the script; print the JSON size, the constraint list, a couple cross-checked veto cells, and the
rep-matrix SVD-match residuals. Write a terse summary to
`.orchestration/runs/WEB-EXPLORER/data_impl_summary.md` ending `DATA-READY` (if JSON produced + SVD
verified) or `DATA-BLOCKED`. Do NOT commit. Use:
`source ~/.bashrc && conda activate ising_bootstrap && export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"`.
NOTE: streaming 2M rows may be slow but memory must stay bounded; if it would exceed ~a few minutes,
process per (r, shard) and aggregate.
