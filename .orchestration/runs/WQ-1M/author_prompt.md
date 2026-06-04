# WQ-1M — SHARDED 1M QUARK-ONLY SCAN (r×M_KK sweep) + ANALYSIS/PLOTS (codex author, gpt-5.x xhigh)

Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing. You are the AUTHOR. First a SHORT plan, then implement. A separate codex reviewer AND an Opus reviewer will INDEPENDENTLY review; BOTH must APPROVE before anything is committed or any big job is submitted. Build on the DUAL-APPROVED quark-only mode (commit 070cdb6) in `scripts/run_full_catalog_scan.py` and the single-node sbatch `scripts/run_wq_quarkonly.sbatch` (commit d5d3e72). DO NOT change constraint physics or the quark/builder math.

## PHYSICS CONTEXT — r
`r` is the MFV up/down weight in the LH-doublet bulk-mass matrix: `C_Q = r·(Y_u Y_u†) + (Y_d Y_d†)` (quarkConstraints/model.py:237), dimensionless, fixed at 0.25 in the current scan. It is NOT an RG scale. The physicist wants to see the quark Yukawa / localization structure at INDEPENDENT r values, plotted SEPARATELY (per-r and per-M_KK), NOT aggregated across scales. So the 1M run must sweep r over a small independent grid, and the analysis must never mix different r (or different M_KK) into one cloud.

## GOAL (≈1,000,000 draws total, quark-only, sharded on SLURM, fast)
1. **r grid:** propose a small set of independent r values (e.g. spanning down-dominated → comparable → up-dominated doublet localization, around the 0.25 default) and JUSTIFY physically. Keep total ≈1M with the locked M_KK grid `1,2,3,5,7,10,15,20,30,50` TeV. E.g. 10 M_KK × 5 r × 20k = 1M (you choose N_r/draws; document).
2. **Record r per output row AND in provenance** (currently absent) so the analysis can separate by r. Minimal serialization change to `scripts/run_full_catalog_scan.py`; default behavior otherwise byte-identical; full-mode unchanged. ALSO serialize the fitted quark Yukawa singular values (up & down) per row if not already present — needed for Yukawa-structure plots (the row already has bulk_c_Q/c_u/c_d, max_abs_quark_yukawa, and the drawn seed; add fitted singular values if cheap and absent).
3. **Sharded SLURM array launcher** (new sbatch + any tiny helper): shard the (r, draw-range) space across array tasks, each running the full M_KK grid for one r-shard with a DISJOINT `--base-seed` so NO two draws ever share a seed. **This is THE correctness risk — PROVE seed-disjointness** (tile seed = base + stride·tile_id, draw seed = tile.seed + draw_idx; choose per-shard base offsets > stride·n_tiles + n_draws). Resume-safe (per-tile summaries). Output under `scan_outputs/wq_quarkonly_1M_<ts>/r<value>/...`. Partition serial_requeue, account randall_lab, generous cores (e.g. 48) + mem; many concurrent tasks OK (the user has authorized aggressive compute). Include a merge step (or document the merge the analysis does).
4. **Analysis + plotting script** `scripts/analyze_wq_quarkonly.py` (new): read the merged/sharded JSONL and produce PNGs + a short markdown report, WITHOUT mixing scales:
   - **M_KK plots:** overall survival fraction vs M_KK; per-constraint survival/veto fraction vs M_KK; a per-constraint constraining-power ranking (which process kills the most points, rigorous vs proxy separated honestly).
   - **Yukawa plots:** distributions of fitted up/down Yukawa singular values and of bulk_c_Q/c_u/c_d (the geometric localization); max_abs_quark_yukawa distribution (perturbativity).
   - **Per-r EVOLUTION (the key ask):** separate panels/curves per r value showing how the fitted Yukawas + c_Q localization + survival evolve as r changes — explicitly NOT aggregated across r, and NOT aggregated across M_KK (facet by both). Label scales clearly.
   - Save to `scan_outputs/.../plots/` + `analysis_report.md`. Use matplotlib (Agg backend, no display).
   - Add a unit test that runs the analysis on the EXISTING `.orchestration/runs/WQ-QUARKONLY/smoke-1k/` output (schema-identical) and asserts it emits the expected PNGs + a non-empty report without error.

## GATES / TESTS
- `python -m pytest tests/ -q` stays green + the new analysis test.
- Validate the analysis on the smoke-1k output NOW (do not wait for the 1M; the orchestrator submits the 1M array after the gate). Report what plots it produced.
- Do NOT submit the 1M SLURM array yourself — the orchestrator submits it after BOTH reviewers APPROVE.

## OUTPUT (≤18 lines)
Short plan; the proposed r grid + total decomposition (M_KK×r×draws=≈1M) + justification; the seed-disjointness proof (formula + the per-shard base-offset you chose); files changed/added; serialization additions (r + fitted Yukawa singular values); the analysis run on smoke-1k (which PNGs + report); pytest counts. End with: WQ-1M-AUTHOR-DONE.
