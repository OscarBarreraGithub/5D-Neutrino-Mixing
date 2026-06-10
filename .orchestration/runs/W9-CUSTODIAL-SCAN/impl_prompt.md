# W9 — IMPLEMENT the DUAL-APPROVED custodial-scan CODE (flag + comparison builder + sbatch)

The plan `.orchestration/runs/W9-CUSTODIAL-SCAN/plan_codex.md` is DUAL-APPROVED. Implement the CODE
portions now (NOT the SLURM run — the orchestrator launches that after this commits). If infeasible,
STOP and write `.orchestration/runs/W9-CUSTODIAL-SCAN/impl_blocker.md`.

NOTE: W7 (`8ca55c8`) and W8 (`e6df1d8`) are now committed and touched point_builder.py /
rs_ew_builder.py / EW001 / T010 / T011 / T014 / rs_ew_couplings.py / oblique_stu.py. Start from the
CURRENT tree; re-read those before editing.

## Implement (per approved plan)
1. **`--ew-model {minimal_rs,custodial_rs_plr}` flag** through `scripts/run_full_catalog_scan.py`:
   add `ew_model` to ScanConfig (default `minimal_rs`); `_config_payload` POPS it when `==minimal_rs`
   so the pinned hash `45e21a07585f7489` and WQ quark-only `c6939cc65d71f86a` stay UNCHANGED; thread
   `cfg.ew_model` into EVERY build site (quark-only `_evaluate_draw` ~465-479, full ~570-586,
   universal-c sanity ~1200-1206 and ~1261-1277) TOGETHER with `model_label=cfg.ew_model` (the
   `_validate_injected_spectrum` co-edit, rs_ew_builder.py:304-305, else it raises). Universal-c
   sanity runs under `cfg.ew_model` (the decided policy).
2. **Comparison builder** (per plan §4): pairs two run roots on identical (r, M_KK, draw_seed),
   normalizing raw `row["seed"]->draw_seed`, `params.M_KK/1000->mkk_tev`, top-level
   `quark_fit_r->r`; emits per-(r,M_KK) survival for BOTH runs + which constraint vetoes; writes the
   documented schema + README so a UI can plot minimal-vs-custodial survival-vs-M_KK per r. Use the
   exact file formats the plan specifies. Put it where the plan says (e.g. a `comparison/` builder script).
3. **Custodial sbatch / run config** reusing the SAME grid+seeds as
   `scan_outputs/wq_quarkonly_1M_20128400` (reuse `scripts/wq_quarkonly_1m_plan.py` unchanged), with
   `--ew-model custodial_rs_plr --quark-only`, output root
   `scan_outputs/wq_quarkonly_1M_custodial_<jobid>/`. Do NOT submit it.

## Tests required (per plan)
both-mode plumbing (quark-only AND full `_evaluate_draw`); worker round-trip
`_config_from_payload(_config_payload(cfg))` preserves ew_model for both modes; minimal byte-identity
(existing PR1 hash test stays green; WQ `c6939cc65d71f86a`); spectrum/ew_model mismatch raises;
comparison builder on a tiny synthetic pair of runs (seed normalization + survival pairing + enum).

## After implementing
```
source ~/.bashrc && conda activate ising_bootstrap && export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
python -m pytest -q 2>&1 | tail -30
```
Write `.orchestration/runs/W9-CUSTODIAL-SCAN/impl_summary.md` (files changed, tests added, pass count,
confirmation both pinned hashes unchanged, the custodial sbatch path) ending `IMPL-READY` (green) or
`IMPL-BLOCKED`. Do NOT commit and do NOT submit any SLURM job.
