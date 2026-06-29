# RS-anarchy follow-up campaign — orchestrator submission guide

> **LEGACY (May 2026 ΔF=2/anarchy follow-up campaign, pre-audit).** These runs
> predate the June 2026 audit and do not reflect the corrected minimal-RS floors
> (typical ~30 TeV ε_K; existence ~18-20 TeV S,T,U; Z→bb ~5 TeV). See
> `docs/FLOOR_SUMMARY.md` and `reports/collaborator_2026-06/CONTENT.md`.

This file lists the exact `sbatch` commands (in the right order) for the
five follow-up runs. **Do not run these blindly**: read the constraints
on each run first.

All scripts assume the repo is checked out at
`/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing` and that the
conda env `ising_bootstrap` is available at
`/n/home09/obarrera/.conda/envs/ising_bootstrap`.

The baseline run lives at:
`scan_outputs/rs_anarchy_20260507T030811/`  (800k draws, 8 tiles, factor-3
PDG match — 153,097 PDG-passing).

---

## Run 1 — PDG-tightness split (no compute, plotting only)

Already executed at code-deliverable time. Verifies:

- `results/figures/quark/rs_anarchy_mkk_min_hist_by_pdg_tightness.{pdf,png}`
- factor-3 sub-ensemble has `153,097` PDG-passing draws (matches baseline `n_pdg_pass`).

Re-run any time:

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing
/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python \
    scripts/rs_anarchy_mkk_min_hist_by_pdg_tightness.py
```

---

## Snapshot the baseline figures BEFORE Run A

Run A regenerates the headline figures with new artifacts. To preserve
the existing baseline-derived figures, copy them aside first:

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing
cp -r results/figures/quark results/figures/quark_baseline_800k
```

Run A's plot commands write into `results/figures/quark/runA/`, so the
baseline figures in `results/figures/quark/` are not overwritten anyway —
the snapshot is belt-and-braces.

---

## Run 3 — c-value sensitivity sweep (4 jobs, parallel)

Each job is `n_draws=200000`, 8 tiles. Each writes to
`scan_outputs/rs_anarchy_run3_<pattern_tag>_<TS>/`. The four submissions:

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing

sbatch --export=ALL,PATTERN_TAG=baseline,\
C_Q="0.63,0.57,0.20",\
C_U="0.66,0.50,-0.50",\
C_D="0.66,0.61,0.55" \
    scripts/run_rs_anarchy_run3.sbatch

sbatch --export=ALL,PATTERN_TAG=qtop_shifted,\
C_Q="0.63,0.57,0.30",\
C_U="0.66,0.50,-0.55",\
C_D="0.66,0.61,0.55" \
    scripts/run_rs_anarchy_run3.sbatch

sbatch --export=ALL,PATTERN_TAG=moreUV,\
C_Q="0.68,0.62,0.25",\
C_U="0.71,0.55,-0.45",\
C_D="0.71,0.66,0.60" \
    scripts/run_rs_anarchy_run3.sbatch

sbatch --export=ALL,PATTERN_TAG=moreIR,\
C_Q="0.58,0.52,0.15",\
C_U="0.61,0.45,-0.55",\
C_D="0.61,0.56,0.50" \
    scripts/run_rs_anarchy_run3.sbatch
```

After all four complete, plot:

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing
/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python \
    scripts/rs_anarchy_mkk_min_hist_by_cvals.py \
    --run baseline=scan_outputs/rs_anarchy_run3_baseline_<TS> \
    --run qtop_shifted=scan_outputs/rs_anarchy_run3_qtop_shifted_<TS> \
    --run moreUV=scan_outputs/rs_anarchy_run3_moreUV_<TS> \
    --run moreIR=scan_outputs/rs_anarchy_run3_moreIR_<TS>
```

(Replace `<TS>` with the actual timestamp — the script prints the path
on completion.)

---

## Run A — High-statistics baseline (1 job)

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing
sbatch scripts/run_rs_anarchy_runA.sbatch
```

After completion (and the snapshot above), regenerate headline figures
with the runA artifacts:

```bash
RUN_A_DIR=scan_outputs/rs_anarchy_runA_<TS>   # fill in the actual TS
/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python \
    scripts/rs_anarchy_mkk_min_hist.py \
    --draws ${RUN_A_DIR}/draws.jsonl \
    --out-dir results/figures/quark/runA \
    --label-tag runA_1Mdraws

/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python \
    scripts/plot_rs_anarchy_summary.py \
    --summary ${RUN_A_DIR}/tile_summary.json \
    --draws ${RUN_A_DIR}/draws.jsonl \
    --out-dir results/figures/quark/runA \
    --label-tag runA_1Mdraws

/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python \
    scripts/rs_anarchy_gate_sensitivity.py \
    --draws ${RUN_A_DIR}/draws.jsonl \
    --summary ${RUN_A_DIR}/tile_summary.json \
    --out-dir results/figures/quark/runA \
    --label-tag runA_1Mdraws
```

**Halt trigger** at the M_KK=3 TeV tile: expected
`n_pdg_pass = 163,600 +/- 2000`. Use:

```bash
/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python -c "
import json
d = json.load(open('${RUN_A_DIR}/tile_summary.json'))
t = next(t for t in d['tiles'] if abs(t['M_KK_GeV'] - 3000) < 1)
print(t['n_pdg_pass'], t['pdg_pass_fraction'])
"
```

If the count is outside `[161,600, 165,600]` (5-sigma window) something
has changed in the math; do not proceed until you understand why.

---

## Run B — Y-prior sensitivity (3 jobs, parallel)

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing
sbatch --export=ALL,PRIOR_TAG=narrow_uniform   scripts/run_rs_anarchy_runB.sbatch
sbatch --export=ALL,PRIOR_TAG=wide_uniform     scripts/run_rs_anarchy_runB.sbatch
sbatch --export=ALL,PRIOR_TAG=gaussian_3sigma  scripts/run_rs_anarchy_runB.sbatch
```

After all three complete, overlay-plot against the baseline:

```bash
/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python \
    scripts/rs_anarchy_mkk_min_hist_by_yprior.py \
    --run baseline=scan_outputs/rs_anarchy_20260507T030811 \
    --run narrow_uniform=scan_outputs/rs_anarchy_runB_narrow_uniform_<TS> \
    --run wide_uniform=scan_outputs/rs_anarchy_runB_wide_uniform_<TS> \
    --run gaussian_3sigma=scan_outputs/rs_anarchy_runB_gaussian_3sigma_<TS>
```

---

## Run C — Tighter PDG gate + g_s* = 4 overlay (1 job)

```bash
cd /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing
sbatch scripts/run_rs_anarchy_runC.sbatch
```

After completion:

```bash
/n/home09/obarrera/.conda/envs/ising_bootstrap/bin/python \
    scripts/rs_anarchy_cfw_comparison.py \
    --run scan_outputs/rs_anarchy_runC_<TS>
```

---

## Suggested submission order

1. (no submit) Run 1 plotter — already executed.
2. (no submit) Snapshot `cp -r results/figures/quark results/figures/quark_baseline_800k`.
3. Submit Run A (longest, 3h) and Run 3's four jobs in parallel.
4. Submit Run B's three jobs (1h each) — can overlap with the above.
5. Submit Run C (1.5h).
6. As each batch completes, run its plotter (commands above).
7. Once Run A finishes, run the regen-figures block; verify the halt
   trigger first.

All five runs are independent at the data-generation level. Plotters
that overlay multiple runs (Run 3 cvals, Run B yprior) require *all*
their input runs to be present.
