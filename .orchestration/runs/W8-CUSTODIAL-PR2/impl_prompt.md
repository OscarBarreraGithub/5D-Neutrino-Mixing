# W8 — IMPLEMENT the DUAL-APPROVED custodial PR2 plan

The plan `.orchestration/runs/W8-CUSTODIAL-PR2/plan_codex.md` is DUAL-APPROVED (codex + Opus, both
re-reviewed; oracles independently reproduced). Implement it EXACTLY — production code + tests. If
the plan proves infeasible, STOP and write the blocker to
`.orchestration/runs/W8-CUSTODIAL-PR2/impl_blocker.md` rather than improvising.

NOTE: W7 (μ→eγ LMFV) just committed (`8ca55c8`) and touched `point_builder.py` + `rs_ew_builder.py`.
Start from the CURRENT tree (git pull of state is already local). Re-read those two files before editing.

## Non-negotiables (from the approved plan)
- `minimal_rs` BYTE-IDENTICAL; PR1 tree behavior unchanged; PR1 pinned hash `45e21a07585f7489`
  preserved; the "15 TeV point survives custodial" benchmark still holds.
- NO fabricated negative ΔT: `top_partner_loop_t_sign=-1` ALONE must raise; a negative ΔT requires a
  finite numeric `top_partner_loop_delta_t_override`. Vertex proxy is never relabeled as a bidoublet T calc.
- Thread `top_partner_loop_delta_t_override` through point_builder, rs_ew_builder, AND the
  `build_rs_ew_couplings` core signature. Replace the second raise at
  `quarkConstraints/rs_ew_couplings.py:491-492` with `ew_model=="minimal_rs"`-gated dispatch.
- Carena δg_L^b added to `z_delta_g_L_d[2,2]` in the repo's additive `(t3 − Q s²)` convention with NO
  spurious s_z/g_Z/g_L^SM factor (the approved normalization map). Do NOT touch `z_delta_g_R_d`.
- T014: keep RH minimal; assert LH off-diagonals zeroed and T014 == RH-minimal contribution; never
  claim total T014 rate is zero.
- Deterministic missing-sign: "compute magnitudes but do not apply" (never raise solely for missing sign).
- Three distinct flags: `top_partner_zbb_loop_numerics_included`, `top_partner_t_loop_numerics_included`,
  and the EW001-facing alias — keep them separate, not one overloaded boolean.
- EW001 passes `delta_t_loop` only when the T-loop flag is True (else 0.0); keep its ValueError path.

## Tests required (per plan)
singlet ΔT/δg oracle; bidoublet-vertex δg (genuinely negative) oracle; T010 pseudo-observable oracle;
T014 RH-minimal oracle; minimal_rs + PR1 byte-identity (hash 45e21a07585f7489); 15 TeV survives;
loop-deferred custodial == PR1 tree-only T; missing-sign computes-but-not-applies; override threading;
honest-omission flags present.

## After implementing
```
source ~/.bashrc && conda activate ising_bootstrap && export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
python -m pytest -q 2>&1 | tail -30
```
Write `.orchestration/runs/W8-CUSTODIAL-PR2/impl_summary.md` (files changed, tests added, pass count,
the oracle numbers, confirmation minimal_rs hash unchanged + 15 TeV survives) ending `IMPL-READY`
(green) or `IMPL-BLOCKED`. Do NOT commit — the orchestrator commits after dual review.
