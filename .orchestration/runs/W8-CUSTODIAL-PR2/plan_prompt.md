# W8 — Custodial PR2 implementation PLAN (codex author)

You are authoring an IMPLEMENTATION PLAN only. Do NOT write production code in this run.
Output the plan to `.orchestration/runs/W8-CUSTODIAL-PR2/plan_codex.md`, end with `PLAN-READY`.

## Context — PR1 is already committed (`71b8453`)
PR1 added `ew_model="custodial_rs_plr"` (default `minimal_rs` BYTE-IDENTICAL). It already:
- zeros ONLY the protected DIAGONAL down-left `z_delta_g_L_d` (off-diagonal Z→bs/bd/sd FCNC +
  T014 UNCHANGED); residual κ_b/L (default 0); oblique keeps c_S, T→−π/(4cW²L), U=0.
- FLAGS but does NOT compute: the one-loop top-partner numerics and custodial-FCNC modeling.

Read these FIRST and ground everything in them:
- `.orchestration/runs/CUSTODIAL-RESEARCH/CORRECTED_PRESCRIPTION.md` (the approved spec)
- `quarkConstraints/rs_ew_couplings.py` (`_apply_custodial_rs_plr_proxy` ~926, build path ~490-740)
- `quarkConstraints/oblique_stu.py` + `flavor_catalog_constraints/physics_adapters/oblique_stu.py`
- `tests/test_rs_ew_custodial_pr1.py`
- Carena-Pontón-Santiago-Wagner hep-ph/0701055 Eqs. 28-30 (one-loop top-partner T and δg_bL).

## Goal (PR2 — turn the PR1 flags into computed, flagged numerics)
1. **One-loop top-partner ΔT and δg_L^b** (Carena et al Eq.28-30): add the leading custodian
   (top-partner) one-loop contributions to T and to the b_L coupling, as an EXPLICIT,
   separately-tagged term layered on the PR1 tree result. Sign is model-dependent (bidoublet
   partners push T negative) — represent the model dependence as metadata, do NOT hard-pick a
   sign silently. At M_KK~2-3 TeV expect ΔT_loop~O(0.1), δg_L^b|loop~1e-3.
2. **Custodial-FCNC modeling**: the residual off-diagonal down-left Z couplings under custodial
   (PR1 left these at minimal). Specify the representation-aware treatment (all-gen bidoublet,
   mass-basis) consistent with the prescription, with an explicit flag for what stays proxy.
3. Keep `minimal_rs` and the PR1 tree behavior byte-identical; PR2 effects only active under
   `custodial_rs_plr` and only via NEW, separately-tagged metadata fields.

## What the plan must contain
- Exact formulas (with the variables already available in `rs_ew_couplings`: a(c) tower, F(c_Q3),
  L=ln(M_Pl/TeV), c_W, v, M_KK convention — keep the physical m_1≈2.45 Λ_IR convention EXPLICIT).
- Files touched + functions; new metadata keys; how EW001/T010/T011 consume the new terms.
- A flag/knob for including-vs-deferring the loop term (must be honestly labelled in diagnostics).
- Tests: numeric oracle for ΔT_loop and δg_L^b|loop at a benchmark M_KK; minimal_rs + PR1-tree
  byte-identity regression; the "15 TeV point survives custodial" PR1 benchmark still holds;
  sign-convention metadata test; honest-omission flags present.
- Isolation: list files; flag overlap with W7 (lepton sector) and W9 (scan harness).

## Constraints
Dual-gate: a second codex + an Opus will review. Cite file:line, give formulas + sources, list
tests. Determinism, graceful degradation, no fabricated physics. End with `PLAN-READY`.
