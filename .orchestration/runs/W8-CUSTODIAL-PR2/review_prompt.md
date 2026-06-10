# W8 plan REVIEW — independent critique

Independently critique the plan at `.orchestration/runs/W8-CUSTODIAL-PR2/plan_codex.md`. Read the
requirements `.orchestration/runs/W8-CUSTODIAL-PR2/plan_prompt.md`, the approved spec
`.orchestration/runs/CUSTODIAL-RESEARCH/CORRECTED_PRESCRIPTION.md`, and verify against the REAL
repo (`quarkConstraints/rs_ew_couplings.py`, `oblique_stu.py`, `tests/test_rs_ew_custodial_pr1.py`).
Do NOT rubber-stamp.

Check specifically:
1. **Loop formulas:** are the Carena et al hep-ph/0701055 Eq.28-30 one-loop top-partner ΔT and
   δg_L^b correctly transcribed, with the model-dependent SIGN represented as metadata (NOT
   silently hard-picked)? Magnitudes at M_KK~2-3 TeV sane (ΔT_loop~O(0.1), δg_L^b|loop~1e-3)?
2. **PR1 invariants preserved:** minimal_rs byte-identical; PR1 tree behavior unchanged; PR2
   effects only under custodial_rs_plr via NEW separately-tagged metadata; the "15 TeV point
   survives custodial" PR1 benchmark still holds.
3. **Custodial-FCNC modeling:** representation-aware (all-gen bidoublet, mass-basis), consistent
   with the prescription, with explicit proxy flags for what stays deferred.
4. **M_KK convention** kept explicit (physical m_1≈2.45 Λ_IR); coefficients applied in matching
   convention.
5. **Tests + honesty:** loop oracle, byte-identity regression, sign-metadata test, honest-omission
   flags; determinism; no fabricated physics.

End with EXACTLY one line: `VERDICT: APPROVE` or `VERDICT: NEEDS-FIXES`, then numbered fixes.
