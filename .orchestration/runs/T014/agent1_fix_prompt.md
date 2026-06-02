# T014 FIX PASS (agent1) — you implemented T014; BOTH reviewers found blockers. Fix them, keep physics otherwise intact, re-run tests. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.

Constraint: T014 — flavor-changing Z decays (down-sector FCNC: Z→bs̄+b̄s, Z→bd̄+b̄d, Z→sd̄+s̄d). SECONDARY/top_higgs_ew/HARD. Your files:
- flavor_catalog_constraints/secondary/top_higgs_ew/T014.py
- flavor_catalog_constraints/physics_adapters/zpole.py (you appended FCNC helpers here)
- tests/constraints/secondary/top_higgs_ew/test_T014.py

## BLOCKER 1 (agent2, PHYSICS) — wrong branching-fraction denominator
The catalog limit B(Z→qq̄') < 2.9e-3 is a fraction of the TOTAL Z width, but your code computes Γ_FCNC / (Γ_had^SM + Γ_FCNC) using only the hadronic width (sum over u,d,s,c,b) in physics_adapters/zpole.py (~:273, :289-290). This makes the HARD veto ~1/0.699 ≈ 1.43× too tight (your threshold corresponds to a true B≈2.0e-3, not 2.9e-3).
**Fix:** use the TOTAL SM Z width in the denominator: B = Γ_FCNC / (Γ_Z^total,SM + Γ_FCNC), where Γ_Z^total,SM ≈ 2.4952 GeV (or whatever the zpole core already exposes as the total Z width — REUSE the existing total-width value from quarkConstraints/zpole.py, do not hardcode a fresh literal if the core has it). Confirm the ratio Γ_had/Γ_Z ≈ 0.699 explains the shift.

## BLOCKER 2 (agent3, CODE/architecture) — physics must live in the core module, adapter must be a thin wrapper
You implemented the FCNC width/proxy physics INSIDE physics_adapters/zpole.py; quarkConstraints/zpole.py has no off-diagonal down-sector helper. This violates the repo boundary (constraints → adapter → core physics; the adapter must be a thin wrapper, real physics lives in quarkConstraints/).
**Fix:** add an APPEND-ONLY off-diagonal down-sector FCNC-Z helper (the width/coupling functions + any dataclass) to quarkConstraints/zpole.py, and reduce physics_adapters/zpole.py to a thin wrapper that calls it. Do NOT modify any pre-existing quarkConstraints/zpole.py function (append only). Keep the existing T010/T011 zpole behavior unchanged.

## SHOULD-FIX (agent2) — effective-coupling diagnostic normalization
The effective-coupling diagnostic (T014.py ~:308 / zpole.py ~:330-352) inherits the same denominator problem (gives sqrt(|δgL|²+|δgR|²)=0.03497 for 2.9e-3; a total-Z denominator in the charge-summed convention gives ≈0.0421). Recompute it consistently with the total-Z-width normalization once Blocker 1 is fixed.

## Rules
- Keep numeric ConstraintResult fields real floats; complex couplings only in diagnostics. Severity HARD. Missing-coupling path stays non-crashing.
- Anchors via load_anchor (loud fail). Keep NEEDS-HUMAN-PHYSICS flag for the RS FCNC-Z matching proxy.
- Tests must cross-check INDEPENDENTLY of the adapter by recomputing from quarkConstraints/zpole.py core. Update any pinned numbers that legitimately change due to the denominator fix; a safe point passes, an excluded point fails.
- Re-run: `python -m pytest tests/constraints/secondary/top_higgs_ew/test_T014.py -q` AND `python -m pytest tests/constraints/ -q`. Report counts.

OUTPUT (≤14 lines): what you changed (file:line), the new total-Z-width BR formula + the corrected effective coupling, where the core helper now lives, pytest counts. End with: FIX-DONE.
