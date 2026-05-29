You are agent3 — CODE reviewer + numerical verifier (codex). Review code/tests of K004 = BR(K+->pi+nunu), which added NEW physics modules. Do NOT rewrite; numbered findings + run things. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.
Files: quarkConstraints/rare_kaon_snd.py (new), flavor_catalog_constraints/physics_adapters/rare_kaon.py (new adapter), flavor_catalog_constraints/primary/kaon/K004.py, tests/constraints/primary/kaon/test_K004.py. Scaffold contract in flavor_catalog_constraints/base.py (numeric ConstraintResult fields real floats; complex in diagnostics).
CHECK (with evidence):
1. ISOLATION: K004 imports no other constraint; changes limited to the new physics module + new adapter + K004.py + its test; `git diff -- flavor_catalog_constraints/physics_adapters/deltaf2.py` empty (must NOT have touched the ΔF=2 adapter); no other constraint modified. Constraint reaches physics ONLY via the new adapter.
2. CONTRACT: numeric result fields real floats; complex in diagnostics; severity HARD; missing-couplings non-crashing.
3. ANCHOR: experimental BR + SM from K004.yaml via load_anchor (NOT hardcoded); loud fail on missing/mismatch (probe).
4. NUMERICAL VERIFICATION: independently confirm the SM-limit BR (~8.5e-11) by calling the new core directly with zero NP; confirm a large-NP point pushes BR outside the budget (fails). Report ACTUAL numbers.
5. DETERMINISM: evaluate() pure. CODE QUALITY: error handling, the new module is self-contained and importable.
6. Run `python -m pytest tests/constraints/ -q`; report counts. Registry smoke.
OUTPUT (<=22 lines): numbered findings (BLOCKER/SHOULD-FIX/NIT) file:line + fix, the actual SM-limit + NP numbers, pytest counts. End: CODE-OK or CODE-NEEDS-FIXES.
