1. SHOULD-FIX tests/constraints/primary/top_higgs_ew/test_T006.py:73 - test recomputes the width with a local formula, not the underlying `quarkConstraints.top_fcnc.gluon_dipole_branching_fraction`; fix by deriving `zeta_L/R` then calling that core evaluator.

Isolation/contract OK: T006 imports no other constraint or direct `quarkConstraints`; tracked diff for T005/top_fcnc/deltaf2/quarkConstraints is empty. Missing-extra path returns non-crashing result.

Independent verifier: predicted `1.118111915877719e-4`, core `1.118111915877719e-4`, abs diff `0`; budget `2.0e-5`, ratio `5.590559579388595`.

Safe/excluded: safe left=0.4 gives pred `1.5522594927586553e-5`, ratio `0.7761297463793275`, pass; excluded left=0.5 gives pred `2.425405457435399e-5`, ratio `1.2127027287176995`, fail.

Anchor probes: missing value_id and mismatched `load_anchor` both raise `AnchorError`. Registry smoke: `T006 top_higgs_ew HARD BR(t -> u g)`.

Pytest: `python -m pytest tests/constraints/ -q` -> `678 passed in 14.95s`.

CODE-NEEDS-FIXES