1. BLOCKER: none.
2. SHOULD-FIX: none.
3. NIT: none.

Isolation: CR014 imports no other constraint; adapter change is append-only at [collider_resonance.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/collider_resonance.py:586). AST check: preexisting_changed=[], preexisting_missing=[], added=[resolve_top_philic_vector_mass_gev, top_philic_vector_four_top_prediction_from_mass_gev].

Contract: numeric result fields are real floats on evaluated paths; missing extras returns non-crashing pass with predicted=None, ratio=None; invalid mass returns non-crashing failure.

Anchor: active YAML row loads as 850.0 GeV = 0.85 TeV, width/m=0.5. Bad value_id probe raises AnchorError loudly.

Independent cross-check vs `quarkConstraints.collider_resonance`: 1.0 TeV -> ratio 0.85 pass=True; 0.9 TeV -> ratio 0.9444444444444444 pass=True; 0.8 TeV -> ratio 1.0625 pass=False. CR014 evaluate matches these numbers exactly.

Determinism: same point evaluated twice produced equal results and did not mutate `point.extras`.

Registry import smoke: `fcc.get("CR014")` returns collider_rs HARD constraint.

Pytest: `python -m pytest tests/constraints/ -q` -> 1054 passed in 34.48s.

CODE-OK