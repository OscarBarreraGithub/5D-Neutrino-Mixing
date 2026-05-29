1. BLOCKER `quarkConstraints/rare_kaon_snd.py:58`: isolation fails; `git diff quarkConstraints/` is non-empty with unrelated rare-kaon edits. Fix: remove/split those changes from B002.
2. BLOCKER `tests/constraints/primary/beauty/test_B002.py:269`: required direct core cross-check is missing; the test reuses the new adapter, not a direct `quarkConstraints.deltaf2` evaluator. Fix: compare against core `evaluate_bd_mixing_with_running` plus direct core complex `M12` construction.
3. SHOULD-FIX `flavor_catalog_constraints/physics_adapters/deltaf2.py:312`: required adapter `bd_mixing_from_wilsons_with_running` still returns `MesonMixingResult`; B002 uses new wrapper at `:342` for phase. Fix: align the adapter contract or B002 source.
4. SHOULD-FIX `flavor_catalog_constraints/primary/beauty/B002.py:195`: beta/penguin values use `load_pdg_block/find_block`, not `load_anchor`. Fix: use `load_anchor(..., value_key=...)` or document the exception.
5. NIT `flavor_catalog_constraints/physics_adapters/deltaf2.py:273`: unrelated D0 complex wrapper is in the same adapter diff. Fix: split from B002.

Cross-check safe `(1e-3+1e-3j)`: `M12=0+2.9852989857064115e-15j GeV`, `phi=0.950163068823887 deg`, `predicted=0.7218829884163276`, `ratio=0.6716246168494265`, `budget=0.01769290183565696`, pass; direct agreement `0`.
Cross-check excluded `(5e-3+5e-3j)`: `M12=0+7.463247464266028e-14j GeV`, `phi=22.52011154516739 deg`, `predicted=0.925739373309186`, `ratio=12.19355509419042`, fail; direct agreement `0`.
Contract/determinism: severity `HARD`; missing extra returns non-crashing `passes=True, predicted=None, ratio=None`; numeric fields are floats, complex M12 stays in diagnostics; repeated eval equal, no coupling mutation.
Anchor/registry: B002 registry smoke passed; missing candidate, missing beta block, and bad numeric beta all raised `AnchorError`.
Pytest: `python -m pytest tests/constraints/ -q` -> `89 passed in 2.70s`.

CODE-NEEDS-FIXES