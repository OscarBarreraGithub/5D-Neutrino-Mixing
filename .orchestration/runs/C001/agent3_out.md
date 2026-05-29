1. BLOCKER tests/constraints/primary/charm/test_C001.py:190: the “direct” cross-check hides a core-vs-YAML budget mismatch. It compares `result.ratio` to `audited.abs_m12_np / expected_budget`, not to `audited.ratio_to_budget`, and `result.budget` to YAML, not `audited.budget`. Fix by either aligning/documenting the override explicitly in-test, or comparing direct ratio/budget and updating the implementation policy.

2. SHOULD-FIX repo status: isolation is not clean at worktree level. Besides C001/test/adapter, status shows unrelated B001/B003/K002 files and orchestration changes. Fix by splitting/cleaning the C001 patch before merge.

3. NIT tests/constraints/primary/charm/test_C001.py:86: anchor tests cover YAML match but not loud failure. My probe confirms `AnchorError` on missing candidate/value key; add that regression test.

Cross-check numbers: safe C001 `pred=1.14472480657973238e-16`, `ratio=3.48895094964868160e-02`, `budget=3.28099999999999987e-15`, pass=True; direct core `pred=1.14472480657973238e-16`, `ratio=3.66311938105514379e-02`, `budget=3.12500000000000009e-15`, pass=True.

Cross-check numbers: excluded C001 `pred=4.57889922631892903e-14`, `ratio=1.39558037985947241e+01`, `budget=3.28099999999999987e-15`, pass=False; direct core `pred=4.57889922631892903e-14`, `ratio=1.46524775242205720e+01`, `budget=3.12500000000000009e-15`, pass=False.

Other checks: `git diff quarkConstraints/` empty; C001 imports no other constraint; numeric result fields are real floats; complex values stay in diagnostics; missing-couplings path is non-crashing; deterministic/purity probe passed.

Registry smoke: 5 constraints registered, C001 registered as charm, 0 import failures.

Pytest: `python -m pytest tests/constraints/ -q` -> 51 passed in 3.35s.

CODE-NEEDS-FIXES