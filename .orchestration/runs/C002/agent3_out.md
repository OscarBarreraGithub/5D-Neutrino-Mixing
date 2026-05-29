1. BLOCKER quarkConstraints/rare_kaon_snd.py:61: isolation fails; `git diff quarkConstraints/` is non-empty, and worktree also has unrelated B002/K005 files. Fix: isolate C002 to its constraint file, test, and append-only deltaf2 adapter wrapper.
2. BLOCKER tests/constraints/primary/charm/test_C002.py:83: CP oracle reuses the same adapter called by C002.py:360, so it is not independent. Fix: compute expected complex M12 via `quarkConstraints.deltaf2._evolve_wilsons` + `compute_m12_np`, with direct `evaluate_d0_mixing_with_running` abs check.
3. SHOULD-FIX flavor_catalog_constraints/primary/charm/C002.py:195: C002 q/phi anchors use `load_pdg_block`/`find_block`, not `load_anchor` as requested. Fix: route value-bearing anchors through `load_anchor` or add a specific tested justification for bespoke parsing.

Contract/probes: numeric result fields are real floats, complex M12 stays in diagnostics, missing extra returns non-crashing pass; missing key/value/bad CI all raise `AnchorError`; registry smoke `C002 charm HARD`.

Direct numerical check: safe M12=`-9.85283569583856e-17+0j` GeV, predicted=`0`, ratio=`0`, budget=`0.04500000000000004`, pass=True.
Direct numerical check: finite M12=`-8.48893839512015e-17-7.67967533192358e-17j` GeV, predicted=`0.023406508174104161`, ratio=`0.52014462609120315`, pass=True.
Direct numerical check: excluded M12=`-3.39557535804806e-14-3.07187013276943e-14j` GeV, predicted=`9.3626032696416637`, ratio=`208.05785043648123`, pass=False.
Determinism: repeated evaluate equal; input arrays unchanged; direct abs-vs-core delta `0.000e+00`.
Pytest: `tests/constraints/primary/charm/test_C002.py -q` = 7 passed; `python -m pytest tests/constraints/ -q` = 89 passed in 3.82s.

CODE-NEEDS-FIXES