1. BLOCKER `flavor_catalog_constraints/primary/beauty/B026.py:263`: mismatched `load_anchor` results are accepted. Probe returned `mismatch_probe_no_failure wrong_block 0.281 0.254`. Fix: after each `load_anchor`, assert `block_key` is the expected candidate and raise `AnchorError` on mismatch.
2. SHOULD-FIX `tests/constraints/primary/beauty/test_B026.py:126`: tests cover missing anchor and bad direct `load_anchor` candidates, but not a monkeypatched mismatched anchor returned to `_load_b026_anchor`. Add that regression test.

Isolation/contract/determinism otherwise OK: B026 imports physics only through `physics_adapters.semileptonic_lfu`; no tracked diffs in B025, `quarkConstraints/semileptonic_lfu.py`, semileptonic adapter, or `deltaf2.py`; missing couplings return a non-crashing HARD result; repeated evaluation is equal and did not mutate inputs.

Numerical cross-check: exp=0.281, SM=0.254, budget=0.0390830459735946. Complex NP shift `(0.09562630065761932+0.01652182157958868j)` gives constraint=core=manual `0.304970170165376`. Safe point `0.2800702727183877`, ratio `0.023788506203955572`, pass=True. Excluded point `0.36592125940580106`, ratio `2.1728413763649788`, pass=False.

Pytest: `python -m pytest tests/constraints/ -q` -> `678 passed in 15.39s`. Registry smoke passed: `B026 beauty HARD`.

CODE-NEEDS-FIXES