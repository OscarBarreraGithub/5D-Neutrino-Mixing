1. No BLOCKER — `base.py:171`: probes raised as intended: NaN ratio `ValueError`, Inf budget `ValueError`, complex `TypeError`, `passes=1` `TypeError`, bad severity `TypeError`; normal result constructs.
2. No BLOCKER — `anchors.py:271`: omitted `expected_*` is no-op; wrong `expected_block_key` raised `AnchorError`.
3. No BLOCKER — `base.py:110`: dict `raw`/`extras` mutation raised `TypeError`; reads still work; `evaluate_all(empty_point)` returned 103 results.
4. No BLOCKER — `registry.py:169`: `reset_for_tests(); discover()` probe: `103 -> 0 -> 103`, `same_ids=True`, `import_failures=0`.
5. No BLOCKER — backward compat spot-check: K001/B001/L001 loaded anchors and evaluated with old omitted-expectation style; all returned `ConstraintResult`, `passes=True`.
6. No BLOCKER — tests are real failure-path tests: `test_contract.py:164`, `210`, `243`, `255`; full suite: `1061 passed in 39.14s`.
7. NIT — `base.py:116`: raw immutability is top-level `Mapping` only; non-mapping mutable raw objects remain mutable, consistent with `raw: Any`.

SCAFFOLD-FIX-OK