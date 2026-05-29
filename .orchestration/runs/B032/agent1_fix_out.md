Changed B032 to virtualize selected list entries through scaffold `load_anchor`, with loud `block_key` mismatch checks; updated B032 tests to spy and force mismatch.

Pytest: `tests/constraints/primary/beauty/test_B032.py -q` -> `9 passed`.
Pytest: `tests/constraints/ -q` -> `770 passed`.