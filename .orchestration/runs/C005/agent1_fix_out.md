Correct value: `1.28111287014103257225e-17` from rare-charm core, C005 adapter, and manual recomputation.

Updated golden to `1.28111287014103257e-17`; targeted small-value checks now use `rel=1e-12, abs=0.0`.

`tests/constraints/primary/charm/test_C005.py -q`: 12 passed.

`tests/constraints/ -q`: 297 passed.