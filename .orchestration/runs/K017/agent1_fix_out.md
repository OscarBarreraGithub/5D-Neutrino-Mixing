Changed `tests/constraints/primary/kaon/test_K017.py`: core cross-check now builds inputs from YAML SM anchor and calls `quarkConstraints.leptonic_tree.evaluate_leptonic_lfu_ratio(...)` directly.

`test_K017.py -q`: 10 passed.
`tests/constraints/ -q`: 724 passed.