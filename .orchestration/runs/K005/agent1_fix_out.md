Changes:
- K005 list-shaped YAML anchors now route selected value blocks through scaffold `load_anchor`; added regression coverage for those calls.
- Added the intentional-location note in the rare-kaon adapter; kept RS NP proxy and `NEEDS-HUMAN-PHYSICS`.

Pytest:
- `python -m pytest tests/constraints/primary/kaon/test_K005.py -q`: 10 passed.
- `python -m pytest tests/constraints/ -q`: 86 passed, 4 failed, all in unrelated untracked `C002` tests; I did not touch C002 due the K005-only scope.