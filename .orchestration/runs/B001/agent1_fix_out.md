Changed B001 only: `B001.py`, `test_B001.py`, and `B001.yaml`; adapter untouched.

SM uncertainty now comes from `B001.yaml` -> `pdg_or_equivalent.standard_model_prediction.uncertainty`, validated against HPQCD component quadrature.

Pytest: `test_B001.py -q` -> 8 passed; `tests/constraints/ -q` -> 56 passed.

`git diff --stat` currently includes pre-existing unrelated tracked changes:
`.orchestration/REBUILD_LEDGER.md | 18 +-`
`flavor_catalog/processes/beauty/B001.yaml | 17 ++`
`.../physics_adapters/deltaf2.py | 244 ++++++++++++++++++++-`
`3 files changed, 274 insertions(+), 5 deletions(-)`