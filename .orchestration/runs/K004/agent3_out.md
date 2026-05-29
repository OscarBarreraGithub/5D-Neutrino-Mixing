1. BLOCKER [K004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/kaon/K004.py:39): anchor loading bypasses `load_anchor` (`load_pdg_block`/`find_block` at line 260). Fix by using/extending `load_anchor` for exp + SM anchors.

2. BLOCKER [B001.yaml](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog/processes/beauty/B001.yaml:119): isolation not clean; unrelated B001 sidecar diff exists, plus untracked `.orchestration/runs/K004/`. Fix by splitting/removing unrelated artifacts before merge.

3. SHOULD-FIX [test_K004.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/kaon/test_K004.py:118): tests hardcode anchor-derived budget/SM constants. Fix by deriving these from `K004.yaml`/anchor object.

Evidence: `deltaf2.py` diff empty; no other constraint code modified; `K004.py` reaches physics only via `physics_adapters.rare_kaon`.

Contract: severity `HARD`; missing couplings return `passes=True`, `predicted=None`, `ratio=None`; numeric result fields are `float`; complex values are diagnostics.

Anchor probe: missing keys raised `AnchorError`; bad units raised `AnchorError`.

Numerics: SM core/evaluate `BR=8.498092727794516e-11`, pass `True`, ratio `0.5961570616599263`. Large NP (`left_sd=5.0e-2`) `BR=2.0131000105410632e-09`, pass `False`, ratio `98.52161534451356`.

Determinism/import/registry: new modules import cleanly; registry smoke `K004 HARD kaon`; repeated `evaluate()` equal `True`.

Pytest: `python -m pytest tests/constraints/ -q` -> `65 passed in 2.52s`.

CODE-NEEDS-FIXES