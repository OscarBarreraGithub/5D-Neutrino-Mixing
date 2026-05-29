Plan executed: studied T001/scaffold/T002.yaml, implemented T002, added tests, ran full constraints.

Physics: `BR(t -> u Z)` uses the shared top-FCNC vector convention through `physics_adapters/top_fcnc.py`; SM is negligible, so T002 is pure-NP vs collider limit.

Reused vs built: reused `quarkConstraints/top_fcnc.py` and the existing adapter; built only [T002.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T002.py) plus [test_T002.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/top_higgs_ew/test_T002.py).

Budget: `6.2e-5`, source `PDG2025:T002:tZu_left` from `flavor_catalog/processes/top_higgs_ew/T002.yaml`.

Validation: SM anchor `8e-17`; independent manual proxy check gives `BR = 4.817352363171973e-7` for the test coupling at `M_KK=3 TeV`.

Gap: RS top-Z FCNC matching remains flagged `NEEDS-HUMAN-PHYSICS` in docstring and diagnostics.

Tests: `python -m pytest tests/constraints/ -q` -> `265 passed`.