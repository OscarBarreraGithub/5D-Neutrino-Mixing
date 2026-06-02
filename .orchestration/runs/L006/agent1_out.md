Plan: studied scaffold/K001/K004/L006 YAML, added L006 proxy adapter + constraint, added focused tests, ran suite.

Physics: `P(M -> Mbar)` is treated as pure NP with SM `0.0`; proxy maps caller-supplied `G_C/G_F` via `P = P_limit * |g/g_limit|^2`, or accepts direct probability.

Built: `flavor_catalog_constraints/physics_adapters/muonium_conversion.py`; constraint imports only this adapter.

Budget: MACS/PSI PDG anchor from `L006.yaml`: `P < 8.3e-11` at 90% C.L.; coupling calibration `G_C/G_F < 0.0030`.

Validation: independent recomputation for `g=1.5e-4+2.0e-4j` gives `P=5.763888888888889e-13`, ratio `0.006944444444444444`.

Gap: full RS Delta L=2 four-lepton matching is not on `ParameterPoint`; diagnostics flag `NEEDS-HUMAN-PHYSICS`.

Files: `L006.py`, `muonium_conversion.py`, `test_L006.py`.

Tests: `test_L006.py` 10 passed; full `python -m pytest tests/constraints/ -q` 937 passed.