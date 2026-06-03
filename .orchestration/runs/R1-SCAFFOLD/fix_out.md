Plan: patch additive guards, pin them with tests, verify full constraint suite.

Changes: [base.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/base.py:110) freezes `ParameterPoint` mappings; [base.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/base.py:171) enforces bool/severity/finite real numerics.
Changes: [anchors.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/anchors.py:181) adds `value_id`/CL fields; [anchors.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/anchors.py:271) adds optional `expected_*` guards.
Changes: [registry.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/registry.py:168) re-registers cached modules after reset; [TEMPLATE.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/TEMPLATE.py:62) demonstrates anchor validation.
Changes: [test_contract.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/test_contract.py:50) adds reset, anchor mismatch, immutability, finite numeric, passes/severity tests.
Probe: NaN ratio -> `ValueError`; wrong block_key -> `AnchorError`; extras mutation -> `TypeError`; reset->discover -> `103->103`.
Pytest: `python -m pytest tests/constraints/test_contract.py -q` -> `20 passed`; `python -m pytest tests/constraints/ -q` -> `1061 passed`.
FIX-DONE