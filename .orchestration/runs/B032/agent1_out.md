Implemented B032 as an INFO-only non-vetoing stub.

- Physics/source: records HFLAV Dec-2025 `Delta A_CP(K pi) = 0.1101` from `B032.yaml`; no hadronic amplitude is computed.
- Machinery: added small adapter [charmless_b_pik.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/charmless_b_pik.py:1), constraint [B032.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/beauty/B032.py:1), and tests [test_B032.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/beauty/test_B032.py:1).
- Budget: full observed `|Delta A_CP(K pi)| = 0.1101`, sourced from `flavor_catalog/processes/beauty/B032.yaml` direct CP anchors.
- Validation: independent YAML recompute gives `Delta A_CP = 11.01%`, sigma `1.2395%`, ratio `1.0`.
- Gaps: diagnostics flag `NEEDS-HUMAN-PHYSICS` for both SM QCDF/SCET hadronics and RS `Delta B=1` penguin matching.
- Tests: `8 passed` for B032; full `python -m pytest tests/constraints/ -q` gives `769 passed`.

Unrelated working tree changes from other tasks were present; I left them untouched.