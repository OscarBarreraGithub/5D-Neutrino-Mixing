Implemented; no physics-number or coupling changes.

- Core retag: [T010.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T010.py:82) lines 82-95, 457-479; [T011.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T011.py:80) lines 80-93, 637-659.
- Split: minimal-complete uses `minimal_rs_tree_veto_ready=True`; custodial/BKT/top-partner are `*_deferred` diagnostics, not top-level `needs_human_physics`.
- Numeric invariance asserted in [test_rs_ew_phase6a_zbb_fermion_mixing.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/test_rs_ew_phase6a_zbb_fermion_mixing.py:271): T010 predicted `0.21328283177386856`, ratio `4.470421078527965`; T011 predicted `0.9383244719075848`, ratio `0.08826504462857736`.
- Routing gates: failing T010 lands in `excluded_by_rigorous`; SM/universal passing point survives; deferred custodial note stays rigorous.
- Test updates: Zbb tests plus harness test at [test_full_catalog_scan_harness.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/test_full_catalog_scan_harness.py:138).
- Full pytest: `1716 passed, 1 skipped in 681.60s`.
- Also added test-env compatibility only: [conftest.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/conftest.py:8), [test_quark_fit.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/test_quark_fit.py:323).
- No commit made; reviewer approvals still required.

ZBB-RETAG-AUTHOR-DONE