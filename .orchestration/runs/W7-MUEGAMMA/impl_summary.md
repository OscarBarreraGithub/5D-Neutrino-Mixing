W7-MUEGAMMA implementation summary

Files changed:
- flavor_catalog_constraints/physics_adapters/lepton.py
- flavor_catalog_constraints/point_builder.py
- flavor_catalog_constraints/primary/charged_lepton/L001.py
- flavor_catalog_constraints/rs_ew_builder.py
- scripts/run_full_catalog_scan.py
- tests/constraints/primary/charged_lepton/test_L001.py
- tests/constraints/test_contract.py
- tests/test_full_catalog_scan_harness.py
- tests/test_mu_to_e_gamma.py
- tests/test_lmfv_lepton_parameters.py

Implemented:
- Added read-only LMFVLeptonParameters carrier and YukawaResult builder.
- Wired lepton_lmfv_parameters into RS-EW extras and declared point extras.
- Switched L001 to require lepton_lmfv_parameters and use C=0.02 only as dipole diagnostic.
- Preserved BR_NP <= br_limit as the only L001 hard-veto gate.
- Added lepton_lmfv_parameters to quark-only forbidden extras without changing ScanConfig or config hashing.

Tests added/updated:
- L001 carrier oracle, absent/invalid extra behavior, BR-pinned pass/fail independence, no proxy flags, malformed carrier handling.
- Raw/core mu->e gamma benchmark oracle.
- LMFV carrier immutability, finiteness/unitarity/consistency guards, and builder physical-M_KK convention split.
- Point-builder declared key coverage.
- Quark-only forbidden-extra coverage, exact serialized-row regression, config hash pins, and nonperturbative lepton skip-before-L001 regression.

Verification:
- Targeted pytest: 65 passed.
- Full pytest with requested LD_LIBRARY_PATH export: 1734 passed, 1 skipped in 922.90s.
- BR oracle: BR_NP = 1.5508276601368708e-10; BR_NP / br_limit = 1033.8851067579139; passes = False.
- Quark-only config hashes unchanged: full=45e21a07585f7489, quark_only=d96cb734f724aedb.
- L001.tex and L001.yaml untouched; sidecar diff empty.

IMPL-READY
