Plan executed: inspect 4a/3a, implement builder/couplings/Wilsons, add gate tests, run full suite.

Changed: [rs_ew_couplings.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rs_ew_couplings.py), [rs_semileptonic_wilsons.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/rs_semileptonic_wilsons.py), [rs_ew_builder.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/rs_ew_builder.py), [point_builder.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/point_builder.py), new [test_rs_ew_phase4a.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/test_rs_ew_phase4a.py).

API added: `RSLeptonMassBasisCouplings`, optional `lepton_yukawa_result` / `lepton_sweep_inputs`, `z_delta_g_L/R_e`, `z_delta_g_L_nu`, full charged/nunu contacts, `lfv_llqq`, `b_to_s_nunu`, `s_to_d_nunu`.

Asserted samples: `z_delta_g_L_e[0,0]=8.5400871097628e-06`, `z_delta_g_L_nu[0,0]=-1.588676075184686e-05`, `lfv s_to_d e_mu C9=0j`, `s_to_d_nunu X_L[0,0]=0.000820388142621831-0.00011326866594393457j`.

Gate checks: universal-c lepton max Z `0.0`; family-universal offdiag LFV `0`; 3a consistency `b_to_s mu C9=0.3046272275031843+0.005803160503768459j`; Majorana/Dirac nunu diff `<1e-18`.

Tests: 4a `6 passed`; 3a+4a `13 passed`; full `1652 passed, 1 skipped`.

Review: Codex APPROVE; independent sub-agent APPROVE. Opus-specific model/CLI was not available in this environment, so I did not label that review as Opus.

P4A-DONE