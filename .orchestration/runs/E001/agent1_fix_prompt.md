You are agent1 fixing E001 (electron EDM). Code review .orchestration/runs/E001/agent3_out.md = CODE-NEEDS-FIXES (BLOCKER + SHOULD-FIX); physics OK. Repo: /n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing.
FIX:
1. (BLOCKER) The anchor loader accepted a SUPERSEDED limit block (superseded_acme_2018) for the experimental value. Assert experimental.block_key == "canonical_limit" (the current JILA limit) and raise AnchorError on any other block; add a monkeypatch test forcing the superseded block and asserting it's rejected.
2. (SHOULD-FIX) The test's expected value is a hardcoded manual HBARC multiply; recompute via the core quarkConstraints.edm (edm_e_cm_from_cp_odd_dipole / evaluate_charged_lepton_edm) instead.
Keep the NP NEEDS-HUMAN-PHYSICS proxy (RS one-loop CP-odd dipole not on ParameterPoint). Run tests/constraints/ -q green. OUTPUT <=8 lines: changes, pytest counts.
