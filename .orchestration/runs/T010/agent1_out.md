Plan: studied scaffold/examples/T010 YAML, built shared Z-pole machinery, implemented T010, added tests, ran suite.

Physics: `R_q = Gamma_q/Gamma_had`, `A_q = (|g_L|^2-|g_R|^2)/(|g_L|^2+|g_R|^2)` from effective Zff couplings. Sources: `T010.yaml` plus its LEP/SLC snapshot.

Reused scaffold `load_anchor`, registry, `ParameterPoint`; built [quarkConstraints/zpole.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/zpole.py), [physics_adapters/zpole.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/zpole.py), [T010.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T010.py), and [test_T010.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/top_higgs_ew/test_T010.py).

Budget: `R_b` `0.0006726812`, `A_b` `0.02000025`, from YAML experimental uncertainties plus LEP/SLC SM-fit uncertainties.

SM validation: `R_b = 0.2156200000`, `A_b = 0.935535`, `A_FB^b = 0.103279`.

RS NP: documented `delta g_b ~ (m_Z/M_KK)^2 Delta overlap` proxy, flagged `NEEDS-HUMAN-PHYSICS` for missing EW KK/Z/Z'/custodial inputs.

Tests: `test_T010.py` 10 passed; `python -m pytest tests/constraints/ -q` 234 passed.

Note: worktree also contains unrelated pre-existing changes/untracked files outside T010.