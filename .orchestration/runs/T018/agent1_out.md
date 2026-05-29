Plan completed: studied scaffold/examples/YAML, built core + adapter, added T018 constraint and tests, then ran the suite.

Physics: `BR(h -> mu tau) = m_h (|Y_mu_tau|^2 + |Y_tau_mu|^2)/(8 pi Gamma_h)`, SM LFV rate `0`. Source limit loaded from `T018.yaml`.

Built new reusable machinery in [higgs_lfv.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/quarkConstraints/higgs_lfv.py), wrapped by [physics_adapters/higgs_lfv.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/higgs_lfv.py). Constraint imports only the adapter.

Budget: `BR < 1.5e-3`, CMS2021 `CMS2021:T018:mutau_limit`; ATLAS/PDG context is also loaded from YAML. Derived effective Yukawa norm limit: `1.1068120958769132e-3`.

NEEDS-HUMAN-PHYSICS is documented in code and diagnostics: the off-diagonal Higgs-Yukawa is a proxy because RS Higgs-Yukawa FCNC inputs are not on `ParameterPoint`.

Files added: [T018.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T018.py), [test_T018.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/top_higgs_ew/test_T018.py), plus the core and adapter above.

Tests: `test_T018.py` 12 passed; `python -m pytest tests/constraints/ -q` 558 passed. Unrelated dirty/untracked files were left untouched.