Plan: read scaffold/CR007/CR010/CR013.yaml, append adapter helper, add CR013 constraint/tests, run suite.

Physics: CR013 compares `m_GKK^(1) = 3.8317059702075125 * Lambda_IR` to the RS-graviton diphoton mass bound.

Reused vs built: reused CR007 spin-2 graviton mapping and collider resonance core; added only a CR013 diphoton adapter helper.

Budget: `4.8 TeV`, `PDG2025:CR013:CMS2024_RSG_diphoton_kMPl_0p1`, CMS 2024/PDG 2025, `k/M_P = 0.1`, 95% CL.

Validation: no SM term; independent core check for `kk_ew_mass_gev=3100` gives `m_G=4.850880 TeV`, ratio `0.989511`.

NEEDS-HUMAN-PHYSICS: sigma*BR, branching surface, `ktilde`, width, acceptance, spin-0 interpretation.

Files changed: [kk_graviton_resonance.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/kk_graviton_resonance.py:236), [CR013.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/collider_rs/CR013.py:1), [test_CR013.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/collider_rs/test_CR013.py:1).

Tests: `python -m pytest tests/constraints/primary/collider_rs/test_CR013.py -q` 8 passed; `python -m pytest tests/constraints/ -q` 1010 passed.

Note: unrelated dirty worktree entries for B013/K021 and other adapters were left untouched.