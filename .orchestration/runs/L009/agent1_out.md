Plan: studied K001/K004/L002 patterns, added tau-specific adapter, added L009 constraint, added tests, ran suite.

Physics: pure-NP `BR(tau -> 3mu)` with SM `0.0`; reuses LFV three-body dipole/contact/interference formula with `initial=tau`, `final=mu`.

Machinery: reused `quarkConstraints/lfv_three_body.py`; built append-only `lfv_three_body_tau` adapter. L002 `mu -> 3e` functions untouched.

Budget: `1.9e-8` from `L009.yaml` `primary_current_limit` (PDG 2025 tau listing / Belle II ADACHI 24R; Belle II source `arXiv:2405.07386`).

Validation: tau/mu dipole conversion factor `0.0022413535782490348`; `BR(tau->mu gamma)=1e-8` gives `BR(tau->3mu)=2.241353578249035e-11`.

NEEDS-HUMAN-PHYSICS: full RS lepton neutral-current, EW KK/Z/Z', tau dipole, and box matching remain proxy inputs and are flagged in diagnostics.

Files: [lfv_three_body_tau.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/physics_adapters/lfv_three_body_tau.py), [L009.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/charged_lepton/L009.py), [test_L009.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/tests/constraints/primary/charged_lepton/test_L009.py).

Pytest: `test_L009.py` 12 passed; charged-lepton 34 passed; full `tests/constraints/` 546 passed. Existing unrelated dirty/untracked worktree entries were left untouched.