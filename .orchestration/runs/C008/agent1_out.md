Plan: studied scaffold, K001/K004, C006/C007, and C008 YAML; implemented C008; added tests; ran suite.

Physics: pure-NP `BR(D+ -> pi+ e±mu∓)`, SM = `0.0`; full-q2 vector/axial LFV `D -> pi` proxy using `C9+C9p` and `C10+C10p`.

Reused/built: reused C007 form factors and C006 LFV e-mu Wilson proxy; built new LFV semileptonic core + adapter. C004-C007 untouched.

Budget: C008 YAML PDG/LHCb 90% CL limits `2.1e-7` and `2.2e-7`; HARD ratio uses the stricter/max charge-mode saturation.

Validation: independent manual integration test gives `1.2409725765054486e-11` for the test proxy point.

NEEDS-HUMAN-PHYSICS: LFV lepton coupling remains an explicit `lepton_mass_basis_couplings` proxy; no full EW KK/Z/Z', scalar/tensor, resonance, or LHCb acceptance matching.

Files: added `quarkConstraints/rare_charm_lfv_semileptonic.py`, adapter, `C008.py`, and `test_C008.py`.

Tests: `test_C008.py` -> 12 passed; `python -m pytest tests/constraints/ -q` -> 452 passed.