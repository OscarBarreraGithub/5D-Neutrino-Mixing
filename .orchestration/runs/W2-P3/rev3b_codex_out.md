1. BLOCKER: none found; adapter rewire matches plan at [T010.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T010.py:499), [T011.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T011.py:676), [T012.py](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/flavor_catalog_constraints/primary/top_higgs_ew/T012.py:364).
2. SHOULD-FIX: none found; T010/T011 keep PARTIAL/NEEDS-HUMAN, T012 removes proxy NEEDS-HUMAN on rigorous path.
3. NIT: tracked code isolation OK; untracked helper plus two untracked `.orchestration` notes exist, no other constraint/core diffs.
4. Zbb recompute: `dgLb=+2.580074998e-3`, `dgRb=-1.636688349e-3`; shifted `gL=-0.4202532583`, `gR=0.0755299783`.
5. 3a pin recompute: `dgRb=-1.660422757e-3` at its identity `c_d=(a_ref,a_ref,0.20)`, 3 TeV point; not inconsistent with 3b sample.
6. Zcc recompute: `dgLc=-4.269875187e-5`, `dgRc=+1.969187886e-3`; shifted `gL=0.3456239679`, `gR=-0.1523641454`.
7. SM-limit: max `|z_delta_g|=0`; T010/T011/T012 all `predicted==sm_prediction`.
8. Degradation: old-style `build_from_quark_couplings` gives all three `passes=True`, `predicted=None`, `ratio=None`, `evaluated=False`, `missing_extra=rs_ew_couplings`.
9. Test-count verdict: legitimate `-6`; proxy pass/fail parametrizations and proxy mass-scaling tests were replaced by rigorous-path, absent-extra, and SM-limit tests at `test_T010.py:241/257/295`, `test_T011.py:286/316/352`, `test_T012.py:194/226/262`.
10. Numeric result fields are real floats/None; complex values stay in diagnostics.
11. `pytest tests/ -q`: `1646 passed, 1 skipped in 825.45s`.
12. P3B-OK.