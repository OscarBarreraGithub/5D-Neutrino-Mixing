1. BLOCKER: Missing KKO units factor. `quarkConstraints/mu_e_conversion.py:367-374` uses tabulated `D,S,V` as dimensionless in a GeV width; KKO tabulates overlaps in units of `m_mu^(5/2)`, so rate needs `m_mu^5`. Current Ti pure-dipole `CR/BR = 314`, correct is `0.00414`; verdict ratios are too tight by `1/m_mu^5 = 7.59e4`.

2. BLOCKER: L005 budget is too loose vs requested/KKO Ti benchmark. `flavor_catalog/processes/charged_lepton/L005.yaml:96` and `L005.py:119-120` use `4.3e-12`; KKO states `6.1e-13`, and Bernstein-Cooper tables give Ti `6.1e-13/6.5e-13`. This hard bound is loose by `4.3e-12 / 6.1e-13 = 7.05`.

3. SHOULD-FIX: Test locks in the wrong normalization. `tests/constraints/primary/charged_lepton/test_L005.py:99-102` independently recomputes the same missing-`m_mu^5` formula, so it cannot catch the unit error. It should assert the Ti dipole benchmark `CR/BR(mu->e gamma) ~= 0.00414`.

4. NIT: DeltaF2-only checks are N/A here: no `M_12`, no CP-imaginary projection, and no QCD-running evaluator should be used for this charged-lepton conversion observable. The amplitude structure itself is the KKO coherent modulus-squared form at `mu_e_conversion.py:315-365`.

Refs: KKO https://arxiv.org/abs/hep-ph/0203110 ; Ti history https://lss.fnal.gov/archive/2013/pub/fermilab-pub-13-259-ppd.pdf

PHYSICS-NEEDS-FIXES