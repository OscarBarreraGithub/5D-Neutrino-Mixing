1. NIT: ΔF=2 amplitude/CP check is not applicable; T001 is `BR(t -> c Z)`, and the code uses `Γ ∝ |X_L|^2+|X_R|^2`, which is the correct vector-coupling decay rate. `quarkConstraints/top_fcnc.py:202`

2. NIT: 2 GeV QCD running check is not applicable; this is an on-shell top EW decay, not a hadronic ΔF=2 Wilson evaluation. No `*_with_running` path is physically expected here. `flavor_catalog_constraints/primary/top_higgs_ew/T001.py:357`

3. NIT: Width normalization is consistent with `L=(g/2cW) Z_mu qbar gamma^mu X t`: `BR/|X|^2 = 0.4896766`, so the active `1.2e-4` limit corresponds to one-chirality `|X_ct| < 1.57e-2`. `quarkConstraints/top_fcnc.py:205`

4. NIT: Budget is defensible for a pure upper-limit constraint: code uses the strongest PDG/ATLAS 95% CL `t -> Zc` entry, `1.2e-4`; SM `1e-14` is only `8.3e-11` of the budget. `T001.py:306`, `T001.yaml:108`

5. NIT: Anchor values match snapshots: PDG `1.3e-4` LH, `1.2e-4` RH, CMS `4.9e-4`, SM `1e-14`; units are branching fractions, while computed partial widths remain GeV. `flavor_catalog/references/T001/pdg_2025_top_fcnc_tzq.txt:16`

6. NIT: Model-dependence is correctly disclosed: ATLAS limits are tensor-SMEFT benchmark limits, while the RS input is only a documented vector `Ztc` proxy flagged `NEEDS-HUMAN-PHYSICS`. `T001.py:15`, `T001.tex:67`

PHYSICS-OK