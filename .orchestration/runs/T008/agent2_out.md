1. NIT: Check 1 OK/no fix: T008 is a BR constraint, not ΔF=2; no M12 real/imag choice applies. It uses Γ(t→uH)=m_t/(32π)(1-m_H²/m_t²)²(|y_L|²+|y_R|²), correct for L=-H ubar(y_L P_L+y_R P_R)t with m_u neglected. quarkConstraints/top_fcnc.py:248, flavor_catalog_constraints/primary/top_higgs_ew/T008.py:377.

2. NIT: Check 2 OK/no fix: QCD running to μ_had=2 GeV is not a physics requirement for this collider BR(t→uH) bound; the verdict correctly uses top_fcnc Higgs-scalar machinery, not any non-running ΔF=2 path. flavor_catalog_constraints/primary/top_higgs_ew/T008.py:384, quarkConstraints/top_fcnc.py:410.

3. NIT: Check 3 OK/no fix: budget is defensible for a pure-NP 95% CL collider upper limit, not a central residual: active limit = 1.9e-4, ratio = BR_NP/1.9e-4. flavor_catalog_constraints/primary/top_higgs_ew/T008.py:327, :384-386; YAML flavor_catalog/processes/top_higgs_ew/T008.yaml:76-80.

4. NIT: Check 4 OK/no fix: anchors match YAML/snapshots/runtime: PDG/CMS diphoton 1.9e-4, ATLAS 2.8e-4/2.6e-4, CMS 2025 7.2e-4/1.9e-4; snapshot lines agree. flavor_catalog/references/T008/pdg2026_top_hu_datablock.txt:18, :23-27; tests passed 9/9.

5. NIT: Check 5 OK/no fix: HARD severity and NEEDS-HUMAN-PHYSICS RS proxy are appropriate; units are dimensionless BR with GeV only in width/input diagnostics, so no GeV budget mismatch. flavor_catalog_constraints/primary/top_higgs_ew/T008.py:336-338, :418-444.

PHYSICS-OK