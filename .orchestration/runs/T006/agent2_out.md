1. NIT: Correct-amplitude check is N/A for T006: this is a ΔF=1 branching fraction, not Δm/CP mixing; code correctly uses Γ(t→ug) ∝ |ζ_L|²+|ζ_R|², not Re/Im M12. quarkConstraints/top_fcnc.py:228,382; T006.py:408-417.

2. NIT: QCD-running-to-μ_had=2 GeV check is N/A; T006 does not use ΔF=2 Wilsons or deltaf2.py. For this top decay the shared core evaluates at m_t with α_s(m_t)=0.108; full RS C_uG running M_KK→m_t is explicitly flagged NEEDS-HUMAN-PHYSICS. T006.py:83-89,408-414; top_fcnc.py:79,238-244.

3. NIT: Budget is defensible: SM ≈3.6e-14 vs active 95% CL limit 2.0e-5, so SM/budget=1.8e-9 and pure-NP ratio is appropriate. Active CMS limit is stronger than ATLAS 6.1e-5 by 3.05×. T006.yaml:138-148,167-180; T006.py:355,415-417.

4. NIT: Anchor numbers match snapshots: ATLAS BR<6.1e-5, |C_uG^ut|/Λ²<0.057 TeV^-2, σ(ug→t)<3.0 pb; CMS BR<2.0e-5, |κ_tug|/Λ<4.1e-3 TeV^-1; SM estimate 4.6e-12×0.0079≈3.6e-14. refs/T006/*.txt:24-33.

5. NIT: One diagnostic string is charm-specific in the shared adapter: “cg->t collider-recast” is emitted under `shared_top_fcnc_proxy_assumption` for T006; correct generic/up wording is qg→t or ug→t. No numeric effect; local T006 diagnostic is correct. physics_adapters/top_fcnc.py:373-379; T006.py:83-89,421-424.

6. NIT: Severity/units are consistent: HARD compares dimensionless NP BR to dimensionless BR limit; contextual C_uG/κ units remain diagnostic only. T006.py:374-376,453-467; tests/constraints/primary/top_higgs_ew/test_T006.py:96-136.

PHYSICS-OK