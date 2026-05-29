1. NIT: ΔF=2 amplitude/QCD-running checks are N/A for T012; this is a flavor-diagonal Z-pole observable, not an \(M_{12}\) constraint. Code uses \(|g_L|^2, |g_R|^2\) for \(R_c,A_c\), which is correct: `quarkConstraints/zpole.py:225-232`, `245-257`.

2. NIT: Correct core path: T012 routes through `zpole_charm.py` to reusable `zpole.py`, not `deltaf2.py`; no inappropriate non-running ΔF=2 path is used: `T012.py:46-52`, `T012.py:316-335`.

3. NIT: SM anchors validate: live evaluator gives \(R_c=0.172100\), \(A_c=0.66758\), matching requested \(R_c≈0.1721\), \(A_c≈0.667\). Calibration and evaluation are at `zpole_charm.py:80-102`, `T012.py:316-321`.

4. NIT: Budgets are defensible from YAML-only experimental errors: \(σ(R_c)=0.0030\), \(σ(A_c)=0.027\); SM-limit pull is \(A_c\) selected at 0.0897σ. No separate SM-fit uncertainty exists in T012 sidecar: `T012.yaml:71-101`, `T012.py:201-224`.

5. NIT: Anchors match snapshots/audit: \(R_c^0=0.1721±0.0030\), \(A_c=0.670±0.027\), \(A_{FB}^{0,c}=0.0707±0.0035\): `pdg_2025_z_boson_charm.txt:8-20`, `factcheck_top_higgs_ew.md:376-381`.

6. NIT: Severity/notes/units are physics-consistent: observables are dimensionless, HARD is appropriate for LEP/SLC Z-pole precision, and RS Zcc proxy is clearly flagged NEEDS-HUMAN-PHYSICS: `T012.py:20-34`, `zpole_charm.py:49-55`.

PHYSICS-OK