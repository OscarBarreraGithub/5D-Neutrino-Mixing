1. BLOCKER: T011 does not use T011 as the physics source of truth: `T011.yaml` is still a `CLOSED-AS-MERGED` stub pointing to T010, and `T011.py` explicitly falls back to T010 anchors. Correct build should carry standalone T011 anchors/snapshots for `A_FB^{0,b}` and `A_b`. `flavor_catalog/processes/top_higgs_ew/T011.yaml:40`, `flavor_catalog_constraints/primary/top_higgs_ew/T011.py:170`, numbers currently load from T010: `A_FB=0.0992±0.0016`, `A_b=0.923±0.020`.

2. NIT: Correct-amplitude check is N/A for T011: this is not a ΔF=2 mass-difference or CP-observable constraint, so there is no `M_12` magnitude/imaginary-part choice. The implemented physics compares Z-pole asymmetry shifts, using `A_b=(|g_L|^2-|g_R|^2)/(|g_L|^2+|g_R|^2)` and `A_FB^0,b=3/4 A_e A_b`, which is correct. `quarkConstraints/zpole.py:225`, `quarkConstraints/zpole.py:235`.

3. NIT: QCD-running check is N/A for T011: no Wilson coefficients or hadronic matrix elements are involved, so no `*_with_running` evaluator should be required. The code correctly routes through the Z-pole evaluator/proxy path, not `deltaf2`. `flavor_catalog_constraints/primary/top_higgs_ew/T011.py:633`, `flavor_catalog_constraints/physics_adapters/zpole.py:148`.

4. NIT: Budget arithmetic is uncertainty-aware and physically defensible as a loose hard-veto envelope, but it is sourced from T010 because of finding 1. Current budgets are `A_FB`: `|0.0992-0.1037| + sqrt(0.0016^2+0.0008^2)=0.00628885`; `A_b`: `|0.923-0.9346| + sqrt(0.020^2+0.0001^2)=0.03160025`. `flavor_catalog_constraints/primary/top_higgs_ew/T011.py:454`.

5. NIT: Anchor spot-checks match the local PDG/LEP-SLC snapshots once T010 fallback is accepted: formula SM gives `A_b=0.935535`, `A_FB^0,b=0.103279`; LEP/SLC SM-fit snapshot gives `0.9346±0.0001`, `0.1037±0.0008`, `A_FB` pull `2.8σ`. `flavor_catalog/references/T010/lepslc_2006_z_resonance.txt:13`.

PHYSICS-NEEDS-FIXES