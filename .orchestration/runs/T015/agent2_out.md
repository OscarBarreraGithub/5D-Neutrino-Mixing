1. NIT (verdict OK): Amplitude/rate is magnitude-squared LFV coupling, not Re/Im; Δm/CP M12 checks are N/A for this Z-pole BR. `quarkConstraints/zpole_lfv.py:155-167` uses `|δg_L|^2+|δg_R|^2` with charge-state factor 2.

2. NIT (verdict N/A): QCD running check is inapplicable; this is not a ΔF=2 Wilson/hadronic observable. T015 correctly routes through `zpole_lfv`, not `deltaf2`: `flavor_catalog_constraints/primary/top_higgs_ew/T015.py:46-52`, `:352-356`.

3. NIT (verdict OK): Budget is the pure-NP one-sided BR limit, not a central residual: CMS observed 95% CL `1.9e-7`, SM set to `0`. Loaded/applied at `T015.py:275-299`, `:354-406`; effective coupling limit is `sqrt(|δg_L|^2+|δg_R|^2) < 5.888e-4`.

4. NIT (verdict OK): Anchor numbers match YAML/snapshots/audit: CMS `1.9 (2.0)e-7`, `138 fb^-1`; PDG/ATLAS `2.62e-7`, charge-summed. See `T015.yaml:77-90`, `:100-112`, snapshots `cms...txt:20-22`, `pdg...txt:15-21`.

5. NIT: Denominator uses default tree-level zpole total width weight `3.649563`, not ΓZ-calibrated `3.760943`; this makes the coupling limit ~1.5% tighter (`5.888e-4` vs `5.977e-4`). Given the documented proxy status, not blocking. `quarkConstraints/zpole.py:73-84`, `quarkConstraints/zpole_lfv.py:119-167`.

6. NIT (verdict OK): HARD severity and NEEDS-HUMAN-PHYSICS caveat are appropriate; missing off-diagonal coupling is explicitly non-vetoing/unevaluated. `T015.py:21-26`, `:306-307`, `:318-340`.

PHYSICS-OK