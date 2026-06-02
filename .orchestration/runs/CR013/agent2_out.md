1. NIT: ΔF=2 amplitude/QCD-running/bag-budget checks are N/A for CR013; code uses a collider mass-lower-bound path, not `M12` or Wilson running. Correct physics here is `m_GKK >= 4.8 TeV`; implemented at `CR013.py:551-557`, core pass test at `quarkConstraints/collider_resonance.py:166-168`.

2. NIT: Adapter docstring still says it “adds only the CR007-specific mass proxy,” but CR013 diphoton helpers are now appended. No numeric effect; update wording to include CR013. `kk_graviton_resonance.py:3-5`, helpers at `kk_graviton_resonance.py:240-300`.

3. NIT: Diagnostic key `exact_ktilde_mass_limits_tev` includes all exact benchmark rows, not only active `ktilde=0.1` rows; it can include `0.05`/`0.01` ATLAS rows. Verdict still uses correct active `4.8 TeV` row. `CR013.py:589-593`; active selection checked at `CR013.py:414-429`.

PHYSICS-OK