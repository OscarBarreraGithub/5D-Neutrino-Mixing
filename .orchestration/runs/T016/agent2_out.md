1. NIT: ΔF=2 checks are N/A for T016: no `M12`, no CP phase observable, no Wilson/QCD-running path. Correct physics is a Z-pole LFV width, using `|δg_L|^2+|δg_R|^2`; see `quarkConstraints/zpole_lfv.py:157-167`, called from `T016.py:357-361`.

2. NIT: Budget is defensible for LFV pure-NP: SM set to `0`, HARD bound is direct `BR(Z→e±τ∓) < 5.0e-6`, not a residual. Loaded from YAML via `T016.py:278-307`; YAML/source values match `PDG/ATLAS = 5.0e-6`, `CMS = 13.8e-6`, expected `11.4e-6` at `T016.yaml:78-159`.

3. NIT: Channel mapping is physics-consistent with the requested proxy: e-tau matrix entry `(0,2)` is used, and charge-summed rate uses factor `2.0`; see `zpole_lfv_etau.py:172-180` and `zpole_lfv_etau.py:265-267`.

4. NIT: Units/severity/diagnostics are OK: dimensionless branching fraction budget, `Severity.HARD`, `sm_prediction=0.0`, and `NEEDS-HUMAN-PHYSICS` proxy warning present; see `T016.py:316-327`, `T016.py:371-428`. Effective coupling limit from the loaded bound is `sqrt(|δg_L|²+|δg_R|²) < 3.02e-3`.

PHYSICS-OK