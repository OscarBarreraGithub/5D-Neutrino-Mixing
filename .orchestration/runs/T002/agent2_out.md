1. NIT: ΔF=2 amplitude/running checks are not applicable to T002; verdict uses `BR(t -> u Z)` from `Gamma ∝ |X_L|^2+|X_R|^2`, not `M12` or CP-imaginary parts. Correct for a top decay. `T002.py:353`, `top_fcnc.py:202`.

2. NIT: No QCD `*_with_running` path is used, and none is physically required for this top-FCNC branching-ratio constraint. The implementation correctly routes through `top_fcnc`, not `deltaf2`. Running effect: N/A. `T002.py:48`, `top_fcnc.py:486`.

3. NIT: Budget is conservative but slightly chirality-tight: active limit is `6.2e-5` LH, while YAML also has RH `6.6e-5`; pure RH points in `6.2e-5 < BR < 6.6e-5` would fail under the stronger LH benchmark. This is documented as policy, not a physics blocker. `T002.py:303`, `T002.yaml:119`, `T002.yaml:129`.

4. NIT: Anchor values match snapshots: PDG/ATLAS `tZu_L < 6.2e-5`, `tZu_R < 6.6e-5`; CMS `<0.022% = 2.2e-4`; SM `8e-17`. Units are dimensionless branching fractions; partial widths are in GeV. `pdg_2025_top_fcnc_tzq.txt:15`, `:17`, `:19`; `aguilar_saavedra_hepph0409342.txt:30`.

5. NIT: Severity/diagnostics are appropriate: SM is negligible vs limit (`8e-17 / 6.2e-5 = 1.3e-12`), verdict is pure NP over collider bound, and RS top-Z proxy is flagged `NEEDS-HUMAN-PHYSICS`. `T002.py:360`, `T002.py:389`, `T002.py:399`.

PHYSICS-OK