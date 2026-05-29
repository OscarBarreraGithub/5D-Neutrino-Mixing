1. NIT: ΔF=2 amplitude/QCD-running checks are N/A for CR011: no \(M_{12}\), CP phase, Wilson coefficients, or `*_with_running` path enters the verdict; code only records/optionally compares a fiducial VBS cross section (`CR011.py:349-352`, `vbs_longitudinal.py:127-174`).

2. NIT: Budget/anchors are physics-consistent: active observed ATLAS 2025 limit is `0.45 fb` at `95% CL`, expected `0.70 fb`; CMS comparator is `1.17 fb`, expected `0.88 fb`, matching YAML and snapshots (`CR011.yaml:98-143`, `atlas_2025...txt:27-29`, `cms_2020...txt:24-26`).

3. NIT: Severity and stub behavior are appropriate: `Severity.INFO`, non-vetoing by contract, no RS \(M_{KK}\), EFT/aQGC, or resonance recast is inferred; dual SM/EFT and RS matching `NEEDS-HUMAN-PHYSICS` flags are present (`CR011.py:338-440`, `vbs_longitudinal.py:23-37`).

4. NIT: Catalog prose notation is mildly misleading: `pp -> W_L W_L W_L W_L` suggests four-W final-state production; correct shorthand is \(V_LV_L \to V_LV_L\) embedded in `pp -> jj W_L^± W_L^±`. Code observable is correct (`CR011.yaml:4`, `CR011.tex:15-18`, `CR011.py:340`).

PHYSICS-OK