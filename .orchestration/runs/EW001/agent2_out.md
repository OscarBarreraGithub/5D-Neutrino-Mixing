1. NIT: No physics defects found. ΔF=2 amplitude/QCD-running checks are N/A for EW001: code uses S,T χ², not M12 or Wilson running (`EW001.py:390`, `oblique_stu.py:187`).

2. NIT: Budget is defensible: YAML S=0.026±0.075, T=0.047±0.066, ρ=0.90, U=0 map to the 2D 95% χ² budget 5.991 (`EW001.yaml:85`, `EW001.py:278`).

3. NIT: RS proxy is clearly flagged and numerically sensible: ΔS=30 v²/MKK² and ΔT=πL/(2cW²)v²/MKK²; at 3.2 TeV, S≈0.178 but T≈0.423, so χ²≈83.46 fails as expected (`oblique_stu.py:17`, `EW001.py:462`).

4. NIT: Severity/units/anchors look consistent: HARD is appropriate for an observed EW precision exclusion contour; S,T,U are dimensionless, MKK/VEV in GeV, and PDG snapshot values match YAML (`base.py:55`, `EW001.py:357`, `pdg_2025_electroweak_stu.txt:17`).

PHYSICS-OK