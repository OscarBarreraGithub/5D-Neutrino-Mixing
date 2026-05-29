1. NIT: ΔF=2 amplitude check is not applicable to C005; the verdict uses the leptonic BR amplitude `|C10 - C10' + ...|^2`, not `Re/Im M12` (`C005.py:432-439`, `rare_charm_dilepton.py:456-475`).

2. NIT: QCD-running check is not applicable; C005 is a ΔC=1 rare leptonic BR and calls the rare-charm adapter, not any `deltaf2` non-running path (`C005.py:55-60`, `C005.py:432-436`). Running effect: no ΔF=2 Wilson evolution enters this verdict.

3. NIT: Budget is defensible: `BR_SD <= 7.9e-8` uses the YAML PDG/Belle 90% CL upper limit without LD subtraction (`C005.yaml:69-96`, `C005.py:437-461`). This is conservative, not too tight.

4. NIT: Anchor numbers match snapshots: PDG/Belle `7.9e-8` at 90% CL and BaBar `1.7e-7` at 90% CL (`pdg2026_d0_ee_pdgLive.txt:17-32`, `belle2010_arxiv1003_2345.txt:21-23`, `babar2012_arxiv1206_5419.txt:21-23`).

5. NIT: Severity/units/notes are physics-consistent: HARD bound, dimensionless branching fraction from GeV inputs, electron mass gives `m_e^2/m_mu^2 = 2.339e-5` and phase-space-corrected `BR_e/BR_mu = 2.354e-5`; RS proxy is marked `NEEDS-HUMAN-PHYSICS` (`C005.py:350-351`, `C005.py:381-408`, `C005.py:462-471`).

PHYSICS-OK