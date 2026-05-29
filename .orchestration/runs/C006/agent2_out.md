1. NIT - ΔF=2 amplitude check is not applicable: C006 is a ΔC=1 LFV branching-ratio bound, not Δm/CP mixing; verdict uses `BR_NP` magnitude only, no `M12` real/imag part. `C006.py:500-524`; formula in `rare_charm_lfv_dilepton.py:357-373`.

2. NIT - QCD-running check is not applicable to this pure V/A semileptonic LFV proxy: no Delta-F=2 `*_with_running` path is used, and C006 correctly routes through `rare_charm_lfv_dilepton`, not non-running `deltaf2`. For C9/C10 V/A semileptonic self-running, expected LO QCD effect is 0. `C006.py:482-487`; `rare_charm_lfv_dilepton.py:232-290`.

3. NIT - Budget is defensible: LFV SM rate is zero for catalog purposes, so pure-NP bound is the observed 90% CL upper limit, not a central residual. Code uses `budget = 1.3e-8` from YAML and vetoes `predicted/budget > 1`. `C006.py:140-143`, `C006.py:516-524`; YAML `C006.yaml:70-82`.

4. NIT - Anchor numbers match snapshots: PDG/LHCb current limit `1.3e-8` at 90% CL; previous Belle `2.6e-7`, BABAR `3.3e-7`. YAML `C006.yaml:76-90,106-121`; snapshots `pdg2026_d0_emu_pdgLive.txt:20-21`, `lhcb2016_arxiv1512_00322.txt:20-21`, `belle2010_arxiv1003_2345.txt:19-22`, `babar2012_arxiv1206_5419.txt:20-23`.

5. NIT - Stale catalog prose: YAML/TeX still say code coverage “NO” / no implementation and “future use” despite live C006 code. This is documentation drift, not a physics-formula blocker. `C006.yaml:139-153`; `C006.tex:75-96`.

PHYSICS-OK