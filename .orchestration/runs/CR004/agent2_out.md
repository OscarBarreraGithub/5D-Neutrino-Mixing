1. NIT: ΔF=2 amplitude check is N/A; CR004 is a collider mass lower-bound recast, and verdict uses `m_limit / m_pred` with pass iff `m_pred >= 1.57 TeV`, not any `M12` component. `CR004.py:412-417`, `collider_resonance.py:166-168`.

2. NIT: QCD running check is N/A; no Wilson coefficients enter CR004, and it imports/uses only `collider_resonance`, not `deltaf2` or non-running ΔF=2 paths. Running effect: none. `CR004.py:46-52`, `CR004.py:417`.

3. NIT: Budget is physically right only for the documented pure `B -> H b` proxy: active limit is `1570 GeV = 1.57 TeV`; it is tighter by `10 GeV` vs pure `Wt` (`1560`) and `30 GeV` vs pure `Zb` (`1540`) if misread as branching-general. Correct full recast is `sigma(pp->B Bbar)*BR_i BR_j`/HEPData contour. `CR004.py:74-79`, `CR004.py:313-320`, `CR004.py:456-461`.

4. NIT: Anchor numbers match YAML/snapshots: `bH 1570`, `bZ 1540`, `tW 1560 GeV`, all 95% CL PDG/CMS entries. `CR004.yaml:83-128`, `pdg2025_bprime_listing.txt:11-13`, `factcheck_collider_rs.md:82-84`.

5. NIT: Units are consistent: YAML stores GeV, code converts to TeV via `/1000`, core compares TeV-to-TeV; test checks `1800 GeV -> 1.8 TeV`, ratio `1.57/1.8`. `CR004.py:103-112`, `CR004.py:327-335`, `test_CR004.py:114-119`.

6. NIT: Metadata wording `observable = "m(B pair -> tW/bZ/bH)"` can be read as pair invariant mass; correct physics is individual bottom-partner mass `m_B`. Verdict/notes use `m_B = M_KK` correctly. `CR004.py:362`, `CR004.py:421-424`, `CR004.py:456-459`.

PHYSICS-OK