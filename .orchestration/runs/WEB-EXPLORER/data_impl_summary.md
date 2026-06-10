# Scan Explorer Data Implementation Summary

- Built `flavor_catalog/website/src/content/scan_explorer.json` (107,294 bytes).
- Streamed rows: minimal=1,000,000, custodial=1,000,000.
- Included 17 constraints: EW001, T010, T011, B003, B004, B011, B012, B013, K001, CR001, CR002, CR003, CR004, CR008, CR010, CR012, CR013.
- Cross-checks against `comparison/constraint_veto_by_r_mkk.csv`:
  - EW001 custodial r=0.05 M_KK=1 TeV: json=1.0, csv=1.0, abs_diff=0.0
  - B003 custodial r=0.1 M_KK=1 TeV: json=0.638556, csv=0.638556, abs_diff=4.63e-07
  - B013 custodial r=0.25 M_KK=1 TeV: json=0.516305, csv=0.516305, abs_diff=1.09e-07
  - B004 custodial r=1.0 M_KK=2 TeV: json=0.571471, csv=0.571471, abs_diff=4.41e-07
  - B003 minimal r=0.1 M_KK=1 TeV: json=0.638556, csv=0.638556, abs_diff=4.63e-07
- Representative matrix SVD residuals:
  - low r (down-localized) seed=202608040007 M_KK=3 TeV: up=0.000e+00, down=0.000e+00
  - r = 0.25 seed=203008040006 M_KK=3 TeV: up=0.000e+00, down=0.000e+00
  - r = 1.0 (up-dominated) seed=203408040006 M_KK=3 TeV: up=0.000e+00, down=0.000e+00

DATA-READY
