1. NIT: ΔF=2 `M_12` real/imag check is N/A for B006; code correctly uses the ΔB=1 `B_q -> mumu` rate, `|lambda_t^d|^2(|P|^2+|S|^2)`, via `b_d` at `B006.py:442-444`, `rare_b_meson.py:119-129`, `rare_b_dilepton.py:445-494`. No wrong real/imaginary `M_12` usage found.

2. NIT: No `*_with_running` ΔF=2 evaluator is expected here. Current proxy has C9 diagnostic-only and C10-only leptonic rate; scalar/pseudoscalar Wilsons are zero, so QCD-running effect on the active verdict is 0 in this implementation. Path: `B006.py:480-486`, `rare_b_dilepton.py:559-616`.

3. NIT: Budget is defensible for a 90%-CL upper-limit veto: active limit `1.5e-10`, SM anchor `1.06e-10`, budget `4.4e-11`, with total `BR <= 1.5e-10` enforced at `B006.py:425-430`. SM sigma `0.09e-10` is loaded but not used in the ratio; a 1σ-down SM-room scale would be `5.3e-11`, but the total-limit veto is the controlling physics.

4. NIT: Anchors match snapshots/YAML: PDG/CMS `BR(B0->mumu)<1.5e-10` 90% CL, LHCb `<2.3e-10`, HFLAV ratio `<8.1e-2`, Bobeth SM `(1.06+-0.09)e-10`; see `B006.yaml:90-135` and refs `pdg_2026_bdmumu.txt:19-35`, `cms_2212_10311_bdmumu.txt:16-20`, `lhcb_2108_09283_bdmumu.txt:15-19`, `bobeth_1311_0903_bsdll_sm.txt:12-18`.

5. NIT: Severity/units/docstrings look physics-consistent: HARD upper-limit treatment, branching-fraction units, `y_d=0` time factor validated as `1.0`, and RS C9/C10 proxy is clearly flagged `NEEDS-HUMAN-PHYSICS` at `B006.py:13-21`, `B006.py:30-35`, `B006.py:561-568`.

PHYSICS-OK