1. BLOCKER: K008 is not a rigorous general `s -> d e+e-` short-distance NP rate. `quarkConstraints/rare_kaon_dilepton.py:502-504` rescales the full ISU `C_dir^e=2.4` term with one axial `Y_NP`; BMS/ISU direct CP depends on separate vector/axial electron Wilsons (`y7V`, `y7A`), and interference is vector-only. Fix with Wilson-level `y7V/y7A`, or relabel as axial-`Y` proxy only.

2. SHOULD-FIX: Mandatory running check fails as written: `K008.py:616-620` calls `klong_pi0ee_direct_cp_from_couplings`, which reaches `compute_rare_kaon_dilepton_wilsons` at `rare_kaon_dilepton.py:557-562`; no `*_with_running` evaluator or `mu_had=2 GeV` evolution is used. If semileptonic QCD running is intentionally negligible, add that explicit diagnostic; current running effect is unquantified.

3. NIT: Correct CP amplitude projection for the implemented proxy: this is not a ΔF=2 `M12` observable, and the K_L direct-CP piece uses `y_eff.imag` at `rare_kaon_dilepton.py:502-503`, not `real` or `abs(M12)`.

4. SHOULD-FIX: CPC diagnostics are inconsistent with the YAML CPC anchor. YAML/PDG gives `BR_CPC = 0.0047e-10 = 4.7e-13` at `K008.yaml:152-162`, but the computed `cpc_branching_fraction` is forced to `0` via ISU `C_CPC^e=0` at `rare_kaon_dilepton.py:516` and `K008.yaml:284-294`. Diagnostic-only, but CPC was requested.

5. NIT: Budget is defensible only for the declared direct-CP proxy, not full `BR(K_L -> pi0 e+e-)`: `K008.py:622-625` compares direct CP to the total 90% CL limit `2.8e-10`. This is loose for direct-only and not a full total-rate exclusion with unknown interference sign.

6. NIT: Anchors and SM validation look numerically consistent: PDG `<2.8e-10`, KTeV `<3.5e-10`, ISU coefficients match YAML/snapshots; SM direct `4.84e-12` lies in `3-6e-12`, constructive total `3.80e-11` lies in `8-45e-12`, and repo `Im(lambda_t)/1e-4=1.420` is within YAML `1.36±0.12`.

PHYSICS-NEEDS-FIXES