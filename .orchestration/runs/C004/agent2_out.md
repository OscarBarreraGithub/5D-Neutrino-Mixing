1. NIT: ΔF=2 M12 checks are inapplicable to C004; this is a ΔC=1 leptonic BR. Code correctly uses the rate amplitude `C10 - C10'` and `|A|^2`, not Re/Im M12: `quarkConstraints/rare_charm_dilepton.py:456`, `:473`.

2. NIT: No `*_with_running` path is used, but for the active C10-only leptonic proxy this is acceptable: C10 has zero LO QCD running; scalar/pseudoscalar slots are zero in the proxy. Verdict path is `C004.py:386` -> adapter `rare_charm_dilepton.py:98`.

3. NIT: Budget is defensible and from YAML, not a central residual: PDG 90% UL `2.1e-9` is used as `BR_SD` budget at `C004.py:133` and `:391`; YAML anchor is `C004.yaml:85-97`. LD VMD `2.7e-5 * 3.5e-8 = 9.45e-13`, only `4.5e-4` of budget, so not subtracting it is harmless and slightly loose.

4. NIT: Anchor numbers match snapshots: PDG `2.1e-9` 90% (`C004.yaml:85-97`), CMS `2.4e-9` 95% (`:98-113`), LHCb `3.1e-9` 90% (`:114-129`), LD minimum `3.0e-13` and VMD inputs (`:130-144`).

5. NIT: SM short-distance is set to exactly zero via `c10_sm=0.0` (`rare_charm_dilepton.py:122`) because YAML has no nonzero SM-SD anchor; this is acceptable as a negligible-SD policy, but it should not be described as a YAML-validated SM SD number.

PHYSICS-OK