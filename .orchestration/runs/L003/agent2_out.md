1. SHOULD-FIX: `Severity.HARD` is used with the Mu2e projected expected UL `6.7e-17`, not an observed Al bound; `Severity.HARD` is documented as an observed experimental veto, while projections are `SOFT`. Correct physics policy: do not hard-veto present points on a future projection, or use an observed/translated bound. `flavor_catalog_constraints/primary/charged_lepton/L003.py:32`, `:314`, `:420`; `flavor_catalog_constraints/base.py:55`; YAML has SINDRUM Au `7.0e-13`, not Al, and Mu2e projection `6.7e-17` (`L003.yaml:74`, `:92`).

2. SHOULD-FIX: Mixed dipole-contact verdict uses the constructive interference envelope when only L001 BR is available, so ambiguous phase points can be hard-failed. Correct physics: with only `BR(mu->e gamma)` the relative dipole/contact phase is unknown; use explicit `A_L/A_R`, report an allowed interval, or fail only if the lower envelope exceeds the limit. `quarkConstraints/mu_e_conversion.py:315`, `:332`, `:337`, `:371`; `L003.py:423`.

3. NIT: ΔF=2 checks are inapplicable: no `M12`, no real-vs-imaginary CP choice, and no QCD `_with_running` evaluator should appear for low-energy mu-e nucleon proxies. The code uses `|left|^2+|right|^2` for the conversion amplitude. `quarkConstraints/mu_e_conversion.py:321`, `:323`, `:346`.

4. NIT: Anchor/input spot-check OK: Mu2e `6.7e-17`, SES `2.8e-17`; SINDRUM Au `7.0e-13`; L001 BR `1.5e-13`; Al overlaps `D=0.0362,Vp=0.0161,Vn=0.0173,Sp=0.0155,Sn=0.0167`, `Gamma_capt=0.7054e6 s^-1 = 4.643e-19 GeV`, matching KKO hep-ph/0203110 and YAML snapshots. `L003.yaml:74`, `:92`; `quarkConstraints/mu_e_conversion.py:139`.

5. NIT: SM=0 and units/proxy notes are otherwise consistent: CR is dimensionless, width/capture are in GeV, and scalar/vector RS matching is flagged `NEEDS-HUMAN-PHYSICS`. `quarkConstraints/mu_e_conversion.py:16`, `:348`; `L003.py:388`, `:425`.

PHYSICS-NEEDS-FIXES.