1. SHOULD-FIX: `a_S` sign ambiguity is not propagated. K010 data fix `|a_S|`, but the verdict uses only `a_S^eff = +|a_S| + Re(lambda_y7V^NP)/1e-4`; correct physics should evaluate/sign-envelope `s_a |a_S| + delta a_S`, `s_a=±1`, or make the positive-branch assumption explicit. `flavor_catalog_constraints/physics_adapters/rare_kaon_dilepton_ks.py:153`, `:157`, `:161`; `|a_S|=1.056`.

2. NIT: Component choice otherwise passes. K010 is not a ΔF=2 `M12` observable; for this CP-conserving K_S vector-amplitude proxy, using `Re(lambda_y7V^NP)` is the right analogue, and it does not reuse the K_L direct-CP imaginary-part formula. `rare_kaon_dilepton_ks.py:153-157`.

3. NIT: QCD-running rubric needs an explicit waiver. No `*_with_running` evaluator is called; diagnostics/tests state semileptonic running is not applied. For color-singlet semileptonic `Q7V/Q7A`, the documented LO multiplicative QCD factor is 1.0, so effect is 0%, but this is not the ΔF=2 `mu_had=2 GeV` path. `quarkConstraints/rare_kaon_dilepton.py:691-698`, `rare_kaon_dilepton_ks.py:183-187`, `tests/constraints/primary/kaon/test_K010.py:230-231`.

4. NIT: Budget is uncertainty-aware and matches YAML: full-region `BR=(5.8^{+2.9}_{-2.4})e-9`, direction-aware denominator `2.9e-9/2.4e-9`, not a bare central residual. This is right if the intended HARD budget is 1σ. `K010.py:376-394`, `:427-438`; `flavor_catalog/processes/kaon/K010.yaml:160-166`.

5. NIT: Anchor numbers match snapshots: partial `3.0e-9` with combined errors `1.513e-9/1.217e-9`, full `5.8e-9 +2.9/-2.4e-9`, `7` events, `0.15` background. `pdg2026_ks_pi0ee_pdgLive.txt:33-45`, `na48_2003_ks_pi0ee_arxiv.txt:16-20`.

6. NIT: Severity/units/notes are internally consistent: `Severity.HARD`, branching-fraction units, GeV only for mass diagnostics, and RS proxy is flagged `NEEDS-HUMAN-PHYSICS`. `K010.py:445-447`, `:547-553`; `rare_kaon_dilepton_ks.py:45-50`.

PHYSICS-NEEDS-FIXES