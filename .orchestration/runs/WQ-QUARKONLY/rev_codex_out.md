Findings: no blocking defects.
Allowlist verdict: keep all current IN ids; keep EW002 OUT because it consumes lepton-dependent `rs_charged_current` G_F/e,mu epsilons.
Borderlines: EW001 IN (mass-only oblique proxy); EW003 IN as SOFT data-only scalar with optional CC diagnostics absent in quark-only; T010/T011 IN and read quark `z_delta_g_L/R_d[2,2]` b entries.
Mode: quark-only skips lepton draw/Yukawa/perturbativity; builds with CC off, fermion-KK on, Higgs-Yukawas off, `lepton_yukawa_result=None`; evaluates only 37 allowlist ids; EW002/non-allowlist ids absent.
Full mode: default config hash surface preserved and full path still uses lepton build + `registry.evaluate_all`; I see no full-mode behavior change.
Determinism: independent 2-run real check matched params, constraints/ratios, veto lists, allowlist, and lepton-sector tag.
Tally/provenance: per-id points/evaluated/active/failed/vetoed serialized; mode=`quark_only`, lepton sector=`dropped (not rigorous)`.
Smoke: 999/1000 evaluated, 0 constraint exceptions, 37/point, post-cache 0.2786 s/eval; ΔF=2 M_KK relaxation covered/passed.
Dominant smoke vetoes are honest: proxy EW001/B011/B012/B013 vs rigorous B003/B004/K001; T010/T011 failures are partial/hard-coverage-gap, not hidden rigorous vetoes.
Pytest: `python -m pytest tests/ -q` passed: 1711 passed, 1 skipped in 769.93s.
WQ-REVIEW: APPROVE