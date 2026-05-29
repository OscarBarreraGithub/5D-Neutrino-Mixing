1. BLOCKER: K001 does not inherit QCD running. `K001.py:193-197` calls `epsilon_k_from_wilsons`, which calls non-running `evaluate_epsilon_k` (`physics_adapters/deltaf2.py:69-77`) on Wilsons matched at `M_KK` (`quarkConstraints/deltaf2.py:321-326,338-340`). Correct physics: evolve to `mu_had=2 GeV` before 2-GeV matrix elements (`quarkConstraints/deltaf2.py:471-475,520-522,829-860`). Spot-check: test couplings at `tests/constraints/primary/kaon/test_K001.py:136` give ratio `0.266` unrun vs `2.845` run.

2. BLOCKER: The HARD budget is too tight as a standalone exclusion. `K001.py:101-103` uses only `|2.228e-3 - 2.161e-3| = 6.7e-5`; but BGS uncertainty is `~0.18e-3` (`K001.yaml:88-92`) and exp uncertainty is `0.011e-3` (`K001.yaml:73-75`). Correct physics: central residual is `(0.067 +/- 0.183)e-3`, so `6.7e-5` is a central diagnostic, not an uncertainty-aware HARD veto. Repo docs require a band (`docs/audits/epsilon_k_sm_decision.md:34-37,96-100`), with loose edge about `3e-4`.

3. NIT: Formula check passes. Core uses `m12_np.imag` and `abs(kappa/(sqrt(2)*Delta_m_K)*Im M12)` (`quarkConstraints/deltaf2.py:787-790`), not `|M12|` or Re(M12). K001 returns absolute prediction and signed Im diagnostic (`K001.py:199-218`).

4. NIT: Anchor values are consistent. YAML exp `2.228(11)e-3` (`K001.yaml:73-76`) matches snapshot (`pdg2026_epsilon_k.txt:7-10`); SM `2.161e-3 = 2.16(18)e-3` (`K001.yaml:87-92`) matches BGS snapshot (`bgs2020...txt:9-13`); FLAG `B_K(2 GeV)=0.5503(66)` and `hat B_K=0.7533(91)` (`K001.yaml:102-104`) match snapshot (`flag2024...txt:8-12`).

5. NIT: HARD severity is conceptually appropriate for epsilon_K because generic RS phases in Im Delta S=2 are strongly constrained; `K001.py:168-170` is fine. The problem is the budget policy, not the severity label.

6. SHOULD-FIX: `kappa_epsilon=0.94` is included (`quarkConstraints/deltaf2.py:655`), but its uncertainty and BGS long/short-distance uncertainties are not propagated; BGS snapshot gives `kappa_epsilon=0.94(2)` and grouped uncertainties (`bgs2020...txt:10-13`).

7. SHOULD-FIX: Units are dimensionally consistent: GeV inputs (`quarkConstraints/deltaf2.py:636-638`), GeV^3 matrix elements (`:719-724`), GeV^-2 Wilsons (`:332-345`), dimensionless epsilon. Caveat: this consistency assumes RG-evolved Wilsons; K001 currently violates that scale matching. Also B4/B5 3-GeV vs 2-GeV caveat is known (`:643-650`).

8. NIT: Wording is misleading. K001 says points must fit inside “central-value room” without the uncertainty-band caveat (`K001.py:16-19`), and the CP test class says “purely imaginary couplings maximize it” (`tests/test_epsilon_k_physics.py:205-207`) though the test correctly explains relative complex phases are needed (`:221-229`).

Overall verdict: PHYSICS-NEEDS-FIXES.