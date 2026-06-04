K001: `Constraint.evaluate` -> `physics_adapters.deltaf2.epsilon_k_*` -> `quarkConstraints.deltaf2.evaluate_epsilon_k_with_running`.
K002: `Constraint.evaluate` -> `physics_adapters.deltaf2.delta_mk_*` -> `quarkConstraints.deltaf2.evaluate_delta_mk_with_running`.
B001: `Constraint.evaluate` -> `physics_adapters.deltaf2.bd_mixing_*` -> `quarkConstraints.deltaf2.evaluate_bd_mixing_with_running`.
B003: `Constraint.evaluate` -> `physics_adapters.deltaf2.bs_mixing_*` -> `quarkConstraints.deltaf2.evaluate_bs_mixing_with_running`.
C001: `Constraint.evaluate` -> `physics_adapters.deltaf2.d0_mixing_*` -> `quarkConstraints.deltaf2.evaluate_d0_mixing_with_running`.
None of the five anchor constraints import or call `quarkConstraints.paper_0710_1869`.
Actual runtime path values: `eps_K^SM=2.161e-3`, `eps_K^exp=2.228e-3`; `Delta m_K=3.4839159e-15 GeV`.
Actual runtime path values: `Delta m_d^SM=3.6e-13 GeV`, `Delta m_d^exp=3.334e-13 GeV`; `Delta m_s^SM=1.17e-11 GeV`, `Delta m_s^exp=1.1693794e-11 GeV`.
Actual runtime path values: `x_D=0.405%`, `Delta m_D=6.562e-15 GeV`; these match the catalog/PDG-scale anchors.
Factor-2 source: Q1 matrix/operator normalization, not the `M12=<H_eff>/(2m)` step.
`deltaf2.py` uses a pre-M12 Q1 kernel `(2/3) f^2 m B`; paper mode uses `<Q1>=(8/3)f^2m^2B` then `/2m`, giving `(4/3)f^2mB`, exactly 2x for K/Bd/Bs/D.
The Wilson Q1 matching coefficient is the same `/6` in both paths, so the two Q1 kernels must not be interchanged.
No live constraint mixes `deltaf2.py` Wilson/M12 with paper-mode hadronic bundles, or vice versa.
Verdict: convention/isolated path difference for the live constraints; document “paper-mode Q1 hadronic normalization is not interchangeable with repo `deltaf2.py` kernels.”
R5-CONVENTION-OK