Computed with `deltaf2.py` M12 kernel and repo LL VLL running to `mu=2 GeV`.

B_s: `2|M12_SM| = 1.139745e-11 GeV` vs PDG/anchor `1.168800e-11` (ratio `0.975`).
B_d: `2|M12_SM| = 3.497365e-13 GeV` vs PDG/anchor `3.334000e-13` (ratio `1.049`).
K: short-distance box through same kernel gives `2.234895e-15 GeV` vs PDG `3.484000e-15`; K is LD/scheme-sensitive, not the clean factor-2 discriminator.

Status: B_s/B_d budgets are anchored/hardcoded (`Delta_m_exp/2`), not computed dynamically.
Same kernel + SM C1 reproduces the B anchored SM scale, not half of it.
Verdict: NP bound `|M12^NP| <= Delta_m_exp/2` is in the correct full convention.

R5-NORM-CORRECT