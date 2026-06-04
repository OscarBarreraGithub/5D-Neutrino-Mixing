Benchmark: `B_s`, `M_KK=3000 GeV`, default `r=0.25` fit; `g_s=1.000542`, `L_bs=-0.00385793+0.000694371i`, `R_bs=5.70773e-05-2.30976e-05i`.

(1) `deltaf2.py` path: `evaluate_bs_mixing_with_running(mu=2 GeV)` gives `2|M12^NP| = 9.04057e-14 GeV`.

(2) Independent formula: Blanke et al. [arXiv:0809.1073](https://arxiv.org/abs/0809.1073), Eq. (4.33), `d -> s`, same low-scale Wilsons and repo `B_s` inputs, gives `Delta_m_Bs^NP = 5.33454e-14 GeV`.

Ratio `(1)/(2) = 1.69472`. This is not the `0.5x` failure mode; an x2 increase would be the wrong fix.

Self-consistency: `deltaf2.py` does not compute an SM box Wilson/Inami-Lim path; it only stores `DELTA_M_BS_SM`, `DELTA_M_BD_SM`, and `EPSILON_K_SM` constants for budgets.

R5-CORRECT