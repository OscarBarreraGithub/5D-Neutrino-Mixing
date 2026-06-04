B003 path: `B003.py:414 "result = bs_mixing_from_wilsons_with_running("`; adapter `deltaf2.py:401 "result = _evaluate_bs_mixing_with_running(wilsons, mu_had=mu_had)"`.
Core path: `quarkConstraints/deltaf2.py:1035 "evolved = _evolve_wilsons(wilsons, mu_had=mu_had)"`; `:1036 "return evaluate_bs_mixing(evolved)"`; `:987 "m12 = compute_m12_np("`.
M12 formation: `deltaf2.py:931 "wilsons.c1_vll * me_vll"` plus `:932 "+ wilsons.c1_vrr * me_vll"`, `:933 "+ wilsons.c4_lr * me_lr4"`, `:934 "+ wilsons.c5_lr * me_lr5"`.
me/C1 quote: `deltaf2.py:910 "me_vll = (2.0 / 3.0) * fp2_mp * B_1"`; RS C1 source is `:342 "c1_vll=left * left * prefactor / 6.0,"`.
Test: calibrated SM `C1=1.645995189308404e-11` so standard `(4/3) f^2 m B` gives `Delta_m_s=1.087e-11`; through `compute_m12_np/evaluate_bs_mixing`, `Delta_m_s=2|M12|=5.435000000000000e-12`.
me_vll check: `(2/3) f^2 m B` is half the standard `M12` per unit C1, `<O1>/(2m)=(4/3) f^2 m B`; no x2 compensator appears in `compute_m12_np`.
Budget check: B003 uses full-convention anchors divided by 2: `B003.py:339 "experimental_m12_gev=experimental_delta_m_gev / 2.0"` and `:342 "sm_m12_gev=sm_delta_m_gev / 2.0"`.
Therefore budget is full while `M12^NP` is half; NP bound is 2x too loose regardless of the closed-form SM Δm path.
R5-NORM-HALF-NEEDS-X2-FIX.