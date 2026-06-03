Updated [.orchestration/runs/W2-P5/plan.md](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/.orchestration/runs/W2-P5/plan.md).

- B025 recorded as PARTIAL; added consensus §4 correction: FULL 25 / PARTIAL +1, aligned with B013/B014 exclusive handling.
- Renamed W-mixing coefficient to `eta_W`; sign pinned from W mass eigenvector with required sign test.
- `a_ref` now reuses shared `spectrum.a(DEFAULT_A_REF_C=0.65)` object.
- Stored fields aligned: `delta_g_W_ud_L/R`, `delta_g_W_lnu_L/R`, `charged_contact_LL`, `delta_G_F_over_G_F`.
- Exact roots specified via `solve_kk("gauge", "NN", exact=True)`, `x_1~2.4505`, zero mode included in `RSEWSpectrum` charged diagonalization.
- Kept K017 HARD `R_K`, K018 SOFT `K_l3 |Vus|`, EW002 SOFT, EW003 data-level+diagnostics.