Plan completed: read approved design/conventions, added one derivation-only note, validated scope.

Wrote: [derivations/rs_ew_gauge_kk_coupling.tex](/n/holylabs/randall_lab/Lab/obarrera/5D-Neutrino-Mixing/derivations/rs_ew_gauge_kk_coupling.tex:1)

Key pins: `Omega_n(c)=int_epsilon^1 dt w0(c,t) chi_n(t)` and `a_N(c)=sum_n (M_KK^2/m_n^2) chi_n(1) Omega_n(c)`, with `M_KK=x_1 Lambda_IR`, `x_1 ~= 2.4048`.

Normalization: `L_Z = g_Z Z_mu fbar gamma^mu(g_L P_L+g_R P_R)f`; stored `z_delta_g_*` are dimensionless additive shifts inside the parentheses, with no extra `g_Z`.

`a_ref` is subtracted after numerical overlap computation; universal `a(c)=a_ref` must give all stored shifts zero.

Phase 2 checks pinned: universal-c/all-`a_ref` gives all `delta g=0`; IR-localized `b_R` gives `a(c_bR)>a(c_lightR)`, `delta g_R^b<0`, about `1e-3` at `M_KK~3 TeV`.

Deferred: closed-form `a(c)`, fermion-KK mixing, custodial/BKT choices, loop dipoles, and Higgs-flavor matching.

Validation: no production/test files changed; `pdflatex -output-directory=/tmp` succeeds.

Refs: [Agashe-Perez-Soni](https://arxiv.org/abs/hep-ph/0408134), [Casagrande et al. 0807.4937](https://arxiv.org/abs/0807.4937), [Buras-Duling-Gori 0905.2318](https://arxiv.org/abs/0905.2318).

PHASE1-DONE.