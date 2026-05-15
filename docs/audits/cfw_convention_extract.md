# CFW 0804.1954 Convention Extract

Audit target: Csaki--Falkowski--Weiler, arXiv:0804.1954v2.  
Local source read: `/tmp/cfw/0804.1954.pdf` converted with `pdftotext` to
`/tmp/cfw/0804.1954.txt` on 2026-05-15.  The PDF has 43 pages.

## Headline numbers

CFW state in the abstract that the ordinary RS anarchic setup needs a KK-gluon
mass of about `21 TeV`, while the pseudo-Goldstone Higgs setup raises the
bound to about `33 TeV` (0804.1954, abstract, PDF p. 1).  In the RS section
they give an analytic estimate `M_G = (22 +/- 6) TeV`, then their scan says
that most generated points with `M_G < 21 TeV` violate the `Im C4K` bound
(0804.1954, section 2, eqs. 2.19--2.21 and Fig. 1 caption, PDF pp. 8--10).
For the pseudo-Goldstone scan, Fig. 7 and the text state that satisfying
`Lambda_Im C4K > 1.6e5 TeV` requires a KK scale of about `33 TeV`
(0804.1954, section 5.3 and Fig. 7 caption, PDF pp. 31--33).

## Effective operator basis

CFW integrate out the KK gluon and define the effective Hamiltonian with three
operator structures relevant for the quoted bounds:

- `C1`: color-singlet left-left vector current,
  `(bar q_L^i gamma_mu q_L^j)(bar q_L^k gamma_mu q_L^l)`.
- `C4`: scalar LR color-singlet contraction,
  `(bar q_R^i q_L^j)(bar q_L^k q_R^l)`.
- `C5`: scalar LR color-mixed contraction,
  `(bar q_R^i_alpha q_L^j_beta)(bar q_L^k_beta q_R^l_alpha)`.

The single-KK matching has the sign pattern `C4 = - g_R g_L / M_G^2` and
`C5 = + g_R g_L / (3 M_G^2)` in their Hamiltonian convention; Appendix A
gives the same tower-summed kaon relation as `C4K < 0` and `C5K = -C4K/3`
(0804.1954, section 2, eq. 2.18, PDF p. 7; Appendix A, eq. A.7, PDF
pp. 36--37).

Comparison to BMU and this repo:

- CFW's `C4/C5` are the conventional scalar-LR operators used by UTfit.
- The BMU LR map recorded in the Wilson audit is
  `Q1_LR^BMU = -2 O5_LR`, `Q2_LR^BMU = O4_LR`.
- The post-audit code uses the same scalar-LR signs as CFW:
  `C4_LR = -LR/M_KK^2`, `C5_LR = LR/(3 M_KK^2)`, then evolves the
  conventional `[O4_LR, O5_LR]` basis with the BMU-mapped LO ADM
  (see `docs/audits/wilson_rg_inventory.md`).

Conclusion: the CFW KK-gluon LR matching signs are compatible with the
post-hole-6 code basis.  The remaining comparison differences are numerical
inputs, coupling conventions, and sampling design, not a C4/C5 sign flip.

## c-pattern and scan design

For the ordinary RS scan behind Fig. 1, CFW randomly choose `1/R0` and brane
Yukawa couplings, then select 500 points that approximately reproduce quark
masses, CKM magnitudes, and the Jarlskog invariant.  The Fig. 1 caption gives
the third-generation localization choices:

- `c_q3 in [0.4, 0.45]`
- `c_u3 in [-0.3, -0.05]`
- `c_d3 = -0.55`
- `|Y_*| in [1, 3]`

(0804.1954, section 2 and Fig. 1 caption, PDF pp. 9--10).

For the pseudo-Goldstone Higgs scan, CFW sample the EWSB-linked parameters
directly rather than drawing fixed-c anarchic Yukawas.  They use
`c_q3 in [0.2, 0.48]`, `1/R0 in [900, 12500] GeV`, `c_-d3 = -0.55`,
`|mtilde_u| in [0.5, 5]`, `|mtilde_d| in [0.5, 2]`, `r in [0, 0.8]`, set
`Mtilde_u = Mtilde_d = 0`, and then solve for `c_u3` so that the Higgs
potential minimum occurs at the chosen `v/f_pi`.  They also require
`N_CFT >= 5` for calculability (0804.1954, section 5.3, PDF pp. 31--32).

This is an inverse-design scan: CFW first pick or solve the localization
variables that control the mass hierarchy, then reject mass/CKM outliers.
Our RUNA/RUNB scans instead fix one `c` pattern and forward-sample anarchic
Yukawa matrices.

## Yukawa convention

CFW write the IR-brane Yukawa term as `v/sqrt(2) * Ytilde` and obtain
`m_SM = v/sqrt(2) f_q Ytilde f_-u,-d` (0804.1954, section 2, eqs. 2.7--2.9,
PDF pp. 5--6).  With their usual electroweak convention `v = 246 GeV`, this
is numerically the same normalization as our `v_GeV = 174 GeV` mass matrix
`M = v_GeV f_Q Y f_f` in `scripts/run_rs_anarchy.py`.

Mapping: `Ytilde_CFW` maps directly to our `Y` up to the standard
`246/sqrt(2) = 174 GeV` convention.  The difference is the prior: CFW's RS
scan quotes `|Y_*| in [1, 3]`, while our default RUNA uses iid
`Re Y_ij, Im Y_ij ~ U(-1.5, 1.5)` with `|Y_ij| >= 0.1`; Run C uses a
larger floor intended to mimic the CFW prior but keeps our fixed-c design.

## Gauge-coupling convention

CFW define the dimensionless color bulk coupling `g_s^*`.  With no boundary
kinetic terms at tree level, they relate it to the four-dimensional QCD
coupling as `g_s^* = g_s(M_G) sqrt(log(R0/R))`, numerically about `6`
(0804.1954, section 2, eq. 2.20 and surrounding text, PDF pp. 8--9).

They then emphasize boundary-kinetic-term dependence:

- Their published `21 TeV` RS bound corresponds to a boundary-term choice in
  which the bare UV term cancels the one-loop running contribution
  (0804.1954, section 2, PDF p. 10).
- With no UV boundary kinetic term at the Planck scale, one-loop QCD running
  lowers `g_s^*` from about `6` to about `3`, reducing the RS bound from
  `21 TeV` to `10.5 TeV` (0804.1954, section 2, eq. 2.21 and PDF p. 10).
- A strongly coupled choice `g_s^* ~ 4 pi` raises the RS bound to about
  `42 TeV` (0804.1954, section 2, PDF p. 10).

For the pseudo-Goldstone Higgs section, the weak SO(5) coupling satisfies
`N_CFT = 16 pi^2/g_*^2`, and small brane kinetic terms imply `g_* ~ 4`
for the SM weak coupling.  The color coupling used for the Fig. 7 headline is
again `g_s^* ~ g_s sqrt(log(R0/R)) ~ 6`; they say the `33 TeV` estimate can be
relaxed by lowering `g_s^*` through QCD running, while a calculable framework
cannot bring it below about `10 TeV` (0804.1954, section 3, eqs. 3.1--3.2,
PDF p. 11; section 5.3, PDF p. 31).

## PDG-match tightness

CFW compare the generated mass eigenvalues, CKM magnitudes, and Jarlskog
invariant to target values at 3 TeV.  They accept a random matrix only when
the maximum relative distance of any constraint is below `30%`
(0804.1954, section 5.3, eqs. 5.26--5.27 and following paragraph,
PDF pp. 32--33).

As a multiplicative gate, a `30%` relative-distance cut is asymmetric:
`0.7 <= prediction/target <= 1.3`.  A symmetric log-factor gate that contains
that interval uses `factor = 1/(1 - 0.30) = 1.4286`.  It is often summarized
as a "factor 1.3" gate on the upper side, but for our stored residuals the
closest forward-only implementation is a symmetric factor-`1.43` cut on
masses, CKM elements, and `J`.

## Hadronic inputs and epsilon_K convention

CFW do not quote explicit `B_K`, `B_K^(4)`, `B_K^(5)`, decay constants, or
`epsilon_K^SM` values in the paper.  Instead, they use UTfit's
model-independent lower limits on `Lambda_F` for the `Delta F = 2`
operators and state that they corrected the UTfit Wilson coefficients from
`Lambda_F` to `3 TeV` using the RG formulae of Buras et al.
(0804.1954, section 2, Table 2, PDF pp. 7--8).  The load-bearing kaon entry is
`Im C4K`, with `Lambda_F = 160e3 TeV`.

Therefore the CFW hadronic convention is not a list of bag constants that can
be copied into `deltaf2.py`; it is an older UTfit black-box constraint table.
Our post-audit code instead uses explicit FLAG 2024 kaon bags
`B_K(2 GeV) = 0.5503`, `B_4^K(3 GeV) = 0.903`, `B_5^K(3 GeV) = 0.691`, and
BGS 2020 `epsilon_K^SM = 2.161e-3` (see
`docs/audits/bag_param_inventory.md` and
`docs/audits/epsilon_k_sm_decision.md`).

For convention matching, the only directly reproducible CFW-era
`epsilon_K` change in our pipeline is to switch the central NP budget back to
the legacy `epsilon_K^SM = 1.81e-3` value that was previously used in this
repo.  CFW itself does not provide an independent `epsilon_K^SM` number.
